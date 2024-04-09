#pragma once

#include "AuxStructure.h"
#include "Utils.h"

// clang-format off
#include "OpenMeshCraft/Utils/DisableWarnings.h"
#include "parallel_hashmap/phmap.h"
#include "OpenMeshCraft/Utils/EnableWarnings.h"
// clang-format on

#include <bitset>
#include <map>
#include <vector>

namespace OMC {

template <typename Traits>
class TriangleSoup
{
public: /* Types **************************************************************/
	using NT         = typename Traits::NT;
	using EPoint     = typename Traits::EPoint;
	using GPoint     = typename Traits::GPoint;
	using IPoint_SSI = typename Traits::IPoint_SSI;
	using IPoint_LPI = typename Traits::IPoint_LPI;
	using IPoint_TPI = typename Traits::IPoint_TPI;

	using OrientOn2D         = typename Traits::OrientOn2D;
	using LessThan3D         = typename Traits::LessThan3D;
	using CollinearPoints3D  = typename Traits::CollinearPoints3D;
	using MaxCompInTriNormal = typename Traits::MaxCompInTriNormal;

	using PntArena = PointArena<Traits>;

	using Tree = Arr_Tree_Intersection<Traits>;

	/* ----- Edge related structures ----- */

	/// <smaller vertex index, larger vertex index>
	using Edge           = UIPair;
	/// map an edge(index pair) to a unique edge index.
	/// 1. map smaller vertex index by std::vector.
	/// 2. map larger vertex index by flat hash map.
	using EdgeMap        = std::vector<phmap::flat_hash_map<index_t, index_t>>;
	/// mutex for reading and writing EdgeMap.
	using EdgeMapMutexes = std::vector<tbb::spin_mutex>;

	/* ----- [c]oplana[r]/[c]olinea[r] (CCr) edges related structures ----- */

	struct CCrEdgeInfo;

	/* ----- seg2tris related structures ----- */

	/// <smaller vertex index, larger vertex index>
	using Seg      = UIPair;
	/// map a segment to its related triangles.
	/// 1. map Seg to (Seg.first+Seg.second) % Seg2Tris.size() to locate the
	/// flat_hash_map in the outer std::vector.
	/// 2. map Seg to Tris in the inner flat_hash_map.
	using Seg2Tris = std::vector<phmap::flat_hash_map<Seg, std::vector<index_t>>>;
	/// mutex for reading and writing Seg2Tris.
	using SegMutexes = std::vector<tbb::spin_mutex>;

	/* ----- edge2pts related structures ----- */

	struct EdgeComparator;
	using Edge2PntsSet = std::set<index_t, EdgeComparator>;

public: /* Constructors *******************************************************/
	TriangleSoup() = default;

	~TriangleSoup() = default;

	void initialize();

public:
	/***** Below data should be set by user ******/

	/// implement vertices
	tbb::concurrent_vector<GPoint *>               vertices;
	/// implement indices of vertices (used in v_map)
	tbb::concurrent_vector<std::atomic<index_t> *> indices;
	// TODO make use of indices, remove IdxArena if possible.

	/// triangles
	std::vector<index_t> triangles;
	/// labels
	std::vector<Label>   tri_labels;

	/// all explicit and implicit points
	std::vector<PntArena> *pnt_arenas = nullptr;
	/// local indices
	std::vector<IdxArena> *idx_arenas = nullptr;

protected:
	/***** Below data should be initialized by calling initialize() ******/

	size_t num_orig_vtxs;
	size_t num_orig_tris;

	tbb::concurrent_vector<UIPair> edges;
	/// map an edge(index pair) to a unique edge index.
	/// 1. map smaller vertex index by std::vector.
	/// 2. map larger vertex index by flat hash map.
	EdgeMap                        edge_map;
	/// mutex for reading and writing EdgeMap.
	EdgeMapMutexes                 edge_mutexes;
	/// triangle edges
	std::vector<index_t>           tri_edges;

public: /* Size queries *******************************************************/
	size_t numVerts() const { return vertices.size(); }
	size_t numEdges() const { return edges.size(); }
	size_t numTris() const { return triangles.size() / 3; }

	size_t numOrigVerts() const { return num_orig_vtxs; }
	size_t numOrigTris() const { return num_orig_tris; }

public: /* Vertices ***********************************************************/
	const GPoint &vert(index_t v_id) const;

	const NT *vertPtr(index_t v_id) const;

	index_t addImplVert(GPoint *pp, std::atomic<index_t> *ip);

public: /* Edges **************************************************************/
	const UIPair &edge(index_t e_id) const;

	index_t getOrAddEdge(index_t v0_id, index_t v1_id);

public: /* Triangles **********************************************************/
	const index_t *tri(index_t t_id) const;

	index_t triVertID(index_t t_id, index_t off) const;

	const GPoint &triVert(index_t t_id, index_t off) const;

	const NT *triVertPtr(index_t t_id, index_t off) const;

	index_t triEdgeID(index_t t_id, index_t off);

	index_t triEdgeID(index_t t_id, index_t off, std::true_type add_if_not_found);

	Label triLabel(index_t t_id) const;

protected:
	/***** Below data are calculated by arrangements  ******/

	bool any_index_fixed;

	// does triangle have intersections
	std::vector<uint8_t> tri_has_intersections;

	// coplanar triangle
	std::vector<tbb::concurrent_vector<index_t>> coplanar_tris;

	// coplanar edge
	std::vector<tbb::concurrent_vector<CCrEdgeInfo>> coplanar_edges;

	// colinear edge
	tbb::concurrent_vector<tbb::concurrent_vector<CCrEdgeInfo>> colinear_edges;

	// map point coordinates to vertex id
	std::unique_ptr<AuxPointMap<Traits>> v_map;

	// store intersection points on triangles
	std::vector<tbb::concurrent_vector<index_t>> tri2pts;

	// store intersection points on edge
	tbb::concurrent_vector<Edge2PntsSet>    edge2pts;
	tbb::concurrent_vector<tbb::spin_mutex> edge2pts_mutex;

	// store contrained segments on triangles
	std::vector<tbb::concurrent_vector<UIPair>> tri2segs;

	// reverse map contrained segments to triangles on where they locate
	Seg2Tris   seg2tris;
	SegMutexes seg_mutexes;

	/// orthogonal plane of a triangle
	/// (only the triangles with intersections will be calculated)
	std::vector<Plane> tri_plane;

	/// orientation of a triangle on its orthogonal plane
	/// (only the triangles with intersections will be calculated)
	std::vector<Sign> tri_orient;

	// coplanar pockets
	phmap::flat_hash_map<std::vector<index_t>, index_t> pockets_map;

public:
	// mutexes
	tbb::spin_mutex new_vertex_mutex;
	tbb::spin_mutex new_edge_mutex;
	tbb::spin_mutex new_tris_mutex;

public: /* Add **************************************************************/
	/* Coplanar and colinear */

	void addCoplanarTriangles(index_t ta, index_t tb);

	void addCoplanarEdge(index_t t_id, index_t e_id, index_t v0_id,
	                     index_t v1_id);

	void addColinearEdge(index_t e0_id, index_t e1_id, index_t v0_id,
	                     index_t v1_id);

	/* Has intersection */

	void setTriangleHasIntersections(index_t tA_id, index_t tB_id);

	/* Intersection points and contrained segments */

	/* -- vertex in triangle -- */

	void addVertexInTriangle(index_t t_id, index_t v_id);

	/* -- vertex in edge -- */

	tbb::spin_mutex &getE2PMutex(index_t e_id);

	index_t findVertexInEdge(index_t e_id, const GPoint &pnt) const;

	void addVertexInEdge(index_t e_id, index_t v_id);

	void fixVertexInEdge(index_t e_id, index_t old_vid, index_t new_vid);

	/* -- segment in triangle -- */

	void addSegmentInTriangle(index_t t_id, const UIPair &seg);

	void addTrianglesInSegment(const UIPair &seg, index_t t_id);

	void addTrianglesInSegment(const UIPair &seg, index_t tA_id, index_t tB_id);

	/* Map point to unique index */

	template <typename GetIndex>
	std::pair<index_t, bool> addVertexInSortedList(const GPoint         *pp,
	                                               std::atomic<index_t> *ip,
	                                               GetIndex              get_idx);

	/* Unique pocket */

	index_t addVisitedPolygonPocket(const std::vector<index_t> &polygon,
	                                index_t                     pos);

public: /* Query ***********************************************************/
	/* Coplanar and colinear */

	const tbb::concurrent_vector<index_t> &coplanarTriangles(index_t t_id) const;

	const tbb::concurrent_vector<CCrEdgeInfo> &coplanarEdges(index_t t_id) const;

	const tbb::concurrent_vector<CCrEdgeInfo> &colinearEdges(index_t e_id) const;

	/* Has intersection */

	bool triangleHasIntersections(index_t t_id) const;

	/* Intersection points and contrained segments */

	const tbb::concurrent_vector<index_t> &trianglePointsList(index_t t_id) const;

	const Edge2PntsSet &edgePointsList(index_t e_id) const;

	const tbb::concurrent_vector<UIPair> &
	triangleSegmentsList(index_t t_id) const;

	const std::vector<index_t> &segmentTrianglesList(const UIPair &seg) const;

	/* Orthogonal plane and orientation on the orthogonal plane */

	Plane triPlane(index_t t_id) const;

	Sign triOrient(index_t t_id) const;

public: /* Modify *********************************************************/
	void buildVMap(Tree *tree);

	void fixAllIndices();

	void removeAllDuplicates();

	void addEndPointsToE2P();

	void calcPlaneAndOrient();
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "TriangleSoup.inl"
#endif