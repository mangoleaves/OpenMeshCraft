#pragma once

#include "Utils.h"
#include "AuxStructure.h"

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
	using NT     = typename Traits::NT;
	using GPoint = typename Traits::GPoint;

	using PntArena = PointArena<Traits>;

	using Tree = Arr_Tree_Intersection<Traits>;

	/// <smaller vertex index, larger vertex index>
	using Edge           = UIPair;
	/// map an edge(index pair) to a unique edge index.
	/// 1. map smaller vertex index by std::vector.
	/// 2. map larger vertex index by flat hash map.
	using EdgeMap        = std::vector<phmap::flat_hash_map<index_t, index_t>>;
	/// mutex for reading and writing EdgeMap.
	using EdgeMapMutexes = std::vector<tbb::spin_mutex>;

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

	/// triangles
	std::vector<index_t> triangles;
	/// labels
	std::vector<Label>   tri_labels;

	/// all explicit and implicit points
	std::vector<PntArena> *pnt_arenas = nullptr;
	/// local indices
	std::vector<IdxArena> *idx_arenas = nullptr;

	/***** Below data should be initialized by calling initialize() ******/

	size_t num_orig_vtxs;
	size_t num_orig_tris;

	std::atomic<index_t> largest_edge_id;
	/// map an edge(index pair) to a unique edge index.
	/// 1. map smaller vertex index by std::vector.
	/// 2. map larger vertex index by flat hash map.
	EdgeMap              edge_map;
	/// mutex for reading and writing EdgeMap.
	EdgeMapMutexes       edge_mutexes;
	/// triangle edges
	std::vector<index_t> tri_edges;

public: /* Size queries *******************************************************/
	size_t numVerts() const { return vertices.size(); }
	size_t numTris() const { return triangles.size() / 3; }

	size_t numOrigVerts() const { return num_orig_vtxs; }
	size_t numOrigTris() const { return num_orig_tris; }

public: /* Vertices ***********************************************************/
	const GPoint &vert(index_t v_id) const;

	const NT *vertPtr(index_t v_id) const;

	index_t addImplVert(GPoint *pp, std::atomic<index_t> *ip);

public: /* Edges **************************************************************/
	index_t getOrAddEdge(index_t v0_id, index_t v1_id);

public: /* Triangles **********************************************************/
	const index_t *tri(index_t t_id) const;

	index_t triVertID(index_t t_id, index_t off) const;

	const GPoint &triVert(index_t t_id, index_t off) const;

	const NT *triVertPtr(index_t t_id, index_t off) const;

	index_t triEdgeID(index_t t_id, index_t off) const;

	Label triLabel(index_t t_id) const;

public:
	/***** Below data are calculated by arrangements  ******/

	// FIXME clear duplication

	std::vector<uint8_t> tri_has_intersections;

	// coplanar related data
	std::vector<tbb::concurrent_vector<index_t>> coplanar_tris;
	std::vector<tbb::concurrent_vector<index_t>> coplanar_edges;

	// map point coordinates to vertex id
	std::unique_ptr<AuxPointMap<Traits>> v_map;

	// store intersection points on triangles
	std::vector<tbb::concurrent_vector<index_t>> tri2pts;

	// store intersection points on edge
	tbb::concurrent_vector<tbb::concurrent_vector<index_t>> edge2pts;

	// store contrained segments on triangles
	std::vector<tbb::concurrent_vector<UIPair>> tri2segs;

	// reverse map contrained segments to triangles on where they locate
	phmap::parallel_flat_hash_map<
	  /*Key*/ UIPair,
	  /*Value*/ tbb::concurrent_vector<index_t>,
	  /*Hash*/ phmap::priv::hash_default_hash<UIPair>,
	  /*Eq*/ phmap::priv::hash_default_eq<UIPair>,
	  /*Alloc*/
	  phmap::priv::Allocator<
	    phmap::priv::Pair<const UIPair, tbb::concurrent_vector<index_t>>>,
	  /*2**N submaps*/ 4,
	  /*Mutex*/ tbb::spin_mutex>
	  seg2tris;

	// coplanar pockets
	phmap::parallel_flat_hash_map<
	  /*Key*/ std::vector<index_t>,
	  /*Value*/ index_t,
	  /*Hash*/ phmap::priv::hash_default_hash<std::vector<index_t>>,
	  /*Eq*/ phmap::priv::hash_default_eq<std::vector<index_t>>,
	  /*Alloc*/
	  phmap::priv::Allocator<
	    phmap::priv::Pair<const std::vector<index_t>, index_t>>,
	  /*2**N submaps*/ 4,
	  /*Mutex*/ tbb::spin_mutex>
	  pockets_map;

	// mutexes
	tbb::spin_mutex new_vertex_mutex;
	tbb::spin_mutex new_tris_mutex;

public: /* Add **************************************************************/
	/* Coplanar */

	void addCoplanarTriangles(index_t ta, index_t tb);

	void addCoplanarEdge(index_t t_id, index_t e_id);

	/* Has intersection */

	void setTriangleHasIntersections(index_t t_id);

	/* Intersection points and contrained segments */

	void addVertexInTriangle(index_t t_id, index_t v_id);

	void addVertexInEdge(index_t e_id, index_t v_id);

	void addSegmentInTriangle(index_t t_id, const UIPair &seg);

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
	/* Coplanar */

	const tbb::concurrent_vector<index_t> &coplanarTriangles(index_t t_id) const;

	const tbb::concurrent_vector<index_t> &coplanarEdges(index_t t_id) const;

	/* Has intersection */

	bool triangleHasIntersections(index_t t_id) const;

	/* Intersection points and contrained segments */

	const tbb::concurrent_vector<index_t> &trianglePointsList(index_t t_id) const;

	const tbb::concurrent_vector<index_t> &edgePointsList(index_t e_id) const;

	const tbb::concurrent_vector<UIPair> &
	triangleSegmentsList(index_t t_id) const;

	const tbb::concurrent_vector<index_t> &
	segmentTrianglesList(const UIPair &seg) const;

public: /* Modify *********************************************************/
	void build_vmap(Tree *tree);

	void remove_all_duplicates();

	template <typename Comparator>
	void sortEdgeList(index_t eid, Comparator comp);
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "TriangleSoup.inl"
#endif