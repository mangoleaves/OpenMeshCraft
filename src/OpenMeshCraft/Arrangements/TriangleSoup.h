#pragma once

#include "Utils.h"

#ifdef OMC_ARR_GLOBAL_POINT_SET
	#include "OpenMeshCraft/Utils/DisableWarnings.h"

	#include "parallel_hashmap/btree.h"

	#include "OpenMeshCraft/Utils/EnableWarnings.h"
#endif

#include <bitset>
#include <map>
#include <vector>

namespace OMC {

#ifdef OMC_ARR_TS_PARA
template <typename T>
using concurrent_vector = tbb::concurrent_vector<T>;
#else // serial, std::vector is more friendly for debug
template <typename T>
using concurrent_vector = std::vector<T>;
#endif

template <typename Traits>
class TriangleSoup
{
public: /* Types **************************************************************/
	using NT         = typename Traits::NT;
	using Bbox       = typename Traits::BoundingBox;
	using EPoint     = typename Traits::EPoint;
	using GPoint     = typename Traits::GPoint;
	using IPoint_SSI = typename Traits::IPoint_SSI;
	using IPoint_LPI = typename Traits::IPoint_LPI;
	using IPoint_TPI = typename Traits::IPoint_TPI;
	using AsGP       = typename Traits::AsGP;
	using AsEP       = typename Traits::AsEP;
	using ToEP       = typename Traits::ToEP;

	using CalcBbox           = typename Traits::CalcBbox;
	using OrientOn2D         = typename Traits::OrientOn2D;
	using LessThan3D         = typename Traits::LessThan3D;
	using LongestAxis        = typename Traits::LongestAxis;
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
	using EdgeMap        = std::vector<phmap::flat_hash_map<Edge, index_t>>;
	/// mutex for reading and writing EdgeMap.
	using EdgeMapMutexes = std::vector<tbb::spin_mutex>;

	/* ----- [c]oplana[r]/[c]olinea[r] (CCr) edges related structures ----- */

	struct CCrEdgeInfo;

	/* ----- seg2tris related structures ----- */

	/// <smaller vertex index, larger vertex index>
	using Segment       = UIPair;
	/// map a segment to its related triangles.
	/// 1. map Seg to (Seg.first+Seg.second) % Seg2Tris.size() to locate the
	/// flat_hash_map in the outer std::vector.
	/// 2. map Seg to Tris in the inner flat_hash_map.
	using SegMap        = std::vector<phmap::flat_hash_map<Segment, index_t>>;
	/// mutex for reading and writing Seg2Tris.
	using SegMapMutexes = std::vector<tbb::spin_mutex>;

	/* ----- edge2pts related structures ----- */

#ifndef OMC_ARR_GLOBAL_POINT_SET
	struct EdgeComparator;
	using Edge2PntsSet =
	  boost::container::flat_set<index_t, EdgeComparator, AuxVector8<index_t>>;
#else
	using Edge2PntsSet = concurrent_vector<index_t>;
#endif

	/* ----- seg2pts related structures ----- */

#ifndef OMC_ARR_GLOBAL_POINT_SET
	struct SegComparator;
	using Seg2PntsSet =
	  boost::container::flat_set<index_t, SegComparator, AuxVector8<index_t>>;
#else
	using Seg2PntsSet = concurrent_vector<index_t>;
#endif

	/* ----- coplnar pockes related structures ----- */

	using PocketsMap = phmap::flat_hash_map<std::vector<index_t>, index_t>;

public: /* Constructors *******************************************************/
	TriangleSoup() = default;

	~TriangleSoup() = default;

	void initialize();

public:
	/***** Below data should be set by user ******/

	/// implement vertices
	concurrent_vector<GPoint *>                  vertices;
	/// implement indices of vertices (used in v_map)
	tbb::concurrent_vector<std::atomic<index_t>> indices;

	/// triangles
	std::vector<index_t> triangles;
	/// labels
	std::vector<Label>   tri_labels;

	/// all explicit and implicit points
	std::vector<PntArena> *pnt_arenas = nullptr;

public:
	/***** Below data should be initialized by calling initialize() ******/

	size_t num_orig_vtxs;
	size_t num_orig_tris;

	concurrent_vector<Edge> edges;
	/// map an edge(unique index pair) to a unique edge index.
	/// 1. outer map: map smaller vertex index to a flat hash map container by
	/// std::vector with limited length to avoid collision in parallel context.
	/// 2. inner map: map edge to its index by flat hash map.
	EdgeMap                 edge_map;
	/// mutex for reading and writing EdgeMap.
	EdgeMapMutexes          edge_mutexes;
	/// edge map size
	const static size_t     edge_map_size = 64;
	/// triangle edges
	std::vector<index_t>    tri_edges;

public: /* Size queries *******************************************************/
	size_t numVerts() const { return vertices.size(); }
	size_t numEdges() const { return edges.size(); }
	size_t numTris() const { return triangles.size() / 3; }

	size_t numOrigVerts() const { return num_orig_vtxs; }
	size_t numOrigTris() const { return num_orig_tris; }

public: /* Vertices ***********************************************************/
	const GPoint &vert(index_t v_id) const;

	const NT *vertPtr(index_t v_id) const;

	index_t addImplVert(GPoint *pp);

public: /* Edges **************************************************************/
	Edge edge(index_t e_id) const;

	index_t getOrAddEdge(index_t v0_id, index_t v1_id);

public: /* Triangles **********************************************************/
	const index_t *tri(index_t t_id) const;

	index_t triVertID(index_t t_id, index_t off) const;

	const GPoint &triVert(index_t t_id, index_t off) const;

	const NT *triVertPtr(index_t t_id, index_t off) const;

	index_t triEdgeID(index_t t_id, index_t off);

	index_t triEdgeID(index_t t_id, index_t off, std::true_type add_if_not_found);

	Label triLabel(index_t t_id) const;

public:
	/***** Below data are calculated by arrangements  ******/

	// does triangle have intersections
	std::vector<uint8_t> tri_has_intersections;

	// coplanar triangle
	std::vector<concurrent_vector<index_t>> coplanar_tris;

	// coplanar edge
	std::vector<concurrent_vector<CCrEdgeInfo>> coplanar_edges;

	// colinear edge
	concurrent_vector<concurrent_vector<CCrEdgeInfo>> colinear_edges;

	// store intersection points on triangles
	std::vector<concurrent_vector<index_t>> tri2pts;

	// store intersection points on edge
	concurrent_vector<Edge2PntsSet>         edge2pts;
	// mutexes for edge2pts (NOTE: must be tbb::concurrent_vector)
	tbb::concurrent_vector<tbb::spin_mutex> edge2pts_mutex;

	// all unique constarined segments
	concurrent_vector<Segment> segments;
	// map segments to their ID
	SegMap                     seg_map;
	// mutexes for seg_map
	SegMapMutexes              seg_mutexes;

	// store contrained segments on triangles
	std::vector<concurrent_vector<index_t>> tri2segs;

	// reverse map contrained segments to triangles on where they locate
	std::vector<concurrent_vector<index_t>> seg2tris;

	// points stored and sorted on segments
	std::vector<Seg2PntsSet>                seg2pts;
	// mutexes for seg2pts
	tbb::concurrent_vector<tbb::spin_mutex> seg2pts_mutex;

	/// orthogonal plane of a triangle
	/// (only the triangles with intersections will be calculated)
	std::vector<Plane> tri_plane;

	/// orientation of a triangle on its orthogonal plane
	/// (only the triangles with intersections will be calculated)
	std::vector<Sign> tri_orient;

	// coplanar pockets (// TODO pocket size? inlined_vector?)
	PocketsMap pockets_map;

	PocketsMap pockets_map_with_tpi;

#ifdef OMC_ARR_GLOBAL_POINT_SET
	struct AuxPoint
	{
		const GPoint *pnt;

		// clang-format off
		AuxPoint(const GPoint *_p) : pnt(_p) {}

		bool operator<(const AuxPoint &rhs) const { return LessThan3D()(*pnt, *rhs.pnt) == Sign::NEGATIVE; }
		bool operator==(const AuxPoint &rhs) const { return LessThan3D()(*pnt, *rhs.pnt) == Sign::ZERO; }
		bool operator>(const AuxPoint &rhs) const { return LessThan3D()(*pnt, *rhs.pnt) == Sign::POSITIVE; }
		// clang-format on
	};

	phmap::btree_map<AuxPoint, index_t> global_point_set;
	tbb::spin_mutex                     new_uniq_point_mutex;
#endif

public:
	// flags
	bool any_index_fixed;

	// mutexes
	tbb::spin_mutex new_vertex_mutex;
	tbb::spin_mutex new_edge_mutex;
	tbb::spin_mutex new_segment_mutex;
	tbb::spin_mutex new_tris_mutex;

#ifdef OMC_ARR_AUX_LPI
	// jolly points
	std::vector<GPoint *> jolly_points;
#endif

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

	index_t getOrAddSegment(const Segment &seg);

	index_t segmentID(const Segment &seg) const;

	Segment segment(index_t seg_id) const;

	void addSegmentInTriangle(index_t t_id, index_t seg_id);

	void addTrianglesInSegment(index_t seg_id, index_t t_id);

	/* -- vertex in segment -- */

	tbb::spin_mutex &getS2PMutex(index_t seg_id);

	index_t findVertexInSeg(index_t seg_id, const GPoint &pnt) const;

	void addVertexInSeg(index_t seg_id, index_t v_id);

	void fixVertexInSeg(index_t seg_id, index_t old_vid, index_t new_vid);

	/* Unique pocket */

	index_t addVisitedPolygonPocket(const std::vector<index_t> &polygon,
	                                index_t pos, bool with_tpi);

#ifdef OMC_ARR_GLOBAL_POINT_SET
	std::pair<index_t, bool> addUniquePoint(GPoint &pnt);
#endif

public: /* Query ***********************************************************/
	/* Coplanar and colinear */

	const concurrent_vector<index_t> &coplanarTriangles(index_t t_id) const;

	const concurrent_vector<CCrEdgeInfo> &coplanarEdges(index_t t_id) const;

	const concurrent_vector<CCrEdgeInfo> &colinearEdges(index_t e_id) const;

	/* Has intersection */

	bool triangleHasIntersections(index_t t_id) const;

	/* Intersection points and contrained segments */

	const concurrent_vector<index_t> &trianglePointsList(index_t t_id) const;

	const Edge2PntsSet &edgePointsList(index_t e_id) const;

	const concurrent_vector<index_t> &triangleSegmentsList(index_t t_id) const;

	const concurrent_vector<index_t> &segmentTrianglesList(index_t seg_id) const;

	/* Coplanar pockets */

	PocketsMap &pocketsMap() { return pockets_map; }

	PocketsMap &pocketsMapWithTPI() { return pockets_map_with_tpi; }

	/* Orthogonal plane and orientation on the orthogonal plane */

	Plane triPlane(index_t t_id) const;

	Sign triOrient(index_t t_id) const;

public: /* Modify *********************************************************/
	void removeDuplicatesBeforeFix();

	void fixColinearEdgesIntersections();

	void fixAllIndices();

	void addEndPointsToE2P();

	void calcOrthogonalPlane();

	void calcTriangleOrient();
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "TriangleSoup.inl"
#endif