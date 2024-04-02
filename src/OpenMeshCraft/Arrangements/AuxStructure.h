#pragma once

#include "TriangleSoup.h"

// clang-format off
#include "OpenMeshCraft/Utils/DisableWarnings.h"
#include "parallel_hashmap/btree.h"
#include "OpenMeshCraft/Utils/EnableWarnings.h"
// clang-format on

#include <mutex>
#include <set>

namespace OMC {

/// @brief Wrap of GPoint*
template <typename Traits>
class AuxPoint
{
public:
	using GPoint     = typename Traits::GPoint;
	using LessThan3D = typename Traits::LessThan3D;

	const GPoint *pt;
	AuxPoint(const GPoint *pt_)
	  : pt(pt_)
	{
	}

	bool operator<(const AuxPoint &rhs) const
	{
		return LessThan3D()(*pt, *(rhs.pt)) == Sign::NEGATIVE;
	}
	bool operator==(const AuxPoint &rhs) const
	{
		return LessThan3D().coincident(*pt, *(rhs.pt));
	}
};

/// @brief Map GPoint to a unique index
template <typename Traits>
class AuxPointMap
{
public:
	using GPoint = typename Traits::GPoint;

public:
	AuxPointMap() = default;
	virtual ~AuxPointMap() {}

public:
	/// @brief Insert point into map.
	/// @param item pair of geometry implementation and input index
	/// @return pair of index (input index if the point is inserted in map first
	/// time, otherwise the previous saved index) and the result of the insert
	/// operation (true if succeed to insert otherwise false).
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) = 0;
};

template <typename Traits>
class AuxPointMap_BTree : public AuxPointMap<Traits>
{
public:
	using GPoint = typename Traits::GPoint;

public:
	AuxPointMap_BTree() = default;
	virtual ~AuxPointMap_BTree() {}

public:
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) override final;

protected:
	phmap::btree_map<AuxPoint<Traits>, std::atomic<index_t> *> map;
};

template <typename Traits>
class AuxPointMap_ConcurrentMap : public AuxPointMap<Traits>
{
public:
	using GPoint = typename Traits::GPoint;

public:
	AuxPointMap_ConcurrentMap() = default;
	virtual ~AuxPointMap_ConcurrentMap() {}

public:
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) override final;

protected:
	tbb::concurrent_map<AuxPoint<Traits>, std::atomic<index_t> *> map;
};

// Forward declaration
template <typename Traits>
class Arr_Tree_Intersection;

template <typename Traits>
class AuxPointMap_Tree : public AuxPointMap<Traits>
{
public:
	using GPoint = typename Traits::GPoint;
	using Tree   = Arr_Tree_Intersection<Traits>;

public:
	AuxPointMap_Tree(Tree *_tree)
	  : tree(_tree)
	{
	}
	virtual ~AuxPointMap_Tree() {}

public:
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) override final;

	template <typename GetIndex>
	std::pair<index_t, bool> insert_F(const GPoint *pp, std::atomic<index_t> *idx,
	                                  GetIndex get_idx);

private:
	Tree *tree;
};

/// @brief Auxiliary structure, storing intersection information for triangle
/// soup.
template <typename Traits>
class AuxiliaryStructure
{
private:
	using GPoint = typename Traits::GPoint;
	using Tree   = Arr_Tree_Intersection<Traits>;

public: /* Constructors *****************************************************/
	AuxiliaryStructure() = default;

	void initialize(const TriangleSoup<Traits> &ts);

	void build_vmap(const TriangleSoup<Traits> &ts, Tree *tree);

	AuxiliaryStructure(AuxiliaryStructure &&rhs) noexcept;

	~AuxiliaryStructure() {}

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
	                                size_t                      pos);

public: /* Query ***********************************************************/
	/* Coplanar */

	const std::vector<size_t> &coplanarTriangles(index_t t_id) const;

	const std::vector<size_t> &sortedCoplanarTriangles(index_t t_id);

	const std::vector<size_t> &coplanarEdges(index_t t_id) const;

	bool triangleHasCoplanars(index_t t_id) const;

	/* Has intersection */

	bool triangleHasIntersections(index_t t_id) const;

	/* Intersection points and contrained segments */

	const std::vector<size_t> &trianglePointsList(index_t t_id) const;

	const std::vector<size_t> &edgePointsList(index_t e_id) const;

	const std::vector<UIPair> &triangleSegmentsList(index_t t_id) const;

	const std::vector<size_t> &segmentTrianglesList(const UIPair &seg) const;

public: /* Modify *********************************************************/
	template <typename Comparator>
	void sortEdgeList(index_t eid, Comparator comp);

public:
	// intersections
	std::vector<UIPair>                                 intersection_list;
	// triangles
	std::vector<uint8_t>                                tri_has_intersections;
	std::vector<std::vector<index_t>>                   coplanar_tris;
	std::vector<std::vector<index_t>>                   coplanar_edges;
	std::vector<uint8_t>                                coplanar_tris_sorted;
	// map point coordinates to vertex id
	std::unique_ptr<AuxPointMap<Traits>>                v_map;
	// store intersection points on triangles and edges
	std::vector<std::vector<index_t>>                   tri2pts;
	std::vector<std::vector<index_t>>                   edge2pts;
	// store contrained segments on triangles
	std::vector<std::vector<UIPair>>                    tri2segs;
	// reverse map contrained segments to triangles on where they locate
	phmap::flat_hash_map<UIPair, std::vector<index_t>>  seg2tris;
	// pockets
	phmap::flat_hash_set<std::vector<index_t>>          visited_pockets;
	phmap::flat_hash_map<std::vector<index_t>, index_t> pockets_map;

	// mutexes
	tbb::spin_mutex new_vertex_mutex;
	tbb::spin_mutex new_tris_mutex;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "AuxStructure.inl"
#endif