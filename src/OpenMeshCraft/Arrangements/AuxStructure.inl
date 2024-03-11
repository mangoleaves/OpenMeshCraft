#pragma once

#include "AuxStructure.h"

namespace OMC {

template <typename Traits>
std::pair<std::atomic<index_t> *, bool>
AuxPointMap_BTree<Traits>::insert(const GPoint *pp, std::atomic<index_t> *idx)
{
	auto ins = map.insert({AuxPoint<Traits>(pp), idx});
	return std::make_pair(
	  // the position of v (pos if first time, or the existed position otherwise)
	  ins.first->second,
	  // the result of the insert operation (true or false)
	  ins.second);
}

template <typename Traits>
std::pair<std::atomic<index_t> *, bool>
AuxPointMap_ConcurrentMap<Traits>::insert(const GPoint         *pp,
                                          std::atomic<index_t> *idx)
{
	auto ins = map.insert({AuxPoint<Traits>(pp), idx});
	return std::make_pair(
	  // the position of v (pos if first time, or the existed position otherwise)
	  ins.first->second,
	  // the result of the insert operation (true or false)
	  ins.second);
}

template <typename Traits>
std::pair<std::atomic<index_t> *, bool>
AuxPointMap_OcTree<Traits>::insert(const GPoint *pp, std::atomic<index_t> *idx)
{
	return tree->insert_point(pp, idx);
}

template <typename Traits>
template <typename GetIndex>
std::pair<index_t, bool> AuxPointMap_OcTree<Traits>::insert_F(
  const GPoint *pp, std::atomic<index_t> *idx, GetIndex get_idx)
{
	return tree->insert_point_F(pp, idx, get_idx);
}

template <typename Traits>
void AuxiliaryStructure<Traits>::initialize(const TriangleSoup<Traits> &ts)
{
	coplanar_tris.resize(ts.numTris());
	coplanar_edges.resize(ts.numTris());
	tri_has_intersections.resize(ts.numTris(), false);
	coplanar_tris_sorted.resize(ts.numTris(), false);

	tri2pts.resize(ts.numTris());
	edge2pts.resize(ts.numEdges());
	tri2segs.resize(ts.numTris());
}

template <typename Traits>
void AuxiliaryStructure<Traits>::build_vmap(const TriangleSoup<Traits> &ts,
                                            Tree                       *tree)
{
	tree->clear_points();
	v_map = std::make_unique<AuxPointMap_OcTree<Traits>>(tree);
	tbb::parallel_for(index_t(0), ts.numVerts(),
	                  [this, &ts](index_t v_id)
	                  { v_map->insert(ts.vertices[v_id], ts.indices[v_id]); });
}

template <typename Traits>
AuxiliaryStructure<Traits>::AuxiliaryStructure(
  AuxiliaryStructure &&rhs) noexcept
{
	intersection_list     = std::move(rhs.intersection_list);
	tri_has_intersections = std::move(rhs.tri_has_intersections);
	coplanar_tris_sorted  = std::move(rhs.coplanar_tris_sorted);
	coplanar_tris         = std::move(rhs.coplanar_tris);
	coplanar_edges        = std::move(rhs.coplanar_edges);
	v_map                 = std::move(rhs.v_map);
	tri2pts               = std::move(rhs.tri2pts);
	edge2pts              = std::move(rhs.edge2pts);
	tri2segs              = std::move(rhs.tri2segs);
	seg2tris              = std::move(rhs.seg2tris);
	visited_pockets       = std::move(rhs.visited_pockets);
	pockets_map           = std::move(rhs.pockets_map);
}

template <typename Traits>
void AuxiliaryStructure<Traits>::addVertexInTriangle(index_t t_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2pts.size(), "out of range");
	auto &points = tri2pts[t_id];
	if (contains(points, v_id))
		return;
	if (points.empty())
		points.reserve(8);
	points.push_back(v_id);
}

template <typename Traits>
void AuxiliaryStructure<Traits>::addVertexInEdge(index_t e_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	auto &points = edge2pts[e_id];
	if (contains(points, v_id))
		return;
	if (points.empty())
		points.reserve(8);
	points.push_back(v_id);
}

template <typename Traits>
void AuxiliaryStructure<Traits>::addSegmentInTriangle(index_t       t_id,
                                                      const UIPair &seg)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2segs.size(), "out of range");
	OMC_EXPENSIVE_ASSERT(isUnique(seg), "segment is not unique");
	auto &segments = tri2segs[t_id];
	if (segments.empty())
		segments.reserve(8);
	if (contains(segments, seg))
		return;
	segments.push_back(seg);
}

template <typename Traits>
void AuxiliaryStructure<Traits>::addTrianglesInSegment(const UIPair &seg,
                                                       index_t       tA_id,
                                                       index_t       tB_id)
{
	OMC_EXPENSIVE_ASSERT(isUnique(seg), "segment is not unique");
	auto &tris = seg2tris[seg];
	if (tris.empty())
		tris.reserve(4);
	if (tA_id == tB_id)
	{
		if (!contains(tris, tA_id))
			tris.push_back(tA_id);
	}
	else
	{
		if (!contains(tris, tA_id))
			tris.push_back(tA_id);
		if (!contains(tris, tB_id))
			tris.push_back(tB_id);
	}
}

template <typename Traits>
void AuxiliaryStructure<Traits>::addCoplanarTriangles(index_t ta, index_t tb)
{
	OMC_EXPENSIVE_ASSERT(ta != tb, "same triangles");
	OMC_EXPENSIVE_ASSERT(ta < coplanar_tris.size() && tb < coplanar_tris.size(),
	                     "out of range");

	if (coplanar_tris[ta].empty())
		coplanar_tris[ta].reserve(8);
	if (coplanar_tris[tb].empty())
		coplanar_tris[tb].reserve(8);
	coplanar_tris[ta].push_back(tb);
	coplanar_tris[tb].push_back(ta);
}

template <typename Traits>
void AuxiliaryStructure<Traits>::addCoplanarEdge(index_t t_id, index_t e_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_edges.size(), "out of range");

	if (coplanar_edges[t_id].empty())
		coplanar_edges[t_id].reserve(8);
	if (contains(coplanar_edges[t_id], e_id))
		return;
	coplanar_edges[t_id].push_back(e_id);
}

template <typename Traits>
const std::vector<size_t> &
AuxiliaryStructure<Traits>::coplanarTriangles(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_tris.size(), "out of range");
	return coplanar_tris[t_id];
}

template <typename Traits>
const std::vector<size_t> &
AuxiliaryStructure<Traits>::sortedCoplanarTriangles(index_t t_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_tris.size(), "out of range");
	if (!coplanar_tris_sorted[t_id])
	{
		std::sort(coplanar_tris[t_id].begin(), coplanar_tris[t_id].end());
		coplanar_tris_sorted[t_id] = true;
	}
	return coplanar_tris[t_id];
}

template <typename Traits>
const std::vector<size_t> &
AuxiliaryStructure<Traits>::coplanarEdges(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_edges.size(), "out of range");
	return coplanar_edges[t_id];
}

template <typename Traits>
bool AuxiliaryStructure<Traits>::triangleHasCoplanars(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_tris.size(), "out of range");
	return (!coplanar_tris[t_id].empty() || !coplanar_edges[t_id].empty());
}

template <typename Traits>
void AuxiliaryStructure<Traits>::setTriangleHasIntersections(index_t t_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri_has_intersections.size(), "out of range");
	tri_has_intersections[t_id] = true;
}

template <typename Traits>
bool AuxiliaryStructure<Traits>::triangleHasIntersections(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri_has_intersections.size(), "out of range");
	return tri_has_intersections[t_id];
}

template <typename Traits>
const std::vector<size_t> &
AuxiliaryStructure<Traits>::trianglePointsList(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2pts.size(), "out of range");
	return tri2pts[t_id];
}

template <typename Traits>
const std::vector<size_t> &
AuxiliaryStructure<Traits>::edgePointsList(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	return edge2pts[e_id];
}

template <typename Traits>
const std::vector<UIPair> &
AuxiliaryStructure<Traits>::triangleSegmentsList(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2segs.size(), "out of range");
	return tri2segs[t_id];
}

template <typename Traits>
const std::vector<size_t> &
AuxiliaryStructure<Traits>::segmentTrianglesList(const UIPair &seg) const
{
	auto res = seg2tris.find(seg);
	OMC_EXPENSIVE_ASSERT(res != seg2tris.end(), "out of range");
	return res->second;
}

template <typename Traits>
template <typename GetIndex>
std::pair<index_t, bool> AuxiliaryStructure<Traits>::addVertexInSortedList(
  const GPoint *pp, std::atomic<index_t> *ip, GetIndex get_idx)
{
	return static_cast<AuxPointMap_OcTree<Traits> *>(v_map.get())
	  ->insert_F(pp, ip, get_idx);
}

template <typename Traits>
index_t AuxiliaryStructure<Traits>::addVisitedPolygonPocket(
  const std::vector<index_t> &polygon, size_t pos)
{
	auto poly_it = pockets_map.insert(std::make_pair(polygon, pos));

	if (poly_it.second)
		return InvalidIndex; // polygon not present yet

	return poly_it.first->second;
}

template <typename Traits>
template <typename Comparator>
void AuxiliaryStructure<Traits>::sortEdgeList(index_t eid, Comparator comp)
{
	std::sort(edge2pts[eid].begin(), edge2pts[eid].end(), comp);
}

} // namespace OMC