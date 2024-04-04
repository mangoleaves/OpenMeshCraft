#pragma once

#include "TriangleSoup.h"

namespace OMC {

template <typename Traits>
void TriangleSoup<Traits>::initialize()
{
	num_orig_vtxs = static_cast<size_t>(vertices.size());
	num_orig_tris = static_cast<size_t>(triangles.size() / 3);

	largest_edge_id.store(0);
	edge_map.resize(numOrigVerts());
	edge_mutexes.resize(numOrigVerts());

	tri_edges.resize(numOrigTris() * 3, InvalidIndex);

	coplanar_tris.resize(numOrigTris());
	coplanar_edges.resize(numOrigTris());
	tri_has_intersections.resize(numOrigTris(), false);
	coplanar_tris_sorted.resize(numOrigTris(), false);

	tri2pts.resize(numOrigTris());
	edge2pts.reserve(numOrigVerts());
	tri2segs.resize(numOrigTris());
}

template <typename Traits>
auto TriangleSoup<Traits>::vert(index_t v_id) const -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	return *vertices[v_id];
}

template <typename Traits>
auto TriangleSoup<Traits>::vertPtr(index_t v_id) const -> const NT *
{
	OMC_EXPENSIVE_ASSERT(v_id < num_orig_vtxs,
	                     "vtx id out of range of original points");
	return vertices[v_id]->data();
}

template <typename Traits>
index_t TriangleSoup<Traits>::addImplVert(GPoint *pp, std::atomic<index_t> *ip)
{
	vertices.push_back(pp);
	indices.push_back(ip);
	return vertices.size() - 1;
}

template <typename Traits>
index_t TriangleSoup<Traits>::getOrAddEdge(index_t v0_id, index_t v1_id)
{
	Edge edge = uniquePair(v0_id, v1_id);

	phmap::flat_hash_map<index_t, index_t> &em = edge_map[edge.first];

	{ // find the edge in edge_map or add the edge to edge_map
		std::lock_guard<tbb::spin_mutex> lock(edge_mutexes[edge.first]);

		auto find_iter = em.find(edge.second);

		if (find_iter == em.end())
		{ // edge does not exist, add it
			index_t edge_id     = largest_edge_id.fetch_add(1);
			auto    insert_iter = em.insert({edge.second, edge_id});
			OMC_EXPENSIVE_ASSERT(insert_iter.second, "fail to add edge.");
			return edge_id;
		}
		else // edge exists, return it
			return find_iter->second;
	}
}

template <typename Traits>
const index_t *TriangleSoup<Traits>::tri(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return &triangles[3 * t_id];
}

template <typename Traits>
index_t TriangleSoup<Traits>::triVertID(index_t t_id, index_t off) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return triangles[3 * t_id + off];
}

template <typename Traits>
auto TriangleSoup<Traits>::triVert(index_t t_id, index_t off) const
  -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vert(triangles[3 * t_id + off]);
}

template <typename Traits>
auto TriangleSoup<Traits>::triVertPtr(index_t t_id, index_t off) const
  -> const NT *
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vertPtr(triangles[3 * t_id + off]);
}

template <typename Traits>
index_t TriangleSoup<Traits>::triEdgeID(index_t t_id, index_t off) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	if (is_valid_idx(tri_edges[3 * t_id + off]))
		return tri_edges[3 * t_id + off];

	index_t e_id = getOrAddEdge(triangles[3 * t_id + off],
	                            triangles[3 * t_id + ((off + 1) % 3)]);

	tri_edges[3 * t_id + off] = e_id;

	return e_id;
}

template <typename Traits>
Label TriangleSoup<Traits>::triLabel(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return tri_labels[t_id];
}

template <typename Traits>
void TriangleSoup<Traits>::build_vmap(Tree *tree)
{
	tree->clear_points();
	v_map = std::make_unique<AuxPointMap_Tree<Traits>>(tree);
	tbb::parallel_for(index_t(0), numVerts(),
	                  [this](index_t v_id)
	                  { v_map->insert(vertices[v_id], indices[v_id]); });
}

template <typename Traits>
void TriangleSoup<Traits>::addVertexInTriangle(index_t t_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2pts.size(), "out of range");
	tbb::concurrent_vector<index_t> &points = tri2pts[t_id];
	if (points.empty())
		points.reserve(8);
	points.push_back(v_id);
}

template <typename Traits>
void TriangleSoup<Traits>::addVertexInEdge(index_t e_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	tbb::concurrent_vector<index_t> &points = edge2pts[e_id];
	if (points.empty())
		points.reserve(8);
	points.push_back(v_id);
}

template <typename Traits>
void TriangleSoup<Traits>::addSegmentInTriangle(index_t t_id, const UIPair &seg)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2segs.size(), "out of range");
	OMC_EXPENSIVE_ASSERT(isUnique(seg), "segment is not unique");
	tbb::concurrent_vector<UIPair> &segments = tri2segs[t_id];
	if (segments.empty())
		segments.reserve(8);
	segments.push_back(seg);
}

template <typename Traits>
void TriangleSoup<Traits>::addTrianglesInSegment(const UIPair &seg,
                                                 index_t tA_id, index_t tB_id)
{
	OMC_EXPENSIVE_ASSERT(isUnique(seg), "segment is not unique");
	tbb::concurrent_vector<index_t> &tris = seg2tris[seg];
	if (tris.empty())
		tris.reserve(4);
	tris.push_back(tA_id);
	tris.push_back(tB_id);
}

template <typename Traits>
void TriangleSoup<Traits>::addCoplanarTriangles(index_t ta, index_t tb)
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
void TriangleSoup<Traits>::addCoplanarEdge(index_t t_id, index_t e_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_edges.size(), "out of range");

	if (coplanar_edges[t_id].empty())
		coplanar_edges[t_id].reserve(8);
	coplanar_edges[t_id].push_back(e_id);
}

template <typename Traits>
const tbb::concurrent_vector<index_t> &
TriangleSoup<Traits>::coplanarTriangles(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_tris.size(), "out of range");
	return coplanar_tris[t_id];
}

template <typename Traits>
const tbb::concurrent_vector<index_t> &
TriangleSoup<Traits>::coplanarEdges(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_edges.size(), "out of range");
	return coplanar_edges[t_id];
}

template <typename Traits>
void TriangleSoup<Traits>::setTriangleHasIntersections(index_t t_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri_has_intersections.size(), "out of range");
	tri_has_intersections[t_id] = true;
}

template <typename Traits>
bool TriangleSoup<Traits>::triangleHasIntersections(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri_has_intersections.size(), "out of range");
	return tri_has_intersections[t_id];
}

template <typename Traits>
const tbb::concurrent_vector<index_t> &
TriangleSoup<Traits>::trianglePointsList(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2pts.size(), "out of range");
	return tri2pts[t_id];
}

template <typename Traits>
const tbb::concurrent_vector<index_t> &
TriangleSoup<Traits>::edgePointsList(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	return edge2pts[e_id];
}

template <typename Traits>
const tbb::concurrent_vector<UIPair> &
TriangleSoup<Traits>::triangleSegmentsList(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2segs.size(), "out of range");
	return tri2segs[t_id];
}

template <typename Traits>
const tbb::concurrent_vector<index_t> &
TriangleSoup<Traits>::segmentTrianglesList(const UIPair &seg) const
{
	auto res = seg2tris.find(seg);
	OMC_EXPENSIVE_ASSERT(res != seg2tris.end(), "out of range");
	return res->second;
}

template <typename Traits>
template <typename GetIndex>
std::pair<index_t, bool> TriangleSoup<Traits>::addVertexInSortedList(
  const GPoint *pp, std::atomic<index_t> *ip, GetIndex get_idx)
{
	return static_cast<AuxPointMap_Tree<Traits> *>(v_map.get())
	  ->insert_F(pp, ip, get_idx);
}

template <typename Traits>
index_t TriangleSoup<Traits>::addVisitedPolygonPocket(
  const std::vector<index_t> &polygon, index_t pos)
{
	auto poly_it = pockets_map.insert(std::make_pair(polygon, pos));

	if (poly_it.second)
		return InvalidIndex; // polygon not present yet

	return poly_it.first->second;
}

template <typename Traits>
void TriangleSoup<Traits>::remove_all_duplicates()
{
	tbb::parallel_for(size_t(0), numOrigTris(),
	                  [this](size_t t_id)
	                  {
		                  remove_duplicates(coplanar_tris[t_id]);
		                  remove_duplicates(coplanar_edges[t_id]);
		                  remove_duplicates(tri2pts[t_id]);
		                  remove_duplicates(tri2segs[t_id]);
	                  });

	tbb::parallel_for_each(edge2pts.begin(), edge2pts.end(),
	                       [this](tbb::concurrent_vector<index_t> &e2p)
	                       { remove_duplicates(e2p); });

	tbb::parallel_for_each(seg2tris.begin(), seg2tris.end(),
	                       [this](auto &s2t) { remove_duplicates(s2t.second); });
}

template <typename Traits>
template <typename Comparator>
void TriangleSoup<Traits>::sortEdgeList(index_t eid, Comparator comp)
{
	std::sort(edge2pts[eid].begin(), edge2pts[eid].end(), comp);
}

} // namespace OMC