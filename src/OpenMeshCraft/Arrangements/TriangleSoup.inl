#pragma once

#include "TriangleSoup.h"

namespace OMC {

template <typename Traits>
void TriangleSoup<Traits>::initialize()
{
	num_orig_vtxs = static_cast<size_t>(vertices.size());
	num_orig_tris = static_cast<size_t>(triangles.size() / 3);

	edges.reserve(numOrigVerts());
	edge_map.resize(numOrigVerts());
	edge_mutexes = std::vector<tbb::spin_mutex>(numOrigVerts());

	tri_edges.resize(numOrigTris() * 3, InvalidIndex);

	tri_has_intersections.resize(numOrigTris(), false);
	coplanar_tris.resize(numOrigTris());
	coplanar_edges.resize(numOrigTris());

	tri2pts.resize(numOrigTris());
	edge2pts.reserve(numOrigVerts());
	tri2segs.resize(numOrigTris());

	seg2tris.resize(tbb::this_task_arena::max_concurrency());
	seg_mutexes =
	  std::vector<tbb::spin_mutex>(tbb::this_task_arena::max_concurrency());
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
const UIPair &TriangleSoup<Traits>::edge(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge index out of range.");
	return edges[e_id];
}

template <typename Traits>
index_t TriangleSoup<Traits>::getOrAddEdge(index_t v0_id, index_t v1_id)
{
	Edge edge = uniquePair(v0_id, v1_id);

	phmap::flat_hash_map<index_t, index_t> &em = edge_map[edge.first];

	// lock until finding and/or adding operation end.
	std::lock_guard<tbb::spin_mutex> lock_edge_map(edge_mutexes[edge.first]);

	auto find_iter = em.find(edge.second);
	if (find_iter != em.end()) // edge exists, just return it.
		return find_iter->second;
	else // edge does not exist, add it
	{
		// 1. add a new edge and related instances.
		typename decltype(edges)::iterator    edges_iter;
		typename decltype(edge2pts)::iterator e2p_iter;
		{
			// lock until all adding operations end.
			std::lock_guard<tbb::spin_mutex> lock_add_edge(new_edge_mutex);

			edges_iter = edges.push_back(edge);
			e2p_iter   = edge2pts.emplace_back();
			OMC_EXPENSIVE_ASSERT(e2p_iter - edge2pts.begin() ==
			                       edges_iter - edges.begin(),
			                     "size mismatch.");
		}
		// 2. get new edge's index
		index_t edge_id     = edges_iter - edges.begin();
		// 3. build map between new edge and its index
		auto    insert_iter = em.insert({edge.second, edge_id});
		OMC_EXPENSIVE_ASSERT(insert_iter.second, "fail to add edge.");
		return edge_id;
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
auto TriangleSoup<Traits>::triVert(index_t t_id,
                                   index_t off) const -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vert(triangles[3 * t_id + off]);
}

template <typename Traits>
auto TriangleSoup<Traits>::triVertPtr(index_t t_id,
                                      index_t off) const -> const NT *
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vertPtr(triangles[3 * t_id + off]);
}

template <typename Traits>
index_t TriangleSoup<Traits>::triEdgeID(index_t t_id, index_t off)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	if (is_valid_idx(tri_edges[3 * t_id + off]))
		// edge index is cached, return it.
		return tri_edges[3 * t_id + off];
	// otherwise check if the edge exists
	index_t v0_id = triangles[3 * t_id + off];
	index_t v1_id = triangles[3 * t_id + ((off + 1) % 3)];
	Edge    edge  = uniquePair(v0_id, v1_id);
	if (edge_map[edge.first].find(edge.second) == edge_map[edge.first].end())
		// the edge does not exist, return invalid index.
		return InvalidIndex;
	// otherwise get the index of edge and return it.
	tri_edges[3 * t_id + off] = edge_map[edge.first].at(edge.second);
	return tri_edges[3 * t_id + off];
}

template <typename Traits>
index_t TriangleSoup<Traits>::triEdgeID(index_t t_id, index_t off,
                                        std::true_type)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	// edge index is cached, return it.
	if (is_valid_idx(tri_edges[3 * t_id + off]))
		return tri_edges[3 * t_id + off];
	// otherwise get or add the edge
	index_t e_id              = getOrAddEdge(triangles[3 * t_id + off],
	                                         triangles[3 * t_id + ((off + 1) % 3)]);
	// get the index and return it
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
void TriangleSoup<Traits>::buildVMap(Tree *tree)
{
	tree->clear_points();
	v_map = std::make_unique<AuxPointMap_Tree<Traits>>(tree);
	tbb::parallel_for(index_t(0), numVerts(), [this](index_t v_id)
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
	index_t outer_map = (seg.first + seg.second) % seg2tris.size();
	{
		std::lock_guard<tbb::spin_mutex> lock(seg_mutexes[outer_map]);
		std::vector<index_t>            &tris = seg2tris[outer_map][seg];
		if (tris.empty())
			tris.reserve(8);
		tris.push_back(tA_id);
		tris.push_back(tB_id);
	}
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
void TriangleSoup<Traits>::setTriangleHasIntersections(index_t tA_id,
                                                       index_t tB_id)
{
	OMC_EXPENSIVE_ASSERT(tA_id < tri_has_intersections.size(), "out of range");
	OMC_EXPENSIVE_ASSERT(tB_id < tri_has_intersections.size(), "out of range");
	tri_has_intersections[tA_id] = true;
	tri_has_intersections[tB_id] = true;
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
const std::vector<index_t> &
TriangleSoup<Traits>::segmentTrianglesList(const UIPair &seg) const
{
	index_t outer_map = (seg.first + seg.second) % seg2tris.size();
	auto    res       = seg2tris[outer_map].find(seg);
	OMC_EXPENSIVE_ASSERT(res != seg2tris.end(), "out of range");
	return res->second;
}

template <typename Traits>
Plane TriangleSoup<Traits>::triPlane(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numOrigTris(), "out of range.");
	OMC_EXPENSIVE_ASSERT(planeToInt(tri_plane[t_id]) >= 0, "invalid plane");
	return tri_plane[t_id];
}

template <typename Traits>
Sign TriangleSoup<Traits>::triOrient(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numOrigTris(), "out of range.");
	OMC_EXPENSIVE_ASSERT(is_sign_reliable(tri_orient[t_id]), "invalid sign");
	return tri_orient[t_id];
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
void TriangleSoup<Traits>::removeAllDuplicates()
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

	tbb::parallel_for_each(
	  seg2tris.begin(), seg2tris.end(),
	  [this](phmap::flat_hash_map<UIPair, std::vector<index_t>> &s2t)
	  {
		  for (auto &iter : s2t)
			  remove_duplicates(iter.second);
	  });
}

template <typename Traits>
void TriangleSoup<Traits>::calcPlaneAndOrient()
{
	tri_plane.resize(numOrigTris(), intToPlane(-1));
	tri_orient.resize(numOrigTris(), Sign::UNCERTAIN);

	auto calc_plane_and_orient = [this](index_t t_id)
	{
		if (triangleHasIntersections(t_id))
		{
			index_t v0_id = triVertID(t_id, 0);
			index_t v1_id = triVertID(t_id, 1);
			index_t v2_id = triVertID(t_id, 2);

			tri_plane[t_id] = intToPlane(
			  MaxCompInTriNormal()(vertPtr(v0_id), vertPtr(v1_id), vertPtr(v2_id)));
			tri_orient[t_id] = OrientOn2D()(vertPtr(v0_id), vertPtr(v1_id),
			                                vertPtr(v2_id), tri_plane[t_id]);
		}
	};

	tbb::parallel_for(size_t(0), numOrigTris(), calc_plane_and_orient);
}

template <typename Traits>
template <typename Comparator>
void TriangleSoup<Traits>::sortEdgeList(index_t eid, Comparator comp)
{
	std::sort(edge2pts[eid].begin(), edge2pts[eid].end(), comp);
}

} // namespace OMC