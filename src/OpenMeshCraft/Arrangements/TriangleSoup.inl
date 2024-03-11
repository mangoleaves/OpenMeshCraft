#pragma once

#include "TriangleSoup.h"

namespace OMC {

template <typename Traits>
void TriangleSoup<Traits>::initialize()
{
	num_orig_vtxs = static_cast<size_t>(vertices.size());
	num_orig_tris = static_cast<size_t>(triangles.size() / 3);

	edges.reserve(numVerts() + numTris());
	edge_map.reserve(numVerts() + numTris());
	tri_planes.resize(numTris());
	tri_orientations.resize(numTris());

	// planes and edges
	auto tmp_edges            = std::vector<Edge>(num_orig_tris * 3);
	// this is done parallely since it is expensive
	auto init_tri_plane_edges = [this, &tmp_edges](index_t t_id)
	{
		index_t v0_id = triVertID(t_id, 0);
		index_t v1_id = triVertID(t_id, 1);
		index_t v2_id = triVertID(t_id, 2);

		tri_planes[t_id] = intToPlane(MaxCompInTriNormal()(
		  vertX(v0_id), vertY(v0_id), vertZ(v0_id), vertX(v1_id), vertY(v1_id),
		  vertZ(v1_id), vertX(v2_id), vertY(v2_id), vertZ(v2_id)));
		tri_orientations[t_id] =
		  OrientOn2D()(vert(v0_id), vert(v1_id), vert(v2_id), tri_planes[t_id]);

		tmp_edges[t_id * 3 + 0] = uniquePair(v0_id, v1_id);
		tmp_edges[t_id * 3 + 1] = uniquePair(v1_id, v2_id);
		tmp_edges[t_id * 3 + 2] = uniquePair(v2_id, v0_id);
	};

	tbb::parallel_for((size_t)0, num_orig_tris, init_tri_plane_edges);

	tbb::parallel_sort(tmp_edges.begin(), tmp_edges.end());

	for (index_t e_id = 0; e_id < tmp_edges.size(); e_id++)
	{
		if (e_id == 0 || tmp_edges[e_id] != tmp_edges[e_id - 1])
			addEdge(tmp_edges[e_id].first, tmp_edges[e_id].second);
	}
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
auto TriangleSoup<Traits>::vertX(index_t v_id) const -> NT
{
	OMC_EXPENSIVE_ASSERT(v_id < num_orig_vtxs,
	                     "vtx id out of range of original points");
	return vertices[v_id]->x();
}

template <typename Traits>
auto TriangleSoup<Traits>::vertY(index_t v_id) const -> NT
{
	OMC_EXPENSIVE_ASSERT(v_id < num_orig_vtxs,
	                     "vtx id out of range of original points");
	return vertices[v_id]->y();
}

template <typename Traits>
auto TriangleSoup<Traits>::vertZ(index_t v_id) const -> NT
{
	OMC_EXPENSIVE_ASSERT(v_id < num_orig_vtxs,
	                     "vtx id out of range of original points");
	return vertices[v_id]->z();
}

template <typename Traits>
index_t TriangleSoup<Traits>::addImplVert(GPoint *pp, std::atomic<index_t> *ip)
{
	vertices.push_back(pp);
	indices.push_back(ip);
	return vertices.size() - 1;
}

template <typename Traits>
auto TriangleSoup<Traits>::edge(index_t e_id) const -> const Edge &
{
	OMC_EXPENSIVE_ASSERT(e_id < edges.size(), "e_id out of range.");
	return edges[e_id];
}

template <typename Traits>
index_t TriangleSoup<Traits>::edgeID(index_t v0_id, index_t v1_id) const
{
	auto it = edge_map.find(uniquePair(v0_id, v1_id));

	if (it == edge_map.end())
		return InvalidIndex;
	return it->second; // edge id
}

template <typename Traits>
auto TriangleSoup<Traits>::edgeVert(index_t e_id, size_t off) const
  -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(e_id < edges.size(), "e_id out of range");
	if (off == 0)
		return vert(edges[e_id].first);
	else
		return vert(edges[e_id].second);
}

template <typename Traits>
auto TriangleSoup<Traits>::edgeVertPtr(index_t e_id, size_t off) const
  -> const NT *
{
	OMC_EXPENSIVE_ASSERT(e_id < edges.size(), "e_id out of range");
	if (off == 0)
		return vertPtr(edges[e_id].first);
	else
		return vertPtr(edges[e_id].second);
}

template <typename Traits>
index_t TriangleSoup<Traits>::edgeOppositeToVert(index_t t_id,
                                                 index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");

	index_t e_id = InvalidIndex;

	if (triVertID(t_id, 0) == v_id)
		e_id = edgeID(triVertID(t_id, 1), triVertID(t_id, 2));
	else if (triVertID(t_id, 1) == v_id)
		e_id = edgeID(triVertID(t_id, 0), triVertID(t_id, 2));
	else if (triVertID(t_id, 2) == v_id)
		e_id = edgeID(triVertID(t_id, 0), triVertID(t_id, 1));

	OMC_EXPENSIVE_ASSERT(e_id != InvalidIndex, "Opposite edge not found");
	return e_id;
}

template <typename Traits>
void TriangleSoup<Traits>::addEdge(index_t v0_id, index_t v1_id)
{
	index_t tmp_id = static_cast<index_t>(edges.size());
	Edge    e      = uniquePair(v0_id, v1_id);

	auto it = edge_map.insert({e, tmp_id});

	if (it.second)
		edges.push_back(e);
}

template <typename Traits>
bool TriangleSoup<Traits>::edgeContainsVert(index_t e_id, index_t v_id) const
{
	return edges[e_id].first == v_id || edges[e_id].second == v_id;
}

template <typename Traits>
const index_t *TriangleSoup<Traits>::tri(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return &triangles[3 * t_id];
}

template <typename Traits>
index_t TriangleSoup<Traits>::triVertID(index_t t_id, size_t off) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return triangles[3 * t_id + off];
}

template <typename Traits>
auto TriangleSoup<Traits>::triVert(index_t t_id, size_t off) const
  -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vert(triangles[3 * t_id + off]);
}

template <typename Traits>
auto TriangleSoup<Traits>::triVertPtr(index_t t_id, size_t off) const
  -> const NT *
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vertPtr(triangles[3 * t_id + off]);
}

template <typename Traits>
index_t TriangleSoup<Traits>::triEdgeID(index_t t_id, size_t off) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	index_t e_id =
	  edgeID(triangles[3 * t_id + off], triangles[3 * t_id + ((off + 1) % 3)]);

	OMC_EXPENSIVE_ASSERT(e_id != InvalidIndex, "no triangle edge found");
	return e_id;
}

template <typename Traits>
Plane TriangleSoup<Traits>::triPlane(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return tri_planes[t_id];
}

template <typename Traits>
Sign TriangleSoup<Traits>::triOrientation(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return tri_orientations[t_id];
}

template <typename Traits>
bool TriangleSoup<Traits>::triContainsVert(index_t t_id, index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "v_id out of range");

	return (triangles[3 * t_id] == v_id || triangles[3 * t_id + 1] == v_id ||
	        triangles[3 * t_id + 2] == v_id);
}

template <typename Traits>
bool TriangleSoup<Traits>::triContainsEdge(const index_t t_id, index_t ev0_id,
                                           index_t ev1_id) const
{
	return (triContainsVert(t_id, ev0_id) && triContainsVert(t_id, ev1_id));
}

template <typename Traits>
Label TriangleSoup<Traits>::triLabel(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return tri_labels[t_id];
}

} // namespace OMC