#pragma once

#include "FastTriMesh.h"
#include "Tree.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"
#include "OpenMeshCraft/Utils/Macros.h"

namespace OMC {

template <typename Traits>
FastTriMesh<Traits>::FastTriMesh()
{
}

template <typename Traits>
FastTriMesh<Traits>::FastTriMesh(const std::vector<GPoint *> &in_verts,
                                 const std::vector<index_t>  &in_tris)
{
	OMC_ASSERT(in_tris.size() % 3 == 0, "triangle error");

	preAllocateSpace(in_verts.size());
	// add vertices
	for (index_t v_id = 0; v_id < in_verts.size(); ++v_id)
	{
		vertices.emplace_back(in_verts[v_id], /*info*/ 0);
	}
	v2t.resize(vertices.size());
	v2e.resize(vertices.size());
	if constexpr (std::is_same_v<AuxVector16<index_t>, std::vector<index_t>>)
	{
		reserve(v2t, 16);
		reserve(v2e, 16);
	}

	// add triangles and (possibly duplicate) edges,
	// i.e., build t2v, e2v, v2t
	triangles.reserve(in_tris.size() / 3);
	edges.reserve(in_tris.size());
	for (index_t t_id = 0; t_id < in_tris.size() / 3; ++t_id)
	{
		index_t tv0_id = in_tris[3 * t_id], tv1_id = in_tris[3 * t_id + 1],
		        tv2_id = in_tris[3 * t_id + 2];
		triangles.emplace_back(tv0_id, tv1_id, tv2_id, /*info*/ 0);
		v2t[tv0_id].push_back(t_id);
		v2t[tv1_id].push_back(t_id);
		v2t[tv2_id].push_back(t_id);

		edges.emplace_back(tv0_id, tv1_id, /*constr*/ false);
		edges.emplace_back(tv1_id, tv2_id, /*constr*/ false);
		edges.emplace_back(tv2_id, tv0_id, /*constr*/ false);
	}

	// build unique edges
	tbb::parallel_sort(edges.begin(), edges.end());
	edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

	// build v2e
	for (index_t e_id = 0; e_id < edges.size(); ++e_id)
	{
		v2e[edges[e_id].verts.first].push_back(e_id);
		v2e[edges[e_id].verts.second].push_back(e_id);
	}

	// build f2e (iTri.edges)
	// clang-format off
	tbb::parallel_for(
	  (size_t)0, triangles.size(),
	  [this](size_t t_id)
	  {
		  auto [tv0_id, tv1_id, tv2_id] = triangles[t_id].verts;
		  index_t e0_id = (std::lower_bound(edges.begin(), edges.end(), iEdge(tv0_id, tv1_id)) - edges.begin());
		  index_t e1_id = (std::lower_bound(edges.begin(), edges.end(), iEdge(tv1_id, tv2_id)) - edges.begin());
		  index_t e2_id = (std::lower_bound(edges.begin(), edges.end(), iEdge(tv2_id, tv0_id)) - edges.begin());
		  triangles[t_id].edges = {e0_id, e1_id, e2_id};
	  });
	// clang-format on

	// build e2t
	e2t.resize(edges.size());
	if constexpr (std::is_same_v<AuxVector4<index_t>, std::vector<index_t>>)
	{
		reserve(e2t, 4);
	}
	for (index_t t_id = 0; t_id < triangles.size(); t_id++)
	{
		e2t[triangles[t_id].edges[0]].push_back(t_id);
		e2t[triangles[t_id].edges[1]].push_back(t_id);
		e2t[triangles[t_id].edges[2]].push_back(t_id);
	}
}

template <typename Traits>
auto FastTriMesh<Traits>::operator=(Self &rhs) noexcept -> Self &
{
	vertices       = rhs.vertices;
	edges          = rhs.edges;
	triangles      = rhs.triangles;
	v2e            = rhs.v2e;
	v2t            = rhs.v2t;
	e2t            = rhs.e2t;
	rev_vtx_map    = rhs.rev_vtx_map;
	triangle_plane = rhs.triangle_plane;
	orientation    = rhs.orientation;

	return *this;
}

template <typename Traits>
auto FastTriMesh<Traits>::operator=(Self &&rhs) noexcept -> Self &
{
	vertices       = std::move(rhs.vertices);
	edges          = std::move(rhs.edges);
	triangles      = std::move(rhs.triangles);
	v2e            = std::move(rhs.v2e);
	v2t            = std::move(rhs.v2t);
	e2t            = std::move(rhs.e2t);
	rev_vtx_map    = std::move(rhs.rev_vtx_map);
	triangle_plane = rhs.triangle_plane;
	orientation    = rhs.orientation;

	return *this;
}

template <typename Traits>
FastTriMesh<Traits>::~FastTriMesh()
{
}

template <typename Traits>
void FastTriMesh<Traits>::initialize(const GPoint *tv0, const GPoint *tv1,
                                     const GPoint *tv2, const index_t *tv_id,
                                     const Plane &ref_p, const Sign &ori)
{
	clear();

	addVert(tv0, tv_id[0]);
	addVert(tv1, tv_id[1]);
	addVert(tv2, tv_id[2]);
	addTri(0, 1, 2);

	triangle_plane = ref_p;
	orientation    = ori;
}

template <typename Traits>
void FastTriMesh<Traits>::preAllocateSpace(size_t estimated_num_verts)
{
	vertices.reserve(estimated_num_verts);
	rev_vtx_map.reserve(vertices.size());
	edges.reserve(vertices.size() / 2);
	triangles.reserve(vertices.size() / 3);
	v2e.reserve(vertices.size());
	v2t.reserve(vertices.size());
	e2t.reserve(edges.size());
}

template <typename Traits>
void FastTriMesh<Traits>::clear()
{
	vertices.clear();
	rev_vtx_map.clear();
	edges.clear();
	triangles.clear();
	v2e.clear();
	v2t.clear();
	e2t.clear();
}

template <typename Traits>
size_t FastTriMesh<Traits>::numVerts() const
{
	return vertices.size();
}

template <typename Traits>
size_t FastTriMesh<Traits>::numEdges() const
{
	return edges.size();
}

template <typename Traits>
size_t FastTriMesh<Traits>::numTriangles() const
{
	return triangles.size();
}

template <typename Traits>
void FastTriMesh<Traits>::resetVerticesInfo()
{
	for (iVertex &v : vertices)
		v.info = 0;
}

template <typename Traits>
void FastTriMesh<Traits>::resetTrianglesInfo()
{
	for (iTri &t : triangles)
		t.info = 0;
}

template <typename Traits>
index_t FastTriMesh<Traits>::vertInfo(const index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	return vertices[v_id].info;
}

template <typename Traits>
void FastTriMesh<Traits>::setVertInfo(const index_t v_id, const index_t info)
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	vertices[v_id].info = info;
}

template <typename Traits>
bool FastTriMesh<Traits>::vertFlag(const index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	return vertices[v_id].flag;
}

template <typename Traits>
void FastTriMesh<Traits>::setVertFlag(const index_t v_id, const bool flag)
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	vertices[v_id].flag = flag;
}

template <typename Traits>
index_t FastTriMesh<Traits>::triInfo(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");

	return triangles[t_id].info;
}

template <typename Traits>
void FastTriMesh<Traits>::setTriInfo(index_t t_id, index_t val)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	triangles[t_id].info = val;
}

template <typename Traits>
index_t FastTriMesh<Traits>::meshInfo() const
{
	return mesh_info;
}

template <typename Traits>
void FastTriMesh<Traits>::setMeshInfo(index_t val)
{
	mesh_info = val;
}

template <typename Traits>
const Label &FastTriMesh<Traits>::triLabel(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");

	return triangles[t_id].label;
}

template <typename Traits>
void FastTriMesh<Traits>::setTriLabel(index_t t_id, const Label &label)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	triangles[t_id].label = label;
}

template <typename Traits>
index_t FastTriMesh<Traits>::vertLocalID(index_t glob_v_id) const
{
	auto it = rev_vtx_map.find(glob_v_id);
	OMC_EXPENSIVE_ASSERT(it != rev_vtx_map.end(),
	                     "vtx id not found in reverse map");
	return it->second;
}

template <typename Traits>
bool FastTriMesh<Traits>::edgeContainsVert(index_t e_id, index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	const auto &ev = edges[e_id].verts;
	return ev.first == v_id || ev.second == v_id;
}

template <typename Traits>
bool FastTriMesh<Traits>::triContainsVert(index_t t_id, index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "tri id out of range");
	const auto &tv = triangles[t_id].verts;
	return tv[0] == v_id || tv[1] == v_id || tv[2] == v_id;
}

template <typename Traits>
const AuxVector4<index_t> &FastTriMesh<Traits>::adjE2T(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	return e2t[e_id];
}

template <typename Traits>
const AuxVector16<index_t> &FastTriMesh<Traits>::adjV2T(index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	return v2t[v_id];
}

template <typename Traits>
const AuxVector16<index_t> &FastTriMesh<Traits>::adjV2E(index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	return v2e[v_id];
}

template <typename Traits>
const std::array<index_t, 3> &FastTriMesh<Traits>::adjT2E(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	return triangles[t_id].edges;
}

template <typename Traits>
index_t FastTriMesh<Traits>::triVertOppositeTo(index_t t_id, index_t v0_id,
                                               index_t v1_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(v0_id != v1_id, "verts are equal");
	OMC_EXPENSIVE_ASSERT(
	  (triContainsVert(t_id, v0_id) && triContainsVert(t_id, v1_id)),
	  "tri doesn't contain vtx");
	// possibly overflow, but finally right
	return triangles[t_id].verts[0] + triangles[t_id].verts[1] +
	       triangles[t_id].verts[2] - v0_id - v1_id;
}

template <typename Traits>
index_t FastTriMesh<Traits>::edgeOppToVert(index_t t_id, index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(triContainsVert(t_id, v_id), "tri doesn't contain vtx");
	index_t     e_id = InvalidIndex;
	const auto &tv   = triangles[t_id].verts;

	if (tv[0] == v_id)
		e_id = triEdgeID(t_id, 1);
	else if (tv[1] == v_id)
		e_id = triEdgeID(t_id, 2);
	else if (tv[2] == v_id)
		e_id = triEdgeID(t_id, 0);

	OMC_EXPENSIVE_ASSERT(is_valid_idx(e_id), "opposite edge not found in tri");
	return e_id;
}

template <typename Traits>
index_t FastTriMesh<Traits>::triOppToEdge(index_t e_id, index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(edgeIsManifold(e_id),
	                     "edge is non-manifold, opp tri is undefined.");
	if (e2t[e_id].size() == 1)
		return InvalidIndex; // boundary edge
	// possibly overflow, but finally right
	return e2t[e_id][0] + e2t[e_id][1] - t_id;
}

template <typename Traits>
index_t FastTriMesh<Traits>::nextVertInTri(index_t t_id,
                                           index_t curr_v_id) const
{
	const auto &tv = triangles[t_id].verts;

	if (tv[0] == curr_v_id)
		return tv[1];
	if (tv[1] == curr_v_id)
		return tv[2];
	if (tv[2] == curr_v_id)
		return tv[0];

	OMC_ASSERT(false, "This should not happen");
	return InvalidIndex; // warning killer
}

template <typename Traits>
index_t FastTriMesh<Traits>::prevVertInTri(index_t t_id,
                                           index_t curr_v_id) const
{
	const auto &tv = triangles[t_id].verts;

	if (tv[0] == curr_v_id)
		return tv[2];
	if (tv[1] == curr_v_id)
		return tv[0];
	if (tv[2] == curr_v_id)
		return tv[1];

	OMC_ASSERT(false, "This should not happen");
	return InvalidIndex; // warning killer
}

template <typename Traits>
index_t FastTriMesh<Traits>::addVert(const GPoint *p, index_t origin_v_id)
{
	index_t v_id = static_cast<index_t>(numVerts());
	vertices.emplace_back(p, origin_v_id);

	v2e.emplace_back();
	v2t.emplace_back();
	if constexpr (std::is_same_v<AuxVector16<index_t>, std::vector<index_t>>)
	{
		v2e.back().reserve(8);
		v2t.back().reserve(8);
	}
	rev_vtx_map[origin_v_id] = v_id;

	return v_id;
}

template <typename Traits>
index_t FastTriMesh<Traits>::addEdge(index_t ev0_id, index_t ev1_id)
{
	OMC_EXPENSIVE_ASSERT((ev0_id < numVerts() && ev1_id < numVerts()),
	                     "vtx id out of range.");
	OMC_EXPENSIVE_ASSERT((ev0_id != ev1_id), "degenerate triangle");

	index_t e_id = edgeID(ev0_id, ev1_id);
	if (is_valid_idx(e_id))
		return e_id;

	e_id = edges.size();
	edges.emplace_back(ev0_id, ev1_id, /*constr*/ false);

	e2t.emplace_back();
	if constexpr (std::is_same_v<AuxVector4<index_t>, std::vector<index_t>>)
	{
		e2t.back().reserve(4);
	}

	v2e[ev0_id].push_back(e_id);
	v2e[ev1_id].push_back(e_id);
	return e_id;
}

template <typename Traits>
index_t FastTriMesh<Traits>::addTri(index_t tv0_id, index_t tv1_id,
                                    index_t tv2_id)
{
	OMC_EXPENSIVE_ASSERT(
	  (tv0_id < numVerts() && tv1_id < numVerts(), tv2_id < numVerts()),
	  "vtx id out of range");
	OMC_EXPENSIVE_ASSERT(
	  (tv0_id != tv1_id && tv0_id != tv2_id && tv1_id != tv2_id),
	  "degenerate triangle");

	index_t t_id = triID(tv0_id, tv1_id, tv2_id);
	if (is_valid_idx(t_id))
		return t_id;

	t_id = numTriangles();

	triangles.emplace_back(tv0_id, tv1_id, tv2_id);

	index_t e0_id          = addEdge(tv0_id, tv1_id);
	index_t e1_id          = addEdge(tv1_id, tv2_id);
	index_t e2_id          = addEdge(tv2_id, tv0_id);
	triangles.back().edges = {e0_id, e1_id, e2_id};

	e2t[e0_id].push_back(t_id);
	e2t[e1_id].push_back(t_id);
	e2t[e2_id].push_back(t_id);
	v2t[tv0_id].push_back(t_id);
	v2t[tv1_id].push_back(t_id);
	v2t[tv2_id].push_back(t_id);
	return t_id;
}

template <typename Traits>
void FastTriMesh<Traits>::removeTri(index_t t_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	// resolve reference in v2t...
	index_t v0_id = triVertID(t_id, 0);
	index_t v1_id = triVertID(t_id, 1);
	index_t v2_id = triVertID(t_id, 2);
	removeFromVec(v2t[v0_id], t_id);
	removeFromVec(v2t[v1_id], t_id);
	removeFromVec(v2t[v2_id], t_id);
	// ...and e2t.
	index_t e0_id = triEdgeID(t_id, 0);
	index_t e1_id = triEdgeID(t_id, 1);
	index_t e2_id = triEdgeID(t_id, 2);
	removeFromVec(e2t[e0_id], t_id);
	removeFromVec(e2t[e1_id], t_id);
	removeFromVec(e2t[e2_id], t_id);
	// find dangling edges and...
	size_t                 n_dangling_edges = 0;
	std::array<index_t, 3> dangling_edges; // higher ids first
	if (e2t[e0_id].empty())
		dangling_edges[n_dangling_edges++] = e0_id;
	if (e2t[e1_id].empty())
		dangling_edges[n_dangling_edges++] = e1_id;
	if (e2t[e2_id].empty())
		dangling_edges[n_dangling_edges++] = e2_id;
	std::sort(dangling_edges.begin(), dangling_edges.begin() + n_dangling_edges,
	          std::greater<index_t>());
	// ...remove reference in v2e and unreferenced dangling edges.
	for (index_t i = 0; i < n_dangling_edges; i++)
	{
		index_t e_id_  = dangling_edges[i];
		index_t v0_id_ = edges[e_id_].verts.first;
		index_t v1_id_ = edges[e_id_].verts.second;
		removeFromVec(v2e[v0_id_], e_id_);
		removeFromVec(v2e[v1_id_], e_id_);
		removeEdgeUnref(e_id_);
	}
	// finnaly remove unreferenced triangle
	removeTriUnref(t_id);
}

template <typename Traits>
template <typename Container>
void FastTriMesh<Traits>::removeTris(const Container &t_ids)
{
	for (index_t t_id : t_ids)
		removeTri(t_id);
}

template <typename Traits>
void FastTriMesh<Traits>::removeEdge(index_t e_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	// remove adjacent triangles and finally remove the dangling edge e_id.
	removeTris(e2t[e_id]);
}

template <typename Traits>
void FastTriMesh<Traits>::splitEdge(const index_t e_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");

	index_t ev0_id = edges[e_id].verts.first;
	index_t ev1_id = edges[e_id].verts.second;

	// for each adjacent triangle
	AuxVector4<index_t> tmp_e2t = e2t[e_id];
	for (index_t t_id : tmp_e2t)
	{
		index_t v_opp = triVertOppositeTo(t_id, ev0_id, ev1_id);
		if (nextVertInTri(t_id, ev0_id) != ev1_id)
			std::swap(ev0_id, ev1_id);
		// split it to two new triangles
		addTri(v_opp, ev0_id, v_id);
		addTri(v_opp, v_id, ev1_id);
	}
	// then remove old triangles and edge
	removeTris(tmp_e2t);
}

template <typename Traits>
void FastTriMesh<Traits>::splitEdge(const index_t e_id, index_t v_id,
                                    SplitTree &tree)
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");

	index_t             ev0_id  = edges[e_id].verts.first;
	index_t             ev1_id  = edges[e_id].verts.second;
	// for each adjacent triangle
	AuxVector4<index_t> tmp_e2t = e2t[e_id];
	for (index_t t_id : tmp_e2t)
	{
		index_t v_opp = triVertOppositeTo(t_id, ev0_id, ev1_id);
		if (nextVertInTri(t_id, ev0_id) != ev1_id)
		{ // if true, ev1->ev0, swap them to ev0->ev1.
			std::swap(ev0_id, ev1_id);
		}
		// split it to two new triangles
		index_t t0_id = addTri(v_opp, ev0_id, v_id);
		index_t t1_id = addTri(v_opp, v_id, ev1_id);
		// add node to tree
		size_t  n0_id = tree.addNode(v_opp, ev0_id, v_id);
		size_t  n1_id = tree.addNode(v_opp, v_id, ev1_id);

		size_t node_id = triInfo(t_id);
		// set split point
		tree.setSplitPoint(node_id, v_id);
		// offset vertices
		tree.offsetVertesInNode(node_id, 3 - triVertOffset(t_id, ev0_id));
		// build parent-children
		tree.addChildren(node_id, n0_id, n1_id);
		// save node info
		setTriInfo(t0_id, n0_id);
		setTriInfo(t1_id, n1_id);
	}
	// remove old triangles and edge
	removeTris(tmp_e2t);
}

template <typename Traits>
void FastTriMesh<Traits>::splitTri(index_t t_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	// split the triangle to three new triangles
	addTri(triVertID(t_id, 0), triVertID(t_id, 1), v_id);
	addTri(triVertID(t_id, 1), triVertID(t_id, 2), v_id);
	addTri(triVertID(t_id, 2), triVertID(t_id, 0), v_id);
	// remove the old triangle
	removeTri(t_id);
}

template <typename Traits>
void FastTriMesh<Traits>::splitTri(index_t t_id, index_t v_id, SplitTree &tree)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");

	index_t node_id = triInfo(t_id);
	index_t v0 = triVertID(t_id, 0), v1 = triVertID(t_id, 1),
	        v2 = triVertID(t_id, 2);

	// split to three new triangles
	index_t t0_id = addTri(v0, v1, v_id);
	index_t t1_id = addTri(v1, v2, v_id);
	index_t t2_id = addTri(v2, v0, v_id);
	// add three nodes to tree
	index_t n0_id = tree.addNode(v0, v1, v_id);
	index_t n1_id = tree.addNode(v1, v2, v_id);
	index_t n2_id = tree.addNode(v2, v0, v_id);
	// set split point
	tree.setSplitPoint(node_id, v_id);
	// build parent-children
	tree.addChildren(node_id, n0_id, n1_id, n2_id);
	// the triangle label contains the node position
	setTriInfo(t0_id, n0_id);
	setTriInfo(t1_id, n1_id);
	setTriInfo(t2_id, n2_id);
	// remove the old triangle
	removeTri(t_id);
}

template <typename Traits>
void FastTriMesh<Traits>::flipTri(index_t t_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	std::swap(triangles[t_id].verts[0], triangles[t_id].verts[2]);
	std::swap(triangles[t_id].edges[0], triangles[t_id].edges[1]);
}

template <typename Traits>
auto FastTriMesh<Traits>::vert(index_t v_id) const -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range.");
	return *vertices[v_id].point;
}

template <typename Traits>
auto FastTriMesh<Traits>::triVert(index_t t_id, size_t off) const
  -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	return *vertices[triangles[t_id].verts[off]].point;
}

template <typename Traits>
index_t FastTriMesh<Traits>::triVertOffset(index_t t_id, index_t v_id) const
{
	const auto &tv = triangles[t_id].verts;

	if (tv[0] == v_id)
		return 0;
	if (tv[1] == v_id)
		return 1;
	if (tv[2] == v_id)
		return 2;

	OMC_ASSERT(false, "This should not happen");
	return InvalidIndex; // warning killer
}

template <typename Traits>
const index_t *FastTriMesh<Traits>::tri(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	return triangles[t_id].verts.data();
}

template <typename Traits>
index_t FastTriMesh<Traits>::triID(index_t tv0_id, index_t tv1_id,
                                   index_t tv2_id) const
{
	OMC_EXPENSIVE_ASSERT(
	  (tv0_id < numVerts() && tv1_id < numVerts() && tv2_id < numVerts()),
	  "vtx id out of range");
	index_t v_id =
	  v2t[tv0_id].size() <= v2t[tv1_id].size()
	    ? (v2t[tv0_id].size() <= v2t[tv2_id].size() ? tv0_id : tv2_id)
	    : (v2t[tv1_id].size() <= v2t[tv2_id].size() ? tv1_id : tv2_id);

	if (v_id == tv0_id)
	{
		for (index_t t_id : v2t[v_id])
			if (triContainsVert(t_id, tv1_id) && triContainsVert(t_id, tv2_id))
				return t_id;
	}
	else if (v_id == tv1_id)
	{
		for (index_t t_id : v2t[v_id])
			if (triContainsVert(t_id, tv0_id) && triContainsVert(t_id, tv2_id))
				return t_id;
	}
	else
	{
		for (index_t t_id : v2t[v_id])
			if (triContainsVert(t_id, tv0_id) && triContainsVert(t_id, tv1_id))
				return t_id;
	}
	return InvalidIndex;
}

template <typename Traits>
index_t FastTriMesh<Traits>::edgeID(index_t ev0_id, index_t ev1_id) const
{
	OMC_EXPENSIVE_ASSERT(ev0_id != ev1_id, "edge with equal endpoints");
	OMC_EXPENSIVE_ASSERT((ev0_id < numVerts() && ev1_id < numVerts()),
	                     "vtx id out of range");
	if (v2e[ev0_id].size() <= v2e[ev1_id].size())
	{
		for (index_t e_id : v2e[ev0_id])
			if (edgeContainsVert(e_id, ev1_id))
				return e_id;
	}
	else
	{
		for (index_t e_id : v2e[ev1_id])
			if (edgeContainsVert(e_id, ev0_id))
				return e_id;
	}
	return InvalidIndex;
}

template <typename Traits>
index_t FastTriMesh<Traits>::triVertID(index_t t_id, size_t off) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(off < 3, "offset out of range.");
	return triangles[t_id].verts[off];
}

template <typename Traits>
index_t FastTriMesh<Traits>::triEdgeID(index_t t_id, index_t off) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	OMC_EXPENSIVE_ASSERT(off < 3, "offset out of range.");
	return triangles[t_id].edges[off];
}

template <typename Traits>
index_t FastTriMesh<Traits>::edgeVertID(index_t e_id, index_t off) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	OMC_EXPENSIVE_ASSERT(off < 2, "offset out of range.");
	return reinterpret_cast<const index_t *>(&edges[e_id].verts)[off];
}

template <typename Traits>
Sign FastTriMesh<Traits>::triOrientation(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range");
	Sign          ori = Sign::ZERO;
	const GPoint &v0 = triVert(t_id, 0), &v1 = triVert(t_id, 1),
	             &v2 = triVert(t_id, 2);
	switch (refPlane())
	{
	case XY:
		ori = OrientOn2D().on_xy(v0, v1, v2);
		break;
	case YZ:
		ori = OrientOn2D().on_yz(v0, v1, v2);
		break;
	default:
		ori = OrientOn2D().on_zx(v0, v1, v2);
		break;
	}
	OMC_ASSERT(ori != Sign::ZERO, "tri orientation is zero.");
	return ori;
}

template <typename Traits>
void FastTriMesh<Traits>::setEdgeConstr(index_t e_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	edges[e_id].constr = true;
}

template <typename Traits>
bool FastTriMesh<Traits>::edgeIsConstr(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	return edges[e_id].constr;
}

template <typename Traits>
bool FastTriMesh<Traits>::edgeIsBoundary(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	return e2t[e_id].size() == 1;
}

template <typename Traits>
bool FastTriMesh<Traits>::edgeIsManifold(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range");
	return e2t[e_id].size() <= 2;
}

template <typename Traits>
Plane FastTriMesh<Traits>::refPlane() const
{
	return triangle_plane;
}

template <typename Traits>
Sign FastTriMesh<Traits>::Orientation() const
{
	return orientation;
}

template <typename Traits>
size_t FastTriMesh<Traits>::vertValence(index_t v_id) const
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	return v2t[v_id].size();
}

template <typename Traits>
index_t FastTriMesh<Traits>::rotateAroundVertex(index_t center_v_id,
                                                index_t start_t_id, size_t step,
                                                bool CCW)
{
	/// @brief Step to next triangle
	/// @retval Valid index of the next triangle.
	/// @retval Invalid index if cross edge is boundary edge.
	auto step2NextTriangle =
	  [this](index_t center_v_id, index_t curr_t_id, bool CCW)
	{
		index_t center_v_off = triVertOffset(curr_t_id, center_v_id);
		index_t cross_e_off  = CCW ? (center_v_off + 2) % 3 : center_v_off;
		index_t cross_e_id   = triEdgeID(curr_t_id, cross_e_off);
		return triOppToEdge(cross_e_id, curr_t_id);
	};

	/// @brief Step multiple triangles, may stop at boundary.
	/// @return Pair of <index, size>, the first index is the current triangle's
	/// index, the second size is the remain steps.
	auto stepMultTriangles = [this, &step2NextTriangle](index_t center_v_id,
	                                                    index_t start_t_id,
	                                                    size_t step, bool CCW)
	{
		index_t curr_t_id = start_t_id;
		size_t  step_cnt  = 0;
		// step along the direction defined by CCW
		while (step_cnt < step)
		{
			index_t next_t_id = step2NextTriangle(center_v_id, curr_t_id, CCW);
			if (!is_valid_idx(next_t_id)) // meet boundary edge
				break;
			curr_t_id = next_t_id;
			step_cnt++;
		}
		if (step_cnt == step) // step to the target triangle
			return std::pair<index_t, size_t>(curr_t_id, 0);
		// else we meet boundary edge, walk back to the other boundary.
		// (we assume that mesh is manifold, thus a vertex at most connects to two
		// boundary edges.)
		while (true)
		{
			index_t next_t_id = step2NextTriangle(center_v_id, curr_t_id, !CCW);
			if (!is_valid_idx(next_t_id)) // meet another boundary edge
				break;
			curr_t_id = next_t_id;
		}
		step_cnt++;
		return std::pair<index_t, size_t>(curr_t_id, step - step_cnt);
	};

	index_t curr_t_id   = start_t_id;
	size_t  remain_step = step;
	while (remain_step > 0)
	{
		std::tie(curr_t_id, remain_step) =
		  stepMultTriangles(center_v_id, curr_t_id, remain_step, CCW);
	}
	return curr_t_id;
}

template <typename Traits>
template <typename Container>
void FastTriMesh<Traits>::removeFromVec(Container &vec, index_t elem)
{
	auto it = std::find(vec.begin(), vec.end(), elem);
	if (it != vec.end())
	{
		*it = vec.back();
		vec.pop_back();
		OMC_EXPENSIVE_ASSERT((std::find(vec.begin(), vec.end(), elem) == vec.end()),
		                     "same elem remains and is not be removed.");
	}
	else
	{
		OMC_ASSERT(false, "should not happen.");
	}
}

template <typename Traits>
template <typename Container>
void FastTriMesh<Traits>::updateInVec(Container &vec, index_t old_id,
                                      index_t new_id)
{
	auto it = std::find(vec.begin(), vec.end(), old_id);
	OMC_EXPENSIVE_ASSERT(it != vec.end(), "can't find old id");
	*it = new_id;
	OMC_EXPENSIVE_ASSERT((std::find(vec.begin(), vec.end(), old_id) == vec.end()),
	                     "old id remains and is not deleted.");
}

template <typename Traits>
void FastTriMesh<Traits>::removeEdgeUnref(index_t e_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge id out of range.");
	index_t last_eid = edges.size() - 1;
	if (e_id != last_eid) // e_id is a middle element
	{
		// gonna move last edge to edges[e_id]

		// but before moving, update v2e...
		updateInVec(v2e[edges[last_eid].verts.first], last_eid, e_id);
		updateInVec(v2e[edges[last_eid].verts.second], last_eid, e_id);
		// ...and t2e (iTri.edges)
		for (index_t t_id : e2t[last_eid])
			updateInVec(triangles[t_id].edges, last_eid, e_id);
		// finnaly update edges and e2t
		edges[e_id] = std::move(edges[last_eid]);
		e2t[e_id]   = std::move(e2t[last_eid]);
	}
	edges.pop_back();
	e2t.pop_back();
}

template <typename Traits>
void FastTriMesh<Traits>::removeTriUnref(index_t t_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTriangles(), "tri id out of range.");
	index_t last_tid = triangles.size() - 1;
	if (t_id != last_tid) // t_id is a middle element
	{
		// gonna move last triangle to triangles[t_id]

		// but before moving, update v2t...
		updateInVec(v2t[triangles[last_tid].verts[0]], last_tid, t_id);
		updateInVec(v2t[triangles[last_tid].verts[1]], last_tid, t_id);
		updateInVec(v2t[triangles[last_tid].verts[2]], last_tid, t_id);
		// ...and e2t.
		updateInVec(e2t[triangles[last_tid].edges[0]], last_tid, t_id);
		updateInVec(e2t[triangles[last_tid].edges[1]], last_tid, t_id);
		updateInVec(e2t[triangles[last_tid].edges[2]], last_tid, t_id);
		// finnaly update triangles
		triangles[t_id] = std::move(triangles[last_tid]);
	}
	triangles.pop_back();
}

} // namespace OMC