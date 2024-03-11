#pragma once

#include "OcTree.h"

#include <algorithm>

namespace OMC {

template <typename Traits>
auto OcTree<Traits>::adj_vertex(index_t vidx, index_t dim, bool dir) -> index_t
{
	return m_adj_vertices[vidx][(dim << 1) + dir];
}

template <typename Traits>
auto OcTree<Traits>::calc_vertex_position(index_t nd_idx, index_t nv_idx)
  -> OrPoint
{
	NodeCRef nd = this->node(nd_idx);
	OrPoint  c  = this->node_center(nd);
	OrPoint  h  = this->node_side_length(nd) * NT(0.5);
	switch (nv_idx)
	{
	case 0:
		return OrPoint(c.x() - h.x(), c.y() - h.y(), c.z() - h.z());
	case 1:
		return OrPoint(c.x() + h.x(), c.y() - h.y(), c.z() - h.z());
	case 2:
		return OrPoint(c.x() - h.x(), c.y() + h.y(), c.z() - h.z());
	case 3:
		return OrPoint(c.x() + h.x(), c.y() + h.y(), c.z() - h.z());
	case 4:
		return OrPoint(c.x() - h.x(), c.y() - h.y(), c.z() + h.z());
	case 5:
		return OrPoint(c.x() + h.x(), c.y() - h.y(), c.z() + h.z());
	case 6:
		return OrPoint(c.x() - h.x(), c.y() + h.y(), c.z() + h.z());
	case 7:
		return OrPoint(c.x() + h.x(), c.y() + h.y(), c.z() + h.z());
	default:
#ifdef OMC_ENABLE_ASSERT
		OMC_ASSERT(false, "vertex index {} out of range.", nv_idx);
#else
		return OrPoint(NT(0.), NT(0.), NT(0.));
#endif
	}
}

template <typename Traits>
void OcTree<Traits>::set_adj_vertices_pair(index_t vl, index_t vh, index_t dim)
{
	m_adj_vertices[vl][(dim << 1) + true]  = vh;
	m_adj_vertices[vh][(dim << 1) + false] = vl;
}

template <typename Traits>
void OcTree<Traits>::build_vertices()
{
	this->m_vertices.clear();
	m_adj_vertices.clear();
	// During building vertices, the VertexPtr stored in Node is temporarily
	// interpreted as the index of vertex. This interpretation is used to random
	// access different containers.

	// each node has 6 vertices located on its 6 faces,
	// we store this relation in face_vertices.
	std::vector<std::array<index_t, 2 * Dimension>> face_vertices;
	// each vertex has 6 vertices located on its 6 adjacent edges,
	// we store this relation in edge_vertices.
	std::vector<std::array<index_t, 2 * Dimension>> edge_vertices;

	// initialize face_vertices and edge_vertices
	std::array<index_t, 2 * Dimension> invalid_indices;
	std::fill(invalid_indices.begin(), invalid_indices.end(), InvalidIndex);
	face_vertices.resize(this->m_nodes.size(), invalid_indices);

	// Given a node and its face, we try to find the face vertex on it.
	auto find_face_vertex = [this, &face_vertices, &edge_vertices,
	                         &invalid_indices](index_t nd_idx, index_t dim,
	                                           bool dir) -> index_t
	{
		NodeCRef nd           = this->node(nd_idx);
		// find the adjacent node
		index_t  adj_node_idx = this->adjacent_node(nd, dim, dir);
		// if the adjacent node exists and has the same depth with this node,
		// we try to find face vertex on the same face in adjacent node.

		dim = dim << 1; // NOTE: dim is multiply by 2.
		if (is_valid_idx(adj_node_idx) &&
		    this->node(adj_node_idx).depth() == nd.depth() &&
		    is_valid_idx(face_vertices[adj_node_idx][dim + (!dir)]))
		{
			face_vertices[nd_idx][dim + dir] =
			  face_vertices[adj_node_idx][dim + (!dir)];
		}
		else
		{
			face_vertices[nd_idx][dim + dir] = this->new_vertex();
			edge_vertices.push_back(invalid_indices);
			m_adj_vertices.push_back(invalid_indices);
		}
		return face_vertices[nd_idx][dim + dir];
	};

	// An edge of a node is given by (lower index, higher index),
	// try to find the middle vertex on this edge if it exists.
	auto find_edge_vertex = [this, &edge_vertices, &invalid_indices](
	                          NodeRef nd, index_t li, index_t hi) -> index_t
	{
		// use xor to get the dimension
		index_t dim = ((li ^ hi) >> 1) << 1;
		index_t vl  = nd.vertex(li);
		index_t vh  = nd.vertex(hi);
		// find at higher direction of lower vertex.
		if (is_valid_idx(edge_vertices[vl][dim + true]))
			return edge_vertices[vl][dim + true];
		// find at lower direction of higher vertex. (does this neccessary?)
		if (is_valid_idx(edge_vertices[vh][dim + false]))
			return edge_vertices[vh][dim + false];

		index_t new_vidx = this->new_vertex();
		// set edge vertices
		edge_vertices.push_back(invalid_indices);
		edge_vertices[vl][dim + true]  = new_vidx;
		edge_vertices[vh][dim + false] = new_vidx;
		// set adjacent relationship
		m_adj_vertices.push_back(invalid_indices);
		m_adj_vertices[vl][dim + true]        = new_vidx;
		m_adj_vertices[vh][dim + false]       = new_vidx;
		m_adj_vertices[new_vidx][dim + false] = vl;
		m_adj_vertices[new_vidx][dim + true]  = vh;
		return new_vidx;
	};

	// build 8 corner vertices for root node.
	for (index_t i = 0; i < Degree; i++)
	{
		this->root_node().vertex(i) = this->new_vertex();
	}
	m_adj_vertices.resize(Degree, invalid_indices);
	edge_vertices.resize(Degree, invalid_indices);

	// traversal nodes to build vertices for them.
	std::queue<index_t> nodes_to_traverse;
	index_t             curr_depth;

	nodes_to_traverse.push(this->root_node_idx());
	curr_depth = 0;

	// Iterate nodes
	while (!nodes_to_traverse.empty())
	{
		index_t nd_idx = nodes_to_traverse.front();
		NodeRef nd     = this->node(nd_idx);
		nodes_to_traverse.pop();

		if (curr_depth != nd.depth())
		{
			OMC_ASSERT(curr_depth == nd.depth() - 1,
			           "depth delta is larger than 1 during traverse.");
			curr_depth = nd.depth();
			// reset all edge_vertices
			std::fill(edge_vertices.begin(), edge_vertices.end(), invalid_indices);
		}

		// When traverse to this node, its corner vertices have been built.
		// If it is leaf node, do nothing.
		if (nd.is_leaf())
			continue;
		// else we build face vertices and edge vertices for this node.

		const auto &nd_vs = nd.vertices();
		index_t     v0    = nd_vs[0];
		index_t     v1    = nd_vs[1];
		index_t     v2    = nd_vs[2];
		index_t     v3    = nd_vs[3];
		index_t     v4    = nd_vs[4];
		index_t     v5    = nd_vs[5];
		index_t     v6    = nd_vs[6];
		index_t     v7    = nd_vs[7];

		// try to find face vertices in adjacent nodes,
		// if not found, we build new vertices and store in face_vertices.
		index_t fv0 = find_face_vertex(nd_idx, /*X*/ 0, false);
		index_t fv1 = find_face_vertex(nd_idx, /*X*/ 0, true);
		index_t fv2 = find_face_vertex(nd_idx, /*Y*/ 1, false);
		index_t fv3 = find_face_vertex(nd_idx, /*Y*/ 1, true);
		index_t fv4 = find_face_vertex(nd_idx, /*Z*/ 2, false);
		index_t fv5 = find_face_vertex(nd_idx, /*Z*/ 2, true);

		// try to find edge vertices,
		// if not found, we build new vertices and store in edge_vertices.
		index_t ev01 = find_edge_vertex(nd, 0, 1);
		index_t ev02 = find_edge_vertex(nd, 0, 2);
		index_t ev04 = find_edge_vertex(nd, 0, 4);
		index_t ev13 = find_edge_vertex(nd, 1, 3);
		index_t ev15 = find_edge_vertex(nd, 1, 5);
		index_t ev23 = find_edge_vertex(nd, 2, 3);
		index_t ev26 = find_edge_vertex(nd, 2, 6);
		index_t ev37 = find_edge_vertex(nd, 3, 7);
		index_t ev45 = find_edge_vertex(nd, 4, 5);
		index_t ev46 = find_edge_vertex(nd, 4, 6);
		index_t ev57 = find_edge_vertex(nd, 5, 7);
		index_t ev67 = find_edge_vertex(nd, 6, 7);

		// build center vertices
		index_t cv = this->new_vertex();
		edge_vertices.push_back(invalid_indices);
		m_adj_vertices.push_back(invalid_indices);

		// children nodes inherit vertices.
		this->node(nd.child(0)).vertices() = {v0,   ev01, fv4, ev02,
		                                      ev04, fv2,  cv,  fv0};
		this->node(nd.child(1)).vertices() = {ev01, v1,   ev13, fv4,
		                                      fv2,  ev15, fv1,  cv};
		this->node(nd.child(2)).vertices() = {ev02, fv4, ev23, v2,
		                                      fv0,  cv,  fv3,  ev26};
		this->node(nd.child(3)).vertices() = {fv4, ev13, v3,   ev23,
		                                      cv,  fv1,  ev37, fv3};
		this->node(nd.child(4)).vertices() = {ev04, fv2,  cv,  fv0,
		                                      v4,   ev45, fv5, ev46};
		this->node(nd.child(5)).vertices() = {fv2,  ev15, fv1,  cv,
		                                      ev45, v5,   ev57, fv5};
		this->node(nd.child(6)).vertices() = {fv0,  cv,  fv3,  ev26,
		                                      ev46, fv5, ev67, v6};
		this->node(nd.child(7)).vertices() = {cv,  fv1,  ev37, fv3,
		                                      fv5, ev57, v7,   ev67};

		// set adjancen relationship between
		// (1) between center vertex and all face vertices
		set_adj_vertices_pair(fv0, cv, 0);
		set_adj_vertices_pair(cv, fv1, 0);
		set_adj_vertices_pair(fv2, cv, 1);
		set_adj_vertices_pair(cv, fv3, 1);
		set_adj_vertices_pair(fv4, cv, 2);
		set_adj_vertices_pair(cv, fv5, 2);
		// (2) between face vertices and edge vertices
		set_adj_vertices_pair(ev04, fv0, /*Y*/ 1);
		set_adj_vertices_pair(fv0, ev26, /*Y*/ 1);
		set_adj_vertices_pair(ev02, fv0, /*Z*/ 2);
		set_adj_vertices_pair(fv0, ev46, /*Z*/ 2);
		set_adj_vertices_pair(ev15, fv1, /*Y*/ 1);
		set_adj_vertices_pair(fv1, ev37, /*Y*/ 1);
		set_adj_vertices_pair(ev13, fv1, /*Z*/ 2);
		set_adj_vertices_pair(fv1, ev57, /*Z*/ 2);

		set_adj_vertices_pair(ev04, fv2, /*X*/ 0);
		set_adj_vertices_pair(fv2, ev15, /*X*/ 0);
		set_adj_vertices_pair(ev01, fv2, /*Z*/ 2);
		set_adj_vertices_pair(fv2, ev45, /*Z*/ 2);
		set_adj_vertices_pair(ev26, fv3, /*X*/ 0);
		set_adj_vertices_pair(fv3, ev37, /*X*/ 0);
		set_adj_vertices_pair(ev23, fv3, /*Z*/ 2);
		set_adj_vertices_pair(fv3, ev67, /*Z*/ 2);

		set_adj_vertices_pair(ev02, fv4, /*X*/ 0);
		set_adj_vertices_pair(fv4, ev13, /*X*/ 0);
		set_adj_vertices_pair(ev01, fv4, /*Y*/ 1);
		set_adj_vertices_pair(fv4, ev23, /*Y*/ 1);
		set_adj_vertices_pair(ev46, fv5, /*X*/ 0);
		set_adj_vertices_pair(fv5, ev57, /*X*/ 0);
		set_adj_vertices_pair(ev45, fv5, /*Y*/ 1);
		set_adj_vertices_pair(fv5, ev67, /*Y*/ 1);

		for (index_t i = 0; i < Degree; i++)
			nodes_to_traverse.push(nd.child(i));
	}

	// calculate positions for all vertices
	std::vector<uint8_t> position_calculated(this->m_vertices.size(), false);
	nodes_to_traverse.push(this->root_node_idx());

	while (!nodes_to_traverse.empty())
	{
		index_t nd_idx = nodes_to_traverse.front();
		NodeRef nd     = this->node(nd_idx);
		nodes_to_traverse.pop();

		if (nd.is_leaf())
		{
			for (index_t i = 0; i < Degree; i++)
			{
				index_t vi = nd.vertex(i);
				if (!position_calculated[vi])
				{
					this->vertex(vi).position() = calc_vertex_position(nd_idx, i);
					position_calculated[vi]     = true;
				}
			}
		}
		else
		{
			for (index_t i = 0; i < Degree; i++)
				nodes_to_traverse.push(nd.child(i));
		}
	}
}

template <typename Traits>
void OcTree<Traits>::calc_box_for_children(NodeRef nd, OrPointCRef c)
{
	OrPointCRef minb = nd.box().min_bound();
	OrPointCRef maxb = nd.box().max_bound();

	// clang-format off
	// child 0(000, -x-y-z)
	this->node(nd.child(0)).box() = Bbox(minb, c);
	// child 1(001, +x-y-z)
	this->node(nd.child(1)).box() = Bbox(OrPoint(c.x(), minb.y(), minb.z()), OrPoint(maxb.x(), c.y(), c.z()));
	// child 2(010, -x+y-z)
	this->node(nd.child(2)).box() = Bbox(OrPoint(minb.x(), c.y(), minb.z()), OrPoint(c.x(), maxb.y(), c.z()));
	// child 3(011, +x+y-z)
	this->node(nd.child(3)).box() = Bbox(OrPoint(c.x(), c.y(), minb.z()), OrPoint(maxb.x(), maxb.y(), c.z()));
	// child 4(100, -x-y+z)
	this->node(nd.child(4)).box() = Bbox(OrPoint(minb.x(), minb.y(), c.z()), OrPoint(c.x(), c.y(), maxb.z()));
	// child 5(101, +x-y+z)
	this->node(nd.child(5)).box() = Bbox(OrPoint(c.x(), minb.y(), c.z()), OrPoint(maxb.x(), c.y(), maxb.z()));
	// child 6(110, -x+y+z)
	this->node(nd.child(6)).box() = Bbox(OrPoint(minb.x(), c.y(), c.z()), OrPoint(c.x(), maxb.y(), maxb.z()));
	// child 7(111, +x+y+z)
	this->node(nd.child(7)).box() = Bbox(c, maxb);
	// clang-format on
}

} // namespace OMC