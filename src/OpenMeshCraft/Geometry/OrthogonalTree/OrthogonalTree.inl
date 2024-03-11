#pragma once

#include "OrthogonalTree.h"

#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

#include "tbb/tbb.h"

#include <algorithm>
#include <deque>
#include <execution>
#include <queue>

namespace OMC {

template <typename Traits>
template <typename Primitives, typename Indices>
void OrthogonalTree<Traits>::insert_primitives(const Primitives &primitives,
                                               const Indices    &indices)
{
	OMC_THROW_INVALID_ARGUMENT_IF(
	  primitives.size() != indices.size(),
	  "size of primitives and indices are different.");

	clear();

	m_boxes.resize(primitives.size());
	tbb::parallel_for(size_t(0), primitives.size(),
	                  [this, &primitives, &indices](size_t i)
	                  {
		                  m_boxes[i].bbox() = m_calc_bbox(primitives[i]);
		                  m_boxes[i].id()   = static_cast<index_t>(indices[i]);
	                  });
}

template <typename Traits>
template <typename Bboxes, typename Indices>
void OrthogonalTree<Traits>::insert_boxes(const Bboxes  &bboxes,
                                          const Indices &indices)
{
	static_assert(std::is_same_v<remove_cvref_t<decltype(*bboxes.begin())>, Bbox>,
	              "Bounding box types are different.");
	OMC_THROW_INVALID_ARGUMENT_IF(
	  bboxes.size() != indices.size(),
	  "size of primitives and indices are different.");

	clear();

	m_boxes.resize(bboxes.size());
	tbb::parallel_for(size_t(0), bboxes.size(),
	                  [this, &bboxes, &indices](size_t i)
	                  {
		                  m_boxes[i].bbox() = bboxes[i];
		                  m_boxes[i].id()   = indices[i];
	                  });
}

template <typename Traits>
void OrthogonalTree<Traits>::construct(bool compact_box, NT enlarge_ratio,
                                       NT dupl_thres, index_t depth_delta)
{
	// save behavior control flags and data for further use or query.
	m_enlarge_ratio = enlarge_ratio;
	m_dupl_thres    = dupl_thres;
	m_depth_delta   = depth_delta;
	// Find the tight bounding box of boxes
	m_bbox          = calc_box_from_boxes()(m_boxes.begin(), m_boxes.end());

	// Find the side length of box
	OrPoint bbox_center = (m_bbox.max_bound() + m_bbox.min_bound()) * NT(0.5);
	OrPoint side_length = m_bbox.max_bound() - m_bbox.min_bound();
	OMC_ASSERT(enlarge_ratio >= NT(1.), "enlarge ratio must be larger than 1.");
	if (!compact_box)
	{
		index_t longest_axis = m_bbox.longest_axis();
		side_length          = OrPoint(side_length[longest_axis]);
	}
	side_length *= enlarge_ratio;

	// Reshape the box to a regular square / cube (all side lengths are equal).
	m_bbox.min_bound() = bbox_center - side_length * NT(0.5);
	m_bbox.max_bound() = bbox_center + side_length * NT(0.5);

	// store the side length at all depth.
	m_side_length_per_depth.clear();
	m_side_length_per_depth.reserve(MaxDepth);
	m_side_length_per_depth.push_back(side_length);
	for (index_t i = 1; i < MaxDepth; i++)
		m_side_length_per_depth.push_back(m_side_length_per_depth.back() * 0.5);

	// store boxes into root node
	m_nodes.clear();
	m_nodes.emplace_back();
	root_node().box() = m_bbox;
	root_node().boxes().resize(m_boxes.size());
	root_node().size() = m_boxes.size();
	std::transform(std::execution::par_unseq, m_boxes.begin(), m_boxes.end(),
	               root_node().boxes().begin(), [](OrBbox &box) { return &box; });

	// initialize assign count if neccessary
	m_assign_cnt.clear();
	m_assign_cnt.resize(m_boxes.size(), 1);

	// Initialize a queue of nodes that need to be split
	// (width-first, lower depth first)
	std::queue<index_t> nodes_to_split;
	nodes_to_split.push(m_root_idx);

	// Split nodes recursively.
	while (!nodes_to_split.empty())
	{
		index_t cur_node_idx = nodes_to_split.front();
		nodes_to_split.pop();

		NodeRef cur_node = node(cur_node_idx);

		// Check if this node needs to be splitted
		// TODO: A split predicate need receive what?
		if (cur_node.depth() < MaxDepth - 1 && m_split_pred(*this, cur_node))
		{
			// Split the node, redistributing its boxes to its children
			split(cur_node_idx);

			// post check
			if (cur_node.dupl_degree() > dupl_thres)
			{ // if duplication degree is too large, collapse it.
				collapse(cur_node_idx);
			}
			else
			{ // process each of its children
				for (index_t i = 0; i < Degree; ++i)
					nodes_to_split.push(cur_node.child(i));
			}
		}
	}

	// refines the tree such that the difference of depth between two
	// neighbor leaves is never more than depth_delta
	grade(dupl_thres, depth_delta);

	if constexpr (EnableVertices)
	{
		build_vertices();
	}
}

template <typename Traits>
void OrthogonalTree<Traits>::clear()
{
	m_nodes.clear();
	m_bbox = Bbox();
	m_side_length_per_depth.clear();
	m_boxes.clear();
	m_assign_cnt.clear();
	m_vertices.clear();
}

template <typename Traits>
void OrthogonalTree<Traits>::grade(NT dupl_thres, index_t depth_delta)
{
	if (root_node().is_leaf())
		return;
	if (depth_delta > MaxDepth)
		return;

	// collect all the leaf nodes
	std::queue<index_t> leaf_nodes;
	for (index_t idx = 0; idx < m_nodes.size(); idx++)
	{
		if (node(idx).is_leaf())
			leaf_nodes.push(idx);
	}

	// Iterate over the nodes
	while (!leaf_nodes.empty())
	{
		index_t cur_node_idx = leaf_nodes.front();
		NodeRef cur_node     = node(cur_node_idx);
		leaf_nodes.pop();

		// Skip current node if it isn't a leaf.
		if (!cur_node.is_leaf())
			continue;

		auto CheckNeigbor = [this, &leaf_nodes, &depth_delta, &dupl_thres,
		                     &cur_node](index_t dim, bool dir) -> void
		{
			index_t neighbor_idx = adjacent_node(cur_node, dim, dir);

			if (!is_valid_idx(neighbor_idx)
			    //^ this neighbor does not exist
			    || !node(neighbor_idx).is_leaf()
			    //^ neighbor is already been split
			    || (cur_node.depth() - node(neighbor_idx).depth()) <= depth_delta
			    //^ different of depth is less than depth_delta
			    || cur_node.depth() >= MaxDepth
			    //^ depth is larger than MaxDepth
			)
				return;

			split(neighbor_idx);

			// post check
			if (node(neighbor_idx).dupl_degree() > dupl_thres)
			{ // if duplication degree is too large, collapse it.
				collapse(neighbor_idx);
			}
			else
			{ // process each of its children
				for (index_t i = 0; i < Degree; i++)
					leaf_nodes.push(node(neighbor_idx).child(i));
			}
		};

		// traverse its neighbors
		for (index_t dim = 0; dim < Dimension; dim++)
		{
			CheckNeigbor(/*dimension*/ dim, /*direction*/ false);
			CheckNeigbor(/*dimension*/ dim, /*direction*/ true);
		}
	}
}

template <typename Traits>
auto OrthogonalTree<Traits>::new_children() -> index_t
{
	index_t first_idx = (index_t)m_nodes.size();
	m_nodes.insert(m_nodes.end(), Degree, Node());
	return first_idx;
}

template <typename Traits>
auto OrthogonalTree<Traits>::new_vertex() -> index_t
{
	index_t idx = (index_t)m_vertices.size();
	m_vertices.emplace_back();
	return idx;
}

template <typename Traits>
void OrthogonalTree<Traits>::child_inherit_parent(
  index_t child_idx, index_t parent_idx, const LocalCoordinates &local_coord)
{
	NodeRef child = node(child_idx);

	child.parent() = parent_idx;
	if (is_valid_idx(parent_idx))
	{
		NodeRef parent = node(parent_idx);
		child.depth()  = parent.depth() + 1;
		for (index_t i = 0; i < Dimension; i++)
		{
			child.global_coordinates()[i] = parent.global_coordinates()[i] << 1;
			child.global_coordinates()[i] += local_coord[i];
		}
	}
	else
	{
		child.depth() = 0;
		for (index_t i = 0; i < Dimension; i++)
			child.global_coordinates()[i] = 0;
	}
}

template <typename Traits>
void OrthogonalTree<Traits>::split(index_t node_idx)
{
	// split node to children node
	NodeRef nd = node(node_idx);

	nd.children() = new_children();
	// <--- This is not thread-safe!
	//  --- put it outside split() and lock it by mutex.

	for (index_t index = 0; index < Degree; index++)
		child_inherit_parent(nd.child(index), node_idx, LocalCoordinates(index));

	// get the center to split boxes
	OrPoint center = node_center(nd);

	// calculate box for children
	calc_box_for_children(nd, center);

	// reassign boxes to children
	assign_boxes(nd, center);

	// clear
	if constexpr (!StoreBoxesInInternalNodes)
	{
		nd.boxes().clear();
		nd.boxes().shrink_to_fit();
	}
}

template <typename Traits>
void OrthogonalTree<Traits>::collapse(index_t node_idx)
{
	OMC_EXPENSIVE_ASSERT(node(node_idx).is_internal(),
	                     "only allow to collapse internal node");
	// get node and its children
	NodeRef nd = node(node_idx);
	index_t ch = nd.children();

	// clang-format off
	OMC_EXPENSIVE_ASSERT_AUX_CODE(for (size_t i = 0; i < Degree; i++) {)
	  OMC_EXPENSIVE_ASSERT(
	    node(ch + i).is_leaf(), "only allow to collapse parent of leaf nodes.");
	OMC_EXPENSIVE_ASSERT_AUX_CODE(});
	// clang-format on

	// collect data in children to this node
	if constexpr (!StoreBoxesInInternalNodes)
	{
		auto &boxes = nd.boxes();
		boxes.clear();
		boxes.reserve(nd.size());
		for (size_t i = 0; i < Degree; i++)
		{
			auto &ch_boxes = node(ch + i).boxes();
			boxes.insert(boxes.end(), ch_boxes.begin(), ch_boxes.end());
		}
		std::sort(boxes.begin(), boxes.end(),
		          [](OrBboxCPtr lhs, OrBboxCPtr rhs)
		          { return lhs->id() < rhs->id(); });
		boxes.erase(std::unique(boxes.begin(), boxes.end(),
		                        [](OrBboxCPtr lhs, OrBboxCPtr rhs)
		                        { return lhs->id() == rhs->id(); }),
		            boxes.end());
	}

	// delete its children
	index_t last_ch = m_nodes.size() - Degree;
	if (ch != last_ch)
	{
		// copy last children to ch
		node(m_nodes[last_ch].parent()).children() = ch;
		for (size_t i = 0; i < Degree; i++)
			node(ch + i) = std::move(node(last_ch + i));
	}
	// erase last Degree children
	m_nodes.erase(m_nodes.begin() + last_ch, m_nodes.end());
	// invalidate child index
	nd.children() = InvalidIndex;
}

template <typename Traits>
void OrthogonalTree<Traits>::assign_boxes(NodeRef nd, OrPointCRef center)
{
	for (index_t i = 0; i < Degree; i++)
	{
		node(nd.child(i)).boxes().clear();
		node(nd.child(i)).boxes().reserve(nd.boxes().size() / Degree);
	}

	auto assign_box = [this, &nd, &center](OrBboxCPtr box_ptr)
	{
		size_t assigned_cnt = 0;
#if 1 // this implementation is possily faster?
		std::pair<std::bitset<Dimension>, std::bitset<Dimension>> res =
		  compare_box_with_center(*box_ptr, center);
		size_t rf = res.first.to_ulong();
		size_t rs = res.second.to_ulong();
		// for each child, check if assign box to it.
		for (size_t i = 0; i < Degree; i++)
		{
			// if box overlaps part of child node.
			if ((((i ^ rf) | (~(i ^ rs))) & (Degree - 1)) == Degree - 1)
			{
				// assign box to this child node.
				node(nd.child(i)).boxes().push_back(box_ptr);
				assigned_cnt += 1;
				m_assign_cnt[box_ptr - m_boxes.data()] += 1;
			}
		}
#else
		// for each child, check if assign box to it.
		for (index_t i = 0; i < Degree; i++)
		{
			if (DoIntersect()(node(nd.child(i)).box(), box_ptr->bbox()))
			{
				// assign box to this child node.
				node(nd.child(i)).boxes().push_back(box_ptr);
				assigned_cnt += 1;
				m_assign_cnt[box_ptr - m_boxes.data()] += 1;
			}
		}
#endif
		OMC_EXPENSIVE_ASSERT(assigned_cnt > 0,
		                     "A box is not assigned to any child.");
		// after being assigned, reduce assign count by 1 because this box
		// is not assigned to parent node any more.
		m_assign_cnt[box_ptr - m_boxes.data()] -= 1;
	};

	std::for_each(nd.boxes().begin(), nd.boxes().end(), assign_box);

	size_t total_assigned_cnt =
	  std::accumulate(nd.boxes().begin(), nd.boxes().end(), size_t(0),
	                  [this](size_t cnt, OrBboxCPtr box_ptr)
	                  { return cnt + m_assign_cnt[box_ptr - m_boxes.data()]; });
	nd.dupl_degree() = (float)total_assigned_cnt / (float)nd.size();

	for (index_t i = 0; i < Degree; i++)
		node(nd.child(i)).size() = node(nd.child(i)).boxes().size();
}

template <typename Traits>
auto OrthogonalTree<Traits>::compare_box_with_center(OrBboxCRef  box,
                                                     OrPointCRef center)
  -> std::pair<std::bitset<Dimension>, std::bitset<Dimension>>
{
	std::pair<std::bitset<Dimension>, std::bitset<Dimension>> res;
	for (index_t i = 0; i < Dimension; i++)
	{
		res.first[i]  = box.min_coord(i) < center[i];
		res.second[i] = box.max_coord(i) >= center[i];
	}
	return res;
}

template <typename Traits>
auto OrthogonalTree<Traits>::node_center(NodeCRef nd) const -> OrPoint
{
	OrPoint side_length = m_side_length_per_depth[nd.depth()];

	// Get the relative position of the center
	OrPoint c;
	for (index_t i = 0; i < Dimension; i++)
		c[i] = nd.global_coordinates()[i] * side_length[i] + side_length[i] * 0.5;

	// Add min bound of box to get the absolute position of the center
	c += m_bbox.min_bound();

	return c;
}

template <typename Traits>
auto OrthogonalTree<Traits>::node_side_length(NodeCRef nd) const -> OrPoint
{
	return m_side_length_per_depth[nd.depth()];
}

template <typename Traits>
std::vector<index_t> OrthogonalTree<Traits>::all_leaf_nodes() const
{
	std::vector<index_t> res;
	res.reserve(m_nodes.size());
	for (size_t i = 0; i < m_nodes.size(); i++)
	{
		if (node(i).is_leaf())
			res.push_back(i);
	}
	return res;
}

template <typename Traits>
auto OrthogonalTree<Traits>::adjacent_node(NodeCRef nd, index_t dim,
                                           bool dir) const -> index_t
{
	// The root node has no adjacent nodes.
	if (nd.is_root())
		return InvalidIndex;

	// Check if this child has the adjacent sibling along the direction
	if (nd.local_coordinates(dim) != dir)
	{
		// This means the adjacent node is a direct sibling, return it.
		return node(nd.parent())
		  .child(nd.local_coordinates().set(dim, dir).to_ulong());
	}

	// Find the parent's neighbor in that direction if it exists
	index_t adj_node_of_parent = adjacent_node(node(nd.parent()), dim, dir);

	// If the parent has no neighbor, then this node doesn't have one
	if (!is_valid_idx(adj_node_of_parent))
		return InvalidIndex;

	// If the parent's adjacent node has no children, then it's this node's
	// adjacent node
	if (node(adj_node_of_parent).is_leaf())
		return adj_node_of_parent;

	// Return the nearest node of the parent by reverse the direction.
	return node(adj_node_of_parent)
	  .child(nd.local_coordinates().set(dim, !dir).to_ulong());
}

template <typename Traits>
bool OrthogonalTree<Traits>::is_topologically_same(NodeCRef lhs,
                                                   NodeCRef rhs) const
{
	// If one node is a leaf, and the other isn't, they're not the same
	if (lhs.is_leaf() != rhs.is_leaf())
		return false;

	// If both nodes are non-leaf nodes
	if (!lhs.is_leaf())
	{
		// Check all the children
		for (index_t i = 0; i < Degree; ++i)
		{
			// If any child cell is different, they're not the same
			if (!is_topologically_same(node(lhs.child(i)), node(rhs.child(i))))
				return false;
		}
	}

	// If both nodes are leaf nodes, they must have same topology.
	// WARN: std::array::operator== until C++20.
	return (lhs.global_coordinates() == rhs.global_coordinates());
}

template <typename Traits>
template <typename PointT>
auto OrthogonalTree<Traits>::locate(const PointT &point) const -> index_t
{
	// Make sure the point is enclosed by the orthtree
	OMC_THROW_INVALID_ARGUMENT_IF(
	  !m_do_intersect(point, m_bbox),
	  "point does not locate in the tree's bounding box.");

	// Start at the root node
	index_t node_for_point = m_root_idx;

	// Descend the tree until reaching a leaf node
	while (!node(node_for_point).is_leaf())
	{
		NodeCRef nd = node(node_for_point);

		// Find the point to split around
		OrPoint center = node_center(nd);

		// Find the index of the correct sub-node
		LocalCoordinates index;
		for (index_t dim = 0; dim < Dimension; dim++)
			// Lower: 0(false), higher: 1(true).
			index[dim] = (point[dim] > center[dim]);

		// Find the correct sub-node of the current node
		node_for_point = nd.child(index.to_ulong());
	}

	// Return the result
	return node_for_point;
}

template <typename Traits>
template <typename TraversalTrait>
void OrthogonalTree<Traits>::traversal(TraversalTrait &traits) const
{
	traversal_node(root_node(), traits);
}

template <typename Traits>
template <typename TraversalTrait>
bool OrthogonalTree<Traits>::traversal_node(NodeCRef        nd,
                                            TraversalTrait &traits) const
{
	bool go_next = true;
	if (traits.do_inter(nd.box()))
	{
		if (nd.is_internal())
		{
			// traversal children
			for (index_t i = 0; i < Degree; i++)
			{
				if (go_next)
					go_next = traversal_node(node(nd.child(i)), traits);
				else
					break;
			}
		}
		else // leaf node
		{
			// traversal saved boxes
			for (OrBboxCPtr box_ptr : nd.boxes())
			{
				if (go_next)
					go_next = traits.intersection(*box_ptr);
				else
					break;
			}
		}
	}
	return go_next;
}

} // namespace OMC