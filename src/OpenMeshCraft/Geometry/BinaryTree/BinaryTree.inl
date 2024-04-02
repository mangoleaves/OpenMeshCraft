#pragma once

#include "BinaryTree.h"

#include <algorithm>
#include <execution>

namespace OMC {

template <typename Traits>
void BinaryTree<Traits>::shallow_copy(const BinaryTree &rhs)
{
	// copy nodes
	m_nodes.clear();
	m_nodes.resize(rhs.m_nodes.size());
	tbb::parallel_for(size_t(0), rhs.m_nodes.size(), [this, &rhs](size_t i)
	                  { m_nodes[i].shallow_copy(rhs.m_nodes[i]); });
	// copy box
	m_bbox              = rhs.m_bbox;
	/// Split predicate.
	m_split_pred        = rhs.m_split_pred;
	/// Shape refine predicate.
	m_shape_refine_pred = rhs.m_shape_refine_pred;
	/// Do intersect
	m_do_intersect      = rhs.m_do_intersect;
	/// Calculate boungding box
	m_calc_bbox         = rhs.m_calc_bbox;
	// copy other attributes
	m_enlarge_ratio     = rhs.m_enlarge_ratio;
}

template <typename Traits>
template <typename Primitives, typename Indices>
void BinaryTree<Traits>::insert_primitives(const Primitives &primitives,
                                           const Indices    &indices)
{
	OMC_THROW_INVALID_ARGUMENT_IF(
	  primitives.size() != indices.size(),
	  "size of primitives and indices are different.");

	clear();

	// m_boxes.reserve((size_t)(primitives.size() * 1.2));
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
void BinaryTree<Traits>::insert_boxes(const Bboxes  &bboxes,
                                      const Indices &indices)
{
	static_assert(std::is_same_v<remove_cvref_t<decltype(*bboxes.begin())>, Bbox>,
	              "Bounding box types are different.");
	OMC_THROW_INVALID_ARGUMENT_IF(
	  bboxes.size() != indices.size(),
	  "size of primitives and indices are different.");

	clear();

	// m_boxes.reserve((size_t)(bboxes.size() * 1.2));
	m_boxes.resize(bboxes.size());
	tbb::parallel_for(size_t(0), bboxes.size(),
	                  [this, &bboxes, &indices](size_t i)
	                  {
		                  m_boxes[i].bbox() = bboxes[i];
		                  m_boxes[i].id()   = indices[i];
	                  });
}

template <typename Traits>
void BinaryTree<Traits>::construct(bool compact_box, NT enlarge_ratio)
{
	// save behavior control flags and data for further use or query.
	m_enlarge_ratio = enlarge_ratio;
	// Find the tight bounding box of input boxes
	m_bbox          = calc_box_from_boxes()(m_boxes.begin(), m_boxes.end());
	// set the tight box to root node
	m_nodes.clear();
	m_nodes.emplace_back();
	root_node().tbox() = m_bbox;

	// Find the side length of box
	TreePoint bbox_center = (m_bbox.max_bound() + m_bbox.min_bound()) * NT(0.5);
	TreePoint side_length = m_bbox.max_bound() - m_bbox.min_bound();
	OMC_ASSERT(enlarge_ratio >= NT(1.), "enlarge ratio must be larger than 1.");
	if (!compact_box)
	{
		index_t longest_axis = m_bbox.longest_axis();
		side_length          = TreePoint(side_length[longest_axis]);
	}
	side_length *= enlarge_ratio;

	// Reshape the box to a regular square / cube (all side lengths are equal).
	m_bbox.min_bound() = bbox_center - side_length * NT(0.5);
	m_bbox.max_bound() = bbox_center + side_length * NT(0.5);

	// store boxes into root node
	root_node().box()   = m_bbox;
	root_node().size()  = m_boxes.size();
	root_node().boxes() = typename Node::BboxesContainer(0, m_boxes.size());

	constexpr size_t NUM_NODES_FOR_PARALLEL = 8;

	// Initialize a queue of nodes that need to be split
	std::deque<index_t> nodes_to_split;
	nodes_to_split.push_back(m_root_idx);

	// sequentially split nodes until there is no node to split
	// or we get enough nodes to split them parallelly.
	while (!nodes_to_split.empty() &&
	       nodes_to_split.size() < NUM_NODES_FOR_PARALLEL)
	{
		index_t cur_node_idx = nodes_to_split.front();
		nodes_to_split.pop_front();
		// handle the task
		NodeRef cur_node = node(cur_node_idx);

		// Check if this node needs to be splitted
		if (cur_node.depth() < MaxDepth && m_split_pred(*this, cur_node))
		{
			// Split the node, redistributing its boxes to its children
			if (split(cur_node_idx))
			{
				// process each of its children
				for (index_t i = 0; i < cur_node.children_size(); ++i)
					nodes_to_split.push_back(cur_node.child(i));
			}
		}
	}

	// if there remain enough nodes, split them parallelly.
	if (!nodes_to_split.empty())
	{
		auto split_node = [this](size_t node_idx)
		{
			std::queue<index_t> local_nodes_to_split;
			local_nodes_to_split.push(node_idx);

			while (!local_nodes_to_split.empty())
			{
				index_t cur_node_idx = local_nodes_to_split.front();
				local_nodes_to_split.pop();
				// handle the task
				NodeRef cur_node = node(cur_node_idx);

				// Check if this node needs to be splitted
				if (cur_node.depth() < MaxDepth && m_split_pred(*this, cur_node))
				{
					// Split the node, redistributing its boxes to its children
					if (split(cur_node_idx))
					{
						// process each of its children
						for (index_t i = 0; i < cur_node.children_size(); ++i)
							local_nodes_to_split.push(cur_node.child(i));
					}
				}
			}
		};

		tbb::parallel_for_each(nodes_to_split.begin(), nodes_to_split.end(),
		                       split_node);
	}
}

template <typename Traits>
void BinaryTree<Traits>::clear()
{
	m_nodes.clear();
	m_bbox = Bbox();
	m_boxes.clear();
}

template <typename Traits>
auto BinaryTree<Traits>::new_children(size_t n_children) -> index_t
{
	std::lock_guard<tbb::spin_mutex> lock(m_new_children_mutex);

	index_t first_idx = (index_t)m_nodes.size();
	m_nodes.grow_by(n_children, Node());
	return first_idx;
}

template <typename Traits>
bool BinaryTree<Traits>::split(index_t node_idx)
{
	NodeRef nd = node(node_idx);

	// get where to split boxes
	index_t split_dim;
	NT      split_coord;

	bool do_split = m_split_manner(*this, nd, split_dim, split_coord);
	if (!do_split)
		return false;

	// update the node itself
	nd.split_dim()   = split_dim;
	nd.split_coord() = split_coord;

	// new children
	index_t children_idx = new_children(Degree);
	nd.children()        = children_idx;

	// update lower child
	NodeRef lower_ch  = node(nd.child(0));
	lower_ch.depth()  = nd.depth() + 1;
	lower_ch.parent() = node_idx;
	lower_ch.box()    = nd.tbox();

	// update higher child
	NodeRef higher_ch  = node(nd.child(1));
	higher_ch.depth()  = nd.depth() + 1;
	higher_ch.parent() = node_idx;
	higher_ch.box()    = nd.tbox();

	// update box of children
	lower_ch.box().max_bound()[split_dim]  = split_coord;
	higher_ch.box().min_bound()[split_dim] = split_coord;

	index_t begin = *nd.boxes().begin(), end = *nd.boxes().end();
	// update boxes for nd and its children
	// box cross split_coord or does not cross.
	auto    lower_begin =
	  std::partition(m_boxes.begin() + begin, m_boxes.begin() + end,
	                 [&](const TreeBbox &b)
	                 {
		                 return b.min_coord(split_dim) < split_coord &&
		                        b.max_coord(split_dim) >= split_coord;
	                 });
	// box in lower or higher.
	auto higher_begin =
	  std::partition(lower_begin, m_boxes.begin() + end, [&](const TreeBbox &b)
	                 { return b.min_coord(split_dim) < split_coord; });

	nd.boxes() =
	  typename Node::BboxesContainer(begin, lower_begin - m_boxes.begin());

	lower_ch.boxes() = typename Node::BboxesContainer(
	  lower_begin - m_boxes.begin(), higher_begin - m_boxes.begin());

	higher_ch.boxes() =
	  typename Node::BboxesContainer(higher_begin - m_boxes.begin(), end);

	lower_ch.size()  = lower_ch.boxes().size();
	higher_ch.size() = higher_ch.boxes().size();

	if (m_split_manner.need_tight_box())
	{
		calc_tbox_for_children(nd);
	}
	else
	{
		lower_ch.tbox()  = lower_ch.box();
		higher_ch.tbox() = higher_ch.box();
	}

	return true;
}

template <typename Traits>
void BinaryTree<Traits>::shape_refine()
{
	// Initialize a queue of nodes that need to be split
	std::queue<index_t> nodes_to_split;
	for (index_t nidx = 0; nidx < m_nodes.size(); nidx++)
	{
		if (node(nidx).is_leaf())
		{
			nodes_to_split.push(nidx);
		}
	}

	// Split nodes
	while (!nodes_to_split.empty())
	{
		index_t cur_node_idx = nodes_to_split.front();
		nodes_to_split.pop();

		NodeRef cur_node = node(cur_node_idx);

		// Check if this node needs to be splitted
		if (cur_node.depth() < MaxDepth
		    /*&& m_shape_refine_pred(*this, cur_node)*/)
		{
			// TODO split the node

			// process each of its children
			for (index_t i = 0; i < cur_node.children_size(); ++i)
				nodes_to_split.push(cur_node.child(i));
		}
	}
}

template <typename Traits>
void BinaryTree<Traits>::calc_tbox_for_children(NodeRef nd)
{
	for (index_t chi = 0; chi < nd.children_size(); chi++)
	{
		NodeRef ch = node(nd.child(chi));
		ch.tbox() =
		  calc_box_from_box_indices()(*this, ch.boxes().begin(), ch.boxes().end());
		ch.tbox().min_bound().maximize(ch.box().min_bound());
		ch.tbox().max_bound().minimize(ch.box().max_bound());
	}
}

/**
 * @brief The actual max depth of this tree (NOT the MaxDepth).
 */
template <typename Traits>
size_t BinaryTree<Traits>::depth() const
{
	size_t dep = 0;
	for (NodeCRef nd : m_nodes)
		if (nd.depth() > dep)
			dep = nd.depth();
	return dep;
}

/**
 * @brief get indices of all nodes and correspond height.
 * @return std::vector<std::pair<index_t, size_t>>
 * <index_t: index of node, size_t: height>
 */
template <typename Traits>
std::vector<std::pair<index_t, size_t>>
BinaryTree<Traits>::all_nodes_with_height() const
{
	std::vector<uint8_t>                    node_visited(m_nodes.size(), false);
	std::vector<size_t>                     node_height(m_nodes.size(), 0);
	std::stack<index_t>                     traversal_node;
	std::vector<std::pair<index_t, size_t>> res;

	// depth first traversal
	index_t nd_idx = root_node_idx();
	traversal_node.push(nd_idx);

	while (!traversal_node.empty())
	{
		// get next node
		nd_idx = traversal_node.top();
		traversal_node.pop();

		if (node(nd_idx).is_leaf())
		{
			// calculate height
			node_height[nd_idx] = 0;
			res.push_back(std::pair<index_t, size_t>(nd_idx, node_height[nd_idx]));
		}
		else // node is internal
		{
			if (node_visited[nd_idx])
			{
				// calculate height
				node_height[nd_idx] = 1 + std::max(node_height[node(nd_idx).child(0)],
				                                   node_height[node(nd_idx).child(1)]);
				res.push_back(std::pair<index_t, size_t>(nd_idx, node_height[nd_idx]));
			}
			else
			{
				// continue traversal
				node_visited[nd_idx] = true;
				traversal_node.push(nd_idx);                // last itself
				traversal_node.push(node(nd_idx).child(1)); // then higher
				traversal_node.push(node(nd_idx).child(0)); // first lower
			}
		}
	}
	return res;
}

template <typename Traits>
std::vector<index_t> BinaryTree<Traits>::all_nodes() const
{
	std::vector<index_t> res;
	res.resize(m_nodes.size());
	std::iota(res.begin(), res.end(), 0);
	return res;
}

template <typename Traits>
std::vector<index_t> BinaryTree<Traits>::all_leaf_nodes() const
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
template <typename TraversalTrait>
void BinaryTree<Traits>::traversal(TraversalTrait &traits) const
{
	traversal_node(root_node(), traits);
}

template <typename Traits>
template <typename TraversalTrait>
bool BinaryTree<Traits>::traversal_node(NodeCRef        nd,
                                        TraversalTrait &traits) const
{
	bool go_next = true;
	if (traits.do_inter(nd.tbox()))
	{
		// traversal saved boxes if they are exist
		for (index_t box_idx : nd.boxes())
		{
			if (go_next)
				go_next = traits.intersection(box(box_idx));
			else
				break;
		}
		// traversal children if they exist
		for (index_t i = 0; i < nd.children_size(); i++)
		{
			if (go_next)
				go_next = traversal_node(node(nd.child(i)), traits);
			else
				break;
		}
	}
	return go_next;
}

} // namespace OMC