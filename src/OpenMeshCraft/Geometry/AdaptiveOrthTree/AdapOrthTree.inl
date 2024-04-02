#pragma once

#include "AdapOrthTree.h"

#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

#include "tbb/tbb.h"

#include <algorithm>
#include <deque>
#include <execution>
#include <queue>

namespace OMC {

template <typename Traits>
void AdapOrthTree<Traits>::shallow_copy(const AdapOrthTree &rhs)
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
void AdapOrthTree<Traits>::insert_primitives(const Primitives &primitives,
                                             const Indices    &indices)
{
	OMC_THROW_INVALID_ARGUMENT_IF(
	  primitives.size() != indices.size(),
	  "size of primitives and indices are different.");

	clear();

	m_boxes.reserve((size_t)(primitives.size() * 1.2));
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
void AdapOrthTree<Traits>::insert_boxes(const Bboxes  &bboxes,
                                        const Indices &indices)
{
	static_assert(std::is_same_v<remove_cvref_t<decltype(*bboxes.begin())>, Bbox>,
	              "Bounding box types are different.");
	OMC_THROW_INVALID_ARGUMENT_IF(
	  bboxes.size() != indices.size(),
	  "size of primitives and indices are different.");

	clear();

	m_boxes.reserve((size_t)(bboxes.size() * 1.2));
	m_boxes.resize(bboxes.size());
	tbb::parallel_for(size_t(0), bboxes.size(),
	                  [this, &bboxes, &indices](size_t i)
	                  {
		                  m_boxes[i].bbox() = bboxes[i];
		                  m_boxes[i].id()   = indices[i];
	                  });
}

template <typename Traits>
void AdapOrthTree<Traits>::construct(bool compact_box, NT enlarge_ratio,
                                     NT adaptive_thres)
{
	// save behavior control flags and data for further use or query.
	m_enlarge_ratio  = enlarge_ratio;
	m_adaptive_thres = adaptive_thres;
	// Find the tight bounding box of input boxes
	m_bbox           = calc_box_from_boxes()(m_boxes.begin(), m_boxes.end());
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
	root_node().box() = m_bbox;
	root_node().boxes().resize(m_boxes.size());
	root_node().size() = m_boxes.size();
	std::iota(root_node().boxes().begin(), root_node().boxes().end(), index_t(0));

	constexpr size_t NUM_NODES_FOR_PARALLEL = 10;

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
void AdapOrthTree<Traits>::clear()
{
	m_nodes.clear();
	m_bbox = Bbox();
	m_boxes.clear();
}

template <typename Traits>
auto AdapOrthTree<Traits>::new_children(size_t n_children) -> index_t
{
	std::lock_guard<tbb::spin_mutex> lock(m_new_children_mutex);

	index_t first_idx = (index_t)m_nodes.size();
	m_nodes.grow_by(n_children, Node());
	return first_idx;
}

template <typename Traits>
bool AdapOrthTree<Traits>::split(index_t node_idx)
{
	// split node to children node
	NodeRef nd = node(node_idx);

	// get the center to split boxes
	TreePoint center = node_center(nd);

	// reassign boxes to children
	std::array<std::vector<index_t>, Degree> assign_res;
	std::array<size_t, Dimension>            lower, higher;
	assign_boxes(nd, center, assign_res, lower, higher);

	// check if boxes are partitionalble on all dimensions.
	std::array<bool, Dimension> partitionable;
	for (size_t i = 0; i < Dimension; i++)
	{
		OMC_ASSERT(lower[i] + higher[i] >= nd.size(), "orphan box");
		partitionable[i] =
		  // check if boxes are not partitionable
		  (((NT)(lower[i] + higher[i] - nd.size()) / (NT)nd.size()) <
		   m_adaptive_thres) &&
		  // check if partition is unbalence
		  (lower[i] > 0 && higher[i] > 0);
	}

	if (std::find(partitionable.begin(), partitionable.end(), true) ==
	    partitionable.end())
		return false; // this node is not partitionable.

	// calculate box for children
	std::array<Bbox, Degree> child_boxes;
	calc_box_for_children(nd, center, child_boxes);

	// this node is not partitionable on some dimensions,
	// collapse child nodes and boxes.
	bool need_moving = std::find(partitionable.begin(), partitionable.end(),
	                             false) != partitionable.end();
	if (need_moving)
	{
		// calculate collapse destination.
		std::array<index_t, Degree> collapse_dest;
		collapse_destination(partitionable, collapse_dest);
		// do collapse
		for (size_t i = 0; i < Degree; i++)
		{
			index_t dest = collapse_dest[i];
			if (dest == i)
				continue;

			std::vector<index_t> &dst = assign_res[dest], &src = assign_res[i];
			std::vector<index_t>  tmp = std::move(dst);

			dst.reserve(tmp.size() + src.size());
			merge_unique(tmp.begin(), tmp.end(), src.begin(), src.end(),
			             std::back_inserter(dst), std::less<index_t>(),
			             std::equal_to<index_t>());
			src.clear();
			src.shrink_to_fit();
			child_boxes[dest] += child_boxes[i];
		}
		// allocate new children
		unsigned int children_size =
		  1u << (unsigned int)std::count_if(
		    partitionable.begin(), partitionable.end(), [](bool b) { return b; });

		index_t children_idx = new_children(children_size);
		nd.children()        = children_idx;
		nd.children_size()   = children_size;
		nd.child_map()       = collapse_dest;

		// really assign box, set size
		std::array<index_t, Degree> sorted_dest = collapse_dest;
		std::sort(sorted_dest.begin(), sorted_dest.end());
		OMC_UNUSED auto end_it =
		  std::unique(sorted_dest.begin(), sorted_dest.end());
		for (size_t i = 0; i < children_size; i++)
		{
			NodeRef ch  = node(nd.child(i));
			ch.depth()  = nd.depth() + 1;
			ch.parent() = node_idx;
			ch.box()    = child_boxes[sorted_dest[i]];
			ch.boxes()  = assign_res[sorted_dest[i]];
			ch.size()   = assign_res[sorted_dest[i]].size();
			for (index_t &dest : nd.child_map())
				if (dest == sorted_dest[i])
					dest = i;
		}
	}
	else
	{
		index_t children_idx = new_children(Degree);
		nd.children()        = children_idx;
		nd.children_size()   = Degree;
		std::iota(nd.child_map().begin(), nd.child_map().end(), 0);
		for (size_t i = 0; i < Degree; i++)
		{
			NodeRef ch  = node(nd.child(i));
			ch.depth()  = nd.depth() + 1;
			ch.parent() = node_idx;
			ch.box()    = child_boxes[i];
			ch.boxes()  = assign_res[i];
			ch.size()   = assign_res[i].size();
		}
	}

	calc_tbox_for_children(nd);

	// clear
	if constexpr (!StoreBoxesInInternalNodes)
	{
		nd.boxes() = std::vector<index_t>();
	}

	return true;
}

template <typename Traits>
void AdapOrthTree<Traits>::collapse(index_t node_idx)
{
	OMC_EXPENSIVE_ASSERT(node(node_idx).is_internal(),
	                     "only allow to collapse internal node");
	// get node and its children
	NodeRef nd = node(node_idx);
	index_t ch = nd.children();

	// clang-format off
	OMC_EXPENSIVE_ASSERT_AUX_CODE(
		for (size_t i = 0; i < nd.children_size(); i++) {)
	  OMC_EXPENSIVE_ASSERT(
	    node(ch + i).is_leaf(), "only allow to collapse parent of leaf nodes.");
	OMC_EXPENSIVE_ASSERT_AUX_CODE(});
	// clang-format on

	// collect data in children to this node
	if constexpr (!StoreBoxesInInternalNodes)
	{
		std::vector<index_t> &boxes = nd.boxes();
		boxes.clear();
		boxes.reserve(nd.size());
		for (size_t i = 0; i < nd.children_size(); i++)
		{
			std::vector<index_t> &ch_boxes = node(ch + i).boxes();
			boxes.insert(boxes.end(), ch_boxes.begin(), ch_boxes.end());
			ch_boxes = std::vector<index_t>();
		}
		std::sort(boxes.begin(), boxes.end(), [](TreeBboxCPtr lhs, TreeBboxCPtr rhs)
		          { return lhs->id() < rhs->id(); });
		boxes.erase(std::unique(boxes.begin(), boxes.end(),
		                        [](TreeBboxCPtr lhs, TreeBboxCPtr rhs)
		                        { return lhs->id() == rhs->id(); }),
		            boxes.end());
	}

	// FIXME TODO delete its children, need a garbage collection strategy?

	// invalidate child index
	nd.children() = InvalidIndex;
}

template <typename Traits>
void AdapOrthTree<Traits>::assign_boxes(
  NodeRef nd, TreePointCRef center,
  std::array<std::vector<index_t>, Degree> &assign_res,
  std::array<size_t, Dimension> &lower, std::array<size_t, Dimension> &higher)
{
	std::for_each(assign_res.begin(), assign_res.end(),
	              [&nd](std::vector<index_t> &ass)
	              {
		              ass.clear();
		              ass.reserve(nd.size() / Degree * 2);
	              });
	std::fill(lower.begin(), lower.end(), 0);
	std::fill(higher.begin(), higher.end(), 0);

	auto assign_box =
	  [this, &nd, &center, &assign_res, &lower, &higher](index_t box_idx)
	{
		std::pair<std::bitset<Dimension>, std::bitset<Dimension>> res =
		  compare_box_with_center(box(box_idx), center);

		for (size_t i = 0; i < Dimension; i++)
		{
			lower[i] += res.first[i];
			higher[i] += res.second[i];
		}

		size_t rf = res.first.to_ulong();
		size_t rs = res.second.to_ulong();
		// for each child, check if assign box to it.
		for (size_t i = 0; i < Degree; i++)
		{
			// if box overlaps part of child node.
			if ((((i ^ rf) | (~(i ^ rs))) & (Degree - 1)) == Degree - 1)
			{
				// assign box to this child node.
				assign_res[i].push_back(box_idx);
			}
		}
	};

	std::for_each(nd.boxes().begin(), nd.boxes().end(), assign_box);
}

template <typename Traits>
auto AdapOrthTree<Traits>::compare_box_with_center(TreeBboxCRef  box,
                                                   TreePointCRef center)
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
void AdapOrthTree<Traits>::shape_refine()
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
		// TODO: A split predicate need receive what?
		std::array<bool, Dimension> partitionable;
		if (cur_node.depth() < MaxDepth &&
		    m_shape_refine_pred(*this, cur_node, partitionable))
		{
			TreePoint center = node_center(cur_node);

			// Split the node on partitionable dimensions.
			// Below codes are similar to split(), I copy most of them here.
			std::array<Bbox, Degree> child_boxes;
			calc_box_for_children(cur_node, center, child_boxes);

			bool need_moving = std::find(partitionable.begin(), partitionable.end(),
			                             false) != partitionable.end();

			if (need_moving)
			{
				// calculate collapse destination.
				std::array<index_t, Degree> collapse_dest;
				collapse_destination(partitionable, collapse_dest);
				// do collapse
				for (size_t i = 0; i < Degree; i++)
				{
					index_t dest = collapse_dest[i];
					if (dest == i)
						continue;
					child_boxes[dest] += child_boxes[i];
				}
				// allocate new children
				unsigned int children_size =
				  1u << (unsigned int)std::count_if(partitionable.begin(),
				                                    partitionable.end(),
				                                    [](bool b) { return b; });

				index_t children_idx     = new_children(children_size);
				cur_node.children()      = children_idx;
				cur_node.children_size() = children_size;
				cur_node.child_map()     = collapse_dest;

				std::array<index_t, Degree> sorted_dest = collapse_dest;
				std::sort(sorted_dest.begin(), sorted_dest.end());
				OMC_UNUSED auto end_it =
				  std::unique(sorted_dest.begin(), sorted_dest.end());
				for (size_t i = 0; i < children_size; i++)
				{
					NodeRef ch  = node(cur_node.child(i));
					ch.depth()  = cur_node.depth() + 1;
					ch.parent() = cur_node_idx;
					ch.box()    = child_boxes[sorted_dest[i]];
					ch.tbox()   = ch.box();
					ch.size()   = 0;
					for (index_t &dest : cur_node.child_map())
						if (dest == sorted_dest[i])
							dest = i;
				}
			}
			else
			{
				index_t children_idx     = new_children(Degree);
				cur_node.children()      = children_idx;
				cur_node.children_size() = Degree;
				std::iota(cur_node.child_map().begin(), cur_node.child_map().end(), 0);
				for (size_t i = 0; i < Degree; i++)
				{
					NodeRef ch  = node(cur_node.child(i));
					ch.depth()  = cur_node.depth() + 1;
					ch.parent() = cur_node_idx;
					ch.box()    = child_boxes[i];
					ch.tbox()   = ch.box();
					ch.size()   = 0;
				}
			}
			// process each of its children
			for (index_t i = 0; i < cur_node.children_size(); ++i)
				nodes_to_split.push(cur_node.child(i));
		}
	}
}

template <typename Traits>
void AdapOrthTree<Traits>::calc_tbox_for_children(NodeRef nd)
{
#if 0
	for (index_t chi = 0; chi < nd.children_size(); chi++)
	{
		NodeRef ch = node(nd.child(chi));
		ch.tbox() =
		  calc_box_from_box_indices()(*this, ch.boxes().begin(), ch.boxes().end());
		ch.tbox().min_bound().maximize(ch.box().min_bound());
		ch.tbox().max_bound().minimize(ch.box().max_bound());
	}
#else
	for (index_t chi = 0; chi < nd.children_size(); chi++)
	{
		NodeRef ch = node(nd.child(chi));
		ch.tbox()  = ch.box();
	}
#endif
}

template <typename Traits>
void AdapOrthTree<Traits>::collapse_destination(
  const std::array<bool, Dimension> &partitionable,
  std::array<index_t, Degree>       &destination)
{
	std::iota(destination.begin(), destination.end(), 0);

	for (int i = Dimension - 1; i >= 0; i--)
	{
		if (partitionable[i])
			continue;
		// collapse on this dimension
		size_t offset = size_t(1) << size_t(i);
		for (size_t j = 0; j < Degree; j++)
			if (destination[j] & offset)
				destination[j] -= offset;
	}
}

template <typename Traits>
auto AdapOrthTree<Traits>::node_center(NodeCRef nd) const -> TreePoint
{
	return (nd.tbox().min_bound() + nd.tbox().max_bound()) * 0.5;
}

template <typename Traits>
std::vector<index_t> AdapOrthTree<Traits>::all_leaf_nodes() const
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
void AdapOrthTree<Traits>::traversal(TraversalTrait &traits) const
{
	traversal_node(root_node(), traits);
}

template <typename Traits>
template <typename TraversalTrait>
bool AdapOrthTree<Traits>::traversal_node(NodeCRef        nd,
                                          TraversalTrait &traits) const
{
	bool go_next = true;
	if (traits.do_inter(nd.tbox()))
	{
		if (nd.is_internal())
		{
			// traversal children
			for (index_t i = 0; i < nd.children_size(); i++)
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
			for (index_t box_idx : nd.boxes())
			{
				if (go_next)
					go_next = traits.intersection(box(box_idx));
				else
					break;
			}
		}
	}
	return go_next;
}

// clang-format off
template <typename Traits>
template <typename Iter, typename OutIter, typename LessPred, typename EqualPred>
void AdapOrthTree<Traits>::merge_unique(Iter b1, Iter e1, Iter b2, Iter e2,
                                        OutIter o, LessPred lp, EqualPred ep)
{
	remove_cvref_t<decltype(*b1)> last;
	if (b1 == e1 || b2 == e2)
	{
		if (b1 != e1) { *o = *b1; last = *b1; ++o; ++b1; }
		else if (b2 != e2) { *o = *b2; last = *b2; ++o; ++b2; }
		else return;
	}
	while (b1 != e1 && b2 != e2)
	{
		if (lp(*b1, *b2))
		{
			if (!ep(*b1, last)) { *o = *b1; last = *b1; ++o; }
			++b1;
		}
		else if (lp(*b2, *b1))
		{
			if (!ep(*b2, last)) { *o = *b2; last = *b2; ++o; }
			++b2;
		}
		else
		{
			if (!ep(*b1, last)) { *o = *b1; last = *b1; ++o; }
			++b1; ++b2;
		}
	}
	while(b1 != e1) { *o = *b1; ++o; ++b1; }
	while(b2 != e2) { *o = *b2; ++o; ++b2; }
}
// clang-format on

} // namespace OMC