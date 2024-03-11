#pragma once

#include "AABBTree.h"

namespace OMC {

template <typename Traits>
AABBTree<Traits>::AABBTree(Prims &&primitives)
{
	insert(std::move(primitives));
	build();
}

template <typename Traits>
AABBTree<Traits>::AABBTree(PrimsIter first, PrimsIter beyond)
{
	insert(first, beyond);
	build();
}

template <typename Traits>
AABBTree<Traits>::~AABBTree()
{
	clear();
}

template <typename Traits>
void AABBTree<Traits>::insert(const Prims &primitives)
{
	clear();
	m_primitives = primitives;
}

template <typename Traits>
void AABBTree<Traits>::insert(Prims &&primitives)
{
	clear();
	m_primitives = std::move(primitives);
}

template <typename Traits>
void AABBTree<Traits>::insert(PrimsIter first, PrimsIter beyond)
{
	clear();
	m_primitives.insert(m_primitives.end(), first, beyond);
}

template <typename Traits>
void AABBTree<Traits>::build()
{
	// clear nodes and avoid wasting memory
	m_nodes.clear();
	m_nodes.shrink_to_fit();

	if (m_primitives.size() > 1)
	{
		// we gonna construct a certain number of nodes
		m_nodes.resize(m_primitives.size() - 1);
		// construct AABB tree.
		expand(m_nodes.data(), m_primitives.begin(), m_primitives.end());
	}
}

template <typename Traits>
void AABBTree<Traits>::clear()
{
	m_nodes.clear();
	m_primitives.clear();
}

template <typename Traits>
template <typename TraversalTrait>
void AABBTree<Traits>::traversal(TraversalTrait &traits) const
{
	switch (size())
	{
	case 0:
		break;
	case 1:
		traits.intersection(m_primitives[0]);
		break;
	default:
		traversal_node(m_nodes.data(), traits, m_primitives.size());
	}
}

template <typename Traits>
void AABBTree<Traits>::expand(NodePtr node, PrimsIter first, PrimsIter beyond)
{
	CalcBbox   calc_bbox;
	SplitPrims split_prims;
	PrimsIter  middle;

	// Calculate box for primitves(first, beyond).
	node->bbox() = calc_bbox(*first);
	for (PrimsIter second = first + 1; second != beyond; second++)
		node->bbox() += calc_bbox(*second);
	// Split primitives to two groups.
	split_prims(first, beyond, middle, node->bbox());

	size_t range = beyond - first;
	switch (range)
	{
	case 2:
		node->left_ptr()  = &(*first);
		node->right_ptr() = &(*(first + 1));
		break;
	case 3:
		node->left_ptr()  = &(*first);
		node->right_ptr() = node + 1;
		expand(node->right_child(), first + 1, beyond);
		break;
	default:
		size_t new_range  = middle - first;
		node->left_ptr()  = node + 1;
		node->right_ptr() = node + new_range;
		expand(node->left_child(), first, middle);
		expand(node->right_child(), middle, beyond);
		break;
	}
}

template <typename Traits>
template <typename TraversalTrait>
bool AABBTree<Traits>::traversal_node(NodeCPtr node, TraversalTrait &trait,
                                      const size_t nb_primitives) const
{
	bool go_next = true;
	switch (nb_primitives)
	{
	case 2:
		go_next = trait.intersection(*node->left_data());
		if (go_next)
			go_next = trait.intersection(*node->right_data());
		return go_next;
	case 3:
		go_next = trait.intersection(*node->left_data());
		if (go_next && trait.do_inter(node->right_child()->bbox()))
			go_next = traversal_node(node->right_child(), trait, 2);
		return go_next;
	default:
		if (trait.do_inter(node->left_child()->bbox()))
		{
			go_next = traversal_node(node->left_child(), trait, nb_primitives / 2);
			if (go_next && trait.do_inter(node->right_child()->bbox()))
				go_next = traversal_node(node->right_child(), trait,
				                         nb_primitives - nb_primitives / 2);
			return go_next;
		}
		else if (trait.do_inter(node->right_child()->bbox()))
		{
			return traversal_node(node->right_child(), trait,
			                      nb_primitives - nb_primitives / 2);
		}
		return true;
	}
}

} // namespace OMC