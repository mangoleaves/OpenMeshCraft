#pragma once

#include "KdTraversal.h"

namespace OMC {

template <typename KdTraits>
OrthogonalNearestSeach<KdTraits>::OrthogonalNearestSeach(NodePtr      tree_root,
                                                         const KdBox &tree_bbox,
                                                         const PointT &query)
  : m_dists(NT(0), NT(0), NT(0))
  , m_query(query)
  , m_square_distance(std::numeric_limits<NT>::max())
{
	compute_nearest_neightbors_orthogonally(
	  tree_root, (m_project_point(tree_bbox, m_query) - m_query).sqrnorm());
}

template <typename KdTraits>
void OrthogonalNearestSeach<KdTraits>::compute_nearest_neightbors_orthogonally(
  NodePtr node, NT sqrdis)
{
	if (node->m_is_leaf)
	{
		LeafPtr leaf_node = static_cast<LeafPtr>(node);
		if (leaf_node->n > 0)
		{
			search_nearest_in_leaf(leaf_node);
		}
	}
	else
	{
		InternalPtr internal_node = static_cast<InternalPtr>(node);
		size_t      cut_dim       = internal_node->m_cut_dim;

		NodePtr best_ch, other_ch;

		NT offset;
		NT val   = m_query[cut_dim];
		NT diff1 = val - internal_node->m_upper_low_val;
		NT diff2 = val - internal_node->m_lower_high_val;

		if (diff1 + diff2 < 0)
		{
			offset   = diff1;
			best_ch  = internal_node->m_lower_ch;
			other_ch = internal_node->m_upper_ch;
		}
		else
		{
			offset   = diff2;
			best_ch  = internal_node->m_upper_ch;
			other_ch = internal_node->m_lower_ch;
		}

		compute_nearest_neightbors_orthogonally(best_ch, sqrdis);
		NT dst           = m_dists[cut_dim];
		NT new_rd        = new_square_distance(sqrdis, dst, offset);
		m_dists[cut_dim] = offset;
		if (new_rd < m_square_distance)
		{
			compute_nearest_neightbors_orthogonally(other_ch, new_rd);
		}
		m_dists[cut_dim] = dst;
	}
}

template <typename KdTraits>
void OrthogonalNearestSeach<KdTraits>::search_nearest_in_leaf(LeafPtr node)
{
	for (KdPointsIter begin = node->data, end = node->data + node->n;
	     begin != end; begin++)
	{
		NT distance = ((*begin).primitive() - m_query).sqrnorm();
		if (distance < m_square_distance)
		{
			m_result          = begin;
			m_square_distance = distance;
		}
	}
}

} // namespace OMC