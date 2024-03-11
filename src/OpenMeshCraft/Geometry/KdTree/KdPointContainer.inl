#pragma once

#include "KdPointContainer.h"

#include <algorithm>

namespace OMC {

template <typename KdTraits>
KdPointContainer<KdTraits>::KdPointContainer(KdPointPtrsIter begin,
                                             KdPointPtrsIter end)
  : m_begin(begin)
  , m_end(end)
{
	if (begin != end)
	{
		m_bbox.update_from_point_pointers(begin, end);
	}
	m_tbox      = m_bbox;
	m_build_dim = m_bbox.longest_axis();
}

template <typename KdTraits>
void KdPointContainer<KdTraits>::split(KdPointContainer &c_low, size_t &sep_dim,
                                       NT &sep_val, bool sliding)
{
	c_low.m_bbox      = m_bbox;
	m_build_dim       = sep_dim;
	c_low.m_build_dim = sep_dim;

	KdPointPtrsIter it = std::partition(
	  m_begin, m_end,
	  [&](KdPointPtr v) { return v->primitive()[sep_dim] < sep_val; });
	// now [begin,it) are lower and [it,end) are upper
	if (sliding) // avoid empty lists
	{
		if (it == m_begin)
		{
			KdPointPtrsIter min_elt = std::min_element(
			  m_begin, m_end,
			  [&](KdPointPtr lhs, KdPointPtr rhs)
			  { return lhs->primitive()[sep_dim] < rhs->primitive()[sep_dim]; });
			if (min_elt != it)
			{
				std::iter_swap(min_elt, it);
			}
			sep_val = (**it)[sep_dim];
			it++;
		}
		if (it == m_end)
		{
			KdPointPtrsIter max_elt = std::max_element(
			  m_begin, m_end,
			  [&](KdPointPtr lhs, KdPointPtr rhs)
			  { return lhs->primitive()[sep_dim] < rhs->primitive()[sep_dim]; });
			it--;
			if (max_elt != it)
			{
				std::iter_swap(max_elt, it);
			}
			sep_val = (**it)[sep_dim];
			it++;
		}
	}

	c_low.m_begin             = m_begin;
	c_low.m_end               = it;
	m_begin                   = it;
	// adjusting boxes
	m_bbox.min_coord(sep_dim) = sep_val;
	m_tbox.update_from_point_pointers(m_begin, m_end);
	c_low.m_bbox.max_coord(sep_dim) = sep_val;
	c_low.m_tbox.update_from_point_pointers(c_low.m_begin, c_low.m_end);
}

} // namespace OMC