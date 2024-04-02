#pragma once

#include "BinaryTree.h"

namespace OMC {

/// @todo Add an option to find first or all intersections.
template <typename Traits>
class Binary_BoxInterTraversal
{
public:
	using NT        = typename Traits::NT;
	using BboxT     = typename Traits::BboxT;
	using TreeBboxT = typename Traits::TreeBboxT;

	using CalcBbox    = typename Traits::CalcBbox;
	using DoIntersect = typename Traits::DoIntersect;

public:
	template <typename QPrimT>
	Binary_BoxInterTraversal(const QPrimT &query)
	{
		m_box_of_query = CalcBbox()(query);
	}

	template <>
	Binary_BoxInterTraversal(const BboxT &query)
	{
		m_box_of_query = query;
	}

	bool intersection(const TreeBboxT &bbox)
	{
		if (DoIntersect()(m_box_of_query, bbox.bbox()))
			m_intersected_ids.push_back(bbox.id());
		return true;
	}

	bool do_inter(const BboxT &bbox) const
	{
		return DoIntersect()(bbox, m_box_of_query);
	}

	const std::vector<index_t> &result() const { return m_intersected_ids; }

private:
	BboxT                m_box_of_query;
	std::vector<index_t> m_intersected_ids;
};

} // namespace OMC