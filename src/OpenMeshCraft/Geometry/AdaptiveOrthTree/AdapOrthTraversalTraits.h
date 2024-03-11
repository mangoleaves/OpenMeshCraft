#pragma once

#include "AdapOrthTree.h"

#include "OpenMeshCraft/Utils/SFINAE.h"

namespace OMC {

/// @todo Add an option to find first or all intersections.
template <typename Traits>
class AdapOrth_BoxInterTraversal
{
public:
	using NT      = typename Traits::NT;
	using BboxT   = typename Traits::BboxT;
	using OrBboxT = typename Traits::OrBboxT;

	using CalcBbox    = typename Traits::CalcBbox;
	using DoIntersect = typename Traits::DoIntersect;

public:
	template <typename QPrimT>
	AdapOrth_BoxInterTraversal(const QPrimT &query)
	{
		m_box_of_query = CalcBbox()(query);
	}

	bool intersection(const OrBboxT &leaf_bbox)
	{
		if (DoIntersect()(m_box_of_query, leaf_bbox.bbox()))
			m_intersected_ids.push_back(leaf_bbox.id());
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