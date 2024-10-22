#pragma once

#include "BoundingBox3T.h"

namespace OMC {

template <typename NT, typename VecT, typename PointT>
inline bool BoundingBox3T<NT, VecT, PointT>::operator==(const BboxT &b) const
{
	return m_min_bound == b.m_min_bound && m_max_bound == b.m_max_bound;
}

template <typename NT, typename VecT, typename PointT>
auto BoundingBox3T<NT, VecT, PointT>::operator+(const BboxT &b) const -> BboxT
{
	BboxT result = *this;
	result.m_min_bound.minimize(b.m_min_bound);
	result.m_max_bound.maximize(b.m_max_bound);
	return result;
}

template <typename NT, typename VecT, typename PointT>
auto BoundingBox3T<NT, VecT, PointT>::operator+=(const BboxT &b) -> BboxT &
{
	m_min_bound.minimize(b.m_min_bound);
	m_max_bound.maximize(b.m_max_bound);
	return *this;
}

template <typename NT, typename VecT, typename PointT>
auto BoundingBox3T<NT, VecT, PointT>::operator+(const PointT &p) const -> BboxT
{
	BboxT result = *this;
	result.m_min_bound.minimize(p);
	result.m_max_bound.maximize(p);
	return result;
}

template <typename NT, typename VecT, typename PointT>
auto BoundingBox3T<NT, VecT, PointT>::operator+=(const PointT &p) -> BboxT &
{
	m_min_bound.minimize(p);
	m_max_bound.maximize(p);
	return *this;
}

template <typename NT, typename VecT, typename PointT>
void BoundingBox3T<NT, VecT, PointT>::enlarge(const VecT &offset)
{
	m_min_bound -= offset;
	m_max_bound += offset;
}

template <typename NT, typename VecT, typename PointT>
size_t BoundingBox3T<NT, VecT, PointT>::longest_axis() const
{
	NT dx = m_max_bound.x() - m_min_bound.x();
	NT dy = m_max_bound.y() - m_min_bound.y();
	NT dz = m_max_bound.z() - m_min_bound.z();

	if (dx >= dy)
	{
		if (dx >= dz)
			return 0;
		else
			return 2; // dz>dx and dx>=dy
	}
	else // dy>dx
	{
		if (dy >= dz)
			return 1;
		else
			return 2; // dz>dy and dy>dx
	}
}

} // namespace OMC