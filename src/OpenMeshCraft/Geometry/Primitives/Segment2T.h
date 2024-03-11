#pragma once

#include "Primitive2.h"

namespace OMC {

/**
 * @brief Segment in 2D that contains two points as start and end.
 * @tparam NT The number type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _PointT>
class Segment2T : public Primitive2<_NT>
{
public:
	using NT       = _NT;
	using PointT   = _PointT;
	using SegmentT = Segment2T<NT, PointT>;

public:
	/**
	 * @brief Construct Segment3T by default.
	 */
	Segment2T() {}

	/**
	 * @brief Construct Segment3T from given \p start and \p end.
	 * @param start start point of segment
	 * @param end end point of segment
	 */
	Segment2T(const PointT &start, const PointT &end)
	  : m_start(start)
	  , m_end(end)
	{
	}

	/// @name Acess the data of Plane
	/// @{
	PointT       &start() { return m_start; }
	const PointT &start() const { return m_start; }
	PointT       &end() { return m_end; }
	const PointT &end() const { return m_end; }
	/// @}

private:
	PointT m_start;
	PointT m_end;
};

} // namespace OMC