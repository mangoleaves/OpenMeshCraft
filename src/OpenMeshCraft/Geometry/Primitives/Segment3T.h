#pragma once

#include "Primitive3.h"

namespace OMC {

/**
 * @brief Segment in 3D that contains two points as start and end.
 * @tparam NT The number type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _PointT>
class Segment3T : public Primitive3<_NT>
{
public:
	using NT       = _NT;
	using PointT   = _PointT;
	using SegmentT = Segment3T<NT, PointT>;

public:
	/**
	 * @brief Construct Segment3T by default.
	 */
	Segment3T() {}

	/**
	 * @brief Construct Segment3T from given \p start and \p end.
	 * @param start start point of segment
	 * @param end end point of segment
	 */
	Segment3T(const PointT &start, const PointT &end)
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