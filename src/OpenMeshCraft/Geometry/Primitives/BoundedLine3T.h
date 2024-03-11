#pragma once

#include "Primitive3.h"

namespace OMC {

/**
 * @brief Bounded Line in 3D.
 * @tparam NT The number type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _PointT>
class BoundedLine3T : public Primitive3<_NT>
{
public:
	using NT     = _NT;
	using PointT = _PointT;

public:
	/**
	 * @brief Construct a bounded Line3T object by default.
	 */
	BoundedLine3T() {}

	/**
	 * @brief Construct a bounded Line3T from given two point \p start and \p end.
	 * Specify whether the line is bounded at the two end points by input
	 * parameters.
	 * @param start start point of the line
	 * @param start_bounded true to bound the line at start point
	 * @param end end point of the line
	 * @param end_bounded true to bound the line at end point
	 */
	explicit BoundedLine3T(const PointT &start, bool start_bounded,
	                       const PointT &end, bool end_bounded) noexcept
	{
		m_start = start;
		m_end   = end;
	}

	/// @name Acess the data of Line
	/// @{
	const PointT &start() const { return m_start; }
	PointT       &start() { return m_start; }

	const PointT &end() const { return m_end; }
	PointT       &end() { return m_end; }

	bool start_bounded() const { return m_start_bounded; }
	bool end_bounded() const { return m_end_bounded; }

	/// @brief bound = true to set the line is bounded at start point
	void start_bounded(bool bound) { m_start_bounded = bound; }
	/// @brief bound = true to set the line is bounded at end point
	void end_bounded(bool bound) { m_end_bounded = bound; }
	/// @}

private:
	/// @brief The first point, aka start point.
	PointT m_start;
	/// @brief The second point, aka end point.
	PointT m_end;
	/// @brief flag to indicate whether the line is bounded at start
	bool   m_start_bounded;
	/// @brief flag to indicate whether the line is bounded at end
	bool   m_end_bounded;
};

} // namespace OMC