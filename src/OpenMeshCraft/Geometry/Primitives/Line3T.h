#pragma once

#include "Primitive3.h"

namespace OMC {

/**
 * @brief Line in 3D that defined by various ways.
 * @tparam NT The number type.
 * @tparam VecT The vector type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _VecT, typename _PointT>
class Line3T : public Primitive3<_NT>
{
public:
	using NT     = _NT;
	using VecT   = _VecT;
	using PointT = _PointT;

public:
	/**
	 * @brief Construct a Line3T object by default.
	 */
	Line3T() {}

	/**
	 * @brief Construct a Line3T from given point \p start and direction \p dir.
	 * @param start start point of the line
	 * @param dir direction of the line
	 */
	explicit Line3T(const PointT &start, const VecT &dir) noexcept
	{
		m_start     = start;
		m_direction = dir;
	}

	/**
	 * @brief Construct a Line3T from given two point \p start and \p end.
	 * Will calculate direction (end-start) and save in object.
	 * @param start start point of the line
	 * @param end end point of the line
	 * @note Possibly lose precision in calculating direction.
	 */
	explicit Line3T(const PointT &start, const PointT &end) noexcept
	{
		m_start     = start;
		m_direction = end - start;
	}

	/// @name Acess the data of Line
	/// @{
	const PointT &start() const { return m_start; }
	PointT       &start() { return m_start; }

	const VecT &direction() const { return m_direction; }
	VecT       &direction() { return m_direction; }
	/// @}

private:
	/// @brief The start point
	PointT m_start;
	/// @brief The direction vector
	VecT   m_direction;
};

} // namespace OMC