#pragma once

#include "Primitive3.h"

namespace OMC {

enum class Ray3TDefinedBy
{
	START_DIR, // start point + direction
	TWO_PNTS   // two end points
};

/**
 * @brief Line in 3D that defined by start point and emit direction.
 * @tparam NT The number type.
 * @tparam VecT The vector type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _VecT, typename _PointT>
class Ray3T : public Primitive3<_NT>
{
public:
	using NT     = _NT;
	using VecT   = _VecT;
	using PointT = _PointT;

public:
	/**
	 * @brief Construct a Ray3T object by default.
	 */
	Ray3T() {}

	/**
	 * @brief Construct a Ray3T from given point \p start and direction \p dir.
	 * @param start start point of the line
	 * @param dir direction of the line
	 */
	explicit Ray3T(const PointT &start, const VecT &dir) noexcept
	{
		m_start     = start;
		m_direction = dir;
	}

	/**
	 * @brief Construct a Ray3T from given two point \p start and \p end.
	 * Calculate the direction (end-start).
	 * @param start start point of the line
	 * @param end end point of the line
	 * @note Possibly lose precision in calculating direction.
	 */
	explicit Ray3T(const PointT &start, const PointT &end) noexcept
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

	void inversed() { m_direction = -m_direction; }

	Ray3T inverse() const
	{
		Ray3T new_vec;
		new_vec.inversed();
		return new_vec;
	}

private:
	/// @brief The start point
	PointT m_start;
	/// @brief The emit direction
	VecT   m_direction;
};

} // namespace OMC