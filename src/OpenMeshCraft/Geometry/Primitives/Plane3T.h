#pragma once

#include "Primitive3.h"

namespace OMC {

/**
 * @brief Plane in 3D that contains a point on the plane and a normal direction
 * of the plane.
 * @tparam NT The number type.
 * @tparam VecT The vector type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _VecT, typename _PointT>
class Plane3T : public Primitive3<_NT>
{
public:
	using NT     = _NT;
	using VecT   = _VecT;
	using PointT = _PointT;
	using PlaneT = Plane3T<NT, VecT, PointT>;

public:
	/**
	 * @brief Construct a Plane3T object by default.
	 */
	Plane3T() {}

	/**
	 * @brief Construct a Plane3T object from given \p point and \p normal.
	 * @param point a point on the plane
	 * @param normal a normal direction of the plane
	 */
	Plane3T(const PointT &point, const VecT &normal)
	  : m_point(point)
	  , m_normal(normal)
	{
	}

	/// @name Acess the data of Plane
	/// @{
	PointT       &point() { return m_point; }
	const PointT &point() const { return m_point; }
	VecT         &normal() { return m_normal; }
	const VecT   &normal() const { return m_normal; }
	/// @}

private:
	/// @brief The point on the plane.
	PointT m_point;
	/// @brief The normal direction of the plane. We do not require it is unit.
	VecT   m_normal;
};

} // namespace OMC