#pragma once

#include "Primitive2.h"

namespace OMC {

/**
 * @brief Sphere in 2D that contains a point as center and and a number as
 * squared radius.
 * @tparam NT The number type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _PointT>
class Sphere2T : public Primitive2<_NT>
{
public:
	using NT      = _NT;
	using PointT  = _PointT;
	using SphereT = Sphere2T<NT, PointT>;

public:
	/**
	 * @brief Construct Sphere2T by default.
	 */
	Sphere2T() {}

	/**
	 * @brief Construct Sphere2T from a point \p center and a number \p
	 * squared_radius .
	 * @param center the center of sphere.
	 * @param squared_radius the squared radius of sphere.
	 */
	Sphere2T(const PointT &center, const NT &squared_radius)
	  : m_center(center)
	  , m_squared_radius(squared_radius)
	{
	}

	/// @name Acess the data of Sphere
	/// @{
	PointT       &center() { return m_center; }
	const PointT &center() const { return m_center; }
	NT           &squared_radius() { return m_squared_radius; }
	const NT     &squared_radius() const { return m_squared_radius; }
	/// @}

private:
	PointT m_center;
	NT     m_squared_radius;
};

} // namespace OMC