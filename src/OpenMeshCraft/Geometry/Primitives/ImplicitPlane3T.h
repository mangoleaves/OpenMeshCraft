#pragma once

#include "Primitive3.h"

namespace OMC {

/**
 * @brief Implicit plane in 3D that represented as implicit function
 * ax+by+cz+d=0.
 * @tparam NT The number type.
 * @tparam VecT The vector type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _VecT, typename _PointT>
class ImplicitPlane3T : public Primitive3<_NT>
{
public:
	using NT             = _NT;
	using VecT           = _VecT;
	using PointT         = _PointT;
	using ImplicitPlaneT = ImplicitPlane3T<NT, VecT, PointT>;

public:
	/**
	 * @brief Construct a ImplicitPlane3T object by default.
	 */
	ImplicitPlane3T() {}

	/**
	 * @brief Construct a ImplicitPlane3T object from given \p point and \p
	 * normal.
	 * @param point a point on the plane
	 * @param normal a normal direction of the plane
	 */
	ImplicitPlane3T(const PointT &point, const VecT &normal)
	{
		from_explicit(point, normal);
	}

	/**
	 * @brief Construct a ImplicitPlane3T object from given \p point and \p
	 * normal.
	 * @param point a point on the plane
	 * @param normal a normal direction of the plane
	 */
	void from_explicit(const PointT &point, const VecT &normal)
	{
		m_a = normal.x();
		m_b = normal.y();
		m_c = normal.z();
		m_d = -point.as_vec().dot(normal);
	}

	VecT normal() { return VecT(m_a, m_b, m_c); }

	const NT &a() const { return m_a; }
	const NT &b() const { return m_b; }
	const NT &c() const { return m_c; }
	const NT &d() const { return m_d; }

	NT &a() { return m_a; }
	NT &b() { return m_b; }
	NT &c() { return m_c; }
	NT &d() { return m_d; }

private:
	/// The plane is represented as ax+by+cz+d=0.
	NT m_a, m_b, m_c, m_d;
};

} // namespace OMC