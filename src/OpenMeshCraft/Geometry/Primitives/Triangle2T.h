#pragma once

#include "Primitive2.h"

namespace OMC {

/**
 * @brief Triangle in 2D that contains three points as its three vertices.
 * @tparam NT The number type.
 * @tparam PointT The point type.
 */
template <typename _NT, typename _PointT>
class Triangle2T : public Primitive2<_NT>
{
public:
	using NT        = _NT;
	using PointT    = _PointT;
	using TriangleT = Triangle2T<NT, PointT>;

public:
	/**
	 * @brief Construct triangle by default.
	 */
	Triangle2T() {}

	/**
	 * @brief Construct triangle from three points.
	 */
	Triangle2T(const PointT &v0, const PointT &v1, const PointT &v2)
	  : m_v0(v0)
	  , m_v1(v1)
	  , m_v2(v2)
	{
	}

	/// @name Access data of triangle.
	/// @{
	/**
	 * @brief Get the i-th vertex.
	 * @warning Is "pointer cast + offset" safe?
	 */
	PointT       &operator[](size_t i) { return (&m_v0)[i]; }
	const PointT &operator[](size_t i) const { return (&m_v0)[i]; }

	PointT       &v0() { return m_v0; }
	const PointT &v0() const { return m_v0; }
	PointT       &v1() { return m_v1; }
	const PointT &v1() const { return m_v1; }
	PointT       &v2() { return m_v2; }
	const PointT &v2() const { return m_v2; }
	/// @}

protected:
	PointT m_v0;
	PointT m_v1;
	PointT m_v2;
};

} // namespace OMC