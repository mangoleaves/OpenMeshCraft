#pragma once

#include "Primitive3.h"
#include "Vector3T.h"

namespace OMC {

/**
 * @brief Point in 3D that contains one Vec3 as its coordinate.
 * @tparam NT Number Type.
 */
template <typename _NT>
class Point3T : public Primitive3<_NT>
{
public:
	using NT   = _NT;
	using VecT = Vec3T<NT>;

public:
	Point3T() noexcept
	  : m_p()
	{
	}

	Point3T(const Point3T &)             = default;
	Point3T(Point3T &&)                  = default;
	Point3T &operator=(const Point3T &v) = default;
	Point3T &operator=(Point3T &&v)      = default;

	/**
	 * @brief Construct a new Point3T object, vectorize it by \p v .
	 */
	explicit Point3T(const NT &v) noexcept
	  : m_p(v)
	{
	}

	/**
	 * @brief Construct a new Point3T object from a pointer to number \p v .
	 * @note We assume the pointer points to at least three numbers and won't
	 * check the validity of the pointer.
	 * @param v the pointer to number.
	 */
	explicit Point3T(const NT *v) noexcept
	  : m_p(v)
	{
	}

	/**
	 * @brief Construct a new Point3T object from VecT.
	 */
	template <typename TT, typename = std::enable_if_t<
	                         std::is_same_v<remove_cvref_t<TT>, VecT>> /*SIFNAE*/>
	Point3T(TT &&v) noexcept
	  : m_p(std::forward<TT>(v))
	{
	}

	/**
	 * @brief Construct a new Point3T object from three numbers.
	 */
	template <
	  typename T1, typename T2, typename T3,
	  typename = std::enable_if_t<std::is_same_v<remove_cvref_t<T1>, NT> &&
	                              std::is_same_v<remove_cvref_t<T2>, NT> &&
	                              std::is_same_v<remove_cvref_t<T3>, NT>>>
	explicit Point3T(T1 &&x, T2 &&y, T3 &&z) noexcept
	  : m_p(std::forward<T1>(x), std::forward<T2>(y), std::forward<T3>(z))
	{
	}

	/// @name Data access
	/// @{
	NT       &x() { return m_p.x(); }
	NT       &y() { return m_p.y(); }
	NT       &z() { return m_p.z(); }
	const NT &x() const { return m_p.x(); }
	const NT &y() const { return m_p.y(); }
	const NT &z() const { return m_p.z(); }

	NT       &operator[](size_t i) { return m_p[i]; }
	const NT &operator[](size_t i) const { return m_p[i]; }

	/// @brief By as_vec, you can get a reference and operate it as VecT.
	/// Any modification to the corresponding VecT is also applied to this
	/// Point3T.
	VecT &as_vec() { return m_p; }

	/// @brief By as_vec, you can get a ference and operate it as VecT.
	const VecT &as_vec() const { return m_p; }

	/// @brief By to_vec, you can get a copy of this Point3T,
	/// whose type is now VecT.
	VecT to_vec() { return m_p; }
	/// @}

	/**
	 * @brief negative
	 * @return Point3T the negative point.
	 */
	Point3T operator-() const { return Point3T(-as_vec()); }

	/**
	 * @brief vector = this point - the other point
	 * @param rhs the other point
	 * @return VecT the result vector between two points
	 */
	VecT operator-(const Point3T &rhs) const { return as_vec() - rhs.as_vec(); }

	/**
	 * @brief new point = this point + another point
	 * @param rhs the given point
	 * @return Point3T the result point.
	 */
	Point3T operator+(const Point3T &rhs) const
	{
		return Point3T(as_vec() + rhs.as_vec());
	}

	/**
	 * @brief new point = this point + given vector
	 * @param rhs the given vector
	 * @return Point3T the result point by adding a vector to this point.
	 */
	Point3T operator+(const VecT &rhs) const { return Point3T(as_vec() + rhs); }

	/**
	 * @brief new point = this point - given vector
	 * @param rhs the given vector
	 * @return Point3T the result point by adding a vector to this point.
	 */
	Point3T operator-(const VecT &rhs) const { return Point3T(as_vec() - rhs); }

	/**
	 * @brief new point = point * scale_factor
	 * @param rhs the scale factor
	 * @return Point3T the scaled new point
	 */
	Point3T operator*(const NT &scale_factor) const
	{
		return Point3T(as_vec() * scale_factor);
	}

	/**
	 * @brief new point = point / scale_factor
	 * @param rhs the scale factor
	 * @return Point3T the scaled new point
	 */
	Point3T operator/(const NT &scale_factor) const
	{
		return Point3T(as_vec() / scale_factor);
	}

	/**
	 * @brief this point += another point
	 * @param rhs the given point
	 * @return Point3T the result point
	 */
	Point3T &operator+=(const Point3T &rhs)
	{
		as_vec() += rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point += given vector
	 * @param rhs the given vector
	 * @return Point3T the result point by adding a vector to this point.
	 */
	Point3T &operator+=(const VecT &rhs)
	{
		as_vec() += rhs;
		return *this;
	}

	/**
	 * @brief this point -= another point
	 * @param rhs the given point
	 * @return Point3T the result point.
	 */
	Point3T &operator-=(const Point3T &rhs)
	{
		as_vec() -= rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point -= given vector
	 * @param rhs the given vector
	 * @return Point3T the result point by adding a vector to this point.
	 */
	Point3T &operator-=(const VecT &rhs)
	{
		as_vec() -= rhs;
		return *this;
	}

	/**
	 * @brief point *= scale_factor
	 * @param rhs the scale factor
	 * @return Point3T the in-place scaled point
	 */
	Point3T operator*=(const NT &scale_factor)
	{
		as_vec() *= scale_factor;
		return *this;
	}

	/**
	 * @brief point /= scale_factor
	 * @param rhs the scale factor
	 * @return Point3T the in-place scaled point
	 */
	Point3T operator/=(const NT &scale_factor)
	{
		as_vec() /= scale_factor;
		return *this;
	}

	/// @brief Minimize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \min(x_i, y_i) \f).
	void minimize(const Point3T &rhs) { m_p.minimize(rhs.m_p); }

	/// @brief Maximize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \max(x_i, y_i) \f).
	void maximize(const Point3T &rhs) { m_p.maximize(rhs.m_p); }

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)                                             \
	inline bool operator op(const Point3T &rhs) const { return m_p op rhs.m_p; } \
	inline bool operator op(const NT &rhs) const { return m_p op rhs; }

	ELEMENT_WISELY_COMPARE(==);
	ELEMENT_WISELY_COMPARE(!=);
	ELEMENT_WISELY_COMPARE(<);
	ELEMENT_WISELY_COMPARE(<=);
	ELEMENT_WISELY_COMPARE(>);
	ELEMENT_WISELY_COMPARE(>=);
#undef ELEMENT_WISELY_COMPARE

private:
	VecT m_p;
};

/**
 * @brief new point = point * scale_factor
 * @param rhs the scale factor
 * @return Point3T the scaled new point
 */
template <typename NT>
inline Point3T<NT> operator*(const NT &scale_factor, const Point3T<NT> &point)
{
	return point * scale_factor;
}

} // namespace OMC