#pragma once

#include "Primitive2.h"
#include "Vector2T.h"

namespace OMC {

/**
 * @brief Point in 2D that contains one Vec2 as its coordinate.
 * @tparam NT Number Type.
 */
template <typename _NT>
class Point2T : public Primitive2<_NT>
{
public:
	using NT   = _NT;
	using VecT = Vec2T<NT>;

public:
	Point2T() noexcept
	  : m_p()
	{
	}

	/**
	 * @brief Construct a new Point2T object, vectorize it by \p v .
	 * @tparam
	 */
	explicit Point2T(const NT &v) noexcept
	  : m_p(v)
	{
	}

	/**
	 * @brief Construct a new Point2T object from a pointer to number \p v .
	 * @note We assume the pointer points to at least two numbers and won't
	 * check the validity of the pointer.
	 * @param v the pointer to number.
	 */
	explicit Point2T(const NT *v) noexcept
	  : m_p(v)
	{
	}

	/**
	 * @brief Construct a new Point2T object from VecT.
	 * @tparam
	 */
	template <typename TT, typename = std::enable_if_t<
	                         std::is_same_v<remove_cvref_t<TT>, VecT>> /*SIFNAE*/>
	explicit Point2T(TT &&v) noexcept
	  : m_p(std::forward<TT>(v))
	{
	}

	/**
	 * @brief Construct a new Point2T object from two numbers.
	 */
	template <
	  typename T1, typename T2,
	  typename = std::enable_if_t<std::is_same_v<remove_cvref_t<T1>, NT> &&
	                              std::is_same_v<remove_cvref_t<T2>, NT>>>
	explicit Point2T(T1 &&x, T2 &&y) noexcept
	  : m_p(std::forward<T1>(x), std::forward<T2>(y))
	{
	}

	/**
	 * @brief Construct a new Point2T by copying from a given Point2T \p v.
	 * @param v The given Point2T.
	 */
	Point2T(const Point2T &v) noexcept
	  : m_p(v.m_p)
	{
	}

	/**
	 * @brief Construct a new Point2T by copying from a given Point2T \p v.
	 * @param v The given Point2T.
	 */
	Point2T(Point2T &&v) noexcept
	  : m_p(std::move(v.m_p))
	{
	}

	/**
	 * @brief Construct a new Point2T by copying from a given Point2T \p v.
	 * @param v The given Point2T.
	 * @return Point2T& return the reference to this object.
	 */
	Point2T &operator=(const Point2T &v) noexcept
	{
		m_p = v.m_p;
		return *this;
	}

	/**
	 * @brief Construct a new Point2T by copying from a given Point2T \p v.
	 * @param v The given Point2T.
	 * @return Point2T& return the reference to this object.
	 */
	Point2T &operator=(Point2T &&v) noexcept
	{
		m_p = std::move(v.m_p);
		return *this;
	}

	/// @name Data access
	/// @{
	NT       &x() { return m_p.x(); }
	NT       &y() { return m_p.y(); }
	const NT &x() const { return m_p.x(); }
	const NT &y() const { return m_p.y(); }

	NT       &operator[](size_t i) { return m_p[i]; }
	const NT &operator[](size_t i) const { return m_p[i]; }

	/// @brief By as_vec, you can get a reference and operate it as VecT.
	/// Any modification to the corresponding VecT is also applied to this
	/// Point2T.
	VecT &as_vec() { return m_p; }

	/// @brief By as_vec, you can get a ference and operate it as VecT.
	const VecT &as_vec() const { return m_p; }

	/// @brief By to_vec, you can get a copy of this Point2T,
	/// whose type is now VecT.
	VecT to_vec() { return m_p; }
	/// @}

	/**
	 * @brief negative
	 * @return Point2T the negative point.
	 */
	Point2T operator-() const { return Point2T(-as_vec()); }

	/**
	 * @brief vector = this point - the other point
	 * @param rhs the other point
	 * @return VecT the result vector between two points
	 */
	VecT operator-(const Point2T &rhs) const { return as_vec() - rhs.as_vec(); }

	/**
	 * @brief new point = this point + another point
	 * @param rhs the given point
	 * @return Point2T the result point.
	 */
	Point2T operator+(const Point2T &rhs) const
	{
		return Point2T(as_vec() + rhs.as_vec());
	}

	/**
	 * @brief new point = this point + given vector
	 * @param rhs the given vector
	 * @return Point2T the result point by adding a vector to this point.
	 */
	Point2T operator+(const VecT &rhs) const { return Point2T(as_vec() + rhs); }

	/**
	 * @brief new point = this point - given vector
	 * @param rhs the given vector
	 * @return Point2T the result point by adding a vector to this point.
	 */
	Point2T operator-(const VecT &rhs) const { return Point2T(as_vec() - rhs); }

	/**
	 * @brief new point = point * scale_factor
	 * @param rhs the scale factor
	 * @return Point2T the scaled new point
	 */
	Point2T operator*(const NT &scale_factor) const
	{
		return Point2T(as_vec() * scale_factor);
	}

	/**
	 * @brief this point += another point
	 * @param rhs the given point
	 * @return Point2T the result point.
	 */
	Point2T &operator+=(const Point2T &rhs)
	{
		as_vec() += rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point += given vector
	 * @param rhs the given vector
	 * @return Point2T the result point by adding a vector to this point.
	 */
	Point2T &operator+=(const VecT &rhs)
	{
		as_vec() += rhs;
		return *this;
	}

	/**
	 * @brief this point -= another point
	 * @param rhs the given point
	 * @return Point2T the result point.
	 */
	Point2T &operator-=(const Point2T &rhs)
	{
		as_vec() -= rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point -= given vector
	 * @param rhs the given vector
	 * @return Point2T the result point by adding a vector to this point.
	 */
	Point2T &operator-=(const VecT &rhs)
	{
		as_vec() -= rhs;
		return *this;
	}

	/**
	 * @brief point *= scale_factor
	 * @param rhs the scale factor
	 * @return Point2T the in-place scaled point
	 */
	Point2T operator*=(const NT &scale_factor)
	{
		as_vec() *= scale_factor;
		return *this;
	}

	/// @brief Minimize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \min(x_i, y_i) \f).
	void minimize(const Point2T &rhs) { m_p.minimize(rhs.m_p); }

	/// @brief Maximize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \max(x_i, y_i) \f).
	void maximize(const Point2T &rhs) { m_p.maximize(rhs.m_p); }

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)                                             \
	inline bool operator op(const Point2T &rhs) const { return m_p op rhs.m_p; } \
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
 * @return Point2T the scaled new point
 */
template <typename NT>
inline Point2T<NT> operator*(const NT &scale_factor, const Point2T<NT> &point)
{
	return point * scale_factor;
}

} // namespace OMC