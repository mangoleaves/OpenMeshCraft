#pragma once

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

#include <cmath>
#include <type_traits>

namespace OMC {

/**
 * @ingroup Geometry
 * @brief Vec4T, containing two numbers of type `T`, is an inexact type of
 * general Vector4T. We use `x`, `y`, `z` and `w` to refer to its four elements
 * in order.
 * @tparam T Number type.
 */
template <typename T>
class Vec4T
{
public:
	static constexpr size_t size_ = 4;
	using NT                      = T;
	using VT                      = Vec4T<T>;
	using value_type              = T;
	using vector_type             = VT;

public:
	/**
	 * @brief Default constructor.
	 */
	Vec4T() noexcept {}

	/**
	 * @brief Construct a new Vec4T object, vectorize it by \p v .
	 */
	explicit Vec4T(const T &v) noexcept { vectorize(v); }

	/**
	 * @brief Construct a new Vec4T object from four numbers.
	 */
	template <typename T1, typename T2, typename T3, typename T4,
	          typename =
	            std::enable_if_t<std::is_constructible_v<T, remove_cvref_t<T1>> &&
	                             std::is_constructible_v<T, remove_cvref_t<T2>> &&
	                             std::is_constructible_v<T, remove_cvref_t<T3>> &&
	                             std::is_constructible_v<T, remove_cvref_t<T4>>>>
	explicit Vec4T(T1 &&x, T2 &&y, T3 &&z, T4 &&w) noexcept
	  : _x(std::forward<T1>(x))
	  , _y(std::forward<T2>(y))
	  , _z(std::forward<T3>(z))
	  , _w(std::forward<T4>(w))
	{
	}

	/**
	 * @brief Construct a new Vec4T object from a pointer to numbers \p v .
	 * @note We assume the pointer points to at least four numbers and won't
	 * check the validity of the pointer.
	 * @param v the pointer to numbers.
	 */
	explicit Vec4T(const T *v) noexcept
	  : _x(v[0])
	  , _y(v[1])
	  , _z(v[2])
	  , _w(v[3])
	{
	}

	/**
	 * @brief Construct a new Vec4T object by copying from a given Vec4T \p v .
	 * @param v The given Vec4T.
	 */
	Vec4T(const VT &v) = default;

	/**
	 * @brief Construct a new Vec4T object by copying from a given Vec4T \p v .
	 * @param v The given Vec4T.
	 */
	Vec4T(VT &&v) = default;

	/**
	 * @brief Construct a new Vec4T object by copying from a given Vec4T \p v .
	 * @param v The given Vec4T.
	 * @return Vec4T& return the reference to this object.
	 */
	VT &operator=(const VT &v) = default;

	/**
	 * @brief Construct a new Vec4T object by copying from a given Vec4T \p v .
	 * @param v The given Vec4T.
	 * @return Vec4T& return the reference to this object.
	 */
	VT &operator=(VT &&v) = default;

	/**
	 * @brief Destroy the Vec4T object.
	 */
	~Vec4T() = default;

	/// @brief Get the reference to `x`
	inline T       &x() { return _x; }
	/// @brief Get the reference to `y`
	inline T       &y() { return _y; }
	/// @brief Get the reference to `z`
	inline T       &z() { return _z; }
	/// @brief Get the reference to `w`
	inline T       &w() { return _w; }
	/// @brief Get `x`
	inline const T &x() const { return _x; }
	/// @brief Get `y`
	inline const T &y() const { return _y; }
	/// @brief Get `z`
	inline const T &z() const { return _z; }
	/// @brief Get `w`
	inline const T &w() const { return _w; }

	/// @brief size(dimention) of vector
	size_t size() const { return size_; }

	/// @brief Get the pointer to the **mutable** elements.
	T       *data() { return (&_x); }
	/// @brief Get the pointer to the **constant** elements.
	const T *data() const { return (&_x); }

	/// @brief Access the element at dimension \p dim by reference.
	/// @warning Is "pointer cast + offset" safe?
	T &operator[](size_t dim)
	{
		OMC_EXPENSIVE_ASSERT(dim < size_, "dim {} out of range {}.", dim, size_);
		return (&_x)[dim];
	}
	/// @brief Access the element at dimension \p dim .
	/// @warning Is "pointer cast + offset" safe?
	const T &operator[](size_t dim) const
	{
		OMC_EXPENSIVE_ASSERT(dim < size_, "dim {} out of range {}.", dim, size_);
		return (&_x)[dim];
	}

	/// @brief Get the opposite of all elements.
	/// @todo If add more unary operators, try using a macro to repelace them.
	inline VT operator-() const { return VT{-_x, -_y, -_z, -_w}; }

	/// This micro defines two functions for a given binary operator.
	/// The functions apply element wisely binary operator to \p this and \p rhs .
#define ELEMENT_WISELY_BINARY_OPERATOR(op)                             \
	inline VT operator op(const VT &rhs) const                           \
	{                                                                    \
		return VT{_x op rhs._x, _y op rhs._y, _z op rhs._z, _w op rhs._w}; \
	}                                                                    \
	inline VT &&operator op(VT && rhs) const                             \
	{                                                                    \
		rhs._x = _x op rhs._x;                                             \
		rhs._y = _y op rhs._y;                                             \
		rhs._z = _z op rhs._z;                                             \
		rhs._w = _w op rhs._w;                                             \
		return std::move(rhs);                                             \
	}                                                                    \
	inline VT operator op(const T &rhs) const                            \
	{                                                                    \
		return VT{_x op rhs, _y op rhs, _z op rhs, _w op rhs};             \
	}

	ELEMENT_WISELY_BINARY_OPERATOR(-);
	ELEMENT_WISELY_BINARY_OPERATOR(+);
	ELEMENT_WISELY_BINARY_OPERATOR(*);
	ELEMENT_WISELY_BINARY_OPERATOR(/);

#undef ELEMENT_WISELY_BINARY_OPERATOR

	/// This micro defines two functions for a given compound assignment operator.
	/// The functions apply element wisely binary operator to \p this and \p rhs
	/// and assign the result to \p this .
#define ELEMENT_WISELY_COMPOUND_ASSIGNMENT(op) \
	inline VT &operator op(const VT & rhs)       \
	{                                            \
		_x op rhs._x;                              \
		_y op rhs._y;                              \
		_z op rhs._z;                              \
		_w op rhs._w;                              \
		return *this;                              \
	}                                            \
	inline VT &operator op(const T & rhs)        \
	{                                            \
		_x op rhs;                                 \
		_y op rhs;                                 \
		_z op rhs;                                 \
		_w op rhs;                                 \
		return *this;                              \
	}

	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(-=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(+=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(*=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(/=);

#undef ELEMENT_WISELY_COMPOUND_ASSIGNMENT

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)                                       \
	inline bool operator op(const VT &rhs) const                           \
	{                                                                      \
		return _x op rhs._x && _y op rhs._y && _z op rhs._z && _w op rhs._w; \
	}                                                                      \
	inline bool operator op(const T &rhs) const                            \
	{                                                                      \
		return _x op rhs && _y op rhs && _z op rhs && _w op rhs;             \
	}

	ELEMENT_WISELY_COMPARE(==);
	ELEMENT_WISELY_COMPARE(<);
	ELEMENT_WISELY_COMPARE(<=);
	ELEMENT_WISELY_COMPARE(>);
	ELEMENT_WISELY_COMPARE(>=);

#undef ELEMENT_WISELY_COMPARE

	inline bool operator!=(const VT &rhs) const { return !(*this == rhs); }
	inline bool operator!=(const T &rhs) const { return !(*this == rhs); }

	// TODO Add all() and any() for elemen-wisely operator.

	/// @brief Check if \p this is less than \p rhs only on dimension \p dim.
	inline bool less_on(size_t dim, const VT &rhs) const
	{
		return operator[](dim) < rhs[dim];
	}

	/// @brief Check if \p this is less than \p rhs only on dimension \p dim.
	inline bool less_on(size_t dim, const T &rhs) const
	{
		return operator[](dim) < rhs;
	}

	/// @brief Dot product of \p this and \p rhs .
	inline T dot(const VT &rhs) const
	{
		return _x * rhs._x + _y * rhs._y + _z * rhs._z + _w * rhs._w;
	}

	/// @brief Dot product of \p this and \p rhs, it is equivalent to \ref dot.
	inline T operator|(const VT &rhs) const { return this->dot(rhs); }

	/// @brief Get the sum of all elements.
	inline T sum() const { return _x + _y + _z + _w; }

	/// @brief Get the squared sum of all elements.
	inline T sqrnorm() const { return _x * _x + _y * _y + _z * _z + _w * _w; }

	/// @brief Minimize \p this by the other vector \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \min(x_i, y_i) \f).
	inline void minimize(const VT &rhs)
	{
		if (rhs._x < _x)
			_x = rhs._x;
		if (rhs._y < _y)
			_y = rhs._y;
		if (rhs._z < _z)
			_z = rhs._z;
		if (rhs._w < _w)
			_w = rhs._w;
	}

	/// @brief Maximize \p this by the other vector \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \max(x_i, y_i) \f).
	inline void maximize(const VT &rhs)
	{
		if (rhs._x > _x)
			_x = rhs._x;
		if (rhs._y > _y)
			_y = rhs._y;
		if (rhs._z > _z)
			_z = rhs._z;
		if (rhs._w > _w)
			_w = rhs._w;
	}

	/// @brief Vectorize \p this by a given number \p rhs ,
	/// i.e., set all elements of \p this to \p rhs.
	inline VT &vectorize(const T &s)
	{
		_x = s;
		_y = s;
		_z = s;
		_w = s;
		return *this;
	}

	/// @brief Normalize \p this inplace and return the reference to \p this .
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline VT &normalize()
	{
		T sn = sqrnorm();
		return sn > 0 ? (*this) /= std::sqrt(sn) : (*this);
	}

	/// @brief Get a copy of normailized \p this .
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline VT normalized() const
	{
		T sn = sqrnorm();
		return sn > 0 ? (*this) / std::sqrt(sn) : (*this);
	}

	/// @brief Get the \f( \ell_2 \f)-norm of \p this .
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline T norm() const
	{
		return std::sqrt(sqrnorm());
	}

	/// @brief Get the \f( \ell_2 \f)-norm of \p this ,
	/// which is equivalent to \ref norm.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline T length() const
	{
		return norm();
	}

	/// @brief Get the element-wisely pow result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline VT pow(const T &exp) const
	{
		return VT{std::pow(_x, exp), std::pow(_y, exp), std::pow(_z, exp),
		          std::pow(_w, exp)};
	}

	/// @brief Get the \f( \ell_p \f)-norm of \p this .
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline T pnorm(const T &exp) const
	{
		T px = std::pow(std::abs(_x), exp);
		T py = std::pow(std::abs(_y), exp);
		T pz = std::pow(std::abs(_z), exp);
		T pw = std::pow(std::abs(_w), exp);
		return std::pow(px + py + pz + pw, 1. / exp);
	}

	/// @brief Check if all elements are finite.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline bool isfinite() const
	{
		return std::isfinite(_x) && std::isfinite(_y) && std::isfinite(_z) &&
		       std::isfinite(_w);
	}

	/// @brief Get element-wisely floor result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT floor() const
	{
		return VT{std::floor(_x), std::floor(_y), std::floor(_z), std::floor(_w)};
	}

	/// @brief Get element-wisely ceil result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT ceil() const
	{
		return VT{std::ceil(_x), std::ceil(_y), std::ceil(_z), std::ceil(_w)};
	}

	/// @brief Get element-wisely round result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT round() const
	{
		return VT{std::round(_x), std::round(_y), std::round(_z), std::round(_w)};
	}

	/// @brief Get element-wisely absolute result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT abs() const
	{
		return VT{std::abs(_x), std::abs(_y), std::abs(_z), std::abs(_w)};
	}

	/// @brief Negate elements in-place.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT &negate()
	{
		_x = -_x;
		_y = -_y;
		_z = -_z;
		_w = -_w;
		return *this;
	}

protected:
	T _x;
	T _y;
	T _z;
	T _w;
}; // namespace Geometry

/******* Below operators are optimized for right values. *******/

template <typename T>
inline Vec4T<T> operator-(Vec4T<T> &&lhs)
{
	lhs.negate();
	return std::move(lhs);
}

/// This micro defines two functions for a given binary operator.
/// The functions apply element wisely binary operator to \p this and \p rhs .
/// They are specially optimized for right values.
#define ELEMENT_WISELY_BINARY_OPERATOR_RV(op, compound_op)           \
	template <typename T>                                              \
	inline Vec4T<T> &&operator op(Vec4T<T> &&lhs, const Vec4T<T> &rhs) \
	{                                                                  \
		return std::move(lhs compound_op rhs);                           \
	}                                                                  \
	template <typename T>                                              \
	inline Vec4T<T> &&operator op(Vec4T<T> &&lhs, Vec4T<T> &&rhs)      \
	{                                                                  \
		return std::move(lhs compound_op rhs);                           \
	}                                                                  \
	template <typename T>                                              \
	inline Vec4T<T> &&operator op(Vec4T<T> &&lhs, const T & rhs)       \
	{                                                                  \
		return std::move(lhs compound_op rhs);                           \
	}

ELEMENT_WISELY_BINARY_OPERATOR_RV(-, -=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(+, +=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(*, *=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(/, /=);

#undef ELEMENT_WISELY_BINARY_OPERATOR_RV

/// @brief Multiply a number \p lhs and a Vec4T \p rhs .
template <typename T>
inline Vec4T<T> operator*(const T &lhs, const Vec4T<T> &rhs)
{
	return rhs * lhs;
}

/// @brief Multiply a number \p lhs and a Vec4T \p rhs .
template <typename T>
inline Vec4T<T> &&operator*(const T &lhs, Vec4T<T> &&rhs)
{
	return std::move(std::move(rhs) * lhs);
}

/// @brief A number \p lhs divede a Vec4T \p rhs element-wisely.
/// result[i] = lhs / rhs[i].
template <typename T>
inline Vec4T<T> operator/(const T &lhs, const Vec4T<T> &rhs)
{
	return Vec4T<T>{lhs / rhs.x(), lhs / rhs.y(), lhs / rhs.z(), lhs / rhs.w()};
}

/// @brief A number \p lhs divede a Vec4T \p rhs element-wisely.
/// result[i] = lhs / rhs[i].
template <typename T>
inline Vec4T<T> &&operator/(const T &lhs, Vec4T<T> &&rhs)
{
	rhs.x() = lhs / rhs.x();
	rhs.y() = lhs / rhs.y();
	rhs.z() = lhs / rhs.z();
	rhs.w() = lhs / rhs.w();
	return std::move(rhs);
}

using Vec4f = Vec4T<float>;
using Vec4d = Vec4T<double>;
using Vec4i = Vec4T<int>;
using Vec4u = Vec4T<unsigned int>;

} // namespace OMC