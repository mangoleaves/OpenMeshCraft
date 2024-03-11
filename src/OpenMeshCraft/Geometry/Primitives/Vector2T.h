#pragma once

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

#include <cmath>
#include <type_traits>

namespace OMC {

/**
 * @ingroup Geometry
 * @brief Vec2d, containing two numbers of type `T`, is an inexact type of
 * general Vector2T. We use `x` to refer to its first element and `y` to refer
 * to its second element.
 * @tparam T Number type.
 */
template <typename T>
class Vec2T
{
public:
	static constexpr size_t size_ = 2;
	using NT                      = T;
	using VT                      = Vec2T<T>;
	using value_type              = T;
	using vector_type             = VT;

public:
	/**
	 * @brief Default constructor.
	 */
	Vec2T() noexcept {}

	/**
	 * @brief Construct a new Vec2T object, vectorize it by \p v .
	 */
	explicit Vec2T(const T &v) noexcept { vectorize(v); }

	/**
	 * @brief Construct a new Vec2T object from two numbers \p x and \p y .
	 * @tparam TT We want to benefit from the std::forward.
	 */
	template <typename T1, typename T2,
	          typename = std::enable_if_t<std::is_constructible_v<T, T1> &&
	                                      std::is_constructible_v<T, T2>>>
	explicit Vec2T(T1 &&x, T2 &&y) noexcept
	  : _x(std::forward<T1>(x))
	  , _y(std::forward<T2>(y))
	{
	}

	/**
	 * @brief Construct a new Vec2T object from a pointer to numbers \p v .
	 * @note We assume the pointer points to at least two numbers and won't
	 * check the validity of the pointer.
	 * @param v the pointer to numbers.
	 */
	explicit Vec2T(const T *v) noexcept
	  : _x(v[0])
	  , _y(v[1])
	{
	}

	/**
	 * @brief Construct a new Vec2T object by copying from a given Vec2T \p v .
	 * @param v The given Vec2T.
	 */
	Vec2T(const VT &v) noexcept
	  : _x(v._x)
	  , _y(v._y)
	{
	}

	/**
	 * @brief Construct a new Vec2T object by copying from a given Vec2T \p v .
	 * @param v The given Vec2T.
	 */
	Vec2T(VT &&v) noexcept
	  : _x(std::move(v._x))
	  , _y(std::move(v._y))
	{
	}

	/**
	 * @brief Construct a new Vec2T object by copying from a given Vec2T \p v .
	 * @param v The given Vec2T.
	 * @return Vec2T& return the reference to this object.
	 */
	VT &operator=(const VT &v) noexcept
	{
		_x = v._x;
		_y = v._y;
		return *this;
	}

	/**
	 * @brief Construct a new Vec2T object by copying from a given Vec2T \p v .
	 * @param v The given Vec2T.
	 * @return Vec2T& return the reference to this object.
	 */
	VT &operator=(VT &&v) noexcept
	{
		_x = std::move(v._x);
		_y = std::move(v._y);
		return *this;
	}

	/**
	 * @brief Destroy the Vec2T object.
	 */
	~Vec2T() noexcept {}

	/// @brief Get the reference to `x`
	inline T       &x() { return _x; }
	/// @brief Get the reference to `y`
	inline T       &y() { return _y; }
	/// @brief Get `x`
	inline const T &x() const { return _x; }
	/// @brief Get `y`
	inline const T &y() const { return _y; }

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
	inline VT operator-() const { return VT{-_x, -_y}; }

	/// This micro defines two functions for a given binary operator.
	/// The functions apply element wisely binary operator to \p this and \p rhs .
#define ELEMENT_WISELY_BINARY_OPERATOR(op)   \
	inline VT operator op(const VT &rhs) const \
	{                                          \
		return VT{_x op rhs._x, _y op rhs._y};   \
	}                                          \
	inline VT &&operator op(VT && rhs) const   \
	{                                          \
		rhs._x = _x op rhs._x;                   \
		rhs._y = _y op rhs._y;                   \
		return std::move(rhs);                   \
	}                                          \
	inline VT operator op(const T &rhs) const { return VT{_x op rhs, _y op rhs}; }

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
		return *this;                              \
	}                                            \
	inline VT &operator op(const T & rhs)        \
	{                                            \
		_x op rhs;                                 \
		_y op rhs;                                 \
		return *this;                              \
	}

	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(-=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(+=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(*=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(/=);

#undef ELEMENT_WISELY_COMPOUND_ASSIGNMENT

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)             \
	inline bool operator op(const VT &rhs) const \
	{                                            \
		return _x op rhs._x && _y op rhs._y;       \
	}                                            \
	inline bool operator op(const T &rhs) const { return _x op rhs && _y op rhs; }

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
	inline T dot(const VT &rhs) const { return _x * rhs._x + _y * rhs._y; }

	/// @brief Dot product of \p this and \p rhs , which is equivalent to \ref dot
	/// .
	inline T operator|(const VT &rhs) const { return this->dot(rhs); }

	/** @brief Calculate the determinant of a 2x2 matrix.
	 * @details The matrix is defined as
	 @verbatim
	    | this.x   this.y |
	    | rhs.x    rhs.y  |
	 @endverbatim
	*/
	inline T determinant(const VT &rhs) const
	{
		return _x * rhs._y - _y * rhs._x;
	}

	/// @brief Get the sum of all elements.
	inline T sum() const { return _x + _y; }

	/// @brief Get the squared sum of all elements.
	inline T sqrnorm() const { return _x * _x + _y * _y; }

	/// @brief Minimize \p this by the other vector \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \min(x_i, y_i) \f).
	inline void minimize(const VT &rhs)
	{
		if (rhs._x < _x)
			_x = rhs._x;
		if (rhs._y < _y)
			_y = rhs._y;
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
	}

	/// @brief Vectorize \p this by a given number \p rhs ,
	/// i.e., set all elements of \p this to \p rhs.
	inline VT &vectorize(const T &s)
	{
		_x = s;
		_y = s;
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

	/// @brief Get the \f( \ell_2 \f)-norm of \p this , which is equivalent to
	/// \ref norm.
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
		return VT{std::pow(_x, exp), std::pow(_y, exp)};
	}

	/// @brief Get the \f( \ell_p \f)-norm of \p this .
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline T pnorm(const T &exp) const
	{
		T px = std::pow(std::abs(_x), exp);
		T py = std::pow(std::abs(_y), exp);
		return std::pow(px + py, 1. / exp);
	}

	/// @brief Check if all elements are finite.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline bool isfinite() const
	{
		return std::isfinite(_x) && std::isfinite(_y);
	}

	/// @brief Get element-wisely floor result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT floor() const
	{
		return VT{std::floor(_x), std::floor(_y)};
	}

	/// @brief Get element-wisely ceil result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT ceil() const
	{
		return VT{std::ceil(_x), std::ceil(_y)};
	}

	/// @brief Get element-wisely round result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT round() const
	{
		return VT{std::round(x()), std::round(y())};
	}

	/// @brief Get element-wisely absolute result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT abs() const
	{
		return VT{std::abs(_x), std::abs(_y)};
	}

	/// @brief Negate elements in-place.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT &negate()
	{
		_x = -_x;
		_y = -_y;
		return *this;
	}

protected:
	T _x;
	T _y;
}; // namespace Geometry

/******* Below operators are optimized for right values. *******/

template <typename T>
inline Vec2T<T> operator-(Vec2T<T> &&lhs)
{
	lhs.negate();
	return std::move(lhs);
}

/// This micro defines two functions for a given binary operator.
/// The functions apply element wisely binary operator to \p this and \p rhs .
/// They are specially optimized for right values.
#define ELEMENT_WISELY_BINARY_OPERATOR_RV(op, compound_op)           \
	template <typename T>                                              \
	inline Vec2T<T> &&operator op(Vec2T<T> &&lhs, const Vec2T<T> &rhs) \
	{                                                                  \
		return std::move(lhs compound_op rhs);                           \
	}                                                                  \
	template <typename T>                                              \
	inline Vec2T<T> &&operator op(Vec2T<T> &&lhs, Vec2T<T> &&rhs)      \
	{                                                                  \
		return std::move(lhs compound_op rhs);                           \
	}                                                                  \
	template <typename T>                                              \
	inline Vec2T<T> &&operator op(Vec2T<T> &&lhs, const T & rhs)       \
	{                                                                  \
		return std::move(lhs compound_op rhs);                           \
	}

ELEMENT_WISELY_BINARY_OPERATOR_RV(-, -=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(+, +=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(*, *=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(/, /=);

#undef ELEMENT_WISELY_BINARY_OPERATOR_RV

/// @brief Multiply a number \p lhs and a Vec2d \p rhs .
template <typename T>
inline Vec2T<T> operator*(const T &lhs, const Vec2T<T> &rhs)
{
	return rhs * lhs;
}

/// @brief Multiply a number \p lhs and a Vec2d \p rhs .
template <typename T>
inline Vec2T<T> &&operator*(const T &lhs, Vec2T<T> &&rhs)
{
	return std::move(rhs) * lhs;
}

/// @brief A number \p lhs divede a Vec3T \p rhs element-wisely.
/// result[i] = lhs / rhs[i].
template <typename T>
inline Vec2T<T> operator/(const T &lhs, const Vec2T<T> &rhs)
{
	return Vec2T<T>{lhs / rhs.x(), lhs / rhs.y()};
}

/// @brief A number \p lhs divede a Vec3T \p rhs element-wisely.
/// result[i] = lhs / rhs[i].
template <typename T>
inline Vec2T<T> &&operator/(const T &lhs, Vec2T<T> &&rhs)
{
	rhs.x() = lhs / rhs.x();
	rhs.y() = lhs / rhs.y();
	return std::move(rhs);
}

using Vec2f = Vec2T<float>;
using Vec2d = Vec2T<double>;
using Vec2i = Vec2T<int>;
using Vec2u = Vec2T<unsigned int>;

} // namespace OMC