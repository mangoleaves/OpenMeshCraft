#pragma once

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

#include <array>
#include <cmath>
#include <type_traits>

namespace OMC {

/**
 * @ingroup Geometry
 * @brief VecXT, containing `X` numbers of type `T`, is an precision-limited
 * type of general VectorXT.
 *
 * @tparam X number of numbers
 * @tparam T type of numbers
 */
template <size_t X, typename T>
class VecXT
{
	static_assert(X >= 2, "VecXT need at least 2 elements.");

public:
	using NT                      = T;
	using VT                      = VecXT<X, T>;
	static constexpr size_t size_ = X;
	using value_type              = T;
	using vector_type             = VT;

private:
	using internal_type = std::array<T, X>;

public:
	VecXT() noexcept {}
	VecXT(const VT &rhs) noexcept
	  : values(rhs.values)
	{
	}
	VecXT(VT &&rhs) noexcept
	  : values(std::move(rhs.values))
	{
	}
	VT &operator=(const VT &rhs) noexcept
	{
		values = rhs.values;
		return *this;
	}
	VT &operator=(VT &&rhs) noexcept
	{
		values = std::move(rhs.values);
		return *this;
	}
	~VecXT() noexcept {}

	explicit inline VecXT(const T *v) noexcept
	{
		const T *first = v;
		for (auto it = values.begin(); it != values.end(); it++, first++)
			*it = *first;
	}

	explicit inline VecXT(const T &v0) noexcept { vectorize(v0); }

	template <typename T1, typename T2>
	explicit inline VecXT(T1 &&v0, T2 &&v1) noexcept
	  : values{std::forward<T1>(v0), std::forward<T2>(v1)}
	{
	}

	template <typename T1, typename T2, typename T3>
	explicit inline VecXT(T1 &&v0, T2 &&v1, T3 &&v2) noexcept
	  : values{std::forward<T1>(v0), std::forward<T2>(v1), std::forward<T3>(v2)}
	{
	}

	template <typename T1, typename T2, typename T3, typename T4>
	explicit inline VecXT(T1 &&v0, T2 &&v1, T3 &&v2, T4 &&v3) noexcept
	  : values{std::forward<T1>(v0), std::forward<T2>(v1), std::forward<T3>(v2),
	           std::forward<T4>(v3)}
	{
	}

	/** @name Data access
	 * */

	/// @{
	inline T &x() { return values[0]; }
	inline T &y() { return values[1]; }

	template <size_t X_ = X, typename = std::enable_if_t<X_ >= 3>> // SFINAE
	inline T &z()
	{
		return values[2];
	}

	template <size_t X_ = X, typename = std::enable_if_t<X_ >= 4>> // SFINAE
	inline T &w()
	{
		return values[3];
	}

	inline const T &x() const { return values[0]; }
	inline const T &y() const { return values[1]; }

	template <size_t X_ = X, typename = std::enable_if_t<X_ >= 3>>
	inline const T &z() const
	{
		return values[2];
	}

	template <size_t X_ = X, typename = std::enable_if_t<X_ >= 4>>
	inline const T &w() const
	{
		return values[3];
	}

	/// @brief size(dimention) of vector
	size_t size() const { return size_; }

	inline T       *data() { return values.data(); }
	inline const T *data() const { return values.data(); }

	inline T &operator[](size_t dim)
	{
		OMC_EXPENSIVE_ASSERT(dim < size_, "dim {} out of range {}.", dim, size_);
		return values[dim];
	}
	inline const T &operator[](size_t dim) const
	{
		OMC_EXPENSIVE_ASSERT(dim < size_, "dim {} out of range {}.", dim, size_);
		return values[dim];
	}
	/// @}

	/** @name Arithmetic operators
	 */
	/// @{

	/// This micro defines two functions for a given binary operator.
	/// The functions apply element wisely binary operator to \p this and \p rhs .
#define ELEMENT_WISE_BINARY_OPERATOR(op)             \
	VT operator op(const VT &rhs) const                \
	{                                                  \
		VT result;                                       \
		for (size_t i = 0; i < X; i++)                   \
			result.values[i] = values[i] op rhs.values[i]; \
		return result;                                   \
	}                                                  \
	VT &&operator op(VT && rhs) const                  \
	{                                                  \
		for (size_t i = 0; i < X; i++)                   \
			rhs.values[i] = values[i] op rhs.values[i];    \
		return std::move(rhs);                           \
	}                                                  \
	VT operator op(const T &rhs) const                 \
	{                                                  \
		VT result;                                       \
		for (size_t i = 0; i < X; i++)                   \
			result.values[i] = values[i] op rhs;           \
		return result;                                   \
	}                                                  \
	VT operator op(T &&rhs) const                      \
	{                                                  \
		VT result;                                       \
		for (size_t i = 0; i < X; i++)                   \
			result.values[i] = values[i] op rhs;           \
		return result;                                   \
	}

	ELEMENT_WISE_BINARY_OPERATOR(-);
	ELEMENT_WISE_BINARY_OPERATOR(+);
	ELEMENT_WISE_BINARY_OPERATOR(*);
	ELEMENT_WISE_BINARY_OPERATOR(/);

#undef ELEMENT_WISE_BINARY_OPERATOR

	/// This micro defines two functions for a given compound assignment operator.
	/// The functions apply element wisely binary operator to \p this and \p rhs
	/// and assign the result to \p this .
#define ELEMENT_WISELY_COMPOUND_ASSIGNMENT(op) \
	VT &operator op(const VT & rhs)              \
	{                                            \
		for (size_t i = 0; i < X; i++)             \
			values[i] op rhs.values[i];              \
		return *this;                              \
	}                                            \
	VT &operator op(const T & rhs)               \
	{                                            \
		for (size_t i = 0; i < X; i++)             \
			values[i] op rhs;                        \
		return *this;                              \
	}

	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(-=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(+=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(*=);
	ELEMENT_WISELY_COMPOUND_ASSIGNMENT(/=);

#undef ELEMENT_WISELY_COMPOUND_ASSIGNMENT

	/// @}

	/** @name Comparison operators
	 */
	/// @{

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op, reverse_op)   \
	bool operator op(const VecXT<X, T> &rhs) const \
	{                                              \
		for (size_t i = 0; i < X; i++)               \
			if (values[i] reverse_op rhs.values[i])    \
				return false;                            \
		return true;                                 \
	}                                              \
	bool operator op(const T &rhs) const           \
	{                                              \
		for (auto &left : values)                    \
			if (left reverse_op rhs)                   \
				return false;                            \
		return true;                                 \
	}

	ELEMENT_WISELY_COMPARE(==, !=);
	ELEMENT_WISELY_COMPARE(!=, ==);
	ELEMENT_WISELY_COMPARE(<, >=);
	ELEMENT_WISELY_COMPARE(<=, >);
	ELEMENT_WISELY_COMPARE(>, <=);
	ELEMENT_WISELY_COMPARE(>=, <);

#undef ELEMENT_WISELY_COMPARE

	bool less_on(size_t dim, const VT &rhs) const
	{
		return (this->operator[](dim) < rhs[dim]);
	}
	bool less_on(size_t dim, const T &rhs) const
	{
		return (this->operator[](dim) < rhs);
	}
	/// @}

	/** @name Batch operators
	 */
	/// @{
	void minimize(const VT &rhs)
	{
		for (size_t i = 0; i < X; i++)
			if (rhs.values[i] < values[i])
				values[i] = rhs.values[i];
	}

	void maximize(const VT &rhs)
	{
		for (size_t i = 0; i < X; i++)
			if (rhs.values[i] > values[i])
				values[i] = rhs.values[i];
	}

	VT &vectorize(const T &s)
	{
		for (auto &elem : values)
			elem = s;
		return *this;
	}
	/// @}

	/** @name Norm
	 */
	/// @{
	T dot(const VT &rhs) const
	{
		T result(0);
		for (size_t i = 0; i < X; i++)
			result += values[i] * rhs.values[i];
		return result;
	}
	T operator|(const VT &rhs) const { return dot(rhs); }
	T sqrnorm() const { return dot(*this); }

	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	T norm() const
	{
		return std::sqrt(sqrnorm());
	}

	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	T length() const
	{
		return norm();
	}

	/// @brief Get the element-wisely pow result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline VT pow(const T &exp) const
	{
		VT result;
		for (size_t i = 0; i < X; i++)
			result.values[i] = std::pow(values[i], exp);
		return result;
	}

	/// @brief Get the \f( \ell_p \f)-norm of \p this .
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline T pnorm(const T &exp) const
	{
		T result(0);
		for (auto &left : values)
			result += std::pow(std::abs(left), exp);
		return std::pow(result, 1. / exp);
	}
	/// @}

	/** @name Numeric
	 */
	/// @brief Check if all elements are finite.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	inline bool isfinite() const
	{
		for (auto &elem : values)
			if (!std::isfinite(elem))
				return false;
		return true;
	}

	/// @brief Get element-wisely floor result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT floor() const
	{
		VT result;
		for (size_t i = 0; i < X; i++)
			result.values[i] = std::floor(values[i]);
		return result;
	}

	/// @brief Get element-wisely ceil result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT ceil() const
	{
		VT result;
		for (size_t i = 0; i < X; i++)
			result.values[i] = std::ceil(values[i]);
		return result;
	}

	/// @brief Get element-wisely round result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT round() const
	{
		VT result;
		for (size_t i = 0; i < X; i++)
			result.values[i] = std::round(values[i]);
		return result;
	}

	/// @brief Get element-wisely absolute result.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT abs() const
	{
		VT result;
		for (size_t i = 0; i < X; i++)
			result.values[i] = std::abs(values[i]);
		return result;
	}

	/// @brief Negate elements in-place.
	template <typename U = T,
	          typename = std::enable_if_t<std::is_floating_point_v<U>>> // SFINAE
	VT &negate()
	{
		for (auto &res : values)
			res = -res;
		return *this;
	}
	/// @}

	template <size_t _X, typename _T>
	friend VecXT<_X, _T> operator/(const _T &lhs, const VecXT<_X, _T> &rhs);
	template <size_t _X, typename _T>
	friend VecXT<_X, _T> &&operator/(const _T &lhs, VecXT<_X, _T> &&rhs);

private:
	internal_type values;
};

/******* Below operators are optimized for right values. *******/

template <size_t X, typename T>
inline VecXT<X, T> operator-(VecXT<X, T> &&lhs)
{
	lhs.negate();
	return std::move(lhs);
}

/// This micro defines two functions for a given binary operator.
/// The functions apply element wisely binary operator to \p this and \p rhs .
/// They are specially optimized for right values.
#define ELEMENT_WISELY_BINARY_OPERATOR_RV(op, compound_op)                    \
	template <size_t X, typename T>                                             \
	inline VecXT<X, T> &&operator op(VecXT<X, T> &&lhs, const VecXT<X, T> &rhs) \
	{                                                                           \
		return std::move(lhs compound_op rhs);                                    \
	}                                                                           \
	template <size_t X, typename T>                                             \
	inline VecXT<X, T> &&operator op(VecXT<X, T> &&lhs, VecXT<X, T> &&rhs)      \
	{                                                                           \
		return std::move(lhs compound_op rhs);                                    \
	}                                                                           \
	template <size_t X, typename T>                                             \
	inline VecXT<X, T> &&operator op(VecXT<X, T> &&lhs, const T & rhs)          \
	{                                                                           \
		return std::move(lhs compound_op rhs);                                    \
	}

ELEMENT_WISELY_BINARY_OPERATOR_RV(-, -=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(+, +=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(*, *=);
ELEMENT_WISELY_BINARY_OPERATOR_RV(/, /=);

#undef ELEMENT_WISELY_BINARY_OPERATOR_RV

/// @brief Multiply a number \p lhs and a Vec3T \p rhs .
template <size_t X, typename T>
inline VecXT<X, T> operator*(const T &lhs, const VecXT<X, T> &rhs)
{
	return rhs * lhs;
}

/// @brief Multiply a number \p lhs and a Vec3T \p rhs .
template <size_t X, typename T>
inline VecXT<X, T> &&operator*(const T &lhs, VecXT<X, T> &&rhs)
{
	return std::move(rhs) * lhs;
}

/// @brief A number \p lhs divede a Vec3T \p rhs element-wisely.
/// result[i] = lhs / rhs[i].
template <size_t X, typename T>
inline VecXT<X, T> operator/(const T &lhs, const VecXT<X, T> &rhs)
{
	VecXT<X, T> result;
	for (size_t i = 0; i < X; i++)
		result.values[i] = lhs / rhs.values[i];
	return result;
}

/// @brief A number \p lhs divede a Vec3T \p rhs element-wisely.
/// result[i] = lhs / rhs[i].
template <size_t X, typename T>
inline VecXT<X, T> &&operator/(const T &lhs, VecXT<X, T> &&rhs)
{
	for (auto &right : rhs.values)
		right = lhs / right;
	return std::move(rhs);
}

// NOTE: To conviniencely traversal multiple ranges, we use boost::combine.
// But does it have the optimal performance? Need test.

} // namespace OMC