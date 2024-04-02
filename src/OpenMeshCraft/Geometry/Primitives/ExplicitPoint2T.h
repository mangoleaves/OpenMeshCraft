#pragma once

#include "GenericPoint2T.h"

namespace OMC {

/**
 * @brief Explicit point in 2D that contains one Vec2 as its coordinate.
 * @tparam IT interval number type.
 * @tparam ET exact number type.
 * @note This type is almost equivalent to Point2T defined in
 * Geometry/Basic/Primitive, the only difference is that this type inherits from
 * GenericPoint2T and will be used in implicit point context.
 */
template <typename IT_, typename ET_>
class ExplicitPoint2T : public GenericPoint2T<IT_, ET_>
{
public: /* types *************************************************************/
	using NT = double;
	using IT = IT_;
	using ET = ET_;

	using VecT = Vec2T<NT>;

	using GP = GenericPoint2T<IT, ET>;

	using PointType = typename GP::PointType;

public: /* Constructors ******************************************************/
	ExplicitPoint2T() noexcept
	  : GP(PointType::Explicit)
	  , m_p()
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint2T object, vectorize it by \p v .
	 * @tparam
	 */
	explicit ExplicitPoint2T(const NT &v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(v)
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint2T object from a pointer to number \p v
	 * @note We assume the pointer points to at least two numbers and won't
	 * check the validity of the pointer.
	 * @param v the pointer to number.
	 */
	explicit ExplicitPoint2T(const NT *v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(v)
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint2T object from VecT.
	 * @tparam
	 */
	template <typename TT, typename = std::enable_if_t<
	                         std::is_same_v<remove_cvref_t<TT>, VecT>> /*SIFNAE*/>
	explicit ExplicitPoint2T(TT &&v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(std::forward<TT>(v))
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint2T object from two numbers.
	 */
	template <
	  typename T1, typename T2,
	  typename = std::enable_if_t<std::is_same_v<remove_cvref_t<T1>, NT> &&
	                              std::is_same_v<remove_cvref_t<T2>, NT>>>
	explicit ExplicitPoint2T(T1 &&x, T2 &&y) noexcept
	  : GP(PointType::Explicit)
	  , m_p(std::forward<T1>(x), std::forward<T2>(y))
	{
	}

	ExplicitPoint2T(const ExplicitPoint2T &v)            = default;
	ExplicitPoint2T(ExplicitPoint2T &&v)                 = default;
	ExplicitPoint2T &operator=(const ExplicitPoint2T &v) = default;
	ExplicitPoint2T &operator=(ExplicitPoint2T &&v)      = default;

	~ExplicitPoint2T() = default;

public: /* Regular operators *************************************************/
	/// @name Data access
	/// @{
	NT       &x() { return m_p.x(); }
	NT       &y() { return m_p.y(); }
	const NT &x() const { return m_p.x(); }
	const NT &y() const { return m_p.y(); }

	NT       &operator[](size_t i) { return m_p[i]; }
	const NT &operator[](size_t i) const { return m_p[i]; }

	/// @brief Get the coordinate pointer.
	const NT *data() const { return m_p.data(); }

	/// @brief By as_vec, you can get a reference and operate it as VecT.
	/// Any modification to the corresponding VecT is also applied to this
	/// ExplicitPoint2T.
	VecT &as_vec() { return m_p; }

	/// @brief By as_vec, you can get a ference and operate it as VecT.
	const VecT &as_vec() const { return m_p; }

	/// @brief By to_vec, you can get a copy of this ExplicitPoint2T,
	/// whose type is now VecT.
	VecT to_vec() { return m_p; }
	/// @}

	/**
	 * @brief negative
	 * @return ExplicitPoint2T the negative point.
	 */
	ExplicitPoint2T operator-() const { return ExplicitPoint2T(-as_vec()); }

	/**
	 * @brief vector = this point - the other point
	 * @param rhs the other point
	 * @return VecT the result vector between two points
	 */
	VecT operator-(const ExplicitPoint2T &rhs) const
	{
		return as_vec() - rhs.as_vec();
	}

	/**
	 * @brief new point = this point + another point
	 * @param rhs the given point
	 * @return ExplicitPoint2T the result point.
	 */
	ExplicitPoint2T operator+(const ExplicitPoint2T &rhs) const
	{
		return ExplicitPoint2T(as_vec() + rhs.as_vec());
	}

	/**
	 * @brief new point = this point + given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint2T the result point by adding a vector to this point.
	 */
	ExplicitPoint2T operator+(const VecT &rhs) const
	{
		return ExplicitPoint2T(as_vec() + rhs);
	}

	/**
	 * @brief new point = this point - given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint2T the result point by adding a vector to this point.
	 */
	ExplicitPoint2T operator-(const VecT &rhs) const
	{
		return ExplicitPoint2T(as_vec() - rhs);
	}

	/**
	 * @brief new point = point * scale_factor
	 * @param rhs the scale factor
	 * @return ExplicitPoint2T the scaled new point
	 */
	ExplicitPoint2T operator*(const NT &scale_factor) const
	{
		return ExplicitPoint2T(as_vec() * scale_factor);
	}

	/**
	 * @brief this point += another point
	 * @param rhs the given point
	 * @return ExplicitPoint2T the result point.
	 */
	ExplicitPoint2T &operator+=(const ExplicitPoint2T &rhs)
	{
		as_vec() += rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point += given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint2T the result point by adding a vector to this point.
	 */
	ExplicitPoint2T &operator+=(const VecT &rhs)
	{
		as_vec() += rhs;
		return *this;
	}

	/**
	 * @brief this point -= another point
	 * @param rhs the given point
	 * @return ExplicitPoint2T the result point.
	 */
	ExplicitPoint2T &operator-=(const ExplicitPoint2T &rhs)
	{
		as_vec() -= rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point -= given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint2T the result point by adding a vector to this point.
	 */
	ExplicitPoint2T &operator-=(const VecT &rhs)
	{
		as_vec() -= rhs;
		return *this;
	}

	/**
	 * @brief point *= scale_factor
	 * @param rhs the scale factor
	 * @return ExplicitPoint2T the in-place scaled point
	 */
	ExplicitPoint2T operator*=(const NT &scale_factor)
	{
		as_vec() *= scale_factor;
		return *this;
	}

	/// @brief Minimize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \min(x_i, y_i) \f).
	void minimize(const ExplicitPoint2T &rhs) { m_p.minimize(rhs.m_p); }

	/// @brief Maximize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \max(x_i, y_i) \f).
	void maximize(const ExplicitPoint2T &rhs) { m_p.maximize(rhs.m_p); }

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)                          \
	inline bool operator op(const ExplicitPoint2T &rhs) const \
	{                                                         \
		return m_p op rhs.m_p;                                  \
	}                                                         \
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
 * @return ExplicitPoint2T the scaled new point
 */
template <typename IT, typename ET>
inline ExplicitPoint2T<IT, ET> operator*(const double &scale_factor,
                                         const ExplicitPoint2T<IT, ET> &point)
{
	return point * scale_factor;
}

} // namespace OMC