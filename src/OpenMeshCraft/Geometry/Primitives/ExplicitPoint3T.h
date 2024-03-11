#pragma once

#include "GenericPoint3T.h"

namespace OMC {

/**
 * @brief Explicit point in 3D that contains one Vec3 as its coordinate.
 * @tparam IT interval number type.
 * @tparam ET exact number type.
 * @note This type is almost equivalent to Point3T defined in
 * Geometry/Basic/Primitive, the only difference is that this type inherits from
 * GenericPoint3T and will be used in implicit point context.
 */
template <typename IT_, typename ET_>
class ExplicitPoint3T : public GenericPoint3T<IT_, ET_>
{
public: /* types *************************************************************/
	using NT = double;
	using IT = IT_;
	using ET = ET_;

	using VecT = Vec3T<NT>;

	using GP = GenericPoint3T<IT, ET>;

	using PointType = typename GP::PointType;

public:
	ExplicitPoint3T() noexcept
	  : GP(PointType::Explicit)
	  , m_p()
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint3T object, vectorize it by \p v .
	 */
	explicit ExplicitPoint3T(const NT &v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(v)
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint3T object from a pointer to number \p v
	 * .
	 * @note We assume the pointer points to at least three numbers and won't
	 * check the validity of the pointer.
	 * @param v the pointer to number.
	 */
	explicit ExplicitPoint3T(const NT *v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(v)
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint3T object from VecT.
	 */
	template <typename TT, typename = std::enable_if_t<
	                         std::is_same_v<remove_cvref_t<TT>, VecT>> /*SIFNAE*/>
	ExplicitPoint3T(TT &&v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(std::forward<TT>(v))
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint3T object from three numbers.
	 */
	template <
	  typename T1, typename T2, typename T3,
	  typename = std::enable_if_t<std::is_same_v<remove_cvref_t<T1>, NT> &&
	                              std::is_same_v<remove_cvref_t<T2>, NT> &&
	                              std::is_same_v<remove_cvref_t<T3>, NT>>>
	explicit ExplicitPoint3T(T1 &&x, T2 &&y, T3 &&z) noexcept
	  : GP(PointType::Explicit)
	  , m_p(std::forward<T1>(x), std::forward<T2>(y), std::forward<T3>(z))
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint3T object by copying from a given
	 * ExplicitPoint3T \p v
	 * .
	 * @param v The given ExplicitPoint3T.
	 */
	ExplicitPoint3T(const ExplicitPoint3T &v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(v.m_p)
	{
	}

	/**
	 * @brief Construct a new ExplicitPoint3T object by copying from a given
	 * ExplicitPoint3T \p v
	 * .
	 * @param v The given ExplicitPoint3T.
	 */
	ExplicitPoint3T(ExplicitPoint3T &&v) noexcept
	  : GP(PointType::Explicit)
	  , m_p(std::move(v.m_p))
	{
	}

	virtual ~ExplicitPoint3T() noexcept {}

	/**
	 * @brief Construct a new ExplicitPoint3T object by copying from a given
	 * ExplicitPoint3T \p v
	 * .
	 * @param v The given ExplicitPoint3T.
	 * @return ExplicitPoint3T& return the reference to this object.
	 */
	ExplicitPoint3T &operator=(const ExplicitPoint3T &v) noexcept
	{
		this->m_point_type = v.m_point_type;
		m_p                = v.m_p;
		return *this;
	}

	/**
	 * @brief Construct a new ExplicitPoint3T object by copying from a given
	 * ExplicitPoint3T \p v
	 * .
	 * @param v The given ExplicitPoint3T.
	 * @return ExplicitPoint3T& return the reference to this object.
	 */
	ExplicitPoint3T &operator=(ExplicitPoint3T &&v) noexcept
	{
		this->m_point_type = v.m_point_type;
		m_p                = std::move(v.m_p);
		return *this;
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

	/// @brief Get the coordinate pointer.
	const NT *data() const { return m_p.data(); }

	/// @brief By as_vec, you can get a reference and operate it as VecT.
	/// Any modification to the corresponding VecT is also applied to this
	/// ExplicitPoint3T.
	VecT &as_vec() { return m_p; }

	/// @brief By as_vec, you can get a ference and operate it as VecT.
	const VecT &as_vec() const { return m_p; }

	/// @brief By to_vec, you can get a copy of this ExplicitPoint3T,
	/// whose type is now VecT.
	VecT to_vec() { return m_p; }
	/// @}

	/**
	 * @brief negative
	 * @return ExplicitPoint3T the negative point.
	 */
	ExplicitPoint3T operator-() const { return ExplicitPoint3T(-as_vec()); }

	/**
	 * @brief vector = this point - the other point
	 * @param rhs the other point
	 * @return VecT the result vector between two points
	 */
	VecT operator-(const ExplicitPoint3T &rhs) const
	{
		return as_vec() - rhs.as_vec();
	}

	/**
	 * @brief new point = this point + another point
	 * @param rhs the given point
	 * @return ExplicitPoint3T the result point.
	 */
	ExplicitPoint3T operator+(const ExplicitPoint3T &rhs) const
	{
		return ExplicitPoint3T(as_vec() + rhs.as_vec());
	}

	/**
	 * @brief new point = this point + given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint3T the result point by adding a vector to this point.
	 */
	ExplicitPoint3T operator+(const VecT &rhs) const
	{
		return ExplicitPoint3T(as_vec() + rhs);
	}

	/**
	 * @brief new point = this point - given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint3T the result point by adding a vector to this point.
	 */
	ExplicitPoint3T operator-(const VecT &rhs) const
	{
		return ExplicitPoint3T(as_vec() - rhs);
	}

	/**
	 * @brief new point = point * scale_factor
	 * @param rhs the scale factor
	 * @return ExplicitPoint3T the scaled new point
	 */
	ExplicitPoint3T operator*(const NT &scale_factor) const
	{
		return ExplicitPoint3T(as_vec() * scale_factor);
	}

	/**
	 * @brief this point += another point
	 * @param rhs the given point
	 * @return ExplicitPoint3T the result point
	 */
	ExplicitPoint3T &operator+=(const ExplicitPoint3T &rhs)
	{
		as_vec() += rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point += given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint3T the result point by adding a vector to this point.
	 */
	ExplicitPoint3T &operator+=(const VecT &rhs)
	{
		as_vec() += rhs;
		return *this;
	}

	/**
	 * @brief this point -= another point
	 * @param rhs the given point
	 * @return ExplicitPoint3T the result point.
	 */
	ExplicitPoint3T &operator-=(const ExplicitPoint3T &rhs)
	{
		as_vec() -= rhs.as_vec();
		return *this;
	}

	/**
	 * @brief this point -= given vector
	 * @param rhs the given vector
	 * @return ExplicitPoint3T the result point by adding a vector to this point.
	 */
	ExplicitPoint3T &operator-=(const VecT &rhs)
	{
		as_vec() -= rhs;
		return *this;
	}

	/**
	 * @brief point *= scale_factor
	 * @param rhs the scale factor
	 * @return ExplicitPoint3T the in-place scaled point
	 */
	ExplicitPoint3T operator*=(const NT &scale_factor)
	{
		as_vec() *= scale_factor;
		return *this;
	}

	/// @brief Minimize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \min(x_i, y_i) \f).
	void minimize(const ExplicitPoint3T &rhs) { m_p.minimize(rhs.m_p); }

	/// @brief Maximize \p this by the other point \p rhs .
	/// For each element \f( x_i \f) of \p this and \f( y_i \f) of \p rhs ,
	/// \f( x_i \f) is set to \f( \max(x_i, y_i) \f).
	void maximize(const ExplicitPoint3T &rhs) { m_p.maximize(rhs.m_p); }

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)                          \
	inline bool operator op(const ExplicitPoint3T &rhs) const \
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
 * @return ExplicitPoint3T the scaled new point
 */
template <typename IT, typename ET>
inline ExplicitPoint3T<IT, ET> operator*(const double &scale_factor,
                                         const ExplicitPoint3T<IT, ET> &point)
{
	return point * scale_factor;
}

} // namespace OMC