#pragma once

#include "Primitive2.h"
#include "Vector2T.h"

#include <cassert>
#include <variant>

namespace OMC {

template <typename IT, typename ET>
class ExplicitPoint2T;

template <typename IT, typename ET>
class ImplicitPoint2T_SSI;

/// @brief the generic point of 2D exact implicit and explicit points
template <typename IT, typename ET>
class GenericPoint2T : public Primitive2<double>
{
public: /* types *************************************************************/
	using FT   = double;
	using VecT = Vec2T<FT>;

	using EP     = ExplicitPoint2T<IT, ET>;
	using IP_SSI = ImplicitPoint2T_SSI<IT, ET>;

public: /* functions about types **********************************************/
	enum class PointType : uint32_t
	{
		Explicit = 0,
		Implicit = 1,
		SSI      = 2
	};

	PointType m_point_type;

	/// @brief get the point type
	inline PointType point_type() const { return m_point_type; };

	bool is_Explicit() const { return point_type() == PointType::Explicit; }
	bool is_Implicit() const { return point_type() > PointType::Implicit; }
	bool is_SSI() const { return point_type() == PointType::SSI; }

	bool has_ssf() const { return true; }

	/// @brief Convert to explicit point, won't check type again.
	EP &EXP()
	{
		OMC_EXPENSIVE_ASSERT(is_Explicit(), "point type mismatch.");
		return *static_cast<EP *>(this);
	}
	/// @brief Convert to explicit point, won't check type again.
	const EP &EXP() const
	{
		OMC_EXPENSIVE_ASSERT(is_Explicit(), "point type mismatch.");
		return *static_cast<const EP *>(this);
	}

	/// @brief Convert to SSI point, won't check type again.
	IP_SSI &SSI()
	{
		OMC_EXPENSIVE_ASSERT(is_SSI(), "point type mismatch.");
		return *static_cast<IP_SSI *>(this);
	}
	/// @brief Convert to SSI point, won't check type again.
	const IP_SSI &SSI() const
	{
		OMC_EXPENSIVE_ASSERT(is_SSI(), "point type mismatch.");
		return *static_cast<const IP_SSI *>(this);
	}

	/// @brief approximate this generic point by an explicit point
	void get_Explicit(EP &e) const
	{
		if (point_type() == PointType::Explicit)
			e = EXP();
		else /*if (point_type() == PointType::SSI)*/
			SSI().get_Explicit(e);
	}

	/// @brief approximate this generic point by an explicit point
	EP to_Explicit() const
	{
		EP e;
		get_Explicit(e);
		return e;
	}

	/// @brief Get approximate values on all dimensions of this generic point.
	void get_coordinates(FT &x, FT &y) const
	{
		EP e_;
		get_Explicit(e_);
		x = e_.x();
		y = e_.y();
	}

public: /* Constructor and Destructor ****************************************/
	GenericPoint2T(PointType pt)
	  : m_point_type(pt)
	{
	}

	virtual ~GenericPoint2T() noexcept {}

	GenericPoint2T(const GenericPoint2T &gp)
	  : m_point_type(gp.m_point_type)
	{
	}

	GenericPoint2T &operator=(const GenericPoint2T &gp)
	{
		m_point_type = gp.m_point_type;
		return *this;
	}

public: /* get lambda values from implicit points ****************************/
	bool getFilteredLambda(FT &lx, FT &ly, FT &d, FT &mv) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		/*if (point_type() == PointType::SSI)*/
		return SSI().getFilteredLambda(lx, ly, d, mv);
	}

	template <typename _IT = IT,
	          typename     = std::enable_if_t<!std::is_void_v<_IT>>>
	bool getIntervalLambda(_IT &lx, _IT &ly, _IT &d) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		/*if (point_type() == PointType::SSI)*/
		return SSI().getIntervalLambda(lx, ly, d);
	}

	template <typename _ET = ET,
	          typename     = std::enable_if_t<!std::is_void_v<_ET>>>
	void getExactLambda(_ET &lx, _ET &ly, _ET &d) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		/*if (point_type() == PointType::SSI)*/
		return SSI().getExactLambda(lx, ly, d);
	}

	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **d,
	                        int &d_len) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		/*if (point_type() == PointType::SSI)*/
		return SSI().getExpansionLambda(lx, lx_len, ly, ly_len, d, d_len);
	}

public:
	/***********************************************************************
	 * Below functions are wrappers of functions from explicit point.
	 * NOT safe!! just for efficiency :(
	 ***********************************************************************/

	/// @brief Access value on x. Assume this is explicit when call the func.
	FT       &x() { return EXP().x(); }
	/// @brief Access value on x. Assume this is explicit when call the func.
	const FT &x() const { return EXP().x(); }

	/// @brief Access value on y. Assume this is explicit when call the func.
	FT       &y() { return EXP().y(); }
	/// @brief Access value on y. Assume this is explicit when call the func.
	const FT &y() const { return EXP().y(); }

	/// @brief Get the value on dimension \p i. Assume this is explicit point.
	FT       &operator[](size_t i) { return EXP()[i]; }
	/// @brief Get the value on dimension \p i. Assume this is explicit point.
	const FT &operator[](size_t i) const { return EXP()[i]; }

	/// @brief Get the coordinate pointer.
	const FT *data() const { return EXP().data(); }

	/// @brief use point as vector
	VecT       &as_vec() { return EXP().as_vec(); }
	/// @brief use point as vector
	const VecT &as_vec() const { return EXP().as_vec(); }

	/// @brief convert point to vector
	VecT to_vec() const { return EXP().to_vec(); }

	/// @brief negative
	GenericPoint2T operator-() const { return GenericPoint2T(-as_vec()); }
	/// @brief calculate the vector from \p rhs to \p this.
	/// vec = this - rhs.
	VecT           operator-(const GenericPoint2T &rhs) const
	{
		return as_vec() - rhs.as_vec();
	}
	/// @brief new point = this point + another point
	GenericPoint2T operator+(const GenericPoint2T &rhs) const
	{
		return GenericPoint2T(as_vec() + rhs.as_vec());
	}
	/// @brief add this point with a vector
	GenericPoint2T operator+(const VecT &rhs) const
	{
		return GenericPoint2T(as_vec() + rhs);
	}
	/// @brief subtract this point with a vector
	GenericPoint2T operator-(const VecT &rhs) const
	{
		return GenericPoint2T(as_vec() - rhs);
	}
	/// @brief this point += another point
	GenericPoint2T &operator+=(const GenericPoint2T &rhs)
	{
		as_vec() += rhs.as_vec();
		return *this;
	}
	/// @brief add this point with a vector
	GenericPoint2T &operator+=(const VecT &rhs)
	{
		as_vec() += rhs;
		return *this;
	}
	/// @brief this point -= another point
	GenericPoint2T &operator-=(const GenericPoint2T &rhs)
	{
		as_vec() -= rhs.as_vec();
		return *this;
	}
	/// @brief subtract this point with a vector
	GenericPoint2T &operator-=(const VecT &rhs)
	{
		as_vec() -= rhs;
		return *this;
	}
	/// @brief multiply a scalar
	GenericPoint2T &operator*=(const FT &scale_factor)
	{
		as_vec() *= scale_factor;
		return *this;
	}
	/// @brief divide a scalar
	GenericPoint2T &operator/=(const FT &scale_factor)
	{
		as_vec() /= scale_factor;
		return *this;
	}

	/// @brief set this point with minimal values between this and another point
	void minimize(const GenericPoint2T &rhs) { as_vec().minimize(rhs.as_vec()); }
	/// @brief set this point with maximal values between this and another point
	void maximize(const GenericPoint2T &rhs) { as_vec().maximize(rhs.as_vec()); }

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)                         \
	inline bool operator op(const GenericPoint2T &rhs) const \
	{                                                        \
		return as_vec() op rhs.as_vec();                       \
	}                                                        \
	inline bool operator op(const FT &rhs) const { return as_vec() op rhs; }

	ELEMENT_WISELY_COMPARE(==);
	ELEMENT_WISELY_COMPARE(!=);
	ELEMENT_WISELY_COMPARE(<);
	ELEMENT_WISELY_COMPARE(<=);
	ELEMENT_WISELY_COMPARE(>);
	ELEMENT_WISELY_COMPARE(>=);
#undef ELEMENT_WISELY_COMPARE

public: /* Global cache lambda values ****************************************/
	static void enable_global_cached_values(size_t s = 0)
	{
		IP_SSI::gcv().enable();
		if (s != 0)
			resize_global_cached_values(s);
	}

	static void disable_global_cached_values() { IP_SSI::gcv().disable(); }

	static void clear_global_cached_values()
	{
		IP_SSI::gcv().clear_cached_values();
	}

	static void resize_global_cached_values(size_t s) { IP_SSI::gcv().resize(s); }

	static bool global_cached_values_enabled()
	{
		return IP_SSI::gcv().is_enabled();
	}
};

} // namespace OMC