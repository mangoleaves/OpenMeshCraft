#pragma once

#include "Primitive3.h"
#include "Vector3T.h"

#include <cassert>
#include <variant>

namespace OMC {

template <typename IT, typename ET>
class ExplicitPoint3T;

template <typename IT, typename ET>
class ImplicitPoint3T_SSI;

template <typename IT, typename ET>
class ImplicitPoint3T_LPI;

template <typename IT, typename ET>
class ImplicitPoint3T_TPI;

template <typename IT, typename ET>
class ImplicitPoint3T_LNC;

/// @brief the generic point of 3D exact implicit and explicit points
template <typename IT, typename ET>
class GenericPoint3T : public Primitive3<double>
{
public: /* types *************************************************************/
	using FT   = double;
	using VecT = Vec3T<FT>;

	using EP     = ExplicitPoint3T<IT, ET>;
	using IP_SSI = ImplicitPoint3T_SSI<IT, ET>;
	using IP_LPI = ImplicitPoint3T_LPI<IT, ET>;
	using IP_TPI = ImplicitPoint3T_TPI<IT, ET>;
	using IP_LNC = ImplicitPoint3T_LNC<IT, ET>;

public: /* functions about types **********************************************/
	enum class PointType : uint32_t
	{ // make sure simpler point has smaller value
		Explicit = 0,
		Implicit = 1,
		SSI      = 2,
		LPI      = 3,
		TPI      = 4,
		LNC      = 5
	};

	PointType m_point_type;

	/// @brief get the point type
	inline PointType point_type() const { return m_point_type; }

	bool is_Explicit() const { return point_type() == PointType::Explicit; }
	bool is_Implicit() const { return point_type() > PointType::Implicit; }
	bool is_SSI() const { return point_type() == PointType::SSI; }
	bool is_LPI() const { return point_type() == PointType::LPI; }
	bool is_TPI() const { return point_type() == PointType::TPI; }
	bool is_LNC() const { return point_type() == PointType::LNC; }

	bool has_ssf() const { return !is_LNC(); }

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

	/// @brief Convert to LPI point, won't check type again.
	IP_LPI &LPI()
	{
		OMC_EXPENSIVE_ASSERT(is_LPI(), "point type mismatch.");
		return *static_cast<IP_LPI *>(this);
	}
	/// @brief Convert to LPI point, won't check type again.
	const IP_LPI &LPI() const
	{
		OMC_EXPENSIVE_ASSERT(is_LPI(), "point type mismatch.");
		return *static_cast<const IP_LPI *>(this);
	}

	/// @brief Convert to TPI point, won't check type again.
	IP_TPI &TPI()
	{
		OMC_EXPENSIVE_ASSERT(is_TPI(), "point type mismatch.");
		return *static_cast<IP_TPI *>(this);
	}
	/// @brief Convert to TPI point, won't check type again.
	const IP_TPI &TPI() const
	{
		OMC_EXPENSIVE_ASSERT(is_TPI(), "point type mismatch.");
		return *static_cast<const IP_TPI *>(this);
	}

	/// @brief Convert to LNC point, won't check type again.
	IP_LNC &LNC()
	{
		OMC_EXPENSIVE_ASSERT(is_LNC(), "point type mismatch.");
		return *static_cast<IP_LNC *>(this);
	}
	/// @brief Convert to LNC point, won't check type again.
	const IP_LNC &LNC() const
	{
		OMC_EXPENSIVE_ASSERT(is_LNC(), "point type mismatch.");
		return *static_cast<const IP_LNC *>(this);
	}

	/// @brief approximate this generic point by an explicit point
	void get_Explicit(EP &e) const
	{
		if (point_type() == PointType::Explicit)
			e = EXP();
		else if (point_type() == PointType::SSI)
			SSI().get_Explicit(e);
		else if (point_type() == PointType::LPI)
			LPI().get_Explicit(e);
		else if (point_type() == PointType::TPI)
			TPI().get_Explicit(e);
#if defined(INDIRECT_PREDICATES)
		else if (point_type() == PointType::LNC)
			LNC().get_Explicit(e);
#endif
	}

	/// @brief approximate this generic point by an explicit point
	EP to_Explicit() const
	{
		EP e;
		get_Explicit(e);
		return e;
	}

	/// @brief Get approximate values on all dimensions of this generic point.
	void get_coordinates(FT &x, FT &y, FT &z) const
	{
		EP e_;
		get_Explicit(e_);
		x = e_.x();
		y = e_.y();
		z = e_.z();
	}

public: /* Constructor and Destructor ***************************************/
	GenericPoint3T(PointType pt)
	  : m_point_type(pt)
	{
	}

	GenericPoint3T(const GenericPoint3T &gp)            = default;
	GenericPoint3T(GenericPoint3T &&gp)                 = default;
	GenericPoint3T &operator=(const GenericPoint3T &gp) = default;
	GenericPoint3T &operator=(GenericPoint3T &&gp)      = default;

protected:
	/// you can't delete a protected and nonvirtual derived class
	/// by delete the pointer to this base class.
	~GenericPoint3T() = default;

public: /* get lambda values from implicit points ****************************/
#if defined(INDIRECT_PREDICATES)
	bool getFilteredLambda(FT &lx, FT &ly, FT &lz, FT &d, FT &mv) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			return SSI().getFilteredLambda(lx, ly, lz, d, mv);
		else if (point_type() == PointType::LPI)
			return LPI().getFilteredLambda(lx, ly, lz, d, mv);
		else if (point_type() == PointType::TPI)
			return TPI().getFilteredLambda(lx, ly, lz, d, mv);
		else /*if (point_type() == PointType::LNC)*/
			return LNC().getFilteredLambda(lx, ly, lz, d, mv);
	}

	template <typename _IT = IT,
	          typename     = std::enable_if_t<!std::is_void_v<_IT>>>
	bool getIntervalLambda(_IT &lx, _IT &ly, _IT &lz, _IT &d) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			return SSI().getIntervalLambda(lx, ly, lz, d);
		else if (point_type() == PointType::LPI)
			return LPI().getIntervalLambda(lx, ly, lz, d);
		else if (point_type() == PointType::TPI)
			return TPI().getIntervalLambda(lx, ly, lz, d);
		else /*if (point_type() == PointType::LNC)*/
			return LNC().getIntervalLambda(lx, ly, lz, d);
	}

	template <typename _ET = ET,
	          typename     = std::enable_if_t<!std::is_void_v<_ET>>>
	void getExactLambda(_ET &lx, _ET &ly, _ET &lz, _ET &d) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			SSI().getExactLambda(lx, ly, lz, d);
		else if (point_type() == PointType::LPI)
			LPI().getExactLambda(lx, ly, lz, d);
		else if (point_type() == PointType::TPI)
			TPI().getExactLambda(lx, ly, lz, d);
		else /*if (point_type() == PointType::LNC)*/
			LNC().getExactLambda(lx, ly, lz, d);
	}

	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **lz,
	                        int &lz_len, FT **d, int &d_len) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			SSI().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len);
		else if (point_type() == PointType::LPI)
			LPI().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len);
		else if (point_type() == PointType::TPI)
			TPI().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len);
		else /*if (point_type() == PointType::LNC)*/
			LNC().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len);
	}
#elif defined(OFFSET_PREDICATES)
	bool getFilteredLambda(FT &lx, FT &ly, FT &lz, FT &d, FT &bx, FT &by, FT &bz,
	                       FT &mv) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			return SSI().getFilteredLambda(lx, ly, lz, d, bx, by, bz, mv);
		else if (point_type() == PointType::LPI)
			return LPI().getFilteredLambda(lx, ly, lz, d, bx, by, bz, mv);
		else /*if (point_type() == PointType::TPI)*/
			return TPI().getFilteredLambda(lx, ly, lz, d, bx, by, bz, mv);
		// else /*if (point_type() == PointType::LNC)*/
		//	return LNC().getFilteredLambda(lx, ly, lz, d, bx, by, bz, mv);
	}

	template <typename _IT = IT,
	          typename     = std::enable_if_t<!std::is_void_v<_IT>>>
	bool getIntervalLambda(_IT &lx, _IT &ly, _IT &lz, _IT &d, _IT &bx, _IT &by,
	                       _IT &bz) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			return SSI().getIntervalLambda(lx, ly, lz, d, bx, by, bz);
		else if (point_type() == PointType::LPI)
			return LPI().getIntervalLambda(lx, ly, lz, d, bx, by, bz);
		else /*if (point_type() == PointType::TPI)*/
			return TPI().getIntervalLambda(lx, ly, lz, d, bx, by, bz);
		// else /*if (point_type() == PointType::LNC)*/
		//	return LNC().getIntervalLambda(lx, ly, lz, d, bx, by, bz);
	}

	template <typename _ET = ET,
	          typename     = std::enable_if_t<!std::is_void_v<_ET>>>
	void getExactLambda(_ET &lx, _ET &ly, _ET &lz, _ET &d, _ET &bx, _ET &by,
	                    _ET &bz) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			SSI().getExactLambda(lx, ly, lz, d, bx, by, bz);
		else if (point_type() == PointType::LPI)
			LPI().getExactLambda(lx, ly, lz, d, bx, by, bz);
		else if (point_type() == PointType::TPI)
			TPI().getExactLambda(lx, ly, lz, d, bx, by, bz);
		// else /*if (point_type() == PointType::LNC)*/
		//	LNC().getExactLambda(lx, ly, lz, d, bx, by, bz);
	}

	void getExpansionLambda(FT **lx, int &lx_len, FT **ly, int &ly_len, FT **lz,
	                        int &lz_len, FT **d, int &d_len, FT &bx, FT &by,
	                        FT &bz) const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		// clang-format off
		if (point_type() == PointType::SSI)
			SSI().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len, bx, by, bz);
		else if (point_type() == PointType::LPI)
			LPI().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len, bx, by, bz);
		else if (point_type() == PointType::TPI)
			TPI().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len, bx, by, bz);
		// else /*if (point_type() == PointType::LNC)*/
		//	LNC().getExpansionLambda(lx, lx_len, ly, ly_len, lz, lz_len, d, d_len, bx, by, bz);
		// clang-format on
	}
#endif

	FT getIndirectMaxVar() const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			return SSI().getIndirectMaxVar();
		else if (point_type() == PointType::LPI)
			return LPI().getIndirectMaxVar();
		else /*if (point_type() == PointType::TPI)*/
			return TPI().getIndirectMaxVar();
	}

	FT getOffsetMaxVar() const
	{
		OMC_EXPENSIVE_ASSERT(point_type() != PointType::Explicit,
		                     "no lambda for explicit point");
		if (point_type() == PointType::SSI)
			return SSI().getOffsetMaxVar();
		else if (point_type() == PointType::LPI)
			return LPI().getOffsetMaxVar();
		else /*if (point_type() == PointType::TPI)*/
			return TPI().getOffsetMaxVar();
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

	/// @brief Access value on z. Assume this is explicit when call the func.
	FT       &z() { return EXP().z(); }
	/// @brief Access value on x. Assume this is explicit when call the func.
	const FT &z() const { return EXP().z(); }

	/// @brief Get the value on dimension \p i. Assume this is explicit point.
	FT       &operator[](size_t i) { return EXP()[i]; }
	/// @brief Get the value on dimension \p i. Assume this is explicit point.
	const FT &operator[](size_t i) const { return EXP()[i]; }

	/// @brief Get the coordinate pointer.
	const FT *data() const { return EXP().data(); }

	/// @brief use point as vector
	VecT &as_vec() { return EXP().as_vec(); }

	/// @brief use point as vector
	const VecT &as_vec() const { return EXP().as_vec(); }

	/// @brief convert point to vector
	VecT to_vec() const { return EXP().to_vec(); }

	/// @brief negative
	GenericPoint3T operator-() const { return GenericPoint3T(-as_vec()); }
	/// @brief calculate the vector from \p rhs to \p this.
	/// vec = this - rhs.
	VecT           operator-(const GenericPoint3T &rhs) const
	{
		return as_vec() - rhs.as_vec();
	}
	/// @brief new point = this point + another point
	GenericPoint3T operator+(const GenericPoint3T &rhs) const
	{
		return GenericPoint3T(as_vec() + rhs.as_vec());
	}
	/// @brief add this point with a vector
	GenericPoint3T operator+(const VecT &rhs) const
	{
		return GenericPoint3T(as_vec() + rhs);
	}
	/// @brief subtract this point with a vector
	GenericPoint3T operator-(const VecT &rhs) const
	{
		return GenericPoint3T(as_vec() - rhs);
	}
	/// @brief this point += another point
	GenericPoint3T &operator+=(const GenericPoint3T &rhs)
	{
		as_vec() += rhs.as_vec();
		return *this;
	}
	/// @brief add this point with a vector
	GenericPoint3T &operator+=(const VecT &rhs)
	{
		as_vec() += rhs;
		return *this;
	}
	/// @brief this point -= another point
	GenericPoint3T &operator-=(const GenericPoint3T &rhs)
	{
		as_vec() -= rhs.as_vec();
		return *this;
	}
	/// @brief subtract this point with a vector
	GenericPoint3T &operator-=(const VecT &rhs)
	{
		as_vec() -= rhs;
		return *this;
	}
	/// @brief multiply a scalar
	GenericPoint3T &operator*=(const FT &scale_factor)
	{
		as_vec() *= scale_factor;
		return *this;
	}
	/// @brief divide a scalar
	GenericPoint3T &operator/=(const FT &scale_factor)
	{
		as_vec() /= scale_factor;
		return *this;
	}

	/// @brief set this point with minimal values between this and another point
	void minimize(const GenericPoint3T &rhs) { as_vec().minimize(rhs.as_vec()); }
	/// @brief set this point with maximal values between this and another point
	void maximize(const GenericPoint3T &rhs) { as_vec().maximize(rhs.as_vec()); }

/// Compare \p this with \p rhs element-wisely.
#define ELEMENT_WISELY_COMPARE(op)                         \
	inline bool operator op(const GenericPoint3T &rhs) const \
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
		IP_LPI::gcv().enable();
		IP_TPI::gcv().enable();
#if defined(INDIRECT_PREDICATES)
		IP_LNC::gcv().enable();
#endif
		if (s != 0)
			resize_global_cached_values(s);
	}

	static void disable_global_cached_values()
	{
		IP_SSI::gcv().disable();
		IP_LPI::gcv().disable();
		IP_TPI::gcv().disable();
#if defined(INDIRECT_PREDICATES)
		IP_LNC::gcv().disable();
#endif
	}

	static void clear_global_cached_values()
	{
		IP_SSI::gcv().clear_cached_values();
		IP_LPI::gcv().clear_cached_values();
		IP_TPI::gcv().clear_cached_values();
#if defined(INDIRECT_PREDICATES)
		IP_LNC::gcv().clear_cached_values();
#endif
	}

	static void resize_global_cached_values(size_t s)
	{
		IP_SSI::gcv().resize(s);
		IP_LPI::gcv().resize(s);
		IP_TPI::gcv().resize(s);
#if defined(INDIRECT_PREDICATES)
		IP_LNC::gcv().resize(s);
#endif
	}

	static bool global_cached_values_enabled()
	{
		return IP_LPI::gcv().is_enabled(); // any one is ok
	}
};

} // namespace OMC