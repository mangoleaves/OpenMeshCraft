#pragma once

#include "Point2T.h"
#include "Point3T.h"

namespace OMC {

template <typename NT>
class AsGenericPoint_Ex
{
public:
	using GP2 = Point2T<NT>;
	using EP2 = Point2T<NT>;

	using GP3 = Point3T<NT>;
	using EP3 = Point3T<NT>;

public:
	// clang-format off
	GP2       &operator()(      EP2& src) { return *static_cast<GP2 *>(&src); }
	const GP2 &operator()(const EP2& src) { return *static_cast<const GP2 *>(&src); }

	GP3       &operator()(      EP3& src) { return *static_cast<GP3 *>(&src); }
	const GP3 &operator()(const EP3& src) { return *static_cast<const GP3 *>(&src); }
	// clang-format on
};

template <typename NT>
class AsExplicitPoint_Ex
{
public:
	using GP2 = Point2T<NT>;
	using EP2 = Point2T<NT>;

	using GP3 = Point3T<NT>;
	using EP3 = Point3T<NT>;

public:
	EP2       &operator()(GP2 &src) { return src; }
	EP3       &operator()(GP3 &src) { return src; }
	const EP2 &operator()(const GP2 &src) { return src; }
	const EP3 &operator()(const GP3 &src) { return src; }

	bool is_explicit(const GP2 &src) { return true; }
	bool is_explicit(const GP3 &src) { return true; }
};

template <typename NT>
class ToExplicitPoint_Ex
{
public:
	using GP2 = Point2T<NT>;
	using EP2 = Point2T<NT>;

	using GP3 = Point3T<NT>;
	using EP3 = Point3T<NT>;

public:
	const EP2 &operator()(const GP2 &src) { return src; }
	const EP3 &operator()(const GP3 &src) { return src; }
};

template <typename NT>
class CreateImplicitSSI2_Ex
{
public:
	using GP2     = Point2T<NT>;
	using EP2     = Point2T<NT>;
	using IP2_SSI = Point2T<NT>;

	IP2_SSI operator()(const GP2 &l11, const GP2 &l12, const GP2 &l21,
	                   const GP2 &l22)
	{
		OMC_THROW_NOT_IMPLEMENTED();
		return IP2_SSI();
	}
};

template <typename NT>
class CreateImplicitSSI3_Ex
{
public:
	using GP3     = Point3T<NT>;
	using EP3     = Point3T<NT>;
	using IP3_SSI = Point3T<NT>;

public:
	IP3_SSI operator()(const GP3 &a, const GP3 &b, const GP3 &p, const GP3 &q,
	                   int plane)
	{
		OMC_THROW_NOT_IMPLEMENTED();
		return IP3_SSI();
	}
};

template <typename NT>
class CreateImplicitLPI_Ex
{
public:
	using GP3     = Point3T<NT>;
	using EP3     = Point3T<NT>;
	using IP3_LPI = Point3T<NT>;

public:
	IP3_LPI operator()(const GP3 &p, const GP3 &q, const GP3 &r, const GP3 &s,
	                   const GP3 &t)
	{
		OMC_THROW_NOT_IMPLEMENTED();
		return IP3_LPI();
	}
};

template <typename NT>
class CreateImplicitTPI_Ex
{
public:
	using GP3     = Point3T<NT>;
	using EP3     = Point3T<NT>;
	using IP3_TPI = Point3T<NT>;

	IP3_TPI operator()(const GP3 &v1, const GP3 &v2, const GP3 &v3, const GP3 &w1,
	                   const GP3 &w2, const GP3 &w3, const GP3 &u1, const GP3 &u2,
	                   const GP3 &u3)
	{
		OMC_THROW_NOT_IMPLEMENTED();
		return IP3_TPI();
	}
};

} // namespace OMC