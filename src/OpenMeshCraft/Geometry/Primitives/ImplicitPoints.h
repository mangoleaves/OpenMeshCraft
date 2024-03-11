#pragma once

#include "GenericPoint2T.h"
#include "GenericPoint3T.h"

#include "ExplicitPoint2T.h"
#include "ExplicitPoint3T.h"

#include "ImplicitPoint2T.h"
#include "ImplicitPoint3T.h"

namespace OMC {

template <typename IT, typename ET>
class AsGenericPoint_Im
{
public:
	using GP2     = GenericPoint2T<IT, ET>;
	// below points inherit from GP2
	using EP2     = ExplicitPoint2T<IT, ET>;
	using IP2_SSI = ImplicitPoint2T_SSI<IT, ET>;

	using GP3     = GenericPoint3T<IT, ET>;
	// below points inherit from GP3
	using EP3     = ExplicitPoint3T<IT, ET>;
	using IP3_SSI = ImplicitPoint3T_SSI<IT, ET>;
	using IP3_LPI = ImplicitPoint3T_LPI<IT, ET>;
	using IP3_TPI = ImplicitPoint3T_TPI<IT, ET>;
	using IP3_LNC = ImplicitPoint3T_LNC<IT, ET>;

public:
	// clang-format off
	GP2       &operator()(      EP2& src) { return *static_cast<GP2 *>(&src); }
	const GP2 &operator()(const EP2& src) { return *static_cast<const GP2 *>(&src); }

	GP2       &operator()(      IP2_SSI& src) { return *static_cast<GP2 *>(&src); }
	const GP2 &operator()(const IP2_SSI& src) { return *static_cast<const GP2 *>(&src); }

	GP3       &operator()(      EP3& src) { return *static_cast<GP3 *>(&src); }
	const GP3 &operator()(const EP3& src) { return *static_cast<const GP3 *>(&src); }

	GP3       &operator()(      IP3_SSI& src) { return *static_cast<GP3 *>(&src); }
	const GP3 &operator()(const IP3_SSI& src) { return *static_cast<const GP3 *>(&src); }

	GP3       &operator()(      IP3_LPI& src) { return *static_cast<GP3 *>(&src); }
	const GP3 &operator()(const IP3_LPI& src) { return *static_cast<const GP3 *>(&src); }

	GP3       &operator()(      IP3_TPI& src) { return *static_cast<GP3 *>(&src); }
	const GP3 &operator()(const IP3_TPI& src) { return *static_cast<const GP3 *>(&src); }

	GP3       &operator()(      IP3_LNC& src) { return *static_cast<GP3 *>(&src); }
	const GP3 &operator()(const IP3_LNC& src) { return *static_cast<const GP3 *>(&src); }
	// clang-format on
};

template <typename IT, typename ET>
class AsExplicitPoint_Im
{
public:
	using GP2 = GenericPoint2T<IT, ET>;
	// below points inherit from GP2
	using EP2 = ExplicitPoint2T<IT, ET>;

	using GP3 = GenericPoint3T<IT, ET>;
	// below points inherit from GP3
	using EP3 = ExplicitPoint3T<IT, ET>;

public:
	EP2       &operator()(GP2 &src) { return src.EXP(); }
	EP3       &operator()(GP3 &src) { return src.EXP(); }
	const EP2 &operator()(const GP2 &src) { return src.EXP(); }
	const EP3 &operator()(const GP3 &src) { return src.EXP(); }

	bool is_explicit(const GP2 &src) { return src.is_Explicit(); }
	bool is_explicit(const GP3 &src) { return src.is_Explicit(); }
};

template <typename IT, typename ET>
class ToExplicitPoint_Im
{
public:
	using GP2     = GenericPoint2T<IT, ET>;
	// below points inherit from GP2
	using EP2     = ExplicitPoint2T<IT, ET>;
	using IP2_SSI = ImplicitPoint2T_SSI<IT, ET>;

	using GP3     = GenericPoint3T<IT, ET>;
	// below points inherit from GP3
	using EP3     = ExplicitPoint3T<IT, ET>;
	using IP3_SSI = ImplicitPoint3T_SSI<IT, ET>;
	using IP3_LPI = ImplicitPoint3T_LPI<IT, ET>;
	using IP3_TPI = ImplicitPoint3T_TPI<IT, ET>;
	using IP3_LNC = ImplicitPoint3T_LNC<IT, ET>;

public:
	EP2 operator()(const GP2 &src) { return src.to_Explicit(); }
	EP3 operator()(const GP3 &src) { return src.to_Explicit(); }

	EP2 operator()(const IP2_SSI &src) { return src.to_Explicit(); }
	EP3 operator()(const IP3_SSI &src) { return src.to_Explicit(); }
	EP3 operator()(const IP3_LPI &src) { return src.to_Explicit(); }
	EP3 operator()(const IP3_TPI &src) { return src.to_Explicit(); }
	EP3 operator()(const IP3_LNC &src) { return src.to_Explicit(); }
};

template <typename IT, typename ET>
class CreateImplicitSSI2_Im
{
public:
	using GP2     = GenericPoint2T<IT, ET>;
	// below points inherit from GP2
	using EP2     = ExplicitPoint2T<IT, ET>;
	using IP2_SSI = ImplicitPoint2T_SSI<IT, ET>;

	IP2_SSI operator()(const GP2 &l11, const GP2 &l12, const GP2 &l21,
	                   const GP2 &l22)
	{
		OMC_EXPENSIVE_ASSERT((l11.is_Explicit() && l12.is_Explicit() &&
		                      l21.is_Explicit() && l22.is_Explicit()),
		                     "Must initialize implicit point by explicit points.");
		return IP2_SSI(l11.EXP(), l12.EXP(), l21.EXP(), l22.EXP());
	}

	IP2_SSI operator()(const EP2 &l11, const EP2 &l12, const EP2 &l21,
	                   const EP2 &l22)
	{
		return IP2_SSI(l11, l12, l21, l22);
	}
};

template <typename IT, typename ET>
class CreateImplicitSSI3_Im
{
public:
	using GP3     = GenericPoint3T<IT, ET>;
	// below points inherit from GP3
	using EP3     = ExplicitPoint3T<IT, ET>;
	using IP3_SSI = ImplicitPoint3T_SSI<IT, ET>;

public:
	IP3_SSI operator()(const GP3 &a, const GP3 &b, const GP3 &p, const GP3 &q,
	                   int plane)
	{
		OMC_EXPENSIVE_ASSERT(p.is_Explicit() && q.is_Explicit() &&
		                       a.is_Explicit() && b.is_Explicit(),
		                     "Must initialize implicit point by explicit points.");
		OMC_EXPENSIVE_ASSERT(plane >= 0 && plane <= 2, "wrong plane.");
		return IP3_SSI(a.EXP(), b.EXP(), p.EXP(), q.EXP(), plane);
	}

	IP3_SSI operator()(const EP3 &a, const EP3 &b, const EP3 &p, const EP3 &q,
	                   int plane)
	{
		OMC_EXPENSIVE_ASSERT(plane >= 0 && plane <= 2, "wrong plane.");
		return IP3_SSI(a, b, p, q, plane);
	}
};

template <typename IT, typename ET>
class CreateImplicitLPI_Im
{
public:
	using GP3     = GenericPoint3T<IT, ET>;
	// below points inherit from GP3
	using EP3     = ExplicitPoint3T<IT, ET>;
	using IP3_LPI = ImplicitPoint3T_LPI<IT, ET>;

public:
	IP3_LPI operator()(const GP3 &p, const GP3 &q, const GP3 &r, const GP3 &s,
	                   const GP3 &t)
	{
		OMC_EXPENSIVE_ASSERT(
		  (p.is_Explicit() && q.is_Explicit() && r.is_Explicit() &&
		   s.is_Explicit() && t.is_Explicit()),
		  "Must initialize implicit point by explicit points.");
		return IP3_LPI(p.EXP(), q.EXP(), r.EXP(), s.EXP(), t.EXP());
	}

	IP3_LPI operator()(const EP3 &p, const EP3 &q, const EP3 &r, const EP3 &s,
	                   const EP3 &t)
	{
		return IP3_LPI(p, q, r, s, t);
	}
};

template <typename IT, typename ET>
class CreateImplicitTPI_Im
{
public:
	using GP3     = GenericPoint3T<IT, ET>;
	// below points inherit from GP3
	using EP3     = ExplicitPoint3T<IT, ET>;
	using IP3_TPI = ImplicitPoint3T_TPI<IT, ET>;

	IP3_TPI operator()(const GP3 &v1, const GP3 &v2, const GP3 &v3, const GP3 &w1,
	                   const GP3 &w2, const GP3 &w3, const GP3 &u1, const GP3 &u2,
	                   const GP3 &u3)
	{
		OMC_EXPENSIVE_ASSERT(
		  (v1.is_Explicit() && v2.is_Explicit() && v3.is_Explicit() &&
		   w1.is_Explicit() && w2.is_Explicit() && w3.is_Explicit() &&
		   u1.is_Explicit() && u2.is_Explicit() && u3.is_Explicit()),
		  "Must initialize implicit point by explicit points.");
		return IP3_TPI(v1.EXP(), v2.EXP(), v3.EXP(), w1.EXP(), w2.EXP(), w3.EXP(),
		               u1.EXP(), u2.EXP(), u3.EXP());
	}

	IP3_TPI operator()(const EP3 &v1, const EP3 &v2, const EP3 &v3, const EP3 &w1,
	                   const EP3 &w2, const EP3 &w3, const EP3 &u1, const EP3 &u2,
	                   const EP3 &u3)
	{
		return IP3_TPI(v1, v2, v3, w1, w2, w3, u1, u2, u3);
	}
};

template <typename IT, typename ET>
class CreateImplicitLNC_Im
{
public:
	using GP3     = GenericPoint3T<IT, ET>;
	// below points inherit from GP3
	using EP3     = ExplicitPoint3T<IT, ET>;
	using IP3_LNC = ImplicitPoint3T_LNC<IT, ET>;

	IP3_LNC operator()(const GP3 &p, const GP3 &q, const double t)
	{
		OMC_EXPENSIVE_ASSERT(
		  (p.is_Explicit() && q.is_Explicit()),
		  "Must initialize implicit point by explicit points.");
		return IP3_LNC(p.EXP(), q.EXP(), t);
	}

	IP3_LNC operator()(const EP3 &p, const EP3 &q, const double t)
	{
		return IP3_LNC(p, q, t);
	}
};

} // namespace OMC