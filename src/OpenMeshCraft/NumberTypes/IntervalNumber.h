#pragma once

#include "FPU.h"
#include "NumberUtils.h"

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

#ifdef OMC_SSE2
	#include "immintrin.h"
#endif

#include <numeric>
#include <utility>

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#undef IN
#undef TRUE
#undef FALSE

namespace OMC {

/// @tparam Protected Indicates whether the FPU round mode is protected.
/// @note Protect is time-consuming!!! Do it outside calculation blocks.
template <typename Protected = std::false_type>
class IntervalNumber
{
public:
	using IN = IntervalNumber<Protected>;

	using Pair = std::pair<double, double>;

	using Protector = FPU_RoundingProtector<std::false_type>;

	union CastedDouble
	{
		double   d;
		uint64_t u;

		inline CastedDouble() {}
		inline CastedDouble(double a)
		  : d(a)
		{
		}
		inline uint64_t is_negative() const { return u >> 63; }
	};

public: /* Constructors ***************************************************/
	IntervalNumber() = default;

	/// @note for efficiency, ni is the negative of real inf.
	explicit IntervalNumber(const double ni, const double s)
#ifdef OMC_SSE2
	  : m_v(_mm_setr_pd(ni, s))
#else
	  : m_ninf(ni)
	  , m_sup(s)
#endif
	{
		// OMC_EXPENSIVE_ASSERT(
		//   std::isfinite(ni) && std::isfinite(s) && (-ni) <= s, "invalid
		//   interval");
	}

#ifdef OMC_SSE2
	explicit IntervalNumber(__m128d v) noexcept
	  : m_v(v)
	{
	}
#endif

	IntervalNumber(const double d)
	  : IN(-d, d)
	{
	}

	IntervalNumber &operator=(const double d)
	{
		*this = IN(-d, d);
		return *this;
	}

	explicit IntervalNumber(const Pair &p)
	  : IN(-p.first, p.second)
	{
	}

public: /* Pre-defined ***************************************************/
	static IN min() // [min, min]
	{
		return IN(std::numeric_limits<double>::min());
	}

	static IN max() // [max, max]
	{
		return IN(std::numeric_limits<double>::max());
	}

	static IN infinity() // [inf, inf]
	{
		return IN(std::numeric_limits<double>::infinity());
	}

	static IN lowest() // [lowest, lowest]
	{
		return IN(std::numeric_limits<double>::lowest());
	}

	static IN largest() // [-inf, inf]
	{
		return IN(std::numeric_limits<double>::infinity(),
		          std::numeric_limits<double>::infinity());
	}

	static IN smallest() // [-min, min]
	{
		return IN(std::numeric_limits<double>::min(),
		          std::numeric_limits<double>::min());
	}

public: /* Access ********************************************************/
#ifdef OMC_SSE2
	double inf() const { return -_mm_cvtsd_f64(m_v); }
	double sup() const { return _mm_cvtsd_f64(_mm_shuffle_pd(m_v, m_v, 1)); }
#else
	double inf() const { return -m_ninf; }
	double sup() const { return m_sup; }
#endif

	Pair pair() const { return Pair(inf(), sup()); }

#ifdef OMC_SSE2
	__m128d simd() const { return m_v; }
#endif

public: /* Query *********************************************************/
	bool is_point() const { return inf() == sup(); }

	bool is_negative() const { return sup() < 0.; }

	bool is_positive() const { return inf() > 0.; }

	bool is_sign_reliable() const { return is_negative() || is_positive(); }

	Sign sign() const
	{
		return is_sign_reliable()
		         ? (is_negative() ? Sign::NEGATIVE : Sign::POSITIVE)
		         : Sign::UNCERTAIN;
	}

	bool is_nan() const { return inf() != inf(); }

#ifdef OMC_SSE2
	bool is_same(const IN &rhs) const
	{
		return _mm_movemask_pd(_mm_cmpeq_pd(m_v, rhs.m_v)) == 3;
	}

	bool do_overlap(const IN &rhs) const
	{
		__m128d m = _mm_set1_pd(-0.);
		__m128d y = _mm_xor_pd((-rhs).m_v, m); // {-ds,di}
		__m128d c = _mm_cmplt_pd(m_v, y);      // {i>ds,s<di}
		return _mm_movemask_pd(c) == 0;
	}
#else
	bool is_same(const IN &rhs) const
	{
		return m_ninf == rhs.m_ninf && m_sup == rhs.m_sup;
	}

	bool do_overlap(const IN &rhs) const
	{
		return !(rhs.inf() > sup() || rhs.sup() < inf());
	}
#endif

public: /* Modifications *****************************************************/
#ifdef OMC_SSE2
	void invert() { m_v = _mm_shuffle_pd(m_v, m_v, 1); }
#else
	void invert() { std::swap(m_ninf, m_sup); }
#endif

public: /* Binary Operators **************************************************/
#ifdef OMC_SSE2
	IN operator+() const { return *this; }
	IN operator-() const { return IN(_mm_shuffle_pd(m_v, m_v, 1)); }

	friend inline IN operator+(const IN &lhs, const IN &rhs)
	{
		InternalProtector _;
		return IN(_mm_add_pd(lhs.m_v, rhs.m_v));
	}
	friend inline IN operator-(const IN &lhs, const IN &rhs)
	{
		InternalProtector _;
		return IN(_mm_add_pd(lhs.m_v, (-rhs).m_v));
	}

	friend inline IN operator*(const IN &lhs, const IN &rhs)
	{
		InternalProtector _;

		// Brutal, compute all products in all directions.
		// The actual winner (by a hair) on recent hardware before removing NaNs.
		__m128d aa = lhs.m_v;                   // {-ai,as}
		__m128d bb = rhs.m_v;                   // {-bi,bs}
		__m128d m  = _mm_set_sd(-0.);           // {-0,+0}
		__m128d m1 = _mm_set1_pd(-0.);          // {-0,-0}
		__m128d ax = _mm_shuffle_pd(aa, aa, 1); // {as,-ai}
		__m128d ap = _mm_xor_pd(ax, m1);        // {-as,ai}
		__m128d bz = _mm_xor_pd(bb, m);         // {bi,bs}
		__m128d c  = _mm_shuffle_pd(bz, bz, 1); // {bs,bi}

		// The multiplications could produce some NaN, with 0 * inf. Replacing it
		// with inf is safe. min(x,y) (the order is essential) returns its second
		// argument when the first is NaN. An IEEE 754-2019 maximum could help.
		__m128d big = largest().m_v;
		__m128d x1  = _mm_mul_pd(aa, bz); // {-ai*bi,as*bs}
		// x1 = _mm_min_pd(x1,big); // no NaN
		__m128d x2  = _mm_mul_pd(aa, c);   // {-ai*bs,as*bi}
		x2          = _mm_min_pd(x2, big); // no NaN
		__m128d x3  = _mm_mul_pd(ap, bz);  // {-as*bi,ai*bs}
		// x3 = _mm_min_pd(x3,big); // no NaN
		__m128d x4  = _mm_mul_pd(ap, c);   // {-as*bs,ai*bi}
		x4          = _mm_min_pd(x4, big); // no NaN

		__m128d y1 = _mm_max_pd(x1, x2);
		__m128d y2 = _mm_max_pd(x3, x4);
		__m128d r  = _mm_max_pd(y1, y2);
		// Alternative with fewer instructions but more dependency
		// __m128d r =
		// _mm_max_pd(x1,_mm_max_pd(x2,_mm_max_pd(x3,_mm_min_pd(x4,big))));
		return IN(r);
	}
	friend inline IN operator*(double lhs, const IN &rhs)
	{
		OMC_ASSERT(std::isfinite(lhs), "invalid number");
		InternalProtector _;

		__m128d aa, bb;
		// clang-format off
		if (lhs < 0) { aa = _mm_set1_pd(-lhs); bb = rhs.m_v; }
		else         { aa = _mm_set1_pd(lhs);  bb = rhs.m_v; }
		// clang-format on

		__m128d r = _mm_mul_pd(aa, bb);
		r         = _mm_min_pd(r, largest().m_v);
		return IN(r);
	}
	friend inline IN operator*(const IN &lhs, double rhs) { return rhs * lhs; }

	friend inline IN operator/(const IN &lhs, const IN &rhs)
	{
		InternalProtector _;
		OMC_ASSERT((rhs != 0.) == Certainty::TRUE, "divide zero");
	#if defined(OMC_AVX)
		// Current fastest
		// if b>0 we want [ai/(ai>0?bs:bi),as/(as>0?bi:bs)]
		// if b<0 we want [as/(as>0?bs:bi),ai/(ai>0?bi:bs)]
		__m128d m  = _mm_set_sd(-0.);
		__m128d aa = lhs.m_v;
		__m128d bb = rhs.m_v;
		int     i  = _mm_movemask_pd(_mm_cmpge_pd(bb, _mm_set1_pd(0.)));
		if (i == 3)
			return largest();                      // bi<=0 && bs>=0
		__m128d ap  = _mm_xor_pd(aa, m);         // {ai, as}
		__m128d ax  = _mm_shuffle_pd(ap, ap, 1); // {as, ai}
		__m128d bp  = _mm_xor_pd(bb, m);         // {bi, bs}
		__m128d bx  = _mm_shuffle_pd(bp, bp, 1); // {bs, bi}
		__m128d num = _mm_blendv_pd(ap, ax, bp); // {(b>0)?ai:as, (b>0)?as:ai}
		__m128d d   = _mm_blendv_pd(bx, bp, num);
		// Can we rearrange things so we need fewer xor?
		__m128d den = _mm_xor_pd(d, m);
		__m128d r   = _mm_div_pd(num, den);
		return IN(r);
	#elif defined(OMC_SSE41)
		// Similar to the multiplication, but slow, because divisions are slow
		// if b>0 we want [-max(-ai/bi,-ai/bs),max(as/bi,as/bs)] {-ai,as}/{bi,bs}
		// {-ai,as}/{bs,bi} if b<0 we want [-max(-as/bi,-as/bs),max(ai/bi,ai/bs)]
		// {-as,ai}/{bi,bs} {-as,ai}/{bs,bi}
		__m128d m  = _mm_set_sd(-0.);
		__m128d m1 = _mm_set1_pd(-0.);
		__m128d aa = lhs.m_v; // {-ai, as}
		__m128d bb = rhs.m_v; // {-bi, bs}
		int     i  = _mm_movemask_pd(_mm_cmpge_pd(bb, _mm_set1_pd(0.)));
		if (i == 3)
			return largest();                      // bi<=0 && bs>=0
		__m128d ap  = _mm_xor_pd(aa, m1);        // {ai, -as}
		__m128d ax  = _mm_shuffle_pd(ap, ap, 1); // {-as, ai}
		__m128d bp  = _mm_xor_pd(bb, m);         // {bi, bs}
		__m128d num = _mm_blendv_pd(aa, ax, bp);
		__m128d bx  = _mm_shuffle_pd(bp, bp, 1); // {bs, bi}
		__m128d d1  = _mm_div_pd(num, bp);
		__m128d d2  = _mm_div_pd(num, bx);
		__m128d r   = _mm_max_pd(d1, d2);
		return IN(r);
	#else
		OMC_ASSERT((rhs != 0.) == Certainty::TRUE, "divide zero");
		InternalProtector _;
		CastedDouble l1(lhs.inf()), h1(lhs.sup()), l2(rhs.inf()), h2(rhs.sup());
		uint64_t     cfg = (l1.is_negative() << 3) + (h1.is_negative() << 2) +
		               (l2.is_negative() << 1) + (h2.is_negative());

		switch (cfg)
		// sign of [lhs.inf(), lhs.sup(), rhs.inf(), rhs.sup()]
		//          1-0+       1-0+       1-0+       1-0+
		{
		case 0: // [0000] -> [++++]
			return IN((-lhs.inf()) / rhs.sup(), lhs.sup() / rhs.inf());
		case 8: // [1000] -> [-+++]
			return IN((-lhs.inf()) / rhs.inf(), lhs.sup() / rhs.inf());
		case 12: // [1100] -> [--++]
			return IN((-lhs.inf()) / rhs.inf(), lhs.sup() / rhs.sup());
		case 3: // [0011] -> [++--]
			return IN(lhs.sup() / (-rhs.sup()), lhs.inf() / rhs.inf());
		case 11: // [1011] -> [-+--]
			return IN(lhs.sup() / (-rhs.sup()), lhs.inf() / rhs.sup());
		case 15: // [1111] -> [----]
			return IN(lhs.sup() / (-rhs.inf()), lhs.inf() / rhs.sup());
		case 2: // [0010] -> [++-+]
			OMC_FALLTHROUGH;
		case 10: // [1010] -> [-+-+]
			OMC_FALLTHROUGH;
		case 14: // [1110] -> [---+]
			return largest();
			// rhs overlaps zero
		};
		OMC_ASSERT(false, "");
	#endif
	}
	friend inline IN operator/(double lhs, const IN &rhs)
	{
		OMC_ASSERT(std::isfinite(lhs), "invalid number");
		OMC_ASSERT((rhs != 0.) == Certainty::TRUE, "divide zero");
		InternalProtector _;
		int i = _mm_movemask_pd(_mm_cmpge_pd(rhs.m_v, _mm_set1_pd(0.)));
		if (i == 3)
			return largest(); // bi<=0 && bs>=0
		__m128d aa, xx;
		if (lhs > 0)
		{
			aa = _mm_set1_pd(-lhs);
			xx = (-rhs).m_v;
		}
		else if (lhs < 0)
		{
			aa = _mm_set1_pd(lhs);
			xx = rhs.simd();
		}
		else
			return IN(0.);

		return IN(_mm_div_pd(aa, xx));
	}
	friend inline IN operator/(const IN &lhs, double rhs)
	{
		OMC_ASSERT(std::isfinite(rhs), "invalid number");
		OMC_ASSERT(rhs != 0., "divide zero");
		Protector _;
		__m128d   aa, bb;
		if (rhs < 0)
		{
			aa = (-lhs).m_v;
			bb = _mm_set1_pd(-rhs);
		}
		else if (rhs == 0)
			return largest();
		// Now rhs > 0
		return IN(_mm_div_pd(aa, bb));
	}
#else
	IN operator+() const { return *this; }
	IN operator-() const { return IN(m_sup, m_ninf); }

	friend inline IN operator+(const IN &lhs, const IN &rhs)
	{
		InternalProtector _;
		return IN(lhs.m_ninf + rhs.m_ninf, lhs.m_sup + rhs.m_sup);
	}
	friend inline IN operator-(const IN &lhs, const IN &rhs)
	{
		InternalProtector _;
		return IN(rhs.m_sup + lhs.m_ninf, lhs.m_sup + rhs.m_ninf);
	}

	friend inline IN operator*(const IN &lhs, const IN &rhs)
	{
		InternalProtector _;

		CastedDouble l1(lhs.m_ninf), h1(lhs.m_sup), l2(rhs.m_ninf), h2(rhs.m_sup);
		uint64_t     cfg = (l1.is_negative() << 3) + (h1.is_negative() << 2) +
		               (l2.is_negative() << 1) + (h2.is_negative());

		switch (cfg)
		// sign of [lhs.inf(), lhs.sup(), rhs.inf(), rhs.sup()]
		//          1+0-       1-0+       1+0-       1-0+
		{
		case 0: // [0000] -> [-+-+]
			double ll, lh, hl, hh;
			ll = lhs.m_ninf * rhs.m_ninf;
			lh = (lhs.m_ninf * rhs.m_sup);
			hl = (lhs.m_sup * rhs.m_ninf);
			hh = lhs.m_sup * rhs.m_sup;
			if (hl > lh)
				lh = hl;
			if (ll > hh)
				hh = ll;
			return IN(lh, hh);
		case 1: // [0001] -> [-+--]
			return IN(lhs.m_sup * rhs.m_ninf, lhs.m_ninf * rhs.m_ninf);
		case 2: // [0010] -> [-+++]
			return IN(lhs.m_ninf * rhs.m_sup, lhs.m_sup * rhs.m_sup);
		case 4: // [0100] -> [---+]
			return IN(lhs.m_ninf * rhs.m_sup, lhs.m_ninf * rhs.m_ninf);
		case 5: // [0101] -> [----]
			return IN((-lhs.m_sup) * rhs.m_sup, lhs.m_ninf * rhs.m_ninf);
		case 6: // [0110] -> [--++]
			return IN(lhs.m_ninf * rhs.m_sup, lhs.m_sup * (-rhs.m_ninf));
		case 8: // [1000] -> [++-+]
			return IN(lhs.m_sup * rhs.m_ninf, lhs.m_sup * rhs.m_sup);
		case 9: // [1001] -> [++--]
			return IN(lhs.m_sup * rhs.m_ninf, (-lhs.m_ninf) * rhs.m_sup);
		case 10: // [1010] -> [++++]
			return IN(lhs.m_ninf * (-rhs.m_ninf), lhs.m_sup * rhs.m_sup);
		};

		OMC_ASSERT(false, "");
	}
	friend inline IN operator*(double lhs, const IN &rhs)
	{
		OMC_ASSERT(std::isfinite(lhs), "invalid number");
		return IN(lhs * rhs.m_ninf, lhs * rhs.m_sup);
	}
	friend inline IN operator*(const IN &lhs, double rhs) { return rhs * lhs; }

	friend inline IN operator/(const IN &lhs, const IN &rhs)
	{
		OMC_ASSERT((rhs != 0.) == Certainty::TRUE, "divide zero");
		InternalProtector _;
		CastedDouble l1(lhs.m_ninf), h1(lhs.m_sup), l2(rhs.m_ninf), h2(rhs.m_sup);
		uint64_t     cfg = (l1.is_negative() << 3) + (h1.is_negative() << 2) +
		               (l2.is_negative() << 1) + (h2.is_negative());

		switch (cfg)
		// sign of [lhs.inf(), lhs.sup(), rhs.inf(), rhs.sup()]
		//          1+0-       1-0+       1+0-       1-0+
		{
		case 2: // [0010] -> [-+++]
			return IN(lhs.m_ninf / (-rhs.m_ninf), lhs.m_sup / (-rhs.m_ninf));
		case 6: // [0110] -> [--++]
			return IN(lhs.m_ninf / (-rhs.m_ninf), lhs.m_sup / rhs.m_sup);
		case 10: // [1010] -> [++++]
			return IN(lhs.m_ninf / rhs.m_sup, lhs.m_sup / (-rhs.m_ninf));
		case 1: // [0001] -> [-+--]
			return IN(lhs.m_sup / (-rhs.m_sup), (-lhs.m_ninf) / rhs.m_sup);
		case 5: // [0101] -> [----]
			return IN(lhs.m_sup / rhs.m_ninf, (-lhs.m_ninf) / rhs.m_sup);
		case 9: // [1001] -> [++--]
			return IN(lhs.m_sup / (-rhs.m_sup), (-lhs.m_ninf) / (-rhs.m_ninf));
		case 0: // [0000] -> [-+-+]
			OMC_FALLTHROUGH;
		case 4: // [0100] -> [---+]
			OMC_FALLTHROUGH;
		case 8:             // [1000] -> [++-+]
			return largest(); // rhs overlaps zero
		};
		OMC_ASSERT(false, "");
	}
	friend inline IN operator/(double lhs, const IN &rhs)
	{
		OMC_ASSERT(std::isfinite(lhs), "invalid number");
		OMC_ASSERT((rhs != 0.) == Certainty::TRUE, "divide zero");
		return IN(lhs) / rhs;
	}
	friend inline IN operator/(const IN &lhs, double rhs)
	{
		OMC_ASSERT(std::isfinite(rhs), "invalid number");
		OMC_ASSERT(rhs != 0., "divide zero");
		return lhs / IN(rhs);
	}
#endif

	IN &operator+=(const IN &rhs) { return *this = *this + rhs; }
	IN &operator-=(const IN &rhs) { return *this = *this - rhs; }
	IN &operator*=(const IN &rhs) { return *this = *this * rhs; }
	IN &operator/=(const IN &rhs) { return *this = *this / rhs; }

public: /* Comparisons */
	friend inline Certainty operator<(const IN &lhs, const IN &rhs)
	{
		if (lhs.sup() < rhs.inf())
			return Certainty::TRUE;
		if (lhs.inf() >= rhs.sup())
			return Certainty::FALSE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator>(const IN &lhs, const IN &rhs)
	{
		return rhs < lhs;
	}
	friend inline Certainty operator<=(const IN &lhs, const IN &rhs)
	{
		if (lhs.sup() <= rhs.inf())
			return Certainty::TRUE;
		if (lhs.inf() > rhs.sup())
			return Certainty::FALSE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator>=(const IN &lhs, const IN &rhs)
	{
		return rhs <= lhs;
	}

	friend inline Certainty operator==(const IN &lhs, const IN &rhs)
	{
		if (rhs.inf() > lhs.sup() || rhs.sup() < lhs.inf())
			return Certainty::FALSE;
		if (rhs.inf() == lhs.sup() && rhs.sup() == lhs.inf())
			return Certainty::TRUE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator!=(const IN &lhs, const IN &rhs)
	{
		return !(lhs == rhs);
	}

	friend inline Certainty operator<(double lhs, const IN &rhs)
	{
		if (lhs < rhs.inf())
			return Certainty::TRUE;
		if (lhs >= rhs.sup())
			return Certainty::FALSE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator>(double lhs, const IN &rhs)
	{
		return rhs < lhs;
	}
	friend inline Certainty operator<=(double lhs, const IN &rhs)
	{
		if (lhs <= rhs.inf())
			return Certainty::TRUE;
		if (lhs > rhs.sup())
			return Certainty::FALSE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator>=(double lhs, const IN &rhs)
	{
		return rhs <= lhs;
	}

	friend inline Certainty operator==(double lhs, const IN &rhs)
	{
		if (rhs.inf() > lhs || rhs.sup() < lhs)
			return Certainty::FALSE;
		if (rhs.is_point())
			return Certainty::TRUE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator!=(double lhs, const IN &rhs)
	{
		if (rhs.inf() > lhs || rhs.sup() < lhs)
			return Certainty::TRUE;
		if (rhs.is_point())
			return Certainty::FALSE;
		return Certainty::UNCERTAIN;
	}

	friend inline Certainty operator<(const IN &lhs, double rhs)
	{
		if (lhs.sup() < rhs)
			return Certainty::TRUE;
		if (lhs.inf() >= rhs)
			return Certainty::FALSE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator>(const IN &lhs, double rhs)
	{
		return rhs < lhs;
	}

	friend inline Certainty operator<=(const IN &lhs, double rhs)
	{
		if (lhs.sup() <= rhs)
			return Certainty::TRUE;
		if (lhs.inf() > rhs)
			return Certainty::FALSE;
		return Certainty::UNCERTAIN;
	}
	friend inline Certainty operator>=(const IN &lhs, double rhs)
	{
		return rhs <= lhs;
	}

	friend inline Certainty operator==(const IN &lhs, double rhs)
	{
		return rhs == lhs;
	}
	friend inline Certainty operator!=(const IN &lhs, double rhs)
	{
		return rhs != lhs;
	}

public: /* Unary Operators ***************************************************/
	IN abs() const
	{
		if (inf() >= 0.0)
			return *this;
		if (sup() <= 0.0)
			return -(*this);
		return IN(-0.0, std::max(-inf(), sup()));
	}

private:
	using InternalProtector = FPU_RoundingProtector<Protected>;

private:
#ifdef OMC_SSE2
	// (-inf, sup)
	// dst[63:0] := -inf, dst [127:64] : = sup
	__m128d m_v;
#else
	// (-inf, sup)
	double m_ninf, m_sup;
#endif
};

template <typename Protected>
class UnaryOperators<IntervalNumber<Protected>>
{
public:
	using NT = IntervalNumber<Protected>;

public:
	static Sign sign(const NT &x)
	{
		if ((x > 0) == Certainty::TRUE)
			return Sign::POSITIVE;
		else if ((x < 0) == Certainty::TRUE)
			return Sign::NEGATIVE;
		else
			return Sign::UNCERTAIN;
	}

	static NT abs(const NT &x) { return x.abs(); }

	static NT negate(const NT &x) { return -x; }
};

} // namespace OMC