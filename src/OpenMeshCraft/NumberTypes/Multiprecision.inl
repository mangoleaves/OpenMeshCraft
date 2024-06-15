#pragma once

#include "Multiprecision.h"

namespace OMC {

void fast_two_sum(const double a, const double b, double &x1, double &x0)
{
	x1 = a + b;
	x0 = b - (x1 - a);
}

void two_sum(const double a, const double b, double &x1, double &x0)
{
	double _av, _bv, _ar, _br;
	x1  = a + b;
	_bv = x1 - a;
	_av = x1 - _bv;
	_br = b - _bv;
	_ar = a - _av;
	x0  = _ar + _br;
}

void two_one_sum_full(const double a1, const double a0, const double b,
                      double &x2, double &x1, double &x0)
{
	double _i;
	// _i = RN(a0 + b), x0 = a0 + b - RN(a0 + b).
	two_sum(a0, b, _i, x0);
	// x2 = RN(a1 + RN(a0 + b)), x1 = a1 + RN(a0 + b) - RN(a1 + RN(a0 + b)).
	two_sum(a1, _i, x2, x1);
}

void two_one_sum_clip(const double a1, const double a0, const double b,
                      double &x1, double &x0)
{
	double _t0, _t1;
	two_sum(a1, b, _t1, _t0);
	_t0 = _t0 + a0;
	fast_two_sum(_t1, _t0, x1, x0);
}

void two_two_sum_full(const double a1, const double a0, const double b1,
                      const double b0, double &x3, double &x2, double &x1,
                      double &x0)
{
	double _j, _0;
	// _j = RN(a1 + RN(a0 + b0))
	// _0 = a1 + RN(a0 + b0) - RN(a1 + RN(a0 + b0))
	// x0 = a0 + b0 - RN(a0 + b0)
	two_one_sum_full(a1, a0, b0, _j, _0, x0);
	// x3 = RN(_j + RN(_0 + b1))
	// x2 = _j + RN(_0 + b1) - RN(_j + RN(_0 + b1))
	// x1 = _0 + b1 - RN(_0 + b1)
	two_one_sum_full(_j, _0, b1, x3, x2, x1);
}

void two_two_sum_clip(const double a1, const double a0, const double b1,
                      const double b0, double &x1, double &x0)
{
	double _t1, _t0, _s1, _s0, _w1, _w0;
	two_sum(a1, b1, _s1, _s0);
	two_sum(a0, b0, _t1, _t0);
	_s0 = _s0 + _t1;
	fast_two_sum(_s1, _s0, _w1, _w0);
	_w0 = _w0 + _t0;
	fast_two_sum(_w1, _w0, x1, x0);
}

void fast_two_sub(const double a, const double b, double &x1, double &x0)
{
	x1 = a - b;
	x0 = (a - x1) - b;
}

void two_sub(const double a, const double b, double &x1, double &x0)
{
	double _bv, _av, _br, _ar;
	x1  = a - b;
	_bv = a - x1;
	_av = x1 + _bv;
	_br = _bv - b;
	_ar = a - _av;
	x0  = _ar + _br;
}

void two_one_sub_full(const double a1, const double a0, const double b,
                      double &x2, double &x1, double &x0)
{
	double _tt;
	two_sub(a0, b, _tt, x0);
	two_sub(a1, _tt, x2, x1);
}

void two_one_sub_clip(const double a1, const double a0, const double b,
                      double &x1, double &x0)
{
	double _t0, _t1;
	two_sub(a1, b, _t1, _t0);
	_t0 = _t0 + a0;
	two_sub(_t1, _t0, x1, x0);
}

void two_two_sub_full(const double a1, const double a0, const double b1,
                      const double b0, double &x3, double &x2, double &x1,
                      double &x0)
{
	double _t1, _t0;
	two_one_sub_full(a1, a0, b0, _t1, _t0, x0);
	two_one_sub_full(_t1, _t0, b1, x3, x2, x1);
}

void two_two_sub_clip(const double a1, const double a0, const double b1,
                      const double b0, double &x1, double &x0)
{
	double _s0, _s1, _t0, _t1, _w0, _w1;
	two_sub(a1, b1, _s1, _s0);
	two_sub(a0, b0, _t1, _t0);
	_s0 = _s0 + _t1;
	fast_two_sub(_s1, _s0, _w1, _w0);
	_w0 = _w0 + _t0;
	fast_two_sub(_w1, _w0, x1, x0);
}

void one_split(double a, double &ah, double &al)
{
	double _c = 1.3421772800000003e+008 * a; // constant is 2^27
	ah        = _c - (_c - a);
	al        = a - ah;
}

void two_prod(const double a, const double b, double &x1, double &x0)
{
#if defined(OMC_AVX2)
	// I am so confused by when std::fma works across different platforms and
	// compilers, so, it may be more straightforward to directly use FMA
	// intrinsics from AVX2.
	__m128d _x1, _x0, _a, _b;
	_a  = _mm_set_sd(a);
	_b  = _mm_set_sd(b);
	_x1 = _mm_mul_sd(_a, _b);
	_x0 = _mm_fmsub_sd(_a, _b, _x1);
	x1  = _mm_cvtsd_f64(_x1);
	x0  = _mm_cvtsd_f64(_x0);
#else
	double _ah, _al, _bh, _bl;
	x1 = a * b;
	one_split(a, _ah, _al);
	one_split(b, _bh, _bl);
	x0 = ((_ah * _bh - x1) + _ah * _bl + _al * _bh) + _al * _bl;
#endif
}

void two_prod_presplit(double a, double b, double bh, double bl, double &x1,
                       double &x0)
{
	double _ah, _al;
	x1 = a * b;
	one_split(a, _ah, _al);
	x0 = ((_ah * bh - x1) + _ah * bl + _al * bh) + _al * bl;
}

void two_prod_2presplit(double a, double ah, double al, double b, double bh,
                        double bl, double &x1, double &x0)
{
	x1 = a * b;
	x0 = ((ah * bh - x1) + ah * bl + al * bh) + al * bl;
}

void one_square(const double a, double &x1, double &x0)
{
#if defined(OMC_AVX2)
	__m128d _x1, _x0, _a;
	_a  = _mm_set_sd(a);
	_x1 = _mm_mul_sd(_a, _a);
	_x0 = _mm_fmsub_sd(_a, _a, _x1);
	x1  = _mm_cvtsd_f64(_x1);
	x0  = _mm_cvtsd_f64(_x0);
#else
	double _ah, _al;
	x1 = a * a;
	one_split(a, _ah, _al);
	x0 = ((_ah * _ah - x1) + _ah * _al + _al * _ah) + _al * _al;
#endif
}

void one_square_presplit(const double a, const double ah, const double al,
                         double &x1, double &x0)
{
	x1 = a * a;
	x0 = ((ah * ah - x1) + ah * al + al * ah) + al * al;
}

void two_one_prod_full(const double a1, const double a0, const double b,
                       double &x3, double &x2, double &x1, double &x0)
{
	double _i, _j, _0, _k;
#if defined(OMC_AVX2)
	two_prod(a0, b, _i, x0);
	two_prod(a1, b, _j, _0);
#else
	double _bh, _bl;
	one_split(b, _bh, _bl);
	two_prod_presplit(a0, b, _bh, _bl, _i, x0);
	two_prod_presplit(a1, b, _bh, _bl, _j, _0);
#endif
	two_sum(_i, _0, _k, x1);
	fast_two_sum(_j, _k, x3, x2);
}

void two_one_prod_clip(const double a1, const double a0, const double b,
                       double &x1, double &x0)
{
#if defined(OMC_AVX2)
	__m128d _a1, _a0, _b, _x1, _x0, _t0, _t1;
	// store into __m128d
	_a1 = _mm_set_sd(a1);
	_a0 = _mm_set_sd(a0);
	_b  = _mm_set_sd(b);
	// expand two_prod(a1, b, _t1, _t0) here.
	_t1 = _mm_mul_sd(_a1, _b);
	_t0 = _mm_fmsub_sd(_a1, _b, _t1);
	// expand fast_two_sum(_t1, _t0, _x1, _x0) here.
	_x1 = _mm_add_sd(_t1, _t0);
	_x0 = _mm_sub_sd(_t0, _mm_sub_sd(_x1, _t1));
	// extract from __m128d
	x1  = _mm_cvtsd_f64(_x1);
	x0  = _mm_cvtsd_f64(_x0);
#else
	double _t0, _t1, _ss;
	two_prod(a1, b, _t1, _t0);
	_ss = a0 * b;
	_t0 = _t0 + _ss;
	fast_two_sum(_t1, _t0, x1, x0);
#endif
}

void two_two_prod_full(const double &a1, const double &a0, const double &b1,
                       const double &b0, double *x)
{
	double _0, _1, _2;
	double _i, _j, _k, _l, _m, _n;
#if defined(OMC_AVX2)
	two_prod(a0, b0, _i, x[0]);
	two_prod(a1, b0, _j, _0);
	two_sum(_i, _0, _k, _1);
	fast_two_sum(_j, _k, _l, _2);
	two_prod(a0, b1, _i, _0);
	two_sum(_1, _0, _k, x[1]);
	two_sum(_2, _k, _j, _1);
	two_sum(_l, _j, _m, _2);
	two_prod(a1, b1, _j, _0);
#else
	double a0h, a0l, a1h, a1l, b0h, b0l, b1h, b1l;
	one_split(a0, a0h, a0l);
	one_split(b0, b0h, b0l);
	one_split(a1, a1h, a1l);
	one_split(b1, b1h, b1l);
	two_prod_2presplit(a0, a0h, a0l, b0, b0h, b0l, _i, x[0]);
	two_prod_2presplit(a1, a1h, a1l, b0, b0h, b0l, _j, _0);
	two_sum(_i, _0, _k, _1);
	fast_two_sum(_j, _k, _l, _2);
	two_prod_2presplit(a0, a0h, a0l, b1, b1h, b1l, _i, _0);
	two_sum(_1, _0, _k, x[1]);
	two_sum(_2, _k, _j, _1);
	two_sum(_l, _j, _m, _2);
	two_prod_2presplit(a1, a1h, a1l, b1, b1h, b1l, _j, _0);
#endif
	two_sum(_i, _0, _n, _0);
	two_sum(_1, _0, _i, x[2]);
	two_sum(_2, _i, _k, _1);
	two_sum(_m, _k, _l, _2);
	two_sum(_j, _n, _k, _0);
	two_sum(_1, _0, _j, x[3]);
	two_sum(_2, _j, _i, _1);
	two_sum(_l, _i, _m, _2);
	two_sum(_1, _k, _i, x[4]);
	two_sum(_2, _i, _k, x[5]);
	two_sum(_m, _k, x[7], x[6]);
}

void two_two_prod_clip(const double &a1, const double &a0, const double &b1,
                       const double &b0, double &x1, double &x0)
{
#if defined(OMC_AVX2)
	__m128d _a1, _a0, _b1, _b0, _x1, _x0, ch, cl1, cl2, cl3, tl0, tl1;
	// store into __m128d
	_a1 = _mm_set_sd(a1);
	_a0 = _mm_set_sd(a0);
	_b1 = _mm_set_sd(b1);
	_b0 = _mm_set_sd(b0);
	// expand two_prod(a1, b1, ch, cl1) here.
	ch  = _mm_mul_sd(_a1, _b1);
	cl1 = _mm_fmsub_sd(_a1, _b1, ch);
	// remaining operations.
	tl0 = _mm_mul_sd(_a0, _b0);
	tl1 = _mm_fmadd_sd(_a1, _b0, tl0);
	cl2 = _mm_fmadd_sd(_a0, _b1, tl1);
	cl3 = _mm_add_sd(cl1, cl2);
	// expand fast_two_sum(ch, cl3) here.
	_x1 = _mm_add_sd(ch, cl3);
	_x0 = _mm_sub_sd(cl3, _mm_sub_sd(_x1, ch));
	// extract from __m128d
	x1  = _mm_cvtsd_f64(_x1);
	x0  = _mm_cvtsd_f64(_x0);
#else
	double ch, cl1, cl2, cl3, tl1, tl2;
	two_prod(a1, b1, ch, cl1);
	tl1 = a1 * b0;
	tl2 = a0 * b1;
	cl2 = tl1 + tl2;
	cl3 = cl1 + cl2;
	fast_two_sum(ch, cl3, x1, x0);
#endif
}

void two_square_full(const double &a1, const double &a0, double *x)
{
	double _0, _1, _2, _3;
	double _i, _j, _k, _l;
	one_square(a0, _i, x[0]);
	_0 = a0 + a0;
#if defined(OMC_AVX2)
	two_prod(_0, a1, _j, _1);
	one_square(a1, _l, _3);
#else
	double a1h, a1l;
	one_split(a1, a1h, a1l);
	two_prod_presplit(_0, a1, a1h, a1l, _j, _1);
	one_square_presplit(a1, a1h, a1l, _l, _3);
#endif
	two_one_sum_full(_j, _1, _i, _k, _2, x[1]);
	two_two_sum_full(_l, _3, _k, _2, x[5], x[4], x[3], x[2]);
}

void two_square_clip(const double &a1, const double &a0, double &x1, double &x0)
{
#if defined(OMC_AVX2)
	__m128d _a1, _a0, _x1, _x0, ch, cl1, cl2, cl3, _2;
	// store into __m128d
	_a1 = _mm_set_sd(a1);
	_a0 = _mm_set_sd(a0);
	_2  = _mm_set_sd(2.);
	// expand one_square(a1, ch, cl1) here.
	ch  = _mm_mul_sd(_a1, _a1);
	cl1 = _mm_fmsub_sd(_a1, _a1, ch);
	// remaining operations.
	cl2 = _mm_mul_sd(_mm_mul_sd(_a0, _a1), _2);
	cl3 = _mm_add_sd(cl1, cl2);
	// expand fast_two_sum(ch, cl3) here.
	_x1 = _mm_add_sd(ch, cl3);
	_x0 = _mm_sub_sd(cl3, _mm_sub_sd(_x1, ch));
	// extract from __m128d
	x1  = _mm_cvtsd_f64(_x1);
	x0  = _mm_cvtsd_f64(_x0);
#else
	double ch, cl1, cl2, cl3, tl1, tl2;
	one_square(a1, ch, cl1);
	cl2 = (a1 * a0) * 2.;
	cl3 = cl1 + cl2;
	fast_two_sum(ch, cl3, x1, x0);
#endif
}

} // namespace OMC