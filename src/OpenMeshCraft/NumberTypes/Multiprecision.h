#pragma once

#include "OpenMeshCraft/Utils/Macros.h"

#if defined(OMC_AVX2)
	#include "immintrin.h"
#endif

#include <cmath>

namespace OMC {

/**********************************************************************/
/*           Multiple Precision Floating-Point Arithtmetic            */
/**********************************************************************/

/*
 * Implementations of basic multi-precision arithmetic due to Dekker, Knuth,
 * Priest, Shewchuk, and others.
 *
 * Reference survay and papers for theories:
 *
 * - Boldo S, Jeannerod C-P, Melquiond G, Muller J-M. Floating-point arithmetic.
 * Acta Numerica. 2023;32:203-290. doi:10.1017/S0962492922000101
 * - M. Joldes, J-M. Muller, V. Popescu (2017): Tight & rigourous error bounds
 * for basic building blocks of double-word arithmetic. ACM Transactions on
 * Mathematical Software, ACM, 44 (2), pp. 1-27.
 *
 * Reference libraries for implementations:
 *
 * - geogram: https://github.com/BrunoLevy/geogram
 * - ggforce: https://ggforce.data-imaginist.com/
 * - IndirectPredicates: https://github.com/MarcoAttene/Indirect_Predicates
 *
 * Notes:
 *
 * 1. We only consider "double" precision floating-point number now.
 * 2. There are "full" and "clip" versions for some routines, where the "full"
 * one is used in expansion arithmetic and the "clip" one is used in
 * double-double arithmetic.
 *
 * Concepts:
 *
 * - Significand precision (p): 24 for float, 53 for double.
 * - Unit rounding error (u): 2^{-p}.
 * - Unit in last place (ulp).
 * - Unit in first place (ufp).
 * - Rounding functions: round a real number to one of two nearest representable
 * floating-point numbers.
 *   - RN : round to near
 * - Double-double number: a double-double number x[x1,x0] satisfies x = x1 + x0
 * and x1 = RN(x1 + x0), where x1 is the high component and x0 is the low
 * component of the double-double number.
 *
 * Supplements:
 *
 * - How to detect FMA?
 * https://stackoverflow.com/questions/16348909/how-do-i-know-if-i-can-compile-with-fma-instruction-sets
 */

/*****************************************************************************
 * Summary routines
 *****************************************************************************/

/**
 * @brief x + y = a + b. Require that a/2 <= b <= 2a.
 * @param[in] a summand
 * @param[in] b summand
 * @param[out] x1 x1 = RN(a + b)
 * @param[out] x0 x0 = (a + b) - RN(a + b)
 */
inline void fast_two_sum(const double a, const double b, double &x1,
                         double &x0);

/**
 * @brief x + y = a + b
 * @param[in] a summand
 * @param[in] b summand
 * @param[out] x1 x1 = RN(a + b)
 * @param[out] x0 x0 = a + b - RN(a + b)
 */
inline void two_sum(const double a, const double b, double &x1, double &x0);

/**
 * @brief x2 + x1 + x0 = a1 + a0 + b
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b summand
 * @param[out] x2 high component of x,
 * x2 = RN(a1 + RN(a0 + b))
 * @param[out] x1 x1 = see details inside function.
 * @param[out] x0 x0 = see details inside function.
 */
inline void two_one_sum_full(const double a1, const double a0, const double b,
                             double &x2, double &x1, double &x0);

/**
 * @brief add a[a1, a0] and b and store the result in x[x1, x0].
 * It is not error-free. The relative error is 2u^2, where u is the unit
 * rounding error.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b summand
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 */
inline void two_one_sum_clip(const double a1, const double a0, const double b,
                             double &x1, double &x0);

/**
 * @brief x3 + x2 + x1 + x0 = a1 + a0 + b1 + b0
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b1 high component of b
 * @param[in] b0 low component of b
 * @param[out] x3 high component of x,
 * x3 = RN(RN(a1+RN(a0+b0)) + RN(a1+RN(a0+b0) - RN(a1+RN(a0+b0))+b1))
 * @param[out] x2 x2 = see details inside function.
 * @param[out] x1 x1 = see details inside function.
 * @param[out] x0 x0 = see details inside function.
 */
inline void two_two_sum_full(const double a1, const double a0, const double b1,
                             const double b0, double &x3, double &x2,
                             double &x1, double &x0);

/**
 * @brief add a[a1, a0] and b[b1, b0] and store the result in x[x1, x0].
 * It is not error-free. The relative error is 3u^2 + 13u^3, where u is the unit
 * rounding error.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b1 high component of b
 * @param[in] b0 low component of b
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 */
inline void two_two_sum_clip(const double a1, const double a0, const double b1,
                             const double b0, double &x1, double &x0);

/*****************************************************************************
 * Subtract routines
 *****************************************************************************/

/**
 * @brief x + y = a - b. Require that a/2 <= -b <= 2a.
 * It is just the dual operation of fast_two_sum on a and -b.
 * @param[in] a minuend
 * @param[in] b subtrahend
 * @param[out] x1 x1 = RN(a - b)
 * @param[out] x0 x0 = (a - b) - RN(a - b)
 */
inline void fast_two_sub(const double a, const double b, double &x1,
                         double &x0);

/**
 * @brief x + y = a - b
 * @param[in] a minuend
 * @param[in] b subtrahend
 * @param[out] x1 x1 = RN(a - b)
 * @param[out] x0 x0 = a - b - RN(a - b)
 */
inline void two_sub(const double a, const double b, double &x1, double &x0);

/**
 * @brief x2 + x1 + x0 = a1 + a0 - b
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b subtrahend
 * @param[out] x2 x2 = see details inside function
 * @param[out] x1 x1 = see details inside function
 * @param[out] x0 x0 = see details inside function
 */
inline void two_one_sub_full(const double a1, const double a0, const double b,
                             double &x2, double &x1, double &x0);

/**
 * @brief subtract b from a[a1, a0] and store the result in x[x1, x0].
 * It is not error-free. The relative error is 2u^2, where u is the unit
 * rounding error.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b subtrahend
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 */
inline void two_one_sub_clip(const double a1, const double a0, const double b,
                             double &x1, double &x0);

/**
 * @brief x3 + x2 + x1 + x0 = (a1 + a0) - (b1 + b0)
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b1 high component of b
 * @param[in] b0 low component of b
 * @param[out] x3 x3 = see details inside function.
 * @param[out] x2 x2 = see details inside function.
 * @param[out] x1 x1 = see details inside function.
 * @param[out] x0 x0 = see details inside function.
 */
inline void two_two_sub_full(const double a1, const double a0, const double b1,
                             const double b0, double &x3, double &x2,
                             double &x1, double &x0);

/**
 * @brief subtract b[b1, b0] from a[a1, a0] and store the result in x[x1, x0].
 * It is not error-free. The relative error is 3u^2 + 13u^3, where u is the unit
 * rounding error.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b1 high component of b
 * @param[in] b0 low component of b
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 */
inline void two_two_sub_clip(const double a1, const double a0, const double b1,
                             const double b0, double &x1, double &x0);

/*****************************************************************************
 * Product routines
 *****************************************************************************/

/**
 * @brief split a to ah and al with smaller significand size, such that
 * a = ah + al and |al| <= ulp(ah)/2.
 * @param[in] a number to be split.
 * @param[out] ah high component
 * @param[out] al low component
 * @details By requiring ah and al to have floor(p/2)-bit (27 bit for double)
 * significands (which is always possible in base-2 FP arithmetic, even when p
 * is odd), and doing the same for another FP number b decomposed into bh and
 * bl, we ensure that xh*yh, xh*yl, xl*yh and xl&yl are FP numbers (and hence
 * are computed exactly with FP multiplications). This is the key to the
 * algorithm of Dekker (1971), used in place of 2Prod when no FMA instruction is
 * available.
 */
inline void one_split(double a, double &ah, double &al);

/**
 * @brief calculate the product of a and b.
 * @param[in] a multiplier
 * @param[in] b multiplier
 * @param[out] x1 high component of result
 * @param[out] x0 low component of result
 */
inline void two_prod(const double a, const double b, double &x1, double &x0);

/**
 * @brief calculate the product of a and b when b is already split.
 * @param[in] a multiplier
 * @param[in] b multiplier
 * @param[in] bh high component of b
 * @param[in] bl low component of b
 * @param[out] x1 high component of result
 * @param[out] x0 low component of result
 */
inline void two_prod_presplit(double a, double b, double bh, double bl,
                              double &x1, double &x0);

/**
 * @brief calculate the product of a and b when b is already split.
 * @param[in] a multiplier
 * @param[in] b multiplier
 * @param[in] ah high component of a
 * @param[in] al low component of a
 * @param[in] bh high component of b
 * @param[in] bl low component of b
 * @param[out] x1 high component of result
 * @param[out] x0 low component of result
 */
inline void two_prod_2presplit(double a, double ah, double al, double b,
                               double bh, double bl, double &x1, double &x0);

/**
 * @brief squares a number and store the results in x[x1,x0]
 * @param[in] a the number to square
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 */
inline void one_square(const double a, double &x1, double &x0);

/**
 * @brief squares a number and store the results in x[x1,x0]
 * @param[in] a the number to square
 * @param[in] ah high component of a
 * @param[in] al low component of a
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 */
inline void one_square_presplit(const double a, const double ah,
                                const double al, double &x1, double &x0);

/**
 * @brief calculate product of a[a1,a0] and b and store the result in an
 * expansion [x3,x2,x1,x0].
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b multiplier
 * @param[out] x3 high component of x
 * @param[out] x2 component of x
 * @param[out] x1 component of x
 * @param[out] x0 component of x
 */
inline void two_one_prod_full(const double a1, const double a0, const double b,
                              double &x3, double &x2, double &x1, double &x0);

/**
 * @brief calculate product of a[a1,a0] and b and store the result in a
 * double-double number [x1,x0].
 * It is not error-free. See relative error in details.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b multiplier
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 * @details The implementation without fma is DWTimesFP2 in Joldes (2017), which
 * has a relative error of 3u^2. The implementation with fma is DWTimesFP3 in
 * Joldes (2017), which has a relative error of 2u^2. u is the unit rounding
 * error.
 */
inline void two_one_prod_clip(const double a1, const double a0, const double b,
                              double &x1, double &x0);

/**
 * @brief Calculate product of two expansions of length 2.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b1 high component of b
 * @param[in] b0 low component of b
 * @param[out] x result expansion, with longest length of 8.
 */
inline void two_two_prod_full(const double &a1, const double &a0,
                              const double &b1, const double &b0, double *x);

/**
 * @brief Calculate product of two expansions of length 2 and store clipped
 * result in a double-double number x[x1,x0].
 * It is not error-free. See relative error in details.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[in] b1 high component of b
 * @param[in] b0 low component of b
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 * @details The implementation without fma is DWTimesDW1 in Joldes (2017), which
 * has a relative error of 7u^2. The implementation with fma is DWTimesDW3 in
 * Joldes (2017), which has a relative error of 5u^2. u is the unit rounding
 * error.
 */
inline void two_two_prod_clip(const double &a1, const double &a0,
                              const double &b1, const double &b0, double &x1,
                              double &x0);

/**
 * @brief Calculate square of an expansion of length 2.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[out] x result expansion, with longest length of 6.
 */
inline void two_square_full(const double &a1, const double &a0, double *x);

/**
 * @brief Calculate square of an expansion of length 2 and store clipped result
 * in a double-double number x[x1,x0].
 * It is not error-free. The relative error is [TODO], where u is the unit
 * rounding error.
 * @param[in] a1 high component of a
 * @param[in] a0 low component of a
 * @param[out] x1 high component of x
 * @param[out] x0 low component of x
 * @details The implementation with/without fma is a modification of DWTimesDW1
 * in Joldes (2017), which originally has a relative error of 7u^2. The actual
 * relative error should be smaller, but I havn't try to give a formally proved
 * bound.
 */
inline void two_square_clip(const double &a1, const double &a0, double &x1,
                            double &x0);

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Multiprecision.inl"
#endif