#pragma once

#include "Multiprecision.h"

namespace OMC {

/////////////////////////////////////////////////////////////////////
//
// 	   E X P A N S I O N   A R I T H M E T I C
//
/////////////////////////////////////////////////////////////////////

// Allocate extra-memory
#define AllocDoubles(n) ((double *)malloc((n) * sizeof(double)))
#define FreeDoubles(p) (free(p))

// An instance of the following must be created to access functions for
// expansion arithmetic
class expansionObject
{
public:
	/****************************************************************************
	 * Wrappers of basic expansion arithmetic,
	 * due to Dekker, Knuth, Priest, Shewchuk, and others.
	 ***************************************************************************/

	/* Sum routines ************************************************************/

	/// @brief Wrapper of two_sum.
	inline void Two_Sum(const double a, const double b, double *x)
	{
		two_sum(a, b, x[1], x[0]);
	}

	/// @brief Wrapper of two_one_sub_full.
	inline void Two_One_Sum(const double *a, const double b, double *x)
	{
		two_one_sum_full(a[1], a[0], b, x[2], x[1], x[0]);
	}

	/// @brief Wrapper of two_two_sum_full.
	inline void Two_Two_Sum(const double *a, const double *b, double *x)
	{
		two_two_sum_full(a[1], a[0], b[1], b[0], x[3], x[2], x[1], x[0]);
	}

	/* Difference routines *****************************************************/

	/// @brief Wrapper of two_sub.
	inline void Two_Diff(const double a, const double b, double *x)
	{
		two_sub(a, b, x[1], x[0]);
	}

	/// @brief Wrapper of two_one_sub_full.
	inline void Two_One_Diff(const double *a, const double b, double *x)
	{
		two_one_sub_full(a[1], a[0], b, x[2], x[1], x[0]);
	}

	/// @brief Wrapper of two_two_sub_full.
	inline void Two_Two_Diff(const double *a, const double *b, double *x)
	{
		two_two_sub_full(a[1], a[0], b[1], b[0], x[3], x[2], x[1], x[0]);
	}

	/// @brief Calculates the second component 'y' of the expansion [x,y] =
	/// [a]-[b] when 'x' is known
	inline void Two_Diff_Back(const double a, const double b, double &x,
	                          double &y)
	{
		double _bv;
		_bv = a - x;
		y   = (a - (x + _bv)) + (_bv - b);
	}

	/// @brief Wrapper of Two_Diff_Back
	inline void Two_Diff_Back(const double a, const double b, double *x)
	{
		Two_Diff_Back(a, b, x[1], x[0]);
	}

	/* Product routines ********************************************************/

	/// @brief Wrapper of two_prod.
	inline void Two_Prod(const double a, const double b, double *x)
	{
		two_prod(a, b, x[1], x[0]);
	}

	/// @brief Wrapper of two_square.
	inline void Square(const double a, double *x) { one_square(a, x[1], x[0]); }

	/// @brief Wrapper of two_one_prod_full.
	inline void Two_One_Prod(const double *a, const double b, double *x)
	{
		two_one_prod_full(a[1], a[0], b, x[3], x[2], x[1], x[0]);
	}

	// [x] = [a1,a0]^2		Squares an expansion of length 2
	// 'x' must be allocated by the caller with 6 components.
	inline void Two_Square(const double &a1, const double &a0, double *x)
	{
		two_square_full(a1, a0, x);
	}

	// [h7,h6,...,h0] = [a1,a0]*[b1,b0]		Calculates the product of two expansions
	// of length two. 'h' must be allocated by the caller with eight components.
	inline void Two_Two_Prod(const double *a, const double *b, double *x)
	{
		two_two_prod_full(a[1], a[0], b[1], b[0], x);
	}

	/****************************************************************************
	 * Unlimited length floating-point expansion arithmetic,
	 * due to Bruno Levy, Macro Attene, and others.
	 ***************************************************************************/

	// [e] = -[e]		Inplace inversion
	inline void Gen_Invert(const int elen, double *e)
	{
		for (int i = 0; i < elen; i++)
			e[i] = -e[i];
	}

	// [h] = [e] + [f]		Sums two expansions and returns number of components of
	// result 'h' must be allocated by the caller with at least elen+flen
	// components.
	inline int Gen_Sum(const int elen, const double *e, const int flen,
	                   const double *f, double *h);

	// Same as above, but 'h' is allocated internally. The caller must still call
	// 'free' to release the memory.
	inline int Gen_Sum_With_Alloc(const int elen, const double *e, const int flen,
	                              const double *f, double **h)
	{
		*h = AllocDoubles(elen + flen);
		return Gen_Sum(elen, e, flen, f, *h);
	}

	// [h] = [e] + [f]		Subtracts two expansions and returns number of
	// components of result 'h' must be allocated by the caller with at least
	// elen+flen components.
	inline int Gen_Diff(const int elen, const double *e, const int flen,
	                    const double *f, double *h);

	// Same as above, but 'h' is allocated internally. The caller must still call
	// 'free' to release the memory.
	inline int Gen_Diff_With_Alloc(const int elen, const double *e,
	                               const int flen, const double *f, double **h)
	{
		*h = AllocDoubles(elen + flen);
		return Gen_Diff(elen, e, flen, f, *h);
	}

	// [h] = [e] * b		Multiplies an expansion by a scalar
	// 'h' must be allocated by the caller with at least elen*2 components.
	inline int Gen_Scale(const int elen, const double *e, const double &b,
	                     double *h);

	// [h] = [e] * 2		Multiplies an expansion by 2
	// 'h' must be allocated by the caller with at least elen components. This is
	// exact up to overflows.
	inline void Double(const int elen, const double *e, double *h)
	{
		for (int i = 0; i < elen; i++)
			h[i] = 2 * e[i];
	}

	// [h] = [e] * 2		Multiplies an expansion by n
	// If 'n' is a power of two, the multiplication is exact
	inline void ExactScale(const int elen, double *e, const double n)
	{
		for (int i = 0; i < elen; i++)
			e[i] *= n;
	}

	// [h] = [a] * [b]
	// 'h' must be allocated by the caller with at least 2*alen*blen components.
	inline int Sub_product(const int alen, const double *a, const int blen,
	                       const double *b, double *h);

	// [h] = [a] * [b]
	// 'h' must be allocated by the caller with at least MAX(2*alen*blen, 8)
	// components.
	inline int Gen_Product(const int alen, const double *a, const int blen,
	                       const double *b, double *h);

	// Same as above, but 'h' is allocated internally. The caller must still call
	// 'free' to release the memory.
	inline int Gen_Product_With_Alloc(const int alen, const double *a,
	                                  const int blen, const double *b,
	                                  double **h);

	// Assume that *h is pre-allocated with hlen doubles.
	// If more elements are required, *h is re-allocated internally.
	// In any case, the function returns the size of the resulting expansion.
	// The caller must verify whether reallocation took place, and possibly call
	// 'free' to release the memory. When reallocation takes place, *h becomes
	// different from its original value.

	inline int Double_With_PreAlloc(const int elen, const double *e, double **h,
	                                const int hlen);

	inline int Gen_Scale_With_PreAlloc(const int elen, const double *e,
	                                   const double &b, double **h,
	                                   const int hlen);

	inline int Gen_Sum_With_PreAlloc(const int elen, const double *e,
	                                 const int flen, const double *f, double **h,
	                                 const int hlen);

	inline int Gen_Diff_With_PreAlloc(const int elen, const double *e,
	                                  const int flen, const double *f, double **h,
	                                  const int hlen);

	inline int Gen_Product_With_PreAlloc(const int alen, const double *a,
	                                     const int blen, const double *b,
	                                     double **h, const int hlen);

	// Approximates the expansion to a double
	inline double To_Double(const int elen, const double *e);
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ExpansionObject.inl"
#endif