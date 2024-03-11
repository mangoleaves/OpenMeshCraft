#pragma once

#include "ExpansionObject.h"

namespace OMC {

void expansionObject::two_Sum(const double a, const double b, double *xy)
{
	Two_Sum(a, b, xy[1], xy[0]);
}

void expansionObject::two_Diff(const double a, const double b, double *xy)
{
	Two_Diff(a, b, xy[1], xy[0]);
}

void expansionObject::Two_Prod(const double a, const double b, double &x,
                               double &y)
{
	x = a * b;
	//			y = fma(a, b, -x);
	Split(a, _ah, _al);
	Split(b, _bh, _bl);
	y = ((_ah * _bh - x) + _ah * _bl + _al * _bh) + _al * _bl;
}

void expansionObject::Square(const double a, double &x, double &y)
{
	x = a * a;
	Split(a, _ah, _al);
	y = (_al * _al) - ((x - (_ah * _ah)) - ((_ah + _ah) * _al));
}

void expansionObject::two_One_Diff(const double a1, const double a0,
                                   const double b, double &x2, double &x1,
                                   double &x0)
{
	Two_One_Diff(a1, a0, b, x2, x1, x0);
}

void expansionObject::Two_One_Prod(const double a1, const double a0,
                                   const double b, double &x3, double &x2,
                                   double &x1, double &x0)
{
	Split(b, _bh, _bl);
	Two_Prod_PreSplit(a0, b, _bh, _bl, _i, x0);
	Two_Prod_PreSplit(a1, b, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _k, x1);
	Quick_Two_Sum(_j, _k, x3, x2);
}

void expansionObject::Two_Two_Sum(const double a1, const double a0,
                                  const double b1, const double b0, double &x3,
                                  double &x2, double &x1, double &x0)
{
	Two_One_Sum(a1, a0, b0, _j, _0, x0);
	Two_One_Sum(_j, _0, b1, x3, x2, x1);
}

void expansionObject::Two_Two_Diff(const double a1, const double a0,
                                   const double b1, const double b0, double &x3,
                                   double &x2, double &x1, double &x0)
{
	Two_One_Diff(a1, a0, b0, _j, _0, x0);
	Two_One_Diff(_j, _0, b1, _u3, x2, x1);
	x3 = _u3;
}

void expansionObject::Two_Two_Prod(const double a1, const double a0,
                                   const double b1, const double b0, double *h)
{
	double _ch, _cl, _m, _n;
	Split(a0, _ah, _al);
	Split(b0, _bh, _bl);
	Two_Product_2Presplit(a0, _ah, _al, b0, _bh, _bl, _i, h[0]);
	Split(a1, _ch, _cl);
	Two_Product_2Presplit(a1, _ch, _cl, b0, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _k, _1);
	Quick_Two_Sum(_j, _k, _l, _2);
	Split(b1, _bh, _bl);
	Two_Product_2Presplit(a0, _ah, _al, b1, _bh, _bl, _i, _0);
	Two_Sum(_1, _0, _k, h[1]);
	Two_Sum(_2, _k, _j, _1);
	Two_Sum(_l, _j, _m, _2);
	Two_Product_2Presplit(a1, _ch, _cl, b1, _bh, _bl, _j, _0);
	Two_Sum(_i, _0, _n, _0);
	Two_Sum(_1, _0, _i, h[2]);
	Two_Sum(_2, _i, _k, _1);
	Two_Sum(_m, _k, _l, _2);
	Two_Sum(_j, _n, _k, _0);
	Two_Sum(_1, _0, _j, h[3]);
	Two_Sum(_2, _j, _i, _1);
	Two_Sum(_l, _i, _m, _2);
	Two_Sum(_1, _k, _i, h[4]);
	Two_Sum(_2, _i, _k, h[5]);
	Two_Sum(_m, _k, h[7], h[6]);
}

int expansionObject::Gen_Sum(const int elen, const double *e, const int flen,
                             const double *f, double *h)
{
	double Q, Qn, hh, en = e[0], fn = f[0];
	int    e_k, f_k, h_k;

	h_k = e_k = f_k = 0;
	if ((fn > en) == (fn > -en))
	{
		Q = en;
		e_k++;
	}
	else
	{
		Q = fn;
		f_k++;
	}

	if ((e_k < elen) && (f_k < flen))
	{
		en = e[e_k];
		fn = f[f_k];
		if ((fn > en) == (fn > -en))
		{
			Quick_Two_Sum(en, Q, Qn, hh);
			e_k++;
		}
		else
		{
			Quick_Two_Sum(fn, Q, Qn, hh);
			f_k++;
		}
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
		while ((e_k < elen) && (f_k < flen))
		{
			en = e[e_k];
			fn = f[f_k];
			if ((fn > en) == (fn > -en))
			{
				Two_Sum(Q, en, Qn, hh);
				e_k++;
			}
			else
			{
				Two_Sum(Q, fn, Qn, hh);
				f_k++;
			}
			Q = Qn;
			if (hh != 0.0)
				h[h_k++] = hh;
		}
	}

	while (e_k < elen)
	{
		en = e[e_k++];
		Two_Sum(Q, en, Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
	}

	while (f_k < flen)
	{
		fn = f[f_k++];
		Two_Sum(Q, fn, Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
	}
	if ((Q != 0.0) || (h_k == 0))
		h[h_k++] = Q;

	return h_k;
}

int expansionObject::Gen_Diff(const int elen, const double *e, const int flen,
                              const double *f, double *h)
{
	double Q, Qn, hh, en = e[0], fn = -f[0];
	int    e_k, f_k, h_k;

	h_k = e_k = f_k = 0;
	if ((fn > en) == (fn > -en))
	{
		Q = en;
		e_k++;
	}
	else
	{
		Q = fn;
		f_k++;
	}

	if ((e_k < elen) && (f_k < flen))
	{
		en = e[e_k];
		fn = -f[f_k];
		if ((fn > en) == (fn > -en))
		{
			Quick_Two_Sum(en, Q, Qn, hh);
			e_k++;
		}
		else
		{
			Quick_Two_Sum(fn, Q, Qn, hh);
			f_k++;
		}
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
		while ((e_k < elen) && (f_k < flen))
		{
			en = e[e_k];
			fn = -f[f_k];
			if ((fn > en) == (fn > -en))
			{
				Two_Sum(Q, en, Qn, hh);
				e_k++;
			}
			else
			{
				Two_Sum(Q, fn, Qn, hh);
				f_k++;
			}
			Q = Qn;
			if (hh != 0.0)
				h[h_k++] = hh;
		}
	}

	while (e_k < elen)
	{
		en = e[e_k++];
		Two_Sum(Q, en, Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
	}

	while (f_k < flen)
	{
		fn = -f[f_k++];
		Two_Sum(Q, fn, Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
	}
	if ((Q != 0.0) || (h_k == 0))
		h[h_k++] = Q;

	return h_k;
}

int expansionObject::Gen_Scale(const int elen, const double *e, const double &b,
                               double *h)
{
	double Q, sum, hh, pr1, pr0, enow;
	int    e_k, h_k;

	Split(b, _bh, _bl);
	Two_Prod_PreSplit(e[0], b, _bh, _bl, Q, hh);
	h_k = 0;
	if (hh != 0)
		h[h_k++] = hh;

	for (e_k = 1; e_k < elen; e_k++)
	{
		enow = e[e_k];
		Two_Prod_PreSplit(enow, b, _bh, _bl, pr1, pr0);
		Two_Sum(Q, pr0, sum, hh);
		if (hh != 0)
			h[h_k++] = hh;
		Quick_Two_Sum(pr1, sum, Q, hh);
		if (hh != 0)
			h[h_k++] = hh;
	}
	if ((Q != 0.0) || (h_k == 0))
		h[h_k++] = Q;

	return h_k;
}

void expansionObject::Two_Square(const double &a1, const double &a0, double *x)
{
	Square(a0, _j, x[0]);
	_0 = a0 + a0;
	Two_Prod(a1, _0, _k, _1);
	Two_One_Sum(_k, _1, _j, _l, _2, x[1]);
	Square(a1, _j, _1);
	Two_Two_Sum(_j, _1, _l, _2, x[5], x[4], x[3], x[2]);
}

int expansionObject::Sub_product(const int alen, const double *a,
                                 const int blen, const double *b, double *h)
{
	if (alen == 1)
		return Gen_Scale(blen, b, a[0], h);
	int     partial = 2 * alen * blen;
	int     allmem  = 2 * (partial + blen);
	double  ph1_p[1024];
	double *ph1   = (allmem > 1024) ? (AllocDoubles(allmem)) : (ph1_p);
	double *ph2   = ph1 + partial;
	double *th    = ph2 + partial;
	double *ph[2] = {ph1, ph2};
	int     first = 0;
	int     phl   = Gen_Scale(blen, b, a[0], ph[0]);

	for (int i = 1; i < alen; i++)
	{
		int thl = Gen_Scale(blen, b, a[i], th);
		first   = i & 1;
		phl     = Gen_Sum(phl, ph[(i + 1) & 1], thl, th, ph[first]);
	}
	if (first)
		for (int i = 0; i < phl; i++)
			h[i] = ph2[i];
	else
		for (int i = 0; i < phl; i++)
			h[i] = ph1[i];
	if (allmem > 1024)
		FreeDoubles(ph1);
	return phl;
}

int expansionObject::Gen_Product(const int alen, const double *a,
                                 const int blen, const double *b, double *h)
{
	if (blen == 1)
		return Gen_Scale(alen, a, b[0], h);
	else if (alen < blen)
		return Sub_product(alen, a, blen, b, h);
	else
		return Sub_product(blen, b, alen, a, h);
}

double expansionObject::To_Double(const int elen, const double *e)
{
	double Q = e[0];
	for (int e_i = 1; e_i < elen; e_i++)
		Q += e[e_i];
	return Q;
}

int expansionObject::Gen_Product_With_Alloc(const int alen, const double *a,
                                            const int blen, const double *b,
                                            double **h)
{
	int h_len = alen * blen * 2;
	if (h_len < 8)
		h_len = 8;
	*h = AllocDoubles(h_len);
	return Gen_Product(alen, a, blen, b, *h);
}

int expansionObject::Double_With_PreAlloc(const int elen, const double *e,
                                          double **h, const int hlen)
{
	int newlen = elen;
	if (hlen < newlen)
		*h = AllocDoubles(newlen);
	// if (hlen < newlen) printf("REALLOC %d bytes\n", newlen);
	Double(elen, e, *h);
	return newlen;
}

int expansionObject::Gen_Scale_With_PreAlloc(const int elen, const double *e,
                                             const double &b, double **h,
                                             const int hlen)
{
	int newlen = elen * 2;
	if (hlen < newlen)
		*h = AllocDoubles(newlen);
	return Gen_Scale(elen, e, b, *h);
}

int expansionObject::Gen_Sum_With_PreAlloc(const int elen, const double *e,
                                           const int flen, const double *f,
                                           double **h, const int hlen)
{
	int newlen = elen + flen;
	if (hlen < newlen)
		*h = AllocDoubles(newlen);
	return Gen_Sum(elen, e, flen, f, *h);
}

int expansionObject::Gen_Diff_With_PreAlloc(const int elen, const double *e,
                                            const int flen, const double *f,
                                            double **h, const int hlen)
{
	int newlen = elen + flen;
	if (hlen < newlen)
		*h = AllocDoubles(newlen);
	return Gen_Diff(elen, e, flen, f, *h);
}

int expansionObject::Gen_Product_With_PreAlloc(const int alen, const double *a,
                                               const int blen, const double *b,
                                               double **h, const int hlen)
{
	int newlen = alen * blen * 2;
	if (hlen < newlen || hlen < 8)
	{
		if (newlen < 8)
			newlen = 8;
		*h = AllocDoubles(newlen);
	}
	return Gen_Product(alen, a, blen, b, *h);
}

} // namespace OMC

// Sums
#undef Quick_Two_Sum
#undef Two_Sum
#undef Two_One_Sum

// Differences
#undef Two_Diff
#undef Two_One_Diff

// Products
#undef Split
#undef Two_Prod_PreSplit
#undef Two_Product_2Presplit