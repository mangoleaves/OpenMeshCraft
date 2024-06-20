#pragma once

#include "ExpansionObject.h"

#include "FPU.h"

namespace OMC {

int expansionObject::Gen_Sum(const int elen, const double *e, const int flen,
                             const double *f, double *h)
{
	double        Q, Qn, hh, s;
	const double *en = e, *fn = f, *elast = e + elen, *flast = f + flen;
	int           h_k = 0;

	Q = (fabs(*fn) > fabs(*en)) ? (*en++) : (*fn++);

	if ((en < elast) && (fn < flast))
	{
		s = (fabs(*fn) > fabs(*en)) ? (*en++) : (*fn++);
		fast_two_sum(s, Q, Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
		while ((en < elast) && (fn < flast))
		{
			s = (fabs(*fn) > fabs(*en)) ? (*en++) : (*fn++);
			two_sum(s, Q, Qn, hh);
			Q = Qn;
			if (hh != 0.0)
				h[h_k++] = hh;
		}
	}

	while (en < elast)
	{
		two_sum(Q, (*en), Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
		en++;
	}

	while (fn < flast)
	{
		two_sum(Q, (*fn), Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
		fn++;
	}
	if ((Q != 0.0) || (h_k == 0))
		h[h_k++] = Q;

	return h_k;
}

int expansionObject::Gen_Diff(const int elen, const double *e, const int flen,
                              const double *f, double *h)
{
	double        Q, Qn, hh, s;
	const double *en = e, *fn = f, *elast = e + elen, *flast = f + flen;
	int           h_k = 0;

	Q = (fabs(*fn) > fabs(*en)) ? (*en++) : (-*fn++);

	if ((en < elast) && (fn < flast))
	{
		s = (fabs(*fn) > fabs(*en)) ? (*en++) : (-*fn++);
		fast_two_sum(s, Q, Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
		while ((en < elast) && (fn < flast))
		{
			s = (fabs(*fn) > fabs(*en)) ? (*en++) : (-*fn++);
			two_sum(s, Q, Qn, hh);
			Q = Qn;
			if (hh != 0.0)
				h[h_k++] = hh;
		}
	}

	while (en < elast)
	{
		two_sum(Q, (*en), Qn, hh);
		Q = Qn;
		if (hh != 0.0)
			h[h_k++] = hh;
		en++;
	}

	while (fn < flast)
	{
		s = *fn++;
		two_sub(Q, s, Qn, hh);
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
	double        Q, sum, hh, pr1, pr0;
	const double *ei = e, *elast = e + elen;

	int k = 0;
	two_prod(*ei, b, Q, hh);
	if (hh != 0)
		h[k++] = hh;

	while (++ei < elast)
	{
		two_prod(*ei, b, pr1, pr0);
		two_sum(Q, pr0, sum, hh);
		if (hh != 0)
			h[k++] = hh;
		fast_two_sum(pr1, sum, Q, hh);
		if (hh != 0)
			h[k++] = hh;
	}
	if ((Q != 0.0) || (k == 0))
		h[k++] = Q;
	return k;
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

int expansionObject::Double_With_PreAlloc(const int elen, const double *e,
                                          double **h, const int hlen)
{
	int newlen = elen;
	if (hlen < newlen)
		*h = AllocDoubles(newlen);
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

void expansionObject::Compress(int &elen, double *e)
{
	double *h = e;

	int    m = elen;
	double Qnew, q;

	int    bottom = m - 1;
	double Q      = e[bottom];

	for (int i = m - 2; i >= 0; --i)
	{
		fast_two_sum(Q, e[i], Qnew, q);
		Q = Qnew;
		if (q != 0.0)
		{
			h[bottom] = Q;
			--bottom;
			Q = q;
		}
	}
	h[bottom] = Q;

	int top = 0;
	for (int i = bottom + 1; i < m; ++i)
	{
		fast_two_sum(h[i], Q, Qnew, q);
		Q = Qnew;
		if (q != 0.0)
		{
			h[top] = q;
			++top;
		}
	}
	h[top] = Q;
	OMC_PRED_PROFILE_SAVE_COMPRS(elen, top + 1);
	elen = top + 1;
}

void expansionObject::CompressIf(int &elen, double *e)
{
	if (elen > compress_thres)
		Compress(elen, e);
}

double expansionObject::To_Double(const int elen, const double *e)
{
	double Q = e[0];
	for (int e_i = 1; e_i < elen; e_i++)
		Q += e[e_i];
	return Q;
}

std::pair<double, double> expansionObject::To_Interval(const int     elen,
                                                       const double *e)
{
	// Optimized expansion-to-interval conversion:
	//
	// Add components starting from the one of largest magnitude
	// Stop as soon as next component is smaller than ulp (and then
	// expand interval by ulp).

	// rounding to +infinity
	FPU_RoundingProtector<std::false_type> P;

	int    l    = elen;
	double ninf = -e[l - 1];
	double sup  = e[l - 1];

	for (int comp_idx = l - 2; comp_idx >= 0; --comp_idx)
	{
		double comp = e[comp_idx];
		sup += comp;
		ninf -= comp;
		if (comp > 0)
		{
			double new_sup = sup + comp;
			if (new_sup == sup)
			{
				sup = std::nextafter(sup, std::numeric_limits<double>::infinity());
				break;
			}
			else
			{
				sup = new_sup;
			}
		}
		else
		{
			// If we stored inf, we would write:
			//  new_inf  =  inf  + comp
			// But we store ninf = -inf, so we write:
			// -new_ninf = -ninf + comp
			// Which means:
			//  new_ninf =  ninf - comp

			double new_ninf = ninf - comp;
			if (new_ninf == ninf)
			{
				ninf = std::nextafter(ninf, std::numeric_limits<double>::infinity());
				break;
			}
			else
			{
				ninf = new_ninf;
			}
		}
	}

	return std::pair<double, double>(ninf, sup);
}

} // namespace OMC