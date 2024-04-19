#pragma once

#include "IndirectPredicateDetailsHand.h"

#include "OpenMeshCraft/NumberTypes/ExpansionObject.h"
#include "OpenMeshCraft/NumberTypes/IntervalNumber.h"
#include "OpenMeshCraft/NumberTypes/LazyNumber.h"

// Uncomment the following to activate overflow/underflow checks
#define CHECK_FOR_XYZERFLOWS

// #define OMC_NO_SHEWCHUK

#ifndef OMC_NO_SHEWCHUK
extern "C"
{
	double orient2d(const double *pa, const double *pb, const double *pc);

	// In this orient3d, generally, pb, pc and pd come from three vertices of a
	// triangle.
	// NOTE In all orient3d in our lib, the first three points come from
	// three vertices of a triangle.
	double orient3d(const double *pa, const double *pb, const double *pc,
	                const double *pd);

	void orient3d_get_minors(const double *pb, const double *pc, const double *pd,
	                         double *minor, double *perm);

	double orient3d_with_cached_minors(const double *pa, const double *pb,
	                                   const double *pc, const double *pd,
	                                   double *minor, double *perm);

	double incircle(const double *pa, const double *pb, const double *pc,
	                const double *pd);

	double insphere(const double *pa, const double *pb, const double *pc,
	                const double *pd, const double *pe);
}
#endif

namespace OMC {

inline Sign orient2d_filtered(double p1x, double p1y, double p2x, double p2y,
                              double p3x, double p3y)
{
	double dl  = (p2x - p1x) * (p3y - p1y);
	double dr  = (p2y - p1y) * (p3x - p1x);
	double det = dl - dr;
	double eb  = 3.3306690738754706e-016 * (fabs(dl) + fabs(dr));
	return filter_sign(det, eb);
}

inline Sign orient2d_expansion(double p1x, double p1y, double p2x, double p2y,
                               double p3x, double p3y)
{
	expansionObject o;
	double          acx[2], acy[2], bcx[2], bcy[2], dtl[2], dtr[2], B[4];
	double          s[2], t[2], u[4], C1[8], C2[12], D[16];
	int             C1l, C2l, Dl;

	acx[1] = (p1x - p3x);
	bcx[1] = (p2x - p3x);
	acy[1] = (p1y - p3y);
	bcy[1] = (p2y - p3y);

	o.Two_Prod(acx[1], bcy[1], dtl);
	o.Two_Prod(acy[1], bcx[1], dtr);
	o.Two_Two_Diff(dtl, dtr, B);

	double dsm = (fabs(dtl[1]) + fabs(dtr[1]));
	double det = o.To_Double(4, B);
	double eb  = 2.2204460492503146e-16 * dsm;
	Dl         = ((det >= eb) - (-det >= eb));
	if (Dl)
		return static_cast<Sign>(Dl);

	o.Two_Diff_Back(p1x, p3x, acx);
	o.Two_Diff_Back(p2x, p3x, bcx);
	o.Two_Diff_Back(p1y, p3y, acy);
	o.Two_Diff_Back(p2y, p3y, bcy);

	if ((acx[0] == 0.0) && (acy[0] == 0.0) && (bcx[0] == 0.0) && (bcy[0] == 0.0))
		return static_cast<Sign>((det > 0) - (det < 0));

	eb = 1.1093356479670487e-31 * dsm + 3.3306690738754706e-16 * fabs(det);
	det +=
	  (acx[1] * bcy[0] + bcy[1] * acx[0]) - (acy[1] * bcx[0] + bcx[1] * acy[0]);
	Dl = ((det >= eb) - (-det >= eb));
	if (Dl)
		return static_cast<Sign>(Dl);

	o.Two_Prod(acx[0], bcy[1], s);
	o.Two_Prod(acy[0], bcx[1], t);
	o.Two_Two_Diff(s, t, u);
	C1l = o.Gen_Sum(4, B, 4, u, C1);

	o.Two_Prod(acx[1], bcy[0], s);
	o.Two_Prod(acy[1], bcx[0], t);
	o.Two_Two_Diff(s, t, u);
	C2l = o.Gen_Sum(C1l, C1, 4, u, C2);

	o.Two_Prod(acx[0], bcy[0], s);
	o.Two_Prod(acy[0], bcx[0], t);
	o.Two_Two_Diff(s, t, u);
	Dl = o.Gen_Sum(C2l, C2, 4, u, D);

	det = D[Dl - 1];
	return static_cast<Sign>((det >= eb) - (-det >= eb));
}

inline Sign orient2d(double p1x, double p1y, double p2x, double p2y, double p3x,
                     double p3y)
{
	Sign ret = orient2d_filtered(p1x, p1y, p2x, p2y, p3x, p3y);
	if (is_sign_reliable(ret))
		return ret;
	return orient2d_expansion(p1x, p1y, p2x, p2y, p3x, p3y);
}

inline Sign orient2d(const double *p1, const double *p2, const double *p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2d(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
#else
	return OMC::sign(::orient2d(p1, p2, p3));
#endif
}

template <typename IT, typename ET>
Sign orient2d(const GenericPoint2T<IT, ET> &p1,
              const GenericPoint2T<IT, ET> &p2,
              const GenericPoint2T<IT, ET> &p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2d(p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y());
#else
	return OMC::sign(::orient2d(p1.data(), p2.data(), p3.data()));
#endif
}

inline Sign orient3d_filtered(double px, double py, double pz, double qx,
                              double qy, double qz, double rx, double ry,
                              double rz, double sx, double sy, double sz)
{
	double fadx, fbdx, fcdx, fady, fbdy, fcdy, fadz, fbdz, fcdz, eb;
	double fbdxcdy, fcdxbdy, fcdxady, fadxcdy, fadxbdy, fbdxady, det;

	fadx = qx - px;
	fbdx = rx - px;
	fcdx = sx - px;
	fady = qy - py;
	fbdy = ry - py;
	fcdy = sy - py;
	fadz = qz - pz;
	fbdz = rz - pz;
	fcdz = sz - pz;

	fbdxcdy = fbdx * fcdy * fadz;
	fcdxbdy = fcdx * fbdy * fadz;
	fcdxady = fcdx * fady * fbdz;
	fadxcdy = fadx * fcdy * fbdz;
	fadxbdy = fadx * fbdy * fcdz;
	fbdxady = fbdx * fady * fcdz;

	det = (fbdxcdy - fcdxbdy) + (fcdxady - fadxcdy) + (fadxbdy - fbdxady);
	eb =
	  7.7715611723761027e-016 * (fabs(fbdxcdy) + fabs(fcdxbdy) + fabs(fcdxady) +
	                             fabs(fadxcdy) + fabs(fadxbdy) + fabs(fbdxady));
	return filter_sign(det, eb);
}

inline void supo3d1(double *c1, double *c2, double *c3, double *c4, double *c5,
                    double *c6, double *a1, double *a2, double &i, double *k1,
                    double *k2, double *k3, double *k4, int &l1, int &l2)
{
	expansionObject o;
	if (c1[0] == 0.0)
	{
		if (c2[0] == 0.0)
		{
			a1[0] = a2[0] = 0.0;
			l1 = l2 = 1;
		}
		else
		{
			i = -c2[0];
			o.Two_Prod(i, c3[1], a1);
			o.Two_Prod(c2[0], c4[1], a2);
			l1 = l2 = 2;
		}
	}
	else
	{
		if (c2[0] == 0.0)
		{
			i = -c1[0];
			o.Two_Prod(c1[0], c5[1], a1);
			o.Two_Prod(i, c6[1], a2);
			l1 = l2 = 2;
		}
		else
		{
			o.Two_Prod(c1[0], c5[1], k1);
			o.Two_Prod(c2[0], c3[1], k2);
			o.Two_Two_Diff(k1, k2, a1);
			o.Two_Prod(c2[0], c4[1], k3);
			o.Two_Prod(c1[0], c6[1], k4);
			o.Two_Two_Diff(k3, k4, a2);
			l1 = l2 = 4;
		}
	}
}

inline void supo3d2(double *c1, double *c2, double *c3, double *c4, double *u,
                    int &fl, double fin[2][192], int &wh, double *c5, double &i,
                    double *c6, double *c7)

{
	expansionObject o;
	if (c1[0] != 0.0)
	{
		if (c2[0] != 0.0)
		{
			o.Two_Prod(c1[0], c2[0], c3);
			o.Two_One_Prod(c3, c4[1], u);
			fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
			wh = !wh;
			if (c4[0] != 0.0)
			{
				o.Two_One_Prod(c3, c4[0], u);
				fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
				wh = !wh;
			}
		}
		if (c5[0] != 0.0)
		{
			i = -c1[0];
			o.Two_Prod(i, c5[0], c6);
			o.Two_One_Prod(c6, c7[1], u);
			fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
			wh = !wh;
			if (c7[0] != 0.0)
			{
				o.Two_One_Prod(c6, c7[0], u);
				fl = o.Gen_Sum(fl, fin[wh], 4, u, fin[!wh]);
				wh = !wh;
			}
		}
	}
}

inline Sign orient3d_expansion(double pdx, double pdy, double pdz, double pax,
                               double pay, double paz, double pbx, double pby,
                               double pbz, double pcx, double pcy, double pcz)
{
	expansionObject o;
	double          eb, det;
	double adx[2], bdx[2], cdx[2], ady[2], bdy[2], cdy[2], adz[2], bdz[2], cdz[2];
	double bdxcdy[2], cdxbdy[2], cdxady[2], adxcdy[2], adxbdy[2], bdxady[2];
	double bc[4], ca[4], ab[4];
	double bdxt_cdy[2], cdxt_bdy[2], cdxt_ady[2];
	double adxt_cdy[2], adxt_bdy[2], bdxt_ady[2];
	double bdyt_cdx[2], cdyt_bdx[2], cdyt_adx[2];
	double adyt_cdx[2], adyt_bdx[2], bdyt_adx[2];
	double bdxt_cdyt[2], cdxt_bdyt[2], cdxt_adyt[2];
	double adxt_cdyt[2], adxt_bdyt[2], bdxt_adyt[2];
	double u[4], v[12], w[16];
	double adet[8], bdet[8], cdet[8], abdet[16];
	double fin[2][192];
	int    wh = 0;
	double at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
	double bct[8], cat[8], abt[8];
	int    alen, blen, clen, finlen, vlen, wlen;
	int    at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
	int    bctlen, catlen, abtlen;
	int    ablen;
	double inv;
	int    ri;

	adx[1] = pax - pdx;
	bdx[1] = pbx - pdx;
	cdx[1] = pcx - pdx;
	ady[1] = pay - pdy;
	bdy[1] = pby - pdy;
	cdy[1] = pcy - pdy;
	adz[1] = paz - pdz;
	bdz[1] = pbz - pdz;
	cdz[1] = pcz - pdz;

	o.Two_Prod(bdx[1], cdy[1], bdxcdy);
	o.Two_Prod(cdx[1], bdy[1], cdxbdy);
	o.Two_Two_Diff(bdxcdy, cdxbdy, bc);
	alen = o.Gen_Scale(4, bc, adz[1], adet);

	o.Two_Prod(cdx[1], ady[1], cdxady);
	o.Two_Prod(adx[1], cdy[1], adxcdy);
	o.Two_Two_Diff(cdxady, adxcdy, ca);
	blen = o.Gen_Scale(4, ca, bdz[1], bdet);

	o.Two_Prod(adx[1], bdy[1], adxbdy);
	o.Two_Prod(bdx[1], ady[1], bdxady);
	o.Two_Two_Diff(adxbdy, bdxady, ab);
	clen = o.Gen_Scale(4, ab, cdz[1], cdet);

	ablen  = o.Gen_Sum(alen, adet, blen, bdet, abdet);
	finlen = o.Gen_Sum(ablen, abdet, clen, cdet, fin[wh]);

	double xx1 = bdxcdy[1] * adz[1];
	double xx2 = cdxbdy[1] * adz[1];
	double yy1 = cdxady[1] * bdz[1];
	double yy2 = adxcdy[1] * bdz[1];
	double zz1 = adxbdy[1] * cdz[1];
	double zz2 = bdxady[1] * cdz[1];
	double pm =
	  fabs(xx1) + fabs(xx2) + fabs(yy1) + fabs(yy2) + fabs(zz1) + fabs(zz2);

	det = o.To_Double(finlen, fin[wh]);
	eb  = 3.3306690738754731e-016 * pm;
	ri  = (det >= eb) - (-det >= eb);
	if (ri)
		return static_cast<Sign>(ri);

	o.Two_Diff_Back(pax, pdx, adx);
	o.Two_Diff_Back(pbx, pdx, bdx);
	o.Two_Diff_Back(pcx, pdx, cdx);
	o.Two_Diff_Back(pay, pdy, ady);
	o.Two_Diff_Back(pby, pdy, bdy);
	o.Two_Diff_Back(pcy, pdy, cdy);
	o.Two_Diff_Back(paz, pdz, adz);
	o.Two_Diff_Back(pbz, pdz, bdz);
	o.Two_Diff_Back(pcz, pdz, cdz);

	if ((adx[0] == 0.0) && (bdx[0] == 0.0) && (cdx[0] == 0.0) &&
	    (ady[0] == 0.0) && (bdy[0] == 0.0) && (cdy[0] == 0.0) &&
	    (adz[0] == 0.0) && (bdz[0] == 0.0) && (cdz[0] == 0.0))
		return static_cast<Sign>((det > 0) - (det < 0));

	eb = 3.2047474274603644e-031 * pm + 1.1102230246251565e-016 * fabs(det);
	det += (adz[1] * ((bdx[1] * cdy[0] + cdy[1] * bdx[0]) -
	                  (bdy[1] * cdx[0] + cdx[1] * bdy[0])) +
	        adz[0] * (bdx[1] * cdy[1] - bdy[1] * cdx[1])) +
	       (bdz[1] * ((cdx[1] * ady[0] + ady[1] * cdx[0]) -
	                  (cdy[1] * adx[0] + adx[1] * cdy[0])) +
	        bdz[0] * (cdx[1] * ady[1] - cdy[1] * adx[1])) +
	       (cdz[1] * ((adx[1] * bdy[0] + bdy[1] * adx[0]) -
	                  (ady[1] * bdx[0] + bdx[1] * ady[0])) +
	        cdz[0] * (adx[1] * bdy[1] - ady[1] * bdx[1]));
	ri = (det >= eb) - (-det >= eb);
	if (ri)
		return static_cast<Sign>(ri);

	// Filters did not work. Compute exactly...
	supo3d1(adx, ady, bdx, cdx, bdy, cdy, at_b, at_c, inv, adxt_bdy, adyt_bdx,
	        adyt_cdx, adxt_cdy, at_blen, at_clen);

	supo3d1(bdx, bdy, cdx, adx, cdy, ady, bt_c, bt_a, inv, bdxt_cdy, bdyt_cdx,
	        bdyt_adx, bdxt_ady, bt_alen, bt_clen);

	supo3d1(cdx, cdy, adx, bdx, ady, bdy, ct_a, ct_b, inv, cdxt_ady, cdyt_adx,
	        cdyt_bdx, cdxt_bdy, ct_alen, ct_blen);

	bctlen = o.Gen_Sum(bt_clen, bt_c, ct_blen, ct_b, bct);
	wlen   = o.Gen_Scale(bctlen, bct, adz[1], w);
	finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[(!wh)]);
	wh     = !wh;

	catlen = o.Gen_Sum(ct_alen, ct_a, at_clen, at_c, cat);
	wlen   = o.Gen_Scale(catlen, cat, bdz[1], w);
	finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[(!wh)]);
	wh     = !wh;

	abtlen = o.Gen_Sum(at_blen, at_b, bt_alen, bt_a, abt);
	wlen   = o.Gen_Scale(abtlen, abt, cdz[1], w);
	finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[(!wh)]);
	wh     = !wh;

	if (adz[0] != 0.0)
	{
		vlen   = o.Gen_Scale(4, bc, adz[0], v);
		finlen = o.Gen_Sum(finlen, fin[wh], vlen, v, fin[(!wh)]);
		wh     = !wh;
	}
	if (bdz[0] != 0.0)
	{
		vlen   = o.Gen_Scale(4, ca, bdz[0], v);
		finlen = o.Gen_Sum(finlen, fin[wh], vlen, v, fin[(!wh)]);
		wh     = !wh;
	}
	if (cdz[0] != 0.0)
	{
		vlen   = o.Gen_Scale(4, ab, cdz[0], v);
		finlen = o.Gen_Sum(finlen, fin[wh], vlen, v, fin[(!wh)]);
		wh     = !wh;
	}

	supo3d2(adx, bdy, adxt_bdyt, cdz, u, finlen, fin, wh, cdy, inv, adxt_cdyt,
	        bdz);
	supo3d2(bdx, cdy, bdxt_cdyt, adz, u, finlen, fin, wh, ady, inv, bdxt_adyt,
	        cdz);
	supo3d2(cdx, ady, cdxt_adyt, bdz, u, finlen, fin, wh, bdy, inv, cdxt_bdyt,
	        adz);

	if (adz[0] != 0.0)
	{
		wlen   = o.Gen_Scale(bctlen, bct, adz[0], w);
		finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[!wh]);
		wh     = !wh;
	}
	if (bdz[0] != 0.0)
	{
		wlen   = o.Gen_Scale(catlen, cat, bdz[0], w);
		finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[!wh]);
		wh     = !wh;
	}
	if (cdz[0] != 0.0)
	{
		wlen   = o.Gen_Scale(abtlen, abt, cdz[0], w);
		finlen = o.Gen_Sum(finlen, fin[wh], wlen, w, fin[!wh]);
		wh     = !wh;
	}

	det = fin[wh][finlen - 1];
	return static_cast<Sign>((det > 0) - (det < 0));
}

inline Sign orient3d(double px, double py, double pz, double qx, double qy,
                     double qz, double rx, double ry, double rz, double sx,
                     double sy, double sz)
{
	Sign ret;
	ret = orient3d_filtered(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
	if (is_sign_reliable(ret))
		return ret;
	return orient3d_expansion(px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz);
}

inline Sign orient3d(const double *p, const double *q, const double *r,
                     const double *s)
{
#ifdef OMC_NO_SHEWCHUK
	return orient3d(p[0], p[1], p[2], q[0], q[1], q[2], r[0], r[1], r[2], s[0],
	                s[1], s[2]);
#else
	return OMC::sign(::orient3d(s, q, r, p));
#endif
}

template <typename IT, typename ET>
Sign orient3d(const GenericPoint3T<IT, ET> &p, const GenericPoint3T<IT, ET> &q,
              const GenericPoint3T<IT, ET> &r, const GenericPoint3T<IT, ET> &s)
{
#ifdef OMC_NO_SHEWCHUK
	return orient3d(p.x(), p.y(), p.z(), q.x(), q.y(), q.z(), r.x(), r.y(), r.z(),
	                s.x(), s.y(), s.z());
#else
	return OMC::sign(::orient3d(s.data(), q.data(), r.data(), p.data()));
#endif
}

inline void orient3d_get_minors(const double *pa, const double *pb,
                                const double *pc, double *minor, double *perm)
{
#ifdef OMC_NO_SHEWCHUK
	const double bax = pb[0] - pa[0];
	const double cax = pc[0] - pa[0];
	const double bay = pb[1] - pa[1];
	const double cay = pc[1] - pa[1];
	const double baz = pb[2] - pa[2];
	const double caz = pc[2] - pa[2];

	const double baxcaz = bax * caz;
	const double caybaz = cay * baz;
	minor[0]            = baxcaz - caybaz;
	perm[0]             = OMC::abs(baxcaz) + OMC::abs(caybaz);

	const double baxcaz = bax * caz;
	const double caxbaz = cax * baz;
	minor[1]            = caxbaz - baxcaz;
	perm[1]             = OMC::abs(caxbaz) + OMC::abs(baxcaz);

	const double baxcay = bax * cay;
	const double caxbax = cax * bax;
	minor[2]            = baxcay - caxbax;
	perm[2]             = OMC::abs(baxcay) + OMC::abs(caxbax);
#else
	::orient3d_get_minors(pb, pc, pa, minor, perm);
#endif
}

inline Sign orient3d_with_cached_minors(const double *pa, const double *pb,
                                        const double *pc, const double *pd,
                                        const double *minor, const double *perm)
{
#ifdef OMC_NO_SHEWCHUK
	const double dax = pd[0] - pa[0];
	const double day = pd[1] - pa[1];
	const double daz = pd[2] - pa[2];

	const double det = dax * minor[0] + day * minor[1] + daz * minor[2];

	const double permanent =
	  OMC::abs(adx) * perm[0] + OMC::abs(ady) * perm[1] + OMC::abs(adz) * perm[2];

	const double errbound = 7.7715611723761027e-16 * permanent;
	if ((det > errbound) || (-det > errbound))
	{
		return OMC::sign(det);
	}

	return orient3d_expansion(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0],
	                          pc[1], pc[2], pd[0], pd[1], pd[2]);
#else
	return OMC::sign(::orient3d_with_cached_minors(
	  pd, pb, pc, pa, const_cast<double *>(minor), const_cast<double *>(perm)));
#endif
}

inline Sign orient2dxy(double p1x, double p1y, double p2x, double p2y,
                       double p3x, double p3y)
{
	Sign ret;
	ret = orient2d_filtered(p1x, p1y, p2x, p2y, p3x, p3y);
	if (is_sign_reliable(ret))
		return ret;
	return orient2d_expansion(p1x, p1y, p2x, p2y, p3x, p3y);
}

inline Sign orient2dxy(const double *p1, const double *p2, const double *p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2d(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]);
#else
	return OMC::sign(::orient2d(p1, p2, p3));
#endif
}

template <typename IT, typename ET>
Sign orient2dxy(const GenericPoint3T<IT, ET> &p1,
                const GenericPoint3T<IT, ET> &p2,
                const GenericPoint3T<IT, ET> &p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2dxy(p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y());
#else
	return OMC::sign(::orient2d(p1.data(), p2.data(), p3.data()));
#endif
}

Sign orient2dyz(double p1y, double p1z, double p2y, double p2z, double p3y,
                double p3z)
{
	Sign ret;
	ret = orient2d_filtered(p1y, p1z, p2y, p2z, p3y, p3z);
	if (is_sign_reliable(ret))
		return ret;
	return orient2d_expansion(p1y, p1z, p2y, p2z, p3y, p3z);
}

inline Sign orient2dyz(const double *p1, const double *p2, const double *p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2dyz(p1[1], p1[2], p2[1], p2[2], p3[1], p3[2]);
#else
	return OMC::sign(::orient2d(p1 + 1, p2 + 1, p3 + 1));
#endif
}

template <typename IT, typename ET>
Sign orient2dyz(const GenericPoint3T<IT, ET> &p1,
                const GenericPoint3T<IT, ET> &p2,
                const GenericPoint3T<IT, ET> &p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2dyz(p1.y(), p1.z(), p2.y(), p2.z(), p3.y(), p3.z());
#else
	return OMC::sign(::orient2d(p1.data() + 1, p2.data() + 1, p3.data() + 1));
#endif
}

Sign orient2dzx(double p1x, double p1z, double p2x, double p2z, double p3x,
                double p3z)
{
	Sign ret;
	ret = orient2d_filtered(p1z, p1x, p2z, p2x, p3z, p3x);
	if (is_sign_reliable(ret))
		return ret;
	return orient2d_expansion(p1z, p1x, p2z, p2x, p3z, p3x);
}

inline Sign orient2dzx(const double *p1, const double *p2, const double *p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2dzx(p1[0], p1[2], p2[0], p2[2], p3[0], p3[2]);
#else
	double _p1[2] = {p1[2], p1[0]};
	double _p2[2] = {p2[2], p2[0]};
	double _p3[2] = {p3[2], p3[0]};
	return OMC::sign(::orient2d(_p1, _p2, _p3));
#endif
}

template <typename IT, typename ET>
Sign orient2dzx(const GenericPoint3T<IT, ET> &p1,
                const GenericPoint3T<IT, ET> &p2,
                const GenericPoint3T<IT, ET> &p3)
{
#ifdef OMC_NO_SHEWCHUK
	return orient2dzx(p1.x(), p1.z(), p2.x(), p2.z(), p3.x(), p3.z());
#else
	double _p1[2] = {p1[2], p1[0]};
	double _p2[2] = {p2[2], p2[0]};
	double _p3[2] = {p3[2], p3[0]};
	return OMC::sign(::orient2d(_p1, _p2, _p3));
#endif
}

inline int maxComponentInTriangleNormal_filtered(double ov1x, double ov1y,
                                                 double ov1z, double ov2x,
                                                 double ov2y, double ov2z,
                                                 double ov3x, double ov3y,
                                                 double ov3z)
{
	double v3x  = ov3x - ov2x;
	double v3y  = ov3y - ov2y;
	double v3z  = ov3z - ov2z;
	double v2x  = ov2x - ov1x;
	double v2y  = ov2y - ov1y;
	double v2z  = ov2z - ov1z;
	double nvx1 = v2y * v3z;
	double nvx2 = v2z * v3y;
	double nvx  = nvx1 - nvx2;
	double nvy1 = v3x * v2z;
	double nvy2 = v3z * v2x;
	double nvy  = nvy1 - nvy2;
	double nvz1 = v2x * v3y;
	double nvz2 = v2y * v3x;
	double nvz  = nvz1 - nvz2;

	double nvxc = fabs(nvx);
	double nvyc = fabs(nvy);
	double nvzc = fabs(nvz);
	double nv   = nvxc;
	int    dim  = 0;
	if (nvyc > nv)
	{
		nv  = nvyc;
		dim = 1;
	}
	if (nvzc > nv)
	{
		nv  = nvzc;
		dim = 2;
	}

	double eps = 8.88720573725927976811e-16, max1, max2;
	if (dim == 0)
	{
		max1 = fabs(v2y) >= fabs(v2z) ? fabs(v2y) : fabs(v2z);
		max2 = fabs(v3y) >= fabs(v3z) ? fabs(v3y) : fabs(v3z);
	}
	else if (dim == 1)
	{
		max1 = fabs(v2x) >= fabs(v2z) ? fabs(v2x) : fabs(v2z);
		max2 = fabs(v3x) >= fabs(v3z) ? fabs(v3x) : fabs(v3z);
	}
	else // dim==2
	{
		max1 = fabs(v2x) >= fabs(v2y) ? fabs(v2x) : fabs(v2y);
		max2 = fabs(v3x) >= fabs(v3y) ? fabs(v3x) : fabs(v3y);
	}
	if (nv > eps * max1 * max2)
		return dim;

	return -1;
}

inline int maxComponentInTriangleNormal_interval(
  IntervalNumber<std::false_type> ov1x, IntervalNumber<std::false_type> ov1y,
  IntervalNumber<std::false_type> ov1z, IntervalNumber<std::false_type> ov2x,
  IntervalNumber<std::false_type> ov2y, IntervalNumber<std::false_type> ov2z,
  IntervalNumber<std::false_type> ov3x, IntervalNumber<std::false_type> ov3y,
  IntervalNumber<std::false_type> ov3z)
{
	using IT = IntervalNumber<std::false_type>;
	IT::Protector P;

	IT v3x  = ov3x - ov2x;
	IT v3y  = ov3y - ov2y;
	IT v3z  = ov3z - ov2z;
	IT v2x  = ov2x - ov1x;
	IT v2y  = ov2y - ov1y;
	IT v2z  = ov2z - ov1z;
	IT nvx1 = v2y * v3z;
	IT nvx2 = v2z * v3y;
	IT nvx  = nvx1 - nvx2;
	IT nvy1 = v3x * v2z;
	IT nvy2 = v3z * v2x;
	IT nvy  = nvy1 - nvy2;
	IT nvz1 = v2x * v3y;
	IT nvz2 = v2y * v3x;
	IT nvz  = nvz1 - nvz2;

	int    dim = -1;
	double nv  = 0.;
	if (nvx.is_sign_reliable())
	{
		nv  = fabs(nvx.inf() + nvx.sup());
		dim = 0;
	}
	if (nvy.is_sign_reliable())
	{
		nv  = std::max(nv, fabs(nvy.inf() + nvy.sup()));
		dim = 1;
	}
	if (nvz.is_sign_reliable())
	{
		nv  = std::max(nv, fabs(nvz.inf() + nvz.sup()));
		dim = 2;
	}

	return dim;
}

inline int maxComponentInTriangleNormal_expansion(double ov1x, double ov1y,
                                                  double ov1z, double ov2x,
                                                  double ov2y, double ov2z,
                                                  double ov3x, double ov3y,
                                                  double ov3z)
{
	expansionObject o;
	double          v3x[2];
	o.two_Diff(ov3x, ov2x, v3x);
	double v3y[2];
	o.two_Diff(ov3y, ov2y, v3y);
	double v3z[2];
	o.two_Diff(ov3z, ov2z, v3z);
	double v2x[2];
	o.two_Diff(ov2x, ov1x, v2x);
	double v2y[2];
	o.two_Diff(ov2y, ov1y, v2y);
	double v2z[2];
	o.two_Diff(ov2z, ov1z, v2z);
	double nvx1[8];
	o.Two_Two_Prod(v2y, v3z, nvx1);
	double nvx2[8];
	o.Two_Two_Prod(v2z, v3y, nvx2);
	double nvx[16];
	int    nvx_len = o.Gen_Diff(8, nvx1, 8, nvx2, nvx);
	double nvy1[8];
	o.Two_Two_Prod(v3x, v2z, nvy1);
	double nvy2[8];
	o.Two_Two_Prod(v3z, v2x, nvy2);
	double nvy[16];
	int    nvy_len = o.Gen_Diff(8, nvy1, 8, nvy2, nvy);
	double nvz1[8];
	o.Two_Two_Prod(v2x, v3y, nvz1);
	double nvz2[8];
	o.Two_Two_Prod(v2y, v3x, nvz2);
	double nvz[16];
	int    nvz_len = o.Gen_Diff(8, nvz1, 8, nvz2, nvz);

	double nvxc = fabs(nvx[nvx_len - 1]);
	double nvyc = fabs(nvy[nvy_len - 1]);
	double nvzc = fabs(nvz[nvz_len - 1]);
	double nv   = nvxc;
	if (nvyc > nv)
		nv = nvyc;
	if (nvzc > nv)
		return 2;
	if (nv == nvxc)
		return 0;
	return 1;
}

inline int maxComponentInTriangleNormal(double ov1x, double ov1y, double ov1z,
                                        double ov2x, double ov2y, double ov2z,
                                        double ov3x, double ov3y, double ov3z)
{
	int ret;
	ret = maxComponentInTriangleNormal_filtered(ov1x, ov1y, ov1z, ov2x, ov2y,
	                                            ov2z, ov3x, ov3y, ov3z);
	if (ret >= 0)
		return ret;
	ret = maxComponentInTriangleNormal_interval(ov1x, ov1y, ov1z, ov2x, ov2y,
	                                            ov2z, ov3x, ov3y, ov3z);
	if (ret >= 0)
		return ret;
	return maxComponentInTriangleNormal_expansion(ov1x, ov1y, ov1z, ov2x, ov2y,
	                                              ov2z, ov3x, ov3y, ov3z);
}

template <typename IT, typename ET>
int longestAxis_IE_filter(const GenericPoint3T<IT, ET> &p1, double bx,
                          double by, double bz, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return -1;

	double tx0 = b1x - bx;
	double tx1 = tx0 * d1;
	double kx  = tx1 + l1x;

	double ty0 = b1y - by;
	double ty1 = ty0 * d1;
	double ky  = ty1 + l1y;

	double tz0 = b1z - bz;
	double tz1 = tz0 * d1;
	double kz  = tz1 + l1z;

	int    dim = 0;
	double k   = fabs(kx);
	if (fabs(ky) > k)
	{
		dim = 1;
		k   = fabs(ky);
	}
	if (fabs(kz) > k)
	{
		dim = 2;
		k   = fabs(kz);
	}

	if (dim == 0)
		max_var = std::max(fabs(tx0), max_var);
	else if (dim == 1)
		max_var = std::max(fabs(ty0), max_var);
	else
		max_var = std::max(fabs(tz0), max_var);

	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.108624468950439e-15;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.2879996548476053e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.680700271819666e-13;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(k, epsilon) != Sign::UNCERTAIN ? dim : -1;
}

template <typename IT, typename ET>
int longestAxis_IE_interval(const GenericPoint3T<IT, ET> &p1, IT bx, IT by,
                            IT bz)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return -1;

	typename IT::Protector P;

	IT tx0 = b1x - bx;
	IT tx1 = tx0 * d1;
	IT kx  = tx1 + l1x;

	IT ty0 = b1y - by;
	IT ty1 = ty0 * d1;
	IT ky  = ty1 + l1y;

	IT tz0 = b1z - bz;
	IT tz1 = tz0 * d1;
	IT kz  = tz1 + l1z;

	int    dim = -1;
	double k   = 0.;
	if (kx.is_sign_reliable())
	{
		k   = fabs(kx.inf() + kx.sup());
		dim = 0;
	}
	if (ky.is_sign_reliable())
	{
		k   = std::max(k, fabs(ky.inf() + ky.sup()));
		dim = 1;
	}
	if (kz.is_sign_reliable())
	{
		k   = std::max(k, fabs(kz.inf() + kz.sup()));
		dim = 2;
	}

	return dim;
}

template <typename IT, typename ET>
int longestAxis_IE_exact(const GenericPoint3T<IT, ET> &p1, ET bx, ET by, ET bz)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);

	ET tx0 = b1x - bx;
	ET tx1 = tx0 * d1;
	ET kx  = tx1 + l1x;

	ET ty0 = b1y - by;
	ET ty1 = ty0 * d1;
	ET ky  = ty1 + l1y;

	ET tz0 = b1z - bz;
	ET tz1 = tz0 * d1;
	ET kz  = tz1 + l1z;

	if (OMC::abs(ky) > OMC::abs(kx))
		return OMC::abs(kz) > ky ? 2 : 1;
	else
		return OMC::abs(kz) > kx ? 2 : 0;
}

template <typename IT, typename ET>
int longestAxis_IE_expansion(const GenericPoint3T<IT, ET> &p1, double bx,
                             double by, double bz)
{
#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	double diff_kx = NAN, diff_ky = NAN, diff_kz = NAN;
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		// x
		double          t0[2];
		o.two_Diff(b1x, bx, t0);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		double k_p[128], *k = k_p;
		int    k_len = o.Gen_Sum_With_PreAlloc(t1_len, t1, l1x_len, l1x, &k, 128);
		diff_kx      = fabs(k[k_len - 1]);
		if (k_p != k)
			FreeDoubles(k);
		if (t1_p != t1)
			FreeDoubles(t1);
		// y
		o.two_Diff(b1y, by, t0);
		t1      = t1_p;
		t1_len  = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		k       = k_p;
		k_len   = o.Gen_Sum_With_PreAlloc(t1_len, t1, l1y_len, l1y, &k, 128);
		diff_ky = fabs(k[k_len - 1]);
		if (k_p != k)
			FreeDoubles(k);
		if (t1_p != t1)
			FreeDoubles(t1);
		// z
		o.two_Diff(b1z, bz, t0);
		t1      = t1_p;
		t1_len  = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		k       = k_p;
		k_len   = o.Gen_Sum_With_PreAlloc(t1_len, t1, l1z_len, l1z, &k, 128);
		diff_kz = fabs(k[k_len - 1]);
		if (k_p != k)
			FreeDoubles(k);
		if (t1_p != t1)
			FreeDoubles(t1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return longestAxis_IE_exact<IT, ET>(p1, bx, by, bz);
#endif

	if (diff_ky > diff_kx)
		return diff_kz > diff_ky ? 2 : 1;
	else
		return diff_kz > diff_kx ? 2 : 0;
}

template <typename IT, typename ET, bool WithSSFilter>
int longestAxis_IE(const GenericPoint3T<IT, ET> &a,
                   const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	int ret;
	if constexpr (WithSSFilter)
	{
		ret = longestAxis_IE_filter<IT, ET>(a, b.x(), b.y(), b.z(), arr);
		if (ret >= 0)
			return ret;
	}
	ret = longestAxis_IE_interval<IT, ET>(a, b.x(), b.y(), b.z());
	if (ret >= 0)
		return ret;
	return longestAxis_IE_expansion<IT, ET>(a, b.x(), b.y(), b.z());
}

template <typename IT, typename ET>
int longestAxis_II_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return -1;

	double tx0 = b1x - b2x;
	double tx1 = tx0 * d1;
	double tx2 = tx1 * d2;
	double tx3 = l1x * d2;
	double tx4 = l2x * d1;
	double tx5 = tx2 + tx3;
	double kx  = tx5 - tx4;

	double ty0 = b1y - b2y;
	double ty1 = ty0 * d1;
	double ty2 = ty1 * d2;
	double ty3 = l1y * d2;
	double ty4 = l2y * d1;
	double ty5 = ty2 + ty3;
	double ky  = ty5 - ty4;

	double tz0 = b1z - b2z;
	double tz1 = tz0 * d1;
	double tz2 = tz1 * d2;
	double tz3 = l1z * d2;
	double tz4 = l2z * d1;
	double tz5 = tz2 + tz3;
	double kz  = tz5 - tz4;

	int    dim = 0;
	double k   = fabs(kx);
	if (fabs(ky) > k)
	{
		dim = 1;
		k   = fabs(ky);
	}
	if (fabs(kz) > k)
	{
		dim = 2;
		k   = fabs(kz);
	}

	if (dim == 0)
		max_var = std::max(fabs(tx0), max_var);
	else if (dim == 1)
		max_var = std::max(fabs(ty0), max_var);
	else
		max_var = std::max(fabs(tz0), max_var);

	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.6431300764452333e-14;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 6.084585960075558e-14;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.6022738691390292e-12;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.20746164403263e-13;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.27253154330999e-12;
	}
	break;
	case PntArr3::TT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.0712852827055069e-10;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(kx, epsilon) != Sign::UNCERTAIN ? dim : -1;
}

template <typename IT, typename ET>
int longestAxis_II_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return -1;

	typename IT::Protector P;

	IT tx0 = b1x - b2x;
	IT tx1 = tx0 * d1;
	IT tx2 = tx1 * d2;
	IT tx3 = l1x * d2;
	IT tx4 = l2x * d1;
	IT tx5 = tx2 + tx3;
	IT kx  = tx5 - tx4;

	IT ty0 = b1y - b2y;
	IT ty1 = ty0 * d1;
	IT ty2 = ty1 * d2;
	IT ty3 = l1y * d2;
	IT ty4 = l2y * d1;
	IT ty5 = ty2 + ty3;
	IT ky  = ty5 - ty4;

	IT tz0 = b1z - b2z;
	IT tz1 = tz0 * d1;
	IT tz2 = tz1 * d2;
	IT tz3 = l1z * d2;
	IT tz4 = l2z * d1;
	IT tz5 = tz2 + tz3;
	IT kz  = tz5 - tz4;

	int    dim = -1;
	double k   = 0.;
	if (kx.is_sign_reliable())
	{
		k   = fabs(kx.inf() + kx.sup());
		dim = 0;
	}
	if (ky.is_sign_reliable())
	{
		k   = std::max(k, fabs(ky.inf() + ky.sup()));
		dim = 1;
	}
	if (kz.is_sign_reliable())
	{
		k   = std::max(k, fabs(kz.inf() + kz.sup()));
		dim = 2;
	}

	return dim;
}

template <typename IT, typename ET>
int longestAxis_II_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);

	ET tx0 = b1x - b2x;
	ET tx1 = tx0 * d1;
	ET tx2 = tx1 * d2;
	ET tx3 = l1x * d2;
	ET tx4 = l2x * d1;
	ET tx5 = tx2 + tx3;
	ET kx  = tx5 - tx4;

	ET ty0 = b1y - b2y;
	ET ty1 = ty0 * d1;
	ET ty2 = ty1 * d2;
	ET ty3 = l1y * d2;
	ET ty4 = l2y * d1;
	ET ty5 = ty2 + ty3;
	ET ky  = ty5 - ty4;

	ET tz0 = b1z - b2z;
	ET tz1 = tz0 * d1;
	ET tz2 = tz1 * d2;
	ET tz3 = l1z * d2;
	ET tz4 = l2z * d1;
	ET tz5 = tz2 + tz3;
	ET kz  = tz5 - tz4;

	if (OMC::abs(ky) > OMC::abs(kx))
		return OMC::abs(kz) > ky ? 2 : 1;
	else
		return OMC::abs(kz) > kx ? 2 : 0;
}

template <typename IT, typename ET>
int longestAxis_II_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2)
{
#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z,
	                   l2x_p[128], *l2x = l2x_p, l2y_p[128], *l2y = l2y_p,
	                   l2z_p[128], *l2z = l2z_p, d2_p[128], *d2 = d2_p, b2x, b2y,
	                   b2z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128, l2x_len = 128,
	    l2y_len = 128, l2z_len = 128, d2_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	double diff_kx = NAN, diff_ky = NAN, diff_kz = NAN;
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		// x
		double          t0[2];
		o.two_Diff(b1x, b2x, t0);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		double t2_p[128], *t2 = t2_p;
		int t2_len = o.Gen_Product_With_PreAlloc(t1_len, t1, d2_len, d2, &t2, 128);
		double t3_p[128], *t3 = t3_p;
		int    t3_len =
		  o.Gen_Product_With_PreAlloc(l1x_len, l1x, d2_len, d2, &t3, 128);
		double t4_p[128], *t4 = t4_p;
		int    t4_len =
		  o.Gen_Product_With_PreAlloc(l2x_len, l2x, d1_len, d1, &t4, 128);
		double t5_p[128], *t5 = t5_p;
		int    t5_len = o.Gen_Sum_With_PreAlloc(t2_len, t2, t3_len, t3, &t5, 128);
		double k_p[128], *k = k_p;
		int    k_len = o.Gen_Diff_With_PreAlloc(t5_len, t5, t4_len, t4, &k, 128);
		diff_kx      = k[k_len - 1];
		if (k_p != k)
			FreeDoubles(k);
		if (t5_p != t5)
			FreeDoubles(t5);
		if (t4_p != t4)
			FreeDoubles(t4);
		if (t3_p != t3)
			FreeDoubles(t3);
		if (t2_p != t2)
			FreeDoubles(t2);
		if (t1_p != t1)
			FreeDoubles(t1);
		// y
		o.two_Diff(b1y, b2y, t0);
		t1      = t1_p;
		t1_len  = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		t2      = t2_p;
		t2_len  = o.Gen_Product_With_PreAlloc(t1_len, t1, d2_len, d2, &t2, 128);
		t3      = t3_p;
		t3_len  = o.Gen_Product_With_PreAlloc(l1y_len, l1y, d2_len, d2, &t3, 128);
		t4      = t4_p;
		t4_len  = o.Gen_Product_With_PreAlloc(l2y_len, l2y, d1_len, d1, &t4, 128);
		t5      = t5_p;
		t5_len  = o.Gen_Sum_With_PreAlloc(t2_len, t2, t3_len, t3, &t5, 128);
		k       = k_p;
		k_len   = o.Gen_Diff_With_PreAlloc(t5_len, t5, t4_len, t4, &k, 128);
		diff_ky = k[k_len - 1];
		if (k_p != k)
			FreeDoubles(k);
		if (t5_p != t5)
			FreeDoubles(t5);
		if (t4_p != t4)
			FreeDoubles(t4);
		if (t3_p != t3)
			FreeDoubles(t3);
		if (t2_p != t2)
			FreeDoubles(t2);
		if (t1_p != t1)
			FreeDoubles(t1);
		// z
		o.two_Diff(b1z, b2z, t0);
		t1      = t1_p;
		t1_len  = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		t2      = t2_p;
		t2_len  = o.Gen_Product_With_PreAlloc(t1_len, t1, d2_len, d2, &t2, 128);
		t3      = t3_p;
		t3_len  = o.Gen_Product_With_PreAlloc(l1z_len, l1z, d2_len, d2, &t3, 128);
		t4      = t4_p;
		t4_len  = o.Gen_Product_With_PreAlloc(l2z_len, l2z, d1_len, d1, &t4, 128);
		t5      = t5_p;
		t5_len  = o.Gen_Sum_With_PreAlloc(t2_len, t2, t3_len, t3, &t5, 128);
		k       = k_p;
		k_len   = o.Gen_Diff_With_PreAlloc(t5_len, t5, t4_len, t4, &k, 128);
		diff_kz = k[k_len - 1];
		if (k_p != k)
			FreeDoubles(k);
		if (t5_p != t5)
			FreeDoubles(t5);
		if (t4_p != t4)
			FreeDoubles(t4);
		if (t3_p != t3)
			FreeDoubles(t3);
		if (t2_p != t2)
			FreeDoubles(t2);
		if (t1_p != t1)
			FreeDoubles(t1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return longestAxis_II_exact<IT, ET>(p1, p2);
#endif

	if (diff_ky > diff_kx)
		return diff_kz > diff_ky ? 2 : 1;
	else
		return diff_kz > diff_kx ? 2 : 0;
}

template <typename IT, typename ET, bool WithSSFilter>
int longestAxis_II(const GenericPoint3T<IT, ET> &a,
                   const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	int ret;
	if constexpr (WithSSFilter)
	{
		ret = longestAxis_II_filtered<IT, ET>(a, b, arr);
		if (ret >= 0)
			return ret;
	}
	ret = longestAxis_II_interval<IT, ET>(a, b);
	if (ret >= 0)
		return ret;
	return longestAxis_II_expansion<IT, ET>(a, b);
}

inline std::array<Sign, 3> lessThanOnAll_EE(double x1, double y1, double z1,
                                            double x2, double y2, double z2)
{
	return std::array<Sign, 3>{static_cast<Sign>((x1 > x2) - (x1 < x2)),
	                           static_cast<Sign>((y1 > y2) - (y1 < y2)),
	                           static_cast<Sign>((z1 > z2) - (z1 < z2))};
}

template <typename IT, typename ET>
std::array<Sign, 3> lessThanOnAll_EE(const GenericPoint3T<IT, ET> &a,
                                     const GenericPoint3T<IT, ET> &b)
{
	return lessThanOnAll_EE(a.x(), a.y(), a.z(), b.x(), b.y(), b.z());
}

template <typename IT, typename ET, bool WithSSFilter>
std::array<Sign, 3> lessThanOnAll_IE(const GenericPoint3T<IT, ET> &p1,
                                     double x2, double y2, double z2,
                                     PntArr3 arr)
{
	Sign retx = lessThanOnX_IE<IT, ET, WithSSFilter>(p1, x2, arr);
	Sign rety = lessThanOnY_IE<IT, ET, WithSSFilter>(p1, y2, arr);
	Sign retz = lessThanOnZ_IE<IT, ET, WithSSFilter>(p1, z2, arr);
	return std::array<Sign, 3>{retx, rety, retz};
}

template <typename IT, typename ET, bool WithSSFilter>
std::array<Sign, 3> lessThanOnAll_IE(const GenericPoint3T<IT, ET> &a,
                                     const GenericPoint3T<IT, ET> &b,
                                     PntArr3                       arr)
{
	return lessThanOnAll_IE<IT, ET, WithSSFilter>(a, b.x(), b.y(), b.z(), arr);
}

template <typename IT, typename ET, bool WithSSFilter>
std::array<Sign, 3> lessThanOnAll_II(const GenericPoint3T<IT, ET> &p1,
                                     const GenericPoint3T<IT, ET> &p2,
                                     PntArr3                       arr)
{
	Sign retx = lessThanOnX_II<IT, ET, WithSSFilter>(p1, p2, arr);
	Sign rety = lessThanOnY_II<IT, ET, WithSSFilter>(p1, p2, arr);
	Sign retz = lessThanOnZ_II<IT, ET, WithSSFilter>(p1, p2, arr);
	return std::array<Sign, 3>{retx, rety, retz};
}

inline Sign lessThan_EE(double x1, double y1, double z1, double x2, double y2,
                        double z2)
{
	int ret = (x1 > x2) - (x1 < x2);
	if (ret)
		return static_cast<Sign>(ret);
	ret = (y1 > y2) - (y1 < y2);
	if (ret)
		return static_cast<Sign>(ret);
	return static_cast<Sign>((z1 > z2) - (z1 < z2));
}

template <typename IT, typename ET>
Sign lessThan_EE(const GenericPoint3T<IT, ET> &a,
                 const GenericPoint3T<IT, ET> &b)
{
	return lessThan_EE(a.x(), a.y(), a.z(), b.x(), b.y(), b.z());
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThan_IE(const GenericPoint3T<IT, ET> &p1, double x2, double y2,
                 double z2, PntArr3 arr)
{
	Sign ret = lessThanOnX_IE<IT, ET, WithSSFilter>(p1, x2, arr);
	if (is_sign_posneg(ret))
		return ret;
	ret = lessThanOnY_IE<IT, ET, WithSSFilter>(p1, y2, arr);
	if (is_sign_posneg(ret))
		return ret;
	return lessThanOnZ_IE<IT, ET, WithSSFilter>(p1, z2, arr);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThan_IE(const GenericPoint3T<IT, ET> &a,
                 const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	return lessThan_IE<IT, ET, WithSSFilter>(a, b.x(), b.y(), b.z(), arr);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThan_II(const GenericPoint3T<IT, ET> &p1,
                 const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	Sign ret;
	ret = lessThanOnX_II<IT, ET, WithSSFilter>(p1, p2, arr);
	if (is_sign_posneg(ret))
		return ret;
	ret = lessThanOnY_II<IT, ET, WithSSFilter>(p1, p2, arr);
	if (is_sign_posneg(ret))
		return ret;
	return lessThanOnZ_II<IT, ET, WithSSFilter>(p1, p2, arr);
}

#define OrientOn2D_LengthThreshold 100

#if defined(INDIRECT_PREDICATES)

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double op3x,
                                double op3y)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p,
	                  l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64],
	                  *d2 = d2_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	bool expansion_calculated = false;
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          a_p[64], *a = a_p;
		int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
		double b_p[64], *b = b_p;
		int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
		double c_p[64], *c = c_p;
		int    c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3y, &c, 64);
		double e_p[64], *e = e_p;
		int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
		double f_p[64], *f = f_p;
		int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
		double g_p[64], *g = g_p;
		int    g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3x, &g, 64);
		double ab_p[64], *ab = ab_p;
		int    ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
		double cd_p[64], *cd = cd_p;
		int    cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1y_len, l1y, &cd, 64);
		double ef_p[64], *ef = ef_p;
		int    ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
		double gh_p[64], *gh = gh_p;
		int    gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1x_len, l1x, &gh, 64);

		double abcd_p[64], *abcd = abcd_p;
		double efgh_p[64], *efgh = efgh_p;
		double L_p[64], *L       = L_p;
		int    abcd_len, efgh_len, L_len;
		if (ab_len * cd_len <= OrientOn2D_LengthThreshold &&
		    ef_len * gh_len <= OrientOn2D_LengthThreshold)
		{
			abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
			efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
			L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);
			return_value         = L[L_len - 1];
			expansion_calculated = true;
		}

		if (L_p != L)
			FreeDoubles(L);
		if (efgh_p != efgh)
			FreeDoubles(efgh);
		if (abcd_p != abcd)
			FreeDoubles(abcd);
		if (gh_p != gh)
			FreeDoubles(gh);
		if (ef_p != ef)
			FreeDoubles(ef);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (g_p != g)
			FreeDoubles(g);
		if (f_p != f)
			FreeDoubles(f);
		if (e_p != e)
			FreeDoubles(e);
		if (c_p != c)
			FreeDoubles(c);
		if (b_p != b)
			FreeDoubles(b);
		if (a_p != a)
			FreeDoubles(a);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dxy_IIE_exact<IT, ET>(p1, p2, op3x, op3y);
	#endif
	if (!expansion_calculated)
		return orientOn2Dxy_IIE_exact<IT, ET>(p1, p2, op3x, op3y);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p,
	                  l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64],
	                  *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64],
	                  *l3y = l3y_p, l3z_p[64], *l3z = l3z_p, d3_p[64], *d3 = d3_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64, l3x_len = 64, l3y_len = 64,
	    l3z_len = 64, d3_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	bool expansion_calculated = false;
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          a_p[64], *a = a_p;
		int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
		double b_p[64], *b = b_p;
		int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
		double c_p[64], *c = c_p;
		int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3y_len, l3y, &c, 64);
		double d_p[64], *d = d_p;
		int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1y_len, l1y, &d, 64);
		double e_p[64], *e = e_p;
		int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
		double f_p[64], *f = f_p;
		int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
		double g_p[64], *g = g_p;
		int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3x_len, l3x, &g, 64);
		double h_p[64], *h = h_p;
		int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1x_len, l1x, &h, 64);
		double ab_p[64], *ab = ab_p;
		int    ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
		double cd_p[64], *cd = cd_p;
		int    cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
		double ef_p[64], *ef = ef_p;
		int    ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
		double gh_p[64], *gh = gh_p;
		int    gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);

		double abcd_p[64], *abcd = abcd_p;
		double efgh_p[64], *efgh = efgh_p;
		double L_p[64], *L       = L_p;
		int    abcd_len, efgh_len, L_len;
		if (ab_len * cd_len <= OrientOn2D_LengthThreshold &&
		    ef_len * gh_len <= OrientOn2D_LengthThreshold)
		{
			abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
			efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
			L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);
			return_value         = L[L_len - 1];
			expansion_calculated = true;
		}

		if (L_p != L)
			FreeDoubles(L);
		if (efgh_p != efgh)
			FreeDoubles(efgh);
		if (abcd_p != abcd)
			FreeDoubles(abcd);
		if (gh_p != gh)
			FreeDoubles(gh);
		if (ef_p != ef)
			FreeDoubles(ef);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (h_p != h)
			FreeDoubles(h);
		if (g_p != g)
			FreeDoubles(g);
		if (f_p != f)
			FreeDoubles(f);
		if (e_p != e)
			FreeDoubles(e);
		if (d_p != d)
			FreeDoubles(d);
		if (c_p != c)
			FreeDoubles(c);
		if (b_p != b)
			FreeDoubles(b);
		if (a_p != a)
			FreeDoubles(a);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (l3z_p != l3z)
			FreeDoubles(l3z);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dxy_III_exact<IT, ET>(p1, p2, p3);
	#endif
	if (!expansion_calculated)
		return orientOn2Dxy_III_exact<IT, ET>(p1, p2, p3);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double op3y,
                                double op3z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p,
	                  l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64],
	                  *d2 = d2_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	bool expansion_calculated = false;
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          a_p[64], *a = a_p;
		int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &a, 64);
		double b_p[64], *b = b_p;
		int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &b, 64);
		double c_p[64], *c = c_p;
		int    c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3z, &c, 64);
		double e_p[64], *e = e_p;
		int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &e, 64);
		double f_p[64], *f = f_p;
		int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &f, 64);
		double g_p[64], *g = g_p;
		int    g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3y, &g, 64);
		double ab_p[64], *ab = ab_p;
		int    ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
		double cd_p[64], *cd = cd_p;
		int    cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1z_len, l1z, &cd, 64);
		double ef_p[64], *ef = ef_p;
		int    ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
		double gh_p[64], *gh = gh_p;
		int    gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1y_len, l1y, &gh, 64);

		double abcd_p[64], *abcd = abcd_p;
		double efgh_p[64], *efgh = efgh_p;
		double L_p[64], *L       = L_p;
		int    abcd_len, efgh_len, L_len;
		if (ab_len * cd_len <= OrientOn2D_LengthThreshold &&
		    ef_len * gh_len <= OrientOn2D_LengthThreshold)
		{
			abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
			efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
			L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);
			return_value         = L[L_len - 1];
			expansion_calculated = true;
		}

		if (L_p != L)
			FreeDoubles(L);
		if (efgh_p != efgh)
			FreeDoubles(efgh);
		if (abcd_p != abcd)
			FreeDoubles(abcd);
		if (gh_p != gh)
			FreeDoubles(gh);
		if (ef_p != ef)
			FreeDoubles(ef);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (g_p != g)
			FreeDoubles(g);
		if (f_p != f)
			FreeDoubles(f);
		if (e_p != e)
			FreeDoubles(e);
		if (c_p != c)
			FreeDoubles(c);
		if (b_p != b)
			FreeDoubles(b);
		if (a_p != a)
			FreeDoubles(a);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dyz_IIE_exact<IT, ET>(p1, p2, op3y, op3z);
	#endif
	if (!expansion_calculated)
		return orientOn2Dyz_IIE_exact<IT, ET>(p1, p2, op3y, op3z);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p,
	                  l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64],
	                  *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64],
	                  *l3y = l3y_p, l3z_p[64], *l3z = l3z_p, d3_p[64], *d3 = d3_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64, l3x_len = 64, l3y_len = 64,
	    l3z_len = 64, d3_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	bool expansion_calculated = false;
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          a_p[64], *a = a_p;
		int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &a, 64);
		double b_p[64], *b = b_p;
		int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &b, 64);
		double c_p[64], *c = c_p;
		int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3z_len, l3z, &c, 64);
		double d_p[64], *d = d_p;
		int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1z_len, l1z, &d, 64);
		double e_p[64], *e = e_p;
		int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &e, 64);
		double f_p[64], *f = f_p;
		int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &f, 64);
		double g_p[64], *g = g_p;
		int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3y_len, l3y, &g, 64);
		double h_p[64], *h = h_p;
		int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1y_len, l1y, &h, 64);
		double ab_p[64], *ab = ab_p;
		int    ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
		double cd_p[64], *cd = cd_p;
		int    cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
		double ef_p[64], *ef = ef_p;
		int    ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
		double gh_p[64], *gh = gh_p;
		int    gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);

		double abcd_p[64], *abcd = abcd_p;
		double efgh_p[64], *efgh = efgh_p;
		double L_p[64], *L       = L_p;
		int    abcd_len, efgh_len, L_len;
		if (ab_len * cd_len <= OrientOn2D_LengthThreshold &&
		    ef_len * gh_len <= OrientOn2D_LengthThreshold)
		{
			abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
			efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
			L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);
			return_value         = L[L_len - 1];
			expansion_calculated = true;
		}

		if (L_p != L)
			FreeDoubles(L);
		if (efgh_p != efgh)
			FreeDoubles(efgh);
		if (abcd_p != abcd)
			FreeDoubles(abcd);
		if (gh_p != gh)
			FreeDoubles(gh);
		if (ef_p != ef)
			FreeDoubles(ef);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (h_p != h)
			FreeDoubles(h);
		if (g_p != g)
			FreeDoubles(g);
		if (f_p != f)
			FreeDoubles(f);
		if (e_p != e)
			FreeDoubles(e);
		if (d_p != d)
			FreeDoubles(d);
		if (c_p != c)
			FreeDoubles(c);
		if (b_p != b)
			FreeDoubles(b);
		if (a_p != a)
			FreeDoubles(a);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (l3z_p != l3z)
			FreeDoubles(l3z);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dyz_III_exact<IT, ET>(p1, p2, p3);
	#endif
	if (!expansion_calculated)
		return orientOn2Dyz_III_exact<IT, ET>(p1, p2, p3);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double op3x,
                                double op3z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p,
	                  l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64],
	                  *d2 = d2_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	bool expansion_calculated = false;
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          a_p[64], *a = a_p;
		int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &a, 64);
		double b_p[64], *b = b_p;
		int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &b, 64);
		double c_p[64], *c = c_p;
		int    c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3x, &c, 64);
		double e_p[64], *e = e_p;
		int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &e, 64);
		double f_p[64], *f = f_p;
		int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &f, 64);
		double g_p[64], *g = g_p;
		int    g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, op3z, &g, 64);
		double ab_p[64], *ab = ab_p;
		int    ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
		double cd_p[64], *cd = cd_p;
		int    cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1x_len, l1x, &cd, 64);
		double ef_p[64], *ef = ef_p;
		int    ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
		double gh_p[64], *gh = gh_p;
		int    gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1z_len, l1z, &gh, 64);

		double abcd_p[64], *abcd = abcd_p;
		double efgh_p[64], *efgh = efgh_p;
		double L_p[64], *L       = L_p;
		int    abcd_len, efgh_len, L_len;
		if (ab_len * cd_len <= OrientOn2D_LengthThreshold &&
		    ef_len * gh_len <= OrientOn2D_LengthThreshold)
		{
			abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
			efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
			L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);
			return_value         = L[L_len - 1];
			expansion_calculated = true;
		}

		if (L_p != L)
			FreeDoubles(L);
		if (efgh_p != efgh)
			FreeDoubles(efgh);
		if (abcd_p != abcd)
			FreeDoubles(abcd);
		if (gh_p != gh)
			FreeDoubles(gh);
		if (ef_p != ef)
			FreeDoubles(ef);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (g_p != g)
			FreeDoubles(g);
		if (f_p != f)
			FreeDoubles(f);
		if (e_p != e)
			FreeDoubles(e);
		if (c_p != c)
			FreeDoubles(c);
		if (b_p != b)
			FreeDoubles(b);
		if (a_p != a)
			FreeDoubles(a);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dzx_IIE_exact<IT, ET>(p1, p2, op3x, op3z);
	#endif
	if (!expansion_calculated)
		return orientOn2Dzx_IIE_exact<IT, ET>(p1, p2, op3x, op3z);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, l2x_p[64], *l2x = l2x_p,
	                  l2y_p[64], *l2y = l2y_p, l2z_p[64], *l2z = l2z_p, d2_p[64],
	                  *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64],
	                  *l3y = l3y_p, l3z_p[64], *l3z = l3z_p, d3_p[64], *d3 = d3_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64, l3x_len = 64, l3y_len = 64,
	    l3z_len = 64, d3_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	bool expansion_calculated = false;
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          a_p[64], *a = a_p;
		int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &a, 64);
		double b_p[64], *b = b_p;
		int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &b, 64);
		double c_p[64], *c = c_p;
		int c_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3x_len, l3x, &c, 64);
		double d_p[64], *d = d_p;
		int d_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1x_len, l1x, &d, 64);
		double e_p[64], *e = e_p;
		int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &e, 64);
		double f_p[64], *f = f_p;
		int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &f, 64);
		double g_p[64], *g = g_p;
		int g_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l3z_len, l3z, &g, 64);
		double h_p[64], *h = h_p;
		int h_len = o.Gen_Product_With_PreAlloc(d3_len, d3, l1z_len, l1z, &h, 64);
		double ab_p[64], *ab = ab_p;
		int    ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
		double cd_p[64], *cd = cd_p;
		int    cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, d_len, d, &cd, 64);
		double ef_p[64], *ef = ef_p;
		int    ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
		double gh_p[64], *gh = gh_p;
		int    gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, h_len, h, &gh, 64);

		double abcd_p[64], *abcd = abcd_p;
		double efgh_p[64], *efgh = efgh_p;
		double L_p[64], *L       = L_p;
		int    abcd_len, efgh_len, L_len;
		if (ab_len * cd_len <= OrientOn2D_LengthThreshold &&
		    ef_len * gh_len <= OrientOn2D_LengthThreshold)
		{
			abcd_len = o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
			efgh_len = o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
			L_len = o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);
			return_value         = L[L_len - 1];
			expansion_calculated = true;
		}

		if (L_p != L)
			FreeDoubles(L);
		if (efgh_p != efgh)
			FreeDoubles(efgh);
		if (abcd_p != abcd)
			FreeDoubles(abcd);
		if (gh_p != gh)
			FreeDoubles(gh);
		if (ef_p != ef)
			FreeDoubles(ef);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (h_p != h)
			FreeDoubles(h);
		if (g_p != g)
			FreeDoubles(g);
		if (f_p != f)
			FreeDoubles(f);
		if (e_p != e)
			FreeDoubles(e);
		if (d_p != d)
			FreeDoubles(d);
		if (c_p != c)
			FreeDoubles(c);
		if (b_p != b)
			FreeDoubles(b);
		if (a_p != a)
			FreeDoubles(a);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (l3z_p != l3z)
			FreeDoubles(l3z);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dzx_III_exact<IT, ET>(p1, p2, p3);
	#endif
	if (!expansion_calculated)
		return orientOn2Dzx_III_exact<IT, ET>(p1, p2, p3);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

#elif defined(OFFSET_PREDICATES)

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double p3x,
                                double p3y)
{
	Sign return_sign = Sign::UNCERTAIN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, b1x, b1y, b1z,
	                  l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64],
	                  *l2z = l2z_p, d2_p[64], *d2 = d2_p, b2x, b2y, b2z;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          b1p3x[2];
		o.two_Diff(b1x, p3x, b1p3x);
		double b1p3y[2];
		o.two_Diff(b1y, p3y, b1p3y);
		double b2p3x[2];
		o.two_Diff(b2x, p3x, b2p3x);
		double b2p3y[2];
		o.two_Diff(b2y, p3y, b2p3y);
		double d1_b1p3x_p[64], *d1_b1p3x = d1_b1p3x_p;
		int    d1_b1p3x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3x, &d1_b1p3x, 64);
		double d1_b1p3y_p[64], *d1_b1p3y = d1_b1p3y_p;
		int    d1_b1p3y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3y, &d1_b1p3y, 64);
		double i1x_p[64], *i1x = i1x_p;
		int    i1x_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3x_len, d1_b1p3x, l1x_len, l1x, &i1x, 64);
		double i1y_p[64], *i1y = i1y_p;
		int    i1y_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3y_len, d1_b1p3y, l1y_len, l1y, &i1y, 64);
		double d2_b2p3x_p[64], *d2_b2p3x = d2_b2p3x_p;
		int    d2_b2p3x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p3x, &d2_b2p3x, 64);
		double d2_b2p3y_p[64], *d2_b2p3y = d2_b2p3y_p;
		int    d2_b2p3y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p3y, &d2_b2p3y, 64);
		double i2x_p[64], *i2x = i2x_p;
		int    i2x_len =
		  o.Gen_Sum_With_PreAlloc(d2_b2p3x_len, d2_b2p3x, l2x_len, l2x, &i2x, 64);
		double i2y_p[64], *i2y = i2y_p;
		int    i2y_len =
		  o.Gen_Sum_With_PreAlloc(d2_b2p3y_len, d2_b2p3y, l2y_len, l2y, &i2y, 64);

		double t0_p[64], *t0 = t0_p;
		int    t0_len;
		double t1_p[64], *t1 = t1_p;
		int    t1_len;
		double det_p[64], *det = det_p;
		int    det_len;

		if (i1x_len * i2y_len <= OrientOn2D_LengthThreshold &&
		    i1y_len * i2x_len <= OrientOn2D_LengthThreshold)
		{
			t0_len = o.Gen_Product_With_PreAlloc(i1x_len, i1x, i2y_len, i2y, &t0, 64);
			t1_len = o.Gen_Product_With_PreAlloc(i1y_len, i1y, i2x_len, i2x, &t1, 64);
			det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 64);
			double return_value = det[det_len - 1];
			if (return_value > 0)
				return_sign = Sign::POSITIVE;
			if (return_value < 0)
				return_sign = Sign::NEGATIVE;
			if (return_value == 0)
				return_sign = Sign::ZERO;
		}
		else
		{
			ET i1x_et = 0, i1y_et = 0, i2x_et = 0, i2y_et = 0;
			// clang-format off
			for (int i = 0; i < i1x_len; i++) i1x_et += i1x[i];
			for (int i = 0; i < i1y_len; i++) i1y_et += i1y[i];
			for (int i = 0; i < i2x_len; i++) i2x_et += i2x[i];
			for (int i = 0; i < i2y_len; i++) i2y_et += i2y[i];
			// clang-format on
			ET t0_et    = i1x_et * i2y_et;
			ET t1_et    = i1y_et * i2x_et;
			ET det_et   = t0_et - t1_et;
			return_sign = OMC::sign(det_et);
		}

		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (i2y_p != i2y)
			FreeDoubles(i2y);
		if (i2x_p != i2x)
			FreeDoubles(i2x);
		if (d2_b2p3y_p != d2_b2p3y)
			FreeDoubles(d2_b2p3y);
		if (d2_b2p3x_p != d2_b2p3x)
			FreeDoubles(d2_b2p3x);
		if (i1y_p != i1y)
			FreeDoubles(i1y);
		if (i1x_p != i1x)
			FreeDoubles(i1x);
		if (d1_b1p3y_p != d1_b1p3y)
			FreeDoubles(d1_b1p3y);
		if (d1_b1p3x_p != d1_b1p3x)
			FreeDoubles(d1_b1p3x);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dxy_IIE_exact<IT, ET>(p1, p2, p3x, p3y);
	#endif

	return return_sign;
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3)
{
	Sign return_sign = Sign::UNCERTAIN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32],
	  *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32], *l1z = l1z_p, d1_p[32],
	  *d1 = d1_p, b1x, b1y, b1z, l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p,
	  l2z_p[32], *l2z = l2z_p, d2_p[32], *d2 = d2_p, b2x, b2y, b2z, l3x_p[32],
	  *l3x = l3x_p, l3y_p[32], *l3y = l3y_p, l3z_p[32], *l3z = l3z_p, d3_p[32],
	  *d3       = d3_p, b3x, b3y, b3z;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32, l3x_len = 32, l3y_len = 32,
	    l3z_len = 32, d3_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len, b3x, b3y, b3z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          b1b3x[2];
		o.two_Diff(b1x, b3x, b1b3x);
		double b1b3y[2];
		o.two_Diff(b1y, b3y, b1b3y);
		double b2b3x[2];
		o.two_Diff(b2x, b3x, b2b3x);
		double b2b3y[2];
		o.two_Diff(b2y, b3y, b2b3y);
		double d1_b1b3x_p[32], *d1_b1b3x = d1_b1b3x_p;
		int    d1_b1b3x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b3x, &d1_b1b3x, 32);
		double d1_b1b3y_p[32], *d1_b1b3y = d1_b1b3y_p;
		int    d1_b1b3y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b3y, &d1_b1b3y, 32);
		double d1_b1b3_l1x_p[32], *d1_b1b3_l1x = d1_b1b3_l1x_p;
		int    d1_b1b3_l1x_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b3x_len, d1_b1b3x, l1x_len, l1x, &d1_b1b3_l1x, 32);
		double d1_b1b3_l1y_p[32], *d1_b1b3_l1y = d1_b1b3_l1y_p;
		int    d1_b1b3_l1y_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b3y_len, d1_b1b3y, l1y_len, l1y, &d1_b1b3_l1y, 32);
		double d3d1_b1b3_l1x_p[32], *d3d1_b1b3_l1x = d3d1_b1b3_l1x_p;
		int    d3d1_b1b3_l1x_len = o.Gen_Product_With_PreAlloc(
      d1_b1b3_l1x_len, d1_b1b3_l1x, d3_len, d3, &d3d1_b1b3_l1x, 32);
		double d3d1_b1b3_l1y_p[32], *d3d1_b1b3_l1y = d3d1_b1b3_l1y_p;
		int    d3d1_b1b3_l1y_len = o.Gen_Product_With_PreAlloc(
      d1_b1b3_l1y_len, d1_b1b3_l1y, d3_len, d3, &d3d1_b1b3_l1y, 32);
		double d2_b2b3x_p[32], *d2_b2b3x = d2_b2b3x_p;
		int    d2_b2b3x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b3x, &d2_b2b3x, 32);
		double d2_b2b3y_p[32], *d2_b2b3y = d2_b2b3y_p;
		int    d2_b2b3y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b3y, &d2_b2b3y, 32);
		double d2_b2b3_l2x_p[32], *d2_b2b3_l2x = d2_b2b3_l2x_p;
		int    d2_b2b3_l2x_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b3x_len, d2_b2b3x, l2x_len, l2x, &d2_b2b3_l2x, 32);
		double d2_b2b3_l2y_p[32], *d2_b2b3_l2y = d2_b2b3_l2y_p;
		int    d2_b2b3_l2y_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b3y_len, d2_b2b3y, l2y_len, l2y, &d2_b2b3_l2y, 32);
		double d3d2_b2b3_l2x_p[32], *d3d2_b2b3_l2x = d3d2_b2b3_l2x_p;
		int    d3d2_b2b3_l2x_len = o.Gen_Product_With_PreAlloc(
      d2_b2b3_l2x_len, d2_b2b3_l2x, d3_len, d3, &d3d2_b2b3_l2x, 32);
		double d3d2_b2b3_l2y_p[32], *d3d2_b2b3_l2y = d3d2_b2b3_l2y_p;
		int    d3d2_b2b3_l2y_len = o.Gen_Product_With_PreAlloc(
      d2_b2b3_l2y_len, d2_b2b3_l2y, d3_len, d3, &d3d2_b2b3_l2y, 32);
		double l3d1x_p[32], *l3d1x = l3d1x_p;
		int    l3d1x_len =
		  o.Gen_Product_With_PreAlloc(l3x_len, l3x, d1_len, d1, &l3d1x, 32);
		double l3d1y_p[32], *l3d1y = l3d1y_p;
		int    l3d1y_len =
		  o.Gen_Product_With_PreAlloc(l3y_len, l3y, d1_len, d1, &l3d1y, 32);
		double l3d2x_p[32], *l3d2x = l3d2x_p;
		int    l3d2x_len =
		  o.Gen_Product_With_PreAlloc(l3x_len, l3x, d2_len, d2, &l3d2x, 32);
		double l3d2y_p[32], *l3d2y = l3d2y_p;
		int    l3d2y_len =
		  o.Gen_Product_With_PreAlloc(l3y_len, l3y, d2_len, d2, &l3d2y, 32);
		double i1x_p[32], *i1x = i1x_p;
		int    i1x_len = o.Gen_Diff_With_PreAlloc(d3d1_b1b3_l1x_len, d3d1_b1b3_l1x,
		                                          l3d1x_len, l3d1x, &i1x, 32);
		double i1y_p[32], *i1y = i1y_p;
		int    i1y_len = o.Gen_Diff_With_PreAlloc(d3d1_b1b3_l1y_len, d3d1_b1b3_l1y,
		                                          l3d1y_len, l3d1y, &i1y, 32);
		double i2x_p[32], *i2x = i2x_p;
		int    i2x_len = o.Gen_Diff_With_PreAlloc(d3d2_b2b3_l2x_len, d3d2_b2b3_l2x,
		                                          l3d2x_len, l3d2x, &i2x, 32);
		double i2y_p[32], *i2y = i2y_p;
		int    i2y_len = o.Gen_Diff_With_PreAlloc(d3d2_b2b3_l2y_len, d3d2_b2b3_l2y,
		                                          l3d2y_len, l3d2y, &i2y, 32);
		double t0_p[32], *t0 = t0_p;
		int    t0_len;
		double t1_p[32], *t1 = t1_p;
		int    t1_len;
		double det_p[32], *det = det_p;
		int    det_len;
		if (i1x_len * i2y_len <= OrientOn2D_LengthThreshold &&
		    i1y_len * i2x_len <= OrientOn2D_LengthThreshold)
		{
			t0_len = o.Gen_Product_With_PreAlloc(i1x_len, i1x, i2y_len, i2y, &t0, 32);
			t1_len = o.Gen_Product_With_PreAlloc(i1y_len, i1y, i2x_len, i2x, &t1, 32);
			det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 32);
			double return_value = det[det_len - 1];
			if (return_value > 0)
				return_sign = Sign::POSITIVE;
			if (return_value < 0)
				return_sign = Sign::NEGATIVE;
			if (return_value == 0)
				return_sign = Sign::ZERO;
		}
		else
		{
			ET i1x_et = 0, i1y_et = 0, i2x_et = 0, i2y_et = 0;
			// clang-format off
			for (int i = 0; i < i1x_len; i++) i1x_et += i1x[i];
			for (int i = 0; i < i1y_len; i++) i1y_et += i1y[i];
			for (int i = 0; i < i2x_len; i++) i2x_et += i2x[i];
			for (int i = 0; i < i2y_len; i++) i2y_et += i2y[i];
			// clang-format on
			ET t0_et    = i1x_et * i2y_et;
			ET t1_et    = i1y_et * i2x_et;
			ET det_et   = t0_et - t1_et;
			return_sign = OMC::sign(det_et);
		}

		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (i2y_p != i2y)
			FreeDoubles(i2y);
		if (i2x_p != i2x)
			FreeDoubles(i2x);
		if (i1y_p != i1y)
			FreeDoubles(i1y);
		if (i1x_p != i1x)
			FreeDoubles(i1x);
		if (l3d2y_p != l3d2y)
			FreeDoubles(l3d2y);
		if (l3d2x_p != l3d2x)
			FreeDoubles(l3d2x);
		if (l3d1y_p != l3d1y)
			FreeDoubles(l3d1y);
		if (l3d1x_p != l3d1x)
			FreeDoubles(l3d1x);
		if (d3d2_b2b3_l2y_p != d3d2_b2b3_l2y)
			FreeDoubles(d3d2_b2b3_l2y);
		if (d3d2_b2b3_l2x_p != d3d2_b2b3_l2x)
			FreeDoubles(d3d2_b2b3_l2x);
		if (d2_b2b3_l2y_p != d2_b2b3_l2y)
			FreeDoubles(d2_b2b3_l2y);
		if (d2_b2b3_l2x_p != d2_b2b3_l2x)
			FreeDoubles(d2_b2b3_l2x);
		if (d2_b2b3y_p != d2_b2b3y)
			FreeDoubles(d2_b2b3y);
		if (d2_b2b3x_p != d2_b2b3x)
			FreeDoubles(d2_b2b3x);
		if (d3d1_b1b3_l1y_p != d3d1_b1b3_l1y)
			FreeDoubles(d3d1_b1b3_l1y);
		if (d3d1_b1b3_l1x_p != d3d1_b1b3_l1x)
			FreeDoubles(d3d1_b1b3_l1x);
		if (d1_b1b3_l1y_p != d1_b1b3_l1y)
			FreeDoubles(d1_b1b3_l1y);
		if (d1_b1b3_l1x_p != d1_b1b3_l1x)
			FreeDoubles(d1_b1b3_l1x);
		if (d1_b1b3y_p != d1_b1b3y)
			FreeDoubles(d1_b1b3y);
		if (d1_b1b3x_p != d1_b1b3x)
			FreeDoubles(d1_b1b3x);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (l3z_p != l3z)
			FreeDoubles(l3z);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dxy_III_exact<IT, ET>(p1, p2, p3);
	#endif

	return return_sign;
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double p3y,
                                double p3z)
{
	Sign return_sign = Sign::UNCERTAIN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, b1x, b1y, b1z,
	                  l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64],
	                  *l2z = l2z_p, d2_p[64], *d2 = d2_p, b2x, b2y, b2z;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          b1p3y[2];
		o.two_Diff(b1y, p3y, b1p3y);
		double b1p3z[2];
		o.two_Diff(b1z, p3z, b1p3z);
		double b2p3y[2];
		o.two_Diff(b2y, p3y, b2p3y);
		double b2p3z[2];
		o.two_Diff(b2z, p3z, b2p3z);
		double d1_b1p3y_p[64], *d1_b1p3y = d1_b1p3y_p;
		int    d1_b1p3y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3y, &d1_b1p3y, 64);
		double d1_b1p3z_p[64], *d1_b1p3z = d1_b1p3z_p;
		int    d1_b1p3z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3z, &d1_b1p3z, 64);
		double i1y_p[64], *i1y = i1y_p;
		int    i1y_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3y_len, d1_b1p3y, l1y_len, l1y, &i1y, 64);
		double i1z_p[64], *i1z = i1z_p;
		int    i1z_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3z_len, d1_b1p3z, l1z_len, l1z, &i1z, 64);
		double d2_b2p3y_p[64], *d2_b2p3y = d2_b2p3y_p;
		int    d2_b2p3y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p3y, &d2_b2p3y, 64);
		double d2_b2p3z_p[64], *d2_b2p3z = d2_b2p3z_p;
		int    d2_b2p3z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p3z, &d2_b2p3z, 64);
		double i2y_p[64], *i2y = i2y_p;
		int    i2y_len =
		  o.Gen_Sum_With_PreAlloc(d2_b2p3y_len, d2_b2p3y, l2y_len, l2y, &i2y, 64);
		double i2z_p[64], *i2z = i2z_p;
		int    i2z_len =
		  o.Gen_Sum_With_PreAlloc(d2_b2p3z_len, d2_b2p3z, l2z_len, l2z, &i2z, 64);

		double t0_p[64], *t0 = t0_p;
		int    t0_len;
		double t1_p[64], *t1 = t1_p;
		int    t1_len;
		double det_p[64], *det = det_p;
		int    det_len;

		if (i1y_len * i2z_len <= OrientOn2D_LengthThreshold &&
		    i1z_len * i2y_len <= OrientOn2D_LengthThreshold)
		{
			t0_len = o.Gen_Product_With_PreAlloc(i1y_len, i1y, i2z_len, i2z, &t0, 64);
			t1_len = o.Gen_Product_With_PreAlloc(i1z_len, i1z, i2y_len, i2y, &t1, 64);
			det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 64);
			double return_value = det[det_len - 1];
			if (return_value > 0)
				return_sign = Sign::POSITIVE;
			if (return_value < 0)
				return_sign = Sign::NEGATIVE;
			if (return_value == 0)
				return_sign = Sign::ZERO;
		}
		else
		{
			ET i1y_et = 0, i1z_et = 0, i2y_et = 0, i2z_et = 0;
			// clang-format off
			for (int i = 0; i < i1y_len; i++) i1y_et += i1y[i];
			for (int i = 0; i < i1z_len; i++) i1z_et += i1z[i];
			for (int i = 0; i < i2y_len; i++) i2y_et += i2y[i];
			for (int i = 0; i < i2z_len; i++) i2z_et += i2z[i];
			// clang-format on
			ET t0_et    = i1y_et * i2z_et;
			ET t1_et    = i1z_et * i2y_et;
			ET det_et   = t0_et - t1_et;
			return_sign = OMC::sign(det_et);
		}

		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (i2z_p != i2z)
			FreeDoubles(i2z);
		if (i2y_p != i2y)
			FreeDoubles(i2y);
		if (d2_b2p3z_p != d2_b2p3z)
			FreeDoubles(d2_b2p3z);
		if (d2_b2p3y_p != d2_b2p3y)
			FreeDoubles(d2_b2p3y);
		if (i1z_p != i1z)
			FreeDoubles(i1z);
		if (i1y_p != i1y)
			FreeDoubles(i1y);
		if (d1_b1p3z_p != d1_b1p3z)
			FreeDoubles(d1_b1p3z);
		if (d1_b1p3y_p != d1_b1p3y)
			FreeDoubles(d1_b1p3y);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dyz_IIE_exact<IT, ET>(p1, p2, p3y, p3z);
	#endif

	return return_sign;
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3)
{
	Sign return_sign = Sign::UNCERTAIN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32],
	  *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32], *l1z = l1z_p, d1_p[32],
	  *d1 = d1_p, b1x, b1y, b1z, l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p,
	  l2z_p[32], *l2z = l2z_p, d2_p[32], *d2 = d2_p, b2x, b2y, b2z, l3x_p[32],
	  *l3x = l3x_p, l3y_p[32], *l3y = l3y_p, l3z_p[32], *l3z = l3z_p, d3_p[32],
	  *d3       = d3_p, b3x, b3y, b3z;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32, l3x_len = 32, l3y_len = 32,
	    l3z_len = 32, d3_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len, b3x, b3y, b3z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          b1b3y[2];
		o.two_Diff(b1y, b3y, b1b3y);
		double b1b3z[2];
		o.two_Diff(b1z, b3z, b1b3z);
		double b2b3y[2];
		o.two_Diff(b2y, b3y, b2b3y);
		double b2b3z[2];
		o.two_Diff(b2z, b3z, b2b3z);
		double d1_b1b3y_p[32], *d1_b1b3y = d1_b1b3y_p;
		int    d1_b1b3y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b3y, &d1_b1b3y, 32);
		double d1_b1b3z_p[32], *d1_b1b3z = d1_b1b3z_p;
		int    d1_b1b3z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b3z, &d1_b1b3z, 32);
		double d1_b1b3_l1y_p[32], *d1_b1b3_l1y = d1_b1b3_l1y_p;
		int    d1_b1b3_l1y_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b3y_len, d1_b1b3y, l1y_len, l1y, &d1_b1b3_l1y, 32);
		double d1_b1b3_l1z_p[32], *d1_b1b3_l1z = d1_b1b3_l1z_p;
		int    d1_b1b3_l1z_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b3z_len, d1_b1b3z, l1z_len, l1z, &d1_b1b3_l1z, 32);
		double d3d1_b1b3_l1y_p[32], *d3d1_b1b3_l1y = d3d1_b1b3_l1y_p;
		int    d3d1_b1b3_l1y_len = o.Gen_Product_With_PreAlloc(
      d1_b1b3_l1y_len, d1_b1b3_l1y, d3_len, d3, &d3d1_b1b3_l1y, 32);
		double d3d1_b1b3_l1z_p[32], *d3d1_b1b3_l1z = d3d1_b1b3_l1z_p;
		int    d3d1_b1b3_l1z_len = o.Gen_Product_With_PreAlloc(
      d1_b1b3_l1z_len, d1_b1b3_l1z, d3_len, d3, &d3d1_b1b3_l1z, 32);
		double d2_b2b3y_p[32], *d2_b2b3y = d2_b2b3y_p;
		int    d2_b2b3y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b3y, &d2_b2b3y, 32);
		double d2_b2b3z_p[32], *d2_b2b3z = d2_b2b3z_p;
		int    d2_b2b3z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b3z, &d2_b2b3z, 32);
		double d2_b2b3_l2y_p[32], *d2_b2b3_l2y = d2_b2b3_l2y_p;
		int    d2_b2b3_l2y_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b3y_len, d2_b2b3y, l2y_len, l2y, &d2_b2b3_l2y, 32);
		double d2_b2b3_l2z_p[32], *d2_b2b3_l2z = d2_b2b3_l2z_p;
		int    d2_b2b3_l2z_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b3z_len, d2_b2b3z, l2z_len, l2z, &d2_b2b3_l2z, 32);
		double d3d2_b2b3_l2y_p[32], *d3d2_b2b3_l2y = d3d2_b2b3_l2y_p;
		int    d3d2_b2b3_l2y_len = o.Gen_Product_With_PreAlloc(
      d2_b2b3_l2y_len, d2_b2b3_l2y, d3_len, d3, &d3d2_b2b3_l2y, 32);
		double d3d2_b2b3_l2z_p[32], *d3d2_b2b3_l2z = d3d2_b2b3_l2z_p;
		int    d3d2_b2b3_l2z_len = o.Gen_Product_With_PreAlloc(
      d2_b2b3_l2z_len, d2_b2b3_l2z, d3_len, d3, &d3d2_b2b3_l2z, 32);
		double l3d1y_p[32], *l3d1y = l3d1y_p;
		int    l3d1y_len =
		  o.Gen_Product_With_PreAlloc(l3y_len, l3y, d1_len, d1, &l3d1y, 32);
		double l3d1z_p[32], *l3d1z = l3d1z_p;
		int    l3d1z_len =
		  o.Gen_Product_With_PreAlloc(l3z_len, l3z, d1_len, d1, &l3d1z, 32);
		double l3d2y_p[32], *l3d2y = l3d2y_p;
		int    l3d2y_len =
		  o.Gen_Product_With_PreAlloc(l3y_len, l3y, d2_len, d2, &l3d2y, 32);
		double l3d2z_p[32], *l3d2z = l3d2z_p;
		int    l3d2z_len =
		  o.Gen_Product_With_PreAlloc(l3z_len, l3z, d2_len, d2, &l3d2z, 32);
		double i1y_p[32], *i1y = i1y_p;
		int    i1y_len = o.Gen_Diff_With_PreAlloc(d3d1_b1b3_l1y_len, d3d1_b1b3_l1y,
		                                          l3d1y_len, l3d1y, &i1y, 32);
		double i1z_p[32], *i1z = i1z_p;
		int    i1z_len = o.Gen_Diff_With_PreAlloc(d3d1_b1b3_l1z_len, d3d1_b1b3_l1z,
		                                          l3d1z_len, l3d1z, &i1z, 32);
		double i2y_p[32], *i2y = i2y_p;
		int    i2y_len = o.Gen_Diff_With_PreAlloc(d3d2_b2b3_l2y_len, d3d2_b2b3_l2y,
		                                          l3d2y_len, l3d2y, &i2y, 32);
		double i2z_p[32], *i2z = i2z_p;
		int    i2z_len = o.Gen_Diff_With_PreAlloc(d3d2_b2b3_l2z_len, d3d2_b2b3_l2z,
		                                          l3d2z_len, l3d2z, &i2z, 32);
		double t0_p[32], *t0 = t0_p;
		int    t0_len;
		double t1_p[32], *t1 = t1_p;
		int    t1_len;
		double det_p[32], *det = det_p;
		int    det_len;

		if (i1y_len * i2z_len <= OrientOn2D_LengthThreshold &&
		    i1z_len * i2y_len <= OrientOn2D_LengthThreshold)
		{
			t0_len = o.Gen_Product_With_PreAlloc(i1y_len, i1y, i2z_len, i2z, &t0, 32);
			t1_len = o.Gen_Product_With_PreAlloc(i1z_len, i1z, i2y_len, i2y, &t1, 32);
			det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 32);
			double return_value = det[det_len - 1];
			if (return_value > 0)
				return_sign = Sign::POSITIVE;
			if (return_value < 0)
				return_sign = Sign::NEGATIVE;
			if (return_value == 0)
				return_sign = Sign::ZERO;
		}
		else
		{
			ET i1y_et = 0, i1z_et = 0, i2y_et = 0, i2z_et = 0;
			// clang-format off
			for (int i = 0; i < i1y_len; i++) i1y_et += i1y[i];
			for (int i = 0; i < i1z_len; i++) i1z_et += i1z[i];
			for (int i = 0; i < i2y_len; i++) i2y_et += i2y[i];
			for (int i = 0; i < i2z_len; i++) i2z_et += i2z[i];
			// clang-format on
			ET t0_et    = i1y_et * i2z_et;
			ET t1_et    = i1z_et * i2y_et;
			ET det_et   = t0_et - t1_et;
			return_sign = OMC::sign(det_et);
		}

		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (i2z_p != i2z)
			FreeDoubles(i2z);
		if (i2y_p != i2y)
			FreeDoubles(i2y);
		if (i1z_p != i1z)
			FreeDoubles(i1z);
		if (i1y_p != i1y)
			FreeDoubles(i1y);
		if (l3d2z_p != l3d2z)
			FreeDoubles(l3d2z);
		if (l3d2y_p != l3d2y)
			FreeDoubles(l3d2y);
		if (l3d1z_p != l3d1z)
			FreeDoubles(l3d1z);
		if (l3d1y_p != l3d1y)
			FreeDoubles(l3d1y);
		if (d3d2_b2b3_l2z_p != d3d2_b2b3_l2z)
			FreeDoubles(d3d2_b2b3_l2z);
		if (d3d2_b2b3_l2y_p != d3d2_b2b3_l2y)
			FreeDoubles(d3d2_b2b3_l2y);
		if (d2_b2b3_l2z_p != d2_b2b3_l2z)
			FreeDoubles(d2_b2b3_l2z);
		if (d2_b2b3_l2y_p != d2_b2b3_l2y)
			FreeDoubles(d2_b2b3_l2y);
		if (d2_b2b3z_p != d2_b2b3z)
			FreeDoubles(d2_b2b3z);
		if (d2_b2b3y_p != d2_b2b3y)
			FreeDoubles(d2_b2b3y);
		if (d3d1_b1b3_l1z_p != d3d1_b1b3_l1z)
			FreeDoubles(d3d1_b1b3_l1z);
		if (d3d1_b1b3_l1y_p != d3d1_b1b3_l1y)
			FreeDoubles(d3d1_b1b3_l1y);
		if (d1_b1b3_l1z_p != d1_b1b3_l1z)
			FreeDoubles(d1_b1b3_l1z);
		if (d1_b1b3_l1y_p != d1_b1b3_l1y)
			FreeDoubles(d1_b1b3_l1y);
		if (d1_b1b3z_p != d1_b1b3z)
			FreeDoubles(d1_b1b3z);
		if (d1_b1b3y_p != d1_b1b3y)
			FreeDoubles(d1_b1b3y);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (l3z_p != l3z)
			FreeDoubles(l3z);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dyz_III_exact<IT, ET>(p1, p2, p3);
	#endif

	return return_sign;
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double p3x,
                                double p3z)
{
	Sign return_sign = Sign::UNCERTAIN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, b1x, b1y, b1z,
	                  l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, l2z_p[64],
	                  *l2z = l2z_p, d2_p[64], *d2 = d2_p, b2x, b2y, b2z;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64, l2x_len = 64,
	    l2y_len = 64, l2z_len = 64, d2_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          b1p3z[2];
		o.two_Diff(b1z, p3z, b1p3z);
		double b1p3x[2];
		o.two_Diff(b1x, p3x, b1p3x);
		double b2p3z[2];
		o.two_Diff(b2z, p3z, b2p3z);
		double b2p3x[2];
		o.two_Diff(b2x, p3x, b2p3x);
		double d1_b1p3z_p[64], *d1_b1p3z = d1_b1p3z_p;
		int    d1_b1p3z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3z, &d1_b1p3z, 64);
		double d1_b1p3x_p[64], *d1_b1p3x = d1_b1p3x_p;
		int    d1_b1p3x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3x, &d1_b1p3x, 64);
		double i1z_p[64], *i1z = i1z_p;
		int    i1z_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3z_len, d1_b1p3z, l1z_len, l1z, &i1z, 64);
		double i1x_p[64], *i1x = i1x_p;
		int    i1x_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3x_len, d1_b1p3x, l1x_len, l1x, &i1x, 64);
		double d2_b2p3z_p[64], *d2_b2p3z = d2_b2p3z_p;
		int    d2_b2p3z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p3z, &d2_b2p3z, 64);
		double d2_b2p3x_p[64], *d2_b2p3x = d2_b2p3x_p;
		int    d2_b2p3x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p3x, &d2_b2p3x, 64);
		double i2z_p[64], *i2z = i2z_p;
		int    i2z_len =
		  o.Gen_Sum_With_PreAlloc(d2_b2p3z_len, d2_b2p3z, l2z_len, l2z, &i2z, 64);
		double i2x_p[64], *i2x = i2x_p;
		int    i2x_len =
		  o.Gen_Sum_With_PreAlloc(d2_b2p3x_len, d2_b2p3x, l2x_len, l2x, &i2x, 64);

		double t0_p[64], *t0 = t0_p;
		int    t0_len;
		double t1_p[64], *t1 = t1_p;
		int    t1_len;
		double det_p[64], *det = det_p;
		int    det_len;

		if (i1z_len * i2x_len <= OrientOn2D_LengthThreshold &&
		    i1x_len * i2z_len <= OrientOn2D_LengthThreshold)
		{
			t0_len = o.Gen_Product_With_PreAlloc(i1z_len, i1z, i2x_len, i2x, &t0, 64);
			t1_len = o.Gen_Product_With_PreAlloc(i1x_len, i1x, i2z_len, i2z, &t1, 64);
			det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 64);
			double return_value = det[det_len - 1];
			if (return_value > 0)
				return_sign = Sign::POSITIVE;
			if (return_value < 0)
				return_sign = Sign::NEGATIVE;
			if (return_value == 0)
				return_sign = Sign::ZERO;
		}
		else
		{
			ET i1x_et = 0, i1z_et = 0, i2x_et = 0, i2z_et = 0;
			// clang-format off
			for (int i = 0; i < i1x_len; i++) i1x_et += i1x[i];
			for (int i = 0; i < i1z_len; i++) i1z_et += i1z[i];
			for (int i = 0; i < i2x_len; i++) i2x_et += i2x[i];
			for (int i = 0; i < i2z_len; i++) i2z_et += i2z[i];
			// clang-format on
			ET t0_et    = i1z_et * i2x_et;
			ET t1_et    = i1x_et * i2z_et;
			ET det_et   = t0_et - t1_et;
			return_sign = OMC::sign(det_et);
		}

		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (i2x_p != i2x)
			FreeDoubles(i2x);
		if (i2z_p != i2z)
			FreeDoubles(i2z);
		if (d2_b2p3x_p != d2_b2p3x)
			FreeDoubles(d2_b2p3x);
		if (d2_b2p3z_p != d2_b2p3z)
			FreeDoubles(d2_b2p3z);
		if (i1x_p != i1x)
			FreeDoubles(i1x);
		if (i1z_p != i1z)
			FreeDoubles(i1z);
		if (d1_b1p3x_p != d1_b1p3x)
			FreeDoubles(d1_b1p3x);
		if (d1_b1p3z_p != d1_b1p3z)
			FreeDoubles(d1_b1p3z);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dzx_IIE_exact<IT, ET>(p1, p2, p3x, p3z);
	#endif

	return return_sign;
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3)
{
	Sign return_sign = Sign::UNCERTAIN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32],
	  *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32], *l1z = l1z_p, d1_p[32],
	  *d1 = d1_p, b1x, b1y, b1z, l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p,
	  l2z_p[32], *l2z = l2z_p, d2_p[32], *d2 = d2_p, b2x, b2y, b2z, l3x_p[32],
	  *l3x = l3x_p, l3y_p[32], *l3y = l3y_p, l3z_p[32], *l3z = l3z_p, d3_p[32],
	  *d3       = d3_p, b3x, b3y, b3z;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32, l3x_len = 32, l3y_len = 32,
	    l3z_len = 32, d3_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len, b3x, b3y, b3z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          b1b3z[2];
		o.two_Diff(b1z, b3z, b1b3z);
		double b1b3x[2];
		o.two_Diff(b1x, b3x, b1b3x);
		double b2b3z[2];
		o.two_Diff(b2z, b3z, b2b3z);
		double b2b3x[2];
		o.two_Diff(b2x, b3x, b2b3x);
		double d1_b1b3z_p[32], *d1_b1b3z = d1_b1b3z_p;
		int    d1_b1b3z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b3z, &d1_b1b3z, 32);
		double d1_b1b3x_p[32], *d1_b1b3x = d1_b1b3x_p;
		int    d1_b1b3x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b3x, &d1_b1b3x, 32);
		double d1_b1b3_l1z_p[32], *d1_b1b3_l1z = d1_b1b3_l1z_p;
		int    d1_b1b3_l1z_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b3z_len, d1_b1b3z, l1z_len, l1z, &d1_b1b3_l1z, 32);
		double d1_b1b3_l1x_p[32], *d1_b1b3_l1x = d1_b1b3_l1x_p;
		int    d1_b1b3_l1x_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b3x_len, d1_b1b3x, l1x_len, l1x, &d1_b1b3_l1x, 32);
		double d3d1_b1b3_l1z_p[32], *d3d1_b1b3_l1z = d3d1_b1b3_l1z_p;
		int    d3d1_b1b3_l1z_len = o.Gen_Product_With_PreAlloc(
      d1_b1b3_l1z_len, d1_b1b3_l1z, d3_len, d3, &d3d1_b1b3_l1z, 32);
		double d3d1_b1b3_l1x_p[32], *d3d1_b1b3_l1x = d3d1_b1b3_l1x_p;
		int    d3d1_b1b3_l1x_len = o.Gen_Product_With_PreAlloc(
      d1_b1b3_l1x_len, d1_b1b3_l1x, d3_len, d3, &d3d1_b1b3_l1x, 32);
		double d2_b2b3z_p[32], *d2_b2b3z = d2_b2b3z_p;
		int    d2_b2b3z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b3z, &d2_b2b3z, 32);
		double d2_b2b3x_p[32], *d2_b2b3x = d2_b2b3x_p;
		int    d2_b2b3x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b3x, &d2_b2b3x, 32);
		double d2_b2b3_l2z_p[32], *d2_b2b3_l2z = d2_b2b3_l2z_p;
		int    d2_b2b3_l2z_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b3z_len, d2_b2b3z, l2z_len, l2z, &d2_b2b3_l2z, 32);
		double d2_b2b3_l2x_p[32], *d2_b2b3_l2x = d2_b2b3_l2x_p;
		int    d2_b2b3_l2x_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b3x_len, d2_b2b3x, l2x_len, l2x, &d2_b2b3_l2x, 32);
		double d3d2_b2b3_l2z_p[32], *d3d2_b2b3_l2z = d3d2_b2b3_l2z_p;
		int    d3d2_b2b3_l2z_len = o.Gen_Product_With_PreAlloc(
      d2_b2b3_l2z_len, d2_b2b3_l2z, d3_len, d3, &d3d2_b2b3_l2z, 32);
		double d3d2_b2b3_l2x_p[32], *d3d2_b2b3_l2x = d3d2_b2b3_l2x_p;
		int    d3d2_b2b3_l2x_len = o.Gen_Product_With_PreAlloc(
      d2_b2b3_l2x_len, d2_b2b3_l2x, d3_len, d3, &d3d2_b2b3_l2x, 32);
		double l3d1z_p[32], *l3d1z = l3d1z_p;
		int    l3d1z_len =
		  o.Gen_Product_With_PreAlloc(l3z_len, l3z, d1_len, d1, &l3d1z, 32);
		double l3d1x_p[32], *l3d1x = l3d1x_p;
		int    l3d1x_len =
		  o.Gen_Product_With_PreAlloc(l3x_len, l3x, d1_len, d1, &l3d1x, 32);
		double l3d2z_p[32], *l3d2z = l3d2z_p;
		int    l3d2z_len =
		  o.Gen_Product_With_PreAlloc(l3z_len, l3z, d2_len, d2, &l3d2z, 32);
		double l3d2x_p[32], *l3d2x = l3d2x_p;
		int    l3d2x_len =
		  o.Gen_Product_With_PreAlloc(l3x_len, l3x, d2_len, d2, &l3d2x, 32);
		double i1z_p[32], *i1z = i1z_p;
		int    i1z_len = o.Gen_Diff_With_PreAlloc(d3d1_b1b3_l1z_len, d3d1_b1b3_l1z,
		                                          l3d1z_len, l3d1z, &i1z, 32);
		double i1x_p[32], *i1x = i1x_p;
		int    i1x_len = o.Gen_Diff_With_PreAlloc(d3d1_b1b3_l1x_len, d3d1_b1b3_l1x,
		                                          l3d1x_len, l3d1x, &i1x, 32);
		double i2z_p[32], *i2z = i2z_p;
		int    i2z_len = o.Gen_Diff_With_PreAlloc(d3d2_b2b3_l2z_len, d3d2_b2b3_l2z,
		                                          l3d2z_len, l3d2z, &i2z, 32);
		double i2x_p[32], *i2x = i2x_p;
		int    i2x_len = o.Gen_Diff_With_PreAlloc(d3d2_b2b3_l2x_len, d3d2_b2b3_l2x,
		                                          l3d2x_len, l3d2x, &i2x, 32);
		double t0_p[32], *t0 = t0_p;
		int    t0_len;
		double t1_p[32], *t1 = t1_p;
		int    t1_len;
		double det_p[32], *det = det_p;
		int    det_len;

		if (i1z_len * i2x_len <= OrientOn2D_LengthThreshold &&
		    i1x_len * i2z_len <= OrientOn2D_LengthThreshold)
		{
			t0_len = o.Gen_Product_With_PreAlloc(i1z_len, i1z, i2x_len, i2x, &t0, 32);
			t1_len = o.Gen_Product_With_PreAlloc(i1x_len, i1x, i2z_len, i2z, &t1, 32);
			det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 32);
			double return_value = det[det_len - 1];
			if (return_value > 0)
				return_sign = Sign::POSITIVE;
			if (return_value < 0)
				return_sign = Sign::NEGATIVE;
			if (return_value == 0)
				return_sign = Sign::ZERO;
		}
		else
		{
			ET i1x_et = 0, i1z_et = 0, i2x_et = 0, i2z_et = 0;
			// clang-format off
			for (int i = 0; i < i1x_len; i++) i1x_et += i1x[i];
			for (int i = 0; i < i1z_len; i++) i1z_et += i1z[i];
			for (int i = 0; i < i2x_len; i++) i2x_et += i2x[i];
			for (int i = 0; i < i2z_len; i++) i2z_et += i2z[i];
			// clang-format on
			ET t0_et    = i1z_et * i2x_et;
			ET t1_et    = i1x_et * i2z_et;
			ET det_et   = t0_et - t1_et;
			return_sign = OMC::sign(det_et);
		}

		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (i2x_p != i2x)
			FreeDoubles(i2x);
		if (i2z_p != i2z)
			FreeDoubles(i2z);
		if (i1x_p != i1x)
			FreeDoubles(i1x);
		if (i1z_p != i1z)
			FreeDoubles(i1z);
		if (l3d2x_p != l3d2x)
			FreeDoubles(l3d2x);
		if (l3d2z_p != l3d2z)
			FreeDoubles(l3d2z);
		if (l3d1x_p != l3d1x)
			FreeDoubles(l3d1x);
		if (l3d1z_p != l3d1z)
			FreeDoubles(l3d1z);
		if (d3d2_b2b3_l2x_p != d3d2_b2b3_l2x)
			FreeDoubles(d3d2_b2b3_l2x);
		if (d3d2_b2b3_l2z_p != d3d2_b2b3_l2z)
			FreeDoubles(d3d2_b2b3_l2z);
		if (d2_b2b3_l2x_p != d2_b2b3_l2x)
			FreeDoubles(d2_b2b3_l2x);
		if (d2_b2b3_l2z_p != d2_b2b3_l2z)
			FreeDoubles(d2_b2b3_l2z);
		if (d2_b2b3x_p != d2_b2b3x)
			FreeDoubles(d2_b2b3x);
		if (d2_b2b3z_p != d2_b2b3z)
			FreeDoubles(d2_b2b3z);
		if (d3d1_b1b3_l1x_p != d3d1_b1b3_l1x)
			FreeDoubles(d3d1_b1b3_l1x);
		if (d3d1_b1b3_l1z_p != d3d1_b1b3_l1z)
			FreeDoubles(d3d1_b1b3_l1z);
		if (d1_b1b3_l1x_p != d1_b1b3_l1x)
			FreeDoubles(d1_b1b3_l1x);
		if (d1_b1b3_l1z_p != d1_b1b3_l1z)
			FreeDoubles(d1_b1b3_l1z);
		if (d1_b1b3x_p != d1_b1b3x)
			FreeDoubles(d1_b1b3x);
		if (d1_b1b3z_p != d1_b1b3z)
			FreeDoubles(d1_b1b3z);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (l1z_p != l1z)
			FreeDoubles(l1z);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (l2z_p != l2z)
			FreeDoubles(l2z);
		if (d2_p != d2)
			FreeDoubles(d2);
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (l3z_p != l3z)
			FreeDoubles(l3z);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dzx_III_exact<IT, ET>(p1, p2, p3);
	#endif

	return return_sign;
}

#endif

} // namespace OMC