#pragma once
#include "ImplicitPointPredicates.h"
#include "IndirectDefinitions.h"

#include "OpenMeshCraft/NumberTypes/ExpansionObject.h"
#include "OpenMeshCraft/NumberTypes/IntervalNumber.h"

namespace OMC {

#if defined(OMC_INDIRECT_PRED)

bool lambda2d_SSI_filtered(double ea1x, double ea1y, double ea2x, double ea2y,
                           double eb1x, double eb1y, double eb2x, double eb2y,
                           double &lambda_x, double &lambda_y,
                           double &lambda_det, double &max_var)
{
	double t1a  = ea1x * ea2y;
	double t1b  = ea2x * ea1y;
	double t1   = t1a - t1b;
	double tx2  = eb1x - eb2x;
	double t3a  = eb1x * eb2y;
	double t3b  = eb2x * eb1y;
	double t3   = t3a - t3b;
	double tx4  = ea1x - ea2x;
	double ty2  = eb1y - eb2y;
	double ty4  = ea1y - ea2y;
	double lxa  = t1 * tx2;
	double lxb  = t3 * tx4;
	lambda_x    = lxa - lxb;
	double lya  = t1 * ty2;
	double lyb  = t3 * ty4;
	lambda_y    = lya - lyb;
	double deta = tx4 * ty2;
	double detb = tx2 * ty4;
	lambda_det  = deta - detb;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(ea1x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ea1y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ea2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ea2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(eb1x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(eb1y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(eb2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(eb2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(tx2)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(tx4)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ty2)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ty4)) > max_var)
		max_var = _tmp_fabs;
	double lambda_det_eps = max_var;
	lambda_det_eps *= lambda_det_eps;
	lambda_det_eps *= 8.881784197001253e-16;

	return ((lambda_det > lambda_det_eps || lambda_det < -lambda_det_eps));
}

bool lambda3d_SSI_filtered(double xa, double ya, double za, double xb,
                           double yb, double zb, double xp, double yp,
                           double xq, double yq, double &lambda_x,
                           double &lambda_y, double &lambda_z, double &lambda_d,
                           double &max_var)
{
	double xba = xb - xa;
	double yba = yb - ya;
	double zba = zb - za;
	double xap = xa - xp;
	double yap = ya - yp;
	double yqp = yq - yp;
	double xqp = xq - xp;
	double c1  = xap * yqp;
	double c2  = xqp * yap;
	double c3  = xba * yqp;
	double c4  = xqp * yba;
	double c12 = c1 - c2;
	lambda_d   = c4 - c3;
	double xf  = xba * c12;
	double xs  = xa * lambda_d;
	lambda_x   = xf + xs;
	double yf  = yba * c12;
	double ys  = ya * lambda_d;
	lambda_y   = yf + ys;
	double zf  = zba * c12;
	double zs  = za * lambda_d;
	lambda_z   = zf + zs;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(xa)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ya)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(za)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xap)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yap)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yqp)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xqp)) > max_var)
		max_var = _tmp_fabs;
	double lambda_d_eps = max_var;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= 8.881784197001252e-16;

	return (lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps);
}

bool lambda3d_LPI_filtered(double px, double py, double pz, double qx,
                           double qy, double qz, double rx, double ry,
                           double rz, double sx, double sy, double sz,
                           double tx, double ty, double tz, double &lambda_d,
                           double &lambda_x, double &lambda_y, double &lambda_z,
                           double &max_var)
{
	double a11   = px - qx;
	double a12   = py - qy;
	double a13   = pz - qz;
	double a21   = sx - rx;
	double a22   = sy - ry;
	double a23   = sz - rz;
	double a31   = tx - rx;
	double a32   = ty - ry;
	double a33   = tz - rz;
	double tv1   = a22 * a33;
	double tv2   = a23 * a32;
	double a2233 = tv1 - tv2;
	double tv3   = a21 * a33;
	double tv4   = a23 * a31;
	double a2133 = tv3 - tv4;
	double tv5   = a21 * a32;
	double tv6   = a22 * a31;
	double a2132 = tv5 - tv6;
	double tv7   = a11 * a2233;
	double tv8   = a12 * a2133;
	double tv9   = a13 * a2132;
	double tt1   = tv7 - tv8;
	lambda_d     = tt1 + tv9;
	double px_rx = px - rx;
	double py_ry = py - ry;
	double pz_rz = pz - rz;
	double tt2   = py_ry * a2133;
	double tt3   = px_rx * a2233;
	double tt4   = pz_rz * a2132;
	double tt5   = tt3 + tt4;
	double n     = tt5 - tt2;
	double ax    = a11 * n;
	double ay    = a12 * n;
	double az    = a13 * n;
	double dpx   = lambda_d * px;
	double dpy   = lambda_d * py;
	double dpz   = lambda_d * pz;
	lambda_x     = dpx - ax;
	lambda_y     = dpy - ay;
	lambda_z     = dpz - az;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(px)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(py)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(pz)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a11)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a12)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a13)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a21)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a22)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a23)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a31)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a32)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(a33)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(px_rx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(py_ry)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(pz_rz)) > max_var)
		max_var = _tmp_fabs;
	double lambda_d_eps = max_var;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= max_var;
	lambda_d_eps *= 4.884981308350689e-15;

	return ((lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps));
}

bool lambda3d_TPI_filtered(double ov1x, double ov1y, double ov1z, double ov2x,
                           double ov2y, double ov2z, double ov3x, double ov3y,
                           double ov3z, double ow1x, double ow1y, double ow1z,
                           double ow2x, double ow2y, double ow2z, double ow3x,
                           double ow3y, double ow3z, double ou1x, double ou1y,
                           double ou1z, double ou2x, double ou2y, double ou2z,
                           double ou3x, double ou3y, double ou3z,
                           double &lambda_x, double &lambda_y, double &lambda_z,
                           double &lambda_d, double &max_var)
{
	double v3x    = ov3x - ov2x;
	double v3y    = ov3y - ov2y;
	double v3z    = ov3z - ov2z;
	double v2x    = ov2x - ov1x;
	double v2y    = ov2y - ov1y;
	double v2z    = ov2z - ov1z;
	double w3x    = ow3x - ow2x;
	double w3y    = ow3y - ow2y;
	double w3z    = ow3z - ow2z;
	double w2x    = ow2x - ow1x;
	double w2y    = ow2y - ow1y;
	double w2z    = ow2z - ow1z;
	double u3x    = ou3x - ou2x;
	double u3y    = ou3y - ou2y;
	double u3z    = ou3z - ou2z;
	double u2x    = ou2x - ou1x;
	double u2y    = ou2y - ou1y;
	double u2z    = ou2z - ou1z;
	double nvx1   = v2y * v3z;
	double nvx2   = v2z * v3y;
	double nvx    = nvx1 - nvx2;
	double nvy1   = v3x * v2z;
	double nvy2   = v3z * v2x;
	double nvy    = nvy1 - nvy2;
	double nvz1   = v2x * v3y;
	double nvz2   = v2y * v3x;
	double nvz    = nvz1 - nvz2;
	double nwx1   = w2y * w3z;
	double nwx2   = w2z * w3y;
	double nwx    = nwx1 - nwx2;
	double nwy1   = w3x * w2z;
	double nwy2   = w3z * w2x;
	double nwy    = nwy1 - nwy2;
	double nwz1   = w2x * w3y;
	double nwz2   = w2y * w3x;
	double nwz    = nwz1 - nwz2;
	double nux1   = u2y * u3z;
	double nux2   = u2z * u3y;
	double nux    = nux1 - nux2;
	double nuy1   = u3x * u2z;
	double nuy2   = u3z * u2x;
	double nuy    = nuy1 - nuy2;
	double nuz1   = u2x * u3y;
	double nuz2   = u2y * u3x;
	double nuz    = nuz1 - nuz2;
	double nwyuz1 = nwy * nuz;
	double nwyuz2 = nwz * nuy;
	double nwyuz  = nwyuz1 - nwyuz2;
	double nwxuz1 = nwx * nuz;
	double nwxuz2 = nwz * nux;
	double nwxuz  = nwxuz1 - nwxuz2;
	double nwxuy1 = nwx * nuy;
	double nwxuy2 = nwy * nux;
	double nwxuy  = nwxuy1 - nwxuy2;
	double nvyuz1 = nvy * nuz;
	double nvyuz2 = nvz * nuy;
	double nvyuz  = nvyuz1 - nvyuz2;
	double nvxuz1 = nvx * nuz;
	double nvxuz2 = nvz * nux;
	double nvxuz  = nvxuz1 - nvxuz2;
	double nvxuy1 = nvx * nuy;
	double nvxuy2 = nvy * nux;
	double nvxuy  = nvxuy1 - nvxuy2;
	double nvywz1 = nvy * nwz;
	double nvywz2 = nvz * nwy;
	double nvywz  = nvywz1 - nvywz2;
	double nvxwz1 = nvx * nwz;
	double nvxwz2 = nvz * nwx;
	double nvxwz  = nvxwz1 - nvxwz2;
	double nvxwy1 = nvx * nwy;
	double nvxwy2 = nvy * nwx;
	double nvxwy  = nvxwy1 - nvxwy2;
	double p1a    = nvx * ov1x;
	double p1b    = nvy * ov1y;
	double p1c    = nvz * ov1z;
	double p1ab   = p1a + p1b;
	double p1     = p1ab + p1c;
	double p2a    = nwx * ow1x;
	double p2b    = nwy * ow1y;
	double p2c    = nwz * ow1z;
	double p2ab   = p2a + p2b;
	double p2     = p2ab + p2c;
	double p3a    = nux * ou1x;
	double p3b    = nuy * ou1y;
	double p3c    = nuz * ou1z;
	double p3ab   = p3a + p3b;
	double p3     = p3ab + p3c;
	double lxa    = p1 * nwyuz;
	double lxb    = p3 * nvywz;
	double lxc    = p2 * nvyuz;
	double lxab   = lxa + lxb;
	lambda_x      = lxab - lxc;
	double lya    = p2 * nvxuz;
	double lyb    = p3 * nvxwz;
	double lyc    = p1 * nwxuz;
	double lybc   = lyc + lyb;
	lambda_y      = lya - lybc;
	double lza    = p3 * nvxwy;
	double lzb    = p1 * nwxuy;
	double lzc    = p2 * nvxuy;
	double lzab   = lza + lzb;
	lambda_z      = lzab - lzc;
	double da     = nvx * nwyuz;
	double db     = nvz * nwxuy;
	double dc     = nvy * nwxuz;
	double dab    = da + db;
	lambda_d      = dab - dc;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(ov1x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ov1y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ov1z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ow1x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ow1y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ow1z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ou1x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ou1y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ou1z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(w3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(w3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(w3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(w2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(w2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(w2z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(u3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(u3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(u3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(u2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(u2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(u2z)) > max_var)
		max_var = _tmp_fabs;
	double lambda_d_eps = max_var;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= max_var;
	lambda_d_eps *= max_var;
	lambda_d_eps *= 8.704148513061234e-14;

	return ((lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps));
}

template <typename IT>
bool lambda2d_SSI_interval(IT ea1x, IT ea1y, IT ea2x, IT ea2y, IT eb1x, IT eb1y,
                           IT eb2x, IT eb2y, IT &lambda_x, IT &lambda_y,
                           IT &lambda_det)
{
	typename IT::Protector P;

	IT t1a(ea1x * ea2y);
	IT t1b(ea2x * ea1y);
	IT t1(t1a - t1b);
	IT tx2(eb1x - eb2x);
	IT t3a(eb1x * eb2y);
	IT t3b(eb2x * eb1y);
	IT t3(t3a - t3b);
	IT tx4(ea1x - ea2x);
	IT ty2(eb1y - eb2y);
	IT ty4(ea1y - ea2y);
	IT lxa(t1 * tx2);
	IT lxb(t3 * tx4);
	lambda_x = lxa - lxb;
	IT lya(t1 * ty2);
	IT lyb(t3 * ty4);
	lambda_y = lya - lyb;
	IT deta(tx4 * ty2);
	IT detb(tx2 * ty4);
	lambda_det = deta - detb;

	return (lambda_det.is_sign_reliable());
}

template <typename ET>
void lambda2d_SSI_exact(const ET &ea1x, const ET &ea1y, const ET &ea2x,
                        const ET &ea2y, const ET &eb1x, const ET &eb1y,
                        const ET &eb2x, const ET &eb2y, ET &lambda_x,
                        ET &lambda_y, ET &lambda_det)
{
	ET t1a = ea1x * ea2y;
	ET t1b = ea2x * ea1y;
	ET t1  = t1a - t1b;
	ET tx2 = eb1x - eb2x;
	ET t3a = eb1x * eb2y;
	ET t3b = eb2x * eb1y;
	ET t3  = t3a - t3b;
	ET tx4 = ea1x - ea2x;
	ET ty2 = eb1y - eb2y;
	ET ty4 = ea1y - ea2y;
	ET lxa = t1 * tx2;
	ET lxb = t3 * tx4;

	lambda_x   = lxa - lxb;
	ET lya     = t1 * ty2;
	ET lyb     = t3 * ty4;
	lambda_y   = lya - lyb;
	ET deta    = tx4 * ty2;
	ET detb    = tx2 * ty4;
	lambda_det = deta - detb;
}

template <typename IT>
bool lambda3d_SSI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xp,
                           IT yp, IT xq, IT yq, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z, IT &lambda_d)
{
	typename IT::Protector P;

	IT xba   = xb - xa;
	IT yba   = yb - ya;
	IT zba   = zb - za;
	IT xap   = xa - xp;
	IT yap   = ya - yp;
	IT yqp   = yq - yp;
	IT xqp   = xq - xp;
	IT c1    = xap * yqp;
	IT c2    = xqp * yap;
	IT c3    = xba * yqp;
	IT c4    = xqp * yba;
	IT c12   = c1 - c2;
	lambda_d = c4 - c3;
	IT xf    = xba * c12;
	IT xs    = xa * lambda_d;
	lambda_x = xf + xs;
	IT yf    = yba * c12;
	IT ys    = ya * lambda_d;
	lambda_y = yf + ys;
	IT zf    = zba * c12;
	IT zs    = za * lambda_d;
	lambda_z = zf + zs;
	return lambda_d.is_sign_reliable();
}

template <typename ET>
void lambda3d_SSI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xp, ET yp,
                        ET xq, ET yq, ET &lambda_x, ET &lambda_y, ET &lambda_z,
                        ET &lambda_d)
{
	ET xba   = xb - xa;
	ET yba   = yb - ya;
	ET zba   = zb - za;
	ET xap   = xa - xp;
	ET yap   = ya - yp;
	ET yqp   = yq - yp;
	ET xqp   = xq - xp;
	ET c1    = xap * yqp;
	ET c2    = xqp * yap;
	ET c3    = xba * yqp;
	ET c4    = xqp * yba;
	ET c12   = c1 - c2;
	lambda_d = c4 - c3;
	ET xf    = xba * c12;
	ET xs    = xa * lambda_d;
	lambda_x = xf + xs;
	ET yf    = yba * c12;
	ET ys    = ya * lambda_d;
	lambda_y = yf + ys;
	ET zf    = zba * c12;
	ET zs    = za * lambda_d;
	lambda_z = zf + zs;
}

template <typename IT>
bool lambda3d_LPI_interval(IT px, IT py, IT pz, IT qx, IT qy, IT qz, IT rx,
                           IT ry, IT rz, IT sx, IT sy, IT sz, IT tx, IT ty,
                           IT tz, IT &lambda_d, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z)
{
	typename IT::Protector P;

	IT a11(px - qx);
	IT a12(py - qy);
	IT a13(pz - qz);
	IT a21(sx - rx);
	IT a22(sy - ry);
	IT a23(sz - rz);
	IT a31(tx - rx);
	IT a32(ty - ry);
	IT a33(tz - rz);
	IT tv1(a22 * a33);
	IT tv2(a23 * a32);
	IT a2233(tv1 - tv2);
	IT tv3(a21 * a33);
	IT tv4(a23 * a31);
	IT a2133(tv3 - tv4);
	IT tv5(a21 * a32);
	IT tv6(a22 * a31);
	IT a2132(tv5 - tv6);
	IT tv7(a11 * a2233);
	IT tv8(a12 * a2133);
	IT tv9(a13 * a2132);
	IT tt1(tv7 - tv8);
	lambda_d = tt1 + tv9;
	IT px_rx(px - rx);
	IT py_ry(py - ry);
	IT pz_rz(pz - rz);
	IT tt2(py_ry * a2133);
	IT tt3(px_rx * a2233);
	IT tt4(pz_rz * a2132);
	IT tt5(tt3 + tt4);
	IT n(tt5 - tt2);
	IT ax(a11 * n);
	IT ay(a12 * n);
	IT az(a13 * n);
	IT dpx(lambda_d * px);
	IT dpy(lambda_d * py);
	IT dpz(lambda_d * pz);
	lambda_x = dpx - ax;
	lambda_y = dpy - ay;
	lambda_z = dpz - az;

	return (lambda_d.is_sign_reliable());
}

template <typename ET>
void lambda3d_LPI_exact(const ET &px, const ET &py, const ET &pz, const ET &qx,
                        const ET &qy, const ET &qz, const ET &rx, const ET &ry,
                        const ET &rz, const ET &sx, const ET &sy, const ET &sz,
                        const ET &tx, const ET &ty, const ET &tz, ET &lambda_d,
                        ET &lambda_x, ET &lambda_y, ET &lambda_z)
{
	ET a11   = px - qx;
	ET a12   = py - qy;
	ET a13   = pz - qz;
	ET a21   = sx - rx;
	ET a22   = sy - ry;
	ET a23   = sz - rz;
	ET a31   = tx - rx;
	ET a32   = ty - ry;
	ET a33   = tz - rz;
	ET tv1   = a22 * a33;
	ET tv2   = a23 * a32;
	ET a2233 = tv1 - tv2;
	ET tv3   = a21 * a33;
	ET tv4   = a23 * a31;
	ET a2133 = tv3 - tv4;
	ET tv5   = a21 * a32;
	ET tv6   = a22 * a31;
	ET a2132 = tv5 - tv6;
	ET tv7   = a11 * a2233;
	ET tv8   = a12 * a2133;
	ET tv9   = a13 * a2132;
	ET tt1   = tv7 - tv8;
	lambda_d = tt1 + tv9;
	ET px_rx = px - rx;
	ET py_ry = py - ry;
	ET pz_rz = pz - rz;
	ET tt2   = py_ry * a2133;
	ET tt3   = px_rx * a2233;
	ET tt4   = pz_rz * a2132;
	ET tt5   = tt3 + tt4;
	ET n     = tt5 - tt2;
	ET ax    = a11 * n;
	ET ay    = a12 * n;
	ET az    = a13 * n;
	ET dpx   = lambda_d * px;
	ET dpy   = lambda_d * py;
	ET dpz   = lambda_d * pz;
	lambda_x = dpx - ax;
	lambda_y = dpy - ay;
	lambda_z = dpz - az;
}

template <typename IT>
bool lambda3d_TPI_interval(IT ov1x, IT ov1y, IT ov1z, IT ov2x, IT ov2y, IT ov2z,
                           IT ov3x, IT ov3y, IT ov3z, IT ow1x, IT ow1y, IT ow1z,
                           IT ow2x, IT ow2y, IT ow2z, IT ow3x, IT ow3y, IT ow3z,
                           IT ou1x, IT ou1y, IT ou1z, IT ou2x, IT ou2y, IT ou2z,
                           IT ou3x, IT ou3y, IT ou3z, IT &lambda_x,
                           IT &lambda_y, IT &lambda_z, IT &lambda_d)
{
	typename IT::Protector P;

	IT v3x(ov3x - ov2x);
	IT v3y(ov3y - ov2y);
	IT v3z(ov3z - ov2z);
	IT v2x(ov2x - ov1x);
	IT v2y(ov2y - ov1y);
	IT v2z(ov2z - ov1z);
	IT w3x(ow3x - ow2x);
	IT w3y(ow3y - ow2y);
	IT w3z(ow3z - ow2z);
	IT w2x(ow2x - ow1x);
	IT w2y(ow2y - ow1y);
	IT w2z(ow2z - ow1z);
	IT u3x(ou3x - ou2x);
	IT u3y(ou3y - ou2y);
	IT u3z(ou3z - ou2z);
	IT u2x(ou2x - ou1x);
	IT u2y(ou2y - ou1y);
	IT u2z(ou2z - ou1z);
	IT nvx1(v2y * v3z);
	IT nvx2(v2z * v3y);
	IT nvx(nvx1 - nvx2);
	IT nvy1(v3x * v2z);
	IT nvy2(v3z * v2x);
	IT nvy(nvy1 - nvy2);
	IT nvz1(v2x * v3y);
	IT nvz2(v2y * v3x);
	IT nvz(nvz1 - nvz2);
	IT nwx1(w2y * w3z);
	IT nwx2(w2z * w3y);
	IT nwx(nwx1 - nwx2);
	IT nwy1(w3x * w2z);
	IT nwy2(w3z * w2x);
	IT nwy(nwy1 - nwy2);
	IT nwz1(w2x * w3y);
	IT nwz2(w2y * w3x);
	IT nwz(nwz1 - nwz2);
	IT nux1(u2y * u3z);
	IT nux2(u2z * u3y);
	IT nux(nux1 - nux2);
	IT nuy1(u3x * u2z);
	IT nuy2(u3z * u2x);
	IT nuy(nuy1 - nuy2);
	IT nuz1(u2x * u3y);
	IT nuz2(u2y * u3x);
	IT nuz(nuz1 - nuz2);
	IT nwyuz1(nwy * nuz);
	IT nwyuz2(nwz * nuy);
	IT nwyuz(nwyuz1 - nwyuz2);
	IT nwxuz1(nwx * nuz);
	IT nwxuz2(nwz * nux);
	IT nwxuz(nwxuz1 - nwxuz2);
	IT nwxuy1(nwx * nuy);
	IT nwxuy2(nwy * nux);
	IT nwxuy(nwxuy1 - nwxuy2);
	IT nvyuz1(nvy * nuz);
	IT nvyuz2(nvz * nuy);
	IT nvyuz(nvyuz1 - nvyuz2);
	IT nvxuz1(nvx * nuz);
	IT nvxuz2(nvz * nux);
	IT nvxuz(nvxuz1 - nvxuz2);
	IT nvxuy1(nvx * nuy);
	IT nvxuy2(nvy * nux);
	IT nvxuy(nvxuy1 - nvxuy2);
	IT nvywz1(nvy * nwz);
	IT nvywz2(nvz * nwy);
	IT nvywz(nvywz1 - nvywz2);
	IT nvxwz1(nvx * nwz);
	IT nvxwz2(nvz * nwx);
	IT nvxwz(nvxwz1 - nvxwz2);
	IT nvxwy1(nvx * nwy);
	IT nvxwy2(nvy * nwx);
	IT nvxwy(nvxwy1 - nvxwy2);
	IT p1a(nvx * ov1x);
	IT p1b(nvy * ov1y);
	IT p1c(nvz * ov1z);
	IT p1ab(p1a + p1b);
	IT p1(p1ab + p1c);
	IT p2a(nwx * ow1x);
	IT p2b(nwy * ow1y);
	IT p2c(nwz * ow1z);
	IT p2ab(p2a + p2b);
	IT p2(p2ab + p2c);
	IT p3a(nux * ou1x);
	IT p3b(nuy * ou1y);
	IT p3c(nuz * ou1z);
	IT p3ab(p3a + p3b);
	IT p3(p3ab + p3c);
	IT lxa(p1 * nwyuz);
	IT lxb(p3 * nvywz);
	IT lxc(p2 * nvyuz);
	IT lxab(lxa + lxb);
	lambda_x = lxab - lxc;
	IT lya(p2 * nvxuz);
	IT lyb(p3 * nvxwz);
	IT lyc(p1 * nwxuz);
	IT lybc(lyc + lyb);
	lambda_y = lya - lybc;
	IT lza(p3 * nvxwy);
	IT lzb(p1 * nwxuy);
	IT lzc(p2 * nvxuy);
	IT lzab(lza + lzb);
	lambda_z = lzab - lzc;
	IT da(nvx * nwyuz);
	IT db(nvz * nwxuy);
	IT dc(nvy * nwxuz);
	IT dab(da + db);
	lambda_d = dab - dc;

	return (lambda_d.is_sign_reliable());
}

template <typename ET>
void lambda3d_TPI_exact(const ET &ov1x, const ET &ov1y, const ET &ov1z,
                        const ET &ov2x, const ET &ov2y, const ET &ov2z,
                        const ET &ov3x, const ET &ov3y, const ET &ov3z,
                        const ET &ow1x, const ET &ow1y, const ET &ow1z,
                        const ET &ow2x, const ET &ow2y, const ET &ow2z,
                        const ET &ow3x, const ET &ow3y, const ET &ow3z,
                        const ET &ou1x, const ET &ou1y, const ET &ou1z,
                        const ET &ou2x, const ET &ou2y, const ET &ou2z,
                        const ET &ou3x, const ET &ou3y, const ET &ou3z,
                        ET &lambda_x, ET &lambda_y, ET &lambda_z, ET &lambda_d)
{
	ET v3x    = ov3x - ov2x;
	ET v3y    = ov3y - ov2y;
	ET v3z    = ov3z - ov2z;
	ET v2x    = ov2x - ov1x;
	ET v2y    = ov2y - ov1y;
	ET v2z    = ov2z - ov1z;
	ET w3x    = ow3x - ow2x;
	ET w3y    = ow3y - ow2y;
	ET w3z    = ow3z - ow2z;
	ET w2x    = ow2x - ow1x;
	ET w2y    = ow2y - ow1y;
	ET w2z    = ow2z - ow1z;
	ET u3x    = ou3x - ou2x;
	ET u3y    = ou3y - ou2y;
	ET u3z    = ou3z - ou2z;
	ET u2x    = ou2x - ou1x;
	ET u2y    = ou2y - ou1y;
	ET u2z    = ou2z - ou1z;
	ET nvx1   = v2y * v3z;
	ET nvx2   = v2z * v3y;
	ET nvx    = nvx1 - nvx2;
	ET nvy1   = v3x * v2z;
	ET nvy2   = v3z * v2x;
	ET nvy    = nvy1 - nvy2;
	ET nvz1   = v2x * v3y;
	ET nvz2   = v2y * v3x;
	ET nvz    = nvz1 - nvz2;
	ET nwx1   = w2y * w3z;
	ET nwx2   = w2z * w3y;
	ET nwx    = nwx1 - nwx2;
	ET nwy1   = w3x * w2z;
	ET nwy2   = w3z * w2x;
	ET nwy    = nwy1 - nwy2;
	ET nwz1   = w2x * w3y;
	ET nwz2   = w2y * w3x;
	ET nwz    = nwz1 - nwz2;
	ET nux1   = u2y * u3z;
	ET nux2   = u2z * u3y;
	ET nux    = nux1 - nux2;
	ET nuy1   = u3x * u2z;
	ET nuy2   = u3z * u2x;
	ET nuy    = nuy1 - nuy2;
	ET nuz1   = u2x * u3y;
	ET nuz2   = u2y * u3x;
	ET nuz    = nuz1 - nuz2;
	ET nwyuz1 = nwy * nuz;
	ET nwyuz2 = nwz * nuy;
	ET nwyuz  = nwyuz1 - nwyuz2;
	ET nwxuz1 = nwx * nuz;
	ET nwxuz2 = nwz * nux;
	ET nwxuz  = nwxuz1 - nwxuz2;
	ET nwxuy1 = nwx * nuy;
	ET nwxuy2 = nwy * nux;
	ET nwxuy  = nwxuy1 - nwxuy2;
	ET nvyuz1 = nvy * nuz;
	ET nvyuz2 = nvz * nuy;
	ET nvyuz  = nvyuz1 - nvyuz2;
	ET nvxuz1 = nvx * nuz;
	ET nvxuz2 = nvz * nux;
	ET nvxuz  = nvxuz1 - nvxuz2;
	ET nvxuy1 = nvx * nuy;
	ET nvxuy2 = nvy * nux;
	ET nvxuy  = nvxuy1 - nvxuy2;
	ET nvywz1 = nvy * nwz;
	ET nvywz2 = nvz * nwy;
	ET nvywz  = nvywz1 - nvywz2;
	ET nvxwz1 = nvx * nwz;
	ET nvxwz2 = nvz * nwx;
	ET nvxwz  = nvxwz1 - nvxwz2;
	ET nvxwy1 = nvx * nwy;
	ET nvxwy2 = nvy * nwx;
	ET nvxwy  = nvxwy1 - nvxwy2;
	ET p1a    = nvx * ov1x;
	ET p1b    = nvy * ov1y;
	ET p1c    = nvz * ov1z;
	ET p1ab   = p1a + p1b;
	ET p1     = p1ab + p1c;
	ET p2a    = nwx * ow1x;
	ET p2b    = nwy * ow1y;
	ET p2c    = nwz * ow1z;
	ET p2ab   = p2a + p2b;
	ET p2     = p2ab + p2c;
	ET p3a    = nux * ou1x;
	ET p3b    = nuy * ou1y;
	ET p3c    = nuz * ou1z;
	ET p3ab   = p3a + p3b;
	ET p3     = p3ab + p3c;
	ET lxa    = p1 * nwyuz;
	ET lxb    = p3 * nvywz;
	ET lxc    = p2 * nvyuz;
	ET lxab   = lxa + lxb;
	lambda_x  = lxab - lxc;
	ET lya    = p2 * nvxuz;
	ET lyb    = p3 * nvxwz;
	ET lyc    = p1 * nwxuz;
	ET lybc   = lyc + lyb;
	lambda_y  = lya - lybc;
	ET lza    = p3 * nvxwy;
	ET lzb    = p1 * nwxuy;
	ET lzc    = p2 * nvxuy;
	ET lzab   = lza + lzb;
	lambda_z  = lzab - lzc;
	ET da     = nvx * nwyuz;
	ET db     = nvz * nwxuy;
	ET dc     = nvy * nwxuz;
	ET dab    = da + db;
	lambda_d  = dab - dc;
}

template <typename IT>
bool lambda3d_LNC_interval(IT px, IT py, IT pz, IT qx, IT qy, IT qz, IT t,
                           IT &lambda_x, IT &lambda_y, IT &lambda_z,
                           IT &lambda_d)
{
	typename IT::Protector P;

	IT vx(px - qx);
	IT vy(py - qy);
	IT vz(pz - qz);
	IT vxt(vx * t);
	IT vyt(vy * t);
	IT vzt(vz * t);
	lambda_x = px - vxt;
	lambda_y = py - vyt;
	lambda_z = pz - vzt;
	lambda_d = 1;

	return true;
}

template <typename ET>
void lambda3d_LNC_exact(ET px, ET py, ET pz, ET qx, ET qy, ET qz, ET t,
                        ET &lambda_x, ET &lambda_y, ET &lambda_z, ET &lambda_d)
{
	ET vx    = px - qx;
	ET vy    = py - qy;
	ET vz    = pz - qz;
	ET vxt   = vx * t;
	ET vyt   = vy * t;
	ET vzt   = vz * t;
	lambda_x = px - vxt;
	lambda_y = py - vyt;
	lambda_z = pz - vzt;
	lambda_d = 1;
}

inline void lambda2d_SSI_expansion(double ea1x, double ea1y, double ea2x,
                                   double ea2y, double eb1x, double eb1y,
                                   double eb2x, double eb2y, double **lambda_x,
                                   int &lambda_x_len, double **lambda_y,
                                   int &lambda_y_len, double **lambda_det,
                                   int &lambda_det_len)
{
	expansionObject o;
	double          t1a[2];
	o.Two_Prod(ea1x, ea2y, t1a);
	double t1b[2];
	o.Two_Prod(ea2x, ea1y, t1b);
	double t1[4];
	o.Two_Two_Diff(t1a, t1b, t1);
	double tx2[2];
	o.two_Diff(eb1x, eb2x, tx2);
	double t3a[2];
	o.Two_Prod(eb1x, eb2y, t3a);
	double t3b[2];
	o.Two_Prod(eb2x, eb1y, t3b);
	double t3[4];
	o.Two_Two_Diff(t3a, t3b, t3);
	double tx4[2];
	o.two_Diff(ea1x, ea2x, tx4);
	double ty2[2];
	o.two_Diff(eb1y, eb2y, ty2);
	double ty4[2];
	o.two_Diff(ea1y, ea2y, ty4);
	double lxa[16];
	int    lxa_len = o.Gen_Product(4, t1, 2, tx2, lxa);
	double lxb[16];
	int    lxb_len = o.Gen_Product(4, t3, 2, tx4, lxb);
	lambda_x_len = o.Gen_Diff_With_PreAlloc(lxa_len, lxa, lxb_len, lxb, lambda_x,
	                                        lambda_x_len);
	double lya[16];
	int    lya_len = o.Gen_Product(4, t1, 2, ty2, lya);
	double lyb[16];
	int    lyb_len = o.Gen_Product(4, t3, 2, ty4, lyb);
	lambda_y_len = o.Gen_Diff_With_PreAlloc(lya_len, lya, lyb_len, lyb, lambda_y,
	                                        lambda_y_len);
	double deta[8];
	int    deta_len = o.Gen_Product(2, tx4, 2, ty2, deta);
	double detb[8];
	int    detb_len = o.Gen_Product(2, tx2, 2, ty4, detb);
	lambda_det_len  = o.Gen_Diff_With_PreAlloc(deta_len, deta, detb_len, detb,
	                                           lambda_det, lambda_det_len);
}

inline void lambda3d_LNC_expansion(double px, double py, double pz, double qx,
                                   double qy, double qz, double t,
                                   double **lambda_x, int &lambda_x_len,
                                   double **lambda_y, int &lambda_y_len,
                                   double **lambda_z, int &lambda_z_len,
                                   double **lambda_d, int &lambda_d_len)
{
	expansionObject o;
	double          vx[2];
	o.two_Diff(px, qx, vx);
	double vy[2];
	o.two_Diff(py, qy, vy);
	double vz[2];
	o.two_Diff(pz, qz, vz);
	double vxt[4];
	o.Two_One_Prod(vx, t, vxt);
	double vyt[4];
	o.Two_One_Prod(vy, t, vyt);
	double vzt[4];
	o.Two_One_Prod(vz, t, vzt);
	lambda_x_len =
	  o.Gen_Diff_With_PreAlloc(1, &px, 4, vxt, lambda_x, lambda_x_len);
	lambda_y_len =
	  o.Gen_Diff_With_PreAlloc(1, &py, 4, vyt, lambda_y, lambda_y_len);
	lambda_z_len =
	  o.Gen_Diff_With_PreAlloc(1, &pz, 4, vzt, lambda_z, lambda_z_len);
	(*lambda_d)[0] = 1;
	lambda_d_len   = 1;
}

inline void lambda3d_SSI_expansion(double xa, double ya, double za, double xb,
                                   double yb, double zb, double xp, double yp,
                                   double xq, double yq, double **lambda_x,
                                   int &lambda_x_len, double **lambda_y,
                                   int &lambda_y_len, double **lambda_z,
                                   int &lambda_z_len, double **lambda_d,
                                   int &lambda_d_len)
{
	expansionObject o;
	double          xba[2];
	o.two_Diff(xb, xa, xba);
	double yba[2];
	o.two_Diff(yb, ya, yba);
	double zba[2];
	o.two_Diff(zb, za, zba);
	double xap[2];
	o.two_Diff(xa, xp, xap);
	double yap[2];
	o.two_Diff(ya, yp, yap);
	double yqp[2];
	o.two_Diff(yq, yp, yqp);
	double xqp[2];
	o.two_Diff(xq, xp, xqp);
	double c1[8];
	int    c1_len = o.Gen_Product(2, xap, 2, yqp, c1);
	double c2[8];
	int    c2_len = o.Gen_Product(2, xqp, 2, yap, c2);
	double c3[8];
	int    c3_len = o.Gen_Product(2, xba, 2, yqp, c3);
	double c4[8];
	int    c4_len = o.Gen_Product(2, xqp, 2, yba, c4);
	double c12[16];
	int    c12_len = o.Gen_Diff(c1_len, c1, c2_len, c2, c12);
	lambda_d_len =
	  o.Gen_Diff_With_PreAlloc(c4_len, c4, c3_len, c3, lambda_d, lambda_d_len);
	double xf[64];
	int    xf_len = o.Gen_Product(2, xba, c12_len, c12, xf);
	double xs[32];
	int    xs_len = o.Gen_Scale(lambda_d_len, *lambda_d, xa, xs);
	lambda_x_len =
	  o.Gen_Sum_With_PreAlloc(xf_len, xf, xs_len, xs, lambda_x, lambda_x_len);
	double yf[64];
	int    yf_len = o.Gen_Product(2, yba, c12_len, c12, yf);
	double ys[32];
	int    ys_len = o.Gen_Scale(lambda_d_len, *lambda_d, ya, ys);
	lambda_y_len =
	  o.Gen_Sum_With_PreAlloc(yf_len, yf, ys_len, ys, lambda_y, lambda_y_len);
	double zf[64];
	int    zf_len = o.Gen_Product(2, zba, c12_len, c12, zf);
	double zs[32];
	int    zs_len = o.Gen_Scale(lambda_d_len, *lambda_d, za, zs);
	lambda_z_len =
	  o.Gen_Sum_With_PreAlloc(zf_len, zf, zs_len, zs, lambda_z, lambda_z_len);
}

inline void lambda3d_LPI_expansion(double px, double py, double pz, double qx,
                                   double qy, double qz, double rx, double ry,
                                   double rz, double sx, double sy, double sz,
                                   double tx, double ty, double tz,
                                   double **lambda_d, int &lambda_d_len,
                                   double **lambda_x, int &lambda_x_len,
                                   double **lambda_y, int &lambda_y_len,
                                   double **lambda_z, int &lambda_z_len)
{
	expansionObject o;
	double          a11[2];
	o.two_Diff(px, qx, a11);
	double a12[2];
	o.two_Diff(py, qy, a12);
	double a13[2];
	o.two_Diff(pz, qz, a13);
	double a21[2];
	o.two_Diff(sx, rx, a21);
	double a22[2];
	o.two_Diff(sy, ry, a22);
	double a23[2];
	o.two_Diff(sz, rz, a23);
	double a31[2];
	o.two_Diff(tx, rx, a31);
	double a32[2];
	o.two_Diff(ty, ry, a32);
	double a33[2];
	o.two_Diff(tz, rz, a33);
	double tv1[8];
	int    tv1_len = o.Gen_Product(2, a22, 2, a33, tv1);
	double tv2[8];
	int    tv2_len = o.Gen_Product(2, a23, 2, a32, tv2);
	double a2233[16];
	int    a2233_len = o.Gen_Diff(tv1_len, tv1, tv2_len, tv2, a2233);
	double tv3[8];
	int    tv3_len = o.Gen_Product(2, a21, 2, a33, tv3);
	double tv4[8];
	int    tv4_len = o.Gen_Product(2, a23, 2, a31, tv4);
	double a2133[16];
	int    a2133_len = o.Gen_Diff(tv3_len, tv3, tv4_len, tv4, a2133);
	double tv5[8];
	int    tv5_len = o.Gen_Product(2, a21, 2, a32, tv5);
	double tv6[8];
	int    tv6_len = o.Gen_Product(2, a22, 2, a31, tv6);
	double a2132[16];
	int    a2132_len = o.Gen_Diff(tv5_len, tv5, tv6_len, tv6, a2132);
	double tv7[64];
	int    tv7_len = o.Gen_Product(2, a11, a2233_len, a2233, tv7);
	double tv8[64];
	int    tv8_len = o.Gen_Product(2, a12, a2133_len, a2133, tv8);
	double tv9[64];
	int    tv9_len = o.Gen_Product(2, a13, a2132_len, a2132, tv9);
	double tt1[128];
	int    tt1_len = o.Gen_Diff(tv7_len, tv7, tv8_len, tv8, tt1);
	lambda_d_len =
	  o.Gen_Sum_With_PreAlloc(tt1_len, tt1, tv9_len, tv9, lambda_d, lambda_d_len);
	double px_rx[2];
	o.two_Diff(px, rx, px_rx);
	double py_ry[2];
	o.two_Diff(py, ry, py_ry);
	double pz_rz[2];
	o.two_Diff(pz, rz, pz_rz);
	double tt2[64];
	int    tt2_len = o.Gen_Product(2, py_ry, a2133_len, a2133, tt2);
	double tt3[64];
	int    tt3_len = o.Gen_Product(2, px_rx, a2233_len, a2233, tt3);
	double tt4[64];
	int    tt4_len = o.Gen_Product(2, pz_rz, a2132_len, a2132, tt4);
	double tt5[128];
	int    tt5_len = o.Gen_Sum(tt3_len, tt3, tt4_len, tt4, tt5);
	double n_p[128], *n = n_p;
	int    n_len = o.Gen_Diff_With_PreAlloc(tt5_len, tt5, tt2_len, tt2, &n, 128);
	double ax_p[128], *ax = ax_p;
	int    ax_len = o.Gen_Product_With_PreAlloc(2, a11, n_len, n, &ax, 128);
	double ay_p[128], *ay = ay_p;
	int    ay_len = o.Gen_Product_With_PreAlloc(2, a12, n_len, n, &ay, 128);
	double az_p[128], *az = az_p;
	int    az_len = o.Gen_Product_With_PreAlloc(2, a13, n_len, n, &az, 128);
	double dpx_p[128], *dpx = dpx_p;
	int    dpx_len =
	  o.Gen_Scale_With_PreAlloc(lambda_d_len, *lambda_d, px, &dpx, 128);
	double dpy_p[128], *dpy = dpy_p;
	int    dpy_len =
	  o.Gen_Scale_With_PreAlloc(lambda_d_len, *lambda_d, py, &dpy, 128);
	double dpz_p[128], *dpz = dpz_p;
	int    dpz_len =
	  o.Gen_Scale_With_PreAlloc(lambda_d_len, *lambda_d, pz, &dpz, 128);
	lambda_x_len =
	  o.Gen_Diff_With_PreAlloc(dpx_len, dpx, ax_len, ax, lambda_x, lambda_x_len);
	lambda_y_len =
	  o.Gen_Diff_With_PreAlloc(dpy_len, dpy, ay_len, ay, lambda_y, lambda_y_len);
	lambda_z_len =
	  o.Gen_Diff_With_PreAlloc(dpz_len, dpz, az_len, az, lambda_z, lambda_z_len);

	if (dpz_p != dpz)
		FreeDoubles(dpz);
	if (dpy_p != dpy)
		FreeDoubles(dpy);
	if (dpx_p != dpx)
		FreeDoubles(dpx);
	if (az_p != az)
		FreeDoubles(az);
	if (ay_p != ay)
		FreeDoubles(ay);
	if (ax_p != ax)
		FreeDoubles(ax);
	if (n_p != n)
		FreeDoubles(n);
}

inline void lambda3d_TPI_expansion(
  double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z,
  double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z,
  double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z,
  double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z,
  double ou3x, double ou3y, double ou3z, double **lambda_x, int &lambda_x_len,
  double **lambda_y, int &lambda_y_len, double **lambda_z, int &lambda_z_len,
  double **lambda_d, int &lambda_d_len)
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
	double w3x[2];
	o.two_Diff(ow3x, ow2x, w3x);
	double w3y[2];
	o.two_Diff(ow3y, ow2y, w3y);
	double w3z[2];
	o.two_Diff(ow3z, ow2z, w3z);
	double w2x[2];
	o.two_Diff(ow2x, ow1x, w2x);
	double w2y[2];
	o.two_Diff(ow2y, ow1y, w2y);
	double w2z[2];
	o.two_Diff(ow2z, ow1z, w2z);
	double u3x[2];
	o.two_Diff(ou3x, ou2x, u3x);
	double u3y[2];
	o.two_Diff(ou3y, ou2y, u3y);
	double u3z[2];
	o.two_Diff(ou3z, ou2z, u3z);
	double u2x[2];
	o.two_Diff(ou2x, ou1x, u2x);
	double u2y[2];
	o.two_Diff(ou2y, ou1y, u2y);
	double u2z[2];
	o.two_Diff(ou2z, ou1z, u2z);
	double nvx1[8];
	int    nvx1_len = o.Gen_Product(2, v2y, 2, v3z, nvx1);
	double nvx2[8];
	int    nvx2_len = o.Gen_Product(2, v2z, 2, v3y, nvx2);
	double nvx[16];
	int    nvx_len = o.Gen_Diff(nvx1_len, nvx1, nvx2_len, nvx2, nvx);
	double nvy1[8];
	int    nvy1_len = o.Gen_Product(2, v3x, 2, v2z, nvy1);
	double nvy2[8];
	int    nvy2_len = o.Gen_Product(2, v3z, 2, v2x, nvy2);
	double nvy[16];
	int    nvy_len = o.Gen_Diff(nvy1_len, nvy1, nvy2_len, nvy2, nvy);
	double nvz1[8];
	int    nvz1_len = o.Gen_Product(2, v2x, 2, v3y, nvz1);
	double nvz2[8];
	int    nvz2_len = o.Gen_Product(2, v2y, 2, v3x, nvz2);
	double nvz[16];
	int    nvz_len = o.Gen_Diff(nvz1_len, nvz1, nvz2_len, nvz2, nvz);
	double nwx1[8];
	int    nwx1_len = o.Gen_Product(2, w2y, 2, w3z, nwx1);
	double nwx2[8];
	int    nwx2_len = o.Gen_Product(2, w2z, 2, w3y, nwx2);
	double nwx[16];
	int    nwx_len = o.Gen_Diff(nwx1_len, nwx1, nwx2_len, nwx2, nwx);
	double nwy1[8];
	int    nwy1_len = o.Gen_Product(2, w3x, 2, w2z, nwy1);
	double nwy2[8];
	int    nwy2_len = o.Gen_Product(2, w3z, 2, w2x, nwy2);
	double nwy[16];
	int    nwy_len = o.Gen_Diff(nwy1_len, nwy1, nwy2_len, nwy2, nwy);
	double nwz1[8];
	int    nwz1_len = o.Gen_Product(2, w2x, 2, w3y, nwz1);
	double nwz2[8];
	int    nwz2_len = o.Gen_Product(2, w2y, 2, w3x, nwz2);
	double nwz[16];
	int    nwz_len = o.Gen_Diff(nwz1_len, nwz1, nwz2_len, nwz2, nwz);
	double nux1[8];
	int    nux1_len = o.Gen_Product(2, u2y, 2, u3z, nux1);
	double nux2[8];
	int    nux2_len = o.Gen_Product(2, u2z, 2, u3y, nux2);
	double nux[16];
	int    nux_len = o.Gen_Diff(nux1_len, nux1, nux2_len, nux2, nux);
	double nuy1[8];
	int    nuy1_len = o.Gen_Product(2, u3x, 2, u2z, nuy1);
	double nuy2[8];
	int    nuy2_len = o.Gen_Product(2, u3z, 2, u2x, nuy2);
	double nuy[16];
	int    nuy_len = o.Gen_Diff(nuy1_len, nuy1, nuy2_len, nuy2, nuy);
	double nuz1[8];
	int    nuz1_len = o.Gen_Product(2, u2x, 2, u3y, nuz1);
	double nuz2[8];
	int    nuz2_len = o.Gen_Product(2, u2y, 2, u3x, nuz2);
	double nuz[16];
	int    nuz_len = o.Gen_Diff(nuz1_len, nuz1, nuz2_len, nuz2, nuz);
	double nwyuz1_p[16], *nwyuz1 = nwyuz1_p;
	int    nwyuz1_len =
	  o.Gen_Product_With_PreAlloc(nwy_len, nwy, nuz_len, nuz, &nwyuz1, 16);
	double nwyuz2_p[16], *nwyuz2 = nwyuz2_p;
	int    nwyuz2_len =
	  o.Gen_Product_With_PreAlloc(nwz_len, nwz, nuy_len, nuy, &nwyuz2, 16);
	double nwyuz_p[16], *nwyuz = nwyuz_p;
	int    nwyuz_len = o.Gen_Diff_With_PreAlloc(nwyuz1_len, nwyuz1, nwyuz2_len,
	                                            nwyuz2, &nwyuz, 16);
	double nwxuz1_p[16], *nwxuz1 = nwxuz1_p;
	int    nwxuz1_len =
	  o.Gen_Product_With_PreAlloc(nwx_len, nwx, nuz_len, nuz, &nwxuz1, 16);
	double nwxuz2_p[16], *nwxuz2 = nwxuz2_p;
	int    nwxuz2_len =
	  o.Gen_Product_With_PreAlloc(nwz_len, nwz, nux_len, nux, &nwxuz2, 16);
	double nwxuz_p[16], *nwxuz = nwxuz_p;
	int    nwxuz_len = o.Gen_Diff_With_PreAlloc(nwxuz1_len, nwxuz1, nwxuz2_len,
	                                            nwxuz2, &nwxuz, 16);
	double nwxuy1_p[16], *nwxuy1 = nwxuy1_p;
	int    nwxuy1_len =
	  o.Gen_Product_With_PreAlloc(nwx_len, nwx, nuy_len, nuy, &nwxuy1, 16);
	double nwxuy2_p[16], *nwxuy2 = nwxuy2_p;
	int    nwxuy2_len =
	  o.Gen_Product_With_PreAlloc(nwy_len, nwy, nux_len, nux, &nwxuy2, 16);
	double nwxuy_p[16], *nwxuy = nwxuy_p;
	int    nwxuy_len = o.Gen_Diff_With_PreAlloc(nwxuy1_len, nwxuy1, nwxuy2_len,
	                                            nwxuy2, &nwxuy, 16);
	double nvyuz1_p[16], *nvyuz1 = nvyuz1_p;
	int    nvyuz1_len =
	  o.Gen_Product_With_PreAlloc(nvy_len, nvy, nuz_len, nuz, &nvyuz1, 16);
	double nvyuz2_p[16], *nvyuz2 = nvyuz2_p;
	int    nvyuz2_len =
	  o.Gen_Product_With_PreAlloc(nvz_len, nvz, nuy_len, nuy, &nvyuz2, 16);
	double nvyuz_p[16], *nvyuz = nvyuz_p;
	int    nvyuz_len = o.Gen_Diff_With_PreAlloc(nvyuz1_len, nvyuz1, nvyuz2_len,
	                                            nvyuz2, &nvyuz, 16);
	double nvxuz1_p[16], *nvxuz1 = nvxuz1_p;
	int    nvxuz1_len =
	  o.Gen_Product_With_PreAlloc(nvx_len, nvx, nuz_len, nuz, &nvxuz1, 16);
	double nvxuz2_p[16], *nvxuz2 = nvxuz2_p;
	int    nvxuz2_len =
	  o.Gen_Product_With_PreAlloc(nvz_len, nvz, nux_len, nux, &nvxuz2, 16);
	double nvxuz_p[16], *nvxuz = nvxuz_p;
	int    nvxuz_len = o.Gen_Diff_With_PreAlloc(nvxuz1_len, nvxuz1, nvxuz2_len,
	                                            nvxuz2, &nvxuz, 16);
	double nvxuy1_p[16], *nvxuy1 = nvxuy1_p;
	int    nvxuy1_len =
	  o.Gen_Product_With_PreAlloc(nvx_len, nvx, nuy_len, nuy, &nvxuy1, 16);
	double nvxuy2_p[16], *nvxuy2 = nvxuy2_p;
	int    nvxuy2_len =
	  o.Gen_Product_With_PreAlloc(nvy_len, nvy, nux_len, nux, &nvxuy2, 16);
	double nvxuy_p[16], *nvxuy = nvxuy_p;
	int    nvxuy_len = o.Gen_Diff_With_PreAlloc(nvxuy1_len, nvxuy1, nvxuy2_len,
	                                            nvxuy2, &nvxuy, 16);
	double nvywz1_p[16], *nvywz1 = nvywz1_p;
	int    nvywz1_len =
	  o.Gen_Product_With_PreAlloc(nvy_len, nvy, nwz_len, nwz, &nvywz1, 16);
	double nvywz2_p[16], *nvywz2 = nvywz2_p;
	int    nvywz2_len =
	  o.Gen_Product_With_PreAlloc(nvz_len, nvz, nwy_len, nwy, &nvywz2, 16);
	double nvywz_p[16], *nvywz = nvywz_p;
	int    nvywz_len = o.Gen_Diff_With_PreAlloc(nvywz1_len, nvywz1, nvywz2_len,
	                                            nvywz2, &nvywz, 16);
	double nvxwz1_p[16], *nvxwz1 = nvxwz1_p;
	int    nvxwz1_len =
	  o.Gen_Product_With_PreAlloc(nvx_len, nvx, nwz_len, nwz, &nvxwz1, 16);
	double nvxwz2_p[16], *nvxwz2 = nvxwz2_p;
	int    nvxwz2_len =
	  o.Gen_Product_With_PreAlloc(nvz_len, nvz, nwx_len, nwx, &nvxwz2, 16);
	double nvxwz_p[16], *nvxwz = nvxwz_p;
	int    nvxwz_len = o.Gen_Diff_With_PreAlloc(nvxwz1_len, nvxwz1, nvxwz2_len,
	                                            nvxwz2, &nvxwz, 16);
	double nvxwy1_p[16], *nvxwy1 = nvxwy1_p;
	int    nvxwy1_len =
	  o.Gen_Product_With_PreAlloc(nvx_len, nvx, nwy_len, nwy, &nvxwy1, 16);
	double nvxwy2_p[16], *nvxwy2 = nvxwy2_p;
	int    nvxwy2_len =
	  o.Gen_Product_With_PreAlloc(nvy_len, nvy, nwx_len, nwx, &nvxwy2, 16);
	double nvxwy_p[16], *nvxwy = nvxwy_p;
	int    nvxwy_len = o.Gen_Diff_With_PreAlloc(nvxwy1_len, nvxwy1, nvxwy2_len,
	                                            nvxwy2, &nvxwy, 16);
	double p1a_p[16], *p1a = p1a_p;
	int    p1a_len = o.Gen_Scale_With_PreAlloc(nvx_len, nvx, ov1x, &p1a, 16);
	double p1b_p[16], *p1b = p1b_p;
	int    p1b_len = o.Gen_Scale_With_PreAlloc(nvy_len, nvy, ov1y, &p1b, 16);
	double p1c_p[16], *p1c = p1c_p;
	int    p1c_len = o.Gen_Scale_With_PreAlloc(nvz_len, nvz, ov1z, &p1c, 16);
	double p1ab_p[16], *p1ab = p1ab_p;
	int p1ab_len = o.Gen_Sum_With_PreAlloc(p1a_len, p1a, p1b_len, p1b, &p1ab, 16);
	double p1_p[16], *p1 = p1_p;
	int p1_len = o.Gen_Sum_With_PreAlloc(p1ab_len, p1ab, p1c_len, p1c, &p1, 16);
	double p2a_p[16], *p2a = p2a_p;
	int    p2a_len = o.Gen_Scale_With_PreAlloc(nwx_len, nwx, ow1x, &p2a, 16);
	double p2b_p[16], *p2b = p2b_p;
	int    p2b_len = o.Gen_Scale_With_PreAlloc(nwy_len, nwy, ow1y, &p2b, 16);
	double p2c_p[16], *p2c = p2c_p;
	int    p2c_len = o.Gen_Scale_With_PreAlloc(nwz_len, nwz, ow1z, &p2c, 16);
	double p2ab_p[16], *p2ab = p2ab_p;
	int p2ab_len = o.Gen_Sum_With_PreAlloc(p2a_len, p2a, p2b_len, p2b, &p2ab, 16);
	double p2_p[16], *p2 = p2_p;
	int p2_len = o.Gen_Sum_With_PreAlloc(p2ab_len, p2ab, p2c_len, p2c, &p2, 16);
	double p3a_p[16], *p3a = p3a_p;
	int    p3a_len = o.Gen_Scale_With_PreAlloc(nux_len, nux, ou1x, &p3a, 16);
	double p3b_p[16], *p3b = p3b_p;
	int    p3b_len = o.Gen_Scale_With_PreAlloc(nuy_len, nuy, ou1y, &p3b, 16);
	double p3c_p[16], *p3c = p3c_p;
	int    p3c_len = o.Gen_Scale_With_PreAlloc(nuz_len, nuz, ou1z, &p3c, 16);
	double p3ab_p[16], *p3ab = p3ab_p;
	int p3ab_len = o.Gen_Sum_With_PreAlloc(p3a_len, p3a, p3b_len, p3b, &p3ab, 16);
	double p3_p[16], *p3 = p3_p;
	int p3_len = o.Gen_Sum_With_PreAlloc(p3ab_len, p3ab, p3c_len, p3c, &p3, 16);
	double lxa_p[16], *lxa = lxa_p;
	int    lxa_len =
	  o.Gen_Product_With_PreAlloc(p1_len, p1, nwyuz_len, nwyuz, &lxa, 16);
	double lxb_p[16], *lxb = lxb_p;
	int    lxb_len =
	  o.Gen_Product_With_PreAlloc(p3_len, p3, nvywz_len, nvywz, &lxb, 16);
	double lxc_p[16], *lxc = lxc_p;
	int    lxc_len =
	  o.Gen_Product_With_PreAlloc(p2_len, p2, nvyuz_len, nvyuz, &lxc, 16);
	double lxab_p[16], *lxab = lxab_p;
	int lxab_len = o.Gen_Sum_With_PreAlloc(lxa_len, lxa, lxb_len, lxb, &lxab, 16);
	lambda_x_len = o.Gen_Diff_With_PreAlloc(lxab_len, lxab, lxc_len, lxc,
	                                        lambda_x, lambda_x_len);
	double lya_p[16], *lya = lya_p;
	int    lya_len =
	  o.Gen_Product_With_PreAlloc(p2_len, p2, nvxuz_len, nvxuz, &lya, 16);
	double lyb_p[16], *lyb = lyb_p;
	int    lyb_len =
	  o.Gen_Product_With_PreAlloc(p3_len, p3, nvxwz_len, nvxwz, &lyb, 16);
	double lyc_p[16], *lyc = lyc_p;
	int    lyc_len =
	  o.Gen_Product_With_PreAlloc(p1_len, p1, nwxuz_len, nwxuz, &lyc, 16);
	double lybc_p[16], *lybc = lybc_p;
	int lybc_len = o.Gen_Sum_With_PreAlloc(lyc_len, lyc, lyb_len, lyb, &lybc, 16);
	lambda_y_len = o.Gen_Diff_With_PreAlloc(lya_len, lya, lybc_len, lybc,
	                                        lambda_y, lambda_y_len);
	double lza_p[16], *lza = lza_p;
	int    lza_len =
	  o.Gen_Product_With_PreAlloc(p3_len, p3, nvxwy_len, nvxwy, &lza, 16);
	double lzb_p[16], *lzb = lzb_p;
	int    lzb_len =
	  o.Gen_Product_With_PreAlloc(p1_len, p1, nwxuy_len, nwxuy, &lzb, 16);
	double lzc_p[16], *lzc = lzc_p;
	int    lzc_len =
	  o.Gen_Product_With_PreAlloc(p2_len, p2, nvxuy_len, nvxuy, &lzc, 16);
	double lzab_p[16], *lzab = lzab_p;
	int lzab_len = o.Gen_Sum_With_PreAlloc(lza_len, lza, lzb_len, lzb, &lzab, 16);
	lambda_z_len = o.Gen_Diff_With_PreAlloc(lzab_len, lzab, lzc_len, lzc,
	                                        lambda_z, lambda_z_len);
	double da_p[16], *da = da_p;
	int    da_len =
	  o.Gen_Product_With_PreAlloc(nvx_len, nvx, nwyuz_len, nwyuz, &da, 16);
	double db_p[16], *db = db_p;
	int    db_len =
	  o.Gen_Product_With_PreAlloc(nvz_len, nvz, nwxuy_len, nwxuy, &db, 16);
	double dc_p[16], *dc = dc_p;
	int    dc_len =
	  o.Gen_Product_With_PreAlloc(nvy_len, nvy, nwxuz_len, nwxuz, &dc, 16);
	double dab_p[16], *dab = dab_p;
	int    dab_len = o.Gen_Sum_With_PreAlloc(da_len, da, db_len, db, &dab, 16);
	lambda_d_len =
	  o.Gen_Diff_With_PreAlloc(dab_len, dab, dc_len, dc, lambda_d, lambda_d_len);

	if (dab_p != dab)
		FreeDoubles(dab);
	if (dc_p != dc)
		FreeDoubles(dc);
	if (db_p != db)
		FreeDoubles(db);
	if (da_p != da)
		FreeDoubles(da);
	if (lzab_p != lzab)
		FreeDoubles(lzab);
	if (lzc_p != lzc)
		FreeDoubles(lzc);
	if (lzb_p != lzb)
		FreeDoubles(lzb);
	if (lza_p != lza)
		FreeDoubles(lza);
	if (lybc_p != lybc)
		FreeDoubles(lybc);
	if (lyc_p != lyc)
		FreeDoubles(lyc);
	if (lyb_p != lyb)
		FreeDoubles(lyb);
	if (lya_p != lya)
		FreeDoubles(lya);
	if (lxab_p != lxab)
		FreeDoubles(lxab);
	if (lxc_p != lxc)
		FreeDoubles(lxc);
	if (lxb_p != lxb)
		FreeDoubles(lxb);
	if (lxa_p != lxa)
		FreeDoubles(lxa);
	if (p3_p != p3)
		FreeDoubles(p3);
	if (p3ab_p != p3ab)
		FreeDoubles(p3ab);
	if (p3c_p != p3c)
		FreeDoubles(p3c);
	if (p3b_p != p3b)
		FreeDoubles(p3b);
	if (p3a_p != p3a)
		FreeDoubles(p3a);
	if (p2_p != p2)
		FreeDoubles(p2);
	if (p2ab_p != p2ab)
		FreeDoubles(p2ab);
	if (p2c_p != p2c)
		FreeDoubles(p2c);
	if (p2b_p != p2b)
		FreeDoubles(p2b);
	if (p2a_p != p2a)
		FreeDoubles(p2a);
	if (p1_p != p1)
		FreeDoubles(p1);
	if (p1ab_p != p1ab)
		FreeDoubles(p1ab);
	if (p1c_p != p1c)
		FreeDoubles(p1c);
	if (p1b_p != p1b)
		FreeDoubles(p1b);
	if (p1a_p != p1a)
		FreeDoubles(p1a);
	if (nvxwy_p != nvxwy)
		FreeDoubles(nvxwy);
	if (nvxwy2_p != nvxwy2)
		FreeDoubles(nvxwy2);
	if (nvxwy1_p != nvxwy1)
		FreeDoubles(nvxwy1);
	if (nvxwz_p != nvxwz)
		FreeDoubles(nvxwz);
	if (nvxwz2_p != nvxwz2)
		FreeDoubles(nvxwz2);
	if (nvxwz1_p != nvxwz1)
		FreeDoubles(nvxwz1);
	if (nvywz_p != nvywz)
		FreeDoubles(nvywz);
	if (nvywz2_p != nvywz2)
		FreeDoubles(nvywz2);
	if (nvywz1_p != nvywz1)
		FreeDoubles(nvywz1);
	if (nvxuy_p != nvxuy)
		FreeDoubles(nvxuy);
	if (nvxuy2_p != nvxuy2)
		FreeDoubles(nvxuy2);
	if (nvxuy1_p != nvxuy1)
		FreeDoubles(nvxuy1);
	if (nvxuz_p != nvxuz)
		FreeDoubles(nvxuz);
	if (nvxuz2_p != nvxuz2)
		FreeDoubles(nvxuz2);
	if (nvxuz1_p != nvxuz1)
		FreeDoubles(nvxuz1);
	if (nvyuz_p != nvyuz)
		FreeDoubles(nvyuz);
	if (nvyuz2_p != nvyuz2)
		FreeDoubles(nvyuz2);
	if (nvyuz1_p != nvyuz1)
		FreeDoubles(nvyuz1);
	if (nwxuy_p != nwxuy)
		FreeDoubles(nwxuy);
	if (nwxuy2_p != nwxuy2)
		FreeDoubles(nwxuy2);
	if (nwxuy1_p != nwxuy1)
		FreeDoubles(nwxuy1);
	if (nwxuz_p != nwxuz)
		FreeDoubles(nwxuz);
	if (nwxuz2_p != nwxuz2)
		FreeDoubles(nwxuz2);
	if (nwxuz1_p != nwxuz1)
		FreeDoubles(nwxuz1);
	if (nwyuz_p != nwyuz)
		FreeDoubles(nwyuz);
	if (nwyuz2_p != nwyuz2)
		FreeDoubles(nwyuz2);
	if (nwyuz1_p != nwyuz1)
		FreeDoubles(nwyuz1);
}

#elif defined(OMC_OFFSET_PRED)

inline bool lambda3d_LPI_filtered(double xa, double ya, double za, double xb,
                                  double yb, double zb, double xo, double yo,
                                  double zo, double xp, double yp, double zp,
                                  double xq, double yq, double zq,
                                  double &lambda_d, double &lambda_x,
                                  double &lambda_y, double &lambda_z,
                                  double &beta_x, double &beta_y,
                                  double &beta_z, double &max_var)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_lpi_filter);
	double xba      = xb - xa;
	double yba      = yb - ya;
	double zba      = zb - za;
	double xpo      = xp - xo;
	double ypo      = yp - yo;
	double zpo      = zp - zo;
	double xqo      = xq - xo;
	double yqo      = yq - yo;
	double zqo      = zq - zo;
	double xoa      = xo - xa;
	double yoa      = yo - ya;
	double zoa      = zo - za;
	double ypo_zqo  = ypo * zqo;
	double zpo_yqo  = zpo * yqo;
	double zpo_xqo  = zpo * xqo;
	double xpo_zqo  = xpo * zqo;
	double xpo_yqo  = xpo * yqo;
	double ypo_xqo  = ypo * xqo;
	double xnor     = ypo_zqo - zpo_yqo;
	double ynor     = zpo_xqo - xpo_zqo;
	double znor     = xpo_yqo - ypo_xqo;
	double x_ba_nor = xba * xnor;
	double y_ba_nor = yba * ynor;
	double z_ba_nor = zba * znor;
	double temp_0   = x_ba_nor + y_ba_nor;
	lambda_d        = temp_0 + z_ba_nor;
	double x_oa_nor = xoa * xnor;
	double y_oa_nor = yoa * ynor;
	double z_oa_nor = zoa * znor;
	double temp_1   = x_oa_nor + y_oa_nor;
	double n        = temp_1 + z_oa_nor;
	lambda_x        = n * xba;
	lambda_y        = n * yba;
	lambda_z        = n * zba;
	beta_x          = xa;
	beta_y          = ya;
	beta_z          = za;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(xba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xpo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ypo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zpo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xqo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yqo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zqo)) > max_var)
		max_var = _tmp_fabs;
	double lambda_d_eps = max_var;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= max_var;
	lambda_d_eps *= 4.884981308350689e-15;

	if ((_tmp_fabs = fabs(xoa)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yoa)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zoa)) > max_var)
		max_var = _tmp_fabs;

	return (lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps);
}

template <typename IT>
bool lambda3d_LPI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xo,
                           IT yo, IT zo, IT xp, IT yp, IT zp, IT xq, IT yq,
                           IT zq, IT &lambda_d, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z, IT &beta_x, IT &beta_y, IT &beta_z)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_lpi_interval);
	typename IT::Protector P;

	IT xba      = xb - xa;
	IT yba      = yb - ya;
	IT zba      = zb - za;
	IT xpo      = xp - xo;
	IT ypo      = yp - yo;
	IT zpo      = zp - zo;
	IT xqo      = xq - xo;
	IT yqo      = yq - yo;
	IT zqo      = zq - zo;
	IT xoa      = xo - xa;
	IT yoa      = yo - ya;
	IT zoa      = zo - za;
	IT ypo_zqo  = ypo * zqo;
	IT zpo_yqo  = zpo * yqo;
	IT zpo_xqo  = zpo * xqo;
	IT xpo_zqo  = xpo * zqo;
	IT xpo_yqo  = xpo * yqo;
	IT ypo_xqo  = ypo * xqo;
	IT xnor     = ypo_zqo - zpo_yqo;
	IT ynor     = zpo_xqo - xpo_zqo;
	IT znor     = xpo_yqo - ypo_xqo;
	IT x_ba_nor = xba * xnor;
	IT y_ba_nor = yba * ynor;
	IT z_ba_nor = zba * znor;
	IT temp_0   = x_ba_nor + y_ba_nor;
	lambda_d    = temp_0 + z_ba_nor;
	IT x_oa_nor = xoa * xnor;
	IT y_oa_nor = yoa * ynor;
	IT z_oa_nor = zoa * znor;
	IT temp_1   = x_oa_nor + y_oa_nor;
	IT n        = temp_1 + z_oa_nor;
	lambda_x    = n * xba;
	lambda_y    = n * yba;
	lambda_z    = n * zba;
	beta_x      = xa;
	beta_y      = ya;
	beta_z      = za;
	return lambda_d.is_sign_reliable();
}

template <typename ET>
void lambda3d_LPI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xo, ET yo,
                        ET zo, ET xp, ET yp, ET zp, ET xq, ET yq, ET zq,
                        ET &lambda_d, ET &lambda_x, ET &lambda_y, ET &lambda_z,
                        ET &beta_x, ET &beta_y, ET &beta_z)
{
	ET xba      = xb - xa;
	ET yba      = yb - ya;
	ET zba      = zb - za;
	ET xpo      = xp - xo;
	ET ypo      = yp - yo;
	ET zpo      = zp - zo;
	ET xqo      = xq - xo;
	ET yqo      = yq - yo;
	ET zqo      = zq - zo;
	ET xoa      = xo - xa;
	ET yoa      = yo - ya;
	ET zoa      = zo - za;
	ET ypo_zqo  = ypo * zqo;
	ET zpo_yqo  = zpo * yqo;
	ET zpo_xqo  = zpo * xqo;
	ET xpo_zqo  = xpo * zqo;
	ET xpo_yqo  = xpo * yqo;
	ET ypo_xqo  = ypo * xqo;
	ET xnor     = ypo_zqo - zpo_yqo;
	ET ynor     = zpo_xqo - xpo_zqo;
	ET znor     = xpo_yqo - ypo_xqo;
	ET x_ba_nor = xba * xnor;
	ET y_ba_nor = yba * ynor;
	ET z_ba_nor = zba * znor;
	ET temp_0   = x_ba_nor + y_ba_nor;
	lambda_d    = temp_0 + z_ba_nor;
	ET x_oa_nor = xoa * xnor;
	ET y_oa_nor = yoa * ynor;
	ET z_oa_nor = zoa * znor;
	ET temp_1   = x_oa_nor + y_oa_nor;
	ET n        = temp_1 + z_oa_nor;
	lambda_x    = n * xba;
	lambda_y    = n * yba;
	lambda_z    = n * zba;
	beta_x      = xa;
	beta_y      = ya;
	beta_z      = za;
}

inline void lambda3d_LPI_expansion(
  double xa, double ya, double za, double xb, double yb, double zb, double xo,
  double yo, double zo, double xp, double yp, double zp, double xq, double yq,
  double zq, double **lambda_d, int &lambda_d_len, double **lambda_x,
  int &lambda_x_len, double **lambda_y, int &lambda_y_len, double **lambda_z,
  int &lambda_z_len, double &beta_x, double &beta_y, double &beta_z)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_lpi_expansion);
	expansionObject o;
	double          xba[2];
	o.two_Diff(xb, xa, xba);
	double yba[2];
	o.two_Diff(yb, ya, yba);
	double zba[2];
	o.two_Diff(zb, za, zba);
	double xpo[2];
	o.two_Diff(xp, xo, xpo);
	double ypo[2];
	o.two_Diff(yp, yo, ypo);
	double zpo[2];
	o.two_Diff(zp, zo, zpo);
	double xqo[2];
	o.two_Diff(xq, xo, xqo);
	double yqo[2];
	o.two_Diff(yq, yo, yqo);
	double zqo[2];
	o.two_Diff(zq, zo, zqo);
	double xoa[2];
	o.two_Diff(xo, xa, xoa);
	double yoa[2];
	o.two_Diff(yo, ya, yoa);
	double zoa[2];
	o.two_Diff(zo, za, zoa);
	double ypo_zqo[8];
	int    ypo_zqo_len = o.Gen_Product(2, ypo, 2, zqo, ypo_zqo);
	double zpo_yqo[8];
	int    zpo_yqo_len = o.Gen_Product(2, zpo, 2, yqo, zpo_yqo);
	double zpo_xqo[8];
	int    zpo_xqo_len = o.Gen_Product(2, zpo, 2, xqo, zpo_xqo);
	double xpo_zqo[8];
	int    xpo_zqo_len = o.Gen_Product(2, xpo, 2, zqo, xpo_zqo);
	double xpo_yqo[8];
	int    xpo_yqo_len = o.Gen_Product(2, xpo, 2, yqo, xpo_yqo);
	double ypo_xqo[8];
	int    ypo_xqo_len = o.Gen_Product(2, ypo, 2, xqo, ypo_xqo);
	double xnor[16];
	int xnor_len = o.Gen_Diff(ypo_zqo_len, ypo_zqo, zpo_yqo_len, zpo_yqo, xnor);
	double ynor[16];
	int ynor_len = o.Gen_Diff(zpo_xqo_len, zpo_xqo, xpo_zqo_len, xpo_zqo, ynor);
	double znor[16];
	int znor_len = o.Gen_Diff(xpo_yqo_len, xpo_yqo, ypo_xqo_len, ypo_xqo, znor);
	double x_ba_nor[64];
	int    x_ba_nor_len = o.Gen_Product(2, xba, xnor_len, xnor, x_ba_nor);
	double y_ba_nor[64];
	int    y_ba_nor_len = o.Gen_Product(2, yba, ynor_len, ynor, y_ba_nor);
	double z_ba_nor[64];
	int    z_ba_nor_len = o.Gen_Product(2, zba, znor_len, znor, z_ba_nor);
	double temp_0[128];
	int    temp_0_len =
	  o.Gen_Sum(x_ba_nor_len, x_ba_nor, y_ba_nor_len, y_ba_nor, temp_0);
	lambda_d_len = o.Gen_Sum_With_PreAlloc(temp_0_len, temp_0, z_ba_nor_len,
	                                       z_ba_nor, lambda_d, lambda_d_len);
	double x_oa_nor[64];
	int    x_oa_nor_len = o.Gen_Product(2, xoa, xnor_len, xnor, x_oa_nor);
	double y_oa_nor[64];
	int    y_oa_nor_len = o.Gen_Product(2, yoa, ynor_len, ynor, y_oa_nor);
	double z_oa_nor[64];
	int    z_oa_nor_len = o.Gen_Product(2, zoa, znor_len, znor, z_oa_nor);
	double temp_1[128];
	int    temp_1_len =
	  o.Gen_Sum(x_oa_nor_len, x_oa_nor, y_oa_nor_len, y_oa_nor, temp_1);
	double n_p[128], *n = n_p;
	int    n_len = o.Gen_Sum_With_PreAlloc(temp_1_len, temp_1, z_oa_nor_len,
	                                       z_oa_nor, &n, 128);
	lambda_x_len =
	  o.Gen_Product_With_PreAlloc(n_len, n, 2, xba, lambda_x, lambda_x_len);
	lambda_y_len =
	  o.Gen_Product_With_PreAlloc(n_len, n, 2, yba, lambda_y, lambda_y_len);
	lambda_z_len =
	  o.Gen_Product_With_PreAlloc(n_len, n, 2, zba, lambda_z, lambda_z_len);
	beta_x = xa;
	beta_y = ya;
	beta_z = za;

	if (n_p != n)
		FreeDoubles(n);
}

inline bool lambda3d_SSI_filtered(double xa, double ya, double za, double xb,
                                  double yb, double zb, double xp, double yp,
                                  double xq, double yq, double &lambda_d,
                                  double &lambda_x, double &lambda_y,
                                  double &lambda_z, double &beta_x,
                                  double &beta_y, double &beta_z,
                                  double &max_var)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_ssi_filter);

	double xap = xa - xp;
	double yap = ya - yp;
	double yqp = yq - yp;
	double xqp = xq - xp;
	double xba = xb - xa;
	double yba = yb - ya;
	double zba = zb - za;
	double c1  = xap * yqp;
	double c2  = xqp * yap;
	double n   = c1 - c2;
	double c3  = xba * yqp;
	double c4  = xqp * yba;
	lambda_d   = c4 - c3;
	lambda_x   = n * xba;
	lambda_y   = n * yba;
	lambda_z   = n * zba;
	beta_x     = xa;
	beta_y     = ya;
	beta_z     = za;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(yqp)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xqp)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yba)) > max_var)
		max_var = _tmp_fabs;
	double lambda_d_eps = max_var;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= 8.881784197001252e-16;

	if ((_tmp_fabs = fabs(xap)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yap)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zba)) > max_var)
		max_var = _tmp_fabs;

	return (lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps);
}

template <typename IT>
bool lambda3d_SSI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xp,
                           IT yp, IT xq, IT yq, IT &lambda_d, IT &lambda_x,
                           IT &lambda_y, IT &lambda_z, IT &beta_x, IT &beta_y,
                           IT &beta_z)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_ssi_interval);
	typename IT::Protector P;

	IT xap   = xa - xp;
	IT yap   = ya - yp;
	IT yqp   = yq - yp;
	IT xqp   = xq - xp;
	IT xba   = xb - xa;
	IT yba   = yb - ya;
	IT zba   = zb - za;
	IT c1    = xap * yqp;
	IT c2    = xqp * yap;
	IT n     = c1 - c2;
	IT c3    = xba * yqp;
	IT c4    = xqp * yba;
	lambda_d = c4 - c3;
	lambda_x = n * xba;
	lambda_y = n * yba;
	lambda_z = n * zba;
	beta_x   = xa;
	beta_y   = ya;
	beta_z   = za;
	return lambda_d.is_sign_reliable();
}

template <typename ET>
void lambda3d_SSI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xp, ET yp,
                        ET xq, ET yq, ET &lambda_d, ET &lambda_x, ET &lambda_y,
                        ET &lambda_z, ET &beta_x, ET &beta_y, ET &beta_z)
{
	ET xap   = xa - xp;
	ET yap   = ya - yp;
	ET yqp   = yq - yp;
	ET xqp   = xq - xp;
	ET xba   = xb - xa;
	ET yba   = yb - ya;
	ET zba   = zb - za;
	ET c1    = xap * yqp;
	ET c2    = xqp * yap;
	ET n     = c1 - c2;
	ET c3    = xba * yqp;
	ET c4    = xqp * yba;
	lambda_d = c4 - c3;
	lambda_x = n * xba;
	lambda_y = n * yba;
	lambda_z = n * zba;
	beta_x   = xa;
	beta_y   = ya;
	beta_z   = za;
}

inline void lambda3d_SSI_expansion(double xa, double ya, double za, double xb,
                                   double yb, double zb, double xp, double yp,
                                   double xq, double yq, double **lambda_d,
                                   int &lambda_d_len, double **lambda_x,
                                   int &lambda_x_len, double **lambda_y,
                                   int &lambda_y_len, double **lambda_z,
                                   int &lambda_z_len, double &beta_x,
                                   double &beta_y, double &beta_z)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_ssi_expansion);
	expansionObject o;
	double          xap[2];
	o.two_Diff(xa, xp, xap);
	double yap[2];
	o.two_Diff(ya, yp, yap);
	double yqp[2];
	o.two_Diff(yq, yp, yqp);
	double xqp[2];
	o.two_Diff(xq, xp, xqp);
	double xba[2];
	o.two_Diff(xb, xa, xba);
	double yba[2];
	o.two_Diff(yb, ya, yba);
	double zba[2];
	o.two_Diff(zb, za, zba);
	double c1[8];
	int    c1_len = o.Gen_Product(2, xap, 2, yqp, c1);
	double c2[8];
	int    c2_len = o.Gen_Product(2, xqp, 2, yap, c2);
	double n[16];
	int    n_len = o.Gen_Diff(c1_len, c1, c2_len, c2, n);
	double c3[8];
	int    c3_len = o.Gen_Product(2, xba, 2, yqp, c3);
	double c4[8];
	int    c4_len = o.Gen_Product(2, xqp, 2, yba, c4);
	lambda_d_len =
	  o.Gen_Diff_With_PreAlloc(c4_len, c4, c3_len, c3, lambda_d, lambda_d_len);
	lambda_x_len =
	  o.Gen_Product_With_PreAlloc(n_len, n, 2, xba, lambda_x, lambda_x_len);
	lambda_y_len =
	  o.Gen_Product_With_PreAlloc(n_len, n, 2, yba, lambda_y, lambda_y_len);
	lambda_z_len =
	  o.Gen_Product_With_PreAlloc(n_len, n, 2, zba, lambda_z, lambda_z_len);
	beta_x = xa;
	beta_y = ya;
	beta_z = za;
}

inline bool lambda3d_TPI_filtered(
  double xa, double ya, double za, double xb, double yb, double zb, double xc,
  double yc, double zc, double xo, double yo, double zo, double xp, double yp,
  double zp, double xq, double yq, double zq, double xr, double yr, double zr,
  double xs, double ys, double zs, double xt, double yt, double zt,
  double &lambda_d, double &lambda_x, double &lambda_y, double &lambda_z,
  double &beta_x, double &beta_y, double &beta_z, double &max_var)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_tpi_filter);
	double xpo            = xp - xo;
	double ypo            = yp - yo;
	double zpo            = zp - zo;
	double xqo            = xq - xo;
	double yqo            = yq - yo;
	double zqo            = zq - zo;
	double xsr            = xs - xr;
	double ysr            = ys - yr;
	double zsr            = zs - zr;
	double xtr            = xt - xr;
	double ytr            = yt - yr;
	double ztr            = zt - zr;
	double ypo_zqo        = ypo * zqo;
	double zpo_yqo        = zpo * yqo;
	double zpo_xqo        = zpo * xqo;
	double xpo_zqo        = xpo * zqo;
	double xpo_yqo        = xpo * yqo;
	double ypo_xqo        = ypo * xqo;
	double x_nor_opq      = ypo_zqo - zpo_yqo;
	double y_nor_opq      = zpo_xqo - xpo_zqo;
	double z_nor_opq      = xpo_yqo - ypo_xqo;
	double ysr_ztr        = ysr * ztr;
	double zsr_ytr        = zsr * ytr;
	double zsr_xtr        = zsr * xtr;
	double xsr_ztr        = xsr * ztr;
	double xsr_ytr        = xsr * ytr;
	double ysr_xtr        = ysr * xtr;
	double x_nor_rst      = ysr_ztr - zsr_ytr;
	double y_nor_rst      = zsr_xtr - xsr_ztr;
	double z_nor_rst      = xsr_ytr - ysr_xtr;
	double xoa            = xo - xa;
	double yoa            = yo - ya;
	double zoa            = zo - za;
	double xra            = xr - xa;
	double yra            = yr - ya;
	double zra            = zr - za;
	double xba            = xb - xa;
	double yba            = yb - ya;
	double zba            = zb - za;
	double xca            = xc - xa;
	double yca            = yc - ya;
	double zca            = zc - za;
	double temp_0         = xoa * x_nor_opq;
	double temp_1         = yoa * y_nor_opq;
	double temp_2         = zoa * z_nor_opq;
	double temp_3         = temp_0 + temp_1;
	double oa_dot_nor_opq = temp_2 + temp_3;
	double temp_4         = xra * x_nor_rst;
	double temp_5         = yra * y_nor_rst;
	double temp_6         = zra * z_nor_rst;
	double temp_7         = temp_4 + temp_5;
	double ra_dot_nor_rst = temp_6 + temp_7;
	double temp_8         = xba * x_nor_opq;
	double temp_9         = yba * y_nor_opq;
	double temp_10        = zba * z_nor_opq;
	double temp_11        = temp_8 + temp_9;
	double ba_dot_nor_opq = temp_10 + temp_11;
	double temp_12        = xca * x_nor_opq;
	double temp_13        = yca * y_nor_opq;
	double temp_14        = zca * z_nor_opq;
	double temp_15        = temp_12 + temp_13;
	double ca_dot_nor_opq = temp_14 + temp_15;
	double temp_16        = xba * x_nor_rst;
	double temp_17        = yba * y_nor_rst;
	double temp_18        = zba * z_nor_rst;
	double temp_19        = temp_16 + temp_17;
	double ba_dot_nor_rst = temp_18 + temp_19;
	double temp_20        = xca * x_nor_rst;
	double temp_21        = yca * y_nor_rst;
	double temp_22        = zca * z_nor_rst;
	double temp_23        = temp_20 + temp_21;
	double ca_dot_nor_rst = temp_22 + temp_23;
	double temp_100       = ba_dot_nor_opq * ca_dot_nor_rst;
	double temp_101       = ca_dot_nor_opq * ba_dot_nor_rst;
	lambda_d              = temp_100 - temp_101;
	double temp_102       = oa_dot_nor_opq * ca_dot_nor_rst;
	double temp_103       = ca_dot_nor_opq * ra_dot_nor_rst;
	double det_sub0       = temp_102 - temp_103;
	double temp_105       = ba_dot_nor_opq * ra_dot_nor_rst;
	double temp_106       = oa_dot_nor_opq * ba_dot_nor_rst;
	double det_sub1       = temp_105 - temp_106;
	double xu             = det_sub0 * xba;
	double yu             = det_sub0 * yba;
	double zu             = det_sub0 * zba;
	double xv             = det_sub1 * xca;
	double yv             = det_sub1 * yca;
	double zv             = det_sub1 * zca;
	lambda_x              = xu + xv;
	lambda_y              = yu + yv;
	lambda_z              = zu + zv;
	beta_x                = xa;
	beta_y                = ya;
	beta_z                = za;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(xpo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ypo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zpo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xqo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yqo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zqo)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xsr)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ysr)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zsr)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xtr)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ytr)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ztr)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xca)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yca)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zca)) > max_var)
		max_var = _tmp_fabs;
	double lambda_d_eps = max_var;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= lambda_d_eps;
	lambda_d_eps *= max_var;
	lambda_d_eps *= max_var;
	lambda_d_eps *= 1.3145040611561864e-13;

	if ((_tmp_fabs = fabs(xoa)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yoa)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zoa)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xra)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yra)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zra)) > max_var)
		max_var = _tmp_fabs;

	return (lambda_d > lambda_d_eps || lambda_d < -lambda_d_eps);
}

template <typename IT>
bool lambda3d_TPI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xc,
                           IT yc, IT zc, IT xo, IT yo, IT zo, IT xp, IT yp,
                           IT zp, IT xq, IT yq, IT zq, IT xr, IT yr, IT zr,
                           IT xs, IT ys, IT zs, IT xt, IT yt, IT zt,
                           IT &lambda_d, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z, IT &beta_x, IT &beta_y, IT &beta_z)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_tpi_interval);
	typename IT::Protector P;

	IT xpo            = xp - xo;
	IT ypo            = yp - yo;
	IT zpo            = zp - zo;
	IT xqo            = xq - xo;
	IT yqo            = yq - yo;
	IT zqo            = zq - zo;
	IT xsr            = xs - xr;
	IT ysr            = ys - yr;
	IT zsr            = zs - zr;
	IT xtr            = xt - xr;
	IT ytr            = yt - yr;
	IT ztr            = zt - zr;
	IT ypo_zqo        = ypo * zqo;
	IT zpo_yqo        = zpo * yqo;
	IT zpo_xqo        = zpo * xqo;
	IT xpo_zqo        = xpo * zqo;
	IT xpo_yqo        = xpo * yqo;
	IT ypo_xqo        = ypo * xqo;
	IT x_nor_opq      = ypo_zqo - zpo_yqo;
	IT y_nor_opq      = zpo_xqo - xpo_zqo;
	IT z_nor_opq      = xpo_yqo - ypo_xqo;
	IT ysr_ztr        = ysr * ztr;
	IT zsr_ytr        = zsr * ytr;
	IT zsr_xtr        = zsr * xtr;
	IT xsr_ztr        = xsr * ztr;
	IT xsr_ytr        = xsr * ytr;
	IT ysr_xtr        = ysr * xtr;
	IT x_nor_rst      = ysr_ztr - zsr_ytr;
	IT y_nor_rst      = zsr_xtr - xsr_ztr;
	IT z_nor_rst      = xsr_ytr - ysr_xtr;
	IT xoa            = xo - xa;
	IT yoa            = yo - ya;
	IT zoa            = zo - za;
	IT xra            = xr - xa;
	IT yra            = yr - ya;
	IT zra            = zr - za;
	IT xba            = xb - xa;
	IT yba            = yb - ya;
	IT zba            = zb - za;
	IT xca            = xc - xa;
	IT yca            = yc - ya;
	IT zca            = zc - za;
	IT temp_0         = xoa * x_nor_opq;
	IT temp_1         = yoa * y_nor_opq;
	IT temp_2         = zoa * z_nor_opq;
	IT temp_3         = temp_0 + temp_1;
	IT oa_dot_nor_opq = temp_2 + temp_3;
	IT temp_4         = xra * x_nor_rst;
	IT temp_5         = yra * y_nor_rst;
	IT temp_6         = zra * z_nor_rst;
	IT temp_7         = temp_4 + temp_5;
	IT ra_dot_nor_rst = temp_6 + temp_7;
	IT temp_8         = xba * x_nor_opq;
	IT temp_9         = yba * y_nor_opq;
	IT temp_10        = zba * z_nor_opq;
	IT temp_11        = temp_8 + temp_9;
	IT ba_dot_nor_opq = temp_10 + temp_11;
	IT temp_12        = xca * x_nor_opq;
	IT temp_13        = yca * y_nor_opq;
	IT temp_14        = zca * z_nor_opq;
	IT temp_15        = temp_12 + temp_13;
	IT ca_dot_nor_opq = temp_14 + temp_15;
	IT temp_16        = xba * x_nor_rst;
	IT temp_17        = yba * y_nor_rst;
	IT temp_18        = zba * z_nor_rst;
	IT temp_19        = temp_16 + temp_17;
	IT ba_dot_nor_rst = temp_18 + temp_19;
	IT temp_20        = xca * x_nor_rst;
	IT temp_21        = yca * y_nor_rst;
	IT temp_22        = zca * z_nor_rst;
	IT temp_23        = temp_20 + temp_21;
	IT ca_dot_nor_rst = temp_22 + temp_23;
	IT temp_100       = ba_dot_nor_opq * ca_dot_nor_rst;
	IT temp_101       = ca_dot_nor_opq * ba_dot_nor_rst;
	lambda_d          = temp_100 - temp_101;
	IT temp_102       = oa_dot_nor_opq * ca_dot_nor_rst;
	IT temp_103       = ca_dot_nor_opq * ra_dot_nor_rst;
	IT det_sub0       = temp_102 - temp_103;
	IT temp_105       = ba_dot_nor_opq * ra_dot_nor_rst;
	IT temp_106       = oa_dot_nor_opq * ba_dot_nor_rst;
	IT det_sub1       = temp_105 - temp_106;
	IT xu             = det_sub0 * xba;
	IT yu             = det_sub0 * yba;
	IT zu             = det_sub0 * zba;
	IT xv             = det_sub1 * xca;
	IT yv             = det_sub1 * yca;
	IT zv             = det_sub1 * zca;
	lambda_x          = xu + xv;
	lambda_y          = yu + yv;
	lambda_z          = zu + zv;
	beta_x            = xa;
	beta_y            = ya;
	beta_z            = za;
	return lambda_d.is_sign_reliable();
}

template <typename ET>
void lambda3d_TPI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xc, ET yc,
                        ET zc, ET xo, ET yo, ET zo, ET xp, ET yp, ET zp, ET xq,
                        ET yq, ET zq, ET xr, ET yr, ET zr, ET xs, ET ys, ET zs,
                        ET xt, ET yt, ET zt, ET &lambda_d, ET &lambda_x,
                        ET &lambda_y, ET &lambda_z, ET &beta_x, ET &beta_y,
                        ET &beta_z)
{
	ET xpo            = xp - xo;
	ET ypo            = yp - yo;
	ET zpo            = zp - zo;
	ET xqo            = xq - xo;
	ET yqo            = yq - yo;
	ET zqo            = zq - zo;
	ET xsr            = xs - xr;
	ET ysr            = ys - yr;
	ET zsr            = zs - zr;
	ET xtr            = xt - xr;
	ET ytr            = yt - yr;
	ET ztr            = zt - zr;
	ET ypo_zqo        = ypo * zqo;
	ET zpo_yqo        = zpo * yqo;
	ET zpo_xqo        = zpo * xqo;
	ET xpo_zqo        = xpo * zqo;
	ET xpo_yqo        = xpo * yqo;
	ET ypo_xqo        = ypo * xqo;
	ET x_nor_opq      = ypo_zqo - zpo_yqo;
	ET y_nor_opq      = zpo_xqo - xpo_zqo;
	ET z_nor_opq      = xpo_yqo - ypo_xqo;
	ET ysr_ztr        = ysr * ztr;
	ET zsr_ytr        = zsr * ytr;
	ET zsr_xtr        = zsr * xtr;
	ET xsr_ztr        = xsr * ztr;
	ET xsr_ytr        = xsr * ytr;
	ET ysr_xtr        = ysr * xtr;
	ET x_nor_rst      = ysr_ztr - zsr_ytr;
	ET y_nor_rst      = zsr_xtr - xsr_ztr;
	ET z_nor_rst      = xsr_ytr - ysr_xtr;
	ET xoa            = xo - xa;
	ET yoa            = yo - ya;
	ET zoa            = zo - za;
	ET xra            = xr - xa;
	ET yra            = yr - ya;
	ET zra            = zr - za;
	ET xba            = xb - xa;
	ET yba            = yb - ya;
	ET zba            = zb - za;
	ET xca            = xc - xa;
	ET yca            = yc - ya;
	ET zca            = zc - za;
	ET temp_0         = xoa * x_nor_opq;
	ET temp_1         = yoa * y_nor_opq;
	ET temp_2         = zoa * z_nor_opq;
	ET temp_3         = temp_0 + temp_1;
	ET oa_dot_nor_opq = temp_2 + temp_3;
	ET temp_4         = xra * x_nor_rst;
	ET temp_5         = yra * y_nor_rst;
	ET temp_6         = zra * z_nor_rst;
	ET temp_7         = temp_4 + temp_5;
	ET ra_dot_nor_rst = temp_6 + temp_7;
	ET temp_8         = xba * x_nor_opq;
	ET temp_9         = yba * y_nor_opq;
	ET temp_10        = zba * z_nor_opq;
	ET temp_11        = temp_8 + temp_9;
	ET ba_dot_nor_opq = temp_10 + temp_11;
	ET temp_12        = xca * x_nor_opq;
	ET temp_13        = yca * y_nor_opq;
	ET temp_14        = zca * z_nor_opq;
	ET temp_15        = temp_12 + temp_13;
	ET ca_dot_nor_opq = temp_14 + temp_15;
	ET temp_16        = xba * x_nor_rst;
	ET temp_17        = yba * y_nor_rst;
	ET temp_18        = zba * z_nor_rst;
	ET temp_19        = temp_16 + temp_17;
	ET ba_dot_nor_rst = temp_18 + temp_19;
	ET temp_20        = xca * x_nor_rst;
	ET temp_21        = yca * y_nor_rst;
	ET temp_22        = zca * z_nor_rst;
	ET temp_23        = temp_20 + temp_21;
	ET ca_dot_nor_rst = temp_22 + temp_23;
	ET temp_100       = ba_dot_nor_opq * ca_dot_nor_rst;
	ET temp_101       = ca_dot_nor_opq * ba_dot_nor_rst;
	lambda_d          = temp_100 - temp_101;
	ET temp_102       = oa_dot_nor_opq * ca_dot_nor_rst;
	ET temp_103       = ca_dot_nor_opq * ra_dot_nor_rst;
	ET det_sub0       = temp_102 - temp_103;
	ET temp_105       = ba_dot_nor_opq * ra_dot_nor_rst;
	ET temp_106       = oa_dot_nor_opq * ba_dot_nor_rst;
	ET det_sub1       = temp_105 - temp_106;
	ET xu             = det_sub0 * xba;
	ET yu             = det_sub0 * yba;
	ET zu             = det_sub0 * zba;
	ET xv             = det_sub1 * xca;
	ET yv             = det_sub1 * yca;
	ET zv             = det_sub1 * zca;
	lambda_x          = xu + xv;
	lambda_y          = yu + yv;
	lambda_z          = zu + zv;
	beta_x            = xa;
	beta_y            = ya;
	beta_z            = za;
}

inline void lambda3d_TPI_expansion(
  double xa, double ya, double za, double xb, double yb, double zb, double xc,
  double yc, double zc, double xo, double yo, double zo, double xp, double yp,
  double zp, double xq, double yq, double zq, double xr, double yr, double zr,
  double xs, double ys, double zs, double xt, double yt, double zt,
  double **lambda_d, int &lambda_d_len, double **lambda_x, int &lambda_x_len,
  double **lambda_y, int &lambda_y_len, double **lambda_z, int &lambda_z_len,
  double &beta_x, double &beta_y, double &beta_z)
{
	OMC_PRED_PROFILE_INC_IP_TOTAL(PredicateNames::_tpi_expansion);
	expansionObject o;
	double          xpo[2];
	o.two_Diff(xp, xo, xpo);
	double ypo[2];
	o.two_Diff(yp, yo, ypo);
	double zpo[2];
	o.two_Diff(zp, zo, zpo);
	double xqo[2];
	o.two_Diff(xq, xo, xqo);
	double yqo[2];
	o.two_Diff(yq, yo, yqo);
	double zqo[2];
	o.two_Diff(zq, zo, zqo);
	double xsr[2];
	o.two_Diff(xs, xr, xsr);
	double ysr[2];
	o.two_Diff(ys, yr, ysr);
	double zsr[2];
	o.two_Diff(zs, zr, zsr);
	double xtr[2];
	o.two_Diff(xt, xr, xtr);
	double ytr[2];
	o.two_Diff(yt, yr, ytr);
	double ztr[2];
	o.two_Diff(zt, zr, ztr);
	double ypo_zqo[8];
	int    ypo_zqo_len = o.Gen_Product(2, ypo, 2, zqo, ypo_zqo);
	double zpo_yqo[8];
	int    zpo_yqo_len = o.Gen_Product(2, zpo, 2, yqo, zpo_yqo);
	double zpo_xqo[8];
	int    zpo_xqo_len = o.Gen_Product(2, zpo, 2, xqo, zpo_xqo);
	double xpo_zqo[8];
	int    xpo_zqo_len = o.Gen_Product(2, xpo, 2, zqo, xpo_zqo);
	double xpo_yqo[8];
	int    xpo_yqo_len = o.Gen_Product(2, xpo, 2, yqo, xpo_yqo);
	double ypo_xqo[8];
	int    ypo_xqo_len = o.Gen_Product(2, ypo, 2, xqo, ypo_xqo);
	double x_nor_opq[16];
	int    x_nor_opq_len =
	  o.Gen_Diff(ypo_zqo_len, ypo_zqo, zpo_yqo_len, zpo_yqo, x_nor_opq);
	double y_nor_opq[16];
	int    y_nor_opq_len =
	  o.Gen_Diff(zpo_xqo_len, zpo_xqo, xpo_zqo_len, xpo_zqo, y_nor_opq);
	double z_nor_opq[16];
	int    z_nor_opq_len =
	  o.Gen_Diff(xpo_yqo_len, xpo_yqo, ypo_xqo_len, ypo_xqo, z_nor_opq);
	double ysr_ztr[8];
	int    ysr_ztr_len = o.Gen_Product(2, ysr, 2, ztr, ysr_ztr);
	double zsr_ytr[8];
	int    zsr_ytr_len = o.Gen_Product(2, zsr, 2, ytr, zsr_ytr);
	double zsr_xtr[8];
	int    zsr_xtr_len = o.Gen_Product(2, zsr, 2, xtr, zsr_xtr);
	double xsr_ztr[8];
	int    xsr_ztr_len = o.Gen_Product(2, xsr, 2, ztr, xsr_ztr);
	double xsr_ytr[8];
	int    xsr_ytr_len = o.Gen_Product(2, xsr, 2, ytr, xsr_ytr);
	double ysr_xtr[8];
	int    ysr_xtr_len = o.Gen_Product(2, ysr, 2, xtr, ysr_xtr);
	double x_nor_rst[16];
	int    x_nor_rst_len =
	  o.Gen_Diff(ysr_ztr_len, ysr_ztr, zsr_ytr_len, zsr_ytr, x_nor_rst);
	double y_nor_rst[16];
	int    y_nor_rst_len =
	  o.Gen_Diff(zsr_xtr_len, zsr_xtr, xsr_ztr_len, xsr_ztr, y_nor_rst);
	double z_nor_rst[16];
	int    z_nor_rst_len =
	  o.Gen_Diff(xsr_ytr_len, xsr_ytr, ysr_xtr_len, ysr_xtr, z_nor_rst);
	double xoa[2];
	o.two_Diff(xo, xa, xoa);
	double yoa[2];
	o.two_Diff(yo, ya, yoa);
	double zoa[2];
	o.two_Diff(zo, za, zoa);
	double xra[2];
	o.two_Diff(xr, xa, xra);
	double yra[2];
	o.two_Diff(yr, ya, yra);
	double zra[2];
	o.two_Diff(zr, za, zra);
	double xba[2];
	o.two_Diff(xb, xa, xba);
	double yba[2];
	o.two_Diff(yb, ya, yba);
	double zba[2];
	o.two_Diff(zb, za, zba);
	double xca[2];
	o.two_Diff(xc, xa, xca);
	double yca[2];
	o.two_Diff(yc, ya, yca);
	double zca[2];
	o.two_Diff(zc, za, zca);
	double temp_0_p[32], *temp_0 = temp_0_p;
	int    temp_0_len =
	  o.Gen_Product_With_PreAlloc(2, xoa, x_nor_opq_len, x_nor_opq, &temp_0, 32);
	double temp_1_p[32], *temp_1 = temp_1_p;
	int    temp_1_len =
	  o.Gen_Product_With_PreAlloc(2, yoa, y_nor_opq_len, y_nor_opq, &temp_1, 32);
	double temp_2_p[32], *temp_2 = temp_2_p;
	int    temp_2_len =
	  o.Gen_Product_With_PreAlloc(2, zoa, z_nor_opq_len, z_nor_opq, &temp_2, 32);
	double temp_3_p[32], *temp_3 = temp_3_p;
	int    temp_3_len = o.Gen_Sum_With_PreAlloc(temp_0_len, temp_0, temp_1_len,
	                                            temp_1, &temp_3, 32);
	double oa_dot_nor_opq_p[32], *oa_dot_nor_opq = oa_dot_nor_opq_p;
	int    oa_dot_nor_opq_len = o.Gen_Sum_With_PreAlloc(
    temp_2_len, temp_2, temp_3_len, temp_3, &oa_dot_nor_opq, 32);
	double temp_4_p[32], *temp_4 = temp_4_p;
	int    temp_4_len =
	  o.Gen_Product_With_PreAlloc(2, xra, x_nor_rst_len, x_nor_rst, &temp_4, 32);
	double temp_5_p[32], *temp_5 = temp_5_p;
	int    temp_5_len =
	  o.Gen_Product_With_PreAlloc(2, yra, y_nor_rst_len, y_nor_rst, &temp_5, 32);
	double temp_6_p[32], *temp_6 = temp_6_p;
	int    temp_6_len =
	  o.Gen_Product_With_PreAlloc(2, zra, z_nor_rst_len, z_nor_rst, &temp_6, 32);
	double temp_7_p[32], *temp_7 = temp_7_p;
	int    temp_7_len = o.Gen_Sum_With_PreAlloc(temp_4_len, temp_4, temp_5_len,
	                                            temp_5, &temp_7, 32);
	double ra_dot_nor_rst_p[32], *ra_dot_nor_rst = ra_dot_nor_rst_p;
	int    ra_dot_nor_rst_len = o.Gen_Sum_With_PreAlloc(
    temp_6_len, temp_6, temp_7_len, temp_7, &ra_dot_nor_rst, 32);
	double temp_8_p[32], *temp_8 = temp_8_p;
	int    temp_8_len =
	  o.Gen_Product_With_PreAlloc(2, xba, x_nor_opq_len, x_nor_opq, &temp_8, 32);
	double temp_9_p[32], *temp_9 = temp_9_p;
	int    temp_9_len =
	  o.Gen_Product_With_PreAlloc(2, yba, y_nor_opq_len, y_nor_opq, &temp_9, 32);
	double temp_10_p[32], *temp_10 = temp_10_p;
	int    temp_10_len =
	  o.Gen_Product_With_PreAlloc(2, zba, z_nor_opq_len, z_nor_opq, &temp_10, 32);
	double temp_11_p[32], *temp_11 = temp_11_p;
	int    temp_11_len = o.Gen_Sum_With_PreAlloc(temp_8_len, temp_8, temp_9_len,
	                                             temp_9, &temp_11, 32);
	double ba_dot_nor_opq_p[32], *ba_dot_nor_opq = ba_dot_nor_opq_p;
	int    ba_dot_nor_opq_len = o.Gen_Sum_With_PreAlloc(
    temp_10_len, temp_10, temp_11_len, temp_11, &ba_dot_nor_opq, 32);
	double temp_12_p[32], *temp_12 = temp_12_p;
	int    temp_12_len =
	  o.Gen_Product_With_PreAlloc(2, xca, x_nor_opq_len, x_nor_opq, &temp_12, 32);
	double temp_13_p[32], *temp_13 = temp_13_p;
	int    temp_13_len =
	  o.Gen_Product_With_PreAlloc(2, yca, y_nor_opq_len, y_nor_opq, &temp_13, 32);
	double temp_14_p[32], *temp_14 = temp_14_p;
	int    temp_14_len =
	  o.Gen_Product_With_PreAlloc(2, zca, z_nor_opq_len, z_nor_opq, &temp_14, 32);
	double temp_15_p[32], *temp_15 = temp_15_p;
	int temp_15_len = o.Gen_Sum_With_PreAlloc(temp_12_len, temp_12, temp_13_len,
	                                          temp_13, &temp_15, 32);
	double ca_dot_nor_opq_p[32], *ca_dot_nor_opq = ca_dot_nor_opq_p;
	int    ca_dot_nor_opq_len = o.Gen_Sum_With_PreAlloc(
    temp_14_len, temp_14, temp_15_len, temp_15, &ca_dot_nor_opq, 32);
	double temp_16_p[32], *temp_16 = temp_16_p;
	int    temp_16_len =
	  o.Gen_Product_With_PreAlloc(2, xba, x_nor_rst_len, x_nor_rst, &temp_16, 32);
	double temp_17_p[32], *temp_17 = temp_17_p;
	int    temp_17_len =
	  o.Gen_Product_With_PreAlloc(2, yba, y_nor_rst_len, y_nor_rst, &temp_17, 32);
	double temp_18_p[32], *temp_18 = temp_18_p;
	int    temp_18_len =
	  o.Gen_Product_With_PreAlloc(2, zba, z_nor_rst_len, z_nor_rst, &temp_18, 32);
	double temp_19_p[32], *temp_19 = temp_19_p;
	int temp_19_len = o.Gen_Sum_With_PreAlloc(temp_16_len, temp_16, temp_17_len,
	                                          temp_17, &temp_19, 32);
	double ba_dot_nor_rst_p[32], *ba_dot_nor_rst = ba_dot_nor_rst_p;
	int    ba_dot_nor_rst_len = o.Gen_Sum_With_PreAlloc(
    temp_18_len, temp_18, temp_19_len, temp_19, &ba_dot_nor_rst, 32);
	double temp_20_p[32], *temp_20 = temp_20_p;
	int    temp_20_len =
	  o.Gen_Product_With_PreAlloc(2, xca, x_nor_rst_len, x_nor_rst, &temp_20, 32);
	double temp_21_p[32], *temp_21 = temp_21_p;
	int    temp_21_len =
	  o.Gen_Product_With_PreAlloc(2, yca, y_nor_rst_len, y_nor_rst, &temp_21, 32);
	double temp_22_p[32], *temp_22 = temp_22_p;
	int    temp_22_len =
	  o.Gen_Product_With_PreAlloc(2, zca, z_nor_rst_len, z_nor_rst, &temp_22, 32);
	double temp_23_p[32], *temp_23 = temp_23_p;
	int temp_23_len = o.Gen_Sum_With_PreAlloc(temp_20_len, temp_20, temp_21_len,
	                                          temp_21, &temp_23, 32);
	double ca_dot_nor_rst_p[32], *ca_dot_nor_rst = ca_dot_nor_rst_p;
	int    ca_dot_nor_rst_len = o.Gen_Sum_With_PreAlloc(
    temp_22_len, temp_22, temp_23_len, temp_23, &ca_dot_nor_rst, 32);
	double temp_100_p[32], *temp_100 = temp_100_p;
	int    temp_100_len = o.Gen_Product_With_PreAlloc(
    ba_dot_nor_opq_len, ba_dot_nor_opq, ca_dot_nor_rst_len, ca_dot_nor_rst,
    &temp_100, 32);
	double temp_101_p[32], *temp_101 = temp_101_p;
	int    temp_101_len = o.Gen_Product_With_PreAlloc(
    ca_dot_nor_opq_len, ca_dot_nor_opq, ba_dot_nor_rst_len, ba_dot_nor_rst,
    &temp_101, 32);
	lambda_d_len = o.Gen_Diff_With_PreAlloc(temp_100_len, temp_100, temp_101_len,
	                                        temp_101, lambda_d, lambda_d_len);
	double temp_102_p[32], *temp_102 = temp_102_p;
	int    temp_102_len = o.Gen_Product_With_PreAlloc(
    oa_dot_nor_opq_len, oa_dot_nor_opq, ca_dot_nor_rst_len, ca_dot_nor_rst,
    &temp_102, 32);
	double temp_103_p[32], *temp_103 = temp_103_p;
	int    temp_103_len = o.Gen_Product_With_PreAlloc(
    ca_dot_nor_opq_len, ca_dot_nor_opq, ra_dot_nor_rst_len, ra_dot_nor_rst,
    &temp_103, 32);
	double det_sub0_p[32], *det_sub0 = det_sub0_p;
	int    det_sub0_len = o.Gen_Diff_With_PreAlloc(
    temp_102_len, temp_102, temp_103_len, temp_103, &det_sub0, 32);
	double temp_105_p[32], *temp_105 = temp_105_p;
	int    temp_105_len = o.Gen_Product_With_PreAlloc(
    ba_dot_nor_opq_len, ba_dot_nor_opq, ra_dot_nor_rst_len, ra_dot_nor_rst,
    &temp_105, 32);
	double temp_106_p[32], *temp_106 = temp_106_p;
	int    temp_106_len = o.Gen_Product_With_PreAlloc(
    oa_dot_nor_opq_len, oa_dot_nor_opq, ba_dot_nor_rst_len, ba_dot_nor_rst,
    &temp_106, 32);
	double det_sub1_p[32], *det_sub1 = det_sub1_p;
	int    det_sub1_len = o.Gen_Diff_With_PreAlloc(
    temp_105_len, temp_105, temp_106_len, temp_106, &det_sub1, 32);
	double xu_p[32], *xu = xu_p;
	int    xu_len =
	  o.Gen_Product_With_PreAlloc(det_sub0_len, det_sub0, 2, xba, &xu, 32);
	double yu_p[32], *yu = yu_p;
	int    yu_len =
	  o.Gen_Product_With_PreAlloc(det_sub0_len, det_sub0, 2, yba, &yu, 32);
	double zu_p[32], *zu = zu_p;
	int    zu_len =
	  o.Gen_Product_With_PreAlloc(det_sub0_len, det_sub0, 2, zba, &zu, 32);
	double xv_p[32], *xv = xv_p;
	int    xv_len =
	  o.Gen_Product_With_PreAlloc(det_sub1_len, det_sub1, 2, xca, &xv, 32);
	double yv_p[32], *yv = yv_p;
	int    yv_len =
	  o.Gen_Product_With_PreAlloc(det_sub1_len, det_sub1, 2, yca, &yv, 32);
	double zv_p[32], *zv = zv_p;
	int    zv_len =
	  o.Gen_Product_With_PreAlloc(det_sub1_len, det_sub1, 2, zca, &zv, 32);
	lambda_x_len =
	  o.Gen_Sum_With_PreAlloc(xu_len, xu, xv_len, xv, lambda_x, lambda_x_len);
	lambda_y_len =
	  o.Gen_Sum_With_PreAlloc(yu_len, yu, yv_len, yv, lambda_y, lambda_y_len);
	lambda_z_len =
	  o.Gen_Sum_With_PreAlloc(zu_len, zu, zv_len, zv, lambda_z, lambda_z_len);
	beta_x = xa;
	beta_y = ya;
	beta_z = za;

	if (zv_p != zv)
		FreeDoubles(zv);
	if (yv_p != yv)
		FreeDoubles(yv);
	if (xv_p != xv)
		FreeDoubles(xv);
	if (zu_p != zu)
		FreeDoubles(zu);
	if (yu_p != yu)
		FreeDoubles(yu);
	if (xu_p != xu)
		FreeDoubles(xu);
	if (det_sub1_p != det_sub1)
		FreeDoubles(det_sub1);
	if (temp_106_p != temp_106)
		FreeDoubles(temp_106);
	if (temp_105_p != temp_105)
		FreeDoubles(temp_105);
	if (det_sub0_p != det_sub0)
		FreeDoubles(det_sub0);
	if (temp_103_p != temp_103)
		FreeDoubles(temp_103);
	if (temp_102_p != temp_102)
		FreeDoubles(temp_102);
	if (temp_101_p != temp_101)
		FreeDoubles(temp_101);
	if (temp_100_p != temp_100)
		FreeDoubles(temp_100);
	if (ca_dot_nor_rst_p != ca_dot_nor_rst)
		FreeDoubles(ca_dot_nor_rst);
	if (temp_23_p != temp_23)
		FreeDoubles(temp_23);
	if (temp_22_p != temp_22)
		FreeDoubles(temp_22);
	if (temp_21_p != temp_21)
		FreeDoubles(temp_21);
	if (temp_20_p != temp_20)
		FreeDoubles(temp_20);
	if (ba_dot_nor_rst_p != ba_dot_nor_rst)
		FreeDoubles(ba_dot_nor_rst);
	if (temp_19_p != temp_19)
		FreeDoubles(temp_19);
	if (temp_18_p != temp_18)
		FreeDoubles(temp_18);
	if (temp_17_p != temp_17)
		FreeDoubles(temp_17);
	if (temp_16_p != temp_16)
		FreeDoubles(temp_16);
	if (ca_dot_nor_opq_p != ca_dot_nor_opq)
		FreeDoubles(ca_dot_nor_opq);
	if (temp_15_p != temp_15)
		FreeDoubles(temp_15);
	if (temp_14_p != temp_14)
		FreeDoubles(temp_14);
	if (temp_13_p != temp_13)
		FreeDoubles(temp_13);
	if (temp_12_p != temp_12)
		FreeDoubles(temp_12);
	if (ba_dot_nor_opq_p != ba_dot_nor_opq)
		FreeDoubles(ba_dot_nor_opq);
	if (temp_11_p != temp_11)
		FreeDoubles(temp_11);
	if (temp_10_p != temp_10)
		FreeDoubles(temp_10);
	if (temp_9_p != temp_9)
		FreeDoubles(temp_9);
	if (temp_8_p != temp_8)
		FreeDoubles(temp_8);
	if (ra_dot_nor_rst_p != ra_dot_nor_rst)
		FreeDoubles(ra_dot_nor_rst);
	if (temp_7_p != temp_7)
		FreeDoubles(temp_7);
	if (temp_6_p != temp_6)
		FreeDoubles(temp_6);
	if (temp_5_p != temp_5)
		FreeDoubles(temp_5);
	if (temp_4_p != temp_4)
		FreeDoubles(temp_4);
	if (oa_dot_nor_opq_p != oa_dot_nor_opq)
		FreeDoubles(oa_dot_nor_opq);
	if (temp_3_p != temp_3)
		FreeDoubles(temp_3);
	if (temp_2_p != temp_2)
		FreeDoubles(temp_2);
	if (temp_1_p != temp_1)
		FreeDoubles(temp_1);
	if (temp_0_p != temp_0)
		FreeDoubles(temp_0);
}

#endif

} // namespace OMC