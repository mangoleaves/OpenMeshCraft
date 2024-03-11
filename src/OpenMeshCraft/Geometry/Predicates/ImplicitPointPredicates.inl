#pragma once
#include "ImplicitPointPredicates.h"

#include "OpenMeshCraft/NumberTypes/ExpansionObject.h"
#include "OpenMeshCraft/NumberTypes/IntervalNumber.h"

namespace OMC {

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

template <typename _IT>
bool lambda2d_SSI_interval(_IT ea1x, _IT ea1y, _IT ea2x, _IT ea2y, _IT eb1x,
                           _IT eb1y, _IT eb2x, _IT eb2y, _IT &lambda_x,
                           _IT &lambda_y, _IT &lambda_det)
{
	typename _IT::Protector P;

	_IT t1a(ea1x * ea2y);
	_IT t1b(ea2x * ea1y);
	_IT t1(t1a - t1b);
	_IT tx2(eb1x - eb2x);
	_IT t3a(eb1x * eb2y);
	_IT t3b(eb2x * eb1y);
	_IT t3(t3a - t3b);
	_IT tx4(ea1x - ea2x);
	_IT ty2(eb1y - eb2y);
	_IT ty4(ea1y - ea2y);
	_IT lxa(t1 * tx2);
	_IT lxb(t3 * tx4);
	lambda_x = lxa - lxb;
	_IT lya(t1 * ty2);
	_IT lyb(t3 * ty4);
	lambda_y = lya - lyb;
	_IT deta(tx4 * ty2);
	_IT detb(tx2 * ty4);
	lambda_det = deta - detb;

	return (lambda_det.is_sign_reliable());
}

template <typename _ET>
void lambda2d_SSI_exact(const _ET &ea1x, const _ET &ea1y, const _ET &ea2x,
                        const _ET &ea2y, const _ET &eb1x, const _ET &eb1y,
                        const _ET &eb2x, const _ET &eb2y, _ET &lambda_x,
                        _ET &lambda_y, _ET &lambda_det)
{
	_ET t1a = ea1x * ea2y;
	_ET t1b = ea2x * ea1y;
	_ET t1  = t1a - t1b;
	_ET tx2 = eb1x - eb2x;
	_ET t3a = eb1x * eb2y;
	_ET t3b = eb2x * eb1y;
	_ET t3  = t3a - t3b;
	_ET tx4 = ea1x - ea2x;
	_ET ty2 = eb1y - eb2y;
	_ET ty4 = ea1y - ea2y;
	_ET lxa = t1 * tx2;
	_ET lxb = t3 * tx4;

	lambda_x   = lxa - lxb;
	_ET lya    = t1 * ty2;
	_ET lyb    = t3 * ty4;
	lambda_y   = lya - lyb;
	_ET deta   = tx4 * ty2;
	_ET detb   = tx2 * ty4;
	lambda_det = deta - detb;
}

template <typename _IT>
bool lambda3d_SSI_interval(_IT xa, _IT ya, _IT za, _IT xb, _IT yb, _IT zb,
                           _IT xp, _IT yp, _IT xq, _IT yq, _IT &lambda_x,
                           _IT &lambda_y, _IT &lambda_z, _IT &lambda_d)
{
	typename _IT::Protector P;

	_IT xba  = xb - xa;
	_IT yba  = yb - ya;
	_IT zba  = zb - za;
	_IT xap  = xa - xp;
	_IT yap  = ya - yp;
	_IT yqp  = yq - yp;
	_IT xqp  = xq - xp;
	_IT c1   = xap * yqp;
	_IT c2   = xqp * yap;
	_IT c3   = xba * yqp;
	_IT c4   = xqp * yba;
	_IT c12  = c1 - c2;
	lambda_d = c4 - c3;
	_IT xf   = xba * c12;
	_IT xs   = xa * lambda_d;
	lambda_x = xf + xs;
	_IT yf   = yba * c12;
	_IT ys   = ya * lambda_d;
	lambda_y = yf + ys;
	_IT zf   = zba * c12;
	_IT zs   = za * lambda_d;
	lambda_z = zf + zs;
	return lambda_d.is_sign_reliable();
}

template <typename _ET>
void lambda3d_SSI_exact(_ET xa, _ET ya, _ET za, _ET xb, _ET yb, _ET zb, _ET xp,
                        _ET yp, _ET xq, _ET yq, _ET &lambda_x, _ET &lambda_y,
                        _ET &lambda_z, _ET &lambda_d)
{
	_ET xba  = xb - xa;
	_ET yba  = yb - ya;
	_ET zba  = zb - za;
	_ET xap  = xa - xp;
	_ET yap  = ya - yp;
	_ET yqp  = yq - yp;
	_ET xqp  = xq - xp;
	_ET c1   = xap * yqp;
	_ET c2   = xqp * yap;
	_ET c3   = xba * yqp;
	_ET c4   = xqp * yba;
	_ET c12  = c1 - c2;
	lambda_d = c4 - c3;
	_ET xf   = xba * c12;
	_ET xs   = xa * lambda_d;
	lambda_x = xf + xs;
	_ET yf   = yba * c12;
	_ET ys   = ya * lambda_d;
	lambda_y = yf + ys;
	_ET zf   = zba * c12;
	_ET zs   = za * lambda_d;
	lambda_z = zf + zs;
}

template <typename _IT>
bool lambda3d_LPI_interval(_IT px, _IT py, _IT pz, _IT qx, _IT qy, _IT qz,
                           _IT rx, _IT ry, _IT rz, _IT sx, _IT sy, _IT sz,
                           _IT tx, _IT ty, _IT tz, _IT &lambda_d, _IT &lambda_x,
                           _IT &lambda_y, _IT &lambda_z)
{
	typename _IT::Protector P;

	_IT a11(px - qx);
	_IT a12(py - qy);
	_IT a13(pz - qz);
	_IT a21(sx - rx);
	_IT a22(sy - ry);
	_IT a23(sz - rz);
	_IT a31(tx - rx);
	_IT a32(ty - ry);
	_IT a33(tz - rz);
	_IT tv1(a22 * a33);
	_IT tv2(a23 * a32);
	_IT a2233(tv1 - tv2);
	_IT tv3(a21 * a33);
	_IT tv4(a23 * a31);
	_IT a2133(tv3 - tv4);
	_IT tv5(a21 * a32);
	_IT tv6(a22 * a31);
	_IT a2132(tv5 - tv6);
	_IT tv7(a11 * a2233);
	_IT tv8(a12 * a2133);
	_IT tv9(a13 * a2132);
	_IT tt1(tv7 - tv8);
	lambda_d = tt1 + tv9;
	_IT px_rx(px - rx);
	_IT py_ry(py - ry);
	_IT pz_rz(pz - rz);
	_IT tt2(py_ry * a2133);
	_IT tt3(px_rx * a2233);
	_IT tt4(pz_rz * a2132);
	_IT tt5(tt3 + tt4);
	_IT n(tt5 - tt2);
	_IT ax(a11 * n);
	_IT ay(a12 * n);
	_IT az(a13 * n);
	_IT dpx(lambda_d * px);
	_IT dpy(lambda_d * py);
	_IT dpz(lambda_d * pz);
	lambda_x = dpx - ax;
	lambda_y = dpy - ay;
	lambda_z = dpz - az;

	return (lambda_d.is_sign_reliable());
}

template <typename _ET>
void lambda3d_LPI_exact(const _ET &px, const _ET &py, const _ET &pz,
                        const _ET &qx, const _ET &qy, const _ET &qz,
                        const _ET &rx, const _ET &ry, const _ET &rz,
                        const _ET &sx, const _ET &sy, const _ET &sz,
                        const _ET &tx, const _ET &ty, const _ET &tz,
                        _ET &lambda_d, _ET &lambda_x, _ET &lambda_y,
                        _ET &lambda_z)
{
	_ET a11   = px - qx;
	_ET a12   = py - qy;
	_ET a13   = pz - qz;
	_ET a21   = sx - rx;
	_ET a22   = sy - ry;
	_ET a23   = sz - rz;
	_ET a31   = tx - rx;
	_ET a32   = ty - ry;
	_ET a33   = tz - rz;
	_ET tv1   = a22 * a33;
	_ET tv2   = a23 * a32;
	_ET a2233 = tv1 - tv2;
	_ET tv3   = a21 * a33;
	_ET tv4   = a23 * a31;
	_ET a2133 = tv3 - tv4;
	_ET tv5   = a21 * a32;
	_ET tv6   = a22 * a31;
	_ET a2132 = tv5 - tv6;
	_ET tv7   = a11 * a2233;
	_ET tv8   = a12 * a2133;
	_ET tv9   = a13 * a2132;
	_ET tt1   = tv7 - tv8;
	lambda_d  = tt1 + tv9;
	_ET px_rx = px - rx;
	_ET py_ry = py - ry;
	_ET pz_rz = pz - rz;
	_ET tt2   = py_ry * a2133;
	_ET tt3   = px_rx * a2233;
	_ET tt4   = pz_rz * a2132;
	_ET tt5   = tt3 + tt4;
	_ET n     = tt5 - tt2;
	_ET ax    = a11 * n;
	_ET ay    = a12 * n;
	_ET az    = a13 * n;
	_ET dpx   = lambda_d * px;
	_ET dpy   = lambda_d * py;
	_ET dpz   = lambda_d * pz;
	lambda_x  = dpx - ax;
	lambda_y  = dpy - ay;
	lambda_z  = dpz - az;
}

template <typename _IT>
bool lambda3d_TPI_interval(_IT ov1x, _IT ov1y, _IT ov1z, _IT ov2x, _IT ov2y,
                           _IT ov2z, _IT ov3x, _IT ov3y, _IT ov3z, _IT ow1x,
                           _IT ow1y, _IT ow1z, _IT ow2x, _IT ow2y, _IT ow2z,
                           _IT ow3x, _IT ow3y, _IT ow3z, _IT ou1x, _IT ou1y,
                           _IT ou1z, _IT ou2x, _IT ou2y, _IT ou2z, _IT ou3x,
                           _IT ou3y, _IT ou3z, _IT &lambda_x, _IT &lambda_y,
                           _IT &lambda_z, _IT &lambda_d)
{
	typename _IT::Protector P;

	_IT v3x(ov3x - ov2x);
	_IT v3y(ov3y - ov2y);
	_IT v3z(ov3z - ov2z);
	_IT v2x(ov2x - ov1x);
	_IT v2y(ov2y - ov1y);
	_IT v2z(ov2z - ov1z);
	_IT w3x(ow3x - ow2x);
	_IT w3y(ow3y - ow2y);
	_IT w3z(ow3z - ow2z);
	_IT w2x(ow2x - ow1x);
	_IT w2y(ow2y - ow1y);
	_IT w2z(ow2z - ow1z);
	_IT u3x(ou3x - ou2x);
	_IT u3y(ou3y - ou2y);
	_IT u3z(ou3z - ou2z);
	_IT u2x(ou2x - ou1x);
	_IT u2y(ou2y - ou1y);
	_IT u2z(ou2z - ou1z);
	_IT nvx1(v2y * v3z);
	_IT nvx2(v2z * v3y);
	_IT nvx(nvx1 - nvx2);
	_IT nvy1(v3x * v2z);
	_IT nvy2(v3z * v2x);
	_IT nvy(nvy1 - nvy2);
	_IT nvz1(v2x * v3y);
	_IT nvz2(v2y * v3x);
	_IT nvz(nvz1 - nvz2);
	_IT nwx1(w2y * w3z);
	_IT nwx2(w2z * w3y);
	_IT nwx(nwx1 - nwx2);
	_IT nwy1(w3x * w2z);
	_IT nwy2(w3z * w2x);
	_IT nwy(nwy1 - nwy2);
	_IT nwz1(w2x * w3y);
	_IT nwz2(w2y * w3x);
	_IT nwz(nwz1 - nwz2);
	_IT nux1(u2y * u3z);
	_IT nux2(u2z * u3y);
	_IT nux(nux1 - nux2);
	_IT nuy1(u3x * u2z);
	_IT nuy2(u3z * u2x);
	_IT nuy(nuy1 - nuy2);
	_IT nuz1(u2x * u3y);
	_IT nuz2(u2y * u3x);
	_IT nuz(nuz1 - nuz2);
	_IT nwyuz1(nwy * nuz);
	_IT nwyuz2(nwz * nuy);
	_IT nwyuz(nwyuz1 - nwyuz2);
	_IT nwxuz1(nwx * nuz);
	_IT nwxuz2(nwz * nux);
	_IT nwxuz(nwxuz1 - nwxuz2);
	_IT nwxuy1(nwx * nuy);
	_IT nwxuy2(nwy * nux);
	_IT nwxuy(nwxuy1 - nwxuy2);
	_IT nvyuz1(nvy * nuz);
	_IT nvyuz2(nvz * nuy);
	_IT nvyuz(nvyuz1 - nvyuz2);
	_IT nvxuz1(nvx * nuz);
	_IT nvxuz2(nvz * nux);
	_IT nvxuz(nvxuz1 - nvxuz2);
	_IT nvxuy1(nvx * nuy);
	_IT nvxuy2(nvy * nux);
	_IT nvxuy(nvxuy1 - nvxuy2);
	_IT nvywz1(nvy * nwz);
	_IT nvywz2(nvz * nwy);
	_IT nvywz(nvywz1 - nvywz2);
	_IT nvxwz1(nvx * nwz);
	_IT nvxwz2(nvz * nwx);
	_IT nvxwz(nvxwz1 - nvxwz2);
	_IT nvxwy1(nvx * nwy);
	_IT nvxwy2(nvy * nwx);
	_IT nvxwy(nvxwy1 - nvxwy2);
	_IT p1a(nvx * ov1x);
	_IT p1b(nvy * ov1y);
	_IT p1c(nvz * ov1z);
	_IT p1ab(p1a + p1b);
	_IT p1(p1ab + p1c);
	_IT p2a(nwx * ow1x);
	_IT p2b(nwy * ow1y);
	_IT p2c(nwz * ow1z);
	_IT p2ab(p2a + p2b);
	_IT p2(p2ab + p2c);
	_IT p3a(nux * ou1x);
	_IT p3b(nuy * ou1y);
	_IT p3c(nuz * ou1z);
	_IT p3ab(p3a + p3b);
	_IT p3(p3ab + p3c);
	_IT lxa(p1 * nwyuz);
	_IT lxb(p3 * nvywz);
	_IT lxc(p2 * nvyuz);
	_IT lxab(lxa + lxb);
	lambda_x = lxab - lxc;
	_IT lya(p2 * nvxuz);
	_IT lyb(p3 * nvxwz);
	_IT lyc(p1 * nwxuz);
	_IT lybc(lyc + lyb);
	lambda_y = lya - lybc;
	_IT lza(p3 * nvxwy);
	_IT lzb(p1 * nwxuy);
	_IT lzc(p2 * nvxuy);
	_IT lzab(lza + lzb);
	lambda_z = lzab - lzc;
	_IT da(nvx * nwyuz);
	_IT db(nvz * nwxuy);
	_IT dc(nvy * nwxuz);
	_IT dab(da + db);
	lambda_d = dab - dc;

	return (lambda_d.is_sign_reliable());
}

template <typename _ET>
void lambda3d_TPI_exact(const _ET &ov1x, const _ET &ov1y, const _ET &ov1z,
                        const _ET &ov2x, const _ET &ov2y, const _ET &ov2z,
                        const _ET &ov3x, const _ET &ov3y, const _ET &ov3z,
                        const _ET &ow1x, const _ET &ow1y, const _ET &ow1z,
                        const _ET &ow2x, const _ET &ow2y, const _ET &ow2z,
                        const _ET &ow3x, const _ET &ow3y, const _ET &ow3z,
                        const _ET &ou1x, const _ET &ou1y, const _ET &ou1z,
                        const _ET &ou2x, const _ET &ou2y, const _ET &ou2z,
                        const _ET &ou3x, const _ET &ou3y, const _ET &ou3z,
                        _ET &lambda_x, _ET &lambda_y, _ET &lambda_z,
                        _ET &lambda_d)
{
	_ET v3x    = ov3x - ov2x;
	_ET v3y    = ov3y - ov2y;
	_ET v3z    = ov3z - ov2z;
	_ET v2x    = ov2x - ov1x;
	_ET v2y    = ov2y - ov1y;
	_ET v2z    = ov2z - ov1z;
	_ET w3x    = ow3x - ow2x;
	_ET w3y    = ow3y - ow2y;
	_ET w3z    = ow3z - ow2z;
	_ET w2x    = ow2x - ow1x;
	_ET w2y    = ow2y - ow1y;
	_ET w2z    = ow2z - ow1z;
	_ET u3x    = ou3x - ou2x;
	_ET u3y    = ou3y - ou2y;
	_ET u3z    = ou3z - ou2z;
	_ET u2x    = ou2x - ou1x;
	_ET u2y    = ou2y - ou1y;
	_ET u2z    = ou2z - ou1z;
	_ET nvx1   = v2y * v3z;
	_ET nvx2   = v2z * v3y;
	_ET nvx    = nvx1 - nvx2;
	_ET nvy1   = v3x * v2z;
	_ET nvy2   = v3z * v2x;
	_ET nvy    = nvy1 - nvy2;
	_ET nvz1   = v2x * v3y;
	_ET nvz2   = v2y * v3x;
	_ET nvz    = nvz1 - nvz2;
	_ET nwx1   = w2y * w3z;
	_ET nwx2   = w2z * w3y;
	_ET nwx    = nwx1 - nwx2;
	_ET nwy1   = w3x * w2z;
	_ET nwy2   = w3z * w2x;
	_ET nwy    = nwy1 - nwy2;
	_ET nwz1   = w2x * w3y;
	_ET nwz2   = w2y * w3x;
	_ET nwz    = nwz1 - nwz2;
	_ET nux1   = u2y * u3z;
	_ET nux2   = u2z * u3y;
	_ET nux    = nux1 - nux2;
	_ET nuy1   = u3x * u2z;
	_ET nuy2   = u3z * u2x;
	_ET nuy    = nuy1 - nuy2;
	_ET nuz1   = u2x * u3y;
	_ET nuz2   = u2y * u3x;
	_ET nuz    = nuz1 - nuz2;
	_ET nwyuz1 = nwy * nuz;
	_ET nwyuz2 = nwz * nuy;
	_ET nwyuz  = nwyuz1 - nwyuz2;
	_ET nwxuz1 = nwx * nuz;
	_ET nwxuz2 = nwz * nux;
	_ET nwxuz  = nwxuz1 - nwxuz2;
	_ET nwxuy1 = nwx * nuy;
	_ET nwxuy2 = nwy * nux;
	_ET nwxuy  = nwxuy1 - nwxuy2;
	_ET nvyuz1 = nvy * nuz;
	_ET nvyuz2 = nvz * nuy;
	_ET nvyuz  = nvyuz1 - nvyuz2;
	_ET nvxuz1 = nvx * nuz;
	_ET nvxuz2 = nvz * nux;
	_ET nvxuz  = nvxuz1 - nvxuz2;
	_ET nvxuy1 = nvx * nuy;
	_ET nvxuy2 = nvy * nux;
	_ET nvxuy  = nvxuy1 - nvxuy2;
	_ET nvywz1 = nvy * nwz;
	_ET nvywz2 = nvz * nwy;
	_ET nvywz  = nvywz1 - nvywz2;
	_ET nvxwz1 = nvx * nwz;
	_ET nvxwz2 = nvz * nwx;
	_ET nvxwz  = nvxwz1 - nvxwz2;
	_ET nvxwy1 = nvx * nwy;
	_ET nvxwy2 = nvy * nwx;
	_ET nvxwy  = nvxwy1 - nvxwy2;
	_ET p1a    = nvx * ov1x;
	_ET p1b    = nvy * ov1y;
	_ET p1c    = nvz * ov1z;
	_ET p1ab   = p1a + p1b;
	_ET p1     = p1ab + p1c;
	_ET p2a    = nwx * ow1x;
	_ET p2b    = nwy * ow1y;
	_ET p2c    = nwz * ow1z;
	_ET p2ab   = p2a + p2b;
	_ET p2     = p2ab + p2c;
	_ET p3a    = nux * ou1x;
	_ET p3b    = nuy * ou1y;
	_ET p3c    = nuz * ou1z;
	_ET p3ab   = p3a + p3b;
	_ET p3     = p3ab + p3c;
	_ET lxa    = p1 * nwyuz;
	_ET lxb    = p3 * nvywz;
	_ET lxc    = p2 * nvyuz;
	_ET lxab   = lxa + lxb;
	lambda_x   = lxab - lxc;
	_ET lya    = p2 * nvxuz;
	_ET lyb    = p3 * nvxwz;
	_ET lyc    = p1 * nwxuz;
	_ET lybc   = lyc + lyb;
	lambda_y   = lya - lybc;
	_ET lza    = p3 * nvxwy;
	_ET lzb    = p1 * nwxuy;
	_ET lzc    = p2 * nvxuy;
	_ET lzab   = lza + lzb;
	lambda_z   = lzab - lzc;
	_ET da     = nvx * nwyuz;
	_ET db     = nvz * nwxuy;
	_ET dc     = nvy * nwxuz;
	_ET dab    = da + db;
	lambda_d   = dab - dc;
}

template <typename _IT>
bool lambda3d_LNC_interval(_IT px, _IT py, _IT pz, _IT qx, _IT qy, _IT qz,
                           _IT t, _IT &lambda_x, _IT &lambda_y, _IT &lambda_z,
                           _IT &lambda_d)
{
	typename _IT::Protector P;

	_IT vx(px - qx);
	_IT vy(py - qy);
	_IT vz(pz - qz);
	_IT vxt(vx * t);
	_IT vyt(vy * t);
	_IT vzt(vz * t);
	lambda_x = px - vxt;
	lambda_y = py - vyt;
	lambda_z = pz - vzt;
	lambda_d = 1;

	return true;
}

template <typename _ET>
void lambda3d_LNC_exact(_ET px, _ET py, _ET pz, _ET qx, _ET qy, _ET qz, _ET t,
                        _ET &lambda_x, _ET &lambda_y, _ET &lambda_z,
                        _ET &lambda_d)
{
	_ET vx   = px - qx;
	_ET vy   = py - qy;
	_ET vz   = pz - qz;
	_ET vxt  = vx * t;
	_ET vyt  = vy * t;
	_ET vzt  = vz * t;
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

} // namespace OMC