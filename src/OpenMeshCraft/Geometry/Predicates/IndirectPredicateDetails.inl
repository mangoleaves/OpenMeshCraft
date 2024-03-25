#pragma once

#include "IndirectPredicateDetailsHand.h"

#include "OpenMeshCraft/NumberTypes/ExpansionObject.h"
#include "OpenMeshCraft/NumberTypes/IntervalNumber.h"
#include "OpenMeshCraft/NumberTypes/LazyNumber.h"

#pragma intrinsic(fabs)

// Uncomment the following to activate overflow/underflow checks
#define CHECK_FOR_XYZERFLOWS

namespace OMC {

inline Sign dotProductSign2D_filtered(double px, double py, double rx,
                                      double ry, double qx, double qy)
{
	double lx = px - qx;
	double ly = py - qy;
	double gx = rx - qx;
	double gy = ry - qy;
	double dx = lx * gx;
	double dy = ly * gy;
	double d  = dx + dy;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gy)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 8.881784197001252e-16;

	return filter_sign(d, epsilon);
}

template <typename IT>
Sign dotProductSign2D_interval(IT px, IT py, IT rx, IT ry, IT qx, IT qy)
{
	typename IT::Protector P;

	IT lx = px - qx;
	IT ly = py - qy;
	IT gx = rx - qx;
	IT gy = ry - qy;
	IT dx = lx * gx;
	IT dy = ly * gy;
	IT d  = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename ET>
Sign dotProductSign2D_exact(ET px, ET py, ET rx, ET ry, ET qx, ET qy)
{
	ET lx = px - qx;
	ET ly = py - qy;
	ET gx = rx - qx;
	ET gy = ry - qy;
	ET dx = lx * gx;
	ET dy = ly * gy;
	ET d  = dx + dy;
	return OMC::sign(d);
}

Sign dotProductSign2D_expansion(double px, double py, double rx, double ry,
                                double qx, double qy)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double ly[2];
	o.two_Diff(py, qy, ly);
	double gx[2];
	o.two_Diff(rx, qx, gx);
	double gy[2];
	o.two_Diff(ry, qy, gy);
	double dx[8];
	int    dx_len = o.Gen_Product(2, lx, 2, gx, dx);
	double dy[8];
	int    dy_len = o.Gen_Product(2, ly, 2, gy, dy);
	double d[16];
	int    d_len = o.Gen_Sum(dx_len, dx, dy_len, dy, d);

	double return_value = d[d_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D(double px, double py, double rx, double ry, double qx,
                      double qy)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = dotProductSign2D_filtered(px, py, rx, ry, qx, qy);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = dotProductSign2D_interval<IT>(px, py, rx, ry, qx, qy);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign2D_expansion(px, py, rx, ry, qx, qy);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D(const GenericPoint2T<IT, ET> &p,
                      const GenericPoint2T<IT, ET> &r,
                      const GenericPoint2T<IT, ET> &q)
{
	return dotProductSign2D<IT, ET, WithSSFilter>(p.x(), p.y(), r.x(), r.y(),
	                                              q.x(), q.y());
}

inline Sign dotProductSign2D4P_filtered(double px, double py, double rx,
                                        double ry, double qx, double qy,
                                        double sx, double sy)
{
	double lx = px - qx;
	double ly = py - qy;
	double gx = rx - sx;
	double gy = ry - sy;
	double dx = lx * gx;
	double dy = ly * gy;
	double d  = dx + dy;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gy)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 8.881784197001252e-16;

	return filter_sign(d, epsilon);
}

template <typename IT>
Sign dotProductSign2D4P_interval(IT px, IT py, IT rx, IT ry, IT qx, IT qy,
                                 IT sx, IT sy)
{
	typename IT::Protector P;

	IT lx = px - qx;
	IT ly = py - qy;
	IT gx = rx - sx;
	IT gy = ry - sy;
	IT dx = lx * gx;
	IT dy = ly * gy;
	IT d  = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename ET>
Sign dotProductSign2D4P_exact(ET px, ET py, ET rx, ET ry, ET qx, ET qy, ET sx,
                              ET sy)
{
	ET lx = px - qx;
	ET ly = py - qy;
	ET gx = rx - sx;
	ET gy = ry - sy;
	ET dx = lx * gx;
	ET dy = ly * gy;
	ET d  = dx + dy;
	return OMC::sign(d);
}

Sign dotProductSign2D4P_expansion(double px, double py, double rx, double ry,
                                  double qx, double qy, double sx, double sy)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double ly[2];
	o.two_Diff(py, qy, ly);
	double gx[2];
	o.two_Diff(rx, sx, gx);
	double gy[2];
	o.two_Diff(ry, sy, gy);
	double dx[8];
	int    dx_len = o.Gen_Product(2, lx, 2, gx, dx);
	double dy[8];
	int    dy_len = o.Gen_Product(2, ly, 2, gy, dy);
	double d[16];
	int    d_len = o.Gen_Sum(dx_len, dx, dy_len, dy, d);

	double return_value = d[d_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D4P(double px, double py, double rx, double ry, double qx,
                        double qy, double sx, double sy)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = dotProductSign2D4P_filtered(px, py, rx, ry, qx, qy, sx, sy);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = dotProductSign2D4P_interval<IT>(px, py, rx, ry, qx, qy, sx, sy);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign2D4P_expansion(px, py, rx, ry, qx, qy, sx, sy);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D4P(const GenericPoint2T<IT, ET> &p,
                        const GenericPoint2T<IT, ET> &r,
                        const GenericPoint2T<IT, ET> &q,
                        const GenericPoint2T<IT, ET> &s)
{
	return dotProductSign2D4P<IT, ET, WithSSFilter>(p.x(), p.y(), r.x(), r.y(),
	                                                q.x(), q.y(), s.x(), s.y());
}

inline Sign dotProductSign3D_filtered(double px, double py, double pz,
                                      double rx, double ry, double rz,
                                      double qx, double qy, double qz)
{
	double lx = px - qx;
	double ly = py - qy;
	double lz = pz - qz;
	double gx = rx - qx;
	double gy = ry - qy;
	double gz = rz - qz;
	double dx = lx * gx;
	double dy = ly * gy;
	double dz = lz * gz;
	double d1 = dx + dy;
	double d  = d1 + dz;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(lz)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gy)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gz)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 1.4432899320127035e-15;

	return filter_sign(d, epsilon);
}

template <typename IT>
Sign dotProductSign3D_interval(IT px, IT py, IT pz, IT rx, IT ry, IT rz, IT qx,
                               IT qy, IT qz)
{
	typename IT::Protector P;

	IT lx = px - qx;
	IT ly = py - qy;
	IT lz = pz - qz;
	IT gx = rx - qx;
	IT gy = ry - qy;
	IT gz = rz - qz;
	IT dx = lx * gx;
	IT dy = ly * gy;
	IT dz = lz * gz;
	IT d1 = dx + dy;
	IT d  = d1 + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename ET>
Sign dotProductSign3D_exact(ET px, ET py, ET pz, ET rx, ET ry, ET rz, ET qx,
                            ET qy, ET qz)
{
	ET lx = px - qx;
	ET ly = py - qy;
	ET lz = pz - qz;
	ET gx = rx - qx;
	ET gy = ry - qy;
	ET gz = rz - qz;
	ET dx = lx * gx;
	ET dy = ly * gy;
	ET dz = lz * gz;
	ET d1 = dx + dy;
	ET d  = d1 + dz;
	return OMC::sign(d);
}

Sign dotProductSign3D_expansion(double px, double py, double pz, double rx,
                                double ry, double rz, double qx, double qy,
                                double qz)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double ly[2];
	o.two_Diff(py, qy, ly);
	double lz[2];
	o.two_Diff(pz, qz, lz);
	double gx[2];
	o.two_Diff(rx, qx, gx);
	double gy[2];
	o.two_Diff(ry, qy, gy);
	double gz[2];
	o.two_Diff(rz, qz, gz);
	double dx[8];
	int    dx_len = o.Gen_Product(2, lx, 2, gx, dx);
	double dy[8];
	int    dy_len = o.Gen_Product(2, ly, 2, gy, dy);
	double dz[8];
	int    dz_len = o.Gen_Product(2, lz, 2, gz, dz);
	double d1[16];
	int    d1_len = o.Gen_Sum(dx_len, dx, dy_len, dy, d1);
	double d[24];
	int    d_len = o.Gen_Sum(d1_len, d1, dz_len, dz, d);

	double return_value = d[d_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D(double px, double py, double pz, double rx, double ry,
                      double rz, double qx, double qy, double qz)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = dotProductSign3D_filtered(px, py, pz, rx, ry, rz, qx, qy, qz);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = dotProductSign3D_interval<IT>(px, py, pz, rx, ry, rz, qx, qy, qz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign3D_expansion(px, py, pz, rx, ry, rz, qx, qy, qz);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D(const GenericPoint3T<IT, ET> &p,
                      const GenericPoint3T<IT, ET> &r,
                      const GenericPoint3T<IT, ET> &q)
{
	return dotProductSign3D<IT, ET, WithSSFilter>(
	  p.x(), p.y(), p.z(), r.x(), r.y(), r.z(), q.x(), q.y(), q.z());
}

inline Sign dotProductSign3D4P_filtered(double px, double py, double pz,
                                        double rx, double ry, double rz,
                                        double qx, double qy, double qz,
                                        double sx, double sy, double sz)
{
	double lx = px - qx;
	double ly = py - qy;
	double lz = pz - qz;
	double gx = rx - sx;
	double gy = ry - sy;
	double gz = rz - sz;
	double dx = lx * gx;
	double dy = ly * gy;
	double dz = lz * gz;
	double d1 = dx + dy;
	double d  = d1 + dz;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(lz)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gy)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gz)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 1.4432899320127035e-15;

	return filter_sign(d, epsilon);
}

template <typename IT>
Sign dotProductSign3D4P_interval(IT px, IT py, IT pz, IT rx, IT ry, IT rz,
                                 IT qx, IT qy, IT qz, IT sx, IT sy, IT sz)
{
	typename IT::Protector P;

	IT lx = px - qx;
	IT ly = py - qy;
	IT lz = pz - qz;
	IT gx = rx - sx;
	IT gy = ry - sy;
	IT gz = rz - sz;
	IT dx = lx * gx;
	IT dy = ly * gy;
	IT dz = lz * gz;
	IT d1 = dx + dy;
	IT d  = d1 + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename ET>
Sign dotProductSign3D4P_exact(ET px, ET py, ET pz, ET rx, ET ry, ET rz, ET qx,
                              ET qy, ET qz, ET sx, ET sy, ET sz)
{
	ET lx = px - qx;
	ET ly = py - qy;
	ET lz = pz - qz;
	ET gx = rx - sx;
	ET gy = ry - sy;
	ET gz = rz - sz;
	ET dx = lx * gx;
	ET dy = ly * gy;
	ET dz = lz * gz;
	ET d1 = dx + dy;
	ET d  = d1 + dz;
	return OMC::sign(d);
}

Sign dotProductSign3D4P_expansion(double px, double py, double pz, double rx,
                                  double ry, double rz, double qx, double qy,
                                  double qz, double sx, double sy, double sz)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double ly[2];
	o.two_Diff(py, qy, ly);
	double lz[2];
	o.two_Diff(pz, qz, lz);
	double gx[2];
	o.two_Diff(rx, sx, gx);
	double gy[2];
	o.two_Diff(ry, sy, gy);
	double gz[2];
	o.two_Diff(rz, sz, gz);
	double dx[8];
	int    dx_len = o.Gen_Product(2, lx, 2, gx, dx);
	double dy[8];
	int    dy_len = o.Gen_Product(2, ly, 2, gy, dy);
	double dz[8];
	int    dz_len = o.Gen_Product(2, lz, 2, gz, dz);
	double d1[16];
	int    d1_len = o.Gen_Sum(dx_len, dx, dy_len, dy, d1);
	double d[24];
	int    d_len = o.Gen_Sum(d1_len, d1, dz_len, dz, d);

	double return_value = d[d_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D4P(double px, double py, double pz, double rx, double ry,
                        double rz, double qx, double qy, double qz, double sx,
                        double sy, double sz)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = dotProductSign3D4P_filtered(px, py, pz, rx, ry, rz, qx, qy, qz, sx,
		                                  sy, sz);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = dotProductSign3D4P_interval<IT>(px, py, pz, rx, ry, rz, qx, qy, qz, sx,
	                                      sy, sz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign3D4P_expansion(px, py, pz, rx, ry, rz, qx, qy, qz, sx,
	                                    sy, sz);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D4P(const GenericPoint3T<IT, ET> &p,
                        const GenericPoint3T<IT, ET> &r,
                        const GenericPoint3T<IT, ET> &q,
                        const GenericPoint3T<IT, ET> &s)
{
	return dotProductSign3D4P<IT, ET, WithSSFilter>(p.x(), p.y(), p.z(), r.x(),
	                                                r.y(), r.z(), q.x(), q.y(),
	                                                q.z(), s.x(), s.y(), s.z());
}

inline Sign dotProductSignOn2Dxy4P_filtered(double px, double py, double rx,
                                            double ry, double qx, double qy,
                                            double sx, double sy)
{
	double lx = px - qx;
	double ly = py - qy;
	double gx = rx - sx;
	double gy = ry - sy;
	double dx = lx * gx;
	double dy = ly * gy;
	double d  = dx + dy;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gy)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 8.881784197001252e-16;

	return filter_sign(d, epsilon);
}

template <typename IT>
Sign dotProductSignOn2Dxy4P_interval(IT px, IT py, IT rx, IT ry, IT qx, IT qy,
                                     IT sx, IT sy)
{
	typename IT::Protector P;

	IT lx = px - qx;
	IT ly = py - qy;
	IT gx = rx - sx;
	IT gy = ry - sy;
	IT dx = lx * gx;
	IT dy = ly * gy;
	IT d  = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename ET>
Sign dotProductSignOn2Dxy4P_exact(ET px, ET py, ET rx, ET ry, ET qx, ET qy,
                                  ET sx, ET sy)
{
	ET lx = px - qx;
	ET ly = py - qy;
	ET gx = rx - sx;
	ET gy = ry - sy;
	ET dx = lx * gx;
	ET dy = ly * gy;
	ET d  = dx + dy;
	return OMC::sign(d);
}

Sign dotProductSignOn2Dxy4P_expansion(double px, double py, double rx,
                                      double ry, double qx, double qy,
                                      double sx, double sy)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double ly[2];
	o.two_Diff(py, qy, ly);
	double gx[2];
	o.two_Diff(rx, sx, gx);
	double gy[2];
	o.two_Diff(ry, sy, gy);
	double dx[8];
	int    dx_len = o.Gen_Product(2, lx, 2, gx, dx);
	double dy[8];
	int    dy_len = o.Gen_Product(2, ly, 2, gy, dy);
	double d[16];
	int    d_len = o.Gen_Sum(dx_len, dx, dy_len, dy, d);

	double return_value = d[d_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dxy4P(double px, double py, double rx, double ry,
                            double qx, double qy, double sx, double sy)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = dotProductSignOn2Dxy4P_filtered(px, py, rx, ry, qx, qy, sx, sy);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = dotProductSignOn2Dxy4P_interval<IT>(px, py, rx, ry, qx, qy, sx, sy);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSignOn2Dxy4P_expansion(px, py, rx, ry, qx, qy, sx, sy);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dxy4P(const GenericPoint3T<IT, ET> &p,
                            const GenericPoint3T<IT, ET> &r,
                            const GenericPoint3T<IT, ET> &q,
                            const GenericPoint3T<IT, ET> &s)
{
	return dotProductSignOn2Dxy4P<IT, ET, WithSSFilter>(
	  p.x(), p.y(), r.x(), r.y(), q.x(), q.y(), s.x(), s.y());
}

inline Sign dotProductSignOn2Dyz4P_filtered(double py, double pz, double ry,
                                            double rz, double qy, double qz,
                                            double sy, double sz)
{
	double ly = py - qy;
	double lz = pz - qz;
	double gy = ry - sy;
	double gz = rz - sz;
	double dy = ly * gy;
	double dz = lz * gz;
	double d  = dy + dz;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(lz)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gy)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gz)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 8.881784197001252e-16;

	return filter_sign(d, epsilon);
}

template <typename IT>
Sign dotProductSignOn2Dyz4P_interval(IT py, IT pz, IT ry, IT rz, IT qy, IT qz,
                                     IT sy, IT sz)
{
	typename IT::Protector P;

	IT ly = py - qy;
	IT lz = pz - qz;
	IT gy = ry - sy;
	IT gz = rz - sz;
	IT dy = ly * gy;
	IT dz = lz * gz;
	IT d  = dy + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename ET>
Sign dotProductSignOn2Dyz4P_exact(ET py, ET pz, ET ry, ET rz, ET qy, ET qz,
                                  ET sy, ET sz)
{
	ET ly = py - qy;
	ET lz = pz - qz;
	ET gy = ry - sy;
	ET gz = rz - sz;
	ET dy = ly * gy;
	ET dz = lz * gz;
	ET d  = dy + dz;
	return OMC::sign(d);
}

Sign dotProductSignOn2Dyz4P_expansion(double py, double pz, double ry,
                                      double rz, double qy, double qz,
                                      double sy, double sz)
{
	expansionObject o;
	double          ly[2];
	o.two_Diff(py, qy, ly);
	double lz[2];
	o.two_Diff(pz, qz, lz);
	double gy[2];
	o.two_Diff(ry, sy, gy);
	double gz[2];
	o.two_Diff(rz, sz, gz);
	double dy[8];
	int    dy_len = o.Gen_Product(2, ly, 2, gy, dy);
	double dz[8];
	int    dz_len = o.Gen_Product(2, lz, 2, gz, dz);
	double d[16];
	int    d_len = o.Gen_Sum(dy_len, dy, dz_len, dz, d);

	double return_value = d[d_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dyz4P(double py, double pz, double ry, double rz,
                            double qy, double qz, double sy, double sz)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = dotProductSignOn2Dyz4P_filtered(py, pz, ry, rz, qy, qz, sy, sz);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = dotProductSignOn2Dyz4P_interval<IT>(py, pz, ry, rz, qy, qz, sy, sz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSignOn2Dyz4P_expansion(py, pz, ry, rz, qy, qz, sy, sz);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dyz4P(const GenericPoint3T<IT, ET> &p,
                            const GenericPoint3T<IT, ET> &r,
                            const GenericPoint3T<IT, ET> &q,
                            const GenericPoint3T<IT, ET> &s)
{
	return dotProductSignOn2Dyz4P<IT, ET, WithSSFilter>(
	  p.y(), p.z(), r.y(), r.z(), q.y(), q.z(), s.y(), s.z());
}

inline Sign dotProductSignOn2Dzx4P_filtered(double px, double pz, double rx,
                                            double rz, double qx, double qz,
                                            double sx, double sz)
{
	double lx = px - qx;
	double lz = pz - qz;
	double gx = rx - sx;
	double gz = rz - sz;
	double dx = lx * gx;
	double dz = lz * gz;
	double d  = dx + dz;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(lz)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(gz)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 8.881784197001252e-16;

	return filter_sign(d, epsilon);
}

template <typename IT>
Sign dotProductSignOn2Dzx4P_interval(IT px, IT pz, IT rx, IT rz, IT qx, IT qz,
                                     IT sx, IT sz)
{
	typename IT::Protector P;

	IT lx = px - qx;
	IT lz = pz - qz;
	IT gx = rx - sx;
	IT gz = rz - sz;
	IT dx = lx * gx;
	IT dz = lz * gz;
	IT d  = dx + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename ET>
Sign dotProductSignOn2Dzx4P_exact(ET px, ET pz, ET rx, ET rz, ET qx, ET qz,
                                  ET sx, ET sz)
{
	ET lx = px - qx;
	ET lz = pz - qz;
	ET gx = rx - sx;
	ET gz = rz - sz;
	ET dx = lx * gx;
	ET dz = lz * gz;
	ET d  = dx + dz;
	return OMC::sign(d);
}

Sign dotProductSignOn2Dzx4P_expansion(double px, double pz, double rx,
                                      double rz, double qx, double qz,
                                      double sx, double sz)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double lz[2];
	o.two_Diff(pz, qz, lz);
	double gx[2];
	o.two_Diff(rx, sx, gx);
	double gz[2];
	o.two_Diff(rz, sz, gz);
	double dx[8];
	int    dx_len = o.Gen_Product(2, lx, 2, gx, dx);
	double dz[8];
	int    dz_len = o.Gen_Product(2, lz, 2, gz, dz);
	double d[16];
	int    d_len = o.Gen_Sum(dx_len, dx, dz_len, dz, d);

	double return_value = d[d_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dzx4P(double px, double pz, double rx, double rz,
                            double qx, double qz, double sx, double sz)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = dotProductSignOn2Dzx4P_filtered(px, pz, rx, rz, qx, qz, sx, sz);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = dotProductSignOn2Dzx4P_interval<IT>(px, pz, rx, rz, qx, qz, sx, sz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSignOn2Dzx4P_expansion(px, pz, rx, rz, qx, qz, sx, sz);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dzx4P(const GenericPoint3T<IT, ET> &p,
                            const GenericPoint3T<IT, ET> &r,
                            const GenericPoint3T<IT, ET> &q,
                            const GenericPoint3T<IT, ET> &s)
{
	return dotProductSignOn2Dzx4P<IT, ET, WithSSFilter>(
	  p.x(), p.z(), r.x(), r.z(), q.x(), q.z(), s.x(), s.z());
}

inline Sign inCircle_filtered(double pax, double pay, double pbx, double pby,
                              double pcx, double pcy, double pdx, double pdy)
{
	double adx    = pax - pdx;
	double ady    = pay - pdy;
	double bdx    = pbx - pdx;
	double bdy    = pby - pdy;
	double cdx    = pcx - pdx;
	double cdy    = pcy - pdy;
	double abdeta = adx * bdy;
	double abdetb = bdx * ady;
	double abdet  = abdeta - abdetb;
	double bcdeta = bdx * cdy;
	double bcdetb = cdx * bdy;
	double bcdet  = bcdeta - bcdetb;
	double cadeta = cdx * ady;
	double cadetb = adx * cdy;
	double cadet  = cadeta - cadetb;
	double alifta = adx * adx;
	double aliftb = ady * ady;
	double alift  = alifta + aliftb;
	double blifta = bdx * bdx;
	double bliftb = bdy * bdy;
	double blift  = blifta + bliftb;
	double clifta = cdx * cdx;
	double cliftb = cdy * cdy;
	double clift  = clifta + cliftb;
	double la     = alift * bcdet;
	double lb     = blift * cadet;
	double lc     = clift * abdet;
	double lab    = la + lb;
	double L      = lab + lc;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(adx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ady)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(bdx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(bdy)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(cdx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(cdy)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= epsilon;
	epsilon *= 1.3766765505351944e-14;

	return filter_sign(L, epsilon);
}

template <typename IT>
Sign inCircle_interval(IT pax, IT pay, IT pbx, IT pby, IT pcx, IT pcy, IT pdx,
                       IT pdy)
{
	typename IT::Protector P;

	IT adx    = pax - pdx;
	IT ady    = pay - pdy;
	IT bdx    = pbx - pdx;
	IT bdy    = pby - pdy;
	IT cdx    = pcx - pdx;
	IT cdy    = pcy - pdy;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT alift  = alifta + aliftb;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT blift  = blifta + bliftb;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT clift  = clifta + cliftb;
	IT la     = alift * bcdet;
	IT lb     = blift * cadet;
	IT lc     = clift * abdet;
	IT lab    = la + lb;
	IT L      = lab + lc;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename ET>
Sign inCircle_exact(ET pax, ET pay, ET pbx, ET pby, ET pcx, ET pcy, ET pdx,
                    ET pdy)
{
	ET adx    = pax - pdx;
	ET ady    = pay - pdy;
	ET bdx    = pbx - pdx;
	ET bdy    = pby - pdy;
	ET cdx    = pcx - pdx;
	ET cdy    = pcy - pdy;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET alift  = alifta + aliftb;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET blift  = blifta + bliftb;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET clift  = clifta + cliftb;
	ET la     = alift * bcdet;
	ET lb     = blift * cadet;
	ET lc     = clift * abdet;
	ET lab    = la + lb;
	ET L      = lab + lc;
	return OMC::sign(L);
}

Sign inCircle_expansion(double pax, double pay, double pbx, double pby,
                        double pcx, double pcy, double pdx, double pdy)
{
	expansionObject o;
	double          adx[2];
	o.two_Diff(pax, pdx, adx);
	double ady[2];
	o.two_Diff(pay, pdy, ady);
	double bdx[2];
	o.two_Diff(pbx, pdx, bdx);
	double bdy[2];
	o.two_Diff(pby, pdy, bdy);
	double cdx[2];
	o.two_Diff(pcx, pdx, cdx);
	double cdy[2];
	o.two_Diff(pcy, pdy, cdy);
	double abdeta[8];
	int    abdeta_len = o.Gen_Product(2, adx, 2, bdy, abdeta);
	double abdetb[8];
	int    abdetb_len = o.Gen_Product(2, bdx, 2, ady, abdetb);
	double abdet[16];
	int    abdet_len = o.Gen_Diff(abdeta_len, abdeta, abdetb_len, abdetb, abdet);
	double bcdeta[8];
	int    bcdeta_len = o.Gen_Product(2, bdx, 2, cdy, bcdeta);
	double bcdetb[8];
	int    bcdetb_len = o.Gen_Product(2, cdx, 2, bdy, bcdetb);
	double bcdet[16];
	int    bcdet_len = o.Gen_Diff(bcdeta_len, bcdeta, bcdetb_len, bcdetb, bcdet);
	double cadeta[8];
	int    cadeta_len = o.Gen_Product(2, cdx, 2, ady, cadeta);
	double cadetb[8];
	int    cadetb_len = o.Gen_Product(2, adx, 2, cdy, cadetb);
	double cadet[16];
	int    cadet_len = o.Gen_Diff(cadeta_len, cadeta, cadetb_len, cadetb, cadet);
	double alifta[8];
	int    alifta_len = o.Gen_Product(2, adx, 2, adx, alifta);
	double aliftb[8];
	int    aliftb_len = o.Gen_Product(2, ady, 2, ady, aliftb);
	double alift[16];
	int    alift_len = o.Gen_Sum(alifta_len, alifta, aliftb_len, aliftb, alift);
	double blifta[8];
	int    blifta_len = o.Gen_Product(2, bdx, 2, bdx, blifta);
	double bliftb[8];
	int    bliftb_len = o.Gen_Product(2, bdy, 2, bdy, bliftb);
	double blift[16];
	int    blift_len = o.Gen_Sum(blifta_len, blifta, bliftb_len, bliftb, blift);
	double clifta[8];
	int    clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
	double cliftb[8];
	int    cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
	double clift[16];
	int    clift_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, clift);
	double la_p[128], *la = la_p;
	int    la_len =
	  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 128);
	double lb_p[128], *lb = lb_p;
	int    lb_len =
	  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 128);
	double lc_p[128], *lc = lc_p;
	int    lc_len =
	  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 128);
	double lab_p[128], *lab = lab_p;
	int    lab_len = o.Gen_Sum_With_PreAlloc(la_len, la, lb_len, lb, &lab, 128);
	double L_p[128], *L = L_p;
	int    L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, lc_len, lc, &L, 128);

	double return_value = L[L_len - 1];
	if (L_p != L)
		FreeDoubles(L);
	if (lab_p != lab)
		FreeDoubles(lab);
	if (lc_p != lc)
		FreeDoubles(lc);
	if (lb_p != lb)
		FreeDoubles(lb);
	if (la_p != la)
		FreeDoubles(la);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign inCircle(double pax, double pay, double pbx, double pby, double pcx,
              double pcy, double pdx, double pdy)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = inCircle_filtered(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = inCircle_interval<IT>(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
	if (is_sign_reliable(ret))
		return ret;
	return inCircle_expansion(pax, pay, pbx, pby, pcx, pcy, pdx, pdy);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign inCircle(const GenericPoint2T<IT, ET> &pa,
              const GenericPoint2T<IT, ET> &pb,
              const GenericPoint2T<IT, ET> &pc,
              const GenericPoint2T<IT, ET> &pd)
{
	return inCircle<IT, ET, WithSSFilter>(pa.x(), pa.y(), pb.x(), pb.y(), pc.x(),
	                                      pc.y(), pd.x(), pd.y());
}

inline Sign inSphere_filtered(double pax, double pay, double paz, double pbx,
                              double pby, double pbz, double pcx, double pcy,
                              double pcz, double pdx, double pdy, double pdz,
                              double pex, double pey, double pez)
{
	double aex    = pax - pex;
	double aey    = pay - pey;
	double aez    = paz - pez;
	double bex    = pbx - pex;
	double bey    = pby - pey;
	double bez    = pbz - pez;
	double cex    = pcx - pex;
	double cey    = pcy - pey;
	double cez    = pcz - pez;
	double dex    = pdx - pex;
	double dey    = pdy - pey;
	double dez    = pdz - pez;
	double aexbey = aex * bey;
	double bexaey = bex * aey;
	double ab     = aexbey - bexaey;
	double bexcey = bex * cey;
	double cexbey = cex * bey;
	double bc     = bexcey - cexbey;
	double cexdey = cex * dey;
	double dexcey = dex * cey;
	double cd     = cexdey - dexcey;
	double dexaey = dex * aey;
	double aexdey = aex * dey;
	double da     = dexaey - aexdey;
	double aexcey = aex * cey;
	double cexaey = cex * aey;
	double ac     = aexcey - cexaey;
	double bexdey = bex * dey;
	double dexbey = dex * bey;
	double bd     = bexdey - dexbey;
	double abc1   = aez * bc;
	double abc2   = bez * ac;
	double abc3   = cez * ab;
	double abc4   = abc1 + abc3;
	double abc    = abc4 - abc2;
	double bcd1   = bez * cd;
	double bcd2   = cez * bd;
	double bcd3   = dez * bc;
	double bcd4   = bcd1 + bcd3;
	double bcd    = bcd4 - bcd2;
	double cda1   = cez * da;
	double cda2   = dez * ac;
	double cda3   = aez * cd;
	double cda4   = cda1 + cda3;
	double cda    = cda4 + cda2;
	double dab1   = dez * ab;
	double dab2   = aez * bd;
	double dab3   = bez * da;
	double dab4   = dab1 + dab3;
	double dab    = dab4 + dab2;
	double al1    = aex * aex;
	double al2    = aey * aey;
	double al3    = aez * aez;
	double al4    = al1 + al2;
	double alift  = al4 + al3;
	double bl1    = bex * bex;
	double bl2    = bey * bey;
	double bl3    = bez * bez;
	double bl4    = bl1 + bl2;
	double blift  = bl4 + bl3;
	double cl1    = cex * cex;
	double cl2    = cey * cey;
	double cl3    = cez * cez;
	double cl4    = cl1 + cl2;
	double clift  = cl4 + cl3;
	double dl1    = dex * dex;
	double dl2    = dey * dey;
	double dl3    = dez * dez;
	double dl4    = dl1 + dl2;
	double dlift  = dl4 + dl3;
	double ds1    = dlift * abc;
	double ds2    = clift * dab;
	double dl     = ds2 - ds1;
	double dr1    = blift * cda;
	double dr2    = alift * bcd;
	double dr     = dr2 - dr1;
	double det    = dl + dr;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(aex)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(aey)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(aez)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(bex)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(bey)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(bez)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(cex)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(cey)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(cez)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(dex)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(dey)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(dez)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= epsilon;
	epsilon *= max_var;
	epsilon *= 1.1457501614131623e-13;

	return filter_sign(det, epsilon);
}

template <typename IT>
Sign inSphere_interval(IT pax, IT pay, IT paz, IT pbx, IT pby, IT pbz, IT pcx,
                       IT pcy, IT pcz, IT pdx, IT pdy, IT pdz, IT pex, IT pey,
                       IT pez)
{
	typename IT::Protector P;

	IT aex    = pax - pex;
	IT aey    = pay - pey;
	IT aez    = paz - pez;
	IT bex    = pbx - pex;
	IT bey    = pby - pey;
	IT bez    = pbz - pez;
	IT cex    = pcx - pex;
	IT cey    = pcy - pey;
	IT cez    = pcz - pez;
	IT dex    = pdx - pex;
	IT dey    = pdy - pey;
	IT dez    = pdz - pez;
	IT aexbey = aex * bey;
	IT bexaey = bex * aey;
	IT ab     = aexbey - bexaey;
	IT bexcey = bex * cey;
	IT cexbey = cex * bey;
	IT bc     = bexcey - cexbey;
	IT cexdey = cex * dey;
	IT dexcey = dex * cey;
	IT cd     = cexdey - dexcey;
	IT dexaey = dex * aey;
	IT aexdey = aex * dey;
	IT da     = dexaey - aexdey;
	IT aexcey = aex * cey;
	IT cexaey = cex * aey;
	IT ac     = aexcey - cexaey;
	IT bexdey = bex * dey;
	IT dexbey = dex * bey;
	IT bd     = bexdey - dexbey;
	IT abc1   = aez * bc;
	IT abc2   = bez * ac;
	IT abc3   = cez * ab;
	IT abc4   = abc1 + abc3;
	IT abc    = abc4 - abc2;
	IT bcd1   = bez * cd;
	IT bcd2   = cez * bd;
	IT bcd3   = dez * bc;
	IT bcd4   = bcd1 + bcd3;
	IT bcd    = bcd4 - bcd2;
	IT cda1   = cez * da;
	IT cda2   = dez * ac;
	IT cda3   = aez * cd;
	IT cda4   = cda1 + cda3;
	IT cda    = cda4 + cda2;
	IT dab1   = dez * ab;
	IT dab2   = aez * bd;
	IT dab3   = bez * da;
	IT dab4   = dab1 + dab3;
	IT dab    = dab4 + dab2;
	IT al1    = aex * aex;
	IT al2    = aey * aey;
	IT al3    = aez * aez;
	IT al4    = al1 + al2;
	IT alift  = al4 + al3;
	IT bl1    = bex * bex;
	IT bl2    = bey * bey;
	IT bl3    = bez * bez;
	IT bl4    = bl1 + bl2;
	IT blift  = bl4 + bl3;
	IT cl1    = cex * cex;
	IT cl2    = cey * cey;
	IT cl3    = cez * cez;
	IT cl4    = cl1 + cl2;
	IT clift  = cl4 + cl3;
	IT dl1    = dex * dex;
	IT dl2    = dey * dey;
	IT dl3    = dez * dez;
	IT dl4    = dl1 + dl2;
	IT dlift  = dl4 + dl3;
	IT ds1    = dlift * abc;
	IT ds2    = clift * dab;
	IT dl     = ds2 - ds1;
	IT dr1    = blift * cda;
	IT dr2    = alift * bcd;
	IT dr     = dr2 - dr1;
	IT det    = dl + dr;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename ET>
Sign inSphere_exact(ET pax, ET pay, ET paz, ET pbx, ET pby, ET pbz, ET pcx,
                    ET pcy, ET pcz, ET pdx, ET pdy, ET pdz, ET pex, ET pey,
                    ET pez)
{
	ET aex    = pax - pex;
	ET aey    = pay - pey;
	ET aez    = paz - pez;
	ET bex    = pbx - pex;
	ET bey    = pby - pey;
	ET bez    = pbz - pez;
	ET cex    = pcx - pex;
	ET cey    = pcy - pey;
	ET cez    = pcz - pez;
	ET dex    = pdx - pex;
	ET dey    = pdy - pey;
	ET dez    = pdz - pez;
	ET aexbey = aex * bey;
	ET bexaey = bex * aey;
	ET ab     = aexbey - bexaey;
	ET bexcey = bex * cey;
	ET cexbey = cex * bey;
	ET bc     = bexcey - cexbey;
	ET cexdey = cex * dey;
	ET dexcey = dex * cey;
	ET cd     = cexdey - dexcey;
	ET dexaey = dex * aey;
	ET aexdey = aex * dey;
	ET da     = dexaey - aexdey;
	ET aexcey = aex * cey;
	ET cexaey = cex * aey;
	ET ac     = aexcey - cexaey;
	ET bexdey = bex * dey;
	ET dexbey = dex * bey;
	ET bd     = bexdey - dexbey;
	ET abc1   = aez * bc;
	ET abc2   = bez * ac;
	ET abc3   = cez * ab;
	ET abc4   = abc1 + abc3;
	ET abc    = abc4 - abc2;
	ET bcd1   = bez * cd;
	ET bcd2   = cez * bd;
	ET bcd3   = dez * bc;
	ET bcd4   = bcd1 + bcd3;
	ET bcd    = bcd4 - bcd2;
	ET cda1   = cez * da;
	ET cda2   = dez * ac;
	ET cda3   = aez * cd;
	ET cda4   = cda1 + cda3;
	ET cda    = cda4 + cda2;
	ET dab1   = dez * ab;
	ET dab2   = aez * bd;
	ET dab3   = bez * da;
	ET dab4   = dab1 + dab3;
	ET dab    = dab4 + dab2;
	ET al1    = aex * aex;
	ET al2    = aey * aey;
	ET al3    = aez * aez;
	ET al4    = al1 + al2;
	ET alift  = al4 + al3;
	ET bl1    = bex * bex;
	ET bl2    = bey * bey;
	ET bl3    = bez * bez;
	ET bl4    = bl1 + bl2;
	ET blift  = bl4 + bl3;
	ET cl1    = cex * cex;
	ET cl2    = cey * cey;
	ET cl3    = cez * cez;
	ET cl4    = cl1 + cl2;
	ET clift  = cl4 + cl3;
	ET dl1    = dex * dex;
	ET dl2    = dey * dey;
	ET dl3    = dez * dez;
	ET dl4    = dl1 + dl2;
	ET dlift  = dl4 + dl3;
	ET ds1    = dlift * abc;
	ET ds2    = clift * dab;
	ET dl     = ds2 - ds1;
	ET dr1    = blift * cda;
	ET dr2    = alift * bcd;
	ET dr     = dr2 - dr1;
	ET det    = dl + dr;
	return OMC::sign(det);
}

Sign inSphere_expansion(double pax, double pay, double paz, double pbx,
                        double pby, double pbz, double pcx, double pcy,
                        double pcz, double pdx, double pdy, double pdz,
                        double pex, double pey, double pez)
{
	expansionObject o;
	double          aex[2];
	o.two_Diff(pax, pex, aex);
	double aey[2];
	o.two_Diff(pay, pey, aey);
	double aez[2];
	o.two_Diff(paz, pez, aez);
	double bex[2];
	o.two_Diff(pbx, pex, bex);
	double bey[2];
	o.two_Diff(pby, pey, bey);
	double bez[2];
	o.two_Diff(pbz, pez, bez);
	double cex[2];
	o.two_Diff(pcx, pex, cex);
	double cey[2];
	o.two_Diff(pcy, pey, cey);
	double cez[2];
	o.two_Diff(pcz, pez, cez);
	double dex[2];
	o.two_Diff(pdx, pex, dex);
	double dey[2];
	o.two_Diff(pdy, pey, dey);
	double dez[2];
	o.two_Diff(pdz, pez, dez);
	double aexbey[8];
	int    aexbey_len = o.Gen_Product(2, aex, 2, bey, aexbey);
	double bexaey[8];
	int    bexaey_len = o.Gen_Product(2, bex, 2, aey, bexaey);
	double ab[16];
	int    ab_len = o.Gen_Diff(aexbey_len, aexbey, bexaey_len, bexaey, ab);
	double bexcey[8];
	int    bexcey_len = o.Gen_Product(2, bex, 2, cey, bexcey);
	double cexbey[8];
	int    cexbey_len = o.Gen_Product(2, cex, 2, bey, cexbey);
	double bc[16];
	int    bc_len = o.Gen_Diff(bexcey_len, bexcey, cexbey_len, cexbey, bc);
	double cexdey[8];
	int    cexdey_len = o.Gen_Product(2, cex, 2, dey, cexdey);
	double dexcey[8];
	int    dexcey_len = o.Gen_Product(2, dex, 2, cey, dexcey);
	double cd[16];
	int    cd_len = o.Gen_Diff(cexdey_len, cexdey, dexcey_len, dexcey, cd);
	double dexaey[8];
	int    dexaey_len = o.Gen_Product(2, dex, 2, aey, dexaey);
	double aexdey[8];
	int    aexdey_len = o.Gen_Product(2, aex, 2, dey, aexdey);
	double da[16];
	int    da_len = o.Gen_Diff(dexaey_len, dexaey, aexdey_len, aexdey, da);
	double aexcey[8];
	int    aexcey_len = o.Gen_Product(2, aex, 2, cey, aexcey);
	double cexaey[8];
	int    cexaey_len = o.Gen_Product(2, cex, 2, aey, cexaey);
	double ac[16];
	int    ac_len = o.Gen_Diff(aexcey_len, aexcey, cexaey_len, cexaey, ac);
	double bexdey[8];
	int    bexdey_len = o.Gen_Product(2, bex, 2, dey, bexdey);
	double dexbey[8];
	int    dexbey_len = o.Gen_Product(2, dex, 2, bey, dexbey);
	double bd[16];
	int    bd_len = o.Gen_Diff(bexdey_len, bexdey, dexbey_len, dexbey, bd);
	double abc1_p[32], *abc1 = abc1_p;
	int    abc1_len = o.Gen_Product_With_PreAlloc(2, aez, bc_len, bc, &abc1, 32);
	double abc2_p[32], *abc2 = abc2_p;
	int    abc2_len = o.Gen_Product_With_PreAlloc(2, bez, ac_len, ac, &abc2, 32);
	double abc3_p[32], *abc3 = abc3_p;
	int    abc3_len = o.Gen_Product_With_PreAlloc(2, cez, ab_len, ab, &abc3, 32);
	double abc4_p[32], *abc4 = abc4_p;
	int    abc4_len =
	  o.Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 32);
	double abc_p[32], *abc = abc_p;
	int    abc_len =
	  o.Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 32);
	double bcd1_p[32], *bcd1 = bcd1_p;
	int    bcd1_len = o.Gen_Product_With_PreAlloc(2, bez, cd_len, cd, &bcd1, 32);
	double bcd2_p[32], *bcd2 = bcd2_p;
	int    bcd2_len = o.Gen_Product_With_PreAlloc(2, cez, bd_len, bd, &bcd2, 32);
	double bcd3_p[32], *bcd3 = bcd3_p;
	int    bcd3_len = o.Gen_Product_With_PreAlloc(2, dez, bc_len, bc, &bcd3, 32);
	double bcd4_p[32], *bcd4 = bcd4_p;
	int    bcd4_len =
	  o.Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 32);
	double bcd_p[32], *bcd = bcd_p;
	int    bcd_len =
	  o.Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 32);
	double cda1_p[32], *cda1 = cda1_p;
	int    cda1_len = o.Gen_Product_With_PreAlloc(2, cez, da_len, da, &cda1, 32);
	double cda2_p[32], *cda2 = cda2_p;
	int    cda2_len = o.Gen_Product_With_PreAlloc(2, dez, ac_len, ac, &cda2, 32);
	double cda3_p[32], *cda3 = cda3_p;
	int    cda3_len = o.Gen_Product_With_PreAlloc(2, aez, cd_len, cd, &cda3, 32);
	double cda4_p[32], *cda4 = cda4_p;
	int    cda4_len =
	  o.Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 32);
	double cda_p[32], *cda = cda_p;
	int    cda_len =
	  o.Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 32);
	double dab1_p[32], *dab1 = dab1_p;
	int    dab1_len = o.Gen_Product_With_PreAlloc(2, dez, ab_len, ab, &dab1, 32);
	double dab2_p[32], *dab2 = dab2_p;
	int    dab2_len = o.Gen_Product_With_PreAlloc(2, aez, bd_len, bd, &dab2, 32);
	double dab3_p[32], *dab3 = dab3_p;
	int    dab3_len = o.Gen_Product_With_PreAlloc(2, bez, da_len, da, &dab3, 32);
	double dab4_p[32], *dab4 = dab4_p;
	int    dab4_len =
	  o.Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 32);
	double dab_p[32], *dab = dab_p;
	int    dab_len =
	  o.Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 32);
	double al1[8];
	int    al1_len = o.Gen_Product(2, aex, 2, aex, al1);
	double al2[8];
	int    al2_len = o.Gen_Product(2, aey, 2, aey, al2);
	double al3[8];
	int    al3_len = o.Gen_Product(2, aez, 2, aez, al3);
	double al4[16];
	int    al4_len = o.Gen_Sum(al1_len, al1, al2_len, al2, al4);
	double alift[24];
	int    alift_len = o.Gen_Sum(al4_len, al4, al3_len, al3, alift);
	double bl1[8];
	int    bl1_len = o.Gen_Product(2, bex, 2, bex, bl1);
	double bl2[8];
	int    bl2_len = o.Gen_Product(2, bey, 2, bey, bl2);
	double bl3[8];
	int    bl3_len = o.Gen_Product(2, bez, 2, bez, bl3);
	double bl4[16];
	int    bl4_len = o.Gen_Sum(bl1_len, bl1, bl2_len, bl2, bl4);
	double blift[24];
	int    blift_len = o.Gen_Sum(bl4_len, bl4, bl3_len, bl3, blift);
	double cl1[8];
	int    cl1_len = o.Gen_Product(2, cex, 2, cex, cl1);
	double cl2[8];
	int    cl2_len = o.Gen_Product(2, cey, 2, cey, cl2);
	double cl3[8];
	int    cl3_len = o.Gen_Product(2, cez, 2, cez, cl3);
	double cl4[16];
	int    cl4_len = o.Gen_Sum(cl1_len, cl1, cl2_len, cl2, cl4);
	double clift[24];
	int    clift_len = o.Gen_Sum(cl4_len, cl4, cl3_len, cl3, clift);
	double dl1[8];
	int    dl1_len = o.Gen_Product(2, dex, 2, dex, dl1);
	double dl2[8];
	int    dl2_len = o.Gen_Product(2, dey, 2, dey, dl2);
	double dl3[8];
	int    dl3_len = o.Gen_Product(2, dez, 2, dez, dl3);
	double dl4[16];
	int    dl4_len = o.Gen_Sum(dl1_len, dl1, dl2_len, dl2, dl4);
	double dlift[24];
	int    dlift_len = o.Gen_Sum(dl4_len, dl4, dl3_len, dl3, dlift);
	double ds1_p[32], *ds1 = ds1_p;
	int    ds1_len =
	  o.Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 32);
	double ds2_p[32], *ds2 = ds2_p;
	int    ds2_len =
	  o.Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 32);
	double dl_p[32], *dl = dl_p;
	int    dl_len = o.Gen_Diff_With_PreAlloc(ds2_len, ds2, ds1_len, ds1, &dl, 32);
	double dr1_p[32], *dr1 = dr1_p;
	int    dr1_len =
	  o.Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1, 32);
	double dr2_p[32], *dr2 = dr2_p;
	int    dr2_len =
	  o.Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 32);
	double dr_p[32], *dr = dr_p;
	int    dr_len = o.Gen_Diff_With_PreAlloc(dr2_len, dr2, dr1_len, dr1, &dr, 32);
	double det_p[32], *det = det_p;
	int    det_len = o.Gen_Sum_With_PreAlloc(dl_len, dl, dr_len, dr, &det, 32);

	double return_value = det[det_len - 1];
	if (det_p != det)
		FreeDoubles(det);
	if (dr_p != dr)
		FreeDoubles(dr);
	if (dr2_p != dr2)
		FreeDoubles(dr2);
	if (dr1_p != dr1)
		FreeDoubles(dr1);
	if (dl_p != dl)
		FreeDoubles(dl);
	if (ds2_p != ds2)
		FreeDoubles(ds2);
	if (ds1_p != ds1)
		FreeDoubles(ds1);
	if (dab_p != dab)
		FreeDoubles(dab);
	if (dab4_p != dab4)
		FreeDoubles(dab4);
	if (dab3_p != dab3)
		FreeDoubles(dab3);
	if (dab2_p != dab2)
		FreeDoubles(dab2);
	if (dab1_p != dab1)
		FreeDoubles(dab1);
	if (cda_p != cda)
		FreeDoubles(cda);
	if (cda4_p != cda4)
		FreeDoubles(cda4);
	if (cda3_p != cda3)
		FreeDoubles(cda3);
	if (cda2_p != cda2)
		FreeDoubles(cda2);
	if (cda1_p != cda1)
		FreeDoubles(cda1);
	if (bcd_p != bcd)
		FreeDoubles(bcd);
	if (bcd4_p != bcd4)
		FreeDoubles(bcd4);
	if (bcd3_p != bcd3)
		FreeDoubles(bcd3);
	if (bcd2_p != bcd2)
		FreeDoubles(bcd2);
	if (bcd1_p != bcd1)
		FreeDoubles(bcd1);
	if (abc_p != abc)
		FreeDoubles(abc);
	if (abc4_p != abc4)
		FreeDoubles(abc4);
	if (abc3_p != abc3)
		FreeDoubles(abc3);
	if (abc2_p != abc2)
		FreeDoubles(abc2);
	if (abc1_p != abc1)
		FreeDoubles(abc1);

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign inSphere(double pax, double pay, double paz, double pbx, double pby,
              double pbz, double pcx, double pcy, double pcz, double pdx,
              double pdy, double pdz, double pex, double pey, double pez)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = inSphere_filtered(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx,
		                        pdy, pdz, pex, pey, pez);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = inSphere_interval<IT>(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx,
	                            pdy, pdz, pex, pey, pez);
	if (is_sign_reliable(ret))
		return ret;
	return inSphere_expansion(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx,
	                          pdy, pdz, pex, pey, pez);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign inSphere(const GenericPoint3T<IT, ET> &pa,
              const GenericPoint3T<IT, ET> &pb,
              const GenericPoint3T<IT, ET> &pc,
              const GenericPoint3T<IT, ET> &pd,
              const GenericPoint3T<IT, ET> &pe)
{
	return inSphere<IT, ET, WithSSFilter>(pa.x(), pa.y(), pa.z(), pb.x(), pb.y(),
	                                      pb.z(), pc.x(), pc.y(), pc.z(), pd.x(),
	                                      pd.y(), pd.z(), pe.x(), pe.y(), pe.z());
}

inline Sign squareDistance2D_filtered(double px, double py, double qx,
                                      double qy, double dis)
{
	double lx      = px - qx;
	double ly      = py - qy;
	double lx2     = lx * lx;
	double ly2     = ly * ly;
	double sqrnorm = lx2 + ly2;
	double diff    = sqrnorm - dis;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(dis)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 1.1102230246251565e-15;

	return filter_sign(diff, epsilon);
}

template <typename IT>
Sign squareDistance2D_interval(IT px, IT py, IT qx, IT qy, IT dis)
{
	typename IT::Protector P;

	IT lx      = px - qx;
	IT ly      = py - qy;
	IT lx2     = lx * lx;
	IT ly2     = ly * ly;
	IT sqrnorm = lx2 + ly2;
	IT diff    = sqrnorm - dis;
	if (!diff.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(diff);
}

template <typename ET>
Sign squareDistance2D_exact(ET px, ET py, ET qx, ET qy, ET dis)
{
	ET lx      = px - qx;
	ET ly      = py - qy;
	ET lx2     = lx * lx;
	ET ly2     = ly * ly;
	ET sqrnorm = lx2 + ly2;
	ET diff    = sqrnorm - dis;
	return OMC::sign(diff);
}

Sign squareDistance2D_expansion(double px, double py, double qx, double qy,
                                double dis)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double ly[2];
	o.two_Diff(py, qy, ly);
	double lx2[8];
	int    lx2_len = o.Gen_Product(2, lx, 2, lx, lx2);
	double ly2[8];
	int    ly2_len = o.Gen_Product(2, ly, 2, ly, ly2);
	double sqrnorm[16];
	int    sqrnorm_len = o.Gen_Sum(lx2_len, lx2, ly2_len, ly2, sqrnorm);
	double diff[17];
	int    diff_len = o.Gen_Diff(sqrnorm_len, sqrnorm, 1, &dis, diff);

	double return_value = diff[diff_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance2D(double px, double py, double qx, double qy, double dis)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = squareDistance2D_filtered(px, py, qx, qy, dis);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = squareDistance2D_interval<IT>(px, py, qx, qy, dis);
	if (is_sign_reliable(ret))
		return ret;
	return squareDistance2D_expansion(px, py, qx, qy, dis);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance2D(const GenericPoint2T<IT, ET> &p,
                      const GenericPoint2T<IT, ET> &q, double dis)
{
	return squareDistance2D<IT, ET, WithSSFilter>(p.x(), p.y(), q.x(), q.y(),
	                                              dis);
}

inline Sign squareDistance3D_filtered(double px, double py, double pz,
                                      double qx, double qy, double qz,
                                      double dis)
{
	double lx        = px - qx;
	double ly        = py - qy;
	double lz        = pz - qz;
	double lx2       = lx * lx;
	double ly2       = ly * ly;
	double lz2       = lz * lz;
	double sqrnormxy = lx2 + ly2;
	double sqrnorm   = sqrnormxy + lz2;
	double diff      = sqrnorm - dis;

	double _tmp_fabs;

	double max_var = 0.0;
	if ((_tmp_fabs = fabs(dis)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(lx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ly)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(lz)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	epsilon *= epsilon;
	epsilon *= 1.887379141862766e-15;

	return filter_sign(diff, epsilon);
}

template <typename IT>
Sign squareDistance3D_interval(IT px, IT py, IT pz, IT qx, IT qy, IT qz, IT dis)
{
	typename IT::Protector P;

	IT lx        = px - qx;
	IT ly        = py - qy;
	IT lz        = pz - qz;
	IT lx2       = lx * lx;
	IT ly2       = ly * ly;
	IT lz2       = lz * lz;
	IT sqrnormxy = lx2 + ly2;
	IT sqrnorm   = sqrnormxy + lz2;
	IT diff      = sqrnorm - dis;
	if (!diff.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(diff);
}

template <typename ET>
Sign squareDistance3D_exact(ET px, ET py, ET pz, ET qx, ET qy, ET qz, ET dis)
{
	ET lx        = px - qx;
	ET ly        = py - qy;
	ET lz        = pz - qz;
	ET lx2       = lx * lx;
	ET ly2       = ly * ly;
	ET lz2       = lz * lz;
	ET sqrnormxy = lx2 + ly2;
	ET sqrnorm   = sqrnormxy + lz2;
	ET diff      = sqrnorm - dis;
	return OMC::sign(diff);
}

Sign squareDistance3D_expansion(double px, double py, double pz, double qx,
                                double qy, double qz, double dis)
{
	expansionObject o;
	double          lx[2];
	o.two_Diff(px, qx, lx);
	double ly[2];
	o.two_Diff(py, qy, ly);
	double lz[2];
	o.two_Diff(pz, qz, lz);
	double lx2[8];
	int    lx2_len = o.Gen_Product(2, lx, 2, lx, lx2);
	double ly2[8];
	int    ly2_len = o.Gen_Product(2, ly, 2, ly, ly2);
	double lz2[8];
	int    lz2_len = o.Gen_Product(2, lz, 2, lz, lz2);
	double sqrnormxy[16];
	int    sqrnormxy_len = o.Gen_Sum(lx2_len, lx2, ly2_len, ly2, sqrnormxy);
	double sqrnorm[24];
	int sqrnorm_len = o.Gen_Sum(sqrnormxy_len, sqrnormxy, lz2_len, lz2, sqrnorm);
	double diff[25];
	int    diff_len = o.Gen_Diff(sqrnorm_len, sqrnorm, 1, &dis, diff);

	double return_value = diff[diff_len - 1];

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance3D(double px, double py, double pz, double qx, double qy,
                      double qz, double dis)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = squareDistance3D_filtered(px, py, pz, qx, qy, qz, dis);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = squareDistance3D_interval<IT>(px, py, pz, qx, qy, qz, dis);
	if (is_sign_reliable(ret))
		return ret;
	return squareDistance3D_expansion(px, py, pz, qx, qy, qz, dis);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance3D(const GenericPoint3T<IT, ET> &p,
                      const GenericPoint3T<IT, ET> &q, double dis)
{
	return squareDistance3D<IT, ET, WithSSFilter>(p.x(), p.y(), p.z(), q.x(),
	                                              q.y(), q.z(), dis);
}

#if defined(INDIRECT_PREDICATES)

template <typename IT, typename ET>
Sign dotProductSign2D_EEI_interval(const GenericPoint2T<IT, ET> &q, IT px,
                                   IT py, IT rx, IT ry)
{
	IT lqx, lqy, dq;
	if (!q.getIntervalLambda(lqx, lqy, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pxq = px * dq;
	IT pyq = py * dq;
	IT rxq = rx * dq;
	IT ryq = ry * dq;
	IT lx  = pxq - lqx;
	IT ly  = pyq - lqy;
	IT gx  = rxq - lqx;
	IT gy  = ryq - lqy;
	IT dx  = lx * gx;
	IT dy  = ly * gy;
	IT d   = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_EEI_exact(const GenericPoint2T<IT, ET> &q, ET px, ET py,
                                ET rx, ET ry)
{
	ET lqx, lqy, dq;
	q.getExactLambda(lqx, lqy, dq);
	ET pxq = px * dq;
	ET pyq = py * dq;
	ET rxq = rx * dq;
	ET ryq = ry * dq;
	ET lx  = pxq - lqx;
	ET ly  = pyq - lqy;
	ET gx  = rxq - lqx;
	ET gy  = ryq - lqy;
	ET dx  = lx * gx;
	ET dy  = ly * gy;
	ET d   = dx + dy;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_EEI_expansion(const GenericPoint2T<IT, ET> &q, double px,
                                    double py, double rx, double ry)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lqx_p[128], *lqx = lqx_p, lqy_p[128], *lqy = lqy_p, dq_p[128],
	                   *dq = dq_p;
	int lqx_len = 128, lqy_len = 128, dq_len = 128;
	q.getExpansionLambda(&lqx, lqx_len, &lqy, lqy_len, &dq, dq_len);
	if ((dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          pxq_p[128], *pxq = pxq_p;
		int    pxq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, px, &pxq, 128);
		double pyq_p[128], *pyq = pyq_p;
		int    pyq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, py, &pyq, 128);
		double rxq_p[128], *rxq = rxq_p;
		int    rxq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, rx, &rxq, 128);
		double ryq_p[128], *ryq = ryq_p;
		int    ryq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, ry, &ryq, 128);
		double lx_p[128], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(pxq_len, pxq, lqx_len, lqx, &lx, 128);
		double ly_p[128], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(pyq_len, pyq, lqy_len, lqy, &ly, 128);
		double gx_p[128], *gx = gx_p;
		int gx_len = o.Gen_Diff_With_PreAlloc(rxq_len, rxq, lqx_len, lqx, &gx, 128);
		double gy_p[128], *gy = gy_p;
		int gy_len = o.Gen_Diff_With_PreAlloc(ryq_len, ryq, lqy_len, lqy, &gy, 128);
		double dx_p[128], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 128);
		double dy_p[128], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 128);
		double d_p[128], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d, 128);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (ryq_p != ryq)
			FreeDoubles(ryq);
		if (rxq_p != rxq)
			FreeDoubles(rxq);
		if (pyq_p != pyq)
			FreeDoubles(pyq);
		if (pxq_p != pxq)
			FreeDoubles(pxq);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lqx_p != lqx)
			FreeDoubles(lqx);
		if (lqy_p != lqy)
			FreeDoubles(lqy);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign2D_EEI_exact<IT, ET>(q, px, py, rx, ry);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign2D_EEI(const GenericPoint2T<IT, ET> &q, double px, double py,
                          double rx, double ry)
{
	Sign ret;
	ret = dotProductSign2D_EEI_interval<IT, ET>(q, px, py, rx, ry);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign2D_EEI_expansion<IT, ET>(q, px, py, rx, ry);
}

template <typename IT, typename ET>
Sign dotProductSign2D_EEI(const GenericPoint2T<IT, ET> &q,
                          const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r)
{
	return dotProductSign2D_EEI<IT, ET>(q, p.x(), p.y(), r.x(), r.y());
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEE_interval(const GenericPoint2T<IT, ET> &p, IT rx,
                                   IT ry, IT qx, IT qy)
{
	IT lpx, lpy, dp;
	if (!p.getIntervalLambda(lpx, lpy, dp))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT qxd = qx * dp;
	IT qyd = qy * dp;
	IT lx  = lpx - qxd;
	IT ly  = lpy - qyd;
	IT gx  = rx - qx;
	IT gy  = ry - qy;
	IT dx  = lx * gx;
	IT dy  = ly * gy;
	IT d   = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEE_exact(const GenericPoint2T<IT, ET> &p, ET rx, ET ry,
                                ET qx, ET qy)
{
	ET lpx, lpy, dp;
	p.getExactLambda(lpx, lpy, dp);
	ET qxd = qx * dp;
	ET qyd = qy * dp;
	ET lx  = lpx - qxd;
	ET ly  = lpy - qyd;
	ET gx  = rx - qx;
	ET gy  = ry - qy;
	ET dx  = lx * gx;
	ET dy  = ly * gy;
	ET d   = dx + dy;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEE_expansion(const GenericPoint2T<IT, ET> &p, double rx,
                                    double ry, double qx, double qy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[128], *lpx = lpx_p, lpy_p[128], *lpy = lpy_p, dp_p[128],
	                   *dp = dp_p;
	int lpx_len = 128, lpy_len = 128, dp_len = 128;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &dp, dp_len);
	if ((dp[dp_len - 1] != 0))
	{
		expansionObject o;
		double          qxd_p[128], *qxd = qxd_p;
		int    qxd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qx, &qxd, 128);
		double qyd_p[128], *qyd = qyd_p;
		int    qyd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qy, &qyd, 128);
		double lx_p[128], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(lpx_len, lpx, qxd_len, qxd, &lx, 128);
		double ly_p[128], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(lpy_len, lpy, qyd_len, qyd, &ly, 128);
		double gx[2];
		o.two_Diff(rx, qx, gx);
		double gy[2];
		o.two_Diff(ry, qy, gy);
		double dx_p[128], *dx = dx_p;
		int    dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, 2, gx, &dx, 128);
		double dy_p[128], *dy = dy_p;
		int    dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, 2, gy, &dy, 128);
		double d_p[128], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d, 128);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (qyd_p != qyd)
			FreeDoubles(qyd);
		if (qxd_p != qxd)
			FreeDoubles(qxd);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign2D_IEE_exact<IT, ET>(p, rx, ry, qx, qy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEE(const GenericPoint2T<IT, ET> &p, double rx, double ry,
                          double qx, double qy)
{
	Sign ret;
	ret = dotProductSign2D_IEE_interval<IT, ET>(p, rx, ry, qx, qy);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign2D_IEE_expansion<IT, ET>(p, rx, ry, qx, qy);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEE(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r,
                          const GenericPoint2T<IT, ET> &q)
{
	return dotProductSign2D_IEE<IT, ET>(p, r.x(), r.y(), q.x(), q.y());
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEI_interval(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &q, IT rx,
                                   IT ry)
{
	IT lpx, lpy, dp, lqx, lqy, dq;
	if (!p.getIntervalLambda(lpx, lpy, dp) || !q.getIntervalLambda(lqx, lqy, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dqp  = dq * dp;
	IT pxq  = lpx * dqp;
	IT pyq  = lpy * dqp;
	IT rxq  = rx * dq;
	IT ryq  = ry * dq;
	IT lqxd = lqx * dp;
	IT lqyd = lqy * dp;
	IT lx   = pxq - lqxd;
	IT ly   = pyq - lqyd;
	IT gx   = rxq - lqx;
	IT gy   = ryq - lqy;
	IT dx   = lx * gx;
	IT dy   = ly * gy;
	IT d    = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEI_exact(const GenericPoint2T<IT, ET> &p,
                                const GenericPoint2T<IT, ET> &q, ET rx, ET ry)
{
	ET lpx, lpy, dp, lqx, lqy, dq;
	p.getExactLambda(lpx, lpy, dp);
	q.getExactLambda(lqx, lqy, dq);
	ET dqp  = dq * dp;
	ET pxq  = lpx * dqp;
	ET pyq  = lpy * dqp;
	ET rxq  = rx * dq;
	ET ryq  = ry * dq;
	ET lqxd = lqx * dp;
	ET lqyd = lqy * dp;
	ET lx   = pxq - lqxd;
	ET ly   = pyq - lqyd;
	ET gx   = rxq - lqx;
	ET gy   = ryq - lqy;
	ET dx   = lx * gx;
	ET dy   = ly * gy;
	ET d    = dx + dy;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEI_expansion(const GenericPoint2T<IT, ET> &p,
                                    const GenericPoint2T<IT, ET> &q, double rx,
                                    double ry)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[64], *lpx = lpx_p, lpy_p[64], *lpy = lpy_p, dp_p[64], *dp = dp_p,
	                  lqx_p[64], *lqx = lqx_p, lqy_p[64], *lqy = lqy_p, dq_p[64],
	                  *dq = dq_p;
	int lpx_len = 64, lpy_len = 64, dp_len = 64, lqx_len = 64, lqy_len = 64,
	    dq_len = 64;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &dp, dp_len);
	q.getExpansionLambda(&lqx, lqx_len, &lqy, lqy_len, &dq, dq_len);
	if ((dp[dp_len - 1] != 0) && (dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          dqp_p[64], *dqp = dqp_p;
		int dqp_len = o.Gen_Product_With_PreAlloc(dq_len, dq, dp_len, dp, &dqp, 64);
		double pxq_p[64], *pxq = pxq_p;
		int    pxq_len =
		  o.Gen_Product_With_PreAlloc(lpx_len, lpx, dqp_len, dqp, &pxq, 64);
		double pyq_p[64], *pyq = pyq_p;
		int    pyq_len =
		  o.Gen_Product_With_PreAlloc(lpy_len, lpy, dqp_len, dqp, &pyq, 64);
		double rxq_p[64], *rxq = rxq_p;
		int    rxq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, rx, &rxq, 64);
		double ryq_p[64], *ryq = ryq_p;
		int    ryq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, ry, &ryq, 64);
		double lqxd_p[64], *lqxd = lqxd_p;
		int    lqxd_len =
		  o.Gen_Product_With_PreAlloc(lqx_len, lqx, dp_len, dp, &lqxd, 64);
		double lqyd_p[64], *lqyd = lqyd_p;
		int    lqyd_len =
		  o.Gen_Product_With_PreAlloc(lqy_len, lqy, dp_len, dp, &lqyd, 64);
		double lx_p[64], *lx = lx_p;
		int    lx_len =
		  o.Gen_Diff_With_PreAlloc(pxq_len, pxq, lqxd_len, lqxd, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int    ly_len =
		  o.Gen_Diff_With_PreAlloc(pyq_len, pyq, lqyd_len, lqyd, &ly, 64);
		double gx_p[64], *gx = gx_p;
		int gx_len = o.Gen_Diff_With_PreAlloc(rxq_len, rxq, lqx_len, lqx, &gx, 64);
		double gy_p[64], *gy = gy_p;
		int gy_len = o.Gen_Diff_With_PreAlloc(ryq_len, ryq, lqy_len, lqy, &gy, 64);
		double dx_p[64], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 64);
		double dy_p[64], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 64);
		double d_p[64], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d, 64);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (lqyd_p != lqyd)
			FreeDoubles(lqyd);
		if (lqxd_p != lqxd)
			FreeDoubles(lqxd);
		if (ryq_p != ryq)
			FreeDoubles(ryq);
		if (rxq_p != rxq)
			FreeDoubles(rxq);
		if (pyq_p != pyq)
			FreeDoubles(pyq);
		if (pxq_p != pxq)
			FreeDoubles(pxq);
		if (dqp_p != dqp)
			FreeDoubles(dqp);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lqx_p != lqx)
			FreeDoubles(lqx);
		if (lqy_p != lqy)
			FreeDoubles(lqy);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign2D_IEI_exact<IT, ET>(p, q, rx, ry);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEI(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &q, double rx, double ry)
{
	Sign ret;
	ret = dotProductSign2D_IEI_interval<IT, ET>(p, q, rx, ry);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign2D_IEI_expansion<IT, ET>(p, q, rx, ry);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IEI(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &q,
                          const GenericPoint2T<IT, ET> &r)
{
	return dotProductSign2D_IEI<IT, ET>(p, q, r.x(), r.y());
}

template <typename IT, typename ET>
Sign dotProductSign2D_IIE_interval(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &r, IT qx,
                                   IT qy)
{
	IT lpx, lpy, dp, lrx, lry, dr;
	if (!p.getIntervalLambda(lpx, lpy, dp) || !r.getIntervalLambda(lrx, lry, dr))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT qxd = qx * dp;
	IT qyd = qy * dp;
	IT lx  = lpx - qxd;
	IT ly  = lpy - qyd;
	IT qxr = qx * dr;
	IT qyr = qy * dr;
	IT gx  = lrx - qxr;
	IT gy  = lry - qyr;
	IT dx  = lx * gx;
	IT dy  = ly * gy;
	IT d   = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IIE_exact(const GenericPoint2T<IT, ET> &p,
                                const GenericPoint2T<IT, ET> &r, ET qx, ET qy)
{
	ET lpx, lpy, dp, lrx, lry, dr;
	p.getExactLambda(lpx, lpy, dp);
	r.getExactLambda(lrx, lry, dr);
	ET qxd = qx * dp;
	ET qyd = qy * dp;
	ET lx  = lpx - qxd;
	ET ly  = lpy - qyd;
	ET qxr = qx * dr;
	ET qyr = qy * dr;
	ET gx  = lrx - qxr;
	ET gy  = lry - qyr;
	ET dx  = lx * gx;
	ET dy  = ly * gy;
	ET d   = dx + dy;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IIE_expansion(const GenericPoint2T<IT, ET> &p,
                                    const GenericPoint2T<IT, ET> &r, double qx,
                                    double qy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[64], *lpx = lpx_p, lpy_p[64], *lpy = lpy_p, dp_p[64], *dp = dp_p,
	                  lrx_p[64], *lrx = lrx_p, lry_p[64], *lry = lry_p, dr_p[64],
	                  *dr = dr_p;
	int lpx_len = 64, lpy_len = 64, dp_len = 64, lrx_len = 64, lry_len = 64,
	    dr_len = 64;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &dp, dp_len);
	r.getExpansionLambda(&lrx, lrx_len, &lry, lry_len, &dr, dr_len);
	if ((dp[dp_len - 1] != 0) && (dr[dr_len - 1] != 0))
	{
		expansionObject o;
		double          qxd_p[64], *qxd = qxd_p;
		int    qxd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qx, &qxd, 64);
		double qyd_p[64], *qyd = qyd_p;
		int    qyd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qy, &qyd, 64);
		double lx_p[64], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(lpx_len, lpx, qxd_len, qxd, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(lpy_len, lpy, qyd_len, qyd, &ly, 64);
		double qxr_p[64], *qxr = qxr_p;
		int    qxr_len = o.Gen_Scale_With_PreAlloc(dr_len, dr, qx, &qxr, 64);
		double qyr_p[64], *qyr = qyr_p;
		int    qyr_len = o.Gen_Scale_With_PreAlloc(dr_len, dr, qy, &qyr, 64);
		double gx_p[64], *gx = gx_p;
		int gx_len = o.Gen_Diff_With_PreAlloc(lrx_len, lrx, qxr_len, qxr, &gx, 64);
		double gy_p[64], *gy = gy_p;
		int gy_len = o.Gen_Diff_With_PreAlloc(lry_len, lry, qyr_len, qyr, &gy, 64);
		double dx_p[64], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 64);
		double dy_p[64], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 64);
		double d_p[64], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d, 64);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (qyr_p != qyr)
			FreeDoubles(qyr);
		if (qxr_p != qxr)
			FreeDoubles(qxr);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (qyd_p != qyd)
			FreeDoubles(qyd);
		if (qxd_p != qxd)
			FreeDoubles(qxd);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lrx_p != lrx)
			FreeDoubles(lrx);
		if (lry_p != lry)
			FreeDoubles(lry);
		if (dr_p != dr)
			FreeDoubles(dr);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign2D_IIE_exact<IT, ET>(p, r, qx, qy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign2D_IIE(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r, double qx, double qy)
{
	Sign ret;
	ret = dotProductSign2D_IIE_interval<IT, ET>(p, r, qx, qy);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign2D_IIE_expansion<IT, ET>(p, r, qx, qy);
}

template <typename IT, typename ET>
Sign dotProductSign2D_IIE(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r,
                          const GenericPoint2T<IT, ET> &q)
{
	return dotProductSign2D_IIE<IT, ET>(p, r, q.x(), q.y());
}

template <typename IT, typename ET>
Sign dotProductSign2D_III_interval(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &r,
                                   const GenericPoint2T<IT, ET> &q)
{
	IT lpx, lpy, dp, lrx, lry, dr, lqx, lqy, dq;
	if (!p.getIntervalLambda(lpx, lpy, dp) ||
	    !r.getIntervalLambda(lrx, lry, dr) || !q.getIntervalLambda(lqx, lqy, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT qxd  = lqx * dp;
	IT qyd  = lqy * dp;
	IT lpxq = lpx * dq;
	IT lpyq = lpy * dq;
	IT lx   = lpxq - qxd;
	IT ly   = lpyq - qyd;
	IT qxr  = lqx * dr;
	IT qyr  = lqy * dr;
	IT lrxq = lrx * dq;
	IT lryq = lry * dq;
	IT gx   = lrxq - qxr;
	IT gy   = lryq - qyr;
	IT dx   = lx * gx;
	IT dy   = ly * gy;
	IT d    = dx + dy;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_III_exact(const GenericPoint2T<IT, ET> &p,
                                const GenericPoint2T<IT, ET> &r,
                                const GenericPoint2T<IT, ET> &q)
{
	ET lpx, lpy, dp, lrx, lry, dr, lqx, lqy, dq;
	p.getExactLambda(lpx, lpy, dp);
	r.getExactLambda(lrx, lry, dr);
	q.getExactLambda(lqx, lqy, dq);
	ET qxd  = lqx * dp;
	ET qyd  = lqy * dp;
	ET lpxq = lpx * dq;
	ET lpyq = lpy * dq;
	ET lx   = lpxq - qxd;
	ET ly   = lpyq - qyd;
	ET qxr  = lqx * dr;
	ET qyr  = lqy * dr;
	ET lrxq = lrx * dq;
	ET lryq = lry * dq;
	ET gx   = lrxq - qxr;
	ET gy   = lryq - qyr;
	ET dx   = lx * gx;
	ET dy   = ly * gy;
	ET d    = dx + dy;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign2D_III_expansion(const GenericPoint2T<IT, ET> &p,
                                    const GenericPoint2T<IT, ET> &r,
                                    const GenericPoint2T<IT, ET> &q)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[64], *lpx = lpx_p, lpy_p[64], *lpy = lpy_p, dp_p[64], *dp = dp_p,
	                  lrx_p[64], *lrx = lrx_p, lry_p[64], *lry = lry_p, dr_p[64],
	                  *dr = dr_p, lqx_p[64], *lqx = lqx_p, lqy_p[64],
	                  *lqy = lqy_p, dq_p[64], *dq = dq_p;
	int lpx_len = 64, lpy_len = 64, dp_len = 64, lrx_len = 64, lry_len = 64,
	    dr_len = 64, lqx_len = 64, lqy_len = 64, dq_len = 64;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &dp, dp_len);
	r.getExpansionLambda(&lrx, lrx_len, &lry, lry_len, &dr, dr_len);
	q.getExpansionLambda(&lqx, lqx_len, &lqy, lqy_len, &dq, dq_len);
	if ((dp[dp_len - 1] != 0) && (dr[dr_len - 1] != 0) && (dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          qxd_p[64], *qxd = qxd_p;
		int             qxd_len =
		  o.Gen_Product_With_PreAlloc(lqx_len, lqx, dp_len, dp, &qxd, 64);
		double qyd_p[64], *qyd = qyd_p;
		int    qyd_len =
		  o.Gen_Product_With_PreAlloc(lqy_len, lqy, dp_len, dp, &qyd, 64);
		double lpxq_p[64], *lpxq = lpxq_p;
		int    lpxq_len =
		  o.Gen_Product_With_PreAlloc(lpx_len, lpx, dq_len, dq, &lpxq, 64);
		double lpyq_p[64], *lpyq = lpyq_p;
		int    lpyq_len =
		  o.Gen_Product_With_PreAlloc(lpy_len, lpy, dq_len, dq, &lpyq, 64);
		double lx_p[64], *lx = lx_p;
		int    lx_len =
		  o.Gen_Diff_With_PreAlloc(lpxq_len, lpxq, qxd_len, qxd, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int    ly_len =
		  o.Gen_Diff_With_PreAlloc(lpyq_len, lpyq, qyd_len, qyd, &ly, 64);
		double qxr_p[64], *qxr = qxr_p;
		int    qxr_len =
		  o.Gen_Product_With_PreAlloc(lqx_len, lqx, dr_len, dr, &qxr, 64);
		double qyr_p[64], *qyr = qyr_p;
		int    qyr_len =
		  o.Gen_Product_With_PreAlloc(lqy_len, lqy, dr_len, dr, &qyr, 64);
		double lrxq_p[64], *lrxq = lrxq_p;
		int    lrxq_len =
		  o.Gen_Product_With_PreAlloc(lrx_len, lrx, dq_len, dq, &lrxq, 64);
		double lryq_p[64], *lryq = lryq_p;
		int    lryq_len =
		  o.Gen_Product_With_PreAlloc(lry_len, lry, dq_len, dq, &lryq, 64);
		double gx_p[64], *gx = gx_p;
		int    gx_len =
		  o.Gen_Diff_With_PreAlloc(lrxq_len, lrxq, qxr_len, qxr, &gx, 64);
		double gy_p[64], *gy = gy_p;
		int    gy_len =
		  o.Gen_Diff_With_PreAlloc(lryq_len, lryq, qyr_len, qyr, &gy, 64);
		double dx_p[64], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 64);
		double dy_p[64], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 64);
		double d_p[64], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d, 64);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (lryq_p != lryq)
			FreeDoubles(lryq);
		if (lrxq_p != lrxq)
			FreeDoubles(lrxq);
		if (qyr_p != qyr)
			FreeDoubles(qyr);
		if (qxr_p != qxr)
			FreeDoubles(qxr);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (lpyq_p != lpyq)
			FreeDoubles(lpyq);
		if (lpxq_p != lpxq)
			FreeDoubles(lpxq);
		if (qyd_p != qyd)
			FreeDoubles(qyd);
		if (qxd_p != qxd)
			FreeDoubles(qxd);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lrx_p != lrx)
			FreeDoubles(lrx);
		if (lry_p != lry)
			FreeDoubles(lry);
		if (dr_p != dr)
			FreeDoubles(dr);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (lqx_p != lqx)
			FreeDoubles(lqx);
		if (lqy_p != lqy)
			FreeDoubles(lqy);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign2D_III_exact<IT, ET>(p, r, q);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign2D_III(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r,
                          const GenericPoint2T<IT, ET> &q)
{
	Sign ret;
	ret = dotProductSign2D_III_interval<IT, ET>(p, r, q);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign2D_III_expansion<IT, ET>(p, r, q);
}

template <typename IT, typename ET>
Sign dotProductSign3D_EEI_interval(const GenericPoint3T<IT, ET> &q, IT px,
                                   IT py, IT pz, IT rx, IT ry, IT rz)
{
	IT lqx, lqy, lqz, dq;
	if (!q.getIntervalLambda(lqx, lqy, lqz, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pxq = px * dq;
	IT pyq = py * dq;
	IT pzq = pz * dq;
	IT rxq = rx * dq;
	IT ryq = ry * dq;
	IT rzq = rz * dq;
	IT lx  = pxq - lqx;
	IT ly  = pyq - lqy;
	IT lz  = pzq - lqz;
	IT gx  = rxq - lqx;
	IT gy  = ryq - lqy;
	IT gz  = rzq - lqz;
	IT dx  = lx * gx;
	IT dy  = ly * gy;
	IT dz  = lz * gz;
	IT d1  = dx + dy;
	IT d   = d1 + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_EEI_exact(const GenericPoint3T<IT, ET> &q, ET px, ET py,
                                ET pz, ET rx, ET ry, ET rz)
{
	ET lqx, lqy, lqz, dq;
	q.getExactLambda(lqx, lqy, lqz, dq);
	ET pxq = px * dq;
	ET pyq = py * dq;
	ET pzq = pz * dq;
	ET rxq = rx * dq;
	ET ryq = ry * dq;
	ET rzq = rz * dq;
	ET lx  = pxq - lqx;
	ET ly  = pyq - lqy;
	ET lz  = pzq - lqz;
	ET gx  = rxq - lqx;
	ET gy  = ryq - lqy;
	ET gz  = rzq - lqz;
	ET dx  = lx * gx;
	ET dy  = ly * gy;
	ET dz  = lz * gz;
	ET d1  = dx + dy;
	ET d   = d1 + dz;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_EEI_expansion(const GenericPoint3T<IT, ET> &q, double px,
                                    double py, double pz, double rx, double ry,
                                    double rz)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lqx_p[64], *lqx = lqx_p, lqy_p[64], *lqy = lqy_p, lqz_p[64],
	                  *lqz = lqz_p, dq_p[64], *dq = dq_p;
	int lqx_len = 64, lqy_len = 64, lqz_len = 64, dq_len = 64;
	q.getExpansionLambda(&lqx, lqx_len, &lqy, lqy_len, &lqz, lqz_len, &dq,
	                     dq_len);
	if ((dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          pxq_p[64], *pxq = pxq_p;
		int    pxq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, px, &pxq, 64);
		double pyq_p[64], *pyq = pyq_p;
		int    pyq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, py, &pyq, 64);
		double pzq_p[64], *pzq = pzq_p;
		int    pzq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, pz, &pzq, 64);
		double rxq_p[64], *rxq = rxq_p;
		int    rxq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, rx, &rxq, 64);
		double ryq_p[64], *ryq = ryq_p;
		int    ryq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, ry, &ryq, 64);
		double rzq_p[64], *rzq = rzq_p;
		int    rzq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, rz, &rzq, 64);
		double lx_p[64], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(pxq_len, pxq, lqx_len, lqx, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(pyq_len, pyq, lqy_len, lqy, &ly, 64);
		double lz_p[64], *lz = lz_p;
		int lz_len = o.Gen_Diff_With_PreAlloc(pzq_len, pzq, lqz_len, lqz, &lz, 64);
		double gx_p[64], *gx = gx_p;
		int gx_len = o.Gen_Diff_With_PreAlloc(rxq_len, rxq, lqx_len, lqx, &gx, 64);
		double gy_p[64], *gy = gy_p;
		int gy_len = o.Gen_Diff_With_PreAlloc(ryq_len, ryq, lqy_len, lqy, &gy, 64);
		double gz_p[64], *gz = gz_p;
		int gz_len = o.Gen_Diff_With_PreAlloc(rzq_len, rzq, lqz_len, lqz, &gz, 64);
		double dx_p[64], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 64);
		double dy_p[64], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 64);
		double dz_p[64], *dz = dz_p;
		int dz_len = o.Gen_Product_With_PreAlloc(lz_len, lz, gz_len, gz, &dz, 64);
		double d1_p[64], *d1 = d1_p;
		int    d1_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d1, 64);
		double d_p[64], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(d1_len, d1, dz_len, dz, &d, 64);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (dz_p != dz)
			FreeDoubles(dz);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gz_p != gz)
			FreeDoubles(gz);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (lz_p != lz)
			FreeDoubles(lz);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (rzq_p != rzq)
			FreeDoubles(rzq);
		if (ryq_p != ryq)
			FreeDoubles(ryq);
		if (rxq_p != rxq)
			FreeDoubles(rxq);
		if (pzq_p != pzq)
			FreeDoubles(pzq);
		if (pyq_p != pyq)
			FreeDoubles(pyq);
		if (pxq_p != pxq)
			FreeDoubles(pxq);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lqx_p != lqx)
			FreeDoubles(lqx);
		if (lqy_p != lqy)
			FreeDoubles(lqy);
		if (lqz_p != lqz)
			FreeDoubles(lqz);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign3D_EEI_exact<IT, ET>(q, px, py, pz, rx, ry, rz);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign3D_EEI(const GenericPoint3T<IT, ET> &q, double px, double py,
                          double pz, double rx, double ry, double rz)
{
	Sign ret;
	ret = dotProductSign3D_EEI_interval<IT, ET>(q, px, py, pz, rx, ry, rz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign3D_EEI_expansion<IT, ET>(q, px, py, pz, rx, ry, rz);
}

template <typename IT, typename ET>
Sign dotProductSign3D_EEI(const GenericPoint3T<IT, ET> &q,
                          const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r)
{
	return dotProductSign3D_EEI<IT, ET>(q, p.x(), p.y(), p.z(), r.x(), r.y(),
	                                    r.z());
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEE_interval(const GenericPoint3T<IT, ET> &p, IT rx,
                                   IT ry, IT rz, IT qx, IT qy, IT qz)
{
	IT lpx, lpy, lpz, dp;
	if (!p.getIntervalLambda(lpx, lpy, lpz, dp))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT qxd = qx * dp;
	IT qyd = qy * dp;
	IT qzd = qz * dp;
	IT lx  = lpx - qxd;
	IT ly  = lpy - qyd;
	IT lz  = lpz - qzd;
	IT gx  = rx - qx;
	IT gy  = ry - qy;
	IT gz  = rz - qz;
	IT dx  = lx * gx;
	IT dy  = ly * gy;
	IT dz  = lz * gz;
	IT d1  = dx + dy;
	IT d   = d1 + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEE_exact(const GenericPoint3T<IT, ET> &p, ET rx, ET ry,
                                ET rz, ET qx, ET qy, ET qz)
{
	ET lpx, lpy, lpz, dp;
	p.getExactLambda(lpx, lpy, lpz, dp);
	ET qxd = qx * dp;
	ET qyd = qy * dp;
	ET qzd = qz * dp;
	ET lx  = lpx - qxd;
	ET ly  = lpy - qyd;
	ET lz  = lpz - qzd;
	ET gx  = rx - qx;
	ET gy  = ry - qy;
	ET gz  = rz - qz;
	ET dx  = lx * gx;
	ET dy  = ly * gy;
	ET dz  = lz * gz;
	ET d1  = dx + dy;
	ET d   = d1 + dz;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEE_expansion(const GenericPoint3T<IT, ET> &p, double rx,
                                    double ry, double rz, double qx, double qy,
                                    double qz)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[128], *lpx = lpx_p, lpy_p[128], *lpy = lpy_p, lpz_p[128],
	                   *lpz = lpz_p, dp_p[128], *dp = dp_p;
	int lpx_len = 128, lpy_len = 128, lpz_len = 128, dp_len = 128;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &lpz, lpz_len, &dp,
	                     dp_len);
	if ((dp[dp_len - 1] != 0))
	{
		expansionObject o;
		double          qxd_p[128], *qxd = qxd_p;
		int    qxd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qx, &qxd, 128);
		double qyd_p[128], *qyd = qyd_p;
		int    qyd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qy, &qyd, 128);
		double qzd_p[128], *qzd = qzd_p;
		int    qzd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qz, &qzd, 128);
		double lx_p[128], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(lpx_len, lpx, qxd_len, qxd, &lx, 128);
		double ly_p[128], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(lpy_len, lpy, qyd_len, qyd, &ly, 128);
		double lz_p[128], *lz = lz_p;
		int lz_len = o.Gen_Diff_With_PreAlloc(lpz_len, lpz, qzd_len, qzd, &lz, 128);
		double gx[2];
		o.two_Diff(rx, qx, gx);
		double gy[2];
		o.two_Diff(ry, qy, gy);
		double gz[2];
		o.two_Diff(rz, qz, gz);
		double dx_p[128], *dx = dx_p;
		int    dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, 2, gx, &dx, 128);
		double dy_p[128], *dy = dy_p;
		int    dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, 2, gy, &dy, 128);
		double dz_p[128], *dz = dz_p;
		int    dz_len = o.Gen_Product_With_PreAlloc(lz_len, lz, 2, gz, &dz, 128);
		double d1_p[128], *d1 = d1_p;
		int    d1_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d1, 128);
		double d_p[128], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(d1_len, d1, dz_len, dz, &d, 128);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (dz_p != dz)
			FreeDoubles(dz);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (lz_p != lz)
			FreeDoubles(lz);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (qzd_p != qzd)
			FreeDoubles(qzd);
		if (qyd_p != qyd)
			FreeDoubles(qyd);
		if (qxd_p != qxd)
			FreeDoubles(qxd);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (lpz_p != lpz)
			FreeDoubles(lpz);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign3D_IEE_exact<IT, ET>(p, rx, ry, rz, qx, qy, qz);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEE(const GenericPoint3T<IT, ET> &p, double rx, double ry,
                          double rz, double qx, double qy, double qz)
{
	Sign ret;
	ret = dotProductSign3D_IEE_interval<IT, ET>(p, rx, ry, rz, qx, qy, qz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign3D_IEE_expansion<IT, ET>(p, rx, ry, rz, qx, qy, qz);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEE(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r,
                          const GenericPoint3T<IT, ET> &q)
{
	return dotProductSign3D_IEE<IT, ET>(p, r.x(), r.y(), r.z(), q.x(), q.y(),
	                                    q.z());
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEI_interval(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &q, IT rx,
                                   IT ry, IT rz)
{
	IT lpx, lpy, lpz, dp, lqx, lqy, lqz, dq;
	if (!p.getIntervalLambda(lpx, lpy, lpz, dp) ||
	    !q.getIntervalLambda(lqx, lqy, lqz, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dqp  = dq * dp;
	IT pxq  = lpx * dqp;
	IT pyq  = lpy * dqp;
	IT pzq  = lpz * dqp;
	IT rxq  = rx * dq;
	IT ryq  = ry * dq;
	IT rzq  = rz * dq;
	IT lqxd = lqx * dp;
	IT lqyd = lqy * dp;
	IT lqzd = lqz * dp;
	IT lx   = pxq - lqxd;
	IT ly   = pyq - lqyd;
	IT lz   = pzq - lqzd;
	IT gx   = rxq - lqx;
	IT gy   = ryq - lqy;
	IT gz   = rzq - lqz;
	IT dx   = lx * gx;
	IT dy   = ly * gy;
	IT dz   = lz * gz;
	IT d1   = dx + dy;
	IT d    = d1 + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEI_exact(const GenericPoint3T<IT, ET> &p,
                                const GenericPoint3T<IT, ET> &q, ET rx, ET ry,
                                ET rz)
{
	ET lpx, lpy, lpz, dp, lqx, lqy, lqz, dq;
	p.getExactLambda(lpx, lpy, lpz, dp);
	q.getExactLambda(lqx, lqy, lqz, dq);
	ET dqp  = dq * dp;
	ET pxq  = lpx * dqp;
	ET pyq  = lpy * dqp;
	ET pzq  = lpz * dqp;
	ET rxq  = rx * dq;
	ET ryq  = ry * dq;
	ET rzq  = rz * dq;
	ET lqxd = lqx * dp;
	ET lqyd = lqy * dp;
	ET lqzd = lqz * dp;
	ET lx   = pxq - lqxd;
	ET ly   = pyq - lqyd;
	ET lz   = pzq - lqzd;
	ET gx   = rxq - lqx;
	ET gy   = ryq - lqy;
	ET gz   = rzq - lqz;
	ET dx   = lx * gx;
	ET dy   = ly * gy;
	ET dz   = lz * gz;
	ET d1   = dx + dy;
	ET d    = d1 + dz;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEI_expansion(const GenericPoint3T<IT, ET> &p,
                                    const GenericPoint3T<IT, ET> &q, double rx,
                                    double ry, double rz)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[64], *lpx = lpx_p, lpy_p[64], *lpy = lpy_p, lpz_p[64],
	                  *lpz = lpz_p, dp_p[64], *dp = dp_p, lqx_p[64], *lqx = lqx_p,
	                  lqy_p[64], *lqy = lqy_p, lqz_p[64], *lqz = lqz_p, dq_p[64],
	                  *dq = dq_p;
	int lpx_len = 64, lpy_len = 64, lpz_len = 64, dp_len = 64, lqx_len = 64,
	    lqy_len = 64, lqz_len = 64, dq_len = 64;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &lpz, lpz_len, &dp,
	                     dp_len);
	q.getExpansionLambda(&lqx, lqx_len, &lqy, lqy_len, &lqz, lqz_len, &dq,
	                     dq_len);
	if ((dp[dp_len - 1] != 0) && (dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          dqp_p[64], *dqp = dqp_p;
		int dqp_len = o.Gen_Product_With_PreAlloc(dq_len, dq, dp_len, dp, &dqp, 64);
		double pxq_p[64], *pxq = pxq_p;
		int    pxq_len =
		  o.Gen_Product_With_PreAlloc(lpx_len, lpx, dqp_len, dqp, &pxq, 64);
		double pyq_p[64], *pyq = pyq_p;
		int    pyq_len =
		  o.Gen_Product_With_PreAlloc(lpy_len, lpy, dqp_len, dqp, &pyq, 64);
		double pzq_p[64], *pzq = pzq_p;
		int    pzq_len =
		  o.Gen_Product_With_PreAlloc(lpz_len, lpz, dqp_len, dqp, &pzq, 64);
		double rxq_p[64], *rxq = rxq_p;
		int    rxq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, rx, &rxq, 64);
		double ryq_p[64], *ryq = ryq_p;
		int    ryq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, ry, &ryq, 64);
		double rzq_p[64], *rzq = rzq_p;
		int    rzq_len = o.Gen_Scale_With_PreAlloc(dq_len, dq, rz, &rzq, 64);
		double lqxd_p[64], *lqxd = lqxd_p;
		int    lqxd_len =
		  o.Gen_Product_With_PreAlloc(lqx_len, lqx, dp_len, dp, &lqxd, 64);
		double lqyd_p[64], *lqyd = lqyd_p;
		int    lqyd_len =
		  o.Gen_Product_With_PreAlloc(lqy_len, lqy, dp_len, dp, &lqyd, 64);
		double lqzd_p[64], *lqzd = lqzd_p;
		int    lqzd_len =
		  o.Gen_Product_With_PreAlloc(lqz_len, lqz, dp_len, dp, &lqzd, 64);
		double lx_p[64], *lx = lx_p;
		int    lx_len =
		  o.Gen_Diff_With_PreAlloc(pxq_len, pxq, lqxd_len, lqxd, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int    ly_len =
		  o.Gen_Diff_With_PreAlloc(pyq_len, pyq, lqyd_len, lqyd, &ly, 64);
		double lz_p[64], *lz = lz_p;
		int    lz_len =
		  o.Gen_Diff_With_PreAlloc(pzq_len, pzq, lqzd_len, lqzd, &lz, 64);
		double gx_p[64], *gx = gx_p;
		int gx_len = o.Gen_Diff_With_PreAlloc(rxq_len, rxq, lqx_len, lqx, &gx, 64);
		double gy_p[64], *gy = gy_p;
		int gy_len = o.Gen_Diff_With_PreAlloc(ryq_len, ryq, lqy_len, lqy, &gy, 64);
		double gz_p[64], *gz = gz_p;
		int gz_len = o.Gen_Diff_With_PreAlloc(rzq_len, rzq, lqz_len, lqz, &gz, 64);
		double dx_p[64], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 64);
		double dy_p[64], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 64);
		double dz_p[64], *dz = dz_p;
		int dz_len = o.Gen_Product_With_PreAlloc(lz_len, lz, gz_len, gz, &dz, 64);
		double d1_p[64], *d1 = d1_p;
		int    d1_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d1, 64);
		double d_p[64], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(d1_len, d1, dz_len, dz, &d, 64);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (dz_p != dz)
			FreeDoubles(dz);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gz_p != gz)
			FreeDoubles(gz);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (lz_p != lz)
			FreeDoubles(lz);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (lqzd_p != lqzd)
			FreeDoubles(lqzd);
		if (lqyd_p != lqyd)
			FreeDoubles(lqyd);
		if (lqxd_p != lqxd)
			FreeDoubles(lqxd);
		if (rzq_p != rzq)
			FreeDoubles(rzq);
		if (ryq_p != ryq)
			FreeDoubles(ryq);
		if (rxq_p != rxq)
			FreeDoubles(rxq);
		if (pzq_p != pzq)
			FreeDoubles(pzq);
		if (pyq_p != pyq)
			FreeDoubles(pyq);
		if (pxq_p != pxq)
			FreeDoubles(pxq);
		if (dqp_p != dqp)
			FreeDoubles(dqp);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (lpz_p != lpz)
			FreeDoubles(lpz);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lqx_p != lqx)
			FreeDoubles(lqx);
		if (lqy_p != lqy)
			FreeDoubles(lqy);
		if (lqz_p != lqz)
			FreeDoubles(lqz);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign3D_IEI_exact<IT, ET>(p, q, rx, ry, rz);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEI(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &q, double rx, double ry,
                          double rz)
{
	Sign ret;
	ret = dotProductSign3D_IEI_interval<IT, ET>(p, q, rx, ry, rz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign3D_IEI_expansion<IT, ET>(p, q, rx, ry, rz);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IEI(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &q,
                          const GenericPoint3T<IT, ET> &r)
{
	return dotProductSign3D_IEI<IT, ET>(p, q, r.x(), r.y(), r.z());
}

template <typename IT, typename ET>
Sign dotProductSign3D_IIE_interval(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &r, IT qx,
                                   IT qy, IT qz)
{
	IT lpx, lpy, lpz, dp, lrx, lry, lrz, dr;
	if (!p.getIntervalLambda(lpx, lpy, lpz, dp) ||
	    !r.getIntervalLambda(lrx, lry, lrz, dr))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT qxd = qx * dp;
	IT qyd = qy * dp;
	IT qzd = qz * dp;
	IT lx  = lpx - qxd;
	IT ly  = lpy - qyd;
	IT lz  = lpz - qzd;
	IT qxr = qx * dr;
	IT qyr = qy * dr;
	IT qzr = qz * dr;
	IT gx  = lrx - qxr;
	IT gy  = lry - qyr;
	IT gz  = lrz - qzr;
	IT dx  = lx * gx;
	IT dy  = ly * gy;
	IT dz  = lz * gz;
	IT d1  = dx + dy;
	IT d   = d1 + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IIE_exact(const GenericPoint3T<IT, ET> &p,
                                const GenericPoint3T<IT, ET> &r, ET qx, ET qy,
                                ET qz)
{
	ET lpx, lpy, lpz, dp, lrx, lry, lrz, dr;
	p.getExactLambda(lpx, lpy, lpz, dp);
	r.getExactLambda(lrx, lry, lrz, dr);
	ET qxd = qx * dp;
	ET qyd = qy * dp;
	ET qzd = qz * dp;
	ET lx  = lpx - qxd;
	ET ly  = lpy - qyd;
	ET lz  = lpz - qzd;
	ET qxr = qx * dr;
	ET qyr = qy * dr;
	ET qzr = qz * dr;
	ET gx  = lrx - qxr;
	ET gy  = lry - qyr;
	ET gz  = lrz - qzr;
	ET dx  = lx * gx;
	ET dy  = ly * gy;
	ET dz  = lz * gz;
	ET d1  = dx + dy;
	ET d   = d1 + dz;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IIE_expansion(const GenericPoint3T<IT, ET> &p,
                                    const GenericPoint3T<IT, ET> &r, double qx,
                                    double qy, double qz)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[64], *lpx = lpx_p, lpy_p[64], *lpy = lpy_p, lpz_p[64],
	                  *lpz = lpz_p, dp_p[64], *dp = dp_p, lrx_p[64], *lrx = lrx_p,
	                  lry_p[64], *lry = lry_p, lrz_p[64], *lrz = lrz_p, dr_p[64],
	                  *dr = dr_p;
	int lpx_len = 64, lpy_len = 64, lpz_len = 64, dp_len = 64, lrx_len = 64,
	    lry_len = 64, lrz_len = 64, dr_len = 64;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &lpz, lpz_len, &dp,
	                     dp_len);
	r.getExpansionLambda(&lrx, lrx_len, &lry, lry_len, &lrz, lrz_len, &dr,
	                     dr_len);
	if ((dp[dp_len - 1] != 0) && (dr[dr_len - 1] != 0))
	{
		expansionObject o;
		double          qxd_p[64], *qxd = qxd_p;
		int    qxd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qx, &qxd, 64);
		double qyd_p[64], *qyd = qyd_p;
		int    qyd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qy, &qyd, 64);
		double qzd_p[64], *qzd = qzd_p;
		int    qzd_len = o.Gen_Scale_With_PreAlloc(dp_len, dp, qz, &qzd, 64);
		double lx_p[64], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(lpx_len, lpx, qxd_len, qxd, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(lpy_len, lpy, qyd_len, qyd, &ly, 64);
		double lz_p[64], *lz = lz_p;
		int lz_len = o.Gen_Diff_With_PreAlloc(lpz_len, lpz, qzd_len, qzd, &lz, 64);
		double qxr_p[64], *qxr = qxr_p;
		int    qxr_len = o.Gen_Scale_With_PreAlloc(dr_len, dr, qx, &qxr, 64);
		double qyr_p[64], *qyr = qyr_p;
		int    qyr_len = o.Gen_Scale_With_PreAlloc(dr_len, dr, qy, &qyr, 64);
		double qzr_p[64], *qzr = qzr_p;
		int    qzr_len = o.Gen_Scale_With_PreAlloc(dr_len, dr, qz, &qzr, 64);
		double gx_p[64], *gx = gx_p;
		int gx_len = o.Gen_Diff_With_PreAlloc(lrx_len, lrx, qxr_len, qxr, &gx, 64);
		double gy_p[64], *gy = gy_p;
		int gy_len = o.Gen_Diff_With_PreAlloc(lry_len, lry, qyr_len, qyr, &gy, 64);
		double gz_p[64], *gz = gz_p;
		int gz_len = o.Gen_Diff_With_PreAlloc(lrz_len, lrz, qzr_len, qzr, &gz, 64);
		double dx_p[64], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 64);
		double dy_p[64], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 64);
		double dz_p[64], *dz = dz_p;
		int dz_len = o.Gen_Product_With_PreAlloc(lz_len, lz, gz_len, gz, &dz, 64);
		double d1_p[64], *d1 = d1_p;
		int    d1_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d1, 64);
		double d_p[64], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(d1_len, d1, dz_len, dz, &d, 64);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (dz_p != dz)
			FreeDoubles(dz);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gz_p != gz)
			FreeDoubles(gz);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (qzr_p != qzr)
			FreeDoubles(qzr);
		if (qyr_p != qyr)
			FreeDoubles(qyr);
		if (qxr_p != qxr)
			FreeDoubles(qxr);
		if (lz_p != lz)
			FreeDoubles(lz);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (qzd_p != qzd)
			FreeDoubles(qzd);
		if (qyd_p != qyd)
			FreeDoubles(qyd);
		if (qxd_p != qxd)
			FreeDoubles(qxd);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (lpz_p != lpz)
			FreeDoubles(lpz);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lrx_p != lrx)
			FreeDoubles(lrx);
		if (lry_p != lry)
			FreeDoubles(lry);
		if (lrz_p != lrz)
			FreeDoubles(lrz);
		if (dr_p != dr)
			FreeDoubles(dr);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign3D_IIE_exact<IT, ET>(p, r, qx, qy, qz);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign3D_IIE(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r, double qx, double qy,
                          double qz)
{
	Sign ret;
	ret = dotProductSign3D_IIE_interval<IT, ET>(p, r, qx, qy, qz);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign3D_IIE_expansion<IT, ET>(p, r, qx, qy, qz);
}

template <typename IT, typename ET>
Sign dotProductSign3D_IIE(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r,
                          const GenericPoint3T<IT, ET> &q)
{
	return dotProductSign3D_IIE<IT, ET>(p, r, q.x(), q.y(), q.z());
}

template <typename IT, typename ET>
Sign dotProductSign3D_III_interval(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &r,
                                   const GenericPoint3T<IT, ET> &q)
{
	IT lpx, lpy, lpz, dp, lrx, lry, lrz, dr, lqx, lqy, lqz, dq;
	if (!p.getIntervalLambda(lpx, lpy, lpz, dp) ||
	    !r.getIntervalLambda(lrx, lry, lrz, dr) ||
	    !q.getIntervalLambda(lqx, lqy, lqz, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT qxd  = lqx * dp;
	IT qyd  = lqy * dp;
	IT qzd  = lqz * dp;
	IT lpxq = lpx * dq;
	IT lpyq = lpy * dq;
	IT lpzq = lpz * dq;
	IT lx   = lpxq - qxd;
	IT ly   = lpyq - qyd;
	IT lz   = lpzq - qzd;
	IT qxr  = lqx * dr;
	IT qyr  = lqy * dr;
	IT qzr  = lqz * dr;
	IT lrxq = lrx * dq;
	IT lryq = lry * dq;
	IT lrzq = lrz * dq;
	IT gx   = lrxq - qxr;
	IT gy   = lryq - qyr;
	IT gz   = lrzq - qzr;
	IT dx   = lx * gx;
	IT dy   = ly * gy;
	IT dz   = lz * gz;
	IT d1   = dx + dy;
	IT d    = d1 + dz;
	if (!d.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_III_exact(const GenericPoint3T<IT, ET> &p,
                                const GenericPoint3T<IT, ET> &r,
                                const GenericPoint3T<IT, ET> &q)
{
	ET lpx, lpy, lpz, dp, lrx, lry, lrz, dr, lqx, lqy, lqz, dq;
	p.getExactLambda(lpx, lpy, lpz, dp);
	r.getExactLambda(lrx, lry, lrz, dr);
	q.getExactLambda(lqx, lqy, lqz, dq);
	ET qxd  = lqx * dp;
	ET qyd  = lqy * dp;
	ET qzd  = lqz * dp;
	ET lpxq = lpx * dq;
	ET lpyq = lpy * dq;
	ET lpzq = lpz * dq;
	ET lx   = lpxq - qxd;
	ET ly   = lpyq - qyd;
	ET lz   = lpzq - qzd;
	ET qxr  = lqx * dr;
	ET qyr  = lqy * dr;
	ET qzr  = lqz * dr;
	ET lrxq = lrx * dq;
	ET lryq = lry * dq;
	ET lrzq = lrz * dq;
	ET gx   = lrxq - qxr;
	ET gy   = lryq - qyr;
	ET gz   = lrzq - qzr;
	ET dx   = lx * gx;
	ET dy   = ly * gy;
	ET dz   = lz * gz;
	ET d1   = dx + dy;
	ET d    = d1 + dz;
	return OMC::sign(d);
}

template <typename IT, typename ET>
Sign dotProductSign3D_III_expansion(const GenericPoint3T<IT, ET> &p,
                                    const GenericPoint3T<IT, ET> &r,
                                    const GenericPoint3T<IT, ET> &q)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double lpx_p[32], *lpx = lpx_p, lpy_p[32], *lpy = lpy_p, lpz_p[32],
	                  *lpz = lpz_p, dp_p[32], *dp = dp_p, lrx_p[32], *lrx = lrx_p,
	                  lry_p[32], *lry = lry_p, lrz_p[32], *lrz = lrz_p, dr_p[32],
	                  *dr = dr_p, lqx_p[32], *lqx = lqx_p, lqy_p[32],
	                  *lqy = lqy_p, lqz_p[32], *lqz = lqz_p, dq_p[32], *dq = dq_p;
	int lpx_len = 32, lpy_len = 32, lpz_len = 32, dp_len = 32, lrx_len = 32,
	    lry_len = 32, lrz_len = 32, dr_len = 32, lqx_len = 32, lqy_len = 32,
	    lqz_len = 32, dq_len = 32;
	p.getExpansionLambda(&lpx, lpx_len, &lpy, lpy_len, &lpz, lpz_len, &dp,
	                     dp_len);
	r.getExpansionLambda(&lrx, lrx_len, &lry, lry_len, &lrz, lrz_len, &dr,
	                     dr_len);
	q.getExpansionLambda(&lqx, lqx_len, &lqy, lqy_len, &lqz, lqz_len, &dq,
	                     dq_len);
	if ((dp[dp_len - 1] != 0) && (dr[dr_len - 1] != 0) && (dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          qxd_p[32], *qxd = qxd_p;
		int             qxd_len =
		  o.Gen_Product_With_PreAlloc(lqx_len, lqx, dp_len, dp, &qxd, 32);
		double qyd_p[32], *qyd = qyd_p;
		int    qyd_len =
		  o.Gen_Product_With_PreAlloc(lqy_len, lqy, dp_len, dp, &qyd, 32);
		double qzd_p[32], *qzd = qzd_p;
		int    qzd_len =
		  o.Gen_Product_With_PreAlloc(lqz_len, lqz, dp_len, dp, &qzd, 32);
		double lpxq_p[32], *lpxq = lpxq_p;
		int    lpxq_len =
		  o.Gen_Product_With_PreAlloc(lpx_len, lpx, dq_len, dq, &lpxq, 32);
		double lpyq_p[32], *lpyq = lpyq_p;
		int    lpyq_len =
		  o.Gen_Product_With_PreAlloc(lpy_len, lpy, dq_len, dq, &lpyq, 32);
		double lpzq_p[32], *lpzq = lpzq_p;
		int    lpzq_len =
		  o.Gen_Product_With_PreAlloc(lpz_len, lpz, dq_len, dq, &lpzq, 32);
		double lx_p[32], *lx = lx_p;
		int    lx_len =
		  o.Gen_Diff_With_PreAlloc(lpxq_len, lpxq, qxd_len, qxd, &lx, 32);
		double ly_p[32], *ly = ly_p;
		int    ly_len =
		  o.Gen_Diff_With_PreAlloc(lpyq_len, lpyq, qyd_len, qyd, &ly, 32);
		double lz_p[32], *lz = lz_p;
		int    lz_len =
		  o.Gen_Diff_With_PreAlloc(lpzq_len, lpzq, qzd_len, qzd, &lz, 32);
		double qxr_p[32], *qxr = qxr_p;
		int    qxr_len =
		  o.Gen_Product_With_PreAlloc(lqx_len, lqx, dr_len, dr, &qxr, 32);
		double qyr_p[32], *qyr = qyr_p;
		int    qyr_len =
		  o.Gen_Product_With_PreAlloc(lqy_len, lqy, dr_len, dr, &qyr, 32);
		double qzr_p[32], *qzr = qzr_p;
		int    qzr_len =
		  o.Gen_Product_With_PreAlloc(lqz_len, lqz, dr_len, dr, &qzr, 32);
		double lrxq_p[32], *lrxq = lrxq_p;
		int    lrxq_len =
		  o.Gen_Product_With_PreAlloc(lrx_len, lrx, dq_len, dq, &lrxq, 32);
		double lryq_p[32], *lryq = lryq_p;
		int    lryq_len =
		  o.Gen_Product_With_PreAlloc(lry_len, lry, dq_len, dq, &lryq, 32);
		double lrzq_p[32], *lrzq = lrzq_p;
		int    lrzq_len =
		  o.Gen_Product_With_PreAlloc(lrz_len, lrz, dq_len, dq, &lrzq, 32);
		double gx_p[32], *gx = gx_p;
		int    gx_len =
		  o.Gen_Diff_With_PreAlloc(lrxq_len, lrxq, qxr_len, qxr, &gx, 32);
		double gy_p[32], *gy = gy_p;
		int    gy_len =
		  o.Gen_Diff_With_PreAlloc(lryq_len, lryq, qyr_len, qyr, &gy, 32);
		double gz_p[32], *gz = gz_p;
		int    gz_len =
		  o.Gen_Diff_With_PreAlloc(lrzq_len, lrzq, qzr_len, qzr, &gz, 32);
		double dx_p[32], *dx = dx_p;
		int dx_len = o.Gen_Product_With_PreAlloc(lx_len, lx, gx_len, gx, &dx, 32);
		double dy_p[32], *dy = dy_p;
		int dy_len = o.Gen_Product_With_PreAlloc(ly_len, ly, gy_len, gy, &dy, 32);
		double dz_p[32], *dz = dz_p;
		int dz_len = o.Gen_Product_With_PreAlloc(lz_len, lz, gz_len, gz, &dz, 32);
		double d1_p[32], *d1 = d1_p;
		int    d1_len = o.Gen_Sum_With_PreAlloc(dx_len, dx, dy_len, dy, &d1, 32);
		double d_p[32], *d = d_p;
		int    d_len = o.Gen_Sum_With_PreAlloc(d1_len, d1, dz_len, dz, &d, 32);

		return_value = d[d_len - 1];
		if (d_p != d)
			FreeDoubles(d);
		if (d1_p != d1)
			FreeDoubles(d1);
		if (dz_p != dz)
			FreeDoubles(dz);
		if (dy_p != dy)
			FreeDoubles(dy);
		if (dx_p != dx)
			FreeDoubles(dx);
		if (gz_p != gz)
			FreeDoubles(gz);
		if (gy_p != gy)
			FreeDoubles(gy);
		if (gx_p != gx)
			FreeDoubles(gx);
		if (lrzq_p != lrzq)
			FreeDoubles(lrzq);
		if (lryq_p != lryq)
			FreeDoubles(lryq);
		if (lrxq_p != lrxq)
			FreeDoubles(lrxq);
		if (qzr_p != qzr)
			FreeDoubles(qzr);
		if (qyr_p != qyr)
			FreeDoubles(qyr);
		if (qxr_p != qxr)
			FreeDoubles(qxr);
		if (lz_p != lz)
			FreeDoubles(lz);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (lpzq_p != lpzq)
			FreeDoubles(lpzq);
		if (lpyq_p != lpyq)
			FreeDoubles(lpyq);
		if (lpxq_p != lpxq)
			FreeDoubles(lpxq);
		if (qzd_p != qzd)
			FreeDoubles(qzd);
		if (qyd_p != qyd)
			FreeDoubles(qyd);
		if (qxd_p != qxd)
			FreeDoubles(qxd);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lpx_p != lpx)
			FreeDoubles(lpx);
		if (lpy_p != lpy)
			FreeDoubles(lpy);
		if (lpz_p != lpz)
			FreeDoubles(lpz);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lrx_p != lrx)
			FreeDoubles(lrx);
		if (lry_p != lry)
			FreeDoubles(lry);
		if (lrz_p != lrz)
			FreeDoubles(lrz);
		if (dr_p != dr)
			FreeDoubles(dr);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (lqx_p != lqx)
			FreeDoubles(lqx);
		if (lqy_p != lqy)
			FreeDoubles(lqy);
		if (lqz_p != lqz)
			FreeDoubles(lqz);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return dotProductSign3D_III_exact<IT, ET>(p, r, q);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign dotProductSign3D_III(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r,
                          const GenericPoint3T<IT, ET> &q)
{
	Sign ret;
	ret = dotProductSign3D_III_interval<IT, ET>(p, r, q);
	if (is_sign_reliable(ret))
		return ret;
	return dotProductSign3D_III_expansion<IT, ET>(p, r, q);
}

template <typename IT, typename ET>
Sign inCirclexy_IEEE_interval(const GenericPoint3T<IT, ET> &p1, IT pbx, IT pby,
                              IT pcx, IT pcy, IT pdx, IT pdy)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pdxt   = pdx * d1;
	IT pdyt   = pdy * d1;
	IT adx    = l1x - pdxt;
	IT ady    = l1y - pdyt;
	IT bdx    = pbx - pdx;
	IT bdy    = pby - pdy;
	IT cdx    = pcx - pdx;
	IT cdy    = pcy - pdy;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT alift  = alifta + aliftb;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT blift  = blifta + bliftb;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT clift  = clifta + cliftb;
	IT la     = alift * bcdet;
	IT lbt    = blift * cadet;
	IT lb     = lbt * d1;
	IT lct    = clift * abdet;
	IT lc     = lct * d1;
	IT lab    = la + lb;
	IT L      = lab + lc;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IEEE_exact(const GenericPoint3T<IT, ET> &p1, ET pbx, ET pby,
                           ET pcx, ET pcy, ET pdx, ET pdy)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET pdxt   = pdx * d1;
	ET pdyt   = pdy * d1;
	ET adx    = l1x - pdxt;
	ET ady    = l1y - pdyt;
	ET bdx    = pbx - pdx;
	ET bdy    = pby - pdy;
	ET cdx    = pcx - pdx;
	ET cdy    = pcy - pdy;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET alift  = alifta + aliftb;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET blift  = blifta + bliftb;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET clift  = clifta + cliftb;
	ET la     = alift * bcdet;
	ET lbt    = blift * cadet;
	ET lb     = lbt * d1;
	ET lct    = clift * abdet;
	ET lc     = lct * d1;
	ET lab    = la + lb;
	ET L      = lab + lc;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IEEE_expansion(const GenericPoint3T<IT, ET> &p1, double pbx,
                               double pby, double pcx, double pcy, double pdx,
                               double pdy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          pdxt_p[64], *pdxt = pdxt_p;
		int    pdxt_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdxt, 64);
		double pdyt_p[64], *pdyt = pdyt_p;
		int    pdyt_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdyt, 64);
		double adx_p[64], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdxt_len, pdxt, &adx, 64);
		double ady_p[64], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdyt_len, pdyt, &ady, 64);
		double bdx[2];
		o.two_Diff(pbx, pdx, bdx);
		double bdy[2];
		o.two_Diff(pby, pdy, bdy);
		double cdx[2];
		o.two_Diff(pcx, pdx, cdx);
		double cdy[2];
		o.two_Diff(pcy, pdy, cdy);
		double abdeta_p[64], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, 2, bdy, &abdeta, 64);
		double abdetb_p[64], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(2, bdx, ady_len, ady, &abdetb, 64);
		double abdet_p[64], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 64);
		double bcdeta[8];
		int    bcdeta_len = o.Gen_Product(2, bdx, 2, cdy, bcdeta);
		double bcdetb[8];
		int    bcdetb_len = o.Gen_Product(2, cdx, 2, bdy, bcdetb);
		double bcdet[16];
		int bcdet_len = o.Gen_Diff(bcdeta_len, bcdeta, bcdetb_len, bcdetb, bcdet);
		double cadeta_p[64], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 64);
		double cadetb_p[64], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 64);
		double cadet_p[64], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 64);
		double alifta_p[64], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 64);
		double aliftb_p[64], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 64);
		double alift_p[64], *alift = alift_p;
		int    alift_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                           aliftb, &alift, 64);
		double blifta[8];
		int    blifta_len = o.Gen_Product(2, bdx, 2, bdx, blifta);
		double bliftb[8];
		int    bliftb_len = o.Gen_Product(2, bdy, 2, bdy, bliftb);
		double blift[16];
		int    blift_len = o.Gen_Sum(blifta_len, blifta, bliftb_len, bliftb, blift);
		double clifta[8];
		int    clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
		double cliftb[8];
		int    cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
		double clift[16];
		int    clift_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, clift);
		double la_p[64], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 64);
		double lbt_p[64], *lbt = lbt_p;
		int    lbt_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lbt, 64);
		double lb_p[64], *lb = lb_p;
		int lb_len = o.Gen_Product_With_PreAlloc(lbt_len, lbt, d1_len, d1, &lb, 64);
		double lct_p[64], *lct = lct_p;
		int    lct_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lct, 64);
		double lc_p[64], *lc = lc_p;
		int lc_len = o.Gen_Product_With_PreAlloc(lct_len, lct, d1_len, d1, &lc, 64);
		double lab_p[64], *lab = lab_p;
		int    lab_len = o.Gen_Sum_With_PreAlloc(la_len, la, lb_len, lb, &lab, 64);
		double L_p[64], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, lc_len, lc, &L, 64);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lct_p != lct)
			FreeDoubles(lct);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (lbt_p != lbt)
			FreeDoubles(lbt);
		if (la_p != la)
			FreeDoubles(la);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (pdyt_p != pdyt)
			FreeDoubles(pdyt);
		if (pdxt_p != pdxt)
			FreeDoubles(pdxt);
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
		return inCirclexy_IEEE_exact<IT, ET>(p1, pbx, pby, pcx, pcy, pdx, pdy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCirclexy_IEEE(const GenericPoint3T<IT, ET> &p1, double pbx, double pby,
                     double pcx, double pcy, double pdx, double pdy)
{
	Sign ret;
	ret = inCirclexy_IEEE_interval<IT, ET>(p1, pbx, pby, pcx, pcy, pdx, pdy);
	if (is_sign_reliable(ret))
		return ret;
	return inCirclexy_IEEE_expansion<IT, ET>(p1, pbx, pby, pcx, pcy, pdx, pdy);
}

template <typename IT, typename ET>
Sign inCirclexy_IEEE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &pb,
                     const GenericPoint3T<IT, ET> &pc,
                     const GenericPoint3T<IT, ET> &pd)
{
	return inCirclexy_IEEE<IT, ET>(p1, pb.x(), pb.y(), pc.x(), pc.y(), pd.x(),
	                               pd.y());
}

template <typename IT, typename ET>
Sign inCirclexy_IIEE_interval(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2, IT pcx, IT pcy,
                              IT pdx, IT pdy)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pdx1   = pdx * d1;
	IT pdy1   = pdy * d1;
	IT adx    = l1x - pdx1;
	IT ady    = l1y - pdy1;
	IT pdx2   = pdx * d2;
	IT pdy2   = pdy * d2;
	IT bdx    = l2x - pdx2;
	IT bdy    = l2y - pdy2;
	IT cdx    = pcx - pdx;
	IT cdy    = pcy - pdy;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT aliftt = alifta + aliftb;
	IT alift  = aliftt * d2;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT blift  = blifta + bliftb;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT cliftt = clifta + cliftb;
	IT clift  = cliftt * d2;
	IT la     = alift * bcdet;
	IT lb     = blift * cadet;
	IT lc     = clift * abdet;
	IT lab    = lc + lb;
	IT lab2   = lab * d1;
	IT L      = lab2 + la;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IIEE_exact(const GenericPoint3T<IT, ET> &p1,
                           const GenericPoint3T<IT, ET> &p2, ET pcx, ET pcy,
                           ET pdx, ET pdy)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET pdx1   = pdx * d1;
	ET pdy1   = pdy * d1;
	ET adx    = l1x - pdx1;
	ET ady    = l1y - pdy1;
	ET pdx2   = pdx * d2;
	ET pdy2   = pdy * d2;
	ET bdx    = l2x - pdx2;
	ET bdy    = l2y - pdy2;
	ET cdx    = pcx - pdx;
	ET cdy    = pcy - pdy;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET aliftt = alifta + aliftb;
	ET alift  = aliftt * d2;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET blift  = blifta + bliftb;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET cliftt = clifta + cliftb;
	ET clift  = cliftt * d2;
	ET la     = alift * bcdet;
	ET lb     = blift * cadet;
	ET lc     = clift * abdet;
	ET lab    = lc + lb;
	ET lab2   = lab * d1;
	ET L      = lab2 + la;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IIEE_expansion(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double pcx,
                               double pcy, double pdx, double pdy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32],
	                  *l1z = l1z_p, d1_p[32], *d1 = d1_p, l2x_p[32], *l2x = l2x_p,
	                  l2y_p[32], *l2y = l2y_p, l2z_p[32], *l2z = l2z_p, d2_p[32],
	                  *d2 = d2_p;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          pdx1_p[32], *pdx1 = pdx1_p;
		int    pdx1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdx1, 32);
		double pdy1_p[32], *pdy1 = pdy1_p;
		int    pdy1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdy1, 32);
		double adx_p[32], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
		double ady_p[32], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
		double pdx2_p[32], *pdx2 = pdx2_p;
		int    pdx2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdx, &pdx2, 32);
		double pdy2_p[32], *pdy2 = pdy2_p;
		int    pdy2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdy, &pdy2, 32);
		double bdx_p[32], *bdx = bdx_p;
		int    bdx_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
		double bdy_p[32], *bdy = bdy_p;
		int    bdy_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
		double cdx[2];
		o.two_Diff(pcx, pdx, cdx);
		double cdy[2];
		o.two_Diff(pcy, pdy, cdy);
		double abdeta_p[32], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
		double abdetb_p[32], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
		double abdet_p[32], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 32);
		double bcdeta_p[32], *bcdeta = bcdeta_p;
		int    bcdeta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, 2, cdy, &bcdeta, 32);
		double bcdetb_p[32], *bcdetb = bcdetb_p;
		int    bcdetb_len =
		  o.Gen_Product_With_PreAlloc(2, cdx, bdy_len, bdy, &bcdetb, 32);
		double bcdet_p[32], *bcdet = bcdet_p;
		int    bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len,
		                                            bcdetb, &bcdet, 32);
		double cadeta_p[32], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 32);
		double cadetb_p[32], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 32);
		double cadet_p[32], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 32);
		double alifta_p[32], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
		double aliftb_p[32], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
		double aliftt_p[32], *aliftt = aliftt_p;
		int    aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                            aliftb, &aliftt, 32);
		double alift_p[32], *alift = alift_p;
		int    alift_len =
		  o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift, 32);
		double blifta_p[32], *blifta = blifta_p;
		int    blifta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
		double bliftb_p[32], *bliftb = bliftb_p;
		int    bliftb_len =
		  o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
		double blift_p[32], *blift = blift_p;
		int    blift_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len,
		                                           bliftb, &blift, 32);
		double clifta[8];
		int    clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
		double cliftb[8];
		int    cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
		double cliftt[16];
		int cliftt_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, cliftt);
		double clift_p[32], *clift = clift_p;
		int    clift_len =
		  o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
		double la_p[32], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
		double lb_p[32], *lb = lb_p;
		int    lb_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
		double lc_p[32], *lc = lc_p;
		int    lc_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
		double lab_p[32], *lab = lab_p;
		int    lab_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab, 32);
		double lab2_p[32], *lab2 = lab2_p;
		int    lab2_len =
		  o.Gen_Product_With_PreAlloc(lab_len, lab, d1_len, d1, &lab2, 32);
		double L_p[32], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab2_len, lab2, la_len, la, &L, 32);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab2_p != lab2)
			FreeDoubles(lab2);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (la_p != la)
			FreeDoubles(la);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bliftb_p != bliftb)
			FreeDoubles(bliftb);
		if (blifta_p != blifta)
			FreeDoubles(blifta);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (aliftt_p != aliftt)
			FreeDoubles(aliftt);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (bcdet_p != bcdet)
			FreeDoubles(bcdet);
		if (bcdetb_p != bcdetb)
			FreeDoubles(bcdetb);
		if (bcdeta_p != bcdeta)
			FreeDoubles(bcdeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (bdy_p != bdy)
			FreeDoubles(bdy);
		if (bdx_p != bdx)
			FreeDoubles(bdx);
		if (pdy2_p != pdy2)
			FreeDoubles(pdy2);
		if (pdx2_p != pdx2)
			FreeDoubles(pdx2);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (pdy1_p != pdy1)
			FreeDoubles(pdy1);
		if (pdx1_p != pdx1)
			FreeDoubles(pdx1);
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
		return inCirclexy_IIEE_exact<IT, ET>(p1, p2, pcx, pcy, pdx, pdy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCirclexy_IIEE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2, double pcx, double pcy,
                     double pdx, double pdy)
{
	Sign ret;
	ret = inCirclexy_IIEE_interval<IT, ET>(p1, p2, pcx, pcy, pdx, pdy);
	if (is_sign_reliable(ret))
		return ret;
	return inCirclexy_IIEE_expansion<IT, ET>(p1, p2, pcx, pcy, pdx, pdy);
}

template <typename IT, typename ET>
Sign inCirclexy_IIEE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &pc,
                     const GenericPoint3T<IT, ET> &pd)
{
	return inCirclexy_IIEE<IT, ET>(p1, p2, pc.x(), pc.y(), pd.x(), pd.y());
}

template <typename IT, typename ET>
Sign inCirclexy_IIIE_interval(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3, IT pdx, IT pdy)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pdx1   = pdx * d1;
	IT pdy1   = pdy * d1;
	IT adx    = l1x - pdx1;
	IT ady    = l1y - pdy1;
	IT pdx2   = pdx * d2;
	IT pdy2   = pdy * d2;
	IT bdx    = l2x - pdx2;
	IT bdy    = l2y - pdy2;
	IT pdx3   = pdx * d3;
	IT pdy3   = pdy * d3;
	IT cdx    = l3x - pdx3;
	IT cdy    = l3y - pdy3;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT aliftt = alifta + aliftb;
	IT alift2 = aliftt * d2;
	IT alift  = alift2 * d3;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT bliftt = blifta + bliftb;
	IT blift  = bliftt * d3;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT cliftt = clifta + cliftb;
	IT clift  = cliftt * d2;
	IT la     = alift * bcdet;
	IT lb     = blift * cadet;
	IT lc     = clift * abdet;
	IT lab2   = lc + lb;
	IT lab    = lab2 * d1;
	IT L      = lab + la;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IIIE_exact(const GenericPoint3T<IT, ET> &p1,
                           const GenericPoint3T<IT, ET> &p2,
                           const GenericPoint3T<IT, ET> &p3, ET pdx, ET pdy)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	ET pdx1   = pdx * d1;
	ET pdy1   = pdy * d1;
	ET adx    = l1x - pdx1;
	ET ady    = l1y - pdy1;
	ET pdx2   = pdx * d2;
	ET pdy2   = pdy * d2;
	ET bdx    = l2x - pdx2;
	ET bdy    = l2y - pdy2;
	ET pdx3   = pdx * d3;
	ET pdy3   = pdy * d3;
	ET cdx    = l3x - pdx3;
	ET cdy    = l3y - pdy3;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET aliftt = alifta + aliftb;
	ET alift2 = aliftt * d2;
	ET alift  = alift2 * d3;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET bliftt = blifta + bliftb;
	ET blift  = bliftt * d3;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET cliftt = clifta + cliftb;
	ET clift  = cliftt * d2;
	ET la     = alift * bcdet;
	ET lb     = blift * cadet;
	ET lc     = clift * abdet;
	ET lab2   = lc + lb;
	ET lab    = lab2 * d1;
	ET L      = lab + la;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IIIE_expansion(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, double pdx,
                               double pdy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32],
	                  *l1z = l1z_p, d1_p[32], *d1 = d1_p, l2x_p[32], *l2x = l2x_p,
	                  l2y_p[32], *l2y = l2y_p, l2z_p[32], *l2z = l2z_p, d2_p[32],
	                  *d2 = d2_p, l3x_p[32], *l3x = l3x_p, l3y_p[32],
	                  *l3y = l3y_p, l3z_p[32], *l3z = l3z_p, d3_p[32], *d3 = d3_p;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32, l3x_len = 32, l3y_len = 32,
	    l3z_len = 32, d3_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          pdx1_p[32], *pdx1 = pdx1_p;
		int    pdx1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdx1, 32);
		double pdy1_p[32], *pdy1 = pdy1_p;
		int    pdy1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdy1, 32);
		double adx_p[32], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
		double ady_p[32], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
		double pdx2_p[32], *pdx2 = pdx2_p;
		int    pdx2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdx, &pdx2, 32);
		double pdy2_p[32], *pdy2 = pdy2_p;
		int    pdy2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdy, &pdy2, 32);
		double bdx_p[32], *bdx = bdx_p;
		int    bdx_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
		double bdy_p[32], *bdy = bdy_p;
		int    bdy_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
		double pdx3_p[32], *pdx3 = pdx3_p;
		int    pdx3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pdx, &pdx3, 32);
		double pdy3_p[32], *pdy3 = pdy3_p;
		int    pdy3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pdy, &pdy3, 32);
		double cdx_p[32], *cdx = cdx_p;
		int    cdx_len =
		  o.Gen_Diff_With_PreAlloc(l3x_len, l3x, pdx3_len, pdx3, &cdx, 32);
		double cdy_p[32], *cdy = cdy_p;
		int    cdy_len =
		  o.Gen_Diff_With_PreAlloc(l3y_len, l3y, pdy3_len, pdy3, &cdy, 32);
		double abdeta_p[32], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
		double abdetb_p[32], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
		double abdet_p[32], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 32);
		double bcdeta_p[32], *bcdeta = bcdeta_p;
		int    bcdeta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 32);
		double bcdetb_p[32], *bcdetb = bcdetb_p;
		int    bcdetb_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 32);
		double bcdet_p[32], *bcdet = bcdet_p;
		int    bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len,
		                                            bcdetb, &bcdet, 32);
		double cadeta_p[32], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 32);
		double cadetb_p[32], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 32);
		double cadet_p[32], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 32);
		double alifta_p[32], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
		double aliftb_p[32], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
		double aliftt_p[32], *aliftt = aliftt_p;
		int    aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                            aliftb, &aliftt, 32);
		double alift2_p[32], *alift2 = alift2_p;
		int    alift2_len =
		  o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 32);
		double alift_p[32], *alift = alift_p;
		int    alift_len =
		  o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 32);
		double blifta_p[32], *blifta = blifta_p;
		int    blifta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
		double bliftb_p[32], *bliftb = bliftb_p;
		int    bliftb_len =
		  o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
		double bliftt_p[32], *bliftt = bliftt_p;
		int    bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len,
		                                            bliftb, &bliftt, 32);
		double blift_p[32], *blift = blift_p;
		int    blift_len =
		  o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 32);
		double clifta_p[32], *clifta = clifta_p;
		int    clifta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 32);
		double cliftb_p[32], *cliftb = cliftb_p;
		int    cliftb_len =
		  o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 32);
		double cliftt_p[32], *cliftt = cliftt_p;
		int    cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len,
		                                            cliftb, &cliftt, 32);
		double clift_p[32], *clift = clift_p;
		int    clift_len =
		  o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
		double la_p[32], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
		double lb_p[32], *lb = lb_p;
		int    lb_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
		double lc_p[32], *lc = lc_p;
		int    lc_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
		double lab2_p[32], *lab2 = lab2_p;
		int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 32);
		double lab_p[32], *lab = lab_p;
		int    lab_len =
		  o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 32);
		double L_p[32], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 32);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lab2_p != lab2)
			FreeDoubles(lab2);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (la_p != la)
			FreeDoubles(la);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (cliftt_p != cliftt)
			FreeDoubles(cliftt);
		if (cliftb_p != cliftb)
			FreeDoubles(cliftb);
		if (clifta_p != clifta)
			FreeDoubles(clifta);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bliftt_p != bliftt)
			FreeDoubles(bliftt);
		if (bliftb_p != bliftb)
			FreeDoubles(bliftb);
		if (blifta_p != blifta)
			FreeDoubles(blifta);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (alift2_p != alift2)
			FreeDoubles(alift2);
		if (aliftt_p != aliftt)
			FreeDoubles(aliftt);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (bcdet_p != bcdet)
			FreeDoubles(bcdet);
		if (bcdetb_p != bcdetb)
			FreeDoubles(bcdetb);
		if (bcdeta_p != bcdeta)
			FreeDoubles(bcdeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (cdy_p != cdy)
			FreeDoubles(cdy);
		if (cdx_p != cdx)
			FreeDoubles(cdx);
		if (pdy3_p != pdy3)
			FreeDoubles(pdy3);
		if (pdx3_p != pdx3)
			FreeDoubles(pdx3);
		if (bdy_p != bdy)
			FreeDoubles(bdy);
		if (bdx_p != bdx)
			FreeDoubles(bdx);
		if (pdy2_p != pdy2)
			FreeDoubles(pdy2);
		if (pdx2_p != pdx2)
			FreeDoubles(pdx2);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (pdy1_p != pdy1)
			FreeDoubles(pdy1);
		if (pdx1_p != pdx1)
			FreeDoubles(pdx1);
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
		return inCirclexy_IIIE_exact<IT, ET>(p1, p2, p3, pdx, pdy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCirclexy_IIIE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &p3, double pdx, double pdy)
{
	Sign ret;
	ret = inCirclexy_IIIE_interval<IT, ET>(p1, p2, p3, pdx, pdy);
	if (is_sign_reliable(ret))
		return ret;
	return inCirclexy_IIIE_expansion<IT, ET>(p1, p2, p3, pdx, pdy);
}

template <typename IT, typename ET>
Sign inCirclexy_IIIE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &p3,
                     const GenericPoint3T<IT, ET> &pd)
{
	return inCirclexy_IIIE<IT, ET>(p1, p2, p3, pd.x(), pd.y());
}

template <typename IT, typename ET>
Sign inCirclexy_IIII_interval(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3,
                              const GenericPoint3T<IT, ET> &p4)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3) ||
	    !p4.getIntervalLambda(l4x, l4y, l4z, d4))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT l1xt   = l1x * d4;
	IT l1yt   = l1y * d4;
	IT l2xt   = l2x * d4;
	IT l2yt   = l2y * d4;
	IT l3xt   = l3x * d4;
	IT l3yt   = l3y * d4;
	IT l4x1   = l4x * d1;
	IT l4y1   = l4y * d1;
	IT adx    = l1xt - l4x1;
	IT ady    = l1yt - l4y1;
	IT l4x2   = l4x * d2;
	IT l4y2   = l4y * d2;
	IT bdx    = l2xt - l4x2;
	IT bdy    = l2yt - l4y2;
	IT l4x3   = l4x * d3;
	IT l4y3   = l4y * d3;
	IT cdx    = l3xt - l4x3;
	IT cdy    = l3yt - l4y3;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT aliftt = alifta + aliftb;
	IT alift2 = aliftt * d2;
	IT alift  = alift2 * d3;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT bliftt = blifta + bliftb;
	IT blift  = bliftt * d3;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT cliftt = clifta + cliftb;
	IT clift  = cliftt * d2;
	IT la     = alift * bcdet;
	IT lb     = blift * cadet;
	IT lc     = clift * abdet;
	IT lab2   = lc + lb;
	IT lab    = lab2 * d1;
	IT L      = lab + la;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IIII_exact(const GenericPoint3T<IT, ET> &p1,
                           const GenericPoint3T<IT, ET> &p2,
                           const GenericPoint3T<IT, ET> &p3,
                           const GenericPoint3T<IT, ET> &p4)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	p4.getExactLambda(l4x, l4y, l4z, d4);
	ET l1xt   = l1x * d4;
	ET l1yt   = l1y * d4;
	ET l2xt   = l2x * d4;
	ET l2yt   = l2y * d4;
	ET l3xt   = l3x * d4;
	ET l3yt   = l3y * d4;
	ET l4x1   = l4x * d1;
	ET l4y1   = l4y * d1;
	ET adx    = l1xt - l4x1;
	ET ady    = l1yt - l4y1;
	ET l4x2   = l4x * d2;
	ET l4y2   = l4y * d2;
	ET bdx    = l2xt - l4x2;
	ET bdy    = l2yt - l4y2;
	ET l4x3   = l4x * d3;
	ET l4y3   = l4y * d3;
	ET cdx    = l3xt - l4x3;
	ET cdy    = l3yt - l4y3;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET aliftt = alifta + aliftb;
	ET alift2 = aliftt * d2;
	ET alift  = alift2 * d3;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET bliftt = blifta + bliftb;
	ET blift  = bliftt * d3;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET cliftt = clifta + cliftb;
	ET clift  = cliftt * d2;
	ET la     = alift * bcdet;
	ET lb     = blift * cadet;
	ET lc     = clift * abdet;
	ET lab2   = lc + lb;
	ET lab    = lab2 * d1;
	ET L      = lab + la;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCirclexy_IIII_expansion(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3,
                               const GenericPoint3T<IT, ET> &p4)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[16],
	  *l1x = l1x_p, l1y_p[16], *l1y = l1y_p, l1z_p[16], *l1z = l1z_p, d1_p[16],
	  *d1 = d1_p, l2x_p[16], *l2x = l2x_p, l2y_p[16], *l2y = l2y_p, l2z_p[16],
	  *l2z = l2z_p, d2_p[16], *d2 = d2_p, l3x_p[16], *l3x = l3x_p, l3y_p[16],
	  *l3y = l3y_p, l3z_p[16], *l3z = l3z_p, d3_p[16], *d3 = d3_p, l4x_p[16],
	  *l4x = l4x_p, l4y_p[16], *l4y = l4y_p, l4z_p[16], *l4z = l4z_p, d4_p[16],
	  *d4       = d4_p;
	int l1x_len = 16, l1y_len = 16, l1z_len = 16, d1_len = 16, l2x_len = 16,
	    l2y_len = 16, l2z_len = 16, d2_len = 16, l3x_len = 16, l3y_len = 16,
	    l3z_len = 16, d3_len = 16, l4x_len = 16, l4y_len = 16, l4z_len = 16,
	    d4_len = 16;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	p4.getExpansionLambda(&l4x, l4x_len, &l4y, l4y_len, &l4z, l4z_len, &d4,
	                      d4_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) &&
	    (d4[d4_len - 1] != 0))
	{
		expansionObject o;
		double          l1xt_p[16], *l1xt = l1xt_p;
		int             l1xt_len =
		  o.Gen_Product_With_PreAlloc(l1x_len, l1x, d4_len, d4, &l1xt, 16);
		double l1yt_p[16], *l1yt = l1yt_p;
		int    l1yt_len =
		  o.Gen_Product_With_PreAlloc(l1y_len, l1y, d4_len, d4, &l1yt, 16);
		double l2xt_p[16], *l2xt = l2xt_p;
		int    l2xt_len =
		  o.Gen_Product_With_PreAlloc(l2x_len, l2x, d4_len, d4, &l2xt, 16);
		double l2yt_p[16], *l2yt = l2yt_p;
		int    l2yt_len =
		  o.Gen_Product_With_PreAlloc(l2y_len, l2y, d4_len, d4, &l2yt, 16);
		double l3xt_p[16], *l3xt = l3xt_p;
		int    l3xt_len =
		  o.Gen_Product_With_PreAlloc(l3x_len, l3x, d4_len, d4, &l3xt, 16);
		double l3yt_p[16], *l3yt = l3yt_p;
		int    l3yt_len =
		  o.Gen_Product_With_PreAlloc(l3y_len, l3y, d4_len, d4, &l3yt, 16);
		double l4x1_p[16], *l4x1 = l4x1_p;
		int    l4x1_len =
		  o.Gen_Product_With_PreAlloc(l4x_len, l4x, d1_len, d1, &l4x1, 16);
		double l4y1_p[16], *l4y1 = l4y1_p;
		int    l4y1_len =
		  o.Gen_Product_With_PreAlloc(l4y_len, l4y, d1_len, d1, &l4y1, 16);
		double adx_p[16], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1xt_len, l1xt, l4x1_len, l4x1, &adx, 16);
		double ady_p[16], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1yt_len, l1yt, l4y1_len, l4y1, &ady, 16);
		double l4x2_p[16], *l4x2 = l4x2_p;
		int    l4x2_len =
		  o.Gen_Product_With_PreAlloc(l4x_len, l4x, d2_len, d2, &l4x2, 16);
		double l4y2_p[16], *l4y2 = l4y2_p;
		int    l4y2_len =
		  o.Gen_Product_With_PreAlloc(l4y_len, l4y, d2_len, d2, &l4y2, 16);
		double bdx_p[16], *bdx = bdx_p;
		int    bdx_len =
		  o.Gen_Diff_With_PreAlloc(l2xt_len, l2xt, l4x2_len, l4x2, &bdx, 16);
		double bdy_p[16], *bdy = bdy_p;
		int    bdy_len =
		  o.Gen_Diff_With_PreAlloc(l2yt_len, l2yt, l4y2_len, l4y2, &bdy, 16);
		double l4x3_p[16], *l4x3 = l4x3_p;
		int    l4x3_len =
		  o.Gen_Product_With_PreAlloc(l4x_len, l4x, d3_len, d3, &l4x3, 16);
		double l4y3_p[16], *l4y3 = l4y3_p;
		int    l4y3_len =
		  o.Gen_Product_With_PreAlloc(l4y_len, l4y, d3_len, d3, &l4y3, 16);
		double cdx_p[16], *cdx = cdx_p;
		int    cdx_len =
		  o.Gen_Diff_With_PreAlloc(l3xt_len, l3xt, l4x3_len, l4x3, &cdx, 16);
		double cdy_p[16], *cdy = cdy_p;
		int    cdy_len =
		  o.Gen_Diff_With_PreAlloc(l3yt_len, l3yt, l4y3_len, l4y3, &cdy, 16);
		double abdeta_p[16], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 16);
		double abdetb_p[16], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 16);
		double abdet_p[16], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 16);
		double bcdeta_p[16], *bcdeta = bcdeta_p;
		int    bcdeta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 16);
		double bcdetb_p[16], *bcdetb = bcdetb_p;
		int    bcdetb_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 16);
		double bcdet_p[16], *bcdet = bcdet_p;
		int    bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len,
		                                            bcdetb, &bcdet, 16);
		double cadeta_p[16], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 16);
		double cadetb_p[16], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 16);
		double cadet_p[16], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 16);
		double alifta_p[16], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 16);
		double aliftb_p[16], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 16);
		double aliftt_p[16], *aliftt = aliftt_p;
		int    aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                            aliftb, &aliftt, 16);
		double alift2_p[16], *alift2 = alift2_p;
		int    alift2_len =
		  o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 16);
		double alift_p[16], *alift = alift_p;
		int    alift_len =
		  o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 16);
		double blifta_p[16], *blifta = blifta_p;
		int    blifta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 16);
		double bliftb_p[16], *bliftb = bliftb_p;
		int    bliftb_len =
		  o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 16);
		double bliftt_p[16], *bliftt = bliftt_p;
		int    bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len,
		                                            bliftb, &bliftt, 16);
		double blift_p[16], *blift = blift_p;
		int    blift_len =
		  o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 16);
		double clifta_p[16], *clifta = clifta_p;
		int    clifta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 16);
		double cliftb_p[16], *cliftb = cliftb_p;
		int    cliftb_len =
		  o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 16);
		double cliftt_p[16], *cliftt = cliftt_p;
		int    cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len,
		                                            cliftb, &cliftt, 16);
		double clift_p[16], *clift = clift_p;
		int    clift_len =
		  o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 16);
		double la_p[16], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 16);
		double lb_p[16], *lb = lb_p;
		int    lb_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 16);
		double lc_p[16], *lc = lc_p;
		int    lc_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 16);
		double lab2_p[16], *lab2 = lab2_p;
		int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 16);
		double lab_p[16], *lab = lab_p;
		int    lab_len =
		  o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 16);
		double L_p[16], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 16);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lab2_p != lab2)
			FreeDoubles(lab2);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (la_p != la)
			FreeDoubles(la);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (cliftt_p != cliftt)
			FreeDoubles(cliftt);
		if (cliftb_p != cliftb)
			FreeDoubles(cliftb);
		if (clifta_p != clifta)
			FreeDoubles(clifta);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bliftt_p != bliftt)
			FreeDoubles(bliftt);
		if (bliftb_p != bliftb)
			FreeDoubles(bliftb);
		if (blifta_p != blifta)
			FreeDoubles(blifta);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (alift2_p != alift2)
			FreeDoubles(alift2);
		if (aliftt_p != aliftt)
			FreeDoubles(aliftt);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (bcdet_p != bcdet)
			FreeDoubles(bcdet);
		if (bcdetb_p != bcdetb)
			FreeDoubles(bcdetb);
		if (bcdeta_p != bcdeta)
			FreeDoubles(bcdeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (cdy_p != cdy)
			FreeDoubles(cdy);
		if (cdx_p != cdx)
			FreeDoubles(cdx);
		if (l4y3_p != l4y3)
			FreeDoubles(l4y3);
		if (l4x3_p != l4x3)
			FreeDoubles(l4x3);
		if (bdy_p != bdy)
			FreeDoubles(bdy);
		if (bdx_p != bdx)
			FreeDoubles(bdx);
		if (l4y2_p != l4y2)
			FreeDoubles(l4y2);
		if (l4x2_p != l4x2)
			FreeDoubles(l4x2);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (l4y1_p != l4y1)
			FreeDoubles(l4y1);
		if (l4x1_p != l4x1)
			FreeDoubles(l4x1);
		if (l3yt_p != l3yt)
			FreeDoubles(l3yt);
		if (l3xt_p != l3xt)
			FreeDoubles(l3xt);
		if (l2yt_p != l2yt)
			FreeDoubles(l2yt);
		if (l2xt_p != l2xt)
			FreeDoubles(l2xt);
		if (l1yt_p != l1yt)
			FreeDoubles(l1yt);
		if (l1xt_p != l1xt)
			FreeDoubles(l1xt);
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

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l4x_p != l4x)
			FreeDoubles(l4x);
		if (l4y_p != l4y)
			FreeDoubles(l4y);
		if (l4z_p != l4z)
			FreeDoubles(l4z);
		if (d4_p != d4)
			FreeDoubles(d4);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return inCirclexy_IIII_exact<IT, ET>(p1, p2, p3, p4);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCirclexy_IIII(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &p3,
                     const GenericPoint3T<IT, ET> &p4)
{
	Sign ret;
	ret = inCirclexy_IIII_interval<IT, ET>(p1, p2, p3, p4);
	if (is_sign_reliable(ret))
		return ret;
	return inCirclexy_IIII_expansion<IT, ET>(p1, p2, p3, p4);
}

template <typename IT, typename ET>
Sign inCircle_IEEE_interval(const GenericPoint2T<IT, ET> &p1, IT pbx, IT pby,
                            IT pcx, IT pcy, IT pdx, IT pdy)
{
	IT l1x, l1y, d1;
	if (!p1.getIntervalLambda(l1x, l1y, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pdxt   = pdx * d1;
	IT pdyt   = pdy * d1;
	IT adx    = l1x - pdxt;
	IT ady    = l1y - pdyt;
	IT bdx    = pbx - pdx;
	IT bdy    = pby - pdy;
	IT cdx    = pcx - pdx;
	IT cdy    = pcy - pdy;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT alift  = alifta + aliftb;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT blift  = blifta + bliftb;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT clift  = clifta + cliftb;
	IT la     = alift * bcdet;
	IT lbt    = blift * cadet;
	IT lb     = lbt * d1;
	IT lct    = clift * abdet;
	IT lc     = lct * d1;
	IT lab    = la + lb;
	IT L      = lab + lc;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IEEE_exact(const GenericPoint2T<IT, ET> &p1, ET pbx, ET pby,
                         ET pcx, ET pcy, ET pdx, ET pdy)
{
	ET l1x, l1y, d1;
	p1.getExactLambda(l1x, l1y, d1);
	ET pdxt   = pdx * d1;
	ET pdyt   = pdy * d1;
	ET adx    = l1x - pdxt;
	ET ady    = l1y - pdyt;
	ET bdx    = pbx - pdx;
	ET bdy    = pby - pdy;
	ET cdx    = pcx - pdx;
	ET cdy    = pcy - pdy;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET alift  = alifta + aliftb;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET blift  = blifta + bliftb;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET clift  = clifta + cliftb;
	ET la     = alift * bcdet;
	ET lbt    = blift * cadet;
	ET lb     = lbt * d1;
	ET lct    = clift * abdet;
	ET lc     = lct * d1;
	ET lab    = la + lb;
	ET L      = lab + lc;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IEEE_expansion(const GenericPoint2T<IT, ET> &p1, double pbx,
                             double pby, double pcx, double pcy, double pdx,
                             double pdy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, d1_p[64], *d1 = d1_p;
	int    l1x_len = 64, l1y_len = 64, d1_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          pdxt_p[64], *pdxt = pdxt_p;
		int    pdxt_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdxt, 64);
		double pdyt_p[64], *pdyt = pdyt_p;
		int    pdyt_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdyt, 64);
		double adx_p[64], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdxt_len, pdxt, &adx, 64);
		double ady_p[64], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdyt_len, pdyt, &ady, 64);
		double bdx[2];
		o.two_Diff(pbx, pdx, bdx);
		double bdy[2];
		o.two_Diff(pby, pdy, bdy);
		double cdx[2];
		o.two_Diff(pcx, pdx, cdx);
		double cdy[2];
		o.two_Diff(pcy, pdy, cdy);
		double abdeta_p[64], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, 2, bdy, &abdeta, 64);
		double abdetb_p[64], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(2, bdx, ady_len, ady, &abdetb, 64);
		double abdet_p[64], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 64);
		double bcdeta[8];
		int    bcdeta_len = o.Gen_Product(2, bdx, 2, cdy, bcdeta);
		double bcdetb[8];
		int    bcdetb_len = o.Gen_Product(2, cdx, 2, bdy, bcdetb);
		double bcdet[16];
		int bcdet_len = o.Gen_Diff(bcdeta_len, bcdeta, bcdetb_len, bcdetb, bcdet);
		double cadeta_p[64], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 64);
		double cadetb_p[64], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 64);
		double cadet_p[64], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 64);
		double alifta_p[64], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 64);
		double aliftb_p[64], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 64);
		double alift_p[64], *alift = alift_p;
		int    alift_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                           aliftb, &alift, 64);
		double blifta[8];
		int    blifta_len = o.Gen_Product(2, bdx, 2, bdx, blifta);
		double bliftb[8];
		int    bliftb_len = o.Gen_Product(2, bdy, 2, bdy, bliftb);
		double blift[16];
		int    blift_len = o.Gen_Sum(blifta_len, blifta, bliftb_len, bliftb, blift);
		double clifta[8];
		int    clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
		double cliftb[8];
		int    cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
		double clift[16];
		int    clift_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, clift);
		double la_p[64], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 64);
		double lbt_p[64], *lbt = lbt_p;
		int    lbt_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lbt, 64);
		double lb_p[64], *lb = lb_p;
		int lb_len = o.Gen_Product_With_PreAlloc(lbt_len, lbt, d1_len, d1, &lb, 64);
		double lct_p[64], *lct = lct_p;
		int    lct_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lct, 64);
		double lc_p[64], *lc = lc_p;
		int lc_len = o.Gen_Product_With_PreAlloc(lct_len, lct, d1_len, d1, &lc, 64);
		double lab_p[64], *lab = lab_p;
		int    lab_len = o.Gen_Sum_With_PreAlloc(la_len, la, lb_len, lb, &lab, 64);
		double L_p[64], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, lc_len, lc, &L, 64);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lct_p != lct)
			FreeDoubles(lct);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (lbt_p != lbt)
			FreeDoubles(lbt);
		if (la_p != la)
			FreeDoubles(la);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (pdyt_p != pdyt)
			FreeDoubles(pdyt);
		if (pdxt_p != pdxt)
			FreeDoubles(pdxt);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return inCircle_IEEE_exact<IT, ET>(p1, pbx, pby, pcx, pcy, pdx, pdy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCircle_IEEE(const GenericPoint2T<IT, ET> &p1, double pbx, double pby,
                   double pcx, double pcy, double pdx, double pdy)
{
	Sign ret;
	ret = inCircle_IEEE_interval<IT, ET>(p1, pbx, pby, pcx, pcy, pdx, pdy);
	if (is_sign_reliable(ret))
		return ret;
	return inCircle_IEEE_expansion<IT, ET>(p1, pbx, pby, pcx, pcy, pdx, pdy);
}

template <typename IT, typename ET>
Sign inCircle_IEEE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &pb,
                   const GenericPoint2T<IT, ET> &pc,
                   const GenericPoint2T<IT, ET> &pd)
{
	return inCircle_IEEE<IT, ET>(p1, pb.x(), pb.y(), pc.x(), pc.y(), pd.x(),
	                             pd.y());
}

template <typename IT, typename ET>
Sign inCircle_IIEE_interval(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2, IT pcx, IT pcy,
                            IT pdx, IT pdy)
{
	IT l1x, l1y, d1, l2x, l2y, d2;
	if (!p1.getIntervalLambda(l1x, l1y, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pdx1   = pdx * d1;
	IT pdy1   = pdy * d1;
	IT adx    = l1x - pdx1;
	IT ady    = l1y - pdy1;
	IT pdx2   = pdx * d2;
	IT pdy2   = pdy * d2;
	IT bdx    = l2x - pdx2;
	IT bdy    = l2y - pdy2;
	IT cdx    = pcx - pdx;
	IT cdy    = pcy - pdy;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT aliftt = alifta + aliftb;
	IT alift  = aliftt * d2;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT blift  = blifta + bliftb;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT cliftt = clifta + cliftb;
	IT clift  = cliftt * d2;
	IT la     = alift * bcdet;
	IT lb     = blift * cadet;
	IT lc     = clift * abdet;
	IT lab    = lc + lb;
	IT lab2   = lab * d1;
	IT L      = lab2 + la;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IIEE_exact(const GenericPoint2T<IT, ET> &p1,
                         const GenericPoint2T<IT, ET> &p2, ET pcx, ET pcy,
                         ET pdx, ET pdy)
{
	ET l1x, l1y, d1, l2x, l2y, d2;
	p1.getExactLambda(l1x, l1y, d1);
	p2.getExactLambda(l2x, l2y, d2);
	ET pdx1   = pdx * d1;
	ET pdy1   = pdy * d1;
	ET adx    = l1x - pdx1;
	ET ady    = l1y - pdy1;
	ET pdx2   = pdx * d2;
	ET pdy2   = pdy * d2;
	ET bdx    = l2x - pdx2;
	ET bdy    = l2y - pdy2;
	ET cdx    = pcx - pdx;
	ET cdy    = pcy - pdy;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET aliftt = alifta + aliftb;
	ET alift  = aliftt * d2;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET blift  = blifta + bliftb;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET cliftt = clifta + cliftb;
	ET clift  = cliftt * d2;
	ET la     = alift * bcdet;
	ET lb     = blift * cadet;
	ET lc     = clift * abdet;
	ET lab    = lc + lb;
	ET lab2   = lab * d1;
	ET L      = lab2 + la;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IIEE_expansion(const GenericPoint2T<IT, ET> &p1,
                             const GenericPoint2T<IT, ET> &p2, double pcx,
                             double pcy, double pdx, double pdy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, d1_p[32], *d1 = d1_p,
	                  l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p, d2_p[32],
	                  *d2 = d2_p;
	int l1x_len = 32, l1y_len = 32, d1_len = 32, l2x_len = 32, l2y_len = 32,
	    d2_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &d2, d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          pdx1_p[32], *pdx1 = pdx1_p;
		int    pdx1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdx1, 32);
		double pdy1_p[32], *pdy1 = pdy1_p;
		int    pdy1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdy1, 32);
		double adx_p[32], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
		double ady_p[32], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
		double pdx2_p[32], *pdx2 = pdx2_p;
		int    pdx2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdx, &pdx2, 32);
		double pdy2_p[32], *pdy2 = pdy2_p;
		int    pdy2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdy, &pdy2, 32);
		double bdx_p[32], *bdx = bdx_p;
		int    bdx_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
		double bdy_p[32], *bdy = bdy_p;
		int    bdy_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
		double cdx[2];
		o.two_Diff(pcx, pdx, cdx);
		double cdy[2];
		o.two_Diff(pcy, pdy, cdy);
		double abdeta_p[32], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
		double abdetb_p[32], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
		double abdet_p[32], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 32);
		double bcdeta_p[32], *bcdeta = bcdeta_p;
		int    bcdeta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, 2, cdy, &bcdeta, 32);
		double bcdetb_p[32], *bcdetb = bcdetb_p;
		int    bcdetb_len =
		  o.Gen_Product_With_PreAlloc(2, cdx, bdy_len, bdy, &bcdetb, 32);
		double bcdet_p[32], *bcdet = bcdet_p;
		int    bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len,
		                                            bcdetb, &bcdet, 32);
		double cadeta_p[32], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(2, cdx, ady_len, ady, &cadeta, 32);
		double cadetb_p[32], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, 2, cdy, &cadetb, 32);
		double cadet_p[32], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 32);
		double alifta_p[32], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
		double aliftb_p[32], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
		double aliftt_p[32], *aliftt = aliftt_p;
		int    aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                            aliftb, &aliftt, 32);
		double alift_p[32], *alift = alift_p;
		int    alift_len =
		  o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift, 32);
		double blifta_p[32], *blifta = blifta_p;
		int    blifta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
		double bliftb_p[32], *bliftb = bliftb_p;
		int    bliftb_len =
		  o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
		double blift_p[32], *blift = blift_p;
		int    blift_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len,
		                                           bliftb, &blift, 32);
		double clifta[8];
		int    clifta_len = o.Gen_Product(2, cdx, 2, cdx, clifta);
		double cliftb[8];
		int    cliftb_len = o.Gen_Product(2, cdy, 2, cdy, cliftb);
		double cliftt[16];
		int cliftt_len = o.Gen_Sum(clifta_len, clifta, cliftb_len, cliftb, cliftt);
		double clift_p[32], *clift = clift_p;
		int    clift_len =
		  o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
		double la_p[32], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
		double lb_p[32], *lb = lb_p;
		int    lb_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
		double lc_p[32], *lc = lc_p;
		int    lc_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
		double lab_p[32], *lab = lab_p;
		int    lab_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab, 32);
		double lab2_p[32], *lab2 = lab2_p;
		int    lab2_len =
		  o.Gen_Product_With_PreAlloc(lab_len, lab, d1_len, d1, &lab2, 32);
		double L_p[32], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab2_len, lab2, la_len, la, &L, 32);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab2_p != lab2)
			FreeDoubles(lab2);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (la_p != la)
			FreeDoubles(la);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bliftb_p != bliftb)
			FreeDoubles(bliftb);
		if (blifta_p != blifta)
			FreeDoubles(blifta);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (aliftt_p != aliftt)
			FreeDoubles(aliftt);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (bcdet_p != bcdet)
			FreeDoubles(bcdet);
		if (bcdetb_p != bcdetb)
			FreeDoubles(bcdetb);
		if (bcdeta_p != bcdeta)
			FreeDoubles(bcdeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (bdy_p != bdy)
			FreeDoubles(bdy);
		if (bdx_p != bdx)
			FreeDoubles(bdx);
		if (pdy2_p != pdy2)
			FreeDoubles(pdy2);
		if (pdx2_p != pdx2)
			FreeDoubles(pdx2);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (pdy1_p != pdy1)
			FreeDoubles(pdy1);
		if (pdx1_p != pdx1)
			FreeDoubles(pdx1);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return inCircle_IIEE_exact<IT, ET>(p1, p2, pcx, pcy, pdx, pdy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCircle_IIEE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2, double pcx, double pcy,
                   double pdx, double pdy)
{
	Sign ret;
	ret = inCircle_IIEE_interval<IT, ET>(p1, p2, pcx, pcy, pdx, pdy);
	if (is_sign_reliable(ret))
		return ret;
	return inCircle_IIEE_expansion<IT, ET>(p1, p2, pcx, pcy, pdx, pdy);
}

template <typename IT, typename ET>
Sign inCircle_IIEE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &pc,
                   const GenericPoint2T<IT, ET> &pd)
{
	return inCircle_IIEE<IT, ET>(p1, p2, pc.x(), pc.y(), pd.x(), pd.y());
}

template <typename IT, typename ET>
Sign inCircle_IIIE_interval(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2,
                            const GenericPoint2T<IT, ET> &p3, IT pdx, IT pdy)
{
	IT l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
	if (!p1.getIntervalLambda(l1x, l1y, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pdx1   = pdx * d1;
	IT pdy1   = pdy * d1;
	IT adx    = l1x - pdx1;
	IT ady    = l1y - pdy1;
	IT pdx2   = pdx * d2;
	IT pdy2   = pdy * d2;
	IT bdx    = l2x - pdx2;
	IT bdy    = l2y - pdy2;
	IT pdx3   = pdx * d3;
	IT pdy3   = pdy * d3;
	IT cdx    = l3x - pdx3;
	IT cdy    = l3y - pdy3;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT aliftt = alifta + aliftb;
	IT alift2 = aliftt * d2;
	IT alift  = alift2 * d3;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT bliftt = blifta + bliftb;
	IT blift  = bliftt * d3;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT cliftt = clifta + cliftb;
	IT clift  = cliftt * d2;
	IT la     = alift * bcdet;
	IT lb     = blift * cadet;
	IT lc     = clift * abdet;
	IT lab2   = lc + lb;
	IT lab    = lab2 * d1;
	IT L      = lab + la;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IIIE_exact(const GenericPoint2T<IT, ET> &p1,
                         const GenericPoint2T<IT, ET> &p2,
                         const GenericPoint2T<IT, ET> &p3, ET pdx, ET pdy)
{
	ET l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
	p1.getExactLambda(l1x, l1y, d1);
	p2.getExactLambda(l2x, l2y, d2);
	p3.getExactLambda(l3x, l3y, d3);
	ET pdx1   = pdx * d1;
	ET pdy1   = pdy * d1;
	ET adx    = l1x - pdx1;
	ET ady    = l1y - pdy1;
	ET pdx2   = pdx * d2;
	ET pdy2   = pdy * d2;
	ET bdx    = l2x - pdx2;
	ET bdy    = l2y - pdy2;
	ET pdx3   = pdx * d3;
	ET pdy3   = pdy * d3;
	ET cdx    = l3x - pdx3;
	ET cdy    = l3y - pdy3;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET aliftt = alifta + aliftb;
	ET alift2 = aliftt * d2;
	ET alift  = alift2 * d3;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET bliftt = blifta + bliftb;
	ET blift  = bliftt * d3;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET cliftt = clifta + cliftb;
	ET clift  = cliftt * d2;
	ET la     = alift * bcdet;
	ET lb     = blift * cadet;
	ET lc     = clift * abdet;
	ET lab2   = lc + lb;
	ET lab    = lab2 * d1;
	ET L      = lab + la;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IIIE_expansion(const GenericPoint2T<IT, ET> &p1,
                             const GenericPoint2T<IT, ET> &p2,
                             const GenericPoint2T<IT, ET> &p3, double pdx,
                             double pdy)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, d1_p[32], *d1 = d1_p,
	                  l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p, d2_p[32],
	                  *d2 = d2_p, l3x_p[32], *l3x = l3x_p, l3y_p[32],
	                  *l3y = l3y_p, d3_p[32], *d3 = d3_p;
	int l1x_len = 32, l1y_len = 32, d1_len = 32, l2x_len = 32, l2y_len = 32,
	    d2_len = 32, l3x_len = 32, l3y_len = 32, d3_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &d2, d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &d3, d3_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          pdx1_p[32], *pdx1 = pdx1_p;
		int    pdx1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdx, &pdx1, 32);
		double pdy1_p[32], *pdy1 = pdy1_p;
		int    pdy1_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pdy, &pdy1, 32);
		double adx_p[32], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pdx1_len, pdx1, &adx, 32);
		double ady_p[32], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, pdy1_len, pdy1, &ady, 32);
		double pdx2_p[32], *pdx2 = pdx2_p;
		int    pdx2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdx, &pdx2, 32);
		double pdy2_p[32], *pdy2 = pdy2_p;
		int    pdy2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pdy, &pdy2, 32);
		double bdx_p[32], *bdx = bdx_p;
		int    bdx_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pdx2_len, pdx2, &bdx, 32);
		double bdy_p[32], *bdy = bdy_p;
		int    bdy_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, pdy2_len, pdy2, &bdy, 32);
		double pdx3_p[32], *pdx3 = pdx3_p;
		int    pdx3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pdx, &pdx3, 32);
		double pdy3_p[32], *pdy3 = pdy3_p;
		int    pdy3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pdy, &pdy3, 32);
		double cdx_p[32], *cdx = cdx_p;
		int    cdx_len =
		  o.Gen_Diff_With_PreAlloc(l3x_len, l3x, pdx3_len, pdx3, &cdx, 32);
		double cdy_p[32], *cdy = cdy_p;
		int    cdy_len =
		  o.Gen_Diff_With_PreAlloc(l3y_len, l3y, pdy3_len, pdy3, &cdy, 32);
		double abdeta_p[32], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
		double abdetb_p[32], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
		double abdet_p[32], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 32);
		double bcdeta_p[32], *bcdeta = bcdeta_p;
		int    bcdeta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 32);
		double bcdetb_p[32], *bcdetb = bcdetb_p;
		int    bcdetb_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 32);
		double bcdet_p[32], *bcdet = bcdet_p;
		int    bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len,
		                                            bcdetb, &bcdet, 32);
		double cadeta_p[32], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 32);
		double cadetb_p[32], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 32);
		double cadet_p[32], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 32);
		double alifta_p[32], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
		double aliftb_p[32], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
		double aliftt_p[32], *aliftt = aliftt_p;
		int    aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                            aliftb, &aliftt, 32);
		double alift2_p[32], *alift2 = alift2_p;
		int    alift2_len =
		  o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 32);
		double alift_p[32], *alift = alift_p;
		int    alift_len =
		  o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 32);
		double blifta_p[32], *blifta = blifta_p;
		int    blifta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
		double bliftb_p[32], *bliftb = bliftb_p;
		int    bliftb_len =
		  o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
		double bliftt_p[32], *bliftt = bliftt_p;
		int    bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len,
		                                            bliftb, &bliftt, 32);
		double blift_p[32], *blift = blift_p;
		int    blift_len =
		  o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 32);
		double clifta_p[32], *clifta = clifta_p;
		int    clifta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 32);
		double cliftb_p[32], *cliftb = cliftb_p;
		int    cliftb_len =
		  o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 32);
		double cliftt_p[32], *cliftt = cliftt_p;
		int    cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len,
		                                            cliftb, &cliftt, 32);
		double clift_p[32], *clift = clift_p;
		int    clift_len =
		  o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
		double la_p[32], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
		double lb_p[32], *lb = lb_p;
		int    lb_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
		double lc_p[32], *lc = lc_p;
		int    lc_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
		double lab2_p[32], *lab2 = lab2_p;
		int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 32);
		double lab_p[32], *lab = lab_p;
		int    lab_len =
		  o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 32);
		double L_p[32], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 32);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lab2_p != lab2)
			FreeDoubles(lab2);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (la_p != la)
			FreeDoubles(la);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (cliftt_p != cliftt)
			FreeDoubles(cliftt);
		if (cliftb_p != cliftb)
			FreeDoubles(cliftb);
		if (clifta_p != clifta)
			FreeDoubles(clifta);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bliftt_p != bliftt)
			FreeDoubles(bliftt);
		if (bliftb_p != bliftb)
			FreeDoubles(bliftb);
		if (blifta_p != blifta)
			FreeDoubles(blifta);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (alift2_p != alift2)
			FreeDoubles(alift2);
		if (aliftt_p != aliftt)
			FreeDoubles(aliftt);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (bcdet_p != bcdet)
			FreeDoubles(bcdet);
		if (bcdetb_p != bcdetb)
			FreeDoubles(bcdetb);
		if (bcdeta_p != bcdeta)
			FreeDoubles(bcdeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (cdy_p != cdy)
			FreeDoubles(cdy);
		if (cdx_p != cdx)
			FreeDoubles(cdx);
		if (pdy3_p != pdy3)
			FreeDoubles(pdy3);
		if (pdx3_p != pdx3)
			FreeDoubles(pdx3);
		if (bdy_p != bdy)
			FreeDoubles(bdy);
		if (bdx_p != bdx)
			FreeDoubles(bdx);
		if (pdy2_p != pdy2)
			FreeDoubles(pdy2);
		if (pdx2_p != pdx2)
			FreeDoubles(pdx2);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (pdy1_p != pdy1)
			FreeDoubles(pdy1);
		if (pdx1_p != pdx1)
			FreeDoubles(pdx1);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return inCircle_IIIE_exact<IT, ET>(p1, p2, p3, pdx, pdy);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCircle_IIIE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &p3, double pdx, double pdy)
{
	Sign ret;
	ret = inCircle_IIIE_interval<IT, ET>(p1, p2, p3, pdx, pdy);
	if (is_sign_reliable(ret))
		return ret;
	return inCircle_IIIE_expansion<IT, ET>(p1, p2, p3, pdx, pdy);
}

template <typename IT, typename ET>
Sign inCircle_IIIE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &p3,
                   const GenericPoint2T<IT, ET> &pd)
{
	return inCircle_IIIE<IT, ET>(p1, p2, p3, pd.x(), pd.y());
}

template <typename IT, typename ET>
Sign inCircle_IIII_interval(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2,
                            const GenericPoint2T<IT, ET> &p3,
                            const GenericPoint2T<IT, ET> &p4)
{
	IT l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, l4x, l4y, d4;
	if (!p1.getIntervalLambda(l1x, l1y, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, d3) ||
	    !p4.getIntervalLambda(l4x, l4y, d4))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT l1xt   = l1x * d4;
	IT l1yt   = l1y * d4;
	IT l2xt   = l2x * d4;
	IT l2yt   = l2y * d4;
	IT l3xt   = l3x * d4;
	IT l3yt   = l3y * d4;
	IT l4x1   = l4x * d1;
	IT l4y1   = l4y * d1;
	IT adx    = l1xt - l4x1;
	IT ady    = l1yt - l4y1;
	IT l4x2   = l4x * d2;
	IT l4y2   = l4y * d2;
	IT bdx    = l2xt - l4x2;
	IT bdy    = l2yt - l4y2;
	IT l4x3   = l4x * d3;
	IT l4y3   = l4y * d3;
	IT cdx    = l3xt - l4x3;
	IT cdy    = l3yt - l4y3;
	IT abdeta = adx * bdy;
	IT abdetb = bdx * ady;
	IT abdet  = abdeta - abdetb;
	IT bcdeta = bdx * cdy;
	IT bcdetb = cdx * bdy;
	IT bcdet  = bcdeta - bcdetb;
	IT cadeta = cdx * ady;
	IT cadetb = adx * cdy;
	IT cadet  = cadeta - cadetb;
	IT alifta = adx * adx;
	IT aliftb = ady * ady;
	IT aliftt = alifta + aliftb;
	IT alift2 = aliftt * d2;
	IT alift  = alift2 * d3;
	IT blifta = bdx * bdx;
	IT bliftb = bdy * bdy;
	IT bliftt = blifta + bliftb;
	IT blift  = bliftt * d3;
	IT clifta = cdx * cdx;
	IT cliftb = cdy * cdy;
	IT cliftt = clifta + cliftb;
	IT clift  = cliftt * d2;
	IT la     = alift * bcdet;
	IT lb     = blift * cadet;
	IT lc     = clift * abdet;
	IT lab2   = lc + lb;
	IT lab    = lab2 * d1;
	IT L      = lab + la;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IIII_exact(const GenericPoint2T<IT, ET> &p1,
                         const GenericPoint2T<IT, ET> &p2,
                         const GenericPoint2T<IT, ET> &p3,
                         const GenericPoint2T<IT, ET> &p4)
{
	ET l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, l4x, l4y, d4;
	p1.getExactLambda(l1x, l1y, d1);
	p2.getExactLambda(l2x, l2y, d2);
	p3.getExactLambda(l3x, l3y, d3);
	p4.getExactLambda(l4x, l4y, d4);
	ET l1xt   = l1x * d4;
	ET l1yt   = l1y * d4;
	ET l2xt   = l2x * d4;
	ET l2yt   = l2y * d4;
	ET l3xt   = l3x * d4;
	ET l3yt   = l3y * d4;
	ET l4x1   = l4x * d1;
	ET l4y1   = l4y * d1;
	ET adx    = l1xt - l4x1;
	ET ady    = l1yt - l4y1;
	ET l4x2   = l4x * d2;
	ET l4y2   = l4y * d2;
	ET bdx    = l2xt - l4x2;
	ET bdy    = l2yt - l4y2;
	ET l4x3   = l4x * d3;
	ET l4y3   = l4y * d3;
	ET cdx    = l3xt - l4x3;
	ET cdy    = l3yt - l4y3;
	ET abdeta = adx * bdy;
	ET abdetb = bdx * ady;
	ET abdet  = abdeta - abdetb;
	ET bcdeta = bdx * cdy;
	ET bcdetb = cdx * bdy;
	ET bcdet  = bcdeta - bcdetb;
	ET cadeta = cdx * ady;
	ET cadetb = adx * cdy;
	ET cadet  = cadeta - cadetb;
	ET alifta = adx * adx;
	ET aliftb = ady * ady;
	ET aliftt = alifta + aliftb;
	ET alift2 = aliftt * d2;
	ET alift  = alift2 * d3;
	ET blifta = bdx * bdx;
	ET bliftb = bdy * bdy;
	ET bliftt = blifta + bliftb;
	ET blift  = bliftt * d3;
	ET clifta = cdx * cdx;
	ET cliftb = cdy * cdy;
	ET cliftt = clifta + cliftb;
	ET clift  = cliftt * d2;
	ET la     = alift * bcdet;
	ET lb     = blift * cadet;
	ET lc     = clift * abdet;
	ET lab2   = lc + lb;
	ET lab    = lab2 * d1;
	ET L      = lab + la;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign inCircle_IIII_expansion(const GenericPoint2T<IT, ET> &p1,
                             const GenericPoint2T<IT, ET> &p2,
                             const GenericPoint2T<IT, ET> &p3,
                             const GenericPoint2T<IT, ET> &p4)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, d1_p[32], *d1 = d1_p,
	                  l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p, d2_p[32],
	                  *d2 = d2_p, l3x_p[32], *l3x = l3x_p, l3y_p[32],
	                  *l3y = l3y_p, d3_p[32], *d3 = d3_p, l4x_p[32], *l4x = l4x_p,
	                  l4y_p[32], *l4y = l4y_p, d4_p[32], *d4 = d4_p;
	int l1x_len = 32, l1y_len = 32, d1_len = 32, l2x_len = 32, l2y_len = 32,
	    d2_len = 32, l3x_len = 32, l3y_len = 32, d3_len = 32, l4x_len = 32,
	    l4y_len = 32, d4_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &d2, d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &d3, d3_len);
	p4.getExpansionLambda(&l4x, l4x_len, &l4y, l4y_len, &d4, d4_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) &&
	    (d4[d4_len - 1] != 0))
	{
		expansionObject o;
		double          l1xt_p[32], *l1xt = l1xt_p;
		int             l1xt_len =
		  o.Gen_Product_With_PreAlloc(l1x_len, l1x, d4_len, d4, &l1xt, 32);
		double l1yt_p[32], *l1yt = l1yt_p;
		int    l1yt_len =
		  o.Gen_Product_With_PreAlloc(l1y_len, l1y, d4_len, d4, &l1yt, 32);
		double l2xt_p[32], *l2xt = l2xt_p;
		int    l2xt_len =
		  o.Gen_Product_With_PreAlloc(l2x_len, l2x, d4_len, d4, &l2xt, 32);
		double l2yt_p[32], *l2yt = l2yt_p;
		int    l2yt_len =
		  o.Gen_Product_With_PreAlloc(l2y_len, l2y, d4_len, d4, &l2yt, 32);
		double l3xt_p[32], *l3xt = l3xt_p;
		int    l3xt_len =
		  o.Gen_Product_With_PreAlloc(l3x_len, l3x, d4_len, d4, &l3xt, 32);
		double l3yt_p[32], *l3yt = l3yt_p;
		int    l3yt_len =
		  o.Gen_Product_With_PreAlloc(l3y_len, l3y, d4_len, d4, &l3yt, 32);
		double l4x1_p[32], *l4x1 = l4x1_p;
		int    l4x1_len =
		  o.Gen_Product_With_PreAlloc(l4x_len, l4x, d1_len, d1, &l4x1, 32);
		double l4y1_p[32], *l4y1 = l4y1_p;
		int    l4y1_len =
		  o.Gen_Product_With_PreAlloc(l4y_len, l4y, d1_len, d1, &l4y1, 32);
		double adx_p[32], *adx = adx_p;
		int    adx_len =
		  o.Gen_Diff_With_PreAlloc(l1xt_len, l1xt, l4x1_len, l4x1, &adx, 32);
		double ady_p[32], *ady = ady_p;
		int    ady_len =
		  o.Gen_Diff_With_PreAlloc(l1yt_len, l1yt, l4y1_len, l4y1, &ady, 32);
		double l4x2_p[32], *l4x2 = l4x2_p;
		int    l4x2_len =
		  o.Gen_Product_With_PreAlloc(l4x_len, l4x, d2_len, d2, &l4x2, 32);
		double l4y2_p[32], *l4y2 = l4y2_p;
		int    l4y2_len =
		  o.Gen_Product_With_PreAlloc(l4y_len, l4y, d2_len, d2, &l4y2, 32);
		double bdx_p[32], *bdx = bdx_p;
		int    bdx_len =
		  o.Gen_Diff_With_PreAlloc(l2xt_len, l2xt, l4x2_len, l4x2, &bdx, 32);
		double bdy_p[32], *bdy = bdy_p;
		int    bdy_len =
		  o.Gen_Diff_With_PreAlloc(l2yt_len, l2yt, l4y2_len, l4y2, &bdy, 32);
		double l4x3_p[32], *l4x3 = l4x3_p;
		int    l4x3_len =
		  o.Gen_Product_With_PreAlloc(l4x_len, l4x, d3_len, d3, &l4x3, 32);
		double l4y3_p[32], *l4y3 = l4y3_p;
		int    l4y3_len =
		  o.Gen_Product_With_PreAlloc(l4y_len, l4y, d3_len, d3, &l4y3, 32);
		double cdx_p[32], *cdx = cdx_p;
		int    cdx_len =
		  o.Gen_Diff_With_PreAlloc(l3xt_len, l3xt, l4x3_len, l4x3, &cdx, 32);
		double cdy_p[32], *cdy = cdy_p;
		int    cdy_len =
		  o.Gen_Diff_With_PreAlloc(l3yt_len, l3yt, l4y3_len, l4y3, &cdy, 32);
		double abdeta_p[32], *abdeta = abdeta_p;
		int    abdeta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, bdy_len, bdy, &abdeta, 32);
		double abdetb_p[32], *abdetb = abdetb_p;
		int    abdetb_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, ady_len, ady, &abdetb, 32);
		double abdet_p[32], *abdet = abdet_p;
		int    abdet_len = o.Gen_Diff_With_PreAlloc(abdeta_len, abdeta, abdetb_len,
		                                            abdetb, &abdet, 32);
		double bcdeta_p[32], *bcdeta = bcdeta_p;
		int    bcdeta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, cdy_len, cdy, &bcdeta, 32);
		double bcdetb_p[32], *bcdetb = bcdetb_p;
		int    bcdetb_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, bdy_len, bdy, &bcdetb, 32);
		double bcdet_p[32], *bcdet = bcdet_p;
		int    bcdet_len = o.Gen_Diff_With_PreAlloc(bcdeta_len, bcdeta, bcdetb_len,
		                                            bcdetb, &bcdet, 32);
		double cadeta_p[32], *cadeta = cadeta_p;
		int    cadeta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, ady_len, ady, &cadeta, 32);
		double cadetb_p[32], *cadetb = cadetb_p;
		int    cadetb_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, cdy_len, cdy, &cadetb, 32);
		double cadet_p[32], *cadet = cadet_p;
		int    cadet_len = o.Gen_Diff_With_PreAlloc(cadeta_len, cadeta, cadetb_len,
		                                            cadetb, &cadet, 32);
		double alifta_p[32], *alifta = alifta_p;
		int    alifta_len =
		  o.Gen_Product_With_PreAlloc(adx_len, adx, adx_len, adx, &alifta, 32);
		double aliftb_p[32], *aliftb = aliftb_p;
		int    aliftb_len =
		  o.Gen_Product_With_PreAlloc(ady_len, ady, ady_len, ady, &aliftb, 32);
		double aliftt_p[32], *aliftt = aliftt_p;
		int    aliftt_len = o.Gen_Sum_With_PreAlloc(alifta_len, alifta, aliftb_len,
		                                            aliftb, &aliftt, 32);
		double alift2_p[32], *alift2 = alift2_p;
		int    alift2_len =
		  o.Gen_Product_With_PreAlloc(aliftt_len, aliftt, d2_len, d2, &alift2, 32);
		double alift_p[32], *alift = alift_p;
		int    alift_len =
		  o.Gen_Product_With_PreAlloc(alift2_len, alift2, d3_len, d3, &alift, 32);
		double blifta_p[32], *blifta = blifta_p;
		int    blifta_len =
		  o.Gen_Product_With_PreAlloc(bdx_len, bdx, bdx_len, bdx, &blifta, 32);
		double bliftb_p[32], *bliftb = bliftb_p;
		int    bliftb_len =
		  o.Gen_Product_With_PreAlloc(bdy_len, bdy, bdy_len, bdy, &bliftb, 32);
		double bliftt_p[32], *bliftt = bliftt_p;
		int    bliftt_len = o.Gen_Sum_With_PreAlloc(blifta_len, blifta, bliftb_len,
		                                            bliftb, &bliftt, 32);
		double blift_p[32], *blift = blift_p;
		int    blift_len =
		  o.Gen_Product_With_PreAlloc(bliftt_len, bliftt, d3_len, d3, &blift, 32);
		double clifta_p[32], *clifta = clifta_p;
		int    clifta_len =
		  o.Gen_Product_With_PreAlloc(cdx_len, cdx, cdx_len, cdx, &clifta, 32);
		double cliftb_p[32], *cliftb = cliftb_p;
		int    cliftb_len =
		  o.Gen_Product_With_PreAlloc(cdy_len, cdy, cdy_len, cdy, &cliftb, 32);
		double cliftt_p[32], *cliftt = cliftt_p;
		int    cliftt_len = o.Gen_Sum_With_PreAlloc(clifta_len, clifta, cliftb_len,
		                                            cliftb, &cliftt, 32);
		double clift_p[32], *clift = clift_p;
		int    clift_len =
		  o.Gen_Product_With_PreAlloc(cliftt_len, cliftt, d2_len, d2, &clift, 32);
		double la_p[32], *la = la_p;
		int    la_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcdet_len, bcdet, &la, 32);
		double lb_p[32], *lb = lb_p;
		int    lb_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cadet_len, cadet, &lb, 32);
		double lc_p[32], *lc = lc_p;
		int    lc_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, abdet_len, abdet, &lc, 32);
		double lab2_p[32], *lab2 = lab2_p;
		int lab2_len = o.Gen_Sum_With_PreAlloc(lc_len, lc, lb_len, lb, &lab2, 32);
		double lab_p[32], *lab = lab_p;
		int    lab_len =
		  o.Gen_Product_With_PreAlloc(lab2_len, lab2, d1_len, d1, &lab, 32);
		double L_p[32], *L = L_p;
		int    L_len = o.Gen_Sum_With_PreAlloc(lab_len, lab, la_len, la, &L, 32);

		return_value = L[L_len - 1];
		if (L_p != L)
			FreeDoubles(L);
		if (lab_p != lab)
			FreeDoubles(lab);
		if (lab2_p != lab2)
			FreeDoubles(lab2);
		if (lc_p != lc)
			FreeDoubles(lc);
		if (lb_p != lb)
			FreeDoubles(lb);
		if (la_p != la)
			FreeDoubles(la);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (cliftt_p != cliftt)
			FreeDoubles(cliftt);
		if (cliftb_p != cliftb)
			FreeDoubles(cliftb);
		if (clifta_p != clifta)
			FreeDoubles(clifta);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bliftt_p != bliftt)
			FreeDoubles(bliftt);
		if (bliftb_p != bliftb)
			FreeDoubles(bliftb);
		if (blifta_p != blifta)
			FreeDoubles(blifta);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (alift2_p != alift2)
			FreeDoubles(alift2);
		if (aliftt_p != aliftt)
			FreeDoubles(aliftt);
		if (aliftb_p != aliftb)
			FreeDoubles(aliftb);
		if (alifta_p != alifta)
			FreeDoubles(alifta);
		if (cadet_p != cadet)
			FreeDoubles(cadet);
		if (cadetb_p != cadetb)
			FreeDoubles(cadetb);
		if (cadeta_p != cadeta)
			FreeDoubles(cadeta);
		if (bcdet_p != bcdet)
			FreeDoubles(bcdet);
		if (bcdetb_p != bcdetb)
			FreeDoubles(bcdetb);
		if (bcdeta_p != bcdeta)
			FreeDoubles(bcdeta);
		if (abdet_p != abdet)
			FreeDoubles(abdet);
		if (abdetb_p != abdetb)
			FreeDoubles(abdetb);
		if (abdeta_p != abdeta)
			FreeDoubles(abdeta);
		if (cdy_p != cdy)
			FreeDoubles(cdy);
		if (cdx_p != cdx)
			FreeDoubles(cdx);
		if (l4y3_p != l4y3)
			FreeDoubles(l4y3);
		if (l4x3_p != l4x3)
			FreeDoubles(l4x3);
		if (bdy_p != bdy)
			FreeDoubles(bdy);
		if (bdx_p != bdx)
			FreeDoubles(bdx);
		if (l4y2_p != l4y2)
			FreeDoubles(l4y2);
		if (l4x2_p != l4x2)
			FreeDoubles(l4x2);
		if (ady_p != ady)
			FreeDoubles(ady);
		if (adx_p != adx)
			FreeDoubles(adx);
		if (l4y1_p != l4y1)
			FreeDoubles(l4y1);
		if (l4x1_p != l4x1)
			FreeDoubles(l4x1);
		if (l3yt_p != l3yt)
			FreeDoubles(l3yt);
		if (l3xt_p != l3xt)
			FreeDoubles(l3xt);
		if (l2yt_p != l2yt)
			FreeDoubles(l2yt);
		if (l2xt_p != l2xt)
			FreeDoubles(l2xt);
		if (l1yt_p != l1yt)
			FreeDoubles(l1yt);
		if (l1xt_p != l1xt)
			FreeDoubles(l1xt);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l4x_p != l4x)
			FreeDoubles(l4x);
		if (l4y_p != l4y)
			FreeDoubles(l4y);
		if (d4_p != d4)
			FreeDoubles(d4);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return inCircle_IIII_exact<IT, ET>(p1, p2, p3, p4);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inCircle_IIII(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &p3,
                   const GenericPoint2T<IT, ET> &p4)
{
	Sign ret;
	ret = inCircle_IIII_interval<IT, ET>(p1, p2, p3, p4);
	if (is_sign_reliable(ret))
		return ret;
	return inCircle_IIII_expansion<IT, ET>(p1, p2, p3, p4);
}

template <typename IT, typename ET>
Sign inSphere_IEEEE_interval(const GenericPoint3T<IT, ET> &p1, IT pbx, IT pby,
                             IT pbz, IT pcx, IT pcy, IT pcz, IT pdx, IT pdy,
                             IT pdz, IT pex, IT pey, IT pez)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pexd   = pex * d1;
	IT peyd   = pey * d1;
	IT pezd   = pez * d1;
	IT aex    = l1x - pexd;
	IT aey    = l1y - peyd;
	IT aez    = l1z - pezd;
	IT bex    = pbx - pex;
	IT bey    = pby - pey;
	IT bez    = pbz - pez;
	IT cex    = pcx - pex;
	IT cey    = pcy - pey;
	IT cez    = pcz - pez;
	IT dex    = pdx - pex;
	IT dey    = pdy - pey;
	IT dez    = pdz - pez;
	IT aexbey = aex * bey;
	IT bexaey = bex * aey;
	IT ab     = aexbey - bexaey;
	IT bexcey = bex * cey;
	IT cexbey = cex * bey;
	IT bc     = bexcey - cexbey;
	IT cexdey = cex * dey;
	IT dexcey = dex * cey;
	IT cd     = cexdey - dexcey;
	IT dexaey = dex * aey;
	IT aexdey = aex * dey;
	IT da     = dexaey - aexdey;
	IT aexcey = aex * cey;
	IT cexaey = cex * aey;
	IT ac     = aexcey - cexaey;
	IT bexdey = bex * dey;
	IT dexbey = dex * bey;
	IT bd     = bexdey - dexbey;
	IT abc1   = aez * bc;
	IT abc2   = bez * ac;
	IT abc3   = cez * ab;
	IT abc4   = abc1 + abc3;
	IT abc    = abc4 - abc2;
	IT bcd1   = bez * cd;
	IT bcd2   = cez * bd;
	IT bcd3   = dez * bc;
	IT bcd4   = bcd1 + bcd3;
	IT bcd    = bcd4 - bcd2;
	IT cda1   = cez * da;
	IT cda2   = dez * ac;
	IT cda3   = aez * cd;
	IT cda4   = cda1 + cda3;
	IT cda    = cda4 + cda2;
	IT dab1   = dez * ab;
	IT dab2   = aez * bd;
	IT dab3   = bez * da;
	IT dab4   = dab1 + dab3;
	IT dab    = dab4 + dab2;
	IT al1    = aex * aex;
	IT al2    = aey * aey;
	IT al3    = aez * aez;
	IT al4    = al1 + al2;
	IT alift  = al4 + al3;
	IT bl1    = bex * bex;
	IT bl2    = bey * bey;
	IT bl3    = bez * bez;
	IT bl4    = bl1 + bl2;
	IT blift  = bl4 + bl3;
	IT cl1    = cex * cex;
	IT cl2    = cey * cey;
	IT cl3    = cez * cez;
	IT cl4    = cl1 + cl2;
	IT clift  = cl4 + cl3;
	IT dl1    = dex * dex;
	IT dl2    = dey * dey;
	IT dl3    = dez * dez;
	IT dl4    = dl1 + dl2;
	IT dlift  = dl4 + dl3;
	IT ds1    = dlift * abc;
	IT ds2    = clift * dab;
	IT dlp    = ds2 - ds1;
	IT dl     = dlp * d1;
	IT dr1p   = blift * cda;
	IT dr1    = dr1p * d1;
	IT dr2    = alift * bcd;
	IT dr     = dr2 - dr1;
	IT det    = dl + dr;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IEEEE_exact(const GenericPoint3T<IT, ET> &p1, ET pbx, ET pby,
                          ET pbz, ET pcx, ET pcy, ET pcz, ET pdx, ET pdy,
                          ET pdz, ET pex, ET pey, ET pez)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET pexd   = pex * d1;
	ET peyd   = pey * d1;
	ET pezd   = pez * d1;
	ET aex    = l1x - pexd;
	ET aey    = l1y - peyd;
	ET aez    = l1z - pezd;
	ET bex    = pbx - pex;
	ET bey    = pby - pey;
	ET bez    = pbz - pez;
	ET cex    = pcx - pex;
	ET cey    = pcy - pey;
	ET cez    = pcz - pez;
	ET dex    = pdx - pex;
	ET dey    = pdy - pey;
	ET dez    = pdz - pez;
	ET aexbey = aex * bey;
	ET bexaey = bex * aey;
	ET ab     = aexbey - bexaey;
	ET bexcey = bex * cey;
	ET cexbey = cex * bey;
	ET bc     = bexcey - cexbey;
	ET cexdey = cex * dey;
	ET dexcey = dex * cey;
	ET cd     = cexdey - dexcey;
	ET dexaey = dex * aey;
	ET aexdey = aex * dey;
	ET da     = dexaey - aexdey;
	ET aexcey = aex * cey;
	ET cexaey = cex * aey;
	ET ac     = aexcey - cexaey;
	ET bexdey = bex * dey;
	ET dexbey = dex * bey;
	ET bd     = bexdey - dexbey;
	ET abc1   = aez * bc;
	ET abc2   = bez * ac;
	ET abc3   = cez * ab;
	ET abc4   = abc1 + abc3;
	ET abc    = abc4 - abc2;
	ET bcd1   = bez * cd;
	ET bcd2   = cez * bd;
	ET bcd3   = dez * bc;
	ET bcd4   = bcd1 + bcd3;
	ET bcd    = bcd4 - bcd2;
	ET cda1   = cez * da;
	ET cda2   = dez * ac;
	ET cda3   = aez * cd;
	ET cda4   = cda1 + cda3;
	ET cda    = cda4 + cda2;
	ET dab1   = dez * ab;
	ET dab2   = aez * bd;
	ET dab3   = bez * da;
	ET dab4   = dab1 + dab3;
	ET dab    = dab4 + dab2;
	ET al1    = aex * aex;
	ET al2    = aey * aey;
	ET al3    = aez * aez;
	ET al4    = al1 + al2;
	ET alift  = al4 + al3;
	ET bl1    = bex * bex;
	ET bl2    = bey * bey;
	ET bl3    = bez * bez;
	ET bl4    = bl1 + bl2;
	ET blift  = bl4 + bl3;
	ET cl1    = cex * cex;
	ET cl2    = cey * cey;
	ET cl3    = cez * cez;
	ET cl4    = cl1 + cl2;
	ET clift  = cl4 + cl3;
	ET dl1    = dex * dex;
	ET dl2    = dey * dey;
	ET dl3    = dez * dez;
	ET dl4    = dl1 + dl2;
	ET dlift  = dl4 + dl3;
	ET ds1    = dlift * abc;
	ET ds2    = clift * dab;
	ET dlp    = ds2 - ds1;
	ET dl     = dlp * d1;
	ET dr1p   = blift * cda;
	ET dr1    = dr1p * d1;
	ET dr2    = alift * bcd;
	ET dr     = dr2 - dr1;
	ET det    = dl + dr;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IEEEE_expansion(const GenericPoint3T<IT, ET> &p1, double pbx,
                              double pby, double pbz, double pcx, double pcy,
                              double pcz, double pdx, double pdy, double pdz,
                              double pex, double pey, double pez)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[16], *l1x = l1x_p, l1y_p[16], *l1y = l1y_p, l1z_p[16],
	                  *l1z = l1z_p, d1_p[16], *d1 = d1_p;
	int l1x_len = 16, l1y_len = 16, l1z_len = 16, d1_len = 16;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          pexd_p[16], *pexd = pexd_p;
		int    pexd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pex, &pexd, 16);
		double peyd_p[16], *peyd = peyd_p;
		int    peyd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pey, &peyd, 16);
		double pezd_p[16], *pezd = pezd_p;
		int    pezd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pez, &pezd, 16);
		double aex_p[16], *aex = aex_p;
		int    aex_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pexd_len, pexd, &aex, 16);
		double aey_p[16], *aey = aey_p;
		int    aey_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, peyd_len, peyd, &aey, 16);
		double aez_p[16], *aez = aez_p;
		int    aez_len =
		  o.Gen_Diff_With_PreAlloc(l1z_len, l1z, pezd_len, pezd, &aez, 16);
		double bex[2];
		o.two_Diff(pbx, pex, bex);
		double bey[2];
		o.two_Diff(pby, pey, bey);
		double bez[2];
		o.two_Diff(pbz, pez, bez);
		double cex[2];
		o.two_Diff(pcx, pex, cex);
		double cey[2];
		o.two_Diff(pcy, pey, cey);
		double cez[2];
		o.two_Diff(pcz, pez, cez);
		double dex[2];
		o.two_Diff(pdx, pex, dex);
		double dey[2];
		o.two_Diff(pdy, pey, dey);
		double dez[2];
		o.two_Diff(pdz, pez, dez);
		double aexbey_p[16], *aexbey = aexbey_p;
		int    aexbey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, 2, bey, &aexbey, 16);
		double bexaey_p[16], *bexaey = bexaey_p;
		int    bexaey_len =
		  o.Gen_Product_With_PreAlloc(2, bex, aey_len, aey, &bexaey, 16);
		double ab_p[16], *ab = ab_p;
		int    ab_len =
		  o.Gen_Diff_With_PreAlloc(aexbey_len, aexbey, bexaey_len, bexaey, &ab, 16);
		double bexcey[8];
		int    bexcey_len = o.Gen_Product(2, bex, 2, cey, bexcey);
		double cexbey[8];
		int    cexbey_len = o.Gen_Product(2, cex, 2, bey, cexbey);
		double bc[16];
		int    bc_len = o.Gen_Diff(bexcey_len, bexcey, cexbey_len, cexbey, bc);
		double cexdey[8];
		int    cexdey_len = o.Gen_Product(2, cex, 2, dey, cexdey);
		double dexcey[8];
		int    dexcey_len = o.Gen_Product(2, dex, 2, cey, dexcey);
		double cd[16];
		int    cd_len = o.Gen_Diff(cexdey_len, cexdey, dexcey_len, dexcey, cd);
		double dexaey_p[16], *dexaey = dexaey_p;
		int    dexaey_len =
		  o.Gen_Product_With_PreAlloc(2, dex, aey_len, aey, &dexaey, 16);
		double aexdey_p[16], *aexdey = aexdey_p;
		int    aexdey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, 2, dey, &aexdey, 16);
		double da_p[16], *da = da_p;
		int    da_len =
		  o.Gen_Diff_With_PreAlloc(dexaey_len, dexaey, aexdey_len, aexdey, &da, 16);
		double aexcey_p[16], *aexcey = aexcey_p;
		int    aexcey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, 2, cey, &aexcey, 16);
		double cexaey_p[16], *cexaey = cexaey_p;
		int    cexaey_len =
		  o.Gen_Product_With_PreAlloc(2, cex, aey_len, aey, &cexaey, 16);
		double ac_p[16], *ac = ac_p;
		int    ac_len =
		  o.Gen_Diff_With_PreAlloc(aexcey_len, aexcey, cexaey_len, cexaey, &ac, 16);
		double bexdey[8];
		int    bexdey_len = o.Gen_Product(2, bex, 2, dey, bexdey);
		double dexbey[8];
		int    dexbey_len = o.Gen_Product(2, dex, 2, bey, dexbey);
		double bd[16];
		int    bd_len = o.Gen_Diff(bexdey_len, bexdey, dexbey_len, dexbey, bd);
		double abc1_p[16], *abc1 = abc1_p;
		int    abc1_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bc_len, bc, &abc1, 16);
		double abc2_p[16], *abc2 = abc2_p;
		int abc2_len = o.Gen_Product_With_PreAlloc(2, bez, ac_len, ac, &abc2, 16);
		double abc3_p[16], *abc3 = abc3_p;
		int abc3_len = o.Gen_Product_With_PreAlloc(2, cez, ab_len, ab, &abc3, 16);
		double abc4_p[16], *abc4 = abc4_p;
		int    abc4_len =
		  o.Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 16);
		double abc_p[16], *abc = abc_p;
		int    abc_len =
		  o.Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 16);
		double bcd1_p[16], *bcd1 = bcd1_p;
		int bcd1_len = o.Gen_Product_With_PreAlloc(2, bez, cd_len, cd, &bcd1, 16);
		double bcd2_p[16], *bcd2 = bcd2_p;
		int bcd2_len = o.Gen_Product_With_PreAlloc(2, cez, bd_len, bd, &bcd2, 16);
		double bcd3_p[16], *bcd3 = bcd3_p;
		int bcd3_len = o.Gen_Product_With_PreAlloc(2, dez, bc_len, bc, &bcd3, 16);
		double bcd4_p[16], *bcd4 = bcd4_p;
		int    bcd4_len =
		  o.Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 16);
		double bcd_p[16], *bcd = bcd_p;
		int    bcd_len =
		  o.Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 16);
		double cda1_p[16], *cda1 = cda1_p;
		int cda1_len = o.Gen_Product_With_PreAlloc(2, cez, da_len, da, &cda1, 16);
		double cda2_p[16], *cda2 = cda2_p;
		int cda2_len = o.Gen_Product_With_PreAlloc(2, dez, ac_len, ac, &cda2, 16);
		double cda3_p[16], *cda3 = cda3_p;
		int    cda3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, cd_len, cd, &cda3, 16);
		double cda4_p[16], *cda4 = cda4_p;
		int    cda4_len =
		  o.Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 16);
		double cda_p[16], *cda = cda_p;
		int    cda_len =
		  o.Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 16);
		double dab1_p[16], *dab1 = dab1_p;
		int dab1_len = o.Gen_Product_With_PreAlloc(2, dez, ab_len, ab, &dab1, 16);
		double dab2_p[16], *dab2 = dab2_p;
		int    dab2_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bd_len, bd, &dab2, 16);
		double dab3_p[16], *dab3 = dab3_p;
		int dab3_len = o.Gen_Product_With_PreAlloc(2, bez, da_len, da, &dab3, 16);
		double dab4_p[16], *dab4 = dab4_p;
		int    dab4_len =
		  o.Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 16);
		double dab_p[16], *dab = dab_p;
		int    dab_len =
		  o.Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 16);
		double al1_p[16], *al1 = al1_p;
		int    al1_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, aex_len, aex, &al1, 16);
		double al2_p[16], *al2 = al2_p;
		int    al2_len =
		  o.Gen_Product_With_PreAlloc(aey_len, aey, aey_len, aey, &al2, 16);
		double al3_p[16], *al3 = al3_p;
		int    al3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, aez_len, aez, &al3, 16);
		double al4_p[16], *al4 = al4_p;
		int al4_len = o.Gen_Sum_With_PreAlloc(al1_len, al1, al2_len, al2, &al4, 16);
		double alift_p[16], *alift = alift_p;
		int    alift_len =
		  o.Gen_Sum_With_PreAlloc(al4_len, al4, al3_len, al3, &alift, 16);
		double bl1[8];
		int    bl1_len = o.Gen_Product(2, bex, 2, bex, bl1);
		double bl2[8];
		int    bl2_len = o.Gen_Product(2, bey, 2, bey, bl2);
		double bl3[8];
		int    bl3_len = o.Gen_Product(2, bez, 2, bez, bl3);
		double bl4[16];
		int    bl4_len = o.Gen_Sum(bl1_len, bl1, bl2_len, bl2, bl4);
		double blift_p[16], *blift = blift_p;
		int    blift_len =
		  o.Gen_Sum_With_PreAlloc(bl4_len, bl4, bl3_len, bl3, &blift, 16);
		double cl1[8];
		int    cl1_len = o.Gen_Product(2, cex, 2, cex, cl1);
		double cl2[8];
		int    cl2_len = o.Gen_Product(2, cey, 2, cey, cl2);
		double cl3[8];
		int    cl3_len = o.Gen_Product(2, cez, 2, cez, cl3);
		double cl4[16];
		int    cl4_len = o.Gen_Sum(cl1_len, cl1, cl2_len, cl2, cl4);
		double clift_p[16], *clift = clift_p;
		int    clift_len =
		  o.Gen_Sum_With_PreAlloc(cl4_len, cl4, cl3_len, cl3, &clift, 16);
		double dl1[8];
		int    dl1_len = o.Gen_Product(2, dex, 2, dex, dl1);
		double dl2[8];
		int    dl2_len = o.Gen_Product(2, dey, 2, dey, dl2);
		double dl3[8];
		int    dl3_len = o.Gen_Product(2, dez, 2, dez, dl3);
		double dl4[16];
		int    dl4_len = o.Gen_Sum(dl1_len, dl1, dl2_len, dl2, dl4);
		double dlift_p[16], *dlift = dlift_p;
		int    dlift_len =
		  o.Gen_Sum_With_PreAlloc(dl4_len, dl4, dl3_len, dl3, &dlift, 16);
		double ds1_p[16], *ds1 = ds1_p;
		int    ds1_len =
		  o.Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 16);
		double ds2_p[16], *ds2 = ds2_p;
		int    ds2_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 16);
		double dlp_p[16], *dlp = dlp_p;
		int    dlp_len =
		  o.Gen_Diff_With_PreAlloc(ds2_len, ds2, ds1_len, ds1, &dlp, 16);
		double dl_p[16], *dl = dl_p;
		int dl_len = o.Gen_Product_With_PreAlloc(dlp_len, dlp, d1_len, d1, &dl, 16);
		double dr1p_p[16], *dr1p = dr1p_p;
		int    dr1p_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1p, 16);
		double dr1_p[16], *dr1 = dr1_p;
		int    dr1_len =
		  o.Gen_Product_With_PreAlloc(dr1p_len, dr1p, d1_len, d1, &dr1, 16);
		double dr2_p[16], *dr2 = dr2_p;
		int    dr2_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 16);
		double dr_p[16], *dr = dr_p;
		int dr_len = o.Gen_Diff_With_PreAlloc(dr2_len, dr2, dr1_len, dr1, &dr, 16);
		double det_p[16], *det = det_p;
		int    det_len = o.Gen_Sum_With_PreAlloc(dl_len, dl, dr_len, dr, &det, 16);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (dr_p != dr)
			FreeDoubles(dr);
		if (dr2_p != dr2)
			FreeDoubles(dr2);
		if (dr1_p != dr1)
			FreeDoubles(dr1);
		if (dr1p_p != dr1p)
			FreeDoubles(dr1p);
		if (dl_p != dl)
			FreeDoubles(dl);
		if (dlp_p != dlp)
			FreeDoubles(dlp);
		if (ds2_p != ds2)
			FreeDoubles(ds2);
		if (ds1_p != ds1)
			FreeDoubles(ds1);
		if (dlift_p != dlift)
			FreeDoubles(dlift);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (al4_p != al4)
			FreeDoubles(al4);
		if (al3_p != al3)
			FreeDoubles(al3);
		if (al2_p != al2)
			FreeDoubles(al2);
		if (al1_p != al1)
			FreeDoubles(al1);
		if (dab_p != dab)
			FreeDoubles(dab);
		if (dab4_p != dab4)
			FreeDoubles(dab4);
		if (dab3_p != dab3)
			FreeDoubles(dab3);
		if (dab2_p != dab2)
			FreeDoubles(dab2);
		if (dab1_p != dab1)
			FreeDoubles(dab1);
		if (cda_p != cda)
			FreeDoubles(cda);
		if (cda4_p != cda4)
			FreeDoubles(cda4);
		if (cda3_p != cda3)
			FreeDoubles(cda3);
		if (cda2_p != cda2)
			FreeDoubles(cda2);
		if (cda1_p != cda1)
			FreeDoubles(cda1);
		if (bcd_p != bcd)
			FreeDoubles(bcd);
		if (bcd4_p != bcd4)
			FreeDoubles(bcd4);
		if (bcd3_p != bcd3)
			FreeDoubles(bcd3);
		if (bcd2_p != bcd2)
			FreeDoubles(bcd2);
		if (bcd1_p != bcd1)
			FreeDoubles(bcd1);
		if (abc_p != abc)
			FreeDoubles(abc);
		if (abc4_p != abc4)
			FreeDoubles(abc4);
		if (abc3_p != abc3)
			FreeDoubles(abc3);
		if (abc2_p != abc2)
			FreeDoubles(abc2);
		if (abc1_p != abc1)
			FreeDoubles(abc1);
		if (ac_p != ac)
			FreeDoubles(ac);
		if (cexaey_p != cexaey)
			FreeDoubles(cexaey);
		if (aexcey_p != aexcey)
			FreeDoubles(aexcey);
		if (da_p != da)
			FreeDoubles(da);
		if (aexdey_p != aexdey)
			FreeDoubles(aexdey);
		if (dexaey_p != dexaey)
			FreeDoubles(dexaey);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (bexaey_p != bexaey)
			FreeDoubles(bexaey);
		if (aexbey_p != aexbey)
			FreeDoubles(aexbey);
		if (aez_p != aez)
			FreeDoubles(aez);
		if (aey_p != aey)
			FreeDoubles(aey);
		if (aex_p != aex)
			FreeDoubles(aex);
		if (pezd_p != pezd)
			FreeDoubles(pezd);
		if (peyd_p != peyd)
			FreeDoubles(peyd);
		if (pexd_p != pexd)
			FreeDoubles(pexd);
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
		return inSphere_IEEEE_exact<IT, ET>(p1, pbx, pby, pbz, pcx, pcy, pcz, pdx,
		                                    pdy, pdz, pex, pey, pez);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inSphere_IEEEE(const GenericPoint3T<IT, ET> &p1, double pbx, double pby,
                    double pbz, double pcx, double pcy, double pcz, double pdx,
                    double pdy, double pdz, double pex, double pey, double pez)
{
	Sign ret;
	ret = inSphere_IEEEE_interval<IT, ET>(p1, pbx, pby, pbz, pcx, pcy, pcz, pdx,
	                                      pdy, pdz, pex, pey, pez);
	if (is_sign_reliable(ret))
		return ret;
	return inSphere_IEEEE_expansion<IT, ET>(p1, pbx, pby, pbz, pcx, pcy, pcz, pdx,
	                                        pdy, pdz, pex, pey, pez);
}

template <typename IT, typename ET>
Sign inSphere_IEEEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &pb,
                    const GenericPoint3T<IT, ET> &pc,
                    const GenericPoint3T<IT, ET> &pd,
                    const GenericPoint3T<IT, ET> &pe)
{
	return inSphere_IEEEE<IT, ET>(p1, pb.x(), pb.y(), pb.z(), pc.x(), pc.y(),
	                              pc.z(), pd.x(), pd.y(), pd.z(), pe.x(), pe.y(),
	                              pe.z());
}

template <typename IT, typename ET>
Sign inSphere_IIEEE_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, IT pcx, IT pcy,
                             IT pcz, IT pdx, IT pdy, IT pdz, IT pex, IT pey,
                             IT pez)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pexd   = pex * d1;
	IT peyd   = pey * d1;
	IT pezd   = pez * d1;
	IT aex    = l1x - pexd;
	IT aey    = l1y - peyd;
	IT aez    = l1z - pezd;
	IT pexd2  = pex * d2;
	IT peyd2  = pey * d2;
	IT pezd2  = pez * d2;
	IT bex    = l2x - pexd2;
	IT bey    = l2y - peyd2;
	IT bez    = l2z - pezd2;
	IT cex    = pcx - pex;
	IT cey    = pcy - pey;
	IT cez    = pcz - pez;
	IT dex    = pdx - pex;
	IT dey    = pdy - pey;
	IT dez    = pdz - pez;
	IT aexbey = aex * bey;
	IT bexaey = bex * aey;
	IT ab     = aexbey - bexaey;
	IT bexcey = bex * cey;
	IT cexbey = cex * bey;
	IT bc     = bexcey - cexbey;
	IT cexdey = cex * dey;
	IT dexcey = dex * cey;
	IT cd     = cexdey - dexcey;
	IT dexaey = dex * aey;
	IT aexdey = aex * dey;
	IT da     = dexaey - aexdey;
	IT aexcey = aex * cey;
	IT cexaey = cex * aey;
	IT ac     = aexcey - cexaey;
	IT bexdey = bex * dey;
	IT dexbey = dex * bey;
	IT bd     = bexdey - dexbey;
	IT abc1   = aez * bc;
	IT abc2   = bez * ac;
	IT abc3   = cez * ab;
	IT abc4   = abc1 + abc3;
	IT abc    = abc4 - abc2;
	IT bcd1   = bez * cd;
	IT bcd2   = cez * bd;
	IT bcd3   = dez * bc;
	IT bcd4   = bcd1 + bcd3;
	IT bcd    = bcd4 - bcd2;
	IT cda1   = cez * da;
	IT cda2   = dez * ac;
	IT cda3   = aez * cd;
	IT cda4   = cda1 + cda3;
	IT cda    = cda4 + cda2;
	IT dab1   = dez * ab;
	IT dab2   = aez * bd;
	IT dab3   = bez * da;
	IT dab4   = dab1 + dab3;
	IT dab    = dab4 + dab2;
	IT al1    = aex * aex;
	IT al2    = aey * aey;
	IT al3    = aez * aez;
	IT al4    = al1 + al2;
	IT alift  = al4 + al3;
	IT bl1    = bex * bex;
	IT bl2    = bey * bey;
	IT bl3    = bez * bez;
	IT bl4    = bl1 + bl2;
	IT blift  = bl4 + bl3;
	IT cl1    = cex * cex;
	IT cl2    = cey * cey;
	IT cl3    = cez * cez;
	IT cl4    = cl1 + cl2;
	IT clift  = cl4 + cl3;
	IT dl1    = dex * dex;
	IT dl2    = dey * dey;
	IT dl3    = dez * dez;
	IT dl4    = dl1 + dl2;
	IT dlift  = dl4 + dl3;
	IT ds1    = dlift * abc;
	IT ds2    = clift * dab;
	IT dl     = ds2 - ds1;
	IT dll    = dl * d1;
	IT dlll   = dll * d2;
	IT dr1    = blift * cda;
	IT dr12   = dr1 * d1;
	IT dr2    = alift * bcd;
	IT dr22   = dr2 * d2;
	IT dr     = dr22 - dr12;
	IT det    = dlll + dr;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIEEE_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2, ET pcx, ET pcy,
                          ET pcz, ET pdx, ET pdy, ET pdz, ET pex, ET pey,
                          ET pez)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET pexd   = pex * d1;
	ET peyd   = pey * d1;
	ET pezd   = pez * d1;
	ET aex    = l1x - pexd;
	ET aey    = l1y - peyd;
	ET aez    = l1z - pezd;
	ET pexd2  = pex * d2;
	ET peyd2  = pey * d2;
	ET pezd2  = pez * d2;
	ET bex    = l2x - pexd2;
	ET bey    = l2y - peyd2;
	ET bez    = l2z - pezd2;
	ET cex    = pcx - pex;
	ET cey    = pcy - pey;
	ET cez    = pcz - pez;
	ET dex    = pdx - pex;
	ET dey    = pdy - pey;
	ET dez    = pdz - pez;
	ET aexbey = aex * bey;
	ET bexaey = bex * aey;
	ET ab     = aexbey - bexaey;
	ET bexcey = bex * cey;
	ET cexbey = cex * bey;
	ET bc     = bexcey - cexbey;
	ET cexdey = cex * dey;
	ET dexcey = dex * cey;
	ET cd     = cexdey - dexcey;
	ET dexaey = dex * aey;
	ET aexdey = aex * dey;
	ET da     = dexaey - aexdey;
	ET aexcey = aex * cey;
	ET cexaey = cex * aey;
	ET ac     = aexcey - cexaey;
	ET bexdey = bex * dey;
	ET dexbey = dex * bey;
	ET bd     = bexdey - dexbey;
	ET abc1   = aez * bc;
	ET abc2   = bez * ac;
	ET abc3   = cez * ab;
	ET abc4   = abc1 + abc3;
	ET abc    = abc4 - abc2;
	ET bcd1   = bez * cd;
	ET bcd2   = cez * bd;
	ET bcd3   = dez * bc;
	ET bcd4   = bcd1 + bcd3;
	ET bcd    = bcd4 - bcd2;
	ET cda1   = cez * da;
	ET cda2   = dez * ac;
	ET cda3   = aez * cd;
	ET cda4   = cda1 + cda3;
	ET cda    = cda4 + cda2;
	ET dab1   = dez * ab;
	ET dab2   = aez * bd;
	ET dab3   = bez * da;
	ET dab4   = dab1 + dab3;
	ET dab    = dab4 + dab2;
	ET al1    = aex * aex;
	ET al2    = aey * aey;
	ET al3    = aez * aez;
	ET al4    = al1 + al2;
	ET alift  = al4 + al3;
	ET bl1    = bex * bex;
	ET bl2    = bey * bey;
	ET bl3    = bez * bez;
	ET bl4    = bl1 + bl2;
	ET blift  = bl4 + bl3;
	ET cl1    = cex * cex;
	ET cl2    = cey * cey;
	ET cl3    = cez * cez;
	ET cl4    = cl1 + cl2;
	ET clift  = cl4 + cl3;
	ET dl1    = dex * dex;
	ET dl2    = dey * dey;
	ET dl3    = dez * dez;
	ET dl4    = dl1 + dl2;
	ET dlift  = dl4 + dl3;
	ET ds1    = dlift * abc;
	ET ds2    = clift * dab;
	ET dl     = ds2 - ds1;
	ET dll    = dl * d1;
	ET dlll   = dll * d2;
	ET dr1    = blift * cda;
	ET dr12   = dr1 * d1;
	ET dr2    = alift * bcd;
	ET dr22   = dr2 * d2;
	ET dr     = dr22 - dr12;
	ET det    = dlll + dr;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIEEE_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2, double pcx,
                              double pcy, double pcz, double pdx, double pdy,
                              double pdz, double pex, double pey, double pez)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[16], *l1x = l1x_p, l1y_p[16], *l1y = l1y_p, l1z_p[16],
	                  *l1z = l1z_p, d1_p[16], *d1 = d1_p, l2x_p[16], *l2x = l2x_p,
	                  l2y_p[16], *l2y = l2y_p, l2z_p[16], *l2z = l2z_p, d2_p[16],
	                  *d2 = d2_p;
	int l1x_len = 16, l1y_len = 16, l1z_len = 16, d1_len = 16, l2x_len = 16,
	    l2y_len = 16, l2z_len = 16, d2_len = 16;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          pexd_p[16], *pexd = pexd_p;
		int    pexd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pex, &pexd, 16);
		double peyd_p[16], *peyd = peyd_p;
		int    peyd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pey, &peyd, 16);
		double pezd_p[16], *pezd = pezd_p;
		int    pezd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pez, &pezd, 16);
		double aex_p[16], *aex = aex_p;
		int    aex_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pexd_len, pexd, &aex, 16);
		double aey_p[16], *aey = aey_p;
		int    aey_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, peyd_len, peyd, &aey, 16);
		double aez_p[16], *aez = aez_p;
		int    aez_len =
		  o.Gen_Diff_With_PreAlloc(l1z_len, l1z, pezd_len, pezd, &aez, 16);
		double pexd2_p[16], *pexd2 = pexd2_p;
		int    pexd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pex, &pexd2, 16);
		double peyd2_p[16], *peyd2 = peyd2_p;
		int    peyd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pey, &peyd2, 16);
		double pezd2_p[16], *pezd2 = pezd2_p;
		int    pezd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pez, &pezd2, 16);
		double bex_p[16], *bex = bex_p;
		int    bex_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pexd2_len, pexd2, &bex, 16);
		double bey_p[16], *bey = bey_p;
		int    bey_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, peyd2_len, peyd2, &bey, 16);
		double bez_p[16], *bez = bez_p;
		int    bez_len =
		  o.Gen_Diff_With_PreAlloc(l2z_len, l2z, pezd2_len, pezd2, &bez, 16);
		double cex[2];
		o.two_Diff(pcx, pex, cex);
		double cey[2];
		o.two_Diff(pcy, pey, cey);
		double cez[2];
		o.two_Diff(pcz, pez, cez);
		double dex[2];
		o.two_Diff(pdx, pex, dex);
		double dey[2];
		o.two_Diff(pdy, pey, dey);
		double dez[2];
		o.two_Diff(pdz, pez, dez);
		double aexbey_p[16], *aexbey = aexbey_p;
		int    aexbey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, bey_len, bey, &aexbey, 16);
		double bexaey_p[16], *bexaey = bexaey_p;
		int    bexaey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, aey_len, aey, &bexaey, 16);
		double ab_p[16], *ab = ab_p;
		int    ab_len =
		  o.Gen_Diff_With_PreAlloc(aexbey_len, aexbey, bexaey_len, bexaey, &ab, 16);
		double bexcey_p[16], *bexcey = bexcey_p;
		int    bexcey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, 2, cey, &bexcey, 16);
		double cexbey_p[16], *cexbey = cexbey_p;
		int    cexbey_len =
		  o.Gen_Product_With_PreAlloc(2, cex, bey_len, bey, &cexbey, 16);
		double bc_p[16], *bc = bc_p;
		int    bc_len =
		  o.Gen_Diff_With_PreAlloc(bexcey_len, bexcey, cexbey_len, cexbey, &bc, 16);
		double cexdey[8];
		int    cexdey_len = o.Gen_Product(2, cex, 2, dey, cexdey);
		double dexcey[8];
		int    dexcey_len = o.Gen_Product(2, dex, 2, cey, dexcey);
		double cd[16];
		int    cd_len = o.Gen_Diff(cexdey_len, cexdey, dexcey_len, dexcey, cd);
		double dexaey_p[16], *dexaey = dexaey_p;
		int    dexaey_len =
		  o.Gen_Product_With_PreAlloc(2, dex, aey_len, aey, &dexaey, 16);
		double aexdey_p[16], *aexdey = aexdey_p;
		int    aexdey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, 2, dey, &aexdey, 16);
		double da_p[16], *da = da_p;
		int    da_len =
		  o.Gen_Diff_With_PreAlloc(dexaey_len, dexaey, aexdey_len, aexdey, &da, 16);
		double aexcey_p[16], *aexcey = aexcey_p;
		int    aexcey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, 2, cey, &aexcey, 16);
		double cexaey_p[16], *cexaey = cexaey_p;
		int    cexaey_len =
		  o.Gen_Product_With_PreAlloc(2, cex, aey_len, aey, &cexaey, 16);
		double ac_p[16], *ac = ac_p;
		int    ac_len =
		  o.Gen_Diff_With_PreAlloc(aexcey_len, aexcey, cexaey_len, cexaey, &ac, 16);
		double bexdey_p[16], *bexdey = bexdey_p;
		int    bexdey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, 2, dey, &bexdey, 16);
		double dexbey_p[16], *dexbey = dexbey_p;
		int    dexbey_len =
		  o.Gen_Product_With_PreAlloc(2, dex, bey_len, bey, &dexbey, 16);
		double bd_p[16], *bd = bd_p;
		int    bd_len =
		  o.Gen_Diff_With_PreAlloc(bexdey_len, bexdey, dexbey_len, dexbey, &bd, 16);
		double abc1_p[16], *abc1 = abc1_p;
		int    abc1_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bc_len, bc, &abc1, 16);
		double abc2_p[16], *abc2 = abc2_p;
		int    abc2_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, ac_len, ac, &abc2, 16);
		double abc3_p[16], *abc3 = abc3_p;
		int abc3_len = o.Gen_Product_With_PreAlloc(2, cez, ab_len, ab, &abc3, 16);
		double abc4_p[16], *abc4 = abc4_p;
		int    abc4_len =
		  o.Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 16);
		double abc_p[16], *abc = abc_p;
		int    abc_len =
		  o.Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 16);
		double bcd1_p[16], *bcd1 = bcd1_p;
		int    bcd1_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, cd_len, cd, &bcd1, 16);
		double bcd2_p[16], *bcd2 = bcd2_p;
		int bcd2_len = o.Gen_Product_With_PreAlloc(2, cez, bd_len, bd, &bcd2, 16);
		double bcd3_p[16], *bcd3 = bcd3_p;
		int bcd3_len = o.Gen_Product_With_PreAlloc(2, dez, bc_len, bc, &bcd3, 16);
		double bcd4_p[16], *bcd4 = bcd4_p;
		int    bcd4_len =
		  o.Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 16);
		double bcd_p[16], *bcd = bcd_p;
		int    bcd_len =
		  o.Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 16);
		double cda1_p[16], *cda1 = cda1_p;
		int cda1_len = o.Gen_Product_With_PreAlloc(2, cez, da_len, da, &cda1, 16);
		double cda2_p[16], *cda2 = cda2_p;
		int cda2_len = o.Gen_Product_With_PreAlloc(2, dez, ac_len, ac, &cda2, 16);
		double cda3_p[16], *cda3 = cda3_p;
		int    cda3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, cd_len, cd, &cda3, 16);
		double cda4_p[16], *cda4 = cda4_p;
		int    cda4_len =
		  o.Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 16);
		double cda_p[16], *cda = cda_p;
		int    cda_len =
		  o.Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 16);
		double dab1_p[16], *dab1 = dab1_p;
		int dab1_len = o.Gen_Product_With_PreAlloc(2, dez, ab_len, ab, &dab1, 16);
		double dab2_p[16], *dab2 = dab2_p;
		int    dab2_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bd_len, bd, &dab2, 16);
		double dab3_p[16], *dab3 = dab3_p;
		int    dab3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, da_len, da, &dab3, 16);
		double dab4_p[16], *dab4 = dab4_p;
		int    dab4_len =
		  o.Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 16);
		double dab_p[16], *dab = dab_p;
		int    dab_len =
		  o.Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 16);
		double al1_p[16], *al1 = al1_p;
		int    al1_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, aex_len, aex, &al1, 16);
		double al2_p[16], *al2 = al2_p;
		int    al2_len =
		  o.Gen_Product_With_PreAlloc(aey_len, aey, aey_len, aey, &al2, 16);
		double al3_p[16], *al3 = al3_p;
		int    al3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, aez_len, aez, &al3, 16);
		double al4_p[16], *al4 = al4_p;
		int al4_len = o.Gen_Sum_With_PreAlloc(al1_len, al1, al2_len, al2, &al4, 16);
		double alift_p[16], *alift = alift_p;
		int    alift_len =
		  o.Gen_Sum_With_PreAlloc(al4_len, al4, al3_len, al3, &alift, 16);
		double bl1_p[16], *bl1 = bl1_p;
		int    bl1_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, bex_len, bex, &bl1, 16);
		double bl2_p[16], *bl2 = bl2_p;
		int    bl2_len =
		  o.Gen_Product_With_PreAlloc(bey_len, bey, bey_len, bey, &bl2, 16);
		double bl3_p[16], *bl3 = bl3_p;
		int    bl3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, bez_len, bez, &bl3, 16);
		double bl4_p[16], *bl4 = bl4_p;
		int bl4_len = o.Gen_Sum_With_PreAlloc(bl1_len, bl1, bl2_len, bl2, &bl4, 16);
		double blift_p[16], *blift = blift_p;
		int    blift_len =
		  o.Gen_Sum_With_PreAlloc(bl4_len, bl4, bl3_len, bl3, &blift, 16);
		double cl1[8];
		int    cl1_len = o.Gen_Product(2, cex, 2, cex, cl1);
		double cl2[8];
		int    cl2_len = o.Gen_Product(2, cey, 2, cey, cl2);
		double cl3[8];
		int    cl3_len = o.Gen_Product(2, cez, 2, cez, cl3);
		double cl4[16];
		int    cl4_len = o.Gen_Sum(cl1_len, cl1, cl2_len, cl2, cl4);
		double clift_p[16], *clift = clift_p;
		int    clift_len =
		  o.Gen_Sum_With_PreAlloc(cl4_len, cl4, cl3_len, cl3, &clift, 16);
		double dl1[8];
		int    dl1_len = o.Gen_Product(2, dex, 2, dex, dl1);
		double dl2[8];
		int    dl2_len = o.Gen_Product(2, dey, 2, dey, dl2);
		double dl3[8];
		int    dl3_len = o.Gen_Product(2, dez, 2, dez, dl3);
		double dl4[16];
		int    dl4_len = o.Gen_Sum(dl1_len, dl1, dl2_len, dl2, dl4);
		double dlift_p[16], *dlift = dlift_p;
		int    dlift_len =
		  o.Gen_Sum_With_PreAlloc(dl4_len, dl4, dl3_len, dl3, &dlift, 16);
		double ds1_p[16], *ds1 = ds1_p;
		int    ds1_len =
		  o.Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 16);
		double ds2_p[16], *ds2 = ds2_p;
		int    ds2_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 16);
		double dl_p[16], *dl = dl_p;
		int dl_len = o.Gen_Diff_With_PreAlloc(ds2_len, ds2, ds1_len, ds1, &dl, 16);
		double dll_p[16], *dll = dll_p;
		int dll_len = o.Gen_Product_With_PreAlloc(dl_len, dl, d1_len, d1, &dll, 16);
		double dlll_p[16], *dlll = dlll_p;
		int    dlll_len =
		  o.Gen_Product_With_PreAlloc(dll_len, dll, d2_len, d2, &dlll, 16);
		double dr1_p[16], *dr1 = dr1_p;
		int    dr1_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1, 16);
		double dr12_p[16], *dr12 = dr12_p;
		int    dr12_len =
		  o.Gen_Product_With_PreAlloc(dr1_len, dr1, d1_len, d1, &dr12, 16);
		double dr2_p[16], *dr2 = dr2_p;
		int    dr2_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 16);
		double dr22_p[16], *dr22 = dr22_p;
		int    dr22_len =
		  o.Gen_Product_With_PreAlloc(dr2_len, dr2, d2_len, d2, &dr22, 16);
		double dr_p[16], *dr = dr_p;
		int    dr_len =
		  o.Gen_Diff_With_PreAlloc(dr22_len, dr22, dr12_len, dr12, &dr, 16);
		double det_p[16], *det = det_p;
		int det_len = o.Gen_Sum_With_PreAlloc(dlll_len, dlll, dr_len, dr, &det, 16);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (dr_p != dr)
			FreeDoubles(dr);
		if (dr22_p != dr22)
			FreeDoubles(dr22);
		if (dr2_p != dr2)
			FreeDoubles(dr2);
		if (dr12_p != dr12)
			FreeDoubles(dr12);
		if (dr1_p != dr1)
			FreeDoubles(dr1);
		if (dlll_p != dlll)
			FreeDoubles(dlll);
		if (dll_p != dll)
			FreeDoubles(dll);
		if (dl_p != dl)
			FreeDoubles(dl);
		if (ds2_p != ds2)
			FreeDoubles(ds2);
		if (ds1_p != ds1)
			FreeDoubles(ds1);
		if (dlift_p != dlift)
			FreeDoubles(dlift);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bl4_p != bl4)
			FreeDoubles(bl4);
		if (bl3_p != bl3)
			FreeDoubles(bl3);
		if (bl2_p != bl2)
			FreeDoubles(bl2);
		if (bl1_p != bl1)
			FreeDoubles(bl1);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (al4_p != al4)
			FreeDoubles(al4);
		if (al3_p != al3)
			FreeDoubles(al3);
		if (al2_p != al2)
			FreeDoubles(al2);
		if (al1_p != al1)
			FreeDoubles(al1);
		if (dab_p != dab)
			FreeDoubles(dab);
		if (dab4_p != dab4)
			FreeDoubles(dab4);
		if (dab3_p != dab3)
			FreeDoubles(dab3);
		if (dab2_p != dab2)
			FreeDoubles(dab2);
		if (dab1_p != dab1)
			FreeDoubles(dab1);
		if (cda_p != cda)
			FreeDoubles(cda);
		if (cda4_p != cda4)
			FreeDoubles(cda4);
		if (cda3_p != cda3)
			FreeDoubles(cda3);
		if (cda2_p != cda2)
			FreeDoubles(cda2);
		if (cda1_p != cda1)
			FreeDoubles(cda1);
		if (bcd_p != bcd)
			FreeDoubles(bcd);
		if (bcd4_p != bcd4)
			FreeDoubles(bcd4);
		if (bcd3_p != bcd3)
			FreeDoubles(bcd3);
		if (bcd2_p != bcd2)
			FreeDoubles(bcd2);
		if (bcd1_p != bcd1)
			FreeDoubles(bcd1);
		if (abc_p != abc)
			FreeDoubles(abc);
		if (abc4_p != abc4)
			FreeDoubles(abc4);
		if (abc3_p != abc3)
			FreeDoubles(abc3);
		if (abc2_p != abc2)
			FreeDoubles(abc2);
		if (abc1_p != abc1)
			FreeDoubles(abc1);
		if (bd_p != bd)
			FreeDoubles(bd);
		if (dexbey_p != dexbey)
			FreeDoubles(dexbey);
		if (bexdey_p != bexdey)
			FreeDoubles(bexdey);
		if (ac_p != ac)
			FreeDoubles(ac);
		if (cexaey_p != cexaey)
			FreeDoubles(cexaey);
		if (aexcey_p != aexcey)
			FreeDoubles(aexcey);
		if (da_p != da)
			FreeDoubles(da);
		if (aexdey_p != aexdey)
			FreeDoubles(aexdey);
		if (dexaey_p != dexaey)
			FreeDoubles(dexaey);
		if (bc_p != bc)
			FreeDoubles(bc);
		if (cexbey_p != cexbey)
			FreeDoubles(cexbey);
		if (bexcey_p != bexcey)
			FreeDoubles(bexcey);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (bexaey_p != bexaey)
			FreeDoubles(bexaey);
		if (aexbey_p != aexbey)
			FreeDoubles(aexbey);
		if (bez_p != bez)
			FreeDoubles(bez);
		if (bey_p != bey)
			FreeDoubles(bey);
		if (bex_p != bex)
			FreeDoubles(bex);
		if (pezd2_p != pezd2)
			FreeDoubles(pezd2);
		if (peyd2_p != peyd2)
			FreeDoubles(peyd2);
		if (pexd2_p != pexd2)
			FreeDoubles(pexd2);
		if (aez_p != aez)
			FreeDoubles(aez);
		if (aey_p != aey)
			FreeDoubles(aey);
		if (aex_p != aex)
			FreeDoubles(aex);
		if (pezd_p != pezd)
			FreeDoubles(pezd);
		if (peyd_p != peyd)
			FreeDoubles(peyd);
		if (pexd_p != pexd)
			FreeDoubles(pexd);
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
		return inSphere_IIEEE_exact<IT, ET>(p1, p2, pcx, pcy, pcz, pdx, pdy, pdz,
		                                    pex, pey, pez);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inSphere_IIEEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, double pcx, double pcy,
                    double pcz, double pdx, double pdy, double pdz, double pex,
                    double pey, double pez)
{
	Sign ret;
	ret = inSphere_IIEEE_interval<IT, ET>(p1, p2, pcx, pcy, pcz, pdx, pdy, pdz,
	                                      pex, pey, pez);
	if (is_sign_reliable(ret))
		return ret;
	return inSphere_IIEEE_expansion<IT, ET>(p1, p2, pcx, pcy, pcz, pdx, pdy, pdz,
	                                        pex, pey, pez);
}

template <typename IT, typename ET>
Sign inSphere_IIEEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &pc,
                    const GenericPoint3T<IT, ET> &pd,
                    const GenericPoint3T<IT, ET> &pe)
{
	return inSphere_IIEEE<IT, ET>(p1, p2, pc.x(), pc.y(), pc.z(), pd.x(), pd.y(),
	                              pd.z(), pe.x(), pe.y(), pe.z());
}

template <typename IT, typename ET>
Sign inSphere_IIIEE_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3, IT pdx, IT pdy,
                             IT pdz, IT pex, IT pey, IT pez)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pexd   = pex * d1;
	IT peyd   = pey * d1;
	IT pezd   = pez * d1;
	IT aex    = l1x - pexd;
	IT aey    = l1y - peyd;
	IT aez    = l1z - pezd;
	IT pexd2  = pex * d2;
	IT peyd2  = pey * d2;
	IT pezd2  = pez * d2;
	IT bex    = l2x - pexd2;
	IT bey    = l2y - peyd2;
	IT bez    = l2z - pezd2;
	IT pexd3  = pex * d3;
	IT peyd3  = pey * d3;
	IT pezd3  = pez * d3;
	IT cex    = l3x - pexd3;
	IT cey    = l3y - peyd3;
	IT cez    = l3z - pezd3;
	IT dex    = pdx - pex;
	IT dey    = pdy - pey;
	IT dez    = pdz - pez;
	IT aexbey = aex * bey;
	IT bexaey = bex * aey;
	IT ab     = aexbey - bexaey;
	IT bexcey = bex * cey;
	IT cexbey = cex * bey;
	IT bc     = bexcey - cexbey;
	IT cexdey = cex * dey;
	IT dexcey = dex * cey;
	IT cd     = cexdey - dexcey;
	IT dexaey = dex * aey;
	IT aexdey = aex * dey;
	IT da     = dexaey - aexdey;
	IT aexcey = aex * cey;
	IT cexaey = cex * aey;
	IT ac     = aexcey - cexaey;
	IT bexdey = bex * dey;
	IT dexbey = dex * bey;
	IT bd     = bexdey - dexbey;
	IT abc1   = aez * bc;
	IT abc2   = bez * ac;
	IT abc3   = cez * ab;
	IT abc4   = abc1 + abc3;
	IT abc    = abc4 - abc2;
	IT bcd1   = bez * cd;
	IT bcd2   = cez * bd;
	IT bcd3   = dez * bc;
	IT bcd4   = bcd1 + bcd3;
	IT bcd    = bcd4 - bcd2;
	IT cda1   = cez * da;
	IT cda2   = dez * ac;
	IT cda3   = aez * cd;
	IT cda4   = cda1 + cda3;
	IT cda    = cda4 + cda2;
	IT dab1   = dez * ab;
	IT dab2   = aez * bd;
	IT dab3   = bez * da;
	IT dab4   = dab1 + dab3;
	IT dab    = dab4 + dab2;
	IT al1    = aex * aex;
	IT al2    = aey * aey;
	IT al3    = aez * aez;
	IT al4    = al1 + al2;
	IT alift  = al4 + al3;
	IT bl1    = bex * bex;
	IT bl2    = bey * bey;
	IT bl3    = bez * bez;
	IT bl4    = bl1 + bl2;
	IT blift  = bl4 + bl3;
	IT cl1    = cex * cex;
	IT cl2    = cey * cey;
	IT cl3    = cez * cez;
	IT cl4    = cl1 + cl2;
	IT clift  = cl4 + cl3;
	IT dl1    = dex * dex;
	IT dl2    = dey * dey;
	IT dl3    = dez * dez;
	IT dl4    = dl1 + dl2;
	IT dlift  = dl4 + dl3;
	IT ds1    = dlift * abc;
	IT ds1n   = ds1 * d3;
	IT ds2    = clift * dab;
	IT dl     = ds2 - ds1n;
	IT dlm    = dl * d1;
	IT dln    = dlm * d2;
	IT dr1    = blift * cda;
	IT dr1n   = dr1 * d1;
	IT dr2    = alift * bcd;
	IT dr2n   = dr2 * d2;
	IT dr     = dr2n - dr1n;
	IT drn    = dr * d3;
	IT det    = dln + drn;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIIEE_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2,
                          const GenericPoint3T<IT, ET> &p3, ET pdx, ET pdy,
                          ET pdz, ET pex, ET pey, ET pez)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	ET pexd   = pex * d1;
	ET peyd   = pey * d1;
	ET pezd   = pez * d1;
	ET aex    = l1x - pexd;
	ET aey    = l1y - peyd;
	ET aez    = l1z - pezd;
	ET pexd2  = pex * d2;
	ET peyd2  = pey * d2;
	ET pezd2  = pez * d2;
	ET bex    = l2x - pexd2;
	ET bey    = l2y - peyd2;
	ET bez    = l2z - pezd2;
	ET pexd3  = pex * d3;
	ET peyd3  = pey * d3;
	ET pezd3  = pez * d3;
	ET cex    = l3x - pexd3;
	ET cey    = l3y - peyd3;
	ET cez    = l3z - pezd3;
	ET dex    = pdx - pex;
	ET dey    = pdy - pey;
	ET dez    = pdz - pez;
	ET aexbey = aex * bey;
	ET bexaey = bex * aey;
	ET ab     = aexbey - bexaey;
	ET bexcey = bex * cey;
	ET cexbey = cex * bey;
	ET bc     = bexcey - cexbey;
	ET cexdey = cex * dey;
	ET dexcey = dex * cey;
	ET cd     = cexdey - dexcey;
	ET dexaey = dex * aey;
	ET aexdey = aex * dey;
	ET da     = dexaey - aexdey;
	ET aexcey = aex * cey;
	ET cexaey = cex * aey;
	ET ac     = aexcey - cexaey;
	ET bexdey = bex * dey;
	ET dexbey = dex * bey;
	ET bd     = bexdey - dexbey;
	ET abc1   = aez * bc;
	ET abc2   = bez * ac;
	ET abc3   = cez * ab;
	ET abc4   = abc1 + abc3;
	ET abc    = abc4 - abc2;
	ET bcd1   = bez * cd;
	ET bcd2   = cez * bd;
	ET bcd3   = dez * bc;
	ET bcd4   = bcd1 + bcd3;
	ET bcd    = bcd4 - bcd2;
	ET cda1   = cez * da;
	ET cda2   = dez * ac;
	ET cda3   = aez * cd;
	ET cda4   = cda1 + cda3;
	ET cda    = cda4 + cda2;
	ET dab1   = dez * ab;
	ET dab2   = aez * bd;
	ET dab3   = bez * da;
	ET dab4   = dab1 + dab3;
	ET dab    = dab4 + dab2;
	ET al1    = aex * aex;
	ET al2    = aey * aey;
	ET al3    = aez * aez;
	ET al4    = al1 + al2;
	ET alift  = al4 + al3;
	ET bl1    = bex * bex;
	ET bl2    = bey * bey;
	ET bl3    = bez * bez;
	ET bl4    = bl1 + bl2;
	ET blift  = bl4 + bl3;
	ET cl1    = cex * cex;
	ET cl2    = cey * cey;
	ET cl3    = cez * cez;
	ET cl4    = cl1 + cl2;
	ET clift  = cl4 + cl3;
	ET dl1    = dex * dex;
	ET dl2    = dey * dey;
	ET dl3    = dez * dez;
	ET dl4    = dl1 + dl2;
	ET dlift  = dl4 + dl3;
	ET ds1    = dlift * abc;
	ET ds1n   = ds1 * d3;
	ET ds2    = clift * dab;
	ET dl     = ds2 - ds1n;
	ET dlm    = dl * d1;
	ET dln    = dlm * d2;
	ET dr1    = blift * cda;
	ET dr1n   = dr1 * d1;
	ET dr2    = alift * bcd;
	ET dr2n   = dr2 * d2;
	ET dr     = dr2n - dr1n;
	ET drn    = dr * d3;
	ET det    = dln + drn;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIIEE_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3, double pdx,
                              double pdy, double pdz, double pex, double pey,
                              double pez)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[16], *l1x = l1x_p, l1y_p[16], *l1y = l1y_p, l1z_p[16],
	                  *l1z = l1z_p, d1_p[16], *d1 = d1_p, l2x_p[16], *l2x = l2x_p,
	                  l2y_p[16], *l2y = l2y_p, l2z_p[16], *l2z = l2z_p, d2_p[16],
	                  *d2 = d2_p, l3x_p[16], *l3x = l3x_p, l3y_p[16],
	                  *l3y = l3y_p, l3z_p[16], *l3z = l3z_p, d3_p[16], *d3 = d3_p;
	int l1x_len = 16, l1y_len = 16, l1z_len = 16, d1_len = 16, l2x_len = 16,
	    l2y_len = 16, l2z_len = 16, d2_len = 16, l3x_len = 16, l3y_len = 16,
	    l3z_len = 16, d3_len = 16;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          pexd_p[16], *pexd = pexd_p;
		int    pexd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pex, &pexd, 16);
		double peyd_p[16], *peyd = peyd_p;
		int    peyd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pey, &peyd, 16);
		double pezd_p[16], *pezd = pezd_p;
		int    pezd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pez, &pezd, 16);
		double aex_p[16], *aex = aex_p;
		int    aex_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pexd_len, pexd, &aex, 16);
		double aey_p[16], *aey = aey_p;
		int    aey_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, peyd_len, peyd, &aey, 16);
		double aez_p[16], *aez = aez_p;
		int    aez_len =
		  o.Gen_Diff_With_PreAlloc(l1z_len, l1z, pezd_len, pezd, &aez, 16);
		double pexd2_p[16], *pexd2 = pexd2_p;
		int    pexd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pex, &pexd2, 16);
		double peyd2_p[16], *peyd2 = peyd2_p;
		int    peyd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pey, &peyd2, 16);
		double pezd2_p[16], *pezd2 = pezd2_p;
		int    pezd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pez, &pezd2, 16);
		double bex_p[16], *bex = bex_p;
		int    bex_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pexd2_len, pexd2, &bex, 16);
		double bey_p[16], *bey = bey_p;
		int    bey_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, peyd2_len, peyd2, &bey, 16);
		double bez_p[16], *bez = bez_p;
		int    bez_len =
		  o.Gen_Diff_With_PreAlloc(l2z_len, l2z, pezd2_len, pezd2, &bez, 16);
		double pexd3_p[16], *pexd3 = pexd3_p;
		int    pexd3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pex, &pexd3, 16);
		double peyd3_p[16], *peyd3 = peyd3_p;
		int    peyd3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pey, &peyd3, 16);
		double pezd3_p[16], *pezd3 = pezd3_p;
		int    pezd3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pez, &pezd3, 16);
		double cex_p[16], *cex = cex_p;
		int    cex_len =
		  o.Gen_Diff_With_PreAlloc(l3x_len, l3x, pexd3_len, pexd3, &cex, 16);
		double cey_p[16], *cey = cey_p;
		int    cey_len =
		  o.Gen_Diff_With_PreAlloc(l3y_len, l3y, peyd3_len, peyd3, &cey, 16);
		double cez_p[16], *cez = cez_p;
		int    cez_len =
		  o.Gen_Diff_With_PreAlloc(l3z_len, l3z, pezd3_len, pezd3, &cez, 16);
		double dex[2];
		o.two_Diff(pdx, pex, dex);
		double dey[2];
		o.two_Diff(pdy, pey, dey);
		double dez[2];
		o.two_Diff(pdz, pez, dez);
		double aexbey_p[16], *aexbey = aexbey_p;
		int    aexbey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, bey_len, bey, &aexbey, 16);
		double bexaey_p[16], *bexaey = bexaey_p;
		int    bexaey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, aey_len, aey, &bexaey, 16);
		double ab_p[16], *ab = ab_p;
		int    ab_len =
		  o.Gen_Diff_With_PreAlloc(aexbey_len, aexbey, bexaey_len, bexaey, &ab, 16);
		double bexcey_p[16], *bexcey = bexcey_p;
		int    bexcey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, cey_len, cey, &bexcey, 16);
		double cexbey_p[16], *cexbey = cexbey_p;
		int    cexbey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, bey_len, bey, &cexbey, 16);
		double bc_p[16], *bc = bc_p;
		int    bc_len =
		  o.Gen_Diff_With_PreAlloc(bexcey_len, bexcey, cexbey_len, cexbey, &bc, 16);
		double cexdey_p[16], *cexdey = cexdey_p;
		int    cexdey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, 2, dey, &cexdey, 16);
		double dexcey_p[16], *dexcey = dexcey_p;
		int    dexcey_len =
		  o.Gen_Product_With_PreAlloc(2, dex, cey_len, cey, &dexcey, 16);
		double cd_p[16], *cd = cd_p;
		int    cd_len =
		  o.Gen_Diff_With_PreAlloc(cexdey_len, cexdey, dexcey_len, dexcey, &cd, 16);
		double dexaey_p[16], *dexaey = dexaey_p;
		int    dexaey_len =
		  o.Gen_Product_With_PreAlloc(2, dex, aey_len, aey, &dexaey, 16);
		double aexdey_p[16], *aexdey = aexdey_p;
		int    aexdey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, 2, dey, &aexdey, 16);
		double da_p[16], *da = da_p;
		int    da_len =
		  o.Gen_Diff_With_PreAlloc(dexaey_len, dexaey, aexdey_len, aexdey, &da, 16);
		double aexcey_p[16], *aexcey = aexcey_p;
		int    aexcey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, cey_len, cey, &aexcey, 16);
		double cexaey_p[16], *cexaey = cexaey_p;
		int    cexaey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, aey_len, aey, &cexaey, 16);
		double ac_p[16], *ac = ac_p;
		int    ac_len =
		  o.Gen_Diff_With_PreAlloc(aexcey_len, aexcey, cexaey_len, cexaey, &ac, 16);
		double bexdey_p[16], *bexdey = bexdey_p;
		int    bexdey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, 2, dey, &bexdey, 16);
		double dexbey_p[16], *dexbey = dexbey_p;
		int    dexbey_len =
		  o.Gen_Product_With_PreAlloc(2, dex, bey_len, bey, &dexbey, 16);
		double bd_p[16], *bd = bd_p;
		int    bd_len =
		  o.Gen_Diff_With_PreAlloc(bexdey_len, bexdey, dexbey_len, dexbey, &bd, 16);
		double abc1_p[16], *abc1 = abc1_p;
		int    abc1_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bc_len, bc, &abc1, 16);
		double abc2_p[16], *abc2 = abc2_p;
		int    abc2_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, ac_len, ac, &abc2, 16);
		double abc3_p[16], *abc3 = abc3_p;
		int    abc3_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, ab_len, ab, &abc3, 16);
		double abc4_p[16], *abc4 = abc4_p;
		int    abc4_len =
		  o.Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 16);
		double abc_p[16], *abc = abc_p;
		int    abc_len =
		  o.Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 16);
		double bcd1_p[16], *bcd1 = bcd1_p;
		int    bcd1_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, cd_len, cd, &bcd1, 16);
		double bcd2_p[16], *bcd2 = bcd2_p;
		int    bcd2_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, bd_len, bd, &bcd2, 16);
		double bcd3_p[16], *bcd3 = bcd3_p;
		int bcd3_len = o.Gen_Product_With_PreAlloc(2, dez, bc_len, bc, &bcd3, 16);
		double bcd4_p[16], *bcd4 = bcd4_p;
		int    bcd4_len =
		  o.Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 16);
		double bcd_p[16], *bcd = bcd_p;
		int    bcd_len =
		  o.Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 16);
		double cda1_p[16], *cda1 = cda1_p;
		int    cda1_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, da_len, da, &cda1, 16);
		double cda2_p[16], *cda2 = cda2_p;
		int cda2_len = o.Gen_Product_With_PreAlloc(2, dez, ac_len, ac, &cda2, 16);
		double cda3_p[16], *cda3 = cda3_p;
		int    cda3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, cd_len, cd, &cda3, 16);
		double cda4_p[16], *cda4 = cda4_p;
		int    cda4_len =
		  o.Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 16);
		double cda_p[16], *cda = cda_p;
		int    cda_len =
		  o.Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 16);
		double dab1_p[16], *dab1 = dab1_p;
		int dab1_len = o.Gen_Product_With_PreAlloc(2, dez, ab_len, ab, &dab1, 16);
		double dab2_p[16], *dab2 = dab2_p;
		int    dab2_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bd_len, bd, &dab2, 16);
		double dab3_p[16], *dab3 = dab3_p;
		int    dab3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, da_len, da, &dab3, 16);
		double dab4_p[16], *dab4 = dab4_p;
		int    dab4_len =
		  o.Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 16);
		double dab_p[16], *dab = dab_p;
		int    dab_len =
		  o.Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 16);
		double al1_p[16], *al1 = al1_p;
		int    al1_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, aex_len, aex, &al1, 16);
		double al2_p[16], *al2 = al2_p;
		int    al2_len =
		  o.Gen_Product_With_PreAlloc(aey_len, aey, aey_len, aey, &al2, 16);
		double al3_p[16], *al3 = al3_p;
		int    al3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, aez_len, aez, &al3, 16);
		double al4_p[16], *al4 = al4_p;
		int al4_len = o.Gen_Sum_With_PreAlloc(al1_len, al1, al2_len, al2, &al4, 16);
		double alift_p[16], *alift = alift_p;
		int    alift_len =
		  o.Gen_Sum_With_PreAlloc(al4_len, al4, al3_len, al3, &alift, 16);
		double bl1_p[16], *bl1 = bl1_p;
		int    bl1_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, bex_len, bex, &bl1, 16);
		double bl2_p[16], *bl2 = bl2_p;
		int    bl2_len =
		  o.Gen_Product_With_PreAlloc(bey_len, bey, bey_len, bey, &bl2, 16);
		double bl3_p[16], *bl3 = bl3_p;
		int    bl3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, bez_len, bez, &bl3, 16);
		double bl4_p[16], *bl4 = bl4_p;
		int bl4_len = o.Gen_Sum_With_PreAlloc(bl1_len, bl1, bl2_len, bl2, &bl4, 16);
		double blift_p[16], *blift = blift_p;
		int    blift_len =
		  o.Gen_Sum_With_PreAlloc(bl4_len, bl4, bl3_len, bl3, &blift, 16);
		double cl1_p[16], *cl1 = cl1_p;
		int    cl1_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, cex_len, cex, &cl1, 16);
		double cl2_p[16], *cl2 = cl2_p;
		int    cl2_len =
		  o.Gen_Product_With_PreAlloc(cey_len, cey, cey_len, cey, &cl2, 16);
		double cl3_p[16], *cl3 = cl3_p;
		int    cl3_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, cez_len, cez, &cl3, 16);
		double cl4_p[16], *cl4 = cl4_p;
		int cl4_len = o.Gen_Sum_With_PreAlloc(cl1_len, cl1, cl2_len, cl2, &cl4, 16);
		double clift_p[16], *clift = clift_p;
		int    clift_len =
		  o.Gen_Sum_With_PreAlloc(cl4_len, cl4, cl3_len, cl3, &clift, 16);
		double dl1[8];
		int    dl1_len = o.Gen_Product(2, dex, 2, dex, dl1);
		double dl2[8];
		int    dl2_len = o.Gen_Product(2, dey, 2, dey, dl2);
		double dl3[8];
		int    dl3_len = o.Gen_Product(2, dez, 2, dez, dl3);
		double dl4[16];
		int    dl4_len = o.Gen_Sum(dl1_len, dl1, dl2_len, dl2, dl4);
		double dlift_p[16], *dlift = dlift_p;
		int    dlift_len =
		  o.Gen_Sum_With_PreAlloc(dl4_len, dl4, dl3_len, dl3, &dlift, 16);
		double ds1_p[16], *ds1 = ds1_p;
		int    ds1_len =
		  o.Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 16);
		double ds1n_p[16], *ds1n = ds1n_p;
		int    ds1n_len =
		  o.Gen_Product_With_PreAlloc(ds1_len, ds1, d3_len, d3, &ds1n, 16);
		double ds2_p[16], *ds2 = ds2_p;
		int    ds2_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 16);
		double dl_p[16], *dl = dl_p;
		int    dl_len =
		  o.Gen_Diff_With_PreAlloc(ds2_len, ds2, ds1n_len, ds1n, &dl, 16);
		double dlm_p[16], *dlm = dlm_p;
		int dlm_len = o.Gen_Product_With_PreAlloc(dl_len, dl, d1_len, d1, &dlm, 16);
		double dln_p[16], *dln = dln_p;
		int    dln_len =
		  o.Gen_Product_With_PreAlloc(dlm_len, dlm, d2_len, d2, &dln, 16);
		double dr1_p[16], *dr1 = dr1_p;
		int    dr1_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1, 16);
		double dr1n_p[16], *dr1n = dr1n_p;
		int    dr1n_len =
		  o.Gen_Product_With_PreAlloc(dr1_len, dr1, d1_len, d1, &dr1n, 16);
		double dr2_p[16], *dr2 = dr2_p;
		int    dr2_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 16);
		double dr2n_p[16], *dr2n = dr2n_p;
		int    dr2n_len =
		  o.Gen_Product_With_PreAlloc(dr2_len, dr2, d2_len, d2, &dr2n, 16);
		double dr_p[16], *dr = dr_p;
		int    dr_len =
		  o.Gen_Diff_With_PreAlloc(dr2n_len, dr2n, dr1n_len, dr1n, &dr, 16);
		double drn_p[16], *drn = drn_p;
		int drn_len = o.Gen_Product_With_PreAlloc(dr_len, dr, d3_len, d3, &drn, 16);
		double det_p[16], *det = det_p;
		int det_len = o.Gen_Sum_With_PreAlloc(dln_len, dln, drn_len, drn, &det, 16);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (drn_p != drn)
			FreeDoubles(drn);
		if (dr_p != dr)
			FreeDoubles(dr);
		if (dr2n_p != dr2n)
			FreeDoubles(dr2n);
		if (dr2_p != dr2)
			FreeDoubles(dr2);
		if (dr1n_p != dr1n)
			FreeDoubles(dr1n);
		if (dr1_p != dr1)
			FreeDoubles(dr1);
		if (dln_p != dln)
			FreeDoubles(dln);
		if (dlm_p != dlm)
			FreeDoubles(dlm);
		if (dl_p != dl)
			FreeDoubles(dl);
		if (ds2_p != ds2)
			FreeDoubles(ds2);
		if (ds1n_p != ds1n)
			FreeDoubles(ds1n);
		if (ds1_p != ds1)
			FreeDoubles(ds1);
		if (dlift_p != dlift)
			FreeDoubles(dlift);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (cl4_p != cl4)
			FreeDoubles(cl4);
		if (cl3_p != cl3)
			FreeDoubles(cl3);
		if (cl2_p != cl2)
			FreeDoubles(cl2);
		if (cl1_p != cl1)
			FreeDoubles(cl1);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bl4_p != bl4)
			FreeDoubles(bl4);
		if (bl3_p != bl3)
			FreeDoubles(bl3);
		if (bl2_p != bl2)
			FreeDoubles(bl2);
		if (bl1_p != bl1)
			FreeDoubles(bl1);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (al4_p != al4)
			FreeDoubles(al4);
		if (al3_p != al3)
			FreeDoubles(al3);
		if (al2_p != al2)
			FreeDoubles(al2);
		if (al1_p != al1)
			FreeDoubles(al1);
		if (dab_p != dab)
			FreeDoubles(dab);
		if (dab4_p != dab4)
			FreeDoubles(dab4);
		if (dab3_p != dab3)
			FreeDoubles(dab3);
		if (dab2_p != dab2)
			FreeDoubles(dab2);
		if (dab1_p != dab1)
			FreeDoubles(dab1);
		if (cda_p != cda)
			FreeDoubles(cda);
		if (cda4_p != cda4)
			FreeDoubles(cda4);
		if (cda3_p != cda3)
			FreeDoubles(cda3);
		if (cda2_p != cda2)
			FreeDoubles(cda2);
		if (cda1_p != cda1)
			FreeDoubles(cda1);
		if (bcd_p != bcd)
			FreeDoubles(bcd);
		if (bcd4_p != bcd4)
			FreeDoubles(bcd4);
		if (bcd3_p != bcd3)
			FreeDoubles(bcd3);
		if (bcd2_p != bcd2)
			FreeDoubles(bcd2);
		if (bcd1_p != bcd1)
			FreeDoubles(bcd1);
		if (abc_p != abc)
			FreeDoubles(abc);
		if (abc4_p != abc4)
			FreeDoubles(abc4);
		if (abc3_p != abc3)
			FreeDoubles(abc3);
		if (abc2_p != abc2)
			FreeDoubles(abc2);
		if (abc1_p != abc1)
			FreeDoubles(abc1);
		if (bd_p != bd)
			FreeDoubles(bd);
		if (dexbey_p != dexbey)
			FreeDoubles(dexbey);
		if (bexdey_p != bexdey)
			FreeDoubles(bexdey);
		if (ac_p != ac)
			FreeDoubles(ac);
		if (cexaey_p != cexaey)
			FreeDoubles(cexaey);
		if (aexcey_p != aexcey)
			FreeDoubles(aexcey);
		if (da_p != da)
			FreeDoubles(da);
		if (aexdey_p != aexdey)
			FreeDoubles(aexdey);
		if (dexaey_p != dexaey)
			FreeDoubles(dexaey);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (dexcey_p != dexcey)
			FreeDoubles(dexcey);
		if (cexdey_p != cexdey)
			FreeDoubles(cexdey);
		if (bc_p != bc)
			FreeDoubles(bc);
		if (cexbey_p != cexbey)
			FreeDoubles(cexbey);
		if (bexcey_p != bexcey)
			FreeDoubles(bexcey);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (bexaey_p != bexaey)
			FreeDoubles(bexaey);
		if (aexbey_p != aexbey)
			FreeDoubles(aexbey);
		if (cez_p != cez)
			FreeDoubles(cez);
		if (cey_p != cey)
			FreeDoubles(cey);
		if (cex_p != cex)
			FreeDoubles(cex);
		if (pezd3_p != pezd3)
			FreeDoubles(pezd3);
		if (peyd3_p != peyd3)
			FreeDoubles(peyd3);
		if (pexd3_p != pexd3)
			FreeDoubles(pexd3);
		if (bez_p != bez)
			FreeDoubles(bez);
		if (bey_p != bey)
			FreeDoubles(bey);
		if (bex_p != bex)
			FreeDoubles(bex);
		if (pezd2_p != pezd2)
			FreeDoubles(pezd2);
		if (peyd2_p != peyd2)
			FreeDoubles(peyd2);
		if (pexd2_p != pexd2)
			FreeDoubles(pexd2);
		if (aez_p != aez)
			FreeDoubles(aez);
		if (aey_p != aey)
			FreeDoubles(aey);
		if (aex_p != aex)
			FreeDoubles(aex);
		if (pezd_p != pezd)
			FreeDoubles(pezd);
		if (peyd_p != peyd)
			FreeDoubles(peyd);
		if (pexd_p != pexd)
			FreeDoubles(pexd);
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
		return inSphere_IIIEE_exact<IT, ET>(p1, p2, p3, pdx, pdy, pdz, pex, pey,
		                                    pez);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inSphere_IIIEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3, double pdx, double pdy,
                    double pdz, double pex, double pey, double pez)
{
	Sign ret;
	ret =
	  inSphere_IIIEE_interval<IT, ET>(p1, p2, p3, pdx, pdy, pdz, pex, pey, pez);
	if (is_sign_reliable(ret))
		return ret;
	return inSphere_IIIEE_expansion<IT, ET>(p1, p2, p3, pdx, pdy, pdz, pex, pey,
	                                        pez);
}

template <typename IT, typename ET>
Sign inSphere_IIIEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &pd,
                    const GenericPoint3T<IT, ET> &pe)
{
	return inSphere_IIIEE<IT, ET>(p1, p2, p3, pd.x(), pd.y(), pd.z(), pe.x(),
	                              pe.y(), pe.z());
}

template <typename IT, typename ET>
Sign inSphere_IIIIE_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3,
                             const GenericPoint3T<IT, ET> &p4, IT pex, IT pey,
                             IT pez)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3) ||
	    !p4.getIntervalLambda(l4x, l4y, l4z, d4))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pexd   = pex * d1;
	IT peyd   = pey * d1;
	IT pezd   = pez * d1;
	IT aex    = l1x - pexd;
	IT aey    = l1y - peyd;
	IT aez    = l1z - pezd;
	IT pexd2  = pex * d2;
	IT peyd2  = pey * d2;
	IT pezd2  = pez * d2;
	IT bex    = l2x - pexd2;
	IT bey    = l2y - peyd2;
	IT bez    = l2z - pezd2;
	IT pexd3  = pex * d3;
	IT peyd3  = pey * d3;
	IT pezd3  = pez * d3;
	IT cex    = l3x - pexd3;
	IT cey    = l3y - peyd3;
	IT cez    = l3z - pezd3;
	IT pexd4  = pex * d4;
	IT peyd4  = pey * d4;
	IT pezd4  = pez * d4;
	IT dex    = l4x - pexd4;
	IT dey    = l4y - peyd4;
	IT dez    = l4z - pezd4;
	IT aexbey = aex * bey;
	IT bexaey = bex * aey;
	IT ab     = aexbey - bexaey;
	IT bexcey = bex * cey;
	IT cexbey = cex * bey;
	IT bc     = bexcey - cexbey;
	IT cexdey = cex * dey;
	IT dexcey = dex * cey;
	IT cd     = cexdey - dexcey;
	IT dexaey = dex * aey;
	IT aexdey = aex * dey;
	IT da     = dexaey - aexdey;
	IT aexcey = aex * cey;
	IT cexaey = cex * aey;
	IT ac     = aexcey - cexaey;
	IT bexdey = bex * dey;
	IT dexbey = dex * bey;
	IT bd     = bexdey - dexbey;
	IT abc1   = aez * bc;
	IT abc2   = bez * ac;
	IT abc3   = cez * ab;
	IT abc4   = abc1 + abc3;
	IT abc    = abc4 - abc2;
	IT bcd1   = bez * cd;
	IT bcd2   = cez * bd;
	IT bcd3   = dez * bc;
	IT bcd4   = bcd1 + bcd3;
	IT bcd    = bcd4 - bcd2;
	IT cda1   = cez * da;
	IT cda2   = dez * ac;
	IT cda3   = aez * cd;
	IT cda4   = cda1 + cda3;
	IT cda    = cda4 + cda2;
	IT dab1   = dez * ab;
	IT dab2   = aez * bd;
	IT dab3   = bez * da;
	IT dab4   = dab1 + dab3;
	IT dab    = dab4 + dab2;
	IT al1    = aex * aex;
	IT al2    = aey * aey;
	IT al3    = aez * aez;
	IT al4    = al1 + al2;
	IT alift  = al4 + al3;
	IT bl1    = bex * bex;
	IT bl2    = bey * bey;
	IT bl3    = bez * bez;
	IT bl4    = bl1 + bl2;
	IT blift  = bl4 + bl3;
	IT cl1    = cex * cex;
	IT cl2    = cey * cey;
	IT cl3    = cez * cez;
	IT cl4    = cl1 + cl2;
	IT clift  = cl4 + cl3;
	IT dl1    = dex * dex;
	IT dl2    = dey * dey;
	IT dl3    = dez * dez;
	IT dl4    = dl1 + dl2;
	IT dlift  = dl4 + dl3;
	IT ds1    = dlift * abc;
	IT ds12   = ds1 * d3;
	IT ds2    = clift * dab;
	IT ds22   = ds2 * d4;
	IT dl     = ds22 - ds12;
	IT dlx1   = dl * d1;
	IT dlx2   = dlx1 * d2;
	IT dr1    = blift * cda;
	IT dr12   = dr1 * d1;
	IT dr2    = alift * bcd;
	IT dr22   = dr2 * d2;
	IT dr     = dr22 - dr12;
	IT drx1   = dr * d3;
	IT drx2   = drx1 * d4;
	IT det    = dlx2 + drx2;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIIIE_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2,
                          const GenericPoint3T<IT, ET> &p3,
                          const GenericPoint3T<IT, ET> &p4, ET pex, ET pey,
                          ET pez)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	p4.getExactLambda(l4x, l4y, l4z, d4);
	ET pexd   = pex * d1;
	ET peyd   = pey * d1;
	ET pezd   = pez * d1;
	ET aex    = l1x - pexd;
	ET aey    = l1y - peyd;
	ET aez    = l1z - pezd;
	ET pexd2  = pex * d2;
	ET peyd2  = pey * d2;
	ET pezd2  = pez * d2;
	ET bex    = l2x - pexd2;
	ET bey    = l2y - peyd2;
	ET bez    = l2z - pezd2;
	ET pexd3  = pex * d3;
	ET peyd3  = pey * d3;
	ET pezd3  = pez * d3;
	ET cex    = l3x - pexd3;
	ET cey    = l3y - peyd3;
	ET cez    = l3z - pezd3;
	ET pexd4  = pex * d4;
	ET peyd4  = pey * d4;
	ET pezd4  = pez * d4;
	ET dex    = l4x - pexd4;
	ET dey    = l4y - peyd4;
	ET dez    = l4z - pezd4;
	ET aexbey = aex * bey;
	ET bexaey = bex * aey;
	ET ab     = aexbey - bexaey;
	ET bexcey = bex * cey;
	ET cexbey = cex * bey;
	ET bc     = bexcey - cexbey;
	ET cexdey = cex * dey;
	ET dexcey = dex * cey;
	ET cd     = cexdey - dexcey;
	ET dexaey = dex * aey;
	ET aexdey = aex * dey;
	ET da     = dexaey - aexdey;
	ET aexcey = aex * cey;
	ET cexaey = cex * aey;
	ET ac     = aexcey - cexaey;
	ET bexdey = bex * dey;
	ET dexbey = dex * bey;
	ET bd     = bexdey - dexbey;
	ET abc1   = aez * bc;
	ET abc2   = bez * ac;
	ET abc3   = cez * ab;
	ET abc4   = abc1 + abc3;
	ET abc    = abc4 - abc2;
	ET bcd1   = bez * cd;
	ET bcd2   = cez * bd;
	ET bcd3   = dez * bc;
	ET bcd4   = bcd1 + bcd3;
	ET bcd    = bcd4 - bcd2;
	ET cda1   = cez * da;
	ET cda2   = dez * ac;
	ET cda3   = aez * cd;
	ET cda4   = cda1 + cda3;
	ET cda    = cda4 + cda2;
	ET dab1   = dez * ab;
	ET dab2   = aez * bd;
	ET dab3   = bez * da;
	ET dab4   = dab1 + dab3;
	ET dab    = dab4 + dab2;
	ET al1    = aex * aex;
	ET al2    = aey * aey;
	ET al3    = aez * aez;
	ET al4    = al1 + al2;
	ET alift  = al4 + al3;
	ET bl1    = bex * bex;
	ET bl2    = bey * bey;
	ET bl3    = bez * bez;
	ET bl4    = bl1 + bl2;
	ET blift  = bl4 + bl3;
	ET cl1    = cex * cex;
	ET cl2    = cey * cey;
	ET cl3    = cez * cez;
	ET cl4    = cl1 + cl2;
	ET clift  = cl4 + cl3;
	ET dl1    = dex * dex;
	ET dl2    = dey * dey;
	ET dl3    = dez * dez;
	ET dl4    = dl1 + dl2;
	ET dlift  = dl4 + dl3;
	ET ds1    = dlift * abc;
	ET ds12   = ds1 * d3;
	ET ds2    = clift * dab;
	ET ds22   = ds2 * d4;
	ET dl     = ds22 - ds12;
	ET dlx1   = dl * d1;
	ET dlx2   = dlx1 * d2;
	ET dr1    = blift * cda;
	ET dr12   = dr1 * d1;
	ET dr2    = alift * bcd;
	ET dr22   = dr2 * d2;
	ET dr     = dr22 - dr12;
	ET drx1   = dr * d3;
	ET drx2   = drx1 * d4;
	ET det    = dlx2 + drx2;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIIIE_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3,
                              const GenericPoint3T<IT, ET> &p4, double pex,
                              double pey, double pez)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[16],
	  *l1x = l1x_p, l1y_p[16], *l1y = l1y_p, l1z_p[16], *l1z = l1z_p, d1_p[16],
	  *d1 = d1_p, l2x_p[16], *l2x = l2x_p, l2y_p[16], *l2y = l2y_p, l2z_p[16],
	  *l2z = l2z_p, d2_p[16], *d2 = d2_p, l3x_p[16], *l3x = l3x_p, l3y_p[16],
	  *l3y = l3y_p, l3z_p[16], *l3z = l3z_p, d3_p[16], *d3 = d3_p, l4x_p[16],
	  *l4x = l4x_p, l4y_p[16], *l4y = l4y_p, l4z_p[16], *l4z = l4z_p, d4_p[16],
	  *d4       = d4_p;
	int l1x_len = 16, l1y_len = 16, l1z_len = 16, d1_len = 16, l2x_len = 16,
	    l2y_len = 16, l2z_len = 16, d2_len = 16, l3x_len = 16, l3y_len = 16,
	    l3z_len = 16, d3_len = 16, l4x_len = 16, l4y_len = 16, l4z_len = 16,
	    d4_len = 16;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	p4.getExpansionLambda(&l4x, l4x_len, &l4y, l4y_len, &l4z, l4z_len, &d4,
	                      d4_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) &&
	    (d4[d4_len - 1] != 0))
	{
		expansionObject o;
		double          pexd_p[16], *pexd = pexd_p;
		int    pexd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pex, &pexd, 16);
		double peyd_p[16], *peyd = peyd_p;
		int    peyd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pey, &peyd, 16);
		double pezd_p[16], *pezd = pezd_p;
		int    pezd_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, pez, &pezd, 16);
		double aex_p[16], *aex = aex_p;
		int    aex_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, pexd_len, pexd, &aex, 16);
		double aey_p[16], *aey = aey_p;
		int    aey_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, peyd_len, peyd, &aey, 16);
		double aez_p[16], *aez = aez_p;
		int    aez_len =
		  o.Gen_Diff_With_PreAlloc(l1z_len, l1z, pezd_len, pezd, &aez, 16);
		double pexd2_p[16], *pexd2 = pexd2_p;
		int    pexd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pex, &pexd2, 16);
		double peyd2_p[16], *peyd2 = peyd2_p;
		int    peyd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pey, &peyd2, 16);
		double pezd2_p[16], *pezd2 = pezd2_p;
		int    pezd2_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, pez, &pezd2, 16);
		double bex_p[16], *bex = bex_p;
		int    bex_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, pexd2_len, pexd2, &bex, 16);
		double bey_p[16], *bey = bey_p;
		int    bey_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, peyd2_len, peyd2, &bey, 16);
		double bez_p[16], *bez = bez_p;
		int    bez_len =
		  o.Gen_Diff_With_PreAlloc(l2z_len, l2z, pezd2_len, pezd2, &bez, 16);
		double pexd3_p[16], *pexd3 = pexd3_p;
		int    pexd3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pex, &pexd3, 16);
		double peyd3_p[16], *peyd3 = peyd3_p;
		int    peyd3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pey, &peyd3, 16);
		double pezd3_p[16], *pezd3 = pezd3_p;
		int    pezd3_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, pez, &pezd3, 16);
		double cex_p[16], *cex = cex_p;
		int    cex_len =
		  o.Gen_Diff_With_PreAlloc(l3x_len, l3x, pexd3_len, pexd3, &cex, 16);
		double cey_p[16], *cey = cey_p;
		int    cey_len =
		  o.Gen_Diff_With_PreAlloc(l3y_len, l3y, peyd3_len, peyd3, &cey, 16);
		double cez_p[16], *cez = cez_p;
		int    cez_len =
		  o.Gen_Diff_With_PreAlloc(l3z_len, l3z, pezd3_len, pezd3, &cez, 16);
		double pexd4_p[16], *pexd4 = pexd4_p;
		int    pexd4_len = o.Gen_Scale_With_PreAlloc(d4_len, d4, pex, &pexd4, 16);
		double peyd4_p[16], *peyd4 = peyd4_p;
		int    peyd4_len = o.Gen_Scale_With_PreAlloc(d4_len, d4, pey, &peyd4, 16);
		double pezd4_p[16], *pezd4 = pezd4_p;
		int    pezd4_len = o.Gen_Scale_With_PreAlloc(d4_len, d4, pez, &pezd4, 16);
		double dex_p[16], *dex = dex_p;
		int    dex_len =
		  o.Gen_Diff_With_PreAlloc(l4x_len, l4x, pexd4_len, pexd4, &dex, 16);
		double dey_p[16], *dey = dey_p;
		int    dey_len =
		  o.Gen_Diff_With_PreAlloc(l4y_len, l4y, peyd4_len, peyd4, &dey, 16);
		double dez_p[16], *dez = dez_p;
		int    dez_len =
		  o.Gen_Diff_With_PreAlloc(l4z_len, l4z, pezd4_len, pezd4, &dez, 16);
		double aexbey_p[16], *aexbey = aexbey_p;
		int    aexbey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, bey_len, bey, &aexbey, 16);
		double bexaey_p[16], *bexaey = bexaey_p;
		int    bexaey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, aey_len, aey, &bexaey, 16);
		double ab_p[16], *ab = ab_p;
		int    ab_len =
		  o.Gen_Diff_With_PreAlloc(aexbey_len, aexbey, bexaey_len, bexaey, &ab, 16);
		double bexcey_p[16], *bexcey = bexcey_p;
		int    bexcey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, cey_len, cey, &bexcey, 16);
		double cexbey_p[16], *cexbey = cexbey_p;
		int    cexbey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, bey_len, bey, &cexbey, 16);
		double bc_p[16], *bc = bc_p;
		int    bc_len =
		  o.Gen_Diff_With_PreAlloc(bexcey_len, bexcey, cexbey_len, cexbey, &bc, 16);
		double cexdey_p[16], *cexdey = cexdey_p;
		int    cexdey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, dey_len, dey, &cexdey, 16);
		double dexcey_p[16], *dexcey = dexcey_p;
		int    dexcey_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, cey_len, cey, &dexcey, 16);
		double cd_p[16], *cd = cd_p;
		int    cd_len =
		  o.Gen_Diff_With_PreAlloc(cexdey_len, cexdey, dexcey_len, dexcey, &cd, 16);
		double dexaey_p[16], *dexaey = dexaey_p;
		int    dexaey_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, aey_len, aey, &dexaey, 16);
		double aexdey_p[16], *aexdey = aexdey_p;
		int    aexdey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, dey_len, dey, &aexdey, 16);
		double da_p[16], *da = da_p;
		int    da_len =
		  o.Gen_Diff_With_PreAlloc(dexaey_len, dexaey, aexdey_len, aexdey, &da, 16);
		double aexcey_p[16], *aexcey = aexcey_p;
		int    aexcey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, cey_len, cey, &aexcey, 16);
		double cexaey_p[16], *cexaey = cexaey_p;
		int    cexaey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, aey_len, aey, &cexaey, 16);
		double ac_p[16], *ac = ac_p;
		int    ac_len =
		  o.Gen_Diff_With_PreAlloc(aexcey_len, aexcey, cexaey_len, cexaey, &ac, 16);
		double bexdey_p[16], *bexdey = bexdey_p;
		int    bexdey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, dey_len, dey, &bexdey, 16);
		double dexbey_p[16], *dexbey = dexbey_p;
		int    dexbey_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, bey_len, bey, &dexbey, 16);
		double bd_p[16], *bd = bd_p;
		int    bd_len =
		  o.Gen_Diff_With_PreAlloc(bexdey_len, bexdey, dexbey_len, dexbey, &bd, 16);
		double abc1_p[16], *abc1 = abc1_p;
		int    abc1_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bc_len, bc, &abc1, 16);
		double abc2_p[16], *abc2 = abc2_p;
		int    abc2_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, ac_len, ac, &abc2, 16);
		double abc3_p[16], *abc3 = abc3_p;
		int    abc3_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, ab_len, ab, &abc3, 16);
		double abc4_p[16], *abc4 = abc4_p;
		int    abc4_len =
		  o.Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 16);
		double abc_p[16], *abc = abc_p;
		int    abc_len =
		  o.Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 16);
		double bcd1_p[16], *bcd1 = bcd1_p;
		int    bcd1_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, cd_len, cd, &bcd1, 16);
		double bcd2_p[16], *bcd2 = bcd2_p;
		int    bcd2_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, bd_len, bd, &bcd2, 16);
		double bcd3_p[16], *bcd3 = bcd3_p;
		int    bcd3_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, bc_len, bc, &bcd3, 16);
		double bcd4_p[16], *bcd4 = bcd4_p;
		int    bcd4_len =
		  o.Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 16);
		double bcd_p[16], *bcd = bcd_p;
		int    bcd_len =
		  o.Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 16);
		double cda1_p[16], *cda1 = cda1_p;
		int    cda1_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, da_len, da, &cda1, 16);
		double cda2_p[16], *cda2 = cda2_p;
		int    cda2_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, ac_len, ac, &cda2, 16);
		double cda3_p[16], *cda3 = cda3_p;
		int    cda3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, cd_len, cd, &cda3, 16);
		double cda4_p[16], *cda4 = cda4_p;
		int    cda4_len =
		  o.Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 16);
		double cda_p[16], *cda = cda_p;
		int    cda_len =
		  o.Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 16);
		double dab1_p[16], *dab1 = dab1_p;
		int    dab1_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, ab_len, ab, &dab1, 16);
		double dab2_p[16], *dab2 = dab2_p;
		int    dab2_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bd_len, bd, &dab2, 16);
		double dab3_p[16], *dab3 = dab3_p;
		int    dab3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, da_len, da, &dab3, 16);
		double dab4_p[16], *dab4 = dab4_p;
		int    dab4_len =
		  o.Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 16);
		double dab_p[16], *dab = dab_p;
		int    dab_len =
		  o.Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 16);
		double al1_p[16], *al1 = al1_p;
		int    al1_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, aex_len, aex, &al1, 16);
		double al2_p[16], *al2 = al2_p;
		int    al2_len =
		  o.Gen_Product_With_PreAlloc(aey_len, aey, aey_len, aey, &al2, 16);
		double al3_p[16], *al3 = al3_p;
		int    al3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, aez_len, aez, &al3, 16);
		double al4_p[16], *al4 = al4_p;
		int al4_len = o.Gen_Sum_With_PreAlloc(al1_len, al1, al2_len, al2, &al4, 16);
		double alift_p[16], *alift = alift_p;
		int    alift_len =
		  o.Gen_Sum_With_PreAlloc(al4_len, al4, al3_len, al3, &alift, 16);
		double bl1_p[16], *bl1 = bl1_p;
		int    bl1_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, bex_len, bex, &bl1, 16);
		double bl2_p[16], *bl2 = bl2_p;
		int    bl2_len =
		  o.Gen_Product_With_PreAlloc(bey_len, bey, bey_len, bey, &bl2, 16);
		double bl3_p[16], *bl3 = bl3_p;
		int    bl3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, bez_len, bez, &bl3, 16);
		double bl4_p[16], *bl4 = bl4_p;
		int bl4_len = o.Gen_Sum_With_PreAlloc(bl1_len, bl1, bl2_len, bl2, &bl4, 16);
		double blift_p[16], *blift = blift_p;
		int    blift_len =
		  o.Gen_Sum_With_PreAlloc(bl4_len, bl4, bl3_len, bl3, &blift, 16);
		double cl1_p[16], *cl1 = cl1_p;
		int    cl1_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, cex_len, cex, &cl1, 16);
		double cl2_p[16], *cl2 = cl2_p;
		int    cl2_len =
		  o.Gen_Product_With_PreAlloc(cey_len, cey, cey_len, cey, &cl2, 16);
		double cl3_p[16], *cl3 = cl3_p;
		int    cl3_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, cez_len, cez, &cl3, 16);
		double cl4_p[16], *cl4 = cl4_p;
		int cl4_len = o.Gen_Sum_With_PreAlloc(cl1_len, cl1, cl2_len, cl2, &cl4, 16);
		double clift_p[16], *clift = clift_p;
		int    clift_len =
		  o.Gen_Sum_With_PreAlloc(cl4_len, cl4, cl3_len, cl3, &clift, 16);
		double dl1_p[16], *dl1 = dl1_p;
		int    dl1_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, dex_len, dex, &dl1, 16);
		double dl2_p[16], *dl2 = dl2_p;
		int    dl2_len =
		  o.Gen_Product_With_PreAlloc(dey_len, dey, dey_len, dey, &dl2, 16);
		double dl3_p[16], *dl3 = dl3_p;
		int    dl3_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, dez_len, dez, &dl3, 16);
		double dl4_p[16], *dl4 = dl4_p;
		int dl4_len = o.Gen_Sum_With_PreAlloc(dl1_len, dl1, dl2_len, dl2, &dl4, 16);
		double dlift_p[16], *dlift = dlift_p;
		int    dlift_len =
		  o.Gen_Sum_With_PreAlloc(dl4_len, dl4, dl3_len, dl3, &dlift, 16);
		double ds1_p[16], *ds1 = ds1_p;
		int    ds1_len =
		  o.Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 16);
		double ds12_p[16], *ds12 = ds12_p;
		int    ds12_len =
		  o.Gen_Product_With_PreAlloc(ds1_len, ds1, d3_len, d3, &ds12, 16);
		double ds2_p[16], *ds2 = ds2_p;
		int    ds2_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 16);
		double ds22_p[16], *ds22 = ds22_p;
		int    ds22_len =
		  o.Gen_Product_With_PreAlloc(ds2_len, ds2, d4_len, d4, &ds22, 16);
		double dl_p[16], *dl = dl_p;
		int    dl_len =
		  o.Gen_Diff_With_PreAlloc(ds22_len, ds22, ds12_len, ds12, &dl, 16);
		double dlx1_p[16], *dlx1 = dlx1_p;
		int    dlx1_len =
		  o.Gen_Product_With_PreAlloc(dl_len, dl, d1_len, d1, &dlx1, 16);
		double dlx2_p[16], *dlx2 = dlx2_p;
		int    dlx2_len =
		  o.Gen_Product_With_PreAlloc(dlx1_len, dlx1, d2_len, d2, &dlx2, 16);
		double dr1_p[16], *dr1 = dr1_p;
		int    dr1_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1, 16);
		double dr12_p[16], *dr12 = dr12_p;
		int    dr12_len =
		  o.Gen_Product_With_PreAlloc(dr1_len, dr1, d1_len, d1, &dr12, 16);
		double dr2_p[16], *dr2 = dr2_p;
		int    dr2_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 16);
		double dr22_p[16], *dr22 = dr22_p;
		int    dr22_len =
		  o.Gen_Product_With_PreAlloc(dr2_len, dr2, d2_len, d2, &dr22, 16);
		double dr_p[16], *dr = dr_p;
		int    dr_len =
		  o.Gen_Diff_With_PreAlloc(dr22_len, dr22, dr12_len, dr12, &dr, 16);
		double drx1_p[16], *drx1 = drx1_p;
		int    drx1_len =
		  o.Gen_Product_With_PreAlloc(dr_len, dr, d3_len, d3, &drx1, 16);
		double drx2_p[16], *drx2 = drx2_p;
		int    drx2_len =
		  o.Gen_Product_With_PreAlloc(drx1_len, drx1, d4_len, d4, &drx2, 16);
		double det_p[16], *det = det_p;
		int    det_len =
		  o.Gen_Sum_With_PreAlloc(dlx2_len, dlx2, drx2_len, drx2, &det, 16);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (drx2_p != drx2)
			FreeDoubles(drx2);
		if (drx1_p != drx1)
			FreeDoubles(drx1);
		if (dr_p != dr)
			FreeDoubles(dr);
		if (dr22_p != dr22)
			FreeDoubles(dr22);
		if (dr2_p != dr2)
			FreeDoubles(dr2);
		if (dr12_p != dr12)
			FreeDoubles(dr12);
		if (dr1_p != dr1)
			FreeDoubles(dr1);
		if (dlx2_p != dlx2)
			FreeDoubles(dlx2);
		if (dlx1_p != dlx1)
			FreeDoubles(dlx1);
		if (dl_p != dl)
			FreeDoubles(dl);
		if (ds22_p != ds22)
			FreeDoubles(ds22);
		if (ds2_p != ds2)
			FreeDoubles(ds2);
		if (ds12_p != ds12)
			FreeDoubles(ds12);
		if (ds1_p != ds1)
			FreeDoubles(ds1);
		if (dlift_p != dlift)
			FreeDoubles(dlift);
		if (dl4_p != dl4)
			FreeDoubles(dl4);
		if (dl3_p != dl3)
			FreeDoubles(dl3);
		if (dl2_p != dl2)
			FreeDoubles(dl2);
		if (dl1_p != dl1)
			FreeDoubles(dl1);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (cl4_p != cl4)
			FreeDoubles(cl4);
		if (cl3_p != cl3)
			FreeDoubles(cl3);
		if (cl2_p != cl2)
			FreeDoubles(cl2);
		if (cl1_p != cl1)
			FreeDoubles(cl1);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bl4_p != bl4)
			FreeDoubles(bl4);
		if (bl3_p != bl3)
			FreeDoubles(bl3);
		if (bl2_p != bl2)
			FreeDoubles(bl2);
		if (bl1_p != bl1)
			FreeDoubles(bl1);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (al4_p != al4)
			FreeDoubles(al4);
		if (al3_p != al3)
			FreeDoubles(al3);
		if (al2_p != al2)
			FreeDoubles(al2);
		if (al1_p != al1)
			FreeDoubles(al1);
		if (dab_p != dab)
			FreeDoubles(dab);
		if (dab4_p != dab4)
			FreeDoubles(dab4);
		if (dab3_p != dab3)
			FreeDoubles(dab3);
		if (dab2_p != dab2)
			FreeDoubles(dab2);
		if (dab1_p != dab1)
			FreeDoubles(dab1);
		if (cda_p != cda)
			FreeDoubles(cda);
		if (cda4_p != cda4)
			FreeDoubles(cda4);
		if (cda3_p != cda3)
			FreeDoubles(cda3);
		if (cda2_p != cda2)
			FreeDoubles(cda2);
		if (cda1_p != cda1)
			FreeDoubles(cda1);
		if (bcd_p != bcd)
			FreeDoubles(bcd);
		if (bcd4_p != bcd4)
			FreeDoubles(bcd4);
		if (bcd3_p != bcd3)
			FreeDoubles(bcd3);
		if (bcd2_p != bcd2)
			FreeDoubles(bcd2);
		if (bcd1_p != bcd1)
			FreeDoubles(bcd1);
		if (abc_p != abc)
			FreeDoubles(abc);
		if (abc4_p != abc4)
			FreeDoubles(abc4);
		if (abc3_p != abc3)
			FreeDoubles(abc3);
		if (abc2_p != abc2)
			FreeDoubles(abc2);
		if (abc1_p != abc1)
			FreeDoubles(abc1);
		if (bd_p != bd)
			FreeDoubles(bd);
		if (dexbey_p != dexbey)
			FreeDoubles(dexbey);
		if (bexdey_p != bexdey)
			FreeDoubles(bexdey);
		if (ac_p != ac)
			FreeDoubles(ac);
		if (cexaey_p != cexaey)
			FreeDoubles(cexaey);
		if (aexcey_p != aexcey)
			FreeDoubles(aexcey);
		if (da_p != da)
			FreeDoubles(da);
		if (aexdey_p != aexdey)
			FreeDoubles(aexdey);
		if (dexaey_p != dexaey)
			FreeDoubles(dexaey);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (dexcey_p != dexcey)
			FreeDoubles(dexcey);
		if (cexdey_p != cexdey)
			FreeDoubles(cexdey);
		if (bc_p != bc)
			FreeDoubles(bc);
		if (cexbey_p != cexbey)
			FreeDoubles(cexbey);
		if (bexcey_p != bexcey)
			FreeDoubles(bexcey);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (bexaey_p != bexaey)
			FreeDoubles(bexaey);
		if (aexbey_p != aexbey)
			FreeDoubles(aexbey);
		if (dez_p != dez)
			FreeDoubles(dez);
		if (dey_p != dey)
			FreeDoubles(dey);
		if (dex_p != dex)
			FreeDoubles(dex);
		if (pezd4_p != pezd4)
			FreeDoubles(pezd4);
		if (peyd4_p != peyd4)
			FreeDoubles(peyd4);
		if (pexd4_p != pexd4)
			FreeDoubles(pexd4);
		if (cez_p != cez)
			FreeDoubles(cez);
		if (cey_p != cey)
			FreeDoubles(cey);
		if (cex_p != cex)
			FreeDoubles(cex);
		if (pezd3_p != pezd3)
			FreeDoubles(pezd3);
		if (peyd3_p != peyd3)
			FreeDoubles(peyd3);
		if (pexd3_p != pexd3)
			FreeDoubles(pexd3);
		if (bez_p != bez)
			FreeDoubles(bez);
		if (bey_p != bey)
			FreeDoubles(bey);
		if (bex_p != bex)
			FreeDoubles(bex);
		if (pezd2_p != pezd2)
			FreeDoubles(pezd2);
		if (peyd2_p != peyd2)
			FreeDoubles(peyd2);
		if (pexd2_p != pexd2)
			FreeDoubles(pexd2);
		if (aez_p != aez)
			FreeDoubles(aez);
		if (aey_p != aey)
			FreeDoubles(aey);
		if (aex_p != aex)
			FreeDoubles(aex);
		if (pezd_p != pezd)
			FreeDoubles(pezd);
		if (peyd_p != peyd)
			FreeDoubles(peyd);
		if (pexd_p != pexd)
			FreeDoubles(pexd);
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

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l4x_p != l4x)
			FreeDoubles(l4x);
		if (l4y_p != l4y)
			FreeDoubles(l4y);
		if (l4z_p != l4z)
			FreeDoubles(l4z);
		if (d4_p != d4)
			FreeDoubles(d4);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return inSphere_IIIIE_exact<IT, ET>(p1, p2, p3, p4, pex, pey, pez);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inSphere_IIIIE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &p4, double pex, double pey,
                    double pez)
{
	Sign ret;
	ret = inSphere_IIIIE_interval<IT, ET>(p1, p2, p3, p4, pex, pey, pez);
	if (is_sign_reliable(ret))
		return ret;
	return inSphere_IIIIE_expansion<IT, ET>(p1, p2, p3, p4, pex, pey, pez);
}

template <typename IT, typename ET>
Sign inSphere_IIIIE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &p4,
                    const GenericPoint3T<IT, ET> &pe)
{
	return inSphere_IIIIE<IT, ET>(p1, p2, p3, p4, pe.x(), pe.y(), pe.z());
}

template <typename IT, typename ET>
Sign inSphere_IIIII_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3,
                             const GenericPoint3T<IT, ET> &p4,
                             const GenericPoint3T<IT, ET> &p5)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4,
	  l5x, l5y, l5z, d5;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3) ||
	    !p4.getIntervalLambda(l4x, l4y, l4z, d4) ||
	    !p5.getIntervalLambda(l5x, l5y, l5z, d5))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT pexd   = l5x * d1;
	IT peyd   = l5y * d1;
	IT pezd   = l5z * d1;
	IT ll1x   = l1x * d5;
	IT ll1y   = l1y * d5;
	IT ll1z   = l1z * d5;
	IT aex    = ll1x - pexd;
	IT aey    = ll1y - peyd;
	IT aez    = ll1z - pezd;
	IT pexd2  = l5x * d2;
	IT peyd2  = l5y * d2;
	IT pezd2  = l5z * d2;
	IT ll2x   = l2x * d5;
	IT ll2y   = l2y * d5;
	IT ll2z   = l2z * d5;
	IT bex    = ll2x - pexd2;
	IT bey    = ll2y - peyd2;
	IT bez    = ll2z - pezd2;
	IT pexd3  = l5x * d3;
	IT peyd3  = l5y * d3;
	IT pezd3  = l5z * d3;
	IT ll3x   = l3x * d5;
	IT ll3y   = l3y * d5;
	IT ll3z   = l3z * d5;
	IT cex    = ll3x - pexd3;
	IT cey    = ll3y - peyd3;
	IT cez    = ll3z - pezd3;
	IT pexd4  = l5x * d4;
	IT peyd4  = l5y * d4;
	IT pezd4  = l5z * d4;
	IT ll4x   = l4x * d5;
	IT ll4y   = l4y * d5;
	IT ll4z   = l4z * d5;
	IT dex    = ll4x - pexd4;
	IT dey    = ll4y - peyd4;
	IT dez    = ll4z - pezd4;
	IT aexbey = aex * bey;
	IT bexaey = bex * aey;
	IT ab     = aexbey - bexaey;
	IT bexcey = bex * cey;
	IT cexbey = cex * bey;
	IT bc     = bexcey - cexbey;
	IT cexdey = cex * dey;
	IT dexcey = dex * cey;
	IT cd     = cexdey - dexcey;
	IT dexaey = dex * aey;
	IT aexdey = aex * dey;
	IT da     = dexaey - aexdey;
	IT aexcey = aex * cey;
	IT cexaey = cex * aey;
	IT ac     = aexcey - cexaey;
	IT bexdey = bex * dey;
	IT dexbey = dex * bey;
	IT bd     = bexdey - dexbey;
	IT abc1   = aez * bc;
	IT abc2   = bez * ac;
	IT abc3   = cez * ab;
	IT abc4   = abc1 + abc3;
	IT abc    = abc4 - abc2;
	IT bcd1   = bez * cd;
	IT bcd2   = cez * bd;
	IT bcd3   = dez * bc;
	IT bcd4   = bcd1 + bcd3;
	IT bcd    = bcd4 - bcd2;
	IT cda1   = cez * da;
	IT cda2   = dez * ac;
	IT cda3   = aez * cd;
	IT cda4   = cda1 + cda3;
	IT cda    = cda4 + cda2;
	IT dab1   = dez * ab;
	IT dab2   = aez * bd;
	IT dab3   = bez * da;
	IT dab4   = dab1 + dab3;
	IT dab    = dab4 + dab2;
	IT al1    = aex * aex;
	IT al2    = aey * aey;
	IT al3    = aez * aez;
	IT al4    = al1 + al2;
	IT alift  = al4 + al3;
	IT bl1    = bex * bex;
	IT bl2    = bey * bey;
	IT bl3    = bez * bez;
	IT bl4    = bl1 + bl2;
	IT blift  = bl4 + bl3;
	IT cl1    = cex * cex;
	IT cl2    = cey * cey;
	IT cl3    = cez * cez;
	IT cl4    = cl1 + cl2;
	IT clift  = cl4 + cl3;
	IT dl1    = dex * dex;
	IT dl2    = dey * dey;
	IT dl3    = dez * dez;
	IT dl4    = dl1 + dl2;
	IT dlift  = dl4 + dl3;
	IT ds1    = dlift * abc;
	IT ds1n   = ds1 * d3;
	IT ds2    = clift * dab;
	IT ds2n   = ds2 * d4;
	IT dl     = ds2n - ds1n;
	IT dla    = dl * d1;
	IT dlb    = dla * d2;
	IT dr1    = blift * cda;
	IT dr1n   = dr1 * d1;
	IT dr2    = alift * bcd;
	IT dr2n   = dr2 * d2;
	IT dr     = dr2n - dr1n;
	IT dra    = dr * d3;
	IT drb    = dra * d4;
	IT det    = dlb + drb;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIIII_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2,
                          const GenericPoint3T<IT, ET> &p3,
                          const GenericPoint3T<IT, ET> &p4,
                          const GenericPoint3T<IT, ET> &p5)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4,
	  l5x, l5y, l5z, d5;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	p4.getExactLambda(l4x, l4y, l4z, d4);
	p5.getExactLambda(l5x, l5y, l5z, d5);
	ET pexd   = l5x * d1;
	ET peyd   = l5y * d1;
	ET pezd   = l5z * d1;
	ET ll1x   = l1x * d5;
	ET ll1y   = l1y * d5;
	ET ll1z   = l1z * d5;
	ET aex    = ll1x - pexd;
	ET aey    = ll1y - peyd;
	ET aez    = ll1z - pezd;
	ET pexd2  = l5x * d2;
	ET peyd2  = l5y * d2;
	ET pezd2  = l5z * d2;
	ET ll2x   = l2x * d5;
	ET ll2y   = l2y * d5;
	ET ll2z   = l2z * d5;
	ET bex    = ll2x - pexd2;
	ET bey    = ll2y - peyd2;
	ET bez    = ll2z - pezd2;
	ET pexd3  = l5x * d3;
	ET peyd3  = l5y * d3;
	ET pezd3  = l5z * d3;
	ET ll3x   = l3x * d5;
	ET ll3y   = l3y * d5;
	ET ll3z   = l3z * d5;
	ET cex    = ll3x - pexd3;
	ET cey    = ll3y - peyd3;
	ET cez    = ll3z - pezd3;
	ET pexd4  = l5x * d4;
	ET peyd4  = l5y * d4;
	ET pezd4  = l5z * d4;
	ET ll4x   = l4x * d5;
	ET ll4y   = l4y * d5;
	ET ll4z   = l4z * d5;
	ET dex    = ll4x - pexd4;
	ET dey    = ll4y - peyd4;
	ET dez    = ll4z - pezd4;
	ET aexbey = aex * bey;
	ET bexaey = bex * aey;
	ET ab     = aexbey - bexaey;
	ET bexcey = bex * cey;
	ET cexbey = cex * bey;
	ET bc     = bexcey - cexbey;
	ET cexdey = cex * dey;
	ET dexcey = dex * cey;
	ET cd     = cexdey - dexcey;
	ET dexaey = dex * aey;
	ET aexdey = aex * dey;
	ET da     = dexaey - aexdey;
	ET aexcey = aex * cey;
	ET cexaey = cex * aey;
	ET ac     = aexcey - cexaey;
	ET bexdey = bex * dey;
	ET dexbey = dex * bey;
	ET bd     = bexdey - dexbey;
	ET abc1   = aez * bc;
	ET abc2   = bez * ac;
	ET abc3   = cez * ab;
	ET abc4   = abc1 + abc3;
	ET abc    = abc4 - abc2;
	ET bcd1   = bez * cd;
	ET bcd2   = cez * bd;
	ET bcd3   = dez * bc;
	ET bcd4   = bcd1 + bcd3;
	ET bcd    = bcd4 - bcd2;
	ET cda1   = cez * da;
	ET cda2   = dez * ac;
	ET cda3   = aez * cd;
	ET cda4   = cda1 + cda3;
	ET cda    = cda4 + cda2;
	ET dab1   = dez * ab;
	ET dab2   = aez * bd;
	ET dab3   = bez * da;
	ET dab4   = dab1 + dab3;
	ET dab    = dab4 + dab2;
	ET al1    = aex * aex;
	ET al2    = aey * aey;
	ET al3    = aez * aez;
	ET al4    = al1 + al2;
	ET alift  = al4 + al3;
	ET bl1    = bex * bex;
	ET bl2    = bey * bey;
	ET bl3    = bez * bez;
	ET bl4    = bl1 + bl2;
	ET blift  = bl4 + bl3;
	ET cl1    = cex * cex;
	ET cl2    = cey * cey;
	ET cl3    = cez * cez;
	ET cl4    = cl1 + cl2;
	ET clift  = cl4 + cl3;
	ET dl1    = dex * dex;
	ET dl2    = dey * dey;
	ET dl3    = dez * dez;
	ET dl4    = dl1 + dl2;
	ET dlift  = dl4 + dl3;
	ET ds1    = dlift * abc;
	ET ds1n   = ds1 * d3;
	ET ds2    = clift * dab;
	ET ds2n   = ds2 * d4;
	ET dl     = ds2n - ds1n;
	ET dla    = dl * d1;
	ET dlb    = dla * d2;
	ET dr1    = blift * cda;
	ET dr1n   = dr1 * d1;
	ET dr2    = alift * bcd;
	ET dr2n   = dr2 * d2;
	ET dr     = dr2n - dr1n;
	ET dra    = dr * d3;
	ET drb    = dra * d4;
	ET det    = dlb + drb;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign inSphere_IIIII_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3,
                              const GenericPoint3T<IT, ET> &p4,
                              const GenericPoint3T<IT, ET> &p5)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[8],
	  *l1x = l1x_p, l1y_p[8], *l1y = l1y_p, l1z_p[8], *l1z = l1z_p, d1_p[8],
	  *d1 = d1_p, l2x_p[8], *l2x = l2x_p, l2y_p[8], *l2y = l2y_p, l2z_p[8],
	  *l2z = l2z_p, d2_p[8], *d2 = d2_p, l3x_p[8], *l3x = l3x_p, l3y_p[8],
	  *l3y = l3y_p, l3z_p[8], *l3z = l3z_p, d3_p[8], *d3 = d3_p, l4x_p[8],
	  *l4x = l4x_p, l4y_p[8], *l4y = l4y_p, l4z_p[8], *l4z = l4z_p, d4_p[8],
	  *d4 = d4_p, l5x_p[8], *l5x = l5x_p, l5y_p[8], *l5y = l5y_p, l5z_p[8],
	  *l5z = l5z_p, d5_p[8], *d5 = d5_p;
	int l1x_len = 8, l1y_len = 8, l1z_len = 8, d1_len = 8, l2x_len = 8,
	    l2y_len = 8, l2z_len = 8, d2_len = 8, l3x_len = 8, l3y_len = 8,
	    l3z_len = 8, d3_len = 8, l4x_len = 8, l4y_len = 8, l4z_len = 8,
	    d4_len = 8, l5x_len = 8, l5y_len = 8, l5z_len = 8, d5_len = 8;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	p4.getExpansionLambda(&l4x, l4x_len, &l4y, l4y_len, &l4z, l4z_len, &d4,
	                      d4_len);
	p5.getExpansionLambda(&l5x, l5x_len, &l5y, l5y_len, &l5z, l5z_len, &d5,
	                      d5_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) &&
	    (d4[d4_len - 1] != 0) && (d5[d5_len - 1] != 0))
	{
		expansionObject o;
		double          pexd_p[8], *pexd = pexd_p;
		int             pexd_len =
		  o.Gen_Product_With_PreAlloc(l5x_len, l5x, d1_len, d1, &pexd, 8);
		double peyd_p[8], *peyd = peyd_p;
		int    peyd_len =
		  o.Gen_Product_With_PreAlloc(l5y_len, l5y, d1_len, d1, &peyd, 8);
		double pezd_p[8], *pezd = pezd_p;
		int    pezd_len =
		  o.Gen_Product_With_PreAlloc(l5z_len, l5z, d1_len, d1, &pezd, 8);
		double ll1x_p[8], *ll1x = ll1x_p;
		int    ll1x_len =
		  o.Gen_Product_With_PreAlloc(l1x_len, l1x, d5_len, d5, &ll1x, 8);
		double ll1y_p[8], *ll1y = ll1y_p;
		int    ll1y_len =
		  o.Gen_Product_With_PreAlloc(l1y_len, l1y, d5_len, d5, &ll1y, 8);
		double ll1z_p[8], *ll1z = ll1z_p;
		int    ll1z_len =
		  o.Gen_Product_With_PreAlloc(l1z_len, l1z, d5_len, d5, &ll1z, 8);
		double aex_p[8], *aex = aex_p;
		int    aex_len =
		  o.Gen_Diff_With_PreAlloc(ll1x_len, ll1x, pexd_len, pexd, &aex, 8);
		double aey_p[8], *aey = aey_p;
		int    aey_len =
		  o.Gen_Diff_With_PreAlloc(ll1y_len, ll1y, peyd_len, peyd, &aey, 8);
		double aez_p[8], *aez = aez_p;
		int    aez_len =
		  o.Gen_Diff_With_PreAlloc(ll1z_len, ll1z, pezd_len, pezd, &aez, 8);
		double pexd2_p[8], *pexd2 = pexd2_p;
		int    pexd2_len =
		  o.Gen_Product_With_PreAlloc(l5x_len, l5x, d2_len, d2, &pexd2, 8);
		double peyd2_p[8], *peyd2 = peyd2_p;
		int    peyd2_len =
		  o.Gen_Product_With_PreAlloc(l5y_len, l5y, d2_len, d2, &peyd2, 8);
		double pezd2_p[8], *pezd2 = pezd2_p;
		int    pezd2_len =
		  o.Gen_Product_With_PreAlloc(l5z_len, l5z, d2_len, d2, &pezd2, 8);
		double ll2x_p[8], *ll2x = ll2x_p;
		int    ll2x_len =
		  o.Gen_Product_With_PreAlloc(l2x_len, l2x, d5_len, d5, &ll2x, 8);
		double ll2y_p[8], *ll2y = ll2y_p;
		int    ll2y_len =
		  o.Gen_Product_With_PreAlloc(l2y_len, l2y, d5_len, d5, &ll2y, 8);
		double ll2z_p[8], *ll2z = ll2z_p;
		int    ll2z_len =
		  o.Gen_Product_With_PreAlloc(l2z_len, l2z, d5_len, d5, &ll2z, 8);
		double bex_p[8], *bex = bex_p;
		int    bex_len =
		  o.Gen_Diff_With_PreAlloc(ll2x_len, ll2x, pexd2_len, pexd2, &bex, 8);
		double bey_p[8], *bey = bey_p;
		int    bey_len =
		  o.Gen_Diff_With_PreAlloc(ll2y_len, ll2y, peyd2_len, peyd2, &bey, 8);
		double bez_p[8], *bez = bez_p;
		int    bez_len =
		  o.Gen_Diff_With_PreAlloc(ll2z_len, ll2z, pezd2_len, pezd2, &bez, 8);
		double pexd3_p[8], *pexd3 = pexd3_p;
		int    pexd3_len =
		  o.Gen_Product_With_PreAlloc(l5x_len, l5x, d3_len, d3, &pexd3, 8);
		double peyd3_p[8], *peyd3 = peyd3_p;
		int    peyd3_len =
		  o.Gen_Product_With_PreAlloc(l5y_len, l5y, d3_len, d3, &peyd3, 8);
		double pezd3_p[8], *pezd3 = pezd3_p;
		int    pezd3_len =
		  o.Gen_Product_With_PreAlloc(l5z_len, l5z, d3_len, d3, &pezd3, 8);
		double ll3x_p[8], *ll3x = ll3x_p;
		int    ll3x_len =
		  o.Gen_Product_With_PreAlloc(l3x_len, l3x, d5_len, d5, &ll3x, 8);
		double ll3y_p[8], *ll3y = ll3y_p;
		int    ll3y_len =
		  o.Gen_Product_With_PreAlloc(l3y_len, l3y, d5_len, d5, &ll3y, 8);
		double ll3z_p[8], *ll3z = ll3z_p;
		int    ll3z_len =
		  o.Gen_Product_With_PreAlloc(l3z_len, l3z, d5_len, d5, &ll3z, 8);
		double cex_p[8], *cex = cex_p;
		int    cex_len =
		  o.Gen_Diff_With_PreAlloc(ll3x_len, ll3x, pexd3_len, pexd3, &cex, 8);
		double cey_p[8], *cey = cey_p;
		int    cey_len =
		  o.Gen_Diff_With_PreAlloc(ll3y_len, ll3y, peyd3_len, peyd3, &cey, 8);
		double cez_p[8], *cez = cez_p;
		int    cez_len =
		  o.Gen_Diff_With_PreAlloc(ll3z_len, ll3z, pezd3_len, pezd3, &cez, 8);
		double pexd4_p[8], *pexd4 = pexd4_p;
		int    pexd4_len =
		  o.Gen_Product_With_PreAlloc(l5x_len, l5x, d4_len, d4, &pexd4, 8);
		double peyd4_p[8], *peyd4 = peyd4_p;
		int    peyd4_len =
		  o.Gen_Product_With_PreAlloc(l5y_len, l5y, d4_len, d4, &peyd4, 8);
		double pezd4_p[8], *pezd4 = pezd4_p;
		int    pezd4_len =
		  o.Gen_Product_With_PreAlloc(l5z_len, l5z, d4_len, d4, &pezd4, 8);
		double ll4x_p[8], *ll4x = ll4x_p;
		int    ll4x_len =
		  o.Gen_Product_With_PreAlloc(l4x_len, l4x, d5_len, d5, &ll4x, 8);
		double ll4y_p[8], *ll4y = ll4y_p;
		int    ll4y_len =
		  o.Gen_Product_With_PreAlloc(l4y_len, l4y, d5_len, d5, &ll4y, 8);
		double ll4z_p[8], *ll4z = ll4z_p;
		int    ll4z_len =
		  o.Gen_Product_With_PreAlloc(l4z_len, l4z, d5_len, d5, &ll4z, 8);
		double dex_p[8], *dex = dex_p;
		int    dex_len =
		  o.Gen_Diff_With_PreAlloc(ll4x_len, ll4x, pexd4_len, pexd4, &dex, 8);
		double dey_p[8], *dey = dey_p;
		int    dey_len =
		  o.Gen_Diff_With_PreAlloc(ll4y_len, ll4y, peyd4_len, peyd4, &dey, 8);
		double dez_p[8], *dez = dez_p;
		int    dez_len =
		  o.Gen_Diff_With_PreAlloc(ll4z_len, ll4z, pezd4_len, pezd4, &dez, 8);
		double aexbey_p[8], *aexbey = aexbey_p;
		int    aexbey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, bey_len, bey, &aexbey, 8);
		double bexaey_p[8], *bexaey = bexaey_p;
		int    bexaey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, aey_len, aey, &bexaey, 8);
		double ab_p[8], *ab = ab_p;
		int    ab_len =
		  o.Gen_Diff_With_PreAlloc(aexbey_len, aexbey, bexaey_len, bexaey, &ab, 8);
		double bexcey_p[8], *bexcey = bexcey_p;
		int    bexcey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, cey_len, cey, &bexcey, 8);
		double cexbey_p[8], *cexbey = cexbey_p;
		int    cexbey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, bey_len, bey, &cexbey, 8);
		double bc_p[8], *bc = bc_p;
		int    bc_len =
		  o.Gen_Diff_With_PreAlloc(bexcey_len, bexcey, cexbey_len, cexbey, &bc, 8);
		double cexdey_p[8], *cexdey = cexdey_p;
		int    cexdey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, dey_len, dey, &cexdey, 8);
		double dexcey_p[8], *dexcey = dexcey_p;
		int    dexcey_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, cey_len, cey, &dexcey, 8);
		double cd_p[8], *cd = cd_p;
		int    cd_len =
		  o.Gen_Diff_With_PreAlloc(cexdey_len, cexdey, dexcey_len, dexcey, &cd, 8);
		double dexaey_p[8], *dexaey = dexaey_p;
		int    dexaey_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, aey_len, aey, &dexaey, 8);
		double aexdey_p[8], *aexdey = aexdey_p;
		int    aexdey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, dey_len, dey, &aexdey, 8);
		double da_p[8], *da = da_p;
		int    da_len =
		  o.Gen_Diff_With_PreAlloc(dexaey_len, dexaey, aexdey_len, aexdey, &da, 8);
		double aexcey_p[8], *aexcey = aexcey_p;
		int    aexcey_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, cey_len, cey, &aexcey, 8);
		double cexaey_p[8], *cexaey = cexaey_p;
		int    cexaey_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, aey_len, aey, &cexaey, 8);
		double ac_p[8], *ac = ac_p;
		int    ac_len =
		  o.Gen_Diff_With_PreAlloc(aexcey_len, aexcey, cexaey_len, cexaey, &ac, 8);
		double bexdey_p[8], *bexdey = bexdey_p;
		int    bexdey_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, dey_len, dey, &bexdey, 8);
		double dexbey_p[8], *dexbey = dexbey_p;
		int    dexbey_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, bey_len, bey, &dexbey, 8);
		double bd_p[8], *bd = bd_p;
		int    bd_len =
		  o.Gen_Diff_With_PreAlloc(bexdey_len, bexdey, dexbey_len, dexbey, &bd, 8);
		double abc1_p[8], *abc1 = abc1_p;
		int    abc1_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bc_len, bc, &abc1, 8);
		double abc2_p[8], *abc2 = abc2_p;
		int    abc2_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, ac_len, ac, &abc2, 8);
		double abc3_p[8], *abc3 = abc3_p;
		int    abc3_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, ab_len, ab, &abc3, 8);
		double abc4_p[8], *abc4 = abc4_p;
		int    abc4_len =
		  o.Gen_Sum_With_PreAlloc(abc1_len, abc1, abc3_len, abc3, &abc4, 8);
		double abc_p[8], *abc = abc_p;
		int    abc_len =
		  o.Gen_Diff_With_PreAlloc(abc4_len, abc4, abc2_len, abc2, &abc, 8);
		double bcd1_p[8], *bcd1 = bcd1_p;
		int    bcd1_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, cd_len, cd, &bcd1, 8);
		double bcd2_p[8], *bcd2 = bcd2_p;
		int    bcd2_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, bd_len, bd, &bcd2, 8);
		double bcd3_p[8], *bcd3 = bcd3_p;
		int    bcd3_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, bc_len, bc, &bcd3, 8);
		double bcd4_p[8], *bcd4 = bcd4_p;
		int    bcd4_len =
		  o.Gen_Sum_With_PreAlloc(bcd1_len, bcd1, bcd3_len, bcd3, &bcd4, 8);
		double bcd_p[8], *bcd = bcd_p;
		int    bcd_len =
		  o.Gen_Diff_With_PreAlloc(bcd4_len, bcd4, bcd2_len, bcd2, &bcd, 8);
		double cda1_p[8], *cda1 = cda1_p;
		int    cda1_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, da_len, da, &cda1, 8);
		double cda2_p[8], *cda2 = cda2_p;
		int    cda2_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, ac_len, ac, &cda2, 8);
		double cda3_p[8], *cda3 = cda3_p;
		int    cda3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, cd_len, cd, &cda3, 8);
		double cda4_p[8], *cda4 = cda4_p;
		int    cda4_len =
		  o.Gen_Sum_With_PreAlloc(cda1_len, cda1, cda3_len, cda3, &cda4, 8);
		double cda_p[8], *cda = cda_p;
		int    cda_len =
		  o.Gen_Sum_With_PreAlloc(cda4_len, cda4, cda2_len, cda2, &cda, 8);
		double dab1_p[8], *dab1 = dab1_p;
		int    dab1_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, ab_len, ab, &dab1, 8);
		double dab2_p[8], *dab2 = dab2_p;
		int    dab2_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, bd_len, bd, &dab2, 8);
		double dab3_p[8], *dab3 = dab3_p;
		int    dab3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, da_len, da, &dab3, 8);
		double dab4_p[8], *dab4 = dab4_p;
		int    dab4_len =
		  o.Gen_Sum_With_PreAlloc(dab1_len, dab1, dab3_len, dab3, &dab4, 8);
		double dab_p[8], *dab = dab_p;
		int    dab_len =
		  o.Gen_Sum_With_PreAlloc(dab4_len, dab4, dab2_len, dab2, &dab, 8);
		double al1_p[8], *al1 = al1_p;
		int    al1_len =
		  o.Gen_Product_With_PreAlloc(aex_len, aex, aex_len, aex, &al1, 8);
		double al2_p[8], *al2 = al2_p;
		int    al2_len =
		  o.Gen_Product_With_PreAlloc(aey_len, aey, aey_len, aey, &al2, 8);
		double al3_p[8], *al3 = al3_p;
		int    al3_len =
		  o.Gen_Product_With_PreAlloc(aez_len, aez, aez_len, aez, &al3, 8);
		double al4_p[8], *al4 = al4_p;
		int al4_len = o.Gen_Sum_With_PreAlloc(al1_len, al1, al2_len, al2, &al4, 8);
		double alift_p[8], *alift = alift_p;
		int    alift_len =
		  o.Gen_Sum_With_PreAlloc(al4_len, al4, al3_len, al3, &alift, 8);
		double bl1_p[8], *bl1 = bl1_p;
		int    bl1_len =
		  o.Gen_Product_With_PreAlloc(bex_len, bex, bex_len, bex, &bl1, 8);
		double bl2_p[8], *bl2 = bl2_p;
		int    bl2_len =
		  o.Gen_Product_With_PreAlloc(bey_len, bey, bey_len, bey, &bl2, 8);
		double bl3_p[8], *bl3 = bl3_p;
		int    bl3_len =
		  o.Gen_Product_With_PreAlloc(bez_len, bez, bez_len, bez, &bl3, 8);
		double bl4_p[8], *bl4 = bl4_p;
		int bl4_len = o.Gen_Sum_With_PreAlloc(bl1_len, bl1, bl2_len, bl2, &bl4, 8);
		double blift_p[8], *blift = blift_p;
		int    blift_len =
		  o.Gen_Sum_With_PreAlloc(bl4_len, bl4, bl3_len, bl3, &blift, 8);
		double cl1_p[8], *cl1 = cl1_p;
		int    cl1_len =
		  o.Gen_Product_With_PreAlloc(cex_len, cex, cex_len, cex, &cl1, 8);
		double cl2_p[8], *cl2 = cl2_p;
		int    cl2_len =
		  o.Gen_Product_With_PreAlloc(cey_len, cey, cey_len, cey, &cl2, 8);
		double cl3_p[8], *cl3 = cl3_p;
		int    cl3_len =
		  o.Gen_Product_With_PreAlloc(cez_len, cez, cez_len, cez, &cl3, 8);
		double cl4_p[8], *cl4 = cl4_p;
		int cl4_len = o.Gen_Sum_With_PreAlloc(cl1_len, cl1, cl2_len, cl2, &cl4, 8);
		double clift_p[8], *clift = clift_p;
		int    clift_len =
		  o.Gen_Sum_With_PreAlloc(cl4_len, cl4, cl3_len, cl3, &clift, 8);
		double dl1_p[8], *dl1 = dl1_p;
		int    dl1_len =
		  o.Gen_Product_With_PreAlloc(dex_len, dex, dex_len, dex, &dl1, 8);
		double dl2_p[8], *dl2 = dl2_p;
		int    dl2_len =
		  o.Gen_Product_With_PreAlloc(dey_len, dey, dey_len, dey, &dl2, 8);
		double dl3_p[8], *dl3 = dl3_p;
		int    dl3_len =
		  o.Gen_Product_With_PreAlloc(dez_len, dez, dez_len, dez, &dl3, 8);
		double dl4_p[8], *dl4 = dl4_p;
		int dl4_len = o.Gen_Sum_With_PreAlloc(dl1_len, dl1, dl2_len, dl2, &dl4, 8);
		double dlift_p[8], *dlift = dlift_p;
		int    dlift_len =
		  o.Gen_Sum_With_PreAlloc(dl4_len, dl4, dl3_len, dl3, &dlift, 8);
		double ds1_p[8], *ds1 = ds1_p;
		int    ds1_len =
		  o.Gen_Product_With_PreAlloc(dlift_len, dlift, abc_len, abc, &ds1, 8);
		double ds1n_p[8], *ds1n = ds1n_p;
		int    ds1n_len =
		  o.Gen_Product_With_PreAlloc(ds1_len, ds1, d3_len, d3, &ds1n, 8);
		double ds2_p[8], *ds2 = ds2_p;
		int    ds2_len =
		  o.Gen_Product_With_PreAlloc(clift_len, clift, dab_len, dab, &ds2, 8);
		double ds2n_p[8], *ds2n = ds2n_p;
		int    ds2n_len =
		  o.Gen_Product_With_PreAlloc(ds2_len, ds2, d4_len, d4, &ds2n, 8);
		double dl_p[8], *dl = dl_p;
		int    dl_len =
		  o.Gen_Diff_With_PreAlloc(ds2n_len, ds2n, ds1n_len, ds1n, &dl, 8);
		double dla_p[8], *dla = dla_p;
		int dla_len = o.Gen_Product_With_PreAlloc(dl_len, dl, d1_len, d1, &dla, 8);
		double dlb_p[8], *dlb = dlb_p;
		int    dlb_len =
		  o.Gen_Product_With_PreAlloc(dla_len, dla, d2_len, d2, &dlb, 8);
		double dr1_p[8], *dr1 = dr1_p;
		int    dr1_len =
		  o.Gen_Product_With_PreAlloc(blift_len, blift, cda_len, cda, &dr1, 8);
		double dr1n_p[8], *dr1n = dr1n_p;
		int    dr1n_len =
		  o.Gen_Product_With_PreAlloc(dr1_len, dr1, d1_len, d1, &dr1n, 8);
		double dr2_p[8], *dr2 = dr2_p;
		int    dr2_len =
		  o.Gen_Product_With_PreAlloc(alift_len, alift, bcd_len, bcd, &dr2, 8);
		double dr2n_p[8], *dr2n = dr2n_p;
		int    dr2n_len =
		  o.Gen_Product_With_PreAlloc(dr2_len, dr2, d2_len, d2, &dr2n, 8);
		double dr_p[8], *dr = dr_p;
		int    dr_len =
		  o.Gen_Diff_With_PreAlloc(dr2n_len, dr2n, dr1n_len, dr1n, &dr, 8);
		double dra_p[8], *dra = dra_p;
		int dra_len = o.Gen_Product_With_PreAlloc(dr_len, dr, d3_len, d3, &dra, 8);
		double drb_p[8], *drb = drb_p;
		int    drb_len =
		  o.Gen_Product_With_PreAlloc(dra_len, dra, d4_len, d4, &drb, 8);
		double det_p[8], *det = det_p;
		int det_len = o.Gen_Sum_With_PreAlloc(dlb_len, dlb, drb_len, drb, &det, 8);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (drb_p != drb)
			FreeDoubles(drb);
		if (dra_p != dra)
			FreeDoubles(dra);
		if (dr_p != dr)
			FreeDoubles(dr);
		if (dr2n_p != dr2n)
			FreeDoubles(dr2n);
		if (dr2_p != dr2)
			FreeDoubles(dr2);
		if (dr1n_p != dr1n)
			FreeDoubles(dr1n);
		if (dr1_p != dr1)
			FreeDoubles(dr1);
		if (dlb_p != dlb)
			FreeDoubles(dlb);
		if (dla_p != dla)
			FreeDoubles(dla);
		if (dl_p != dl)
			FreeDoubles(dl);
		if (ds2n_p != ds2n)
			FreeDoubles(ds2n);
		if (ds2_p != ds2)
			FreeDoubles(ds2);
		if (ds1n_p != ds1n)
			FreeDoubles(ds1n);
		if (ds1_p != ds1)
			FreeDoubles(ds1);
		if (dlift_p != dlift)
			FreeDoubles(dlift);
		if (dl4_p != dl4)
			FreeDoubles(dl4);
		if (dl3_p != dl3)
			FreeDoubles(dl3);
		if (dl2_p != dl2)
			FreeDoubles(dl2);
		if (dl1_p != dl1)
			FreeDoubles(dl1);
		if (clift_p != clift)
			FreeDoubles(clift);
		if (cl4_p != cl4)
			FreeDoubles(cl4);
		if (cl3_p != cl3)
			FreeDoubles(cl3);
		if (cl2_p != cl2)
			FreeDoubles(cl2);
		if (cl1_p != cl1)
			FreeDoubles(cl1);
		if (blift_p != blift)
			FreeDoubles(blift);
		if (bl4_p != bl4)
			FreeDoubles(bl4);
		if (bl3_p != bl3)
			FreeDoubles(bl3);
		if (bl2_p != bl2)
			FreeDoubles(bl2);
		if (bl1_p != bl1)
			FreeDoubles(bl1);
		if (alift_p != alift)
			FreeDoubles(alift);
		if (al4_p != al4)
			FreeDoubles(al4);
		if (al3_p != al3)
			FreeDoubles(al3);
		if (al2_p != al2)
			FreeDoubles(al2);
		if (al1_p != al1)
			FreeDoubles(al1);
		if (dab_p != dab)
			FreeDoubles(dab);
		if (dab4_p != dab4)
			FreeDoubles(dab4);
		if (dab3_p != dab3)
			FreeDoubles(dab3);
		if (dab2_p != dab2)
			FreeDoubles(dab2);
		if (dab1_p != dab1)
			FreeDoubles(dab1);
		if (cda_p != cda)
			FreeDoubles(cda);
		if (cda4_p != cda4)
			FreeDoubles(cda4);
		if (cda3_p != cda3)
			FreeDoubles(cda3);
		if (cda2_p != cda2)
			FreeDoubles(cda2);
		if (cda1_p != cda1)
			FreeDoubles(cda1);
		if (bcd_p != bcd)
			FreeDoubles(bcd);
		if (bcd4_p != bcd4)
			FreeDoubles(bcd4);
		if (bcd3_p != bcd3)
			FreeDoubles(bcd3);
		if (bcd2_p != bcd2)
			FreeDoubles(bcd2);
		if (bcd1_p != bcd1)
			FreeDoubles(bcd1);
		if (abc_p != abc)
			FreeDoubles(abc);
		if (abc4_p != abc4)
			FreeDoubles(abc4);
		if (abc3_p != abc3)
			FreeDoubles(abc3);
		if (abc2_p != abc2)
			FreeDoubles(abc2);
		if (abc1_p != abc1)
			FreeDoubles(abc1);
		if (bd_p != bd)
			FreeDoubles(bd);
		if (dexbey_p != dexbey)
			FreeDoubles(dexbey);
		if (bexdey_p != bexdey)
			FreeDoubles(bexdey);
		if (ac_p != ac)
			FreeDoubles(ac);
		if (cexaey_p != cexaey)
			FreeDoubles(cexaey);
		if (aexcey_p != aexcey)
			FreeDoubles(aexcey);
		if (da_p != da)
			FreeDoubles(da);
		if (aexdey_p != aexdey)
			FreeDoubles(aexdey);
		if (dexaey_p != dexaey)
			FreeDoubles(dexaey);
		if (cd_p != cd)
			FreeDoubles(cd);
		if (dexcey_p != dexcey)
			FreeDoubles(dexcey);
		if (cexdey_p != cexdey)
			FreeDoubles(cexdey);
		if (bc_p != bc)
			FreeDoubles(bc);
		if (cexbey_p != cexbey)
			FreeDoubles(cexbey);
		if (bexcey_p != bexcey)
			FreeDoubles(bexcey);
		if (ab_p != ab)
			FreeDoubles(ab);
		if (bexaey_p != bexaey)
			FreeDoubles(bexaey);
		if (aexbey_p != aexbey)
			FreeDoubles(aexbey);
		if (dez_p != dez)
			FreeDoubles(dez);
		if (dey_p != dey)
			FreeDoubles(dey);
		if (dex_p != dex)
			FreeDoubles(dex);
		if (ll4z_p != ll4z)
			FreeDoubles(ll4z);
		if (ll4y_p != ll4y)
			FreeDoubles(ll4y);
		if (ll4x_p != ll4x)
			FreeDoubles(ll4x);
		if (pezd4_p != pezd4)
			FreeDoubles(pezd4);
		if (peyd4_p != peyd4)
			FreeDoubles(peyd4);
		if (pexd4_p != pexd4)
			FreeDoubles(pexd4);
		if (cez_p != cez)
			FreeDoubles(cez);
		if (cey_p != cey)
			FreeDoubles(cey);
		if (cex_p != cex)
			FreeDoubles(cex);
		if (ll3z_p != ll3z)
			FreeDoubles(ll3z);
		if (ll3y_p != ll3y)
			FreeDoubles(ll3y);
		if (ll3x_p != ll3x)
			FreeDoubles(ll3x);
		if (pezd3_p != pezd3)
			FreeDoubles(pezd3);
		if (peyd3_p != peyd3)
			FreeDoubles(peyd3);
		if (pexd3_p != pexd3)
			FreeDoubles(pexd3);
		if (bez_p != bez)
			FreeDoubles(bez);
		if (bey_p != bey)
			FreeDoubles(bey);
		if (bex_p != bex)
			FreeDoubles(bex);
		if (ll2z_p != ll2z)
			FreeDoubles(ll2z);
		if (ll2y_p != ll2y)
			FreeDoubles(ll2y);
		if (ll2x_p != ll2x)
			FreeDoubles(ll2x);
		if (pezd2_p != pezd2)
			FreeDoubles(pezd2);
		if (peyd2_p != peyd2)
			FreeDoubles(peyd2);
		if (pexd2_p != pexd2)
			FreeDoubles(pexd2);
		if (aez_p != aez)
			FreeDoubles(aez);
		if (aey_p != aey)
			FreeDoubles(aey);
		if (aex_p != aex)
			FreeDoubles(aex);
		if (ll1z_p != ll1z)
			FreeDoubles(ll1z);
		if (ll1y_p != ll1y)
			FreeDoubles(ll1y);
		if (ll1x_p != ll1x)
			FreeDoubles(ll1x);
		if (pezd_p != pezd)
			FreeDoubles(pezd);
		if (peyd_p != peyd)
			FreeDoubles(peyd);
		if (pexd_p != pexd)
			FreeDoubles(pexd);
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

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l4x_p != l4x)
			FreeDoubles(l4x);
		if (l4y_p != l4y)
			FreeDoubles(l4y);
		if (l4z_p != l4z)
			FreeDoubles(l4z);
		if (d4_p != d4)
			FreeDoubles(d4);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l5x_p != l5x)
			FreeDoubles(l5x);
		if (l5y_p != l5y)
			FreeDoubles(l5y);
		if (l5z_p != l5z)
			FreeDoubles(l5z);
		if (d5_p != d5)
			FreeDoubles(d5);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return inSphere_IIIII_exact<IT, ET>(p1, p2, p3, p4, p5);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign inSphere_IIIII(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &p4,
                    const GenericPoint3T<IT, ET> &p5)
{
	Sign ret;
	ret = inSphere_IIIII_interval<IT, ET>(p1, p2, p3, p4, p5);
	if (is_sign_reliable(ret))
		return ret;
	return inSphere_IIIII_expansion<IT, ET>(p1, p2, p3, p4, p5);
}

template <typename IT, typename ET>
Sign lessThanOnX_IE_filtered(const GenericPoint3T<IT, ET> &p1, double bx,
                             PntArr3 arr)
{
	double l1x, l1y, l1z, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var))
		return Sign::UNCERTAIN;

	double dbx = bx * d1;
	double kx  = l1x - dbx;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(bx)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 4.440892098500627e-15;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.9322976378688418e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.9802709739245137e-13;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(kx, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnX_IE_interval(const GenericPoint3T<IT, ET> &p1, IT bx)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dbx = bx * d1;
	IT kx  = l1x - dbx;
	if (!kx.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_IE_exact(const GenericPoint3T<IT, ET> &p1, ET bx)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET dbx = bx * d1;
	ET kx  = l1x - dbx;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_IE_expansion(const GenericPoint3T<IT, ET> &p1, double bx)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          dbx_p[128], *dbx = dbx_p;
		int    dbx_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, bx, &dbx, 128);
		double kx_p[128], *kx = kx_p;
		int kx_len = o.Gen_Diff_With_PreAlloc(l1x_len, l1x, dbx_len, dbx, &kx, 128);

		return_value = kx[kx_len - 1];
		if (kx_p != kx)
			FreeDoubles(kx);
		if (dbx_p != dbx)
			FreeDoubles(dbx);
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
		return lessThanOnX_IE_exact<IT, ET>(p1, bx);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_IE(const GenericPoint3T<IT, ET> &p1, double bx, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnX_IE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnX_IE_filtered<IT, ET>(p1, bx, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnX_IE, arr);
	}
	ret = lessThanOnX_IE_interval<IT, ET>(p1, bx);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnX_IE, arr);
	return lessThanOnX_IE_expansion<IT, ET>(p1, bx);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	return lessThanOnX_IE<IT, ET, WithSSFilter>(p1, b.x(), arr);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var))
		return Sign::UNCERTAIN;

	double k1      = d2 * l1x;
	double k2      = d1 * l2x;
	double kx      = k1 - k2;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 2.2204460492503147e-14;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 7.906088977938107e-14;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.272171604171927e-12;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.922887626377606e-13;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.321380059346694e-12;
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
		epsilon *= 5.5041415869539175e-11;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(kx, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT k1 = d2 * l1x;
	IT k2 = d1 * l2x;
	IT kx = k1 - k2;
	if (!kx.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET k1 = d2 * l1x;
	ET k2 = d1 * l2x;
	ET kx = k1 - k2;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128],
	                   *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128],
	                   *l2z = l2z_p, d2_p[128], *d2 = d2_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128, l2x_len = 128,
	    l2y_len = 128, l2z_len = 128, d2_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          k1_p[128], *k1 = k1_p;
		int             k1_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &k1, 128);
		double k2_p[128], *k2 = k2_p;
		int    k2_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &k2, 128);
		double kx_p[128], *kx = kx_p;
		int    kx_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kx, 128);

		return_value = kx[kx_len - 1];
		if (kx_p != kx)
			FreeDoubles(kx);
		if (k2_p != k2)
			FreeDoubles(k2);
		if (k1_p != k1)
			FreeDoubles(k1);
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
		return lessThanOnX_II_exact<IT, ET>(p1, p2);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnX_II, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnX_II_filtered<IT, ET>(p1, p2, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnX_II, arr);
	}
	ret = lessThanOnX_II_interval<IT, ET>(p1, p2);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnX_II, arr);
	return lessThanOnX_II_expansion<IT, ET>(p1, p2);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_filtered(const GenericPoint3T<IT, ET> &p1, double by,
                             PntArr3 arr)
{
	double l1x, l1y, l1z, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var))
		return Sign::UNCERTAIN;

	double dby = by * d1;
	double ky  = l1y - dby;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(by)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 4.440892098500627e-15;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.9322976378688418e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.9802709739245137e-13;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(ky, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_interval(const GenericPoint3T<IT, ET> &p1, IT by)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dby = by * d1;
	IT ky  = l1y - dby;
	if (!ky.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_exact(const GenericPoint3T<IT, ET> &p1, ET by)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET dby = by * d1;
	ET ky  = l1y - dby;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_expansion(const GenericPoint3T<IT, ET> &p1, double by)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          dby_p[128], *dby = dby_p;
		int    dby_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, by, &dby, 128);
		double ky_p[128], *ky = ky_p;
		int ky_len = o.Gen_Diff_With_PreAlloc(l1y_len, l1y, dby_len, dby, &ky, 128);

		return_value = ky[ky_len - 1];
		if (ky_p != ky)
			FreeDoubles(ky);
		if (dby_p != dby)
			FreeDoubles(dby);
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
		return lessThanOnY_IE_exact<IT, ET>(p1, by);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_IE(const GenericPoint3T<IT, ET> &p1, double by, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnY_IE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnY_IE_filtered<IT, ET>(p1, by, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnY_IE, arr);
	}
	ret = lessThanOnY_IE_interval<IT, ET>(p1, by);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnY_IE, arr);
	return lessThanOnY_IE_expansion<IT, ET>(p1, by);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	return lessThanOnY_IE<IT, ET, WithSSFilter>(p1, b.y(), arr);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var))
		return Sign::UNCERTAIN;

	double k1      = d2 * l1y;
	double k2      = d1 * l2y;
	double ky      = k1 - k2;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 2.2204460492503147e-14;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 7.906088977938107e-14;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.272171604171927e-12;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.922887626377606e-13;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.321380059346694e-12;
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
		epsilon *= 5.5041415869539175e-11;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(ky, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT k1 = d2 * l1y;
	IT k2 = d1 * l2y;
	IT ky = k1 - k2;
	if (!ky.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET k1 = d2 * l1y;
	ET k2 = d1 * l2y;
	ET ky = k1 - k2;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128],
	                   *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128],
	                   *l2z = l2z_p, d2_p[128], *d2 = d2_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128, l2x_len = 128,
	    l2y_len = 128, l2z_len = 128, d2_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          k1_p[128], *k1 = k1_p;
		int             k1_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &k1, 128);
		double k2_p[128], *k2 = k2_p;
		int    k2_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &k2, 128);
		double ky_p[128], *ky = ky_p;
		int    ky_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &ky, 128);

		return_value = ky[ky_len - 1];
		if (ky_p != ky)
			FreeDoubles(ky);
		if (k2_p != k2)
			FreeDoubles(k2);
		if (k1_p != k1)
			FreeDoubles(k1);
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
		return lessThanOnY_II_exact<IT, ET>(p1, p2);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnY_II, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnY_II_filtered<IT, ET>(p1, p2, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnY_II, arr);
	}
	ret = lessThanOnY_II_interval<IT, ET>(p1, p2);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnY_II, arr);
	return lessThanOnY_II_expansion<IT, ET>(p1, p2);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_filtered(const GenericPoint3T<IT, ET> &p1, double bz,
                             PntArr3 arr)
{
	double l1x, l1y, l1z, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var))
		return Sign::UNCERTAIN;

	double dbz = bz * d1;
	double kz  = l1z - dbz;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(bz)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 4.440892098500627e-15;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.9322976378688418e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.9802709739245137e-13;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(kz, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_interval(const GenericPoint3T<IT, ET> &p1, IT bz)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dbz = bz * d1;
	IT kz  = l1z - dbz;
	if (!kz.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_exact(const GenericPoint3T<IT, ET> &p1, ET bz)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET dbz = bz * d1;
	ET kz  = l1z - dbz;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_expansion(const GenericPoint3T<IT, ET> &p1, double bz)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          dbz_p[128], *dbz = dbz_p;
		int    dbz_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, bz, &dbz, 128);
		double kz_p[128], *kz = kz_p;
		int kz_len = o.Gen_Diff_With_PreAlloc(l1z_len, l1z, dbz_len, dbz, &kz, 128);

		return_value = kz[kz_len - 1];
		if (kz_p != kz)
			FreeDoubles(kz);
		if (dbz_p != dbz)
			FreeDoubles(dbz);
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
		return lessThanOnZ_IE_exact<IT, ET>(p1, bz);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_IE(const GenericPoint3T<IT, ET> &p1, double bz, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnZ_IE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnZ_IE_filtered<IT, ET>(p1, bz, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnZ_IE, arr);
	}
	ret = lessThanOnZ_IE_interval<IT, ET>(p1, bz);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnZ_IE, arr);
	return lessThanOnZ_IE_expansion<IT, ET>(p1, bz);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	return lessThanOnZ_IE<IT, ET, WithSSFilter>(p1, b.z(), arr);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var))
		return Sign::UNCERTAIN;

	double k1      = d2 * l1z;
	double k2      = d1 * l2z;
	double kz      = k1 - k2;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 2.2204460492503147e-14;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 7.906088977938107e-14;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.272171604171927e-12;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.922887626377606e-13;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.321380059346694e-12;
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
		epsilon *= 5.5041415869539175e-11;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(kz, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT k1 = d2 * l1z;
	IT k2 = d1 * l2z;
	IT kz = k1 - k2;
	if (!kz.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET k1 = d2 * l1z;
	ET k2 = d1 * l2z;
	ET kz = k1 - k2;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, l2x_p[128],
	                   *l2x = l2x_p, l2y_p[128], *l2y = l2y_p, l2z_p[128],
	                   *l2z = l2z_p, d2_p[128], *d2 = d2_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128, l2x_len = 128,
	    l2y_len = 128, l2z_len = 128, d2_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          k1_p[128], *k1 = k1_p;
		int             k1_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l1z_len, l1z, &k1, 128);
		double k2_p[128], *k2 = k2_p;
		int    k2_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l2z_len, l2z, &k2, 128);
		double kz_p[128], *kz = kz_p;
		int    kz_len = o.Gen_Diff_With_PreAlloc(k1_len, k1, k2_len, k2, &kz, 128);

		return_value = kz[kz_len - 1];
		if (kz_p != kz)
			FreeDoubles(kz);
		if (k2_p != k2)
			FreeDoubles(k2);
		if (k1_p != k1)
			FreeDoubles(k1);
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
		return lessThanOnZ_II_exact<IT, ET>(p1, p2);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnZ_II, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnZ_II_filtered<IT, ET>(p1, p2, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnZ_II, arr);
	}
	ret = lessThanOnZ_II_interval<IT, ET>(p1, p2);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnZ_II, arr);
	return lessThanOnZ_II_expansion<IT, ET>(p1, p2);
}

template <typename IT, typename ET>
Sign orient2D_IEE_filtered(const GenericPoint2T<IT, ET> &p1, double p2x,
                           double p2y, double p3x, double p3y, PntArr2 arr)
{
	double l1x, l1y, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, d1, max_var))
		return Sign::UNCERTAIN;

	double t1x = p2y - p3y;
	double t1y = p3x - p2x;
	double e2  = l1x * t1x;
	double e3  = l1y * t1y;
	double e   = e2 + e3;
	double pr1 = p2x * p3y;
	double pr2 = p2y * p3x;
	double pr  = pr1 - pr2;
	double dpr = d1 * pr;
	double det = dpr + e;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(p2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1y)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr2::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.1104398650596543e-14;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orient2D_IEE_interval(const GenericPoint2T<IT, ET> &p1, IT p2x, IT p2y,
                           IT p3x, IT p3y)
{
	IT l1x, l1y, d1;
	if (!p1.getIntervalLambda(l1x, l1y, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t1x = p2y - p3y;
	IT t1y = p3x - p2x;
	IT e2  = l1x * t1x;
	IT e3  = l1y * t1y;
	IT e   = e2 + e3;
	IT pr1 = p2x * p3y;
	IT pr2 = p2y * p3x;
	IT pr  = pr1 - pr2;
	IT dpr = d1 * pr;
	IT det = dpr + e;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orient2D_IEE_exact(const GenericPoint2T<IT, ET> &p1, ET p2x, ET p2y,
                        ET p3x, ET p3y)
{
	ET l1x, l1y, d1;
	p1.getExactLambda(l1x, l1y, d1);
	ET t1x = p2y - p3y;
	ET t1y = p3x - p2x;
	ET e2  = l1x * t1x;
	ET e3  = l1y * t1y;
	ET e   = e2 + e3;
	ET pr1 = p2x * p3y;
	ET pr2 = p2y * p3x;
	ET pr  = pr1 - pr2;
	ET dpr = d1 * pr;
	ET det = dpr + e;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orient2D_IEE_expansion(const GenericPoint2T<IT, ET> &p1, double p2x,
                            double p2y, double p3x, double p3y)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, d1_p[128],
	                   *d1 = d1_p;
	int l1x_len = 128, l1y_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          t1x[2];
		o.two_Diff(p2y, p3y, t1x);
		double t1y[2];
		o.two_Diff(p3x, p2x, t1y);
		double e2_p[128], *e2 = e2_p;
		int    e2_len = o.Gen_Product_With_PreAlloc(l1x_len, l1x, 2, t1x, &e2, 128);
		double e3_p[128], *e3 = e3_p;
		int    e3_len = o.Gen_Product_With_PreAlloc(l1y_len, l1y, 2, t1y, &e3, 128);
		double e_p[128], *e = e_p;
		int    e_len = o.Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
		double pr1[2];
		o.Two_Prod(p2x, p3y, pr1);
		double pr2[2];
		o.Two_Prod(p2y, p3x, pr2);
		double pr[4];
		o.Two_Two_Diff(pr1, pr2, pr);
		double dpr_p[128], *dpr = dpr_p;
		int    dpr_len = o.Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
		double det_p[128], *det = det_p;
		int    det_len = o.Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (dpr_p != dpr)
			FreeDoubles(dpr);
		if (e_p != e)
			FreeDoubles(e);
		if (e3_p != e3)
			FreeDoubles(e3);
		if (e2_p != e2)
			FreeDoubles(e2);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orient2D_IEE_exact<IT, ET>(p1, p2x, p2y, p3x, p3y);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IEE(const GenericPoint2T<IT, ET> &p1, double p2x, double p2y,
                  double p3x, double p3y, PntArr2 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient2D_IEE_filtered<IT, ET>(p1, p2x, p2y, p3x, p3y, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient2D_IEE_interval<IT, ET>(p1, p2x, p2y, p3x, p3y);
	if (is_sign_reliable(ret))
		return ret;
	return orient2D_IEE_expansion<IT, ET>(p1, p2x, p2y, p3x, p3y);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IEE(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2,
                  const GenericPoint2T<IT, ET> &p3, PntArr2 arr)
{
	return orient2D_IEE<IT, ET, WithSSFilter>(p1, p2.x(), p2.y(), p3.x(), p3.y(),
	                                          arr);
}

template <typename IT, typename ET>
Sign orient2D_IIE_filtered(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2, double p3x,
                           double p3y, PntArr2 arr)
{
	double l1x, l1y, d1, l2x, l2y, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, d2, max_var))
		return Sign::UNCERTAIN;

	double a    = d1 * l2x;
	double b    = d2 * l1x;
	double c    = d1 * p3y;
	double e    = d1 * l2y;
	double f    = d2 * l1y;
	double g    = d1 * p3x;
	double ab   = a - b;
	double cd   = c - l1y;
	double ef   = e - f;
	double gh   = g - l1x;
	double abcd = ab * cd;
	double efgh = ef * gh;
	double L    = abcd - efgh;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3y)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr2::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 3.837902218251096e-13;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orient2D_IIE_interval(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2, IT p3x, IT p3y)
{
	IT l1x, l1y, d1, l2x, l2y, d2;
	if (!p1.getIntervalLambda(l1x, l1y, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2x;
	IT b    = d2 * l1x;
	IT c    = d1 * p3y;
	IT e    = d1 * l2y;
	IT f    = d2 * l1y;
	IT g    = d1 * p3x;
	IT ab   = a - b;
	IT cd   = c - l1y;
	IT ef   = e - f;
	IT gh   = g - l1x;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orient2D_IIE_exact(const GenericPoint2T<IT, ET> &p1,
                        const GenericPoint2T<IT, ET> &p2, ET p3x, ET p3y)
{
	ET l1x, l1y, d1, l2x, l2y, d2;
	p1.getExactLambda(l1x, l1y, d1);
	p2.getExactLambda(l2x, l2y, d2);
	ET a    = d1 * l2x;
	ET b    = d2 * l1x;
	ET c    = d1 * p3y;
	ET e    = d1 * l2y;
	ET f    = d2 * l1y;
	ET g    = d1 * p3x;
	ET ab   = a - b;
	ET cd   = c - l1y;
	ET ef   = e - f;
	ET gh   = g - l1x;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orient2D_IIE_expansion(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2, double p3x,
                            double p3y)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, d1_p[64], *d1 = d1_p,
	                  l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, d2_p[64],
	                  *d2 = d2_p;
	int l1x_len = 64, l1y_len = 64, d1_len = 64, l2x_len = 64, l2y_len = 64,
	    d2_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &d2, d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          a_p[64], *a = a_p;
		int a_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2x_len, l2x, &a, 64);
		double b_p[64], *b = b_p;
		int b_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1x_len, l1x, &b, 64);
		double c_p[64], *c = c_p;
		int    c_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3y, &c, 64);
		double e_p[64], *e = e_p;
		int e_len = o.Gen_Product_With_PreAlloc(d1_len, d1, l2y_len, l2y, &e, 64);
		double f_p[64], *f = f_p;
		int f_len = o.Gen_Product_With_PreAlloc(d2_len, d2, l1y_len, l1y, &f, 64);
		double g_p[64], *g = g_p;
		int    g_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p3x, &g, 64);
		double ab_p[64], *ab = ab_p;
		int    ab_len = o.Gen_Diff_With_PreAlloc(a_len, a, b_len, b, &ab, 64);
		double cd_p[64], *cd = cd_p;
		int    cd_len = o.Gen_Diff_With_PreAlloc(c_len, c, l1y_len, l1y, &cd, 64);
		double ef_p[64], *ef = ef_p;
		int    ef_len = o.Gen_Diff_With_PreAlloc(e_len, e, f_len, f, &ef, 64);
		double gh_p[64], *gh = gh_p;
		int    gh_len = o.Gen_Diff_With_PreAlloc(g_len, g, l1x_len, l1x, &gh, 64);
		double abcd_p[64], *abcd = abcd_p;
		int    abcd_len =
		  o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
		double efgh_p[64], *efgh = efgh_p;
		int    efgh_len =
		  o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
		double L_p[64], *L = L_p;
		int    L_len =
		  o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

		return_value = L[L_len - 1];
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

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orient2D_IIE_exact<IT, ET>(p1, p2, p3x, p3y);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IIE(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2, double p3x, double p3y,
                  PntArr2 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient2D_IIE_filtered<IT, ET>(p1, p2, p3x, p3y, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient2D_IIE_interval<IT, ET>(p1, p2, p3x, p3y);
	if (is_sign_reliable(ret))
		return ret;
	return orient2D_IIE_expansion<IT, ET>(p1, p2, p3x, p3y);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IIE(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2,
                  const GenericPoint2T<IT, ET> &p3, PntArr2 arr)
{
	return orient2D_IIE<IT, ET, WithSSFilter>(p1, p2, p3.x(), p3.y(), arr);
}

template <typename IT, typename ET>
Sign orient2D_III_filtered(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2,
                           const GenericPoint2T<IT, ET> &p3, PntArr2 arr)
{
	double l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, d2, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, d3, max_var))
		return Sign::UNCERTAIN;

	double a       = d1 * l2x;
	double b       = d2 * l1x;
	double c       = d1 * l3y;
	double d       = d3 * l1y;
	double e       = d1 * l2y;
	double f       = d2 * l1y;
	double g       = d1 * l3x;
	double h       = d3 * l1x;
	double ab      = a - b;
	double cd      = c - d;
	double ef      = e - f;
	double gh      = g - h;
	double abcd    = ab * cd;
	double efgh    = ef * gh;
	double L       = abcd - efgh;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr2::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.3645751195667837e-12;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orient2D_III_interval(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2,
                           const GenericPoint2T<IT, ET> &p3)
{
	IT l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
	if (!p1.getIntervalLambda(l1x, l1y, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2x;
	IT b    = d2 * l1x;
	IT c    = d1 * l3y;
	IT d    = d3 * l1y;
	IT e    = d1 * l2y;
	IT f    = d2 * l1y;
	IT g    = d1 * l3x;
	IT h    = d3 * l1x;
	IT ab   = a - b;
	IT cd   = c - d;
	IT ef   = e - f;
	IT gh   = g - h;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orient2D_III_exact(const GenericPoint2T<IT, ET> &p1,
                        const GenericPoint2T<IT, ET> &p2,
                        const GenericPoint2T<IT, ET> &p3)
{
	ET l1x, l1y, d1, l2x, l2y, d2, l3x, l3y, d3;
	p1.getExactLambda(l1x, l1y, d1);
	p2.getExactLambda(l2x, l2y, d2);
	p3.getExactLambda(l3x, l3y, d3);
	ET a    = d1 * l2x;
	ET b    = d2 * l1x;
	ET c    = d1 * l3y;
	ET d    = d3 * l1y;
	ET e    = d1 * l2y;
	ET f    = d2 * l1y;
	ET g    = d1 * l3x;
	ET h    = d3 * l1x;
	ET ab   = a - b;
	ET cd   = c - d;
	ET ef   = e - f;
	ET gh   = g - h;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orient2D_III_expansion(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2,
                            const GenericPoint2T<IT, ET> &p3)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, d1_p[64], *d1 = d1_p,
	                  l2x_p[64], *l2x = l2x_p, l2y_p[64], *l2y = l2y_p, d2_p[64],
	                  *d2 = d2_p, l3x_p[64], *l3x = l3x_p, l3y_p[64],
	                  *l3y = l3y_p, d3_p[64], *d3 = d3_p;
	int l1x_len = 64, l1y_len = 64, d1_len = 64, l2x_len = 64, l2y_len = 64,
	    d2_len = 64, l3x_len = 64, l3y_len = 64, d3_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &d2, d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &d3, d3_len);
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
		int    abcd_len =
		  o.Gen_Product_With_PreAlloc(ab_len, ab, cd_len, cd, &abcd, 64);
		double efgh_p[64], *efgh = efgh_p;
		int    efgh_len =
		  o.Gen_Product_With_PreAlloc(ef_len, ef, gh_len, gh, &efgh, 64);
		double L_p[64], *L = L_p;
		int    L_len =
		  o.Gen_Diff_With_PreAlloc(abcd_len, abcd, efgh_len, efgh, &L, 64);

		return_value = L[L_len - 1];
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

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l1x_p != l1x)
			FreeDoubles(l1x);
		if (l1y_p != l1y)
			FreeDoubles(l1y);
		if (d1_p != d1)
			FreeDoubles(d1);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l2x_p != l2x)
			FreeDoubles(l2x);
		if (l2y_p != l2y)
			FreeDoubles(l2y);
		if (d2_p != d2)
			FreeDoubles(d2);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (l3x_p != l3x)
			FreeDoubles(l3x);
		if (l3y_p != l3y)
			FreeDoubles(l3y);
		if (d3_p != d3)
			FreeDoubles(d3);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orient2D_III_exact<IT, ET>(p1, p2, p3);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_III(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2,
                  const GenericPoint2T<IT, ET> &p3, PntArr2 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient2D_III_filtered<IT, ET>(p1, p2, p3, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient2D_III_interval<IT, ET>(p1, p2, p3);
	if (is_sign_reliable(ret))
		return ret;
	return orient2D_III_expansion<IT, ET>(p1, p2, p3);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_filtered(const GenericPoint3T<IT, ET> &p1, double ax,
                            double ay, double az, double bx, double by,
                            double bz, double cx, double cy, double cz,
                            PntArr3 arr)
{
	double l1x, l1y, l1z, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var))
		return Sign::UNCERTAIN;

	double dcx   = d1 * cx;
	double dcy   = d1 * cy;
	double dcz   = d1 * cz;
	double ix_cx = l1x - dcx;
	double iy_cy = l1y - dcy;
	double ax_cx = ax - cx;
	double ay_cy = ay - cy;
	double az_cz = az - cz;
	double iz_cz = l1z - dcz;
	double bx_cx = bx - cx;
	double by_cy = by - cy;
	double bz_cz = bz - cz;
	double tmc_a = ix_cx * ay_cy;
	double tmc_b = iy_cy * ax_cx;
	double m01   = tmc_a - tmc_b;
	double tmi_a = ix_cx * az_cz;
	double tmi_b = iz_cz * ax_cx;
	double m02   = tmi_a - tmi_b;
	double tma_a = iy_cy * az_cz;
	double tma_b = iz_cz * ay_cy;
	double m12   = tma_a - tma_b;
	double mt1   = m01 * bz_cz;
	double mt2   = m02 * by_cy;
	double mt3   = m12 * bx_cx;
	double mtt   = mt2 - mt1;
	double m012  = mtt - mt3;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(cx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(cy)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(cz)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ax_cx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(ay_cy)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(az_cz)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(bx_cx)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(by_cy)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(bz_cz)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 4.7961634663806794e-14;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.8610395342844048e-13;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.0702836106844058e-12;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_interval(const GenericPoint3T<IT, ET> &p1, IT ax, IT ay,
                            IT az, IT bx, IT by, IT bz, IT cx, IT cy, IT cz)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dcx   = d1 * cx;
	IT dcy   = d1 * cy;
	IT dcz   = d1 * cz;
	IT ix_cx = l1x - dcx;
	IT iy_cy = l1y - dcy;
	IT ax_cx = ax - cx;
	IT ay_cy = ay - cy;
	IT az_cz = az - cz;
	IT iz_cz = l1z - dcz;
	IT bx_cx = bx - cx;
	IT by_cy = by - cy;
	IT bz_cz = bz - cz;
	IT tmc_a = ix_cx * ay_cy;
	IT tmc_b = iy_cy * ax_cx;
	IT m01   = tmc_a - tmc_b;
	IT tmi_a = ix_cx * az_cz;
	IT tmi_b = iz_cz * ax_cx;
	IT m02   = tmi_a - tmi_b;
	IT tma_a = iy_cy * az_cz;
	IT tma_b = iz_cz * ay_cy;
	IT m12   = tma_a - tma_b;
	IT mt1   = m01 * bz_cz;
	IT mt2   = m02 * by_cy;
	IT mt3   = m12 * bx_cx;
	IT mtt   = mt2 - mt1;
	IT m012  = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_exact(const GenericPoint3T<IT, ET> &p1, ET ax, ET ay, ET az,
                         ET bx, ET by, ET bz, ET cx, ET cy, ET cz)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET dcx   = d1 * cx;
	ET dcy   = d1 * cy;
	ET dcz   = d1 * cz;
	ET ix_cx = l1x - dcx;
	ET iy_cy = l1y - dcy;
	ET ax_cx = ax - cx;
	ET ay_cy = ay - cy;
	ET az_cz = az - cz;
	ET iz_cz = l1z - dcz;
	ET bx_cx = bx - cx;
	ET by_cy = by - cy;
	ET bz_cz = bz - cz;
	ET tmc_a = ix_cx * ay_cy;
	ET tmc_b = iy_cy * ax_cx;
	ET m01   = tmc_a - tmc_b;
	ET tmi_a = ix_cx * az_cz;
	ET tmi_b = iz_cz * ax_cx;
	ET m02   = tmi_a - tmi_b;
	ET tma_a = iy_cy * az_cz;
	ET tma_b = iz_cz * ay_cy;
	ET m12   = tma_a - tma_b;
	ET mt1   = m01 * bz_cz;
	ET mt2   = m02 * by_cy;
	ET mt3   = m12 * bx_cx;
	ET mtt   = mt2 - mt1;
	ET m012  = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_expansion(const GenericPoint3T<IT, ET> &p1, double ax,
                             double ay, double az, double bx, double by,
                             double bz, double cx, double cy, double cz)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          dcx_p[64], *dcx = dcx_p;
		int    dcx_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cx, &dcx, 64);
		double dcy_p[64], *dcy = dcy_p;
		int    dcy_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cy, &dcy, 64);
		double dcz_p[64], *dcz = dcz_p;
		int    dcz_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, cz, &dcz, 64);
		double ix_cx_p[64], *ix_cx = ix_cx_p;
		int    ix_cx_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, dcx_len, dcx, &ix_cx, 64);
		double iy_cy_p[64], *iy_cy = iy_cy_p;
		int    iy_cy_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, dcy_len, dcy, &iy_cy, 64);
		double ax_cx[2];
		o.two_Diff(ax, cx, ax_cx);
		double ay_cy[2];
		o.two_Diff(ay, cy, ay_cy);
		double az_cz[2];
		o.two_Diff(az, cz, az_cz);
		double iz_cz_p[64], *iz_cz = iz_cz_p;
		int    iz_cz_len =
		  o.Gen_Diff_With_PreAlloc(l1z_len, l1z, dcz_len, dcz, &iz_cz, 64);
		double bx_cx[2];
		o.two_Diff(bx, cx, bx_cx);
		double by_cy[2];
		o.two_Diff(by, cy, by_cy);
		double bz_cz[2];
		o.two_Diff(bz, cz, bz_cz);
		double tmc_a_p[64], *tmc_a = tmc_a_p;
		int    tmc_a_len =
		  o.Gen_Product_With_PreAlloc(ix_cx_len, ix_cx, 2, ay_cy, &tmc_a, 64);
		double tmc_b_p[64], *tmc_b = tmc_b_p;
		int    tmc_b_len =
		  o.Gen_Product_With_PreAlloc(iy_cy_len, iy_cy, 2, ax_cx, &tmc_b, 64);
		double m01_p[64], *m01 = m01_p;
		int    m01_len =
		  o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 64);
		double tmi_a_p[64], *tmi_a = tmi_a_p;
		int    tmi_a_len =
		  o.Gen_Product_With_PreAlloc(ix_cx_len, ix_cx, 2, az_cz, &tmi_a, 64);
		double tmi_b_p[64], *tmi_b = tmi_b_p;
		int    tmi_b_len =
		  o.Gen_Product_With_PreAlloc(iz_cz_len, iz_cz, 2, ax_cx, &tmi_b, 64);
		double m02_p[64], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 64);
		double tma_a_p[64], *tma_a = tma_a_p;
		int    tma_a_len =
		  o.Gen_Product_With_PreAlloc(iy_cy_len, iy_cy, 2, az_cz, &tma_a, 64);
		double tma_b_p[64], *tma_b = tma_b_p;
		int    tma_b_len =
		  o.Gen_Product_With_PreAlloc(iz_cz_len, iz_cz, 2, ay_cy, &tma_b, 64);
		double m12_p[64], *m12 = m12_p;
		int    m12_len =
		  o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 64);
		double mt1_p[64], *mt1 = mt1_p;
		int mt1_len = o.Gen_Product_With_PreAlloc(m01_len, m01, 2, bz_cz, &mt1, 64);
		double mt2_p[64], *mt2 = mt2_p;
		int mt2_len = o.Gen_Product_With_PreAlloc(m02_len, m02, 2, by_cy, &mt2, 64);
		double mt3_p[64], *mt3 = mt3_p;
		int mt3_len = o.Gen_Product_With_PreAlloc(m12_len, m12, 2, bx_cx, &mt3, 64);
		double mtt_p[64], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 64);
		double m012_p[64], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 64);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tma_b_p != tma_b)
			FreeDoubles(tma_b);
		if (tma_a_p != tma_a)
			FreeDoubles(tma_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmc_b_p != tmc_b)
			FreeDoubles(tmc_b);
		if (tmc_a_p != tmc_a)
			FreeDoubles(tmc_a);
		if (iz_cz_p != iz_cz)
			FreeDoubles(iz_cz);
		if (iy_cy_p != iy_cy)
			FreeDoubles(iy_cy);
		if (ix_cx_p != ix_cx)
			FreeDoubles(ix_cx);
		if (dcz_p != dcz)
			FreeDoubles(dcz);
		if (dcy_p != dcy)
			FreeDoubles(dcy);
		if (dcx_p != dcx)
			FreeDoubles(dcx);
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
		return orient3D_IEEE_exact<IT, ET>(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IEEE(const GenericPoint3T<IT, ET> &p1, double ax, double ay,
                   double az, double bx, double by, double bz, double cx,
                   double cy, double cz, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient3D_IEEE_filtered<IT, ET>(p1, ax, ay, az, bx, by, bz, cx, cy, cz,
		                                     arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IEEE_interval<IT, ET>(p1, ax, ay, az, bx, by, bz, cx, cy, cz);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IEEE_expansion<IT, ET>(p1, ax, ay, az, bx, by, bz, cx, cy,
	                                       cz);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IEEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &a,
                   const GenericPoint3T<IT, ET> &b,
                   const GenericPoint3T<IT, ET> &c, PntArr3 arr)
{
	return orient3D_IEEE<IT, ET, WithSSFilter>(
	  p1, a.x(), a.y(), a.z(), b.x(), b.y(), b.z(), c.x(), c.y(), c.z(), arr);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, double p3x,
                            double p3y, double p3z, double p4x, double p4y,
                            double p4z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var))
		return Sign::UNCERTAIN;

	double d1p4x = d1 * p4x;
	double d1p4y = d1 * p4y;
	double d1p4z = d1 * p4z;
	double d2p4x = d2 * p4x;
	double d2p4y = d2 * p4y;
	double d2p4z = d2 * p4z;
	double p1p4x = l1x - d1p4x;
	double p1p4y = l1y - d1p4y;
	double p1p4z = l1z - d1p4z;
	double p2p4x = l2x - d2p4x;
	double p2p4y = l2y - d2p4y;
	double p2p4z = l2z - d2p4z;
	double p3p4x = p3x - p4x;
	double p3p4y = p3y - p4y;
	double p3p4z = p3z - p4z;
	double tmc_a = p1p4x * p2p4y;
	double tmc_b = p1p4y * p2p4x;
	double m01   = tmc_a - tmc_b;
	double tmi_a = p1p4x * p2p4z;
	double tmi_b = p1p4z * p2p4x;
	double m02   = tmi_a - tmi_b;
	double tma_a = p1p4y * p2p4z;
	double tma_b = p1p4z * p2p4y;
	double m12   = tma_a - tma_b;
	double mt1   = m01 * p3p4z;
	double mt2   = m02 * p3p4y;
	double mt3   = m12 * p3p4x;
	double mtt   = mt2 - mt1;
	double m012  = mtt - mt3;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3p4z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.3609560407276204e-13;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.4603487689945929e-12;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.316812713898393e-11;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 5.12855469897434e-12;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 7.437036403379365e-11;
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
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.0361982383244646e-09;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, IT p3x, IT p3y,
                            IT p3z, IT p4x, IT p4y, IT p4z)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT d1p4x = d1 * p4x;
	IT d1p4y = d1 * p4y;
	IT d1p4z = d1 * p4z;
	IT d2p4x = d2 * p4x;
	IT d2p4y = d2 * p4y;
	IT d2p4z = d2 * p4z;
	IT p1p4x = l1x - d1p4x;
	IT p1p4y = l1y - d1p4y;
	IT p1p4z = l1z - d1p4z;
	IT p2p4x = l2x - d2p4x;
	IT p2p4y = l2y - d2p4y;
	IT p2p4z = l2z - d2p4z;
	IT p3p4x = p3x - p4x;
	IT p3p4y = p3y - p4y;
	IT p3p4z = p3z - p4z;
	IT tmc_a = p1p4x * p2p4y;
	IT tmc_b = p1p4y * p2p4x;
	IT m01   = tmc_a - tmc_b;
	IT tmi_a = p1p4x * p2p4z;
	IT tmi_b = p1p4z * p2p4x;
	IT m02   = tmi_a - tmi_b;
	IT tma_a = p1p4y * p2p4z;
	IT tma_b = p1p4z * p2p4y;
	IT m12   = tma_a - tma_b;
	IT mt1   = m01 * p3p4z;
	IT mt2   = m02 * p3p4y;
	IT mt3   = m12 * p3p4x;
	IT mtt   = mt2 - mt1;
	IT m012  = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2, ET p3x, ET p3y,
                         ET p3z, ET p4x, ET p4y, ET p4z)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET d1p4x = d1 * p4x;
	ET d1p4y = d1 * p4y;
	ET d1p4z = d1 * p4z;
	ET d2p4x = d2 * p4x;
	ET d2p4y = d2 * p4y;
	ET d2p4z = d2 * p4z;
	ET p1p4x = l1x - d1p4x;
	ET p1p4y = l1y - d1p4y;
	ET p1p4z = l1z - d1p4z;
	ET p2p4x = l2x - d2p4x;
	ET p2p4y = l2y - d2p4y;
	ET p2p4z = l2z - d2p4z;
	ET p3p4x = p3x - p4x;
	ET p3p4y = p3y - p4y;
	ET p3p4z = p3z - p4z;
	ET tmc_a = p1p4x * p2p4y;
	ET tmc_b = p1p4y * p2p4x;
	ET m01   = tmc_a - tmc_b;
	ET tmi_a = p1p4x * p2p4z;
	ET tmi_b = p1p4z * p2p4x;
	ET m02   = tmi_a - tmi_b;
	ET tma_a = p1p4y * p2p4z;
	ET tma_b = p1p4z * p2p4y;
	ET m12   = tma_a - tma_b;
	ET mt1   = m01 * p3p4z;
	ET mt2   = m02 * p3p4y;
	ET mt3   = m12 * p3p4x;
	ET mtt   = mt2 - mt1;
	ET m012  = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, double p3x,
                             double p3y, double p3z, double p4x, double p4y,
                             double p4z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32],
	                  *l1z = l1z_p, d1_p[32], *d1 = d1_p, l2x_p[32], *l2x = l2x_p,
	                  l2y_p[32], *l2y = l2y_p, l2z_p[32], *l2z = l2z_p, d2_p[32],
	                  *d2 = d2_p;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          d1p4x_p[32], *d1p4x = d1p4x_p;
		int    d1p4x_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p4x, &d1p4x, 32);
		double d1p4y_p[32], *d1p4y = d1p4y_p;
		int    d1p4y_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p4y, &d1p4y, 32);
		double d1p4z_p[32], *d1p4z = d1p4z_p;
		int    d1p4z_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p4z, &d1p4z, 32);
		double d2p4x_p[32], *d2p4x = d2p4x_p;
		int    d2p4x_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, p4x, &d2p4x, 32);
		double d2p4y_p[32], *d2p4y = d2p4y_p;
		int    d2p4y_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, p4y, &d2p4y, 32);
		double d2p4z_p[32], *d2p4z = d2p4z_p;
		int    d2p4z_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, p4z, &d2p4z, 32);
		double p1p4x_p[32], *p1p4x = p1p4x_p;
		int    p1p4x_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, d1p4x_len, d1p4x, &p1p4x, 32);
		double p1p4y_p[32], *p1p4y = p1p4y_p;
		int    p1p4y_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, d1p4y_len, d1p4y, &p1p4y, 32);
		double p1p4z_p[32], *p1p4z = p1p4z_p;
		int    p1p4z_len =
		  o.Gen_Diff_With_PreAlloc(l1z_len, l1z, d1p4z_len, d1p4z, &p1p4z, 32);
		double p2p4x_p[32], *p2p4x = p2p4x_p;
		int    p2p4x_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, d2p4x_len, d2p4x, &p2p4x, 32);
		double p2p4y_p[32], *p2p4y = p2p4y_p;
		int    p2p4y_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, d2p4y_len, d2p4y, &p2p4y, 32);
		double p2p4z_p[32], *p2p4z = p2p4z_p;
		int    p2p4z_len =
		  o.Gen_Diff_With_PreAlloc(l2z_len, l2z, d2p4z_len, d2p4z, &p2p4z, 32);
		double p3p4x[2];
		o.two_Diff(p3x, p4x, p3p4x);
		double p3p4y[2];
		o.two_Diff(p3y, p4y, p3p4y);
		double p3p4z[2];
		o.two_Diff(p3z, p4z, p3p4z);
		double tmc_a_p[32], *tmc_a = tmc_a_p;
		int    tmc_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4y_len,
		                                               p2p4y, &tmc_a, 32);
		double tmc_b_p[32], *tmc_b = tmc_b_p;
		int    tmc_b_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4x_len,
		                                               p2p4x, &tmc_b, 32);
		double m01_p[32], *m01 = m01_p;
		int    m01_len =
		  o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 32);
		double tmi_a_p[32], *tmi_a = tmi_a_p;
		int    tmi_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4z_len,
		                                               p2p4z, &tmi_a, 32);
		double tmi_b_p[32], *tmi_b = tmi_b_p;
		int    tmi_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4x_len,
		                                               p2p4x, &tmi_b, 32);
		double m02_p[32], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 32);
		double tma_a_p[32], *tma_a = tma_a_p;
		int    tma_a_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4z_len,
		                                               p2p4z, &tma_a, 32);
		double tma_b_p[32], *tma_b = tma_b_p;
		int    tma_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4y_len,
		                                               p2p4y, &tma_b, 32);
		double m12_p[32], *m12 = m12_p;
		int    m12_len =
		  o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 32);
		double mt1_p[32], *mt1 = mt1_p;
		int mt1_len = o.Gen_Product_With_PreAlloc(m01_len, m01, 2, p3p4z, &mt1, 32);
		double mt2_p[32], *mt2 = mt2_p;
		int mt2_len = o.Gen_Product_With_PreAlloc(m02_len, m02, 2, p3p4y, &mt2, 32);
		double mt3_p[32], *mt3 = mt3_p;
		int mt3_len = o.Gen_Product_With_PreAlloc(m12_len, m12, 2, p3p4x, &mt3, 32);
		double mtt_p[32], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 32);
		double m012_p[32], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 32);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tma_b_p != tma_b)
			FreeDoubles(tma_b);
		if (tma_a_p != tma_a)
			FreeDoubles(tma_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmc_b_p != tmc_b)
			FreeDoubles(tmc_b);
		if (tmc_a_p != tmc_a)
			FreeDoubles(tmc_a);
		if (p2p4z_p != p2p4z)
			FreeDoubles(p2p4z);
		if (p2p4y_p != p2p4y)
			FreeDoubles(p2p4y);
		if (p2p4x_p != p2p4x)
			FreeDoubles(p2p4x);
		if (p1p4z_p != p1p4z)
			FreeDoubles(p1p4z);
		if (p1p4y_p != p1p4y)
			FreeDoubles(p1p4y);
		if (p1p4x_p != p1p4x)
			FreeDoubles(p1p4x);
		if (d2p4z_p != d2p4z)
			FreeDoubles(d2p4z);
		if (d2p4y_p != d2p4y)
			FreeDoubles(d2p4y);
		if (d2p4x_p != d2p4x)
			FreeDoubles(d2p4x);
		if (d1p4z_p != d1p4z)
			FreeDoubles(d1p4z);
		if (d1p4y_p != d1p4y)
			FreeDoubles(d1p4y);
		if (d1p4x_p != d1p4x)
			FreeDoubles(d1p4x);
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
		return orient3D_IIEE_exact<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2, double p3x, double p3y,
                   double p3z, double p4x, double p4y, double p4z, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret =
		  orient3D_IIEE_filtered<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IIEE_interval<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IIEE_expansion<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	return orient3D_IIEE<IT, ET, WithSSFilter>(p1, p2, p3.x(), p3.y(), p3.z(),
	                                           p4.x(), p4.y(), p4.z(), arr);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3, double p4x,
                            double p4y, double p4z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var))
		return Sign::UNCERTAIN;

	double d1p4x = d1 * p4x;
	double d1p4y = d1 * p4y;
	double d1p4z = d1 * p4z;
	double d2p4x = d2 * p4x;
	double d2p4y = d2 * p4y;
	double d2p4z = d2 * p4z;
	double d3p4x = d3 * p4x;
	double d3p4y = d3 * p4y;
	double d3p4z = d3 * p4z;
	double p1p4x = l1x - d1p4x;
	double p1p4y = l1y - d1p4y;
	double p1p4z = l1z - d1p4z;
	double p2p4x = l2x - d2p4x;
	double p2p4y = l2y - d2p4y;
	double p2p4z = l2z - d2p4z;
	double p3p4x = l3x - d3p4x;
	double p3p4y = l3y - d3p4y;
	double p3p4z = l3z - d3p4z;
	double tmc_a = p1p4x * p2p4y;
	double tmc_b = p1p4y * p2p4x;
	double m01   = tmc_a - tmc_b;
	double tmi_a = p1p4x * p2p4z;
	double tmi_b = p1p4z * p2p4x;
	double m02   = tmi_a - tmi_b;
	double tma_a = p1p4y * p2p4z;
	double tma_b = p1p4z * p2p4y;
	double m12   = tma_a - tma_b;
	double mt1   = m01 * p3p4z;
	double mt2   = m02 * p3p4y;
	double mt3   = m12 * p3p4x;
	double mtt   = mt2 - mt1;
	double m012  = mtt - mt3;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p4z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.3892888495756848e-12;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1492353076125757e-11;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.6373565003835534e-10;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.725680701449641e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.254050666758681e-10;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 7.008449820489246e-09;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.2701613979343482e-10;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7060943907631994e-09;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.211968919141341e-08;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.8087548287203607e-07;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3, IT p4x, IT p4y,
                            IT p4z)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT d1p4x = d1 * p4x;
	IT d1p4y = d1 * p4y;
	IT d1p4z = d1 * p4z;
	IT d2p4x = d2 * p4x;
	IT d2p4y = d2 * p4y;
	IT d2p4z = d2 * p4z;
	IT d3p4x = d3 * p4x;
	IT d3p4y = d3 * p4y;
	IT d3p4z = d3 * p4z;
	IT p1p4x = l1x - d1p4x;
	IT p1p4y = l1y - d1p4y;
	IT p1p4z = l1z - d1p4z;
	IT p2p4x = l2x - d2p4x;
	IT p2p4y = l2y - d2p4y;
	IT p2p4z = l2z - d2p4z;
	IT p3p4x = l3x - d3p4x;
	IT p3p4y = l3y - d3p4y;
	IT p3p4z = l3z - d3p4z;
	IT tmc_a = p1p4x * p2p4y;
	IT tmc_b = p1p4y * p2p4x;
	IT m01   = tmc_a - tmc_b;
	IT tmi_a = p1p4x * p2p4z;
	IT tmi_b = p1p4z * p2p4x;
	IT m02   = tmi_a - tmi_b;
	IT tma_a = p1p4y * p2p4z;
	IT tma_b = p1p4z * p2p4y;
	IT m12   = tma_a - tma_b;
	IT mt1   = m01 * p3p4z;
	IT mt2   = m02 * p3p4y;
	IT mt3   = m12 * p3p4x;
	IT mtt   = mt2 - mt1;
	IT m012  = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2,
                         const GenericPoint3T<IT, ET> &p3, ET p4x, ET p4y,
                         ET p4z)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	ET d1p4x = d1 * p4x;
	ET d1p4y = d1 * p4y;
	ET d1p4z = d1 * p4z;
	ET d2p4x = d2 * p4x;
	ET d2p4y = d2 * p4y;
	ET d2p4z = d2 * p4z;
	ET d3p4x = d3 * p4x;
	ET d3p4y = d3 * p4y;
	ET d3p4z = d3 * p4z;
	ET p1p4x = l1x - d1p4x;
	ET p1p4y = l1y - d1p4y;
	ET p1p4z = l1z - d1p4z;
	ET p2p4x = l2x - d2p4x;
	ET p2p4y = l2y - d2p4y;
	ET p2p4z = l2z - d2p4z;
	ET p3p4x = l3x - d3p4x;
	ET p3p4y = l3y - d3p4y;
	ET p3p4z = l3z - d3p4z;
	ET tmc_a = p1p4x * p2p4y;
	ET tmc_b = p1p4y * p2p4x;
	ET m01   = tmc_a - tmc_b;
	ET tmi_a = p1p4x * p2p4z;
	ET tmi_b = p1p4z * p2p4x;
	ET m02   = tmi_a - tmi_b;
	ET tma_a = p1p4y * p2p4z;
	ET tma_b = p1p4z * p2p4y;
	ET m12   = tma_a - tma_b;
	ET mt1   = m01 * p3p4z;
	ET mt2   = m02 * p3p4y;
	ET mt3   = m12 * p3p4x;
	ET mtt   = mt2 - mt1;
	ET m012  = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3, double p4x,
                             double p4y, double p4z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32],
	                  *l1z = l1z_p, d1_p[32], *d1 = d1_p, l2x_p[32], *l2x = l2x_p,
	                  l2y_p[32], *l2y = l2y_p, l2z_p[32], *l2z = l2z_p, d2_p[32],
	                  *d2 = d2_p, l3x_p[32], *l3x = l3x_p, l3y_p[32],
	                  *l3y = l3y_p, l3z_p[32], *l3z = l3z_p, d3_p[32], *d3 = d3_p;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32, l3x_len = 32, l3y_len = 32,
	    l3z_len = 32, d3_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0))
	{
		expansionObject o;
		double          d1p4x_p[32], *d1p4x = d1p4x_p;
		int    d1p4x_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p4x, &d1p4x, 32);
		double d1p4y_p[32], *d1p4y = d1p4y_p;
		int    d1p4y_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p4y, &d1p4y, 32);
		double d1p4z_p[32], *d1p4z = d1p4z_p;
		int    d1p4z_len = o.Gen_Scale_With_PreAlloc(d1_len, d1, p4z, &d1p4z, 32);
		double d2p4x_p[32], *d2p4x = d2p4x_p;
		int    d2p4x_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, p4x, &d2p4x, 32);
		double d2p4y_p[32], *d2p4y = d2p4y_p;
		int    d2p4y_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, p4y, &d2p4y, 32);
		double d2p4z_p[32], *d2p4z = d2p4z_p;
		int    d2p4z_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, p4z, &d2p4z, 32);
		double d3p4x_p[32], *d3p4x = d3p4x_p;
		int    d3p4x_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, p4x, &d3p4x, 32);
		double d3p4y_p[32], *d3p4y = d3p4y_p;
		int    d3p4y_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, p4y, &d3p4y, 32);
		double d3p4z_p[32], *d3p4z = d3p4z_p;
		int    d3p4z_len = o.Gen_Scale_With_PreAlloc(d3_len, d3, p4z, &d3p4z, 32);
		double p1p4x_p[32], *p1p4x = p1p4x_p;
		int    p1p4x_len =
		  o.Gen_Diff_With_PreAlloc(l1x_len, l1x, d1p4x_len, d1p4x, &p1p4x, 32);
		double p1p4y_p[32], *p1p4y = p1p4y_p;
		int    p1p4y_len =
		  o.Gen_Diff_With_PreAlloc(l1y_len, l1y, d1p4y_len, d1p4y, &p1p4y, 32);
		double p1p4z_p[32], *p1p4z = p1p4z_p;
		int    p1p4z_len =
		  o.Gen_Diff_With_PreAlloc(l1z_len, l1z, d1p4z_len, d1p4z, &p1p4z, 32);
		double p2p4x_p[32], *p2p4x = p2p4x_p;
		int    p2p4x_len =
		  o.Gen_Diff_With_PreAlloc(l2x_len, l2x, d2p4x_len, d2p4x, &p2p4x, 32);
		double p2p4y_p[32], *p2p4y = p2p4y_p;
		int    p2p4y_len =
		  o.Gen_Diff_With_PreAlloc(l2y_len, l2y, d2p4y_len, d2p4y, &p2p4y, 32);
		double p2p4z_p[32], *p2p4z = p2p4z_p;
		int    p2p4z_len =
		  o.Gen_Diff_With_PreAlloc(l2z_len, l2z, d2p4z_len, d2p4z, &p2p4z, 32);
		double p3p4x_p[32], *p3p4x = p3p4x_p;
		int    p3p4x_len =
		  o.Gen_Diff_With_PreAlloc(l3x_len, l3x, d3p4x_len, d3p4x, &p3p4x, 32);
		double p3p4y_p[32], *p3p4y = p3p4y_p;
		int    p3p4y_len =
		  o.Gen_Diff_With_PreAlloc(l3y_len, l3y, d3p4y_len, d3p4y, &p3p4y, 32);
		double p3p4z_p[32], *p3p4z = p3p4z_p;
		int    p3p4z_len =
		  o.Gen_Diff_With_PreAlloc(l3z_len, l3z, d3p4z_len, d3p4z, &p3p4z, 32);
		double tmc_a_p[32], *tmc_a = tmc_a_p;
		int    tmc_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4y_len,
		                                               p2p4y, &tmc_a, 32);
		double tmc_b_p[32], *tmc_b = tmc_b_p;
		int    tmc_b_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4x_len,
		                                               p2p4x, &tmc_b, 32);
		double m01_p[32], *m01 = m01_p;
		int    m01_len =
		  o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 32);
		double tmi_a_p[32], *tmi_a = tmi_a_p;
		int    tmi_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4z_len,
		                                               p2p4z, &tmi_a, 32);
		double tmi_b_p[32], *tmi_b = tmi_b_p;
		int    tmi_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4x_len,
		                                               p2p4x, &tmi_b, 32);
		double m02_p[32], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 32);
		double tma_a_p[32], *tma_a = tma_a_p;
		int    tma_a_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4z_len,
		                                               p2p4z, &tma_a, 32);
		double tma_b_p[32], *tma_b = tma_b_p;
		int    tma_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4y_len,
		                                               p2p4y, &tma_b, 32);
		double m12_p[32], *m12 = m12_p;
		int    m12_len =
		  o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 32);
		double mt1_p[32], *mt1 = mt1_p;
		int    mt1_len =
		  o.Gen_Product_With_PreAlloc(m01_len, m01, p3p4z_len, p3p4z, &mt1, 32);
		double mt2_p[32], *mt2 = mt2_p;
		int    mt2_len =
		  o.Gen_Product_With_PreAlloc(m02_len, m02, p3p4y_len, p3p4y, &mt2, 32);
		double mt3_p[32], *mt3 = mt3_p;
		int    mt3_len =
		  o.Gen_Product_With_PreAlloc(m12_len, m12, p3p4x_len, p3p4x, &mt3, 32);
		double mtt_p[32], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 32);
		double m012_p[32], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 32);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tma_b_p != tma_b)
			FreeDoubles(tma_b);
		if (tma_a_p != tma_a)
			FreeDoubles(tma_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmc_b_p != tmc_b)
			FreeDoubles(tmc_b);
		if (tmc_a_p != tmc_a)
			FreeDoubles(tmc_a);
		if (p3p4z_p != p3p4z)
			FreeDoubles(p3p4z);
		if (p3p4y_p != p3p4y)
			FreeDoubles(p3p4y);
		if (p3p4x_p != p3p4x)
			FreeDoubles(p3p4x);
		if (p2p4z_p != p2p4z)
			FreeDoubles(p2p4z);
		if (p2p4y_p != p2p4y)
			FreeDoubles(p2p4y);
		if (p2p4x_p != p2p4x)
			FreeDoubles(p2p4x);
		if (p1p4z_p != p1p4z)
			FreeDoubles(p1p4z);
		if (p1p4y_p != p1p4y)
			FreeDoubles(p1p4y);
		if (p1p4x_p != p1p4x)
			FreeDoubles(p1p4x);
		if (d3p4z_p != d3p4z)
			FreeDoubles(d3p4z);
		if (d3p4y_p != d3p4y)
			FreeDoubles(d3p4y);
		if (d3p4x_p != d3p4x)
			FreeDoubles(d3p4x);
		if (d2p4z_p != d2p4z)
			FreeDoubles(d2p4z);
		if (d2p4y_p != d2p4y)
			FreeDoubles(d2p4y);
		if (d2p4x_p != d2p4x)
			FreeDoubles(d2p4x);
		if (d1p4z_p != d1p4z)
			FreeDoubles(d1p4z);
		if (d1p4y_p != d1p4y)
			FreeDoubles(d1p4y);
		if (d1p4x_p != d1p4x)
			FreeDoubles(d1p4x);
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
		return orient3D_IIIE_exact<IT, ET>(p1, p2, p3, p4x, p4y, p4z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIIE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3, double p4x, double p4y,
                   double p4z, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient3D_IIIE_filtered<IT, ET>(p1, p2, p3, p4x, p4y, p4z, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IIIE_interval<IT, ET>(p1, p2, p3, p4x, p4y, p4z);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IIIE_expansion<IT, ET>(p1, p2, p3, p4x, p4y, p4z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIIE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	return orient3D_IIIE<IT, ET, WithSSFilter>(p1, p2, p3, p4.x(), p4.y(), p4.z(),
	                                           arr);
}

template <typename IT, typename ET>
Sign orient3D_IIII_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3,
                            const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z,
	  d4, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var) ||
	    !p4.getFilteredLambda(l4x, l4y, l4z, d4, max_var))
		return Sign::UNCERTAIN;

	double d1p4x   = d1 * l4x;
	double d1p4y   = d1 * l4y;
	double d1p4z   = d1 * l4z;
	double d2p4x   = d2 * l4x;
	double d2p4y   = d2 * l4y;
	double d2p4z   = d2 * l4z;
	double d3p4x   = d3 * l4x;
	double d3p4y   = d3 * l4y;
	double d3p4z   = d3 * l4z;
	double d4l1x   = d4 * l1x;
	double d4l1y   = d4 * l1y;
	double d4l1z   = d4 * l1z;
	double d4l2x   = d4 * l2x;
	double d4l2y   = d4 * l2y;
	double d4l2z   = d4 * l2z;
	double d4l3x   = d4 * l3x;
	double d4l3y   = d4 * l3y;
	double d4l3z   = d4 * l3z;
	double p1p4x   = d4l1x - d1p4x;
	double p1p4y   = d4l1y - d1p4y;
	double p1p4z   = d4l1z - d1p4z;
	double p2p4x   = d4l2x - d2p4x;
	double p2p4y   = d4l2y - d2p4y;
	double p2p4z   = d4l2z - d2p4z;
	double p3p4x   = d4l3x - d3p4x;
	double p3p4y   = d4l3y - d3p4y;
	double p3p4z   = d4l3z - d3p4z;
	double tmc_a   = p1p4x * p2p4y;
	double tmc_b   = p1p4y * p2p4x;
	double m01     = tmc_a - tmc_b;
	double tmi_a   = p1p4x * p2p4z;
	double tmi_b   = p1p4z * p2p4x;
	double m02     = tmi_a - tmi_b;
	double tma_a   = p1p4y * p2p4z;
	double tma_b   = p1p4z * p2p4y;
	double m12     = tma_a - tma_b;
	double mt1     = m01 * p3p4z;
	double mt2     = m02 * p3p4y;
	double mt3     = m12 * p3p4x;
	double mtt     = mt2 - mt1;
	double m012    = mtt - mt3;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1414158507250298e-10;
	}
	break;
	case PntArr3::SSSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.540747695041109e-09;
	}
	break;
	case PntArr3::SSST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.493879559281876e-06;
	}
	break;
	case PntArr3::SSLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1398581634125613e-08;
	}
	break;
	case PntArr3::SSLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7194465556258534e-05;
	}
	break;
	case PntArr3::SSTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.00018409301731026268;
	}
	break;
	case PntArr3::SLLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.591358010979228e-08;
	}
	break;
	case PntArr3::SLLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.37498169812528e-05;
	}
	break;
	case PntArr3::SLTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 0.0005752701378583036;
	}
	break;
	case PntArr3::STTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.00605601238839889;
	}
	break;
	case PntArr3::LLLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1643036135211634e-07;
	}
	break;
	case PntArr3::LLLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.00016759783762410226;
	}
	break;
	case PntArr3::LLTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 0.0017707331971905868;
	}
	break;
	case PntArr3::LTTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.01883943108077826;
	}
	break;
	case PntArr3::TTTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.1952243033447331;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IIII_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3,
                            const GenericPoint3T<IT, ET> &p4)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3) ||
	    !p4.getIntervalLambda(l4x, l4y, l4z, d4))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT d1p4x = d1 * l4x;
	IT d1p4y = d1 * l4y;
	IT d1p4z = d1 * l4z;
	IT d2p4x = d2 * l4x;
	IT d2p4y = d2 * l4y;
	IT d2p4z = d2 * l4z;
	IT d3p4x = d3 * l4x;
	IT d3p4y = d3 * l4y;
	IT d3p4z = d3 * l4z;
	IT d4l1x = d4 * l1x;
	IT d4l1y = d4 * l1y;
	IT d4l1z = d4 * l1z;
	IT d4l2x = d4 * l2x;
	IT d4l2y = d4 * l2y;
	IT d4l2z = d4 * l2z;
	IT d4l3x = d4 * l3x;
	IT d4l3y = d4 * l3y;
	IT d4l3z = d4 * l3z;
	IT p1p4x = d4l1x - d1p4x;
	IT p1p4y = d4l1y - d1p4y;
	IT p1p4z = d4l1z - d1p4z;
	IT p2p4x = d4l2x - d2p4x;
	IT p2p4y = d4l2y - d2p4y;
	IT p2p4z = d4l2z - d2p4z;
	IT p3p4x = d4l3x - d3p4x;
	IT p3p4y = d4l3y - d3p4y;
	IT p3p4z = d4l3z - d3p4z;
	IT tmc_a = p1p4x * p2p4y;
	IT tmc_b = p1p4y * p2p4x;
	IT m01   = tmc_a - tmc_b;
	IT tmi_a = p1p4x * p2p4z;
	IT tmi_b = p1p4z * p2p4x;
	IT m02   = tmi_a - tmi_b;
	IT tma_a = p1p4y * p2p4z;
	IT tma_b = p1p4z * p2p4y;
	IT m12   = tma_a - tma_b;
	IT mt1   = m01 * p3p4z;
	IT mt2   = m02 * p3p4y;
	IT mt3   = m12 * p3p4x;
	IT mtt   = mt2 - mt1;
	IT m012  = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIII_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2,
                         const GenericPoint3T<IT, ET> &p3,
                         const GenericPoint3T<IT, ET> &p4)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, l4x, l4y, l4z, d4;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	p4.getExactLambda(l4x, l4y, l4z, d4);
	ET d1p4x = d1 * l4x;
	ET d1p4y = d1 * l4y;
	ET d1p4z = d1 * l4z;
	ET d2p4x = d2 * l4x;
	ET d2p4y = d2 * l4y;
	ET d2p4z = d2 * l4z;
	ET d3p4x = d3 * l4x;
	ET d3p4y = d3 * l4y;
	ET d3p4z = d3 * l4z;
	ET d4l1x = d4 * l1x;
	ET d4l1y = d4 * l1y;
	ET d4l1z = d4 * l1z;
	ET d4l2x = d4 * l2x;
	ET d4l2y = d4 * l2y;
	ET d4l2z = d4 * l2z;
	ET d4l3x = d4 * l3x;
	ET d4l3y = d4 * l3y;
	ET d4l3z = d4 * l3z;
	ET p1p4x = d4l1x - d1p4x;
	ET p1p4y = d4l1y - d1p4y;
	ET p1p4z = d4l1z - d1p4z;
	ET p2p4x = d4l2x - d2p4x;
	ET p2p4y = d4l2y - d2p4y;
	ET p2p4z = d4l2z - d2p4z;
	ET p3p4x = d4l3x - d3p4x;
	ET p3p4y = d4l3y - d3p4y;
	ET p3p4z = d4l3z - d3p4z;
	ET tmc_a = p1p4x * p2p4y;
	ET tmc_b = p1p4y * p2p4x;
	ET m01   = tmc_a - tmc_b;
	ET tmi_a = p1p4x * p2p4z;
	ET tmi_b = p1p4z * p2p4x;
	ET m02   = tmi_a - tmi_b;
	ET tma_a = p1p4y * p2p4z;
	ET tma_b = p1p4z * p2p4y;
	ET m12   = tma_a - tma_b;
	ET mt1   = m01 * p3p4z;
	ET mt2   = m02 * p3p4y;
	ET mt3   = m12 * p3p4x;
	ET mtt   = mt2 - mt1;
	ET m012  = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIII_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3,
                             const GenericPoint3T<IT, ET> &p4)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32],
	  *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32], *l1z = l1z_p, d1_p[32],
	  *d1 = d1_p, l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p, l2z_p[32],
	  *l2z = l2z_p, d2_p[32], *d2 = d2_p, l3x_p[32], *l3x = l3x_p, l3y_p[32],
	  *l3y = l3y_p, l3z_p[32], *l3z = l3z_p, d3_p[32], *d3 = d3_p, l4x_p[32],
	  *l4x = l4x_p, l4y_p[32], *l4y = l4y_p, l4z_p[32], *l4z = l4z_p, d4_p[32],
	  *d4       = d4_p;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32, l3x_len = 32, l3y_len = 32,
	    l3z_len = 32, d3_len = 32, l4x_len = 32, l4y_len = 32, l4z_len = 32,
	    d4_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len);
	p4.getExpansionLambda(&l4x, l4x_len, &l4y, l4y_len, &l4z, l4z_len, &d4,
	                      d4_len);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) &&
	    (d4[d4_len - 1] != 0))
	{
		expansionObject o;
		double          d1p4x_p[32], *d1p4x = d1p4x_p;
		int             d1p4x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l4x_len, l4x, &d1p4x, 32);
		double d1p4y_p[32], *d1p4y = d1p4y_p;
		int    d1p4y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l4y_len, l4y, &d1p4y, 32);
		double d1p4z_p[32], *d1p4z = d1p4z_p;
		int    d1p4z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l4z_len, l4z, &d1p4z, 32);
		double d2p4x_p[32], *d2p4x = d2p4x_p;
		int    d2p4x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l4x_len, l4x, &d2p4x, 32);
		double d2p4y_p[32], *d2p4y = d2p4y_p;
		int    d2p4y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l4y_len, l4y, &d2p4y, 32);
		double d2p4z_p[32], *d2p4z = d2p4z_p;
		int    d2p4z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l4z_len, l4z, &d2p4z, 32);
		double d3p4x_p[32], *d3p4x = d3p4x_p;
		int    d3p4x_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, l4x_len, l4x, &d3p4x, 32);
		double d3p4y_p[32], *d3p4y = d3p4y_p;
		int    d3p4y_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, l4y_len, l4y, &d3p4y, 32);
		double d3p4z_p[32], *d3p4z = d3p4z_p;
		int    d3p4z_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, l4z_len, l4z, &d3p4z, 32);
		double d4l1x_p[32], *d4l1x = d4l1x_p;
		int    d4l1x_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l1x_len, l1x, &d4l1x, 32);
		double d4l1y_p[32], *d4l1y = d4l1y_p;
		int    d4l1y_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l1y_len, l1y, &d4l1y, 32);
		double d4l1z_p[32], *d4l1z = d4l1z_p;
		int    d4l1z_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l1z_len, l1z, &d4l1z, 32);
		double d4l2x_p[32], *d4l2x = d4l2x_p;
		int    d4l2x_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l2x_len, l2x, &d4l2x, 32);
		double d4l2y_p[32], *d4l2y = d4l2y_p;
		int    d4l2y_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l2y_len, l2y, &d4l2y, 32);
		double d4l2z_p[32], *d4l2z = d4l2z_p;
		int    d4l2z_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l2z_len, l2z, &d4l2z, 32);
		double d4l3x_p[32], *d4l3x = d4l3x_p;
		int    d4l3x_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l3x_len, l3x, &d4l3x, 32);
		double d4l3y_p[32], *d4l3y = d4l3y_p;
		int    d4l3y_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l3y_len, l3y, &d4l3y, 32);
		double d4l3z_p[32], *d4l3z = d4l3z_p;
		int    d4l3z_len =
		  o.Gen_Product_With_PreAlloc(d4_len, d4, l3z_len, l3z, &d4l3z, 32);
		double p1p4x_p[32], *p1p4x = p1p4x_p;
		int    p1p4x_len =
		  o.Gen_Diff_With_PreAlloc(d4l1x_len, d4l1x, d1p4x_len, d1p4x, &p1p4x, 32);
		double p1p4y_p[32], *p1p4y = p1p4y_p;
		int    p1p4y_len =
		  o.Gen_Diff_With_PreAlloc(d4l1y_len, d4l1y, d1p4y_len, d1p4y, &p1p4y, 32);
		double p1p4z_p[32], *p1p4z = p1p4z_p;
		int    p1p4z_len =
		  o.Gen_Diff_With_PreAlloc(d4l1z_len, d4l1z, d1p4z_len, d1p4z, &p1p4z, 32);
		double p2p4x_p[32], *p2p4x = p2p4x_p;
		int    p2p4x_len =
		  o.Gen_Diff_With_PreAlloc(d4l2x_len, d4l2x, d2p4x_len, d2p4x, &p2p4x, 32);
		double p2p4y_p[32], *p2p4y = p2p4y_p;
		int    p2p4y_len =
		  o.Gen_Diff_With_PreAlloc(d4l2y_len, d4l2y, d2p4y_len, d2p4y, &p2p4y, 32);
		double p2p4z_p[32], *p2p4z = p2p4z_p;
		int    p2p4z_len =
		  o.Gen_Diff_With_PreAlloc(d4l2z_len, d4l2z, d2p4z_len, d2p4z, &p2p4z, 32);
		double p3p4x_p[32], *p3p4x = p3p4x_p;
		int    p3p4x_len =
		  o.Gen_Diff_With_PreAlloc(d4l3x_len, d4l3x, d3p4x_len, d3p4x, &p3p4x, 32);
		double p3p4y_p[32], *p3p4y = p3p4y_p;
		int    p3p4y_len =
		  o.Gen_Diff_With_PreAlloc(d4l3y_len, d4l3y, d3p4y_len, d3p4y, &p3p4y, 32);
		double p3p4z_p[32], *p3p4z = p3p4z_p;
		int    p3p4z_len =
		  o.Gen_Diff_With_PreAlloc(d4l3z_len, d4l3z, d3p4z_len, d3p4z, &p3p4z, 32);
		double tmc_a_p[32], *tmc_a = tmc_a_p;
		int    tmc_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4y_len,
		                                               p2p4y, &tmc_a, 32);
		double tmc_b_p[32], *tmc_b = tmc_b_p;
		int    tmc_b_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4x_len,
		                                               p2p4x, &tmc_b, 32);
		double m01_p[32], *m01 = m01_p;
		int    m01_len =
		  o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 32);
		double tmi_a_p[32], *tmi_a = tmi_a_p;
		int    tmi_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4z_len,
		                                               p2p4z, &tmi_a, 32);
		double tmi_b_p[32], *tmi_b = tmi_b_p;
		int    tmi_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4x_len,
		                                               p2p4x, &tmi_b, 32);
		double m02_p[32], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 32);
		double tma_a_p[32], *tma_a = tma_a_p;
		int    tma_a_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4z_len,
		                                               p2p4z, &tma_a, 32);
		double tma_b_p[32], *tma_b = tma_b_p;
		int    tma_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4y_len,
		                                               p2p4y, &tma_b, 32);
		double m12_p[32], *m12 = m12_p;
		int    m12_len =
		  o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 32);
		double mt1_p[32], *mt1 = mt1_p;
		int    mt1_len =
		  o.Gen_Product_With_PreAlloc(m01_len, m01, p3p4z_len, p3p4z, &mt1, 32);
		double mt2_p[32], *mt2 = mt2_p;
		int    mt2_len =
		  o.Gen_Product_With_PreAlloc(m02_len, m02, p3p4y_len, p3p4y, &mt2, 32);
		double mt3_p[32], *mt3 = mt3_p;
		int    mt3_len =
		  o.Gen_Product_With_PreAlloc(m12_len, m12, p3p4x_len, p3p4x, &mt3, 32);
		double mtt_p[32], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 32);
		double m012_p[32], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 32);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tma_b_p != tma_b)
			FreeDoubles(tma_b);
		if (tma_a_p != tma_a)
			FreeDoubles(tma_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmc_b_p != tmc_b)
			FreeDoubles(tmc_b);
		if (tmc_a_p != tmc_a)
			FreeDoubles(tmc_a);
		if (p3p4z_p != p3p4z)
			FreeDoubles(p3p4z);
		if (p3p4y_p != p3p4y)
			FreeDoubles(p3p4y);
		if (p3p4x_p != p3p4x)
			FreeDoubles(p3p4x);
		if (p2p4z_p != p2p4z)
			FreeDoubles(p2p4z);
		if (p2p4y_p != p2p4y)
			FreeDoubles(p2p4y);
		if (p2p4x_p != p2p4x)
			FreeDoubles(p2p4x);
		if (p1p4z_p != p1p4z)
			FreeDoubles(p1p4z);
		if (p1p4y_p != p1p4y)
			FreeDoubles(p1p4y);
		if (p1p4x_p != p1p4x)
			FreeDoubles(p1p4x);
		if (d4l3z_p != d4l3z)
			FreeDoubles(d4l3z);
		if (d4l3y_p != d4l3y)
			FreeDoubles(d4l3y);
		if (d4l3x_p != d4l3x)
			FreeDoubles(d4l3x);
		if (d4l2z_p != d4l2z)
			FreeDoubles(d4l2z);
		if (d4l2y_p != d4l2y)
			FreeDoubles(d4l2y);
		if (d4l2x_p != d4l2x)
			FreeDoubles(d4l2x);
		if (d4l1z_p != d4l1z)
			FreeDoubles(d4l1z);
		if (d4l1y_p != d4l1y)
			FreeDoubles(d4l1y);
		if (d4l1x_p != d4l1x)
			FreeDoubles(d4l1x);
		if (d3p4z_p != d3p4z)
			FreeDoubles(d3p4z);
		if (d3p4y_p != d3p4y)
			FreeDoubles(d3p4y);
		if (d3p4x_p != d3p4x)
			FreeDoubles(d3p4x);
		if (d2p4z_p != d2p4z)
			FreeDoubles(d2p4z);
		if (d2p4y_p != d2p4y)
			FreeDoubles(d2p4y);
		if (d2p4x_p != d2p4x)
			FreeDoubles(d2p4x);
		if (d1p4z_p != d1p4z)
			FreeDoubles(d1p4z);
		if (d1p4y_p != d1p4y)
			FreeDoubles(d1p4y);
		if (d1p4x_p != d1p4x)
			FreeDoubles(d1p4x);
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

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (l4x_p != l4x)
			FreeDoubles(l4x);
		if (l4y_p != l4y)
			FreeDoubles(l4y);
		if (l4z_p != l4z)
			FreeDoubles(l4z);
		if (d4_p != d4)
			FreeDoubles(d4);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orient3D_IIII_exact<IT, ET>(p1, p2, p3, p4);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIII(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient3D_IIII_filtered<IT, ET>(p1, p2, p3, p4, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IIII_interval<IT, ET>(p1, p2, p3, p4);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IIII_expansion<IT, ET>(p1, p2, p3, p4);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2x,
                               double p2y, double p3x, double p3y, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var))
		return Sign::UNCERTAIN;

	double t1x = p2y - p3y;
	double t1y = p3x - p2x;
	double e2  = l1x * t1x;
	double e3  = l1y * t1y;
	double e   = e2 + e3;
	double pr1 = p2x * p3y;
	double pr2 = p2y * p3x;
	double pr  = pr1 - pr2;
	double dpr = d1 * pr;
	double det = dpr + e;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(p2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1y)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.2434497875801758e-14;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 4.75277369543781e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 9.061883188277186e-13;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2x, IT p2y,
                               IT p3x, IT p3y)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t1x = p2y - p3y;
	IT t1y = p3x - p2x;
	IT e2  = l1x * t1x;
	IT e3  = l1y * t1y;
	IT e   = e2 + e3;
	IT pr1 = p2x * p3y;
	IT pr2 = p2y * p3x;
	IT pr  = pr1 - pr2;
	IT dpr = d1 * pr;
	IT det = dpr + e;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2x, ET p2y,
                            ET p3x, ET p3y)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET t1x = p2y - p3y;
	ET t1y = p3x - p2x;
	ET e2  = l1x * t1x;
	ET e3  = l1y * t1y;
	ET e   = e2 + e3;
	ET pr1 = p2x * p3y;
	ET pr2 = p2y * p3x;
	ET pr  = pr1 - pr2;
	ET dpr = d1 * pr;
	ET det = dpr + e;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2x,
                                double p2y, double p3x, double p3y)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          t1x[2];
		o.two_Diff(p2y, p3y, t1x);
		double t1y[2];
		o.two_Diff(p3x, p2x, t1y);
		double e2_p[128], *e2 = e2_p;
		int    e2_len = o.Gen_Product_With_PreAlloc(l1x_len, l1x, 2, t1x, &e2, 128);
		double e3_p[128], *e3 = e3_p;
		int    e3_len = o.Gen_Product_With_PreAlloc(l1y_len, l1y, 2, t1y, &e3, 128);
		double e_p[128], *e = e_p;
		int    e_len = o.Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
		double pr1[2];
		o.Two_Prod(p2x, p3y, pr1);
		double pr2[2];
		o.Two_Prod(p2y, p3x, pr2);
		double pr[4];
		o.Two_Two_Diff(pr1, pr2, pr);
		double dpr_p[128], *dpr = dpr_p;
		int    dpr_len = o.Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
		double det_p[128], *det = det_p;
		int    det_len = o.Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (dpr_p != dpr)
			FreeDoubles(dpr);
		if (e_p != e)
			FreeDoubles(e);
		if (e3_p != e3)
			FreeDoubles(e3);
		if (e2_p != e2)
			FreeDoubles(e2);
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
		return orientOn2Dxy_IEE_exact<IT, ET>(p1, p2x, p2y, p3x, p3y);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IEE(const GenericPoint3T<IT, ET> &p1, double p2x, double p2y,
                      double p3x, double p3y, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dxy_IEE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dxy_IEE_filtered<IT, ET>(p1, p2x, p2y, p3x, p3y, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dxy_IEE, arr);
	}
	ret = orientOn2Dxy_IEE_interval<IT, ET>(p1, p2x, p2y, p3x, p3y);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dxy_IEE, arr);
	return orientOn2Dxy_IEE_expansion<IT, ET>(p1, p2x, p2y, p3x, p3y);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dxy_IEE<IT, ET, WithSSFilter>(p1, p2.x(), p2.y(), p3.x(),
	                                              p3.y(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double op3x,
                               double op3y, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var))
		return Sign::UNCERTAIN;

	double a    = d1 * l2x;
	double b    = d2 * l1x;
	double c    = d1 * op3y;
	double e    = d1 * l2y;
	double f    = d2 * l1y;
	double g    = d1 * op3x;
	double ab   = a - b;
	double cd   = c - l1y;
	double ef   = e - f;
	double gh   = g - l1x;
	double abcd = ab * cd;
	double efgh = ef * gh;
	double L    = abcd - efgh;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(op3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(op3y)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 4.3698378249246237e-13;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.488743156530251e-12;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0438810366396684e-11;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.6996907353794606e-11;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.1849581172128747e-10;
	}
	break;
	case PntArr3::TT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.3071879457225134e-08;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT op3x,
                               IT op3y)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2x;
	IT b    = d2 * l1x;
	IT c    = d1 * op3y;
	IT e    = d1 * l2y;
	IT f    = d2 * l1y;
	IT g    = d1 * op3x;
	IT ab   = a - b;
	IT cd   = c - l1y;
	IT ef   = e - f;
	IT gh   = g - l1x;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET op3x, ET op3y)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET a    = d1 * l2x;
	ET b    = d2 * l1x;
	ET c    = d1 * op3y;
	ET e    = d1 * l2y;
	ET f    = d2 * l1y;
	ET g    = d1 * op3x;
	ET ab   = a - b;
	ET cd   = c - l1y;
	ET ef   = e - f;
	ET gh   = g - l1x;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double op3x,
                      double op3y, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dxy_IIE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dxy_IIE_filtered<IT, ET>(p1, p2, op3x, op3y, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dxy_IIE, arr);
	}
	ret = orientOn2Dxy_IIE_interval<IT, ET>(p1, p2, op3x, op3y);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dxy_IIE, arr);
	return orientOn2Dxy_IIE_expansion<IT, ET>(p1, p2, op3x, op3y);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &op3, PntArr3 arr)
{
	return orientOn2Dxy_IIE<IT, ET, WithSSFilter>(p1, p2, op3.x(), op3.y(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var))
		return Sign::UNCERTAIN;

	double a       = d1 * l2x;
	double b       = d2 * l1x;
	double c       = d1 * l3y;
	double d       = d3 * l1y;
	double e       = d1 * l2y;
	double f       = d2 * l1y;
	double g       = d1 * l3x;
	double h       = d3 * l1x;
	double ab      = a - b;
	double cd      = c - d;
	double ef      = e - f;
	double gh      = g - h;
	double abcd    = ab * cd;
	double efgh    = ef * gh;
	double L       = abcd - efgh;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.53477230924182e-12;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.8889503556637374e-12;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 6.384476280985069e-11;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.6089185539414123e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0530288580289475e-10;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.5007771409235755e-09;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7563428489353397e-10;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 2.144556754402072e-09;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.5356810429144788e-08;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.103174776697444e-06;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2x;
	IT b    = d2 * l1x;
	IT c    = d1 * l3y;
	IT d    = d3 * l1y;
	IT e    = d1 * l2y;
	IT f    = d2 * l1y;
	IT g    = d1 * l3x;
	IT h    = d3 * l1x;
	IT ab   = a - b;
	IT cd   = c - d;
	IT ef   = e - f;
	IT gh   = g - h;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	ET a    = d1 * l2x;
	ET b    = d2 * l1x;
	ET c    = d1 * l3y;
	ET d    = d3 * l1y;
	ET e    = d1 * l2y;
	ET f    = d2 * l1y;
	ET g    = d1 * l3x;
	ET h    = d3 * l1x;
	ET ab   = a - b;
	ET cd   = c - d;
	ET ef   = e - f;
	ET gh   = g - h;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dxy_III, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dxy_III_filtered<IT, ET>(p1, p2, p3, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dxy_III, arr);
	}
	ret = orientOn2Dxy_III_interval<IT, ET>(p1, p2, p3);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dxy_III, arr);
	return orientOn2Dxy_III_expansion<IT, ET>(p1, p2, p3);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2y,
                               double p2z, double p3y, double p3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var))
		return Sign::UNCERTAIN;

	double t1y = p2z - p3z;
	double t1z = p3y - p2y;
	double e2  = l1y * t1y;
	double e3  = l1z * t1z;
	double e   = e2 + e3;
	double pr1 = p2y * p3z;
	double pr2 = p2z * p3y;
	double pr  = pr1 - pr2;
	double dpr = d1 * pr;
	double det = dpr + e;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(p2y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.2434497875801758e-14;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 4.75277369543781e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 9.061883188277186e-13;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2y, IT p2z,
                               IT p3y, IT p3z)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t1y = p2z - p3z;
	IT t1z = p3y - p2y;
	IT e2  = l1y * t1y;
	IT e3  = l1z * t1z;
	IT e   = e2 + e3;
	IT pr1 = p2y * p3z;
	IT pr2 = p2z * p3y;
	IT pr  = pr1 - pr2;
	IT dpr = d1 * pr;
	IT det = dpr + e;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2y, ET p2z,
                            ET p3y, ET p3z)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET t1y = p2z - p3z;
	ET t1z = p3y - p2y;
	ET e2  = l1y * t1y;
	ET e3  = l1z * t1z;
	ET e   = e2 + e3;
	ET pr1 = p2y * p3z;
	ET pr2 = p2z * p3y;
	ET pr  = pr1 - pr2;
	ET dpr = d1 * pr;
	ET det = dpr + e;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2y,
                                double p2z, double p3y, double p3z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          t1y[2];
		o.two_Diff(p2z, p3z, t1y);
		double t1z[2];
		o.two_Diff(p3y, p2y, t1z);
		double e2_p[128], *e2 = e2_p;
		int    e2_len = o.Gen_Product_With_PreAlloc(l1y_len, l1y, 2, t1y, &e2, 128);
		double e3_p[128], *e3 = e3_p;
		int    e3_len = o.Gen_Product_With_PreAlloc(l1z_len, l1z, 2, t1z, &e3, 128);
		double e_p[128], *e = e_p;
		int    e_len = o.Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
		double pr1[2];
		o.Two_Prod(p2y, p3z, pr1);
		double pr2[2];
		o.Two_Prod(p2z, p3y, pr2);
		double pr[4];
		o.Two_Two_Diff(pr1, pr2, pr);
		double dpr_p[128], *dpr = dpr_p;
		int    dpr_len = o.Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
		double det_p[128], *det = det_p;
		int    det_len = o.Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (dpr_p != dpr)
			FreeDoubles(dpr);
		if (e_p != e)
			FreeDoubles(e);
		if (e3_p != e3)
			FreeDoubles(e3);
		if (e2_p != e2)
			FreeDoubles(e2);
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
		return orientOn2Dyz_IEE_exact<IT, ET>(p1, p2y, p2z, p3y, p3z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IEE(const GenericPoint3T<IT, ET> &p1, double p2y, double p2z,
                      double p3y, double p3z, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dyz_IEE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dyz_IEE_filtered<IT, ET>(p1, p2y, p2z, p3y, p3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dyz_IEE, arr);
	}
	ret = orientOn2Dyz_IEE_interval<IT, ET>(p1, p2y, p2z, p3y, p3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dyz_IEE, arr);
	return orientOn2Dyz_IEE_expansion<IT, ET>(p1, p2y, p2z, p3y, p3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dyz_IEE<IT, ET, WithSSFilter>(p1, p2.y(), p2.z(), p3.y(),
	                                              p3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double op3y,
                               double op3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var))
		return Sign::UNCERTAIN;

	double a    = d1 * l2y;
	double b    = d2 * l1y;
	double c    = d1 * op3z;
	double e    = d1 * l2z;
	double f    = d2 * l1z;
	double g    = d1 * op3y;
	double ab   = a - b;
	double cd   = c - l1z;
	double ef   = e - f;
	double gh   = g - l1y;
	double abcd = ab * cd;
	double efgh = ef * gh;
	double L    = abcd - efgh;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(op3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(op3z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 4.3698378249246237e-13;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.488743156530251e-12;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0438810366396684e-11;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.6996907353794606e-11;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.1849581172128747e-10;
	}
	break;
	case PntArr3::TT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.3071879457225134e-08;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT op3y,
                               IT op3z)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2y;
	IT b    = d2 * l1y;
	IT c    = d1 * op3z;
	IT e    = d1 * l2z;
	IT f    = d2 * l1z;
	IT g    = d1 * op3y;
	IT ab   = a - b;
	IT cd   = c - l1z;
	IT ef   = e - f;
	IT gh   = g - l1y;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET op3y, ET op3z)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET a    = d1 * l2y;
	ET b    = d2 * l1y;
	ET c    = d1 * op3z;
	ET e    = d1 * l2z;
	ET f    = d2 * l1z;
	ET g    = d1 * op3y;
	ET ab   = a - b;
	ET cd   = c - l1z;
	ET ef   = e - f;
	ET gh   = g - l1y;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double op3y,
                      double op3z, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dyz_IIE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dyz_IIE_filtered<IT, ET>(p1, p2, op3y, op3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dyz_IIE, arr);
	}
	ret = orientOn2Dyz_IIE_interval<IT, ET>(p1, p2, op3y, op3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dyz_IIE, arr);
	return orientOn2Dyz_IIE_expansion<IT, ET>(p1, p2, op3y, op3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &op3, PntArr3 arr)
{
	return orientOn2Dyz_IIE<IT, ET, WithSSFilter>(p1, p2, op3.y(), op3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var))
		return Sign::UNCERTAIN;

	double a       = d1 * l2y;
	double b       = d2 * l1y;
	double c       = d1 * l3z;
	double d       = d3 * l1z;
	double e       = d1 * l2z;
	double f       = d2 * l1z;
	double g       = d1 * l3y;
	double h       = d3 * l1y;
	double ab      = a - b;
	double cd      = c - d;
	double ef      = e - f;
	double gh      = g - h;
	double abcd    = ab * cd;
	double efgh    = ef * gh;
	double L       = abcd - efgh;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.53477230924182e-12;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.8889503556637374e-12;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 6.384476280985069e-11;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.6089185539414123e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0530288580289475e-10;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.5007771409235755e-09;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7563428489353397e-10;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 2.144556754402072e-09;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.5356810429144788e-08;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.103174776697444e-06;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2y;
	IT b    = d2 * l1y;
	IT c    = d1 * l3z;
	IT d    = d3 * l1z;
	IT e    = d1 * l2z;
	IT f    = d2 * l1z;
	IT g    = d1 * l3y;
	IT h    = d3 * l1y;
	IT ab   = a - b;
	IT cd   = c - d;
	IT ef   = e - f;
	IT gh   = g - h;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	ET a    = d1 * l2y;
	ET b    = d2 * l1y;
	ET c    = d1 * l3z;
	ET d    = d3 * l1z;
	ET e    = d1 * l2z;
	ET f    = d2 * l1z;
	ET g    = d1 * l3y;
	ET h    = d3 * l1y;
	ET ab   = a - b;
	ET cd   = c - d;
	ET ef   = e - f;
	ET gh   = g - h;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dyz_III, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dyz_III_filtered<IT, ET>(p1, p2, p3, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dyz_III, arr);
	}
	ret = orientOn2Dyz_III_interval<IT, ET>(p1, p2, p3);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dyz_III, arr);
	return orientOn2Dyz_III_expansion<IT, ET>(p1, p2, p3);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2x,
                               double p2z, double p3x, double p3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var))
		return Sign::UNCERTAIN;

	double t1z = p2x - p3x;
	double t1x = p3z - p2z;
	double e2  = l1z * t1z;
	double e3  = l1x * t1x;
	double e   = e2 + e3;
	double pr1 = p2z * p3x;
	double pr2 = p2x * p3z;
	double pr  = pr1 - pr2;
	double dpr = d1 * pr;
	double det = dpr + e;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(p2x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(t1x)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.2434497875801758e-14;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 4.75277369543781e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 9.061883188277186e-13;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2x, IT p2z,
                               IT p3x, IT p3z)
{
	IT l1x, l1y, l1z, d1;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t1z = p2x - p3x;
	IT t1x = p3z - p2z;
	IT e2  = l1z * t1z;
	IT e3  = l1x * t1x;
	IT e   = e2 + e3;
	IT pr1 = p2z * p3x;
	IT pr2 = p2x * p3z;
	IT pr  = pr1 - pr2;
	IT dpr = d1 * pr;
	IT det = dpr + e;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2x, ET p2z,
                            ET p3x, ET p3z)
{
	ET l1x, l1y, l1z, d1;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	ET t1z = p2x - p3x;
	ET t1x = p3z - p2z;
	ET e2  = l1z * t1z;
	ET e3  = l1x * t1x;
	ET e   = e2 + e3;
	ET pr1 = p2z * p3x;
	ET pr2 = p2x * p3z;
	ET pr  = pr1 - pr2;
	ET dpr = d1 * pr;
	ET det = dpr + e;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2x,
                                double p2z, double p3x, double p3z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          t1z[2];
		o.two_Diff(p2x, p3x, t1z);
		double t1x[2];
		o.two_Diff(p3z, p2z, t1x);
		double e2_p[128], *e2 = e2_p;
		int    e2_len = o.Gen_Product_With_PreAlloc(l1z_len, l1z, 2, t1z, &e2, 128);
		double e3_p[128], *e3 = e3_p;
		int    e3_len = o.Gen_Product_With_PreAlloc(l1x_len, l1x, 2, t1x, &e3, 128);
		double e_p[128], *e = e_p;
		int    e_len = o.Gen_Sum_With_PreAlloc(e2_len, e2, e3_len, e3, &e, 128);
		double pr1[2];
		o.Two_Prod(p2z, p3x, pr1);
		double pr2[2];
		o.Two_Prod(p2x, p3z, pr2);
		double pr[4];
		o.Two_Two_Diff(pr1, pr2, pr);
		double dpr_p[128], *dpr = dpr_p;
		int    dpr_len = o.Gen_Product_With_PreAlloc(d1_len, d1, 4, pr, &dpr, 128);
		double det_p[128], *det = det_p;
		int    det_len = o.Gen_Sum_With_PreAlloc(dpr_len, dpr, e_len, e, &det, 128);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (dpr_p != dpr)
			FreeDoubles(dpr);
		if (e_p != e)
			FreeDoubles(e);
		if (e3_p != e3)
			FreeDoubles(e3);
		if (e2_p != e2)
			FreeDoubles(e2);
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
		return orientOn2Dzx_IEE_exact<IT, ET>(p1, p2x, p2z, p3x, p3z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IEE(const GenericPoint3T<IT, ET> &p1, double p2x, double p2z,
                      double p3x, double p3z, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dzx_IEE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dzx_IEE_filtered<IT, ET>(p1, p2x, p2z, p3x, p3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dzx_IEE, arr);
	}
	ret = orientOn2Dzx_IEE_interval<IT, ET>(p1, p2x, p2z, p3x, p3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dzx_IEE, arr);
	return orientOn2Dzx_IEE_expansion<IT, ET>(p1, p2x, p2z, p3x, p3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dzx_IEE<IT, ET, WithSSFilter>(p1, p2.x(), p2.z(), p3.x(),
	                                              p3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double op3x,
                               double op3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var))
		return Sign::UNCERTAIN;

	double a    = d1 * l2z;
	double b    = d2 * l1z;
	double c    = d1 * op3x;
	double e    = d1 * l2x;
	double f    = d2 * l1x;
	double g    = d1 * op3z;
	double ab   = a - b;
	double cd   = c - l1x;
	double ef   = e - f;
	double gh   = g - l1z;
	double abcd = ab * cd;
	double efgh = ef * gh;
	double L    = abcd - efgh;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(op3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(op3z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 4.3698378249246237e-13;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.488743156530251e-12;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0438810366396684e-11;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.6996907353794606e-11;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.1849581172128747e-10;
	}
	break;
	case PntArr3::TT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.3071879457225134e-08;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT op3x,
                               IT op3z)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2z;
	IT b    = d2 * l1z;
	IT c    = d1 * op3x;
	IT e    = d1 * l2x;
	IT f    = d2 * l1x;
	IT g    = d1 * op3z;
	IT ab   = a - b;
	IT cd   = c - l1x;
	IT ef   = e - f;
	IT gh   = g - l1z;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET op3x, ET op3z)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	ET a    = d1 * l2z;
	ET b    = d2 * l1z;
	ET c    = d1 * op3x;
	ET e    = d1 * l2x;
	ET f    = d2 * l1x;
	ET g    = d1 * op3z;
	ET ab   = a - b;
	ET cd   = c - l1x;
	ET ef   = e - f;
	ET gh   = g - l1z;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double op3x,
                      double op3z, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dzx_IIE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dzx_IIE_filtered<IT, ET>(p1, p2, op3x, op3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dzx_IIE, arr);
	}
	ret = orientOn2Dzx_IIE_interval<IT, ET>(p1, p2, op3x, op3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dzx_IIE, arr);
	return orientOn2Dzx_IIE_expansion<IT, ET>(p1, p2, op3x, op3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &op3, PntArr3 arr)
{
	return orientOn2Dzx_IIE<IT, ET, WithSSFilter>(p1, p2, op3.x(), op3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, max_var))
		return Sign::UNCERTAIN;

	double a       = d1 * l2z;
	double b       = d2 * l1z;
	double c       = d1 * l3x;
	double d       = d3 * l1x;
	double e       = d1 * l2x;
	double f       = d2 * l1x;
	double g       = d1 * l3z;
	double h       = d3 * l1z;
	double ab      = a - b;
	double cd      = c - d;
	double ef      = e - f;
	double gh      = g - h;
	double abcd    = ab * cd;
	double efgh    = ef * gh;
	double L       = abcd - efgh;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.53477230924182e-12;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.8889503556637374e-12;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 6.384476280985069e-11;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.6089185539414123e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0530288580289475e-10;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.5007771409235755e-09;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7563428489353397e-10;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 2.144556754402072e-09;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.5356810429144788e-08;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.103174776697444e-06;
	}
	break;
	default:
		return Sign::UNCERTAIN; // OMC_EXIT("Unsopported points
		                        // arrangement.");
	}

	return filter_sign(L, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3)
{
	IT l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT a    = d1 * l2z;
	IT b    = d2 * l1z;
	IT c    = d1 * l3x;
	IT d    = d3 * l1x;
	IT e    = d1 * l2x;
	IT f    = d2 * l1x;
	IT g    = d1 * l3z;
	IT h    = d3 * l1z;
	IT ab   = a - b;
	IT cd   = c - d;
	IT ef   = e - f;
	IT gh   = g - h;
	IT abcd = ab * cd;
	IT efgh = ef * gh;
	IT L    = abcd - efgh;
	if (!L.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(L);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3)
{
	ET l1x, l1y, l1z, d1, l2x, l2y, l2z, d2, l3x, l3y, l3z, d3;
	p1.getExactLambda(l1x, l1y, l1z, d1);
	p2.getExactLambda(l2x, l2y, l2z, d2);
	p3.getExactLambda(l3x, l3y, l3z, d3);
	ET a    = d1 * l2z;
	ET b    = d2 * l1z;
	ET c    = d1 * l3x;
	ET d    = d3 * l1x;
	ET e    = d1 * l2x;
	ET f    = d2 * l1x;
	ET g    = d1 * l3z;
	ET h    = d3 * l1z;
	ET ab   = a - b;
	ET cd   = c - d;
	ET ef   = e - f;
	ET gh   = g - h;
	ET abcd = ab * cd;
	ET efgh = ef * gh;
	ET L    = abcd - efgh;
	return OMC::sign(L);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dzx_III, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dzx_III_filtered<IT, ET>(p1, p2, p3, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dzx_III, arr);
	}
	ret = orientOn2Dzx_III_interval<IT, ET>(p1, p2, p3);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dzx_III, arr);
	return orientOn2Dzx_III_expansion<IT, ET>(p1, p2, p3);
}

template <typename IT, typename ET>
Sign squareDistance2D_IE_interval(const GenericPoint2T<IT, ET> &p, IT qx, IT qy,
                                  IT dis)
{
	IT px, py, d;
	if (!p.getIntervalLambda(px, py, d))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dqx     = d * qx;
	IT dqy     = d * qy;
	IT lx      = px - dqx;
	IT ly      = py - dqy;
	IT lx2     = lx * lx;
	IT ly2     = ly * ly;
	IT sqrnorm = lx2 + ly2;
	IT d2      = d * d;
	IT d2dis   = d2 * dis;
	IT diff    = sqrnorm - d2dis;
	if (!diff.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance2D_IE_exact(const GenericPoint2T<IT, ET> &p, ET qx, ET qy,
                               ET dis)
{
	ET px, py, d;
	p.getExactLambda(px, py, d);
	ET dqx     = d * qx;
	ET dqy     = d * qy;
	ET lx      = px - dqx;
	ET ly      = py - dqy;
	ET lx2     = lx * lx;
	ET ly2     = ly * ly;
	ET sqrnorm = lx2 + ly2;
	ET d2      = d * d;
	ET d2dis   = d2 * dis;
	ET diff    = sqrnorm - d2dis;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance2D_IE_expansion(const GenericPoint2T<IT, ET> &p, double qx,
                                   double qy, double dis)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double px_p[128], *px = px_p, py_p[128], *py = py_p, d_p[128], *d = d_p;
	int    px_len = 128, py_len = 128, d_len = 128;
	p.getExpansionLambda(&px, px_len, &py, py_len, &d, d_len);
	if ((d[d_len - 1] != 0))
	{
		expansionObject o;
		double          dqx_p[128], *dqx = dqx_p;
		int    dqx_len = o.Gen_Scale_With_PreAlloc(d_len, d, qx, &dqx, 128);
		double dqy_p[128], *dqy = dqy_p;
		int    dqy_len = o.Gen_Scale_With_PreAlloc(d_len, d, qy, &dqy, 128);
		double lx_p[128], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(px_len, px, dqx_len, dqx, &lx, 128);
		double ly_p[128], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(py_len, py, dqy_len, dqy, &ly, 128);
		double lx2_p[128], *lx2 = lx2_p;
		int    lx2_len =
		  o.Gen_Product_With_PreAlloc(lx_len, lx, lx_len, lx, &lx2, 128);
		double ly2_p[128], *ly2 = ly2_p;
		int    ly2_len =
		  o.Gen_Product_With_PreAlloc(ly_len, ly, ly_len, ly, &ly2, 128);
		double sqrnorm_p[128], *sqrnorm = sqrnorm_p;
		int    sqrnorm_len =
		  o.Gen_Sum_With_PreAlloc(lx2_len, lx2, ly2_len, ly2, &sqrnorm, 128);
		double d2_p[128], *d2 = d2_p;
		int    d2_len = o.Gen_Product_With_PreAlloc(d_len, d, d_len, d, &d2, 128);
		double d2dis_p[128], *d2dis = d2dis_p;
		int    d2dis_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, dis, &d2dis, 128);
		double diff_p[128], *diff = diff_p;
		int    diff_len = o.Gen_Diff_With_PreAlloc(sqrnorm_len, sqrnorm, d2dis_len,
		                                           d2dis, &diff, 128);

		return_value = diff[diff_len - 1];
		if (diff_p != diff)
			FreeDoubles(diff);
		if (d2dis_p != d2dis)
			FreeDoubles(d2dis);
		if (d2_p != d2)
			FreeDoubles(d2);
		if (sqrnorm_p != sqrnorm)
			FreeDoubles(sqrnorm);
		if (ly2_p != ly2)
			FreeDoubles(ly2);
		if (lx2_p != lx2)
			FreeDoubles(lx2);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (dqy_p != dqy)
			FreeDoubles(dqy);
		if (dqx_p != dqx)
			FreeDoubles(dqx);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (px_p != px)
			FreeDoubles(px);
		if (py_p != py)
			FreeDoubles(py);
		if (d_p != d)
			FreeDoubles(d);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return squareDistance2D_IE_exact<IT, ET>(p, qx, qy, dis);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign squareDistance2D_IE(const GenericPoint2T<IT, ET> &p, double qx, double qy,
                         double dis)
{
	Sign ret;
	ret = squareDistance2D_IE_interval<IT, ET>(p, qx, qy, dis);
	if (is_sign_reliable(ret))
		return ret;
	return squareDistance2D_IE_expansion<IT, ET>(p, qx, qy, dis);
}

template <typename IT, typename ET>
Sign squareDistance2D_IE(const GenericPoint2T<IT, ET> &p,
                         const GenericPoint2T<IT, ET> &q, double dis)
{
	return squareDistance2D_IE<IT, ET>(p, q.x(), q.y(), dis);
}

template <typename IT, typename ET>
Sign squareDistance2D_II_interval(const GenericPoint2T<IT, ET> &p,
                                  const GenericPoint2T<IT, ET> &q, IT dis)
{
	IT px, py, dp, qx, qy, dq;
	if (!p.getIntervalLambda(px, py, dp) || !q.getIntervalLambda(qx, qy, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dqx     = dp * qx;
	IT dqy     = dp * qy;
	IT dpx     = dq * px;
	IT dpy     = dq * py;
	IT lx      = dpx - dqx;
	IT ly      = dpy - dqy;
	IT lx2     = lx * lx;
	IT ly2     = ly * ly;
	IT sqrnorm = lx2 + ly2;
	IT dp2     = dp * dp;
	IT dq2     = dq * dq;
	IT d2      = dp2 * dq2;
	IT d2dis   = d2 * dis;
	IT diff    = sqrnorm - d2dis;
	if (!diff.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance2D_II_exact(const GenericPoint2T<IT, ET> &p,
                               const GenericPoint2T<IT, ET> &q, ET dis)
{
	ET px, py, dp, qx, qy, dq;
	p.getExactLambda(px, py, dp);
	q.getExactLambda(qx, qy, dq);
	ET dqx     = dp * qx;
	ET dqy     = dp * qy;
	ET dpx     = dq * px;
	ET dpy     = dq * py;
	ET lx      = dpx - dqx;
	ET ly      = dpy - dqy;
	ET lx2     = lx * lx;
	ET ly2     = ly * ly;
	ET sqrnorm = lx2 + ly2;
	ET dp2     = dp * dp;
	ET dq2     = dq * dq;
	ET d2      = dp2 * dq2;
	ET d2dis   = d2 * dis;
	ET diff    = sqrnorm - d2dis;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance2D_II_expansion(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &q, double dis)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double px_p[64], *px = px_p, py_p[64], *py = py_p, dp_p[64], *dp = dp_p,
	                 qx_p[64], *qx = qx_p, qy_p[64], *qy = qy_p, dq_p[64],
	                 *dq = dq_p;
	int px_len = 64, py_len = 64, dp_len = 64, qx_len = 64, qy_len = 64,
	    dq_len = 64;
	p.getExpansionLambda(&px, px_len, &py, py_len, &dp, dp_len);
	q.getExpansionLambda(&qx, qx_len, &qy, qy_len, &dq, dq_len);
	if ((dp[dp_len - 1] != 0) && (dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          dqx_p[64], *dqx = dqx_p;
		int dqx_len = o.Gen_Product_With_PreAlloc(dp_len, dp, qx_len, qx, &dqx, 64);
		double dqy_p[64], *dqy = dqy_p;
		int dqy_len = o.Gen_Product_With_PreAlloc(dp_len, dp, qy_len, qy, &dqy, 64);
		double dpx_p[64], *dpx = dpx_p;
		int dpx_len = o.Gen_Product_With_PreAlloc(dq_len, dq, px_len, px, &dpx, 64);
		double dpy_p[64], *dpy = dpy_p;
		int dpy_len = o.Gen_Product_With_PreAlloc(dq_len, dq, py_len, py, &dpy, 64);
		double lx_p[64], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(dpx_len, dpx, dqx_len, dqx, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(dpy_len, dpy, dqy_len, dqy, &ly, 64);
		double lx2_p[64], *lx2 = lx2_p;
		int lx2_len = o.Gen_Product_With_PreAlloc(lx_len, lx, lx_len, lx, &lx2, 64);
		double ly2_p[64], *ly2 = ly2_p;
		int ly2_len = o.Gen_Product_With_PreAlloc(ly_len, ly, ly_len, ly, &ly2, 64);
		double sqrnorm_p[64], *sqrnorm = sqrnorm_p;
		int    sqrnorm_len =
		  o.Gen_Sum_With_PreAlloc(lx2_len, lx2, ly2_len, ly2, &sqrnorm, 64);
		double dp2_p[64], *dp2 = dp2_p;
		int dp2_len = o.Gen_Product_With_PreAlloc(dp_len, dp, dp_len, dp, &dp2, 64);
		double dq2_p[64], *dq2 = dq2_p;
		int dq2_len = o.Gen_Product_With_PreAlloc(dq_len, dq, dq_len, dq, &dq2, 64);
		double d2_p[64], *d2 = d2_p;
		int    d2_len =
		  o.Gen_Product_With_PreAlloc(dp2_len, dp2, dq2_len, dq2, &d2, 64);
		double d2dis_p[64], *d2dis = d2dis_p;
		int    d2dis_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, dis, &d2dis, 64);
		double diff_p[64], *diff = diff_p;
		int    diff_len = o.Gen_Diff_With_PreAlloc(sqrnorm_len, sqrnorm, d2dis_len,
		                                           d2dis, &diff, 64);

		return_value = diff[diff_len - 1];
		if (diff_p != diff)
			FreeDoubles(diff);
		if (d2dis_p != d2dis)
			FreeDoubles(d2dis);
		if (d2_p != d2)
			FreeDoubles(d2);
		if (dq2_p != dq2)
			FreeDoubles(dq2);
		if (dp2_p != dp2)
			FreeDoubles(dp2);
		if (sqrnorm_p != sqrnorm)
			FreeDoubles(sqrnorm);
		if (ly2_p != ly2)
			FreeDoubles(ly2);
		if (lx2_p != lx2)
			FreeDoubles(lx2);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (dpy_p != dpy)
			FreeDoubles(dpy);
		if (dpx_p != dpx)
			FreeDoubles(dpx);
		if (dqy_p != dqy)
			FreeDoubles(dqy);
		if (dqx_p != dqx)
			FreeDoubles(dqx);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (px_p != px)
			FreeDoubles(px);
		if (py_p != py)
			FreeDoubles(py);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint2T<IT, ET>::global_cached_values_enabled())
	{
		if (qx_p != qx)
			FreeDoubles(qx);
		if (qy_p != qy)
			FreeDoubles(qy);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return squareDistance2D_II_exact<IT, ET>(p, q, dis);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign squareDistance2D_II(const GenericPoint2T<IT, ET> &p,
                         const GenericPoint2T<IT, ET> &q, double dis)
{
	Sign ret;
	ret = squareDistance2D_II_interval<IT, ET>(p, q, dis);
	if (is_sign_reliable(ret))
		return ret;
	return squareDistance2D_II_expansion<IT, ET>(p, q, dis);
}

template <typename IT, typename ET>
Sign squareDistance3D_IE_interval(const GenericPoint3T<IT, ET> &p, IT qx, IT qy,
                                  IT qz, IT dis)
{
	IT px, py, pz, d;
	if (!p.getIntervalLambda(px, py, pz, d))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dqx       = d * qx;
	IT dqy       = d * qy;
	IT dqz       = d * qz;
	IT lx        = px - dqx;
	IT ly        = py - dqy;
	IT lz        = pz - dqz;
	IT lx2       = lx * lx;
	IT ly2       = ly * ly;
	IT lz2       = lz * lz;
	IT sqrnormxy = lx2 + ly2;
	IT sqrnorm   = sqrnormxy + lz2;
	IT d2        = d * d;
	IT d2dis     = d2 * dis;
	IT diff      = sqrnorm - d2dis;
	if (!diff.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance3D_IE_exact(const GenericPoint3T<IT, ET> &p, ET qx, ET qy,
                               ET qz, ET dis)
{
	ET px, py, pz, d;
	p.getExactLambda(px, py, pz, d);
	ET dqx       = d * qx;
	ET dqy       = d * qy;
	ET dqz       = d * qz;
	ET lx        = px - dqx;
	ET ly        = py - dqy;
	ET lz        = pz - dqz;
	ET lx2       = lx * lx;
	ET ly2       = ly * ly;
	ET lz2       = lz * lz;
	ET sqrnormxy = lx2 + ly2;
	ET sqrnorm   = sqrnormxy + lz2;
	ET d2        = d * d;
	ET d2dis     = d2 * dis;
	ET diff      = sqrnorm - d2dis;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance3D_IE_expansion(const GenericPoint3T<IT, ET> &p, double qx,
                                   double qy, double qz, double dis)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double px_p[64], *px = px_p, py_p[64], *py = py_p, pz_p[64], *pz = pz_p,
	                 d_p[64], *d = d_p;
	int px_len = 64, py_len = 64, pz_len = 64, d_len = 64;
	p.getExpansionLambda(&px, px_len, &py, py_len, &pz, pz_len, &d, d_len);
	if ((d[d_len - 1] != 0))
	{
		expansionObject o;
		double          dqx_p[64], *dqx = dqx_p;
		int             dqx_len = o.Gen_Scale_With_PreAlloc(d_len, d, qx, &dqx, 64);
		double          dqy_p[64], *dqy = dqy_p;
		int             dqy_len = o.Gen_Scale_With_PreAlloc(d_len, d, qy, &dqy, 64);
		double          dqz_p[64], *dqz = dqz_p;
		int             dqz_len = o.Gen_Scale_With_PreAlloc(d_len, d, qz, &dqz, 64);
		double          lx_p[64], *lx = lx_p;
		int    lx_len = o.Gen_Diff_With_PreAlloc(px_len, px, dqx_len, dqx, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int    ly_len = o.Gen_Diff_With_PreAlloc(py_len, py, dqy_len, dqy, &ly, 64);
		double lz_p[64], *lz = lz_p;
		int    lz_len = o.Gen_Diff_With_PreAlloc(pz_len, pz, dqz_len, dqz, &lz, 64);
		double lx2_p[64], *lx2 = lx2_p;
		int lx2_len = o.Gen_Product_With_PreAlloc(lx_len, lx, lx_len, lx, &lx2, 64);
		double ly2_p[64], *ly2 = ly2_p;
		int ly2_len = o.Gen_Product_With_PreAlloc(ly_len, ly, ly_len, ly, &ly2, 64);
		double lz2_p[64], *lz2 = lz2_p;
		int lz2_len = o.Gen_Product_With_PreAlloc(lz_len, lz, lz_len, lz, &lz2, 64);
		double sqrnormxy_p[64], *sqrnormxy = sqrnormxy_p;
		int    sqrnormxy_len =
		  o.Gen_Sum_With_PreAlloc(lx2_len, lx2, ly2_len, ly2, &sqrnormxy, 64);
		double sqrnorm_p[64], *sqrnorm = sqrnorm_p;
		int sqrnorm_len = o.Gen_Sum_With_PreAlloc(sqrnormxy_len, sqrnormxy, lz2_len,
		                                          lz2, &sqrnorm, 64);
		double d2_p[64], *d2 = d2_p;
		int    d2_len = o.Gen_Product_With_PreAlloc(d_len, d, d_len, d, &d2, 64);
		double d2dis_p[64], *d2dis = d2dis_p;
		int    d2dis_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, dis, &d2dis, 64);
		double diff_p[64], *diff = diff_p;
		int    diff_len = o.Gen_Diff_With_PreAlloc(sqrnorm_len, sqrnorm, d2dis_len,
		                                           d2dis, &diff, 64);

		return_value = diff[diff_len - 1];
		if (diff_p != diff)
			FreeDoubles(diff);
		if (d2dis_p != d2dis)
			FreeDoubles(d2dis);
		if (d2_p != d2)
			FreeDoubles(d2);
		if (sqrnorm_p != sqrnorm)
			FreeDoubles(sqrnorm);
		if (sqrnormxy_p != sqrnormxy)
			FreeDoubles(sqrnormxy);
		if (lz2_p != lz2)
			FreeDoubles(lz2);
		if (ly2_p != ly2)
			FreeDoubles(ly2);
		if (lx2_p != lx2)
			FreeDoubles(lx2);
		if (lz_p != lz)
			FreeDoubles(lz);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (dqz_p != dqz)
			FreeDoubles(dqz);
		if (dqy_p != dqy)
			FreeDoubles(dqy);
		if (dqx_p != dqx)
			FreeDoubles(dqx);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (px_p != px)
			FreeDoubles(px);
		if (py_p != py)
			FreeDoubles(py);
		if (pz_p != pz)
			FreeDoubles(pz);
		if (d_p != d)
			FreeDoubles(d);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return squareDistance3D_IE_exact<IT, ET>(p, qx, qy, qz, dis);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign squareDistance3D_IE(const GenericPoint3T<IT, ET> &p, double qx, double qy,
                         double qz, double dis)
{
	Sign ret;
	ret = squareDistance3D_IE_interval<IT, ET>(p, qx, qy, qz, dis);
	if (is_sign_reliable(ret))
		return ret;
	return squareDistance3D_IE_expansion<IT, ET>(p, qx, qy, qz, dis);
}

template <typename IT, typename ET>
Sign squareDistance3D_IE(const GenericPoint3T<IT, ET> &p,
                         const GenericPoint3T<IT, ET> &q, double dis)
{
	return squareDistance3D_IE<IT, ET>(p, q.x(), q.y(), q.z(), dis);
}

template <typename IT, typename ET>
Sign squareDistance3D_II_interval(const GenericPoint3T<IT, ET> &p,
                                  const GenericPoint3T<IT, ET> &q, IT dis)
{
	IT px, py, pz, dp, qx, qy, qz, dq;
	if (!p.getIntervalLambda(px, py, pz, dp) ||
	    !q.getIntervalLambda(qx, qy, qz, dq))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT dqx       = dp * qx;
	IT dqy       = dp * qy;
	IT dqz       = dp * qz;
	IT dpx       = dq * qx;
	IT dpy       = dq * qy;
	IT dpz       = dq * qz;
	IT lx        = dpx - dqx;
	IT ly        = dpy - dqy;
	IT lz        = dpz - dqz;
	IT lx2       = lx * lx;
	IT ly2       = ly * ly;
	IT lz2       = lz * lz;
	IT sqrnormxy = lx2 + ly2;
	IT sqrnorm   = sqrnormxy + lz2;
	IT dp2       = dp * dp;
	IT dq2       = dq * dq;
	IT d2        = dp2 * dq2;
	IT d2dis     = d2 * dis;
	IT diff      = sqrnorm - d2dis;
	if (!diff.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance3D_II_exact(const GenericPoint3T<IT, ET> &p,
                               const GenericPoint3T<IT, ET> &q, ET dis)
{
	ET px, py, pz, dp, qx, qy, qz, dq;
	p.getExactLambda(px, py, pz, dp);
	q.getExactLambda(qx, qy, qz, dq);
	ET dqx       = dp * qx;
	ET dqy       = dp * qy;
	ET dqz       = dp * qz;
	ET dpx       = dq * qx;
	ET dpy       = dq * qy;
	ET dpz       = dq * qz;
	ET lx        = dpx - dqx;
	ET ly        = dpy - dqy;
	ET lz        = dpz - dqz;
	ET lx2       = lx * lx;
	ET ly2       = ly * ly;
	ET lz2       = lz * lz;
	ET sqrnormxy = lx2 + ly2;
	ET sqrnorm   = sqrnormxy + lz2;
	ET dp2       = dp * dp;
	ET dq2       = dq * dq;
	ET d2        = dp2 * dq2;
	ET d2dis     = d2 * dis;
	ET diff      = sqrnorm - d2dis;
	return OMC::sign(diff);
}

template <typename IT, typename ET>
Sign squareDistance3D_II_expansion(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &q, double dis)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double px_p[64], *px = px_p, py_p[64], *py = py_p, pz_p[64], *pz = pz_p,
	                 dp_p[64], *dp = dp_p, qx_p[64], *qx = qx_p, qy_p[64],
	                 *qy = qy_p, qz_p[64], *qz = qz_p, dq_p[64], *dq = dq_p;
	int px_len = 64, py_len = 64, pz_len = 64, dp_len = 64, qx_len = 64,
	    qy_len = 64, qz_len = 64, dq_len = 64;
	p.getExpansionLambda(&px, px_len, &py, py_len, &pz, pz_len, &dp, dp_len);
	q.getExpansionLambda(&qx, qx_len, &qy, qy_len, &qz, qz_len, &dq, dq_len);
	if ((dp[dp_len - 1] != 0) && (dq[dq_len - 1] != 0))
	{
		expansionObject o;
		double          dqx_p[64], *dqx = dqx_p;
		int dqx_len = o.Gen_Product_With_PreAlloc(dp_len, dp, qx_len, qx, &dqx, 64);
		double dqy_p[64], *dqy = dqy_p;
		int dqy_len = o.Gen_Product_With_PreAlloc(dp_len, dp, qy_len, qy, &dqy, 64);
		double dqz_p[64], *dqz = dqz_p;
		int dqz_len = o.Gen_Product_With_PreAlloc(dp_len, dp, qz_len, qz, &dqz, 64);
		double dpx_p[64], *dpx = dpx_p;
		int dpx_len = o.Gen_Product_With_PreAlloc(dq_len, dq, qx_len, qx, &dpx, 64);
		double dpy_p[64], *dpy = dpy_p;
		int dpy_len = o.Gen_Product_With_PreAlloc(dq_len, dq, qy_len, qy, &dpy, 64);
		double dpz_p[64], *dpz = dpz_p;
		int dpz_len = o.Gen_Product_With_PreAlloc(dq_len, dq, qz_len, qz, &dpz, 64);
		double lx_p[64], *lx = lx_p;
		int lx_len = o.Gen_Diff_With_PreAlloc(dpx_len, dpx, dqx_len, dqx, &lx, 64);
		double ly_p[64], *ly = ly_p;
		int ly_len = o.Gen_Diff_With_PreAlloc(dpy_len, dpy, dqy_len, dqy, &ly, 64);
		double lz_p[64], *lz = lz_p;
		int lz_len = o.Gen_Diff_With_PreAlloc(dpz_len, dpz, dqz_len, dqz, &lz, 64);
		double lx2_p[64], *lx2 = lx2_p;
		int lx2_len = o.Gen_Product_With_PreAlloc(lx_len, lx, lx_len, lx, &lx2, 64);
		double ly2_p[64], *ly2 = ly2_p;
		int ly2_len = o.Gen_Product_With_PreAlloc(ly_len, ly, ly_len, ly, &ly2, 64);
		double lz2_p[64], *lz2 = lz2_p;
		int lz2_len = o.Gen_Product_With_PreAlloc(lz_len, lz, lz_len, lz, &lz2, 64);
		double sqrnormxy_p[64], *sqrnormxy = sqrnormxy_p;
		int    sqrnormxy_len =
		  o.Gen_Sum_With_PreAlloc(lx2_len, lx2, ly2_len, ly2, &sqrnormxy, 64);
		double sqrnorm_p[64], *sqrnorm = sqrnorm_p;
		int sqrnorm_len = o.Gen_Sum_With_PreAlloc(sqrnormxy_len, sqrnormxy, lz2_len,
		                                          lz2, &sqrnorm, 64);
		double dp2_p[64], *dp2 = dp2_p;
		int dp2_len = o.Gen_Product_With_PreAlloc(dp_len, dp, dp_len, dp, &dp2, 64);
		double dq2_p[64], *dq2 = dq2_p;
		int dq2_len = o.Gen_Product_With_PreAlloc(dq_len, dq, dq_len, dq, &dq2, 64);
		double d2_p[64], *d2 = d2_p;
		int    d2_len =
		  o.Gen_Product_With_PreAlloc(dp2_len, dp2, dq2_len, dq2, &d2, 64);
		double d2dis_p[64], *d2dis = d2dis_p;
		int    d2dis_len = o.Gen_Scale_With_PreAlloc(d2_len, d2, dis, &d2dis, 64);
		double diff_p[64], *diff = diff_p;
		int    diff_len = o.Gen_Diff_With_PreAlloc(sqrnorm_len, sqrnorm, d2dis_len,
		                                           d2dis, &diff, 64);

		return_value = diff[diff_len - 1];
		if (diff_p != diff)
			FreeDoubles(diff);
		if (d2dis_p != d2dis)
			FreeDoubles(d2dis);
		if (d2_p != d2)
			FreeDoubles(d2);
		if (dq2_p != dq2)
			FreeDoubles(dq2);
		if (dp2_p != dp2)
			FreeDoubles(dp2);
		if (sqrnorm_p != sqrnorm)
			FreeDoubles(sqrnorm);
		if (sqrnormxy_p != sqrnormxy)
			FreeDoubles(sqrnormxy);
		if (lz2_p != lz2)
			FreeDoubles(lz2);
		if (ly2_p != ly2)
			FreeDoubles(ly2);
		if (lx2_p != lx2)
			FreeDoubles(lx2);
		if (lz_p != lz)
			FreeDoubles(lz);
		if (ly_p != ly)
			FreeDoubles(ly);
		if (lx_p != lx)
			FreeDoubles(lx);
		if (dpz_p != dpz)
			FreeDoubles(dpz);
		if (dpy_p != dpy)
			FreeDoubles(dpy);
		if (dpx_p != dpx)
			FreeDoubles(dpx);
		if (dqz_p != dqz)
			FreeDoubles(dqz);
		if (dqy_p != dqy)
			FreeDoubles(dqy);
		if (dqx_p != dqx)
			FreeDoubles(dqx);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (px_p != px)
			FreeDoubles(px);
		if (py_p != py)
			FreeDoubles(py);
		if (pz_p != pz)
			FreeDoubles(pz);
		if (dp_p != dp)
			FreeDoubles(dp);
	}

	if (!GenericPoint3T<IT, ET>::global_cached_values_enabled())
	{
		if (qx_p != qx)
			FreeDoubles(qx);
		if (qy_p != qy)
			FreeDoubles(qy);
		if (qz_p != qz)
			FreeDoubles(qz);
		if (dq_p != dq)
			FreeDoubles(dq);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return squareDistance3D_II_exact<IT, ET>(p, q, dis);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET>
Sign squareDistance3D_II(const GenericPoint3T<IT, ET> &p,
                         const GenericPoint3T<IT, ET> &q, double dis)
{
	Sign ret;
	ret = squareDistance3D_II_interval<IT, ET>(p, q, dis);
	if (is_sign_reliable(ret))
		return ret;
	return squareDistance3D_II_expansion<IT, ET>(p, q, dis);
}

#elif defined(OFFSET_PREDICATES)

template <typename IT, typename ET>
Sign lessThanOnX_IE_filtered(const GenericPoint3T<IT, ET> &p1, double bx,
                             PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return Sign::UNCERTAIN;

	double t0 = b1x - bx;
	double t1 = t0 * d1;
	double kx = t1 + l1x;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(t0)) > max_var)
		max_var = _tmp_fabs;
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

	return filter_sign(kx, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnX_IE_interval(const GenericPoint3T<IT, ET> &p1, IT bx)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t0 = b1x - bx;
	IT t1 = t0 * d1;
	IT kx = t1 + l1x;
	if (!kx.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_IE_exact(const GenericPoint3T<IT, ET> &p1, ET bx)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	ET t0 = b1x - bx;
	ET t1 = t0 * d1;
	ET kx = t1 + l1x;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_IE_expansion(const GenericPoint3T<IT, ET> &p1, double bx)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          t0[2];
		o.two_Diff(b1x, bx, t0);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		double kx_p[128], *kx = kx_p;
		int    kx_len = o.Gen_Sum_With_PreAlloc(t1_len, t1, l1x_len, l1x, &kx, 128);

		return_value = kx[kx_len - 1];
		if (kx_p != kx)
			FreeDoubles(kx);
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
		return lessThanOnX_IE_exact<IT, ET>(p1, bx);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_IE(const GenericPoint3T<IT, ET> &p1, double bx, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnX_IE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnX_IE_filtered<IT, ET>(p1, bx, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnX_IE, arr);
	}
	ret = lessThanOnX_IE_interval<IT, ET>(p1, bx);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnX_IE, arr);
	return lessThanOnX_IE_expansion<IT, ET>(p1, bx);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	return lessThanOnX_IE<IT, ET, WithSSFilter>(p1, b.x(), arr);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return Sign::UNCERTAIN;

	double t0 = b1x - b2x;
	double t1 = t0 * d1;
	double t2 = t1 * d2;
	double t3 = l1x * d2;
	double t4 = l2x * d1;
	double t5 = t2 + t3;
	double kx = t5 - t4;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(t0)) > max_var)
		max_var = _tmp_fabs;
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

	return filter_sign(kx, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t0 = b1x - b2x;
	IT t1 = t0 * d1;
	IT t2 = t1 * d2;
	IT t3 = l1x * d2;
	IT t4 = l2x * d1;
	IT t5 = t2 + t3;
	IT kx = t5 - t4;
	if (!kx.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	ET t0 = b1x - b2x;
	ET t1 = t0 * d1;
	ET t2 = t1 * d2;
	ET t3 = l1x * d2;
	ET t4 = l2x * d1;
	ET t5 = t2 + t3;
	ET kx = t5 - t4;
	return OMC::sign(kx);
}

template <typename IT, typename ET>
Sign lessThanOnX_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2)
{
	double return_value = NAN;
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
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
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
		double kx_p[128], *kx = kx_p;
		int    kx_len = o.Gen_Diff_With_PreAlloc(t5_len, t5, t4_len, t4, &kx, 128);

		return_value = kx[kx_len - 1];
		if (kx_p != kx)
			FreeDoubles(kx);
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
		return lessThanOnX_II_exact<IT, ET>(p1, p2);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnX_II, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnX_II_filtered<IT, ET>(p1, p2, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnX_II, arr);
	}
	ret = lessThanOnX_II_interval<IT, ET>(p1, p2);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnX_II, arr);
	return lessThanOnX_II_expansion<IT, ET>(p1, p2);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_filtered(const GenericPoint3T<IT, ET> &p1, double by,
                             PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return Sign::UNCERTAIN;

	double t0 = b1y - by;
	double t1 = t0 * d1;
	double ky = t1 + l1y;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(t0)) > max_var)
		max_var = _tmp_fabs;
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

	return filter_sign(ky, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_interval(const GenericPoint3T<IT, ET> &p1, IT by)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t0 = b1y - by;
	IT t1 = t0 * d1;
	IT ky = t1 + l1y;
	if (!ky.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_exact(const GenericPoint3T<IT, ET> &p1, ET by)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	ET t0 = b1y - by;
	ET t1 = t0 * d1;
	ET ky = t1 + l1y;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_IE_expansion(const GenericPoint3T<IT, ET> &p1, double by)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          t0[2];
		o.two_Diff(b1y, by, t0);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		double ky_p[128], *ky = ky_p;
		int    ky_len = o.Gen_Sum_With_PreAlloc(t1_len, t1, l1y_len, l1y, &ky, 128);

		return_value = ky[ky_len - 1];
		if (ky_p != ky)
			FreeDoubles(ky);
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
		return lessThanOnY_IE_exact<IT, ET>(p1, by);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_IE(const GenericPoint3T<IT, ET> &p1, double by, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnY_IE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnY_IE_filtered<IT, ET>(p1, by, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnY_IE, arr);
	}
	ret = lessThanOnY_IE_interval<IT, ET>(p1, by);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnY_IE, arr);
	return lessThanOnY_IE_expansion<IT, ET>(p1, by);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	return lessThanOnY_IE<IT, ET, WithSSFilter>(p1, b.y(), arr);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return Sign::UNCERTAIN;

	double t0 = b1y - b2y;
	double t1 = t0 * d1;
	double t2 = t1 * d2;
	double t3 = l1y * d2;
	double t4 = l2y * d1;
	double t5 = t2 + t3;
	double ky = t5 - t4;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(t0)) > max_var)
		max_var = _tmp_fabs;
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

	return filter_sign(ky, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t0 = b1y - b2y;
	IT t1 = t0 * d1;
	IT t2 = t1 * d2;
	IT t3 = l1y * d2;
	IT t4 = l2y * d1;
	IT t5 = t2 + t3;
	IT ky = t5 - t4;
	if (!ky.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	ET t0 = b1y - b2y;
	ET t1 = t0 * d1;
	ET t2 = t1 * d2;
	ET t3 = l1y * d2;
	ET t4 = l2y * d1;
	ET t5 = t2 + t3;
	ET ky = t5 - t4;
	return OMC::sign(ky);
}

template <typename IT, typename ET>
Sign lessThanOnY_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2)
{
	double return_value = NAN;
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
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          t0[2];
		o.two_Diff(b1y, b2y, t0);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		double t2_p[128], *t2 = t2_p;
		int t2_len = o.Gen_Product_With_PreAlloc(t1_len, t1, d2_len, d2, &t2, 128);
		double t3_p[128], *t3 = t3_p;
		int    t3_len =
		  o.Gen_Product_With_PreAlloc(l1y_len, l1y, d2_len, d2, &t3, 128);
		double t4_p[128], *t4 = t4_p;
		int    t4_len =
		  o.Gen_Product_With_PreAlloc(l2y_len, l2y, d1_len, d1, &t4, 128);
		double t5_p[128], *t5 = t5_p;
		int    t5_len = o.Gen_Sum_With_PreAlloc(t2_len, t2, t3_len, t3, &t5, 128);
		double ky_p[128], *ky = ky_p;
		int    ky_len = o.Gen_Diff_With_PreAlloc(t5_len, t5, t4_len, t4, &ky, 128);

		return_value = ky[ky_len - 1];
		if (ky_p != ky)
			FreeDoubles(ky);
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
		return lessThanOnY_II_exact<IT, ET>(p1, p2);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnY_II, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnY_II_filtered<IT, ET>(p1, p2, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnY_II, arr);
	}
	ret = lessThanOnY_II_interval<IT, ET>(p1, p2);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnY_II, arr);
	return lessThanOnY_II_expansion<IT, ET>(p1, p2);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_filtered(const GenericPoint3T<IT, ET> &p1, double bz,
                             PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return Sign::UNCERTAIN;

	double t0 = b1z - bz;
	double t1 = t0 * d1;
	double kz = t1 + l1z;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(t0)) > max_var)
		max_var = _tmp_fabs;
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

	return filter_sign(kz, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_interval(const GenericPoint3T<IT, ET> &p1, IT bz)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t0 = b1z - bz;
	IT t1 = t0 * d1;
	IT kz = t1 + l1z;
	if (!kz.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_exact(const GenericPoint3T<IT, ET> &p1, ET bz)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	ET t0 = b1z - bz;
	ET t1 = t0 * d1;
	ET kz = t1 + l1z;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_IE_expansion(const GenericPoint3T<IT, ET> &p1, double bz)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          t0[2];
		o.two_Diff(b1z, bz, t0);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		double kz_p[128], *kz = kz_p;
		int    kz_len = o.Gen_Sum_With_PreAlloc(t1_len, t1, l1z_len, l1z, &kz, 128);

		return_value = kz[kz_len - 1];
		if (kz_p != kz)
			FreeDoubles(kz);
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
		return lessThanOnZ_IE_exact<IT, ET>(p1, bz);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_IE(const GenericPoint3T<IT, ET> &p1, double bz, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnZ_IE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnZ_IE_filtered<IT, ET>(p1, bz, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnZ_IE, arr);
	}
	ret = lessThanOnZ_IE_interval<IT, ET>(p1, bz);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnZ_IE, arr);
	return lessThanOnZ_IE_expansion<IT, ET>(p1, bz);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr)
{
	return lessThanOnZ_IE<IT, ET, WithSSFilter>(p1, b.z(), arr);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return Sign::UNCERTAIN;

	double t0 = b1z - b2z;
	double t1 = t0 * d1;
	double t2 = t1 * d2;
	double t3 = l1z * d2;
	double t4 = l2z * d1;
	double t5 = t2 + t3;
	double kz = t5 - t4;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(t0)) > max_var)
		max_var = _tmp_fabs;
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

	return filter_sign(kz, epsilon);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT t0 = b1z - b2z;
	IT t1 = t0 * d1;
	IT t2 = t1 * d2;
	IT t3 = l1z * d2;
	IT t4 = l2z * d1;
	IT t5 = t2 + t3;
	IT kz = t5 - t4;
	if (!kz.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	ET t0 = b1z - b2z;
	ET t1 = t0 * d1;
	ET t2 = t1 * d2;
	ET t3 = l1z * d2;
	ET t4 = l2z * d1;
	ET t5 = t2 + t3;
	ET kz = t5 - t4;
	return OMC::sign(kz);
}

template <typename IT, typename ET>
Sign lessThanOnZ_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2)
{
	double return_value = NAN;
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
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          t0[2];
		o.two_Diff(b1z, b2z, t0);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(2, t0, d1_len, d1, &t1, 128);
		double t2_p[128], *t2 = t2_p;
		int t2_len = o.Gen_Product_With_PreAlloc(t1_len, t1, d2_len, d2, &t2, 128);
		double t3_p[128], *t3 = t3_p;
		int    t3_len =
		  o.Gen_Product_With_PreAlloc(l1z_len, l1z, d2_len, d2, &t3, 128);
		double t4_p[128], *t4 = t4_p;
		int    t4_len =
		  o.Gen_Product_With_PreAlloc(l2z_len, l2z, d1_len, d1, &t4, 128);
		double t5_p[128], *t5 = t5_p;
		int    t5_len = o.Gen_Sum_With_PreAlloc(t2_len, t2, t3_len, t3, &t5, 128);
		double kz_p[128], *kz = kz_p;
		int    kz_len = o.Gen_Diff_With_PreAlloc(t5_len, t5, t4_len, t4, &kz, 128);

		return_value = kz[kz_len - 1];
		if (kz_p != kz)
			FreeDoubles(kz);
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
		return lessThanOnZ_II_exact<IT, ET>(p1, p2);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_lessThanOnZ_II, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = lessThanOnZ_II_filtered<IT, ET>(p1, p2, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_lessThanOnZ_II, arr);
	}
	ret = lessThanOnZ_II_interval<IT, ET>(p1, p2);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_lessThanOnZ_II, arr);
	return lessThanOnZ_II_expansion<IT, ET>(p1, p2);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2x,
                            double p2y, double p2z, double p3x, double p3y,
                            double p3z, double p4x, double p4y, double p4z,
                            PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return Sign::UNCERTAIN;

	double b1p4x    = b1x - p4x;
	double b1p4y    = b1y - p4y;
	double b1p4z    = b1z - p4z;
	double d1_b1p4x = d1 * b1p4x;
	double d1_b1p4y = d1 * b1p4y;
	double d1_b1p4z = d1 * b1p4z;
	double p1x_p4x  = d1_b1p4x + l1x;
	double p1y_p4y  = d1_b1p4y + l1y;
	double p1z_p4z  = d1_b1p4z + l1z;
	double p2x_p4x  = p2x - p4x;
	double p2y_p4y  = p2y - p4y;
	double p2z_p4z  = p2z - p4z;
	double tmp4_a   = p1x_p4x * p2y_p4y;
	double tmp4_b   = p1y_p4y * p2x_p4x;
	double m01      = tmp4_a - tmp4_b;
	double tmi_a    = p1x_p4x * p2z_p4z;
	double tmi_b    = p1z_p4z * p2x_p4x;
	double m02      = tmi_a - tmi_b;
	double tmp2_a   = p1y_p4y * p2z_p4z;
	double tmp2_b   = p1z_p4z * p2y_p4y;
	double m12      = tmp2_a - tmp2_b;
	double p3x_p4x  = p3x - p4x;
	double p3y_p4y  = p3y - p4y;
	double p3z_p4z  = p3z - p4z;
	double mt1      = m01 * p3z_p4z;
	double mt2      = m02 * p3y_p4y;
	double mt3      = m12 * p3x_p4x;
	double mtt      = mt2 - mt1;
	double m012     = mtt - mt3;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2x_p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2y_p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2z_p4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3x_p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3y_p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3z_p4z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.552713678800503e-14;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.199125434364624e-13;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.5225156125307035e-12;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2x, IT p2y,
                            IT p2z, IT p3x, IT p3y, IT p3z, IT p4x, IT p4y,
                            IT p4z)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p4x    = b1x - p4x;
	IT b1p4y    = b1y - p4y;
	IT b1p4z    = b1z - p4z;
	IT d1_b1p4x = d1 * b1p4x;
	IT d1_b1p4y = d1 * b1p4y;
	IT d1_b1p4z = d1 * b1p4z;
	IT p1x_p4x  = d1_b1p4x + l1x;
	IT p1y_p4y  = d1_b1p4y + l1y;
	IT p1z_p4z  = d1_b1p4z + l1z;
	IT p2x_p4x  = p2x - p4x;
	IT p2y_p4y  = p2y - p4y;
	IT p2z_p4z  = p2z - p4z;
	IT tmp4_a   = p1x_p4x * p2y_p4y;
	IT tmp4_b   = p1y_p4y * p2x_p4x;
	IT m01      = tmp4_a - tmp4_b;
	IT tmi_a    = p1x_p4x * p2z_p4z;
	IT tmi_b    = p1z_p4z * p2x_p4x;
	IT m02      = tmi_a - tmi_b;
	IT tmp2_a   = p1y_p4y * p2z_p4z;
	IT tmp2_b   = p1z_p4z * p2y_p4y;
	IT m12      = tmp2_a - tmp2_b;
	IT p3x_p4x  = p3x - p4x;
	IT p3y_p4y  = p3y - p4y;
	IT p3z_p4z  = p3z - p4z;
	IT mt1      = m01 * p3z_p4z;
	IT mt2      = m02 * p3y_p4y;
	IT mt3      = m12 * p3x_p4x;
	IT mtt      = mt2 - mt1;
	IT m012     = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2x, ET p2y,
                         ET p2z, ET p3x, ET p3y, ET p3z, ET p4x, ET p4y, ET p4z)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	ET b1p4x    = b1x - p4x;
	ET b1p4y    = b1y - p4y;
	ET b1p4z    = b1z - p4z;
	ET d1_b1p4x = d1 * b1p4x;
	ET d1_b1p4y = d1 * b1p4y;
	ET d1_b1p4z = d1 * b1p4z;
	ET p1x_p4x  = d1_b1p4x + l1x;
	ET p1y_p4y  = d1_b1p4y + l1y;
	ET p1z_p4z  = d1_b1p4z + l1z;
	ET p2x_p4x  = p2x - p4x;
	ET p2y_p4y  = p2y - p4y;
	ET p2z_p4z  = p2z - p4z;
	ET tmp4_a   = p1x_p4x * p2y_p4y;
	ET tmp4_b   = p1y_p4y * p2x_p4x;
	ET m01      = tmp4_a - tmp4_b;
	ET tmi_a    = p1x_p4x * p2z_p4z;
	ET tmi_b    = p1z_p4z * p2x_p4x;
	ET m02      = tmi_a - tmi_b;
	ET tmp2_a   = p1y_p4y * p2z_p4z;
	ET tmp2_b   = p1z_p4z * p2y_p4y;
	ET m12      = tmp2_a - tmp2_b;
	ET p3x_p4x  = p3x - p4x;
	ET p3y_p4y  = p3y - p4y;
	ET p3z_p4z  = p3z - p4z;
	ET mt1      = m01 * p3z_p4z;
	ET mt2      = m02 * p3y_p4y;
	ET mt3      = m12 * p3x_p4x;
	ET mtt      = mt2 - mt1;
	ET m012     = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IEEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2x,
                             double p2y, double p2z, double p3x, double p3y,
                             double p3z, double p4x, double p4y, double p4z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[64], *l1x = l1x_p, l1y_p[64], *l1y = l1y_p, l1z_p[64],
	                  *l1z = l1z_p, d1_p[64], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 64, l1y_len = 64, l1z_len = 64, d1_len = 64;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          b1p4x[2];
		o.two_Diff(b1x, p4x, b1p4x);
		double b1p4y[2];
		o.two_Diff(b1y, p4y, b1p4y);
		double b1p4z[2];
		o.two_Diff(b1z, p4z, b1p4z);
		double d1_b1p4x_p[64], *d1_b1p4x = d1_b1p4x_p;
		int    d1_b1p4x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4x, &d1_b1p4x, 64);
		double d1_b1p4y_p[64], *d1_b1p4y = d1_b1p4y_p;
		int    d1_b1p4y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4y, &d1_b1p4y, 64);
		double d1_b1p4z_p[64], *d1_b1p4z = d1_b1p4z_p;
		int    d1_b1p4z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4z, &d1_b1p4z, 64);
		double p1x_p4x_p[64], *p1x_p4x = p1x_p4x_p;
		int p1x_p4x_len = o.Gen_Sum_With_PreAlloc(d1_b1p4x_len, d1_b1p4x, l1x_len,
		                                          l1x, &p1x_p4x, 64);
		double p1y_p4y_p[64], *p1y_p4y = p1y_p4y_p;
		int p1y_p4y_len = o.Gen_Sum_With_PreAlloc(d1_b1p4y_len, d1_b1p4y, l1y_len,
		                                          l1y, &p1y_p4y, 64);
		double p1z_p4z_p[64], *p1z_p4z = p1z_p4z_p;
		int p1z_p4z_len = o.Gen_Sum_With_PreAlloc(d1_b1p4z_len, d1_b1p4z, l1z_len,
		                                          l1z, &p1z_p4z, 64);
		double p2x_p4x[2];
		o.two_Diff(p2x, p4x, p2x_p4x);
		double p2y_p4y[2];
		o.two_Diff(p2y, p4y, p2y_p4y);
		double p2z_p4z[2];
		o.two_Diff(p2z, p4z, p2z_p4z);
		double tmp4_a_p[64], *tmp4_a = tmp4_a_p;
		int    tmp4_a_len = o.Gen_Product_With_PreAlloc(p1x_p4x_len, p1x_p4x, 2,
		                                                p2y_p4y, &tmp4_a, 64);
		double tmp4_b_p[64], *tmp4_b = tmp4_b_p;
		int    tmp4_b_len = o.Gen_Product_With_PreAlloc(p1y_p4y_len, p1y_p4y, 2,
		                                                p2x_p4x, &tmp4_b, 64);
		double m01_p[64], *m01 = m01_p;
		int    m01_len = o.Gen_Diff_With_PreAlloc(tmp4_a_len, tmp4_a, tmp4_b_len,
		                                          tmp4_b, &m01, 64);
		double tmi_a_p[64], *tmi_a = tmi_a_p;
		int    tmi_a_len =
		  o.Gen_Product_With_PreAlloc(p1x_p4x_len, p1x_p4x, 2, p2z_p4z, &tmi_a, 64);
		double tmi_b_p[64], *tmi_b = tmi_b_p;
		int    tmi_b_len =
		  o.Gen_Product_With_PreAlloc(p1z_p4z_len, p1z_p4z, 2, p2x_p4x, &tmi_b, 64);
		double m02_p[64], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 64);
		double tmp2_a_p[64], *tmp2_a = tmp2_a_p;
		int    tmp2_a_len = o.Gen_Product_With_PreAlloc(p1y_p4y_len, p1y_p4y, 2,
		                                                p2z_p4z, &tmp2_a, 64);
		double tmp2_b_p[64], *tmp2_b = tmp2_b_p;
		int    tmp2_b_len = o.Gen_Product_With_PreAlloc(p1z_p4z_len, p1z_p4z, 2,
		                                                p2y_p4y, &tmp2_b, 64);
		double m12_p[64], *m12 = m12_p;
		int    m12_len = o.Gen_Diff_With_PreAlloc(tmp2_a_len, tmp2_a, tmp2_b_len,
		                                          tmp2_b, &m12, 64);
		double p3x_p4x[2];
		o.two_Diff(p3x, p4x, p3x_p4x);
		double p3y_p4y[2];
		o.two_Diff(p3y, p4y, p3y_p4y);
		double p3z_p4z[2];
		o.two_Diff(p3z, p4z, p3z_p4z);
		double mt1_p[64], *mt1 = mt1_p;
		int    mt1_len =
		  o.Gen_Product_With_PreAlloc(m01_len, m01, 2, p3z_p4z, &mt1, 64);
		double mt2_p[64], *mt2 = mt2_p;
		int    mt2_len =
		  o.Gen_Product_With_PreAlloc(m02_len, m02, 2, p3y_p4y, &mt2, 64);
		double mt3_p[64], *mt3 = mt3_p;
		int    mt3_len =
		  o.Gen_Product_With_PreAlloc(m12_len, m12, 2, p3x_p4x, &mt3, 64);
		double mtt_p[64], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 64);
		double m012_p[64], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 64);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tmp2_b_p != tmp2_b)
			FreeDoubles(tmp2_b);
		if (tmp2_a_p != tmp2_a)
			FreeDoubles(tmp2_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmp4_b_p != tmp4_b)
			FreeDoubles(tmp4_b);
		if (tmp4_a_p != tmp4_a)
			FreeDoubles(tmp4_a);
		if (p1z_p4z_p != p1z_p4z)
			FreeDoubles(p1z_p4z);
		if (p1y_p4y_p != p1y_p4y)
			FreeDoubles(p1y_p4y);
		if (p1x_p4x_p != p1x_p4x)
			FreeDoubles(p1x_p4x);
		if (d1_b1p4z_p != d1_b1p4z)
			FreeDoubles(d1_b1p4z);
		if (d1_b1p4y_p != d1_b1p4y)
			FreeDoubles(d1_b1p4y);
		if (d1_b1p4x_p != d1_b1p4x)
			FreeDoubles(d1_b1p4x);
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
		return orient3D_IEEE_exact<IT, ET>(p1, p2x, p2y, p2z, p3x, p3y, p3z, p4x,
		                                   p4y, p4z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IEEE(const GenericPoint3T<IT, ET> &p1, double p2x, double p2y,
                   double p2z, double p3x, double p3y, double p3z, double p4x,
                   double p4y, double p4z, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient3D_IEEE_filtered<IT, ET>(p1, p2x, p2y, p2z, p3x, p3y, p3z, p4x,
		                                     p4y, p4z, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IEEE_interval<IT, ET>(p1, p2x, p2y, p2z, p3x, p3y, p3z, p4x,
	                                     p4y, p4z);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IEEE_expansion<IT, ET>(p1, p2x, p2y, p2z, p3x, p3y, p3z, p4x,
	                                       p4y, p4z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IEEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	return orient3D_IEEE<IT, ET, WithSSFilter>(p1, p2.x(), p2.y(), p2.z(), p3.x(),
	                                           p3.y(), p3.z(), p4.x(), p4.y(),
	                                           p4.z(), arr);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, double p3x,
                            double p3y, double p3z, double p4x, double p4y,
                            double p4z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return Sign::UNCERTAIN;

	double b1p4x    = b1x - p4x;
	double b1p4y    = b1y - p4y;
	double b1p4z    = b1z - p4z;
	double b2p4x    = b2x - p4x;
	double b2p4y    = b2y - p4y;
	double b2p4z    = b2z - p4z;
	double d1_b1p4x = d1 * b1p4x;
	double d1_b1p4y = d1 * b1p4y;
	double d1_b1p4z = d1 * b1p4z;
	double d2_b2p4x = d2 * b2p4x;
	double d2_b2p4y = d2 * b2p4y;
	double d2_b2p4z = d2 * b2p4z;
	double p1p4x    = l1x + d1_b1p4x;
	double p1p4y    = l1y + d1_b1p4y;
	double p1p4z    = l1z + d1_b1p4z;
	double p2p4x    = l2x + d2_b2p4x;
	double p2p4y    = l2y + d2_b2p4y;
	double p2p4z    = l2z + d2_b2p4z;
	double p3p4x    = p3x - p4x;
	double p3p4y    = p3y - p4y;
	double p3p4z    = p3z - p4z;
	double tmc_a    = p1p4x * p2p4y;
	double tmc_b    = p1p4y * p2p4x;
	double m01      = tmc_a - tmc_b;
	double tmi_a    = p1p4x * p2p4z;
	double tmi_b    = p1p4z * p2p4x;
	double m02      = tmi_a - tmi_b;
	double tma_a    = p1p4y * p2p4z;
	double tma_b    = p1p4z * p2p4y;
	double m12      = tma_a - tma_b;
	double mt1      = m01 * p3p4z;
	double mt2      = m02 * p3p4y;
	double mt3      = m12 * p3p4x;
	double mtt      = mt2 - mt1;
	double m012     = mtt - mt3;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p3p4z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0605739337042936e-13;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 6.714967124010773e-13;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7543300145916922e-11;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 2.3201249949034204e-12;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.8031851893347754e-11;
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
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.3534560139305596e-09;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, IT p3x, IT p3y,
                            IT p3z, IT p4x, IT p4y, IT p4z)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p4x    = b1x - p4x;
	IT b1p4y    = b1y - p4y;
	IT b1p4z    = b1z - p4z;
	IT b2p4x    = b2x - p4x;
	IT b2p4y    = b2y - p4y;
	IT b2p4z    = b2z - p4z;
	IT d1_b1p4x = d1 * b1p4x;
	IT d1_b1p4y = d1 * b1p4y;
	IT d1_b1p4z = d1 * b1p4z;
	IT d2_b2p4x = d2 * b2p4x;
	IT d2_b2p4y = d2 * b2p4y;
	IT d2_b2p4z = d2 * b2p4z;
	IT p1p4x    = l1x + d1_b1p4x;
	IT p1p4y    = l1y + d1_b1p4y;
	IT p1p4z    = l1z + d1_b1p4z;
	IT p2p4x    = l2x + d2_b2p4x;
	IT p2p4y    = l2y + d2_b2p4y;
	IT p2p4z    = l2z + d2_b2p4z;
	IT p3p4x    = p3x - p4x;
	IT p3p4y    = p3y - p4y;
	IT p3p4z    = p3z - p4z;
	IT tmc_a    = p1p4x * p2p4y;
	IT tmc_b    = p1p4y * p2p4x;
	IT m01      = tmc_a - tmc_b;
	IT tmi_a    = p1p4x * p2p4z;
	IT tmi_b    = p1p4z * p2p4x;
	IT m02      = tmi_a - tmi_b;
	IT tma_a    = p1p4y * p2p4z;
	IT tma_b    = p1p4z * p2p4y;
	IT m12      = tma_a - tma_b;
	IT mt1      = m01 * p3p4z;
	IT mt2      = m02 * p3p4y;
	IT mt3      = m12 * p3p4x;
	IT mtt      = mt2 - mt1;
	IT m012     = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2, ET p3x, ET p3y,
                         ET p3z, ET p4x, ET p4y, ET p4z)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	ET b1p4x    = b1x - p4x;
	ET b1p4y    = b1y - p4y;
	ET b1p4z    = b1z - p4z;
	ET b2p4x    = b2x - p4x;
	ET b2p4y    = b2y - p4y;
	ET b2p4z    = b2z - p4z;
	ET d1_b1p4x = d1 * b1p4x;
	ET d1_b1p4y = d1 * b1p4y;
	ET d1_b1p4z = d1 * b1p4z;
	ET d2_b2p4x = d2 * b2p4x;
	ET d2_b2p4y = d2 * b2p4y;
	ET d2_b2p4z = d2 * b2p4z;
	ET p1p4x    = l1x + d1_b1p4x;
	ET p1p4y    = l1y + d1_b1p4y;
	ET p1p4z    = l1z + d1_b1p4z;
	ET p2p4x    = l2x + d2_b2p4x;
	ET p2p4y    = l2y + d2_b2p4y;
	ET p2p4z    = l2z + d2_b2p4z;
	ET p3p4x    = p3x - p4x;
	ET p3p4y    = p3y - p4y;
	ET p3p4z    = p3z - p4z;
	ET tmc_a    = p1p4x * p2p4y;
	ET tmc_b    = p1p4y * p2p4x;
	ET m01      = tmc_a - tmc_b;
	ET tmi_a    = p1p4x * p2p4z;
	ET tmi_b    = p1p4z * p2p4x;
	ET m02      = tmi_a - tmi_b;
	ET tma_a    = p1p4y * p2p4z;
	ET tma_b    = p1p4z * p2p4y;
	ET m12      = tma_a - tma_b;
	ET mt1      = m01 * p3p4z;
	ET mt2      = m02 * p3p4y;
	ET mt3      = m12 * p3p4x;
	ET mtt      = mt2 - mt1;
	ET m012     = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIEE_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, double p3x,
                             double p3y, double p3z, double p4x, double p4y,
                             double p4z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[32], *l1x = l1x_p, l1y_p[32], *l1y = l1y_p, l1z_p[32],
	                  *l1z = l1z_p, d1_p[32], *d1 = d1_p, b1x, b1y, b1z,
	                  l2x_p[32], *l2x = l2x_p, l2y_p[32], *l2y = l2y_p, l2z_p[32],
	                  *l2z = l2z_p, d2_p[32], *d2 = d2_p, b2x, b2y, b2z;
	int l1x_len = 32, l1y_len = 32, l1z_len = 32, d1_len = 32, l2x_len = 32,
	    l2y_len = 32, l2z_len = 32, d2_len = 32;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0))
	{
		expansionObject o;
		double          b1p4x[2];
		o.two_Diff(b1x, p4x, b1p4x);
		double b1p4y[2];
		o.two_Diff(b1y, p4y, b1p4y);
		double b1p4z[2];
		o.two_Diff(b1z, p4z, b1p4z);
		double b2p4x[2];
		o.two_Diff(b2x, p4x, b2p4x);
		double b2p4y[2];
		o.two_Diff(b2y, p4y, b2p4y);
		double b2p4z[2];
		o.two_Diff(b2z, p4z, b2p4z);
		double d1_b1p4x_p[32], *d1_b1p4x = d1_b1p4x_p;
		int    d1_b1p4x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4x, &d1_b1p4x, 32);
		double d1_b1p4y_p[32], *d1_b1p4y = d1_b1p4y_p;
		int    d1_b1p4y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4y, &d1_b1p4y, 32);
		double d1_b1p4z_p[32], *d1_b1p4z = d1_b1p4z_p;
		int    d1_b1p4z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4z, &d1_b1p4z, 32);
		double d2_b2p4x_p[32], *d2_b2p4x = d2_b2p4x_p;
		int    d2_b2p4x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p4x, &d2_b2p4x, 32);
		double d2_b2p4y_p[32], *d2_b2p4y = d2_b2p4y_p;
		int    d2_b2p4y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p4y, &d2_b2p4y, 32);
		double d2_b2p4z_p[32], *d2_b2p4z = d2_b2p4z_p;
		int    d2_b2p4z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p4z, &d2_b2p4z, 32);
		double p1p4x_p[32], *p1p4x = p1p4x_p;
		int    p1p4x_len =
		  o.Gen_Sum_With_PreAlloc(l1x_len, l1x, d1_b1p4x_len, d1_b1p4x, &p1p4x, 32);
		double p1p4y_p[32], *p1p4y = p1p4y_p;
		int    p1p4y_len =
		  o.Gen_Sum_With_PreAlloc(l1y_len, l1y, d1_b1p4y_len, d1_b1p4y, &p1p4y, 32);
		double p1p4z_p[32], *p1p4z = p1p4z_p;
		int    p1p4z_len =
		  o.Gen_Sum_With_PreAlloc(l1z_len, l1z, d1_b1p4z_len, d1_b1p4z, &p1p4z, 32);
		double p2p4x_p[32], *p2p4x = p2p4x_p;
		int    p2p4x_len =
		  o.Gen_Sum_With_PreAlloc(l2x_len, l2x, d2_b2p4x_len, d2_b2p4x, &p2p4x, 32);
		double p2p4y_p[32], *p2p4y = p2p4y_p;
		int    p2p4y_len =
		  o.Gen_Sum_With_PreAlloc(l2y_len, l2y, d2_b2p4y_len, d2_b2p4y, &p2p4y, 32);
		double p2p4z_p[32], *p2p4z = p2p4z_p;
		int    p2p4z_len =
		  o.Gen_Sum_With_PreAlloc(l2z_len, l2z, d2_b2p4z_len, d2_b2p4z, &p2p4z, 32);
		double p3p4x[2];
		o.two_Diff(p3x, p4x, p3p4x);
		double p3p4y[2];
		o.two_Diff(p3y, p4y, p3p4y);
		double p3p4z[2];
		o.two_Diff(p3z, p4z, p3p4z);
		double tmc_a_p[32], *tmc_a = tmc_a_p;
		int    tmc_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4y_len,
		                                               p2p4y, &tmc_a, 32);
		double tmc_b_p[32], *tmc_b = tmc_b_p;
		int    tmc_b_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4x_len,
		                                               p2p4x, &tmc_b, 32);
		double m01_p[32], *m01 = m01_p;
		int    m01_len =
		  o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 32);
		double tmi_a_p[32], *tmi_a = tmi_a_p;
		int    tmi_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4z_len,
		                                               p2p4z, &tmi_a, 32);
		double tmi_b_p[32], *tmi_b = tmi_b_p;
		int    tmi_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4x_len,
		                                               p2p4x, &tmi_b, 32);
		double m02_p[32], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 32);
		double tma_a_p[32], *tma_a = tma_a_p;
		int    tma_a_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4z_len,
		                                               p2p4z, &tma_a, 32);
		double tma_b_p[32], *tma_b = tma_b_p;
		int    tma_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4y_len,
		                                               p2p4y, &tma_b, 32);
		double m12_p[32], *m12 = m12_p;
		int    m12_len =
		  o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 32);
		double mt1_p[32], *mt1 = mt1_p;
		int mt1_len = o.Gen_Product_With_PreAlloc(m01_len, m01, 2, p3p4z, &mt1, 32);
		double mt2_p[32], *mt2 = mt2_p;
		int mt2_len = o.Gen_Product_With_PreAlloc(m02_len, m02, 2, p3p4y, &mt2, 32);
		double mt3_p[32], *mt3 = mt3_p;
		int mt3_len = o.Gen_Product_With_PreAlloc(m12_len, m12, 2, p3p4x, &mt3, 32);
		double mtt_p[32], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 32);
		double m012_p[32], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 32);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tma_b_p != tma_b)
			FreeDoubles(tma_b);
		if (tma_a_p != tma_a)
			FreeDoubles(tma_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmc_b_p != tmc_b)
			FreeDoubles(tmc_b);
		if (tmc_a_p != tmc_a)
			FreeDoubles(tmc_a);
		if (p2p4z_p != p2p4z)
			FreeDoubles(p2p4z);
		if (p2p4y_p != p2p4y)
			FreeDoubles(p2p4y);
		if (p2p4x_p != p2p4x)
			FreeDoubles(p2p4x);
		if (p1p4z_p != p1p4z)
			FreeDoubles(p1p4z);
		if (p1p4y_p != p1p4y)
			FreeDoubles(p1p4y);
		if (p1p4x_p != p1p4x)
			FreeDoubles(p1p4x);
		if (d2_b2p4z_p != d2_b2p4z)
			FreeDoubles(d2_b2p4z);
		if (d2_b2p4y_p != d2_b2p4y)
			FreeDoubles(d2_b2p4y);
		if (d2_b2p4x_p != d2_b2p4x)
			FreeDoubles(d2_b2p4x);
		if (d1_b1p4z_p != d1_b1p4z)
			FreeDoubles(d1_b1p4z);
		if (d1_b1p4y_p != d1_b1p4y)
			FreeDoubles(d1_b1p4y);
		if (d1_b1p4x_p != d1_b1p4x)
			FreeDoubles(d1_b1p4x);
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
		return orient3D_IIEE_exact<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2, double p3x, double p3y,
                   double p3z, double p4x, double p4y, double p4z, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret =
		  orient3D_IIEE_filtered<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IIEE_interval<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IIEE_expansion<IT, ET>(p1, p2, p3x, p3y, p3z, p4x, p4y, p4z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	return orient3D_IIEE<IT, ET, WithSSFilter>(p1, p2, p3.x(), p3.y(), p3.z(),
	                                           p4.x(), p4.y(), p4.z(), arr);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3, double p4x,
                            double p4y, double p4z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var))
		return Sign::UNCERTAIN;

	double b1p4x    = b1x - p4x;
	double b1p4y    = b1y - p4y;
	double b1p4z    = b1z - p4z;
	double b2p4x    = b2x - p4x;
	double b2p4y    = b2y - p4y;
	double b2p4z    = b2z - p4z;
	double b3p4x    = b3x - p4x;
	double b3p4y    = b3y - p4y;
	double b3p4z    = b3z - p4z;
	double d1_b1p4x = d1 * b1p4x;
	double d1_b1p4y = d1 * b1p4y;
	double d1_b1p4z = d1 * b1p4z;
	double d2_b2p4x = d2 * b2p4x;
	double d2_b2p4y = d2 * b2p4y;
	double d2_b2p4z = d2 * b2p4z;
	double d3_b3p4x = d3 * b3p4x;
	double d3_b3p4y = d3 * b3p4y;
	double d3_b3p4z = d3 * b3p4z;
	double p1p4x    = l1x + d1_b1p4x;
	double p1p4y    = l1y + d1_b1p4y;
	double p1p4z    = l1z + d1_b1p4z;
	double p2p4x    = l2x + d2_b2p4x;
	double p2p4y    = l2y + d2_b2p4y;
	double p2p4z    = l2z + d2_b2p4z;
	double p3p4x    = l3x + d3_b3p4x;
	double p3p4y    = l3y + d3_b3p4y;
	double p3p4z    = l3z + d3_b3p4z;
	double tmc_a    = p1p4x * p2p4y;
	double tmc_b    = p1p4y * p2p4x;
	double m01      = tmc_a - tmc_b;
	double tmi_a    = p1p4x * p2p4z;
	double tmi_b    = p1p4z * p2p4x;
	double m02      = tmi_a - tmi_b;
	double tma_a    = p1p4y * p2p4z;
	double tma_b    = p1p4z * p2p4y;
	double m12      = tma_a - tma_b;
	double mt1      = m01 * p3p4z;
	double mt2      = m02 * p3p4y;
	double mt3      = m12 * p3p4x;
	double mtt      = mt2 - mt1;
	double m012     = mtt - mt3;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b3p4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b3p4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b3p4z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 1.0800249583553551e-12;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.538638132516435e-12;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 8.586198418925086e-11;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1412128186894e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.7186095735709623e-10;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 6.132438556960583e-09;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.747767721162821e-11;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 8.746082608146503e-10;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.956868465879957e-08;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.155571105002337e-07;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3, IT p4x, IT p4y,
                            IT p4z)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p4x    = b1x - p4x;
	IT b1p4y    = b1y - p4y;
	IT b1p4z    = b1z - p4z;
	IT b2p4x    = b2x - p4x;
	IT b2p4y    = b2y - p4y;
	IT b2p4z    = b2z - p4z;
	IT b3p4x    = b3x - p4x;
	IT b3p4y    = b3y - p4y;
	IT b3p4z    = b3z - p4z;
	IT d1_b1p4x = d1 * b1p4x;
	IT d1_b1p4y = d1 * b1p4y;
	IT d1_b1p4z = d1 * b1p4z;
	IT d2_b2p4x = d2 * b2p4x;
	IT d2_b2p4y = d2 * b2p4y;
	IT d2_b2p4z = d2 * b2p4z;
	IT d3_b3p4x = d3 * b3p4x;
	IT d3_b3p4y = d3 * b3p4y;
	IT d3_b3p4z = d3 * b3p4z;
	IT p1p4x    = l1x + d1_b1p4x;
	IT p1p4y    = l1y + d1_b1p4y;
	IT p1p4z    = l1z + d1_b1p4z;
	IT p2p4x    = l2x + d2_b2p4x;
	IT p2p4y    = l2y + d2_b2p4y;
	IT p2p4z    = l2z + d2_b2p4z;
	IT p3p4x    = l3x + d3_b3p4x;
	IT p3p4y    = l3y + d3_b3p4y;
	IT p3p4z    = l3z + d3_b3p4z;
	IT tmc_a    = p1p4x * p2p4y;
	IT tmc_b    = p1p4y * p2p4x;
	IT m01      = tmc_a - tmc_b;
	IT tmi_a    = p1p4x * p2p4z;
	IT tmi_b    = p1p4z * p2p4x;
	IT m02      = tmi_a - tmi_b;
	IT tma_a    = p1p4y * p2p4z;
	IT tma_b    = p1p4z * p2p4y;
	IT m12      = tma_a - tma_b;
	IT mt1      = m01 * p3p4z;
	IT mt2      = m02 * p3p4y;
	IT mt3      = m12 * p3p4x;
	IT mtt      = mt2 - mt1;
	IT m012     = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2,
                         const GenericPoint3T<IT, ET> &p3, ET p4x, ET p4y,
                         ET p4z)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	p3.getExactLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z);
	ET b1p4x    = b1x - p4x;
	ET b1p4y    = b1y - p4y;
	ET b1p4z    = b1z - p4z;
	ET b2p4x    = b2x - p4x;
	ET b2p4y    = b2y - p4y;
	ET b2p4z    = b2z - p4z;
	ET b3p4x    = b3x - p4x;
	ET b3p4y    = b3y - p4y;
	ET b3p4z    = b3z - p4z;
	ET d1_b1p4x = d1 * b1p4x;
	ET d1_b1p4y = d1 * b1p4y;
	ET d1_b1p4z = d1 * b1p4z;
	ET d2_b2p4x = d2 * b2p4x;
	ET d2_b2p4y = d2 * b2p4y;
	ET d2_b2p4z = d2 * b2p4z;
	ET d3_b3p4x = d3 * b3p4x;
	ET d3_b3p4y = d3 * b3p4y;
	ET d3_b3p4z = d3 * b3p4z;
	ET p1p4x    = l1x + d1_b1p4x;
	ET p1p4y    = l1y + d1_b1p4y;
	ET p1p4z    = l1z + d1_b1p4z;
	ET p2p4x    = l2x + d2_b2p4x;
	ET p2p4y    = l2y + d2_b2p4y;
	ET p2p4z    = l2z + d2_b2p4z;
	ET p3p4x    = l3x + d3_b3p4x;
	ET p3p4y    = l3y + d3_b3p4y;
	ET p3p4z    = l3z + d3_b3p4z;
	ET tmc_a    = p1p4x * p2p4y;
	ET tmc_b    = p1p4y * p2p4x;
	ET m01      = tmc_a - tmc_b;
	ET tmi_a    = p1p4x * p2p4z;
	ET tmi_b    = p1p4z * p2p4x;
	ET m02      = tmi_a - tmi_b;
	ET tma_a    = p1p4y * p2p4z;
	ET tma_b    = p1p4z * p2p4y;
	ET m12      = tma_a - tma_b;
	ET mt1      = m01 * p3p4z;
	ET mt2      = m02 * p3p4y;
	ET mt3      = m12 * p3p4x;
	ET mtt      = mt2 - mt1;
	ET m012     = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIIE_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3, double p4x,
                             double p4y, double p4z)
{
	double return_value = NAN;
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
		double          b1p4x[2];
		o.two_Diff(b1x, p4x, b1p4x);
		double b1p4y[2];
		o.two_Diff(b1y, p4y, b1p4y);
		double b1p4z[2];
		o.two_Diff(b1z, p4z, b1p4z);
		double b2p4x[2];
		o.two_Diff(b2x, p4x, b2p4x);
		double b2p4y[2];
		o.two_Diff(b2y, p4y, b2p4y);
		double b2p4z[2];
		o.two_Diff(b2z, p4z, b2p4z);
		double b3p4x[2];
		o.two_Diff(b3x, p4x, b3p4x);
		double b3p4y[2];
		o.two_Diff(b3y, p4y, b3p4y);
		double b3p4z[2];
		o.two_Diff(b3z, p4z, b3p4z);
		double d1_b1p4x_p[32], *d1_b1p4x = d1_b1p4x_p;
		int    d1_b1p4x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4x, &d1_b1p4x, 32);
		double d1_b1p4y_p[32], *d1_b1p4y = d1_b1p4y_p;
		int    d1_b1p4y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4y, &d1_b1p4y, 32);
		double d1_b1p4z_p[32], *d1_b1p4z = d1_b1p4z_p;
		int    d1_b1p4z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p4z, &d1_b1p4z, 32);
		double d2_b2p4x_p[32], *d2_b2p4x = d2_b2p4x_p;
		int    d2_b2p4x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p4x, &d2_b2p4x, 32);
		double d2_b2p4y_p[32], *d2_b2p4y = d2_b2p4y_p;
		int    d2_b2p4y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p4y, &d2_b2p4y, 32);
		double d2_b2p4z_p[32], *d2_b2p4z = d2_b2p4z_p;
		int    d2_b2p4z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2p4z, &d2_b2p4z, 32);
		double d3_b3p4x_p[32], *d3_b3p4x = d3_b3p4x_p;
		int    d3_b3p4x_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, 2, b3p4x, &d3_b3p4x, 32);
		double d3_b3p4y_p[32], *d3_b3p4y = d3_b3p4y_p;
		int    d3_b3p4y_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, 2, b3p4y, &d3_b3p4y, 32);
		double d3_b3p4z_p[32], *d3_b3p4z = d3_b3p4z_p;
		int    d3_b3p4z_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, 2, b3p4z, &d3_b3p4z, 32);
		double p1p4x_p[32], *p1p4x = p1p4x_p;
		int    p1p4x_len =
		  o.Gen_Sum_With_PreAlloc(l1x_len, l1x, d1_b1p4x_len, d1_b1p4x, &p1p4x, 32);
		double p1p4y_p[32], *p1p4y = p1p4y_p;
		int    p1p4y_len =
		  o.Gen_Sum_With_PreAlloc(l1y_len, l1y, d1_b1p4y_len, d1_b1p4y, &p1p4y, 32);
		double p1p4z_p[32], *p1p4z = p1p4z_p;
		int    p1p4z_len =
		  o.Gen_Sum_With_PreAlloc(l1z_len, l1z, d1_b1p4z_len, d1_b1p4z, &p1p4z, 32);
		double p2p4x_p[32], *p2p4x = p2p4x_p;
		int    p2p4x_len =
		  o.Gen_Sum_With_PreAlloc(l2x_len, l2x, d2_b2p4x_len, d2_b2p4x, &p2p4x, 32);
		double p2p4y_p[32], *p2p4y = p2p4y_p;
		int    p2p4y_len =
		  o.Gen_Sum_With_PreAlloc(l2y_len, l2y, d2_b2p4y_len, d2_b2p4y, &p2p4y, 32);
		double p2p4z_p[32], *p2p4z = p2p4z_p;
		int    p2p4z_len =
		  o.Gen_Sum_With_PreAlloc(l2z_len, l2z, d2_b2p4z_len, d2_b2p4z, &p2p4z, 32);
		double p3p4x_p[32], *p3p4x = p3p4x_p;
		int    p3p4x_len =
		  o.Gen_Sum_With_PreAlloc(l3x_len, l3x, d3_b3p4x_len, d3_b3p4x, &p3p4x, 32);
		double p3p4y_p[32], *p3p4y = p3p4y_p;
		int    p3p4y_len =
		  o.Gen_Sum_With_PreAlloc(l3y_len, l3y, d3_b3p4y_len, d3_b3p4y, &p3p4y, 32);
		double p3p4z_p[32], *p3p4z = p3p4z_p;
		int    p3p4z_len =
		  o.Gen_Sum_With_PreAlloc(l3z_len, l3z, d3_b3p4z_len, d3_b3p4z, &p3p4z, 32);
		double tmc_a_p[32], *tmc_a = tmc_a_p;
		int    tmc_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4y_len,
		                                               p2p4y, &tmc_a, 32);
		double tmc_b_p[32], *tmc_b = tmc_b_p;
		int    tmc_b_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4x_len,
		                                               p2p4x, &tmc_b, 32);
		double m01_p[32], *m01 = m01_p;
		int    m01_len =
		  o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 32);
		double tmi_a_p[32], *tmi_a = tmi_a_p;
		int    tmi_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4z_len,
		                                               p2p4z, &tmi_a, 32);
		double tmi_b_p[32], *tmi_b = tmi_b_p;
		int    tmi_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4x_len,
		                                               p2p4x, &tmi_b, 32);
		double m02_p[32], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 32);
		double tma_a_p[32], *tma_a = tma_a_p;
		int    tma_a_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4z_len,
		                                               p2p4z, &tma_a, 32);
		double tma_b_p[32], *tma_b = tma_b_p;
		int    tma_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4y_len,
		                                               p2p4y, &tma_b, 32);
		double m12_p[32], *m12 = m12_p;
		int    m12_len =
		  o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 32);
		double mt1_p[32], *mt1 = mt1_p;
		int    mt1_len =
		  o.Gen_Product_With_PreAlloc(m01_len, m01, p3p4z_len, p3p4z, &mt1, 32);
		double mt2_p[32], *mt2 = mt2_p;
		int    mt2_len =
		  o.Gen_Product_With_PreAlloc(m02_len, m02, p3p4y_len, p3p4y, &mt2, 32);
		double mt3_p[32], *mt3 = mt3_p;
		int    mt3_len =
		  o.Gen_Product_With_PreAlloc(m12_len, m12, p3p4x_len, p3p4x, &mt3, 32);
		double mtt_p[32], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 32);
		double m012_p[32], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 32);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tma_b_p != tma_b)
			FreeDoubles(tma_b);
		if (tma_a_p != tma_a)
			FreeDoubles(tma_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmc_b_p != tmc_b)
			FreeDoubles(tmc_b);
		if (tmc_a_p != tmc_a)
			FreeDoubles(tmc_a);
		if (p3p4z_p != p3p4z)
			FreeDoubles(p3p4z);
		if (p3p4y_p != p3p4y)
			FreeDoubles(p3p4y);
		if (p3p4x_p != p3p4x)
			FreeDoubles(p3p4x);
		if (p2p4z_p != p2p4z)
			FreeDoubles(p2p4z);
		if (p2p4y_p != p2p4y)
			FreeDoubles(p2p4y);
		if (p2p4x_p != p2p4x)
			FreeDoubles(p2p4x);
		if (p1p4z_p != p1p4z)
			FreeDoubles(p1p4z);
		if (p1p4y_p != p1p4y)
			FreeDoubles(p1p4y);
		if (p1p4x_p != p1p4x)
			FreeDoubles(p1p4x);
		if (d3_b3p4z_p != d3_b3p4z)
			FreeDoubles(d3_b3p4z);
		if (d3_b3p4y_p != d3_b3p4y)
			FreeDoubles(d3_b3p4y);
		if (d3_b3p4x_p != d3_b3p4x)
			FreeDoubles(d3_b3p4x);
		if (d2_b2p4z_p != d2_b2p4z)
			FreeDoubles(d2_b2p4z);
		if (d2_b2p4y_p != d2_b2p4y)
			FreeDoubles(d2_b2p4y);
		if (d2_b2p4x_p != d2_b2p4x)
			FreeDoubles(d2_b2p4x);
		if (d1_b1p4z_p != d1_b1p4z)
			FreeDoubles(d1_b1p4z);
		if (d1_b1p4y_p != d1_b1p4y)
			FreeDoubles(d1_b1p4y);
		if (d1_b1p4x_p != d1_b1p4x)
			FreeDoubles(d1_b1p4x);
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
		return orient3D_IIIE_exact<IT, ET>(p1, p2, p3, p4x, p4y, p4z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIIE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3, double p4x, double p4y,
                   double p4z, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient3D_IIIE_filtered<IT, ET>(p1, p2, p3, p4x, p4y, p4z, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IIIE_interval<IT, ET>(p1, p2, p3, p4x, p4y, p4z);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IIIE_expansion<IT, ET>(p1, p2, p3, p4x, p4y, p4z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIIE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	return orient3D_IIIE<IT, ET, WithSSFilter>(p1, p2, p3, p4.x(), p4.y(), p4.z(),
	                                           arr);
}

template <typename IT, typename ET>
Sign orient3D_IIII_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3,
                            const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  l3x, l3y, l3z, d3, b3x, b3y, b3z, l4x, l4y, l4z, d4, b4x, b4y, b4z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var) ||
	    !p4.getFilteredLambda(l4x, l4y, l4z, d4, b4x, b4y, b4z, max_var))
		return Sign::UNCERTAIN;

	double d1p4x         = d1 * l4x;
	double d1p4y         = d1 * l4y;
	double d1p4z         = d1 * l4z;
	double d2p4x         = d2 * l4x;
	double d2p4y         = d2 * l4y;
	double d2p4z         = d2 * l4z;
	double d3p4x         = d3 * l4x;
	double d3p4y         = d3 * l4y;
	double d3p4z         = d3 * l4z;
	double b1b4x         = b1x - b4x;
	double b1b4y         = b1y - b4y;
	double b1b4z         = b1z - b4z;
	double b2b4x         = b2x - b4x;
	double b2b4y         = b2y - b4y;
	double b2b4z         = b2z - b4z;
	double b3b4x         = b3x - b4x;
	double b3b4y         = b3y - b4y;
	double b3b4z         = b3z - b4z;
	double d1_b1b4x      = d1 * b1b4x;
	double d1_b1b4y      = d1 * b1b4y;
	double d1_b1b4z      = d1 * b1b4z;
	double d2_b2b4x      = d2 * b2b4x;
	double d2_b2b4y      = d2 * b2b4y;
	double d2_b2b4z      = d2 * b2b4z;
	double d3_b3b4x      = d3 * b3b4x;
	double d3_b3b4y      = d3 * b3b4y;
	double d3_b3b4z      = d3 * b3b4z;
	double d1_b1b4_l1x   = d1_b1b4x + l1x;
	double d1_b1b4_l1y   = d1_b1b4y + l1y;
	double d1_b1b4_l1z   = d1_b1b4z + l1z;
	double d2_b2b4_l2x   = d2_b2b4x + l2x;
	double d2_b2b4_l2y   = d2_b2b4y + l2y;
	double d2_b2b4_l2z   = d2_b2b4z + l2z;
	double d3_b3b4_l3x   = d3_b3b4x + l3x;
	double d3_b3b4_l3y   = d3_b3b4y + l3y;
	double d3_b3b4_l3z   = d3_b3b4z + l3z;
	double d4d1_b1b4_l1x = d4 * d1_b1b4_l1x;
	double d4d1_b1b4_l1y = d4 * d1_b1b4_l1y;
	double d4d1_b1b4_l1z = d4 * d1_b1b4_l1z;
	double d4d2_b2b4_l2x = d4 * d2_b2b4_l2x;
	double d4d2_b2b4_l2y = d4 * d2_b2b4_l2y;
	double d4d2_b2b4_l2z = d4 * d2_b2b4_l2z;
	double d4d3_b3b4_l3x = d4 * d3_b3b4_l3x;
	double d4d3_b3b4_l3y = d4 * d3_b3b4_l3y;
	double d4d3_b3b4_l3z = d4 * d3_b3b4_l3z;
	double p1p4x         = d4d1_b1b4_l1x - d1p4x;
	double p1p4y         = d4d1_b1b4_l1y - d1p4y;
	double p1p4z         = d4d1_b1b4_l1z - d1p4z;
	double p2p4x         = d4d2_b2b4_l2x - d2p4x;
	double p2p4y         = d4d2_b2b4_l2y - d2p4y;
	double p2p4z         = d4d2_b2b4_l2z - d2p4z;
	double p3p4x         = d4d3_b3b4_l3x - d3p4x;
	double p3p4y         = d4d3_b3b4_l3y - d3p4y;
	double p3p4z         = d4d3_b3b4_l3z - d3p4z;
	double tmc_a         = p1p4x * p2p4y;
	double tmc_b         = p1p4y * p2p4x;
	double m01           = tmc_a - tmc_b;
	double tmi_a         = p1p4x * p2p4z;
	double tmi_b         = p1p4z * p2p4x;
	double m02           = tmi_a - tmi_b;
	double tma_a         = p1p4y * p2p4z;
	double tma_b         = p1p4z * p2p4y;
	double m12           = tma_a - tma_b;
	double mt1           = m01 * p3p4z;
	double mt2           = m02 * p3p4y;
	double mt3           = m12 * p3p4x;
	double mtt           = mt2 - mt1;
	double m012          = mtt - mt3;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b4z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b3b4x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b3b4y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b3b4z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.668265773943881e-11;
	}
	break;
	case PntArr3::SSSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.554727971364587e-09;
	}
	break;
	case PntArr3::SSST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.0059331543743702e-05;
	}
	break;
	case PntArr3::SSLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.934037456383351e-09;
	}
	break;
	case PntArr3::SSLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.112779300096265e-05;
	}
	break;
	case PntArr3::SSTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.0005226793582551237;
	}
	break;
	case PntArr3::SLLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.535519549011886e-08;
	}
	break;
	case PntArr3::SLLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 9.51822422479163e-05;
	}
	break;
	case PntArr3::SLTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 0.0015946633955721441;
	}
	break;
	case PntArr3::STTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.026444148039445593;
	}
	break;
	case PntArr3::LLLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 4.8781668676723156e-08;
	}
	break;
	case PntArr3::LLLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.0002952721038056467;
	}
	break;
	case PntArr3::LLTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 0.004930109261224001;
	}
	break;
	case PntArr3::LTTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 0.08154614984406648;
	}
	break;
	case PntArr3::TTTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.327523517713434;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(m012, epsilon);
}

template <typename IT, typename ET>
Sign orient3D_IIII_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3,
                            const GenericPoint3T<IT, ET> &p4)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z, l4x, l4y, l4z, d4, b4x, b4y, b4z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z) ||
	    !p4.getIntervalLambda(l4x, l4y, l4z, d4, b4x, b4y, b4z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT d1p4x         = d1 * l4x;
	IT d1p4y         = d1 * l4y;
	IT d1p4z         = d1 * l4z;
	IT d2p4x         = d2 * l4x;
	IT d2p4y         = d2 * l4y;
	IT d2p4z         = d2 * l4z;
	IT d3p4x         = d3 * l4x;
	IT d3p4y         = d3 * l4y;
	IT d3p4z         = d3 * l4z;
	IT b1b4x         = b1x - b4x;
	IT b1b4y         = b1y - b4y;
	IT b1b4z         = b1z - b4z;
	IT b2b4x         = b2x - b4x;
	IT b2b4y         = b2y - b4y;
	IT b2b4z         = b2z - b4z;
	IT b3b4x         = b3x - b4x;
	IT b3b4y         = b3y - b4y;
	IT b3b4z         = b3z - b4z;
	IT d1_b1b4x      = d1 * b1b4x;
	IT d1_b1b4y      = d1 * b1b4y;
	IT d1_b1b4z      = d1 * b1b4z;
	IT d2_b2b4x      = d2 * b2b4x;
	IT d2_b2b4y      = d2 * b2b4y;
	IT d2_b2b4z      = d2 * b2b4z;
	IT d3_b3b4x      = d3 * b3b4x;
	IT d3_b3b4y      = d3 * b3b4y;
	IT d3_b3b4z      = d3 * b3b4z;
	IT d1_b1b4_l1x   = d1_b1b4x + l1x;
	IT d1_b1b4_l1y   = d1_b1b4y + l1y;
	IT d1_b1b4_l1z   = d1_b1b4z + l1z;
	IT d2_b2b4_l2x   = d2_b2b4x + l2x;
	IT d2_b2b4_l2y   = d2_b2b4y + l2y;
	IT d2_b2b4_l2z   = d2_b2b4z + l2z;
	IT d3_b3b4_l3x   = d3_b3b4x + l3x;
	IT d3_b3b4_l3y   = d3_b3b4y + l3y;
	IT d3_b3b4_l3z   = d3_b3b4z + l3z;
	IT d4d1_b1b4_l1x = d4 * d1_b1b4_l1x;
	IT d4d1_b1b4_l1y = d4 * d1_b1b4_l1y;
	IT d4d1_b1b4_l1z = d4 * d1_b1b4_l1z;
	IT d4d2_b2b4_l2x = d4 * d2_b2b4_l2x;
	IT d4d2_b2b4_l2y = d4 * d2_b2b4_l2y;
	IT d4d2_b2b4_l2z = d4 * d2_b2b4_l2z;
	IT d4d3_b3b4_l3x = d4 * d3_b3b4_l3x;
	IT d4d3_b3b4_l3y = d4 * d3_b3b4_l3y;
	IT d4d3_b3b4_l3z = d4 * d3_b3b4_l3z;
	IT p1p4x         = d4d1_b1b4_l1x - d1p4x;
	IT p1p4y         = d4d1_b1b4_l1y - d1p4y;
	IT p1p4z         = d4d1_b1b4_l1z - d1p4z;
	IT p2p4x         = d4d2_b2b4_l2x - d2p4x;
	IT p2p4y         = d4d2_b2b4_l2y - d2p4y;
	IT p2p4z         = d4d2_b2b4_l2z - d2p4z;
	IT p3p4x         = d4d3_b3b4_l3x - d3p4x;
	IT p3p4y         = d4d3_b3b4_l3y - d3p4y;
	IT p3p4z         = d4d3_b3b4_l3z - d3p4z;
	IT tmc_a         = p1p4x * p2p4y;
	IT tmc_b         = p1p4y * p2p4x;
	IT m01           = tmc_a - tmc_b;
	IT tmi_a         = p1p4x * p2p4z;
	IT tmi_b         = p1p4z * p2p4x;
	IT m02           = tmi_a - tmi_b;
	IT tma_a         = p1p4y * p2p4z;
	IT tma_b         = p1p4z * p2p4y;
	IT m12           = tma_a - tma_b;
	IT mt1           = m01 * p3p4z;
	IT mt2           = m02 * p3p4y;
	IT mt3           = m12 * p3p4x;
	IT mtt           = mt2 - mt1;
	IT m012          = mtt - mt3;
	if (!m012.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIII_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2,
                         const GenericPoint3T<IT, ET> &p3,
                         const GenericPoint3T<IT, ET> &p4)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z, l4x, l4y, l4z, d4, b4x, b4y, b4z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	p3.getExactLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z);
	p4.getExactLambda(l4x, l4y, l4z, d4, b4x, b4y, b4z);
	ET d1p4x         = d1 * l4x;
	ET d1p4y         = d1 * l4y;
	ET d1p4z         = d1 * l4z;
	ET d2p4x         = d2 * l4x;
	ET d2p4y         = d2 * l4y;
	ET d2p4z         = d2 * l4z;
	ET d3p4x         = d3 * l4x;
	ET d3p4y         = d3 * l4y;
	ET d3p4z         = d3 * l4z;
	ET b1b4x         = b1x - b4x;
	ET b1b4y         = b1y - b4y;
	ET b1b4z         = b1z - b4z;
	ET b2b4x         = b2x - b4x;
	ET b2b4y         = b2y - b4y;
	ET b2b4z         = b2z - b4z;
	ET b3b4x         = b3x - b4x;
	ET b3b4y         = b3y - b4y;
	ET b3b4z         = b3z - b4z;
	ET d1_b1b4x      = d1 * b1b4x;
	ET d1_b1b4y      = d1 * b1b4y;
	ET d1_b1b4z      = d1 * b1b4z;
	ET d2_b2b4x      = d2 * b2b4x;
	ET d2_b2b4y      = d2 * b2b4y;
	ET d2_b2b4z      = d2 * b2b4z;
	ET d3_b3b4x      = d3 * b3b4x;
	ET d3_b3b4y      = d3 * b3b4y;
	ET d3_b3b4z      = d3 * b3b4z;
	ET d1_b1b4_l1x   = d1_b1b4x + l1x;
	ET d1_b1b4_l1y   = d1_b1b4y + l1y;
	ET d1_b1b4_l1z   = d1_b1b4z + l1z;
	ET d2_b2b4_l2x   = d2_b2b4x + l2x;
	ET d2_b2b4_l2y   = d2_b2b4y + l2y;
	ET d2_b2b4_l2z   = d2_b2b4z + l2z;
	ET d3_b3b4_l3x   = d3_b3b4x + l3x;
	ET d3_b3b4_l3y   = d3_b3b4y + l3y;
	ET d3_b3b4_l3z   = d3_b3b4z + l3z;
	ET d4d1_b1b4_l1x = d4 * d1_b1b4_l1x;
	ET d4d1_b1b4_l1y = d4 * d1_b1b4_l1y;
	ET d4d1_b1b4_l1z = d4 * d1_b1b4_l1z;
	ET d4d2_b2b4_l2x = d4 * d2_b2b4_l2x;
	ET d4d2_b2b4_l2y = d4 * d2_b2b4_l2y;
	ET d4d2_b2b4_l2z = d4 * d2_b2b4_l2z;
	ET d4d3_b3b4_l3x = d4 * d3_b3b4_l3x;
	ET d4d3_b3b4_l3y = d4 * d3_b3b4_l3y;
	ET d4d3_b3b4_l3z = d4 * d3_b3b4_l3z;
	ET p1p4x         = d4d1_b1b4_l1x - d1p4x;
	ET p1p4y         = d4d1_b1b4_l1y - d1p4y;
	ET p1p4z         = d4d1_b1b4_l1z - d1p4z;
	ET p2p4x         = d4d2_b2b4_l2x - d2p4x;
	ET p2p4y         = d4d2_b2b4_l2y - d2p4y;
	ET p2p4z         = d4d2_b2b4_l2z - d2p4z;
	ET p3p4x         = d4d3_b3b4_l3x - d3p4x;
	ET p3p4y         = d4d3_b3b4_l3y - d3p4y;
	ET p3p4z         = d4d3_b3b4_l3z - d3p4z;
	ET tmc_a         = p1p4x * p2p4y;
	ET tmc_b         = p1p4y * p2p4x;
	ET m01           = tmc_a - tmc_b;
	ET tmi_a         = p1p4x * p2p4z;
	ET tmi_b         = p1p4z * p2p4x;
	ET m02           = tmi_a - tmi_b;
	ET tma_a         = p1p4y * p2p4z;
	ET tma_b         = p1p4z * p2p4y;
	ET m12           = tma_a - tma_b;
	ET mt1           = m01 * p3p4z;
	ET mt2           = m02 * p3p4y;
	ET mt3           = m12 * p3p4x;
	ET mtt           = mt2 - mt1;
	ET m012          = mtt - mt3;
	return OMC::sign(m012);
}

template <typename IT, typename ET>
Sign orient3D_IIII_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3,
                             const GenericPoint3T<IT, ET> &p4)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[16],
	  *l1x = l1x_p, l1y_p[16], *l1y = l1y_p, l1z_p[16], *l1z = l1z_p, d1_p[16],
	  *d1 = d1_p, b1x, b1y, b1z, l2x_p[16], *l2x = l2x_p, l2y_p[16], *l2y = l2y_p,
	  l2z_p[16], *l2z = l2z_p, d2_p[16], *d2 = d2_p, b2x, b2y, b2z, l3x_p[16],
	  *l3x = l3x_p, l3y_p[16], *l3y = l3y_p, l3z_p[16], *l3z = l3z_p, d3_p[16],
	  *d3 = d3_p, b3x, b3y, b3z, l4x_p[16], *l4x = l4x_p, l4y_p[16], *l4y = l4y_p,
	  l4z_p[16], *l4z = l4z_p, d4_p[16], *d4 = d4_p, b4x, b4y, b4z;
	int l1x_len = 16, l1y_len = 16, l1z_len = 16, d1_len = 16, l2x_len = 16,
	    l2y_len = 16, l2z_len = 16, d2_len = 16, l3x_len = 16, l3y_len = 16,
	    l3z_len = 16, d3_len = 16, l4x_len = 16, l4y_len = 16, l4z_len = 16,
	    d4_len = 16;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	p2.getExpansionLambda(&l2x, l2x_len, &l2y, l2y_len, &l2z, l2z_len, &d2,
	                      d2_len, b2x, b2y, b2z);
	p3.getExpansionLambda(&l3x, l3x_len, &l3y, l3y_len, &l3z, l3z_len, &d3,
	                      d3_len, b3x, b3y, b3z);
	p4.getExpansionLambda(&l4x, l4x_len, &l4y, l4y_len, &l4z, l4z_len, &d4,
	                      d4_len, b4x, b4y, b4z);
	if ((d1[d1_len - 1] != 0) && (d2[d2_len - 1] != 0) && (d3[d3_len - 1] != 0) &&
	    (d4[d4_len - 1] != 0))
	{
		expansionObject o;
		double          d1p4x_p[16], *d1p4x = d1p4x_p;
		int             d1p4x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l4x_len, l4x, &d1p4x, 16);
		double d1p4y_p[16], *d1p4y = d1p4y_p;
		int    d1p4y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l4y_len, l4y, &d1p4y, 16);
		double d1p4z_p[16], *d1p4z = d1p4z_p;
		int    d1p4z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, l4z_len, l4z, &d1p4z, 16);
		double d2p4x_p[16], *d2p4x = d2p4x_p;
		int    d2p4x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l4x_len, l4x, &d2p4x, 16);
		double d2p4y_p[16], *d2p4y = d2p4y_p;
		int    d2p4y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l4y_len, l4y, &d2p4y, 16);
		double d2p4z_p[16], *d2p4z = d2p4z_p;
		int    d2p4z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, l4z_len, l4z, &d2p4z, 16);
		double d3p4x_p[16], *d3p4x = d3p4x_p;
		int    d3p4x_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, l4x_len, l4x, &d3p4x, 16);
		double d3p4y_p[16], *d3p4y = d3p4y_p;
		int    d3p4y_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, l4y_len, l4y, &d3p4y, 16);
		double d3p4z_p[16], *d3p4z = d3p4z_p;
		int    d3p4z_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, l4z_len, l4z, &d3p4z, 16);
		double b1b4x[2];
		o.two_Diff(b1x, b4x, b1b4x);
		double b1b4y[2];
		o.two_Diff(b1y, b4y, b1b4y);
		double b1b4z[2];
		o.two_Diff(b1z, b4z, b1b4z);
		double b2b4x[2];
		o.two_Diff(b2x, b4x, b2b4x);
		double b2b4y[2];
		o.two_Diff(b2y, b4y, b2b4y);
		double b2b4z[2];
		o.two_Diff(b2z, b4z, b2b4z);
		double b3b4x[2];
		o.two_Diff(b3x, b4x, b3b4x);
		double b3b4y[2];
		o.two_Diff(b3y, b4y, b3b4y);
		double b3b4z[2];
		o.two_Diff(b3z, b4z, b3b4z);
		double d1_b1b4x_p[16], *d1_b1b4x = d1_b1b4x_p;
		int    d1_b1b4x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b4x, &d1_b1b4x, 16);
		double d1_b1b4y_p[16], *d1_b1b4y = d1_b1b4y_p;
		int    d1_b1b4y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b4y, &d1_b1b4y, 16);
		double d1_b1b4z_p[16], *d1_b1b4z = d1_b1b4z_p;
		int    d1_b1b4z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1b4z, &d1_b1b4z, 16);
		double d2_b2b4x_p[16], *d2_b2b4x = d2_b2b4x_p;
		int    d2_b2b4x_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b4x, &d2_b2b4x, 16);
		double d2_b2b4y_p[16], *d2_b2b4y = d2_b2b4y_p;
		int    d2_b2b4y_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b4y, &d2_b2b4y, 16);
		double d2_b2b4z_p[16], *d2_b2b4z = d2_b2b4z_p;
		int    d2_b2b4z_len =
		  o.Gen_Product_With_PreAlloc(d2_len, d2, 2, b2b4z, &d2_b2b4z, 16);
		double d3_b3b4x_p[16], *d3_b3b4x = d3_b3b4x_p;
		int    d3_b3b4x_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, 2, b3b4x, &d3_b3b4x, 16);
		double d3_b3b4y_p[16], *d3_b3b4y = d3_b3b4y_p;
		int    d3_b3b4y_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, 2, b3b4y, &d3_b3b4y, 16);
		double d3_b3b4z_p[16], *d3_b3b4z = d3_b3b4z_p;
		int    d3_b3b4z_len =
		  o.Gen_Product_With_PreAlloc(d3_len, d3, 2, b3b4z, &d3_b3b4z, 16);
		double d1_b1b4_l1x_p[16], *d1_b1b4_l1x = d1_b1b4_l1x_p;
		int    d1_b1b4_l1x_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b4x_len, d1_b1b4x, l1x_len, l1x, &d1_b1b4_l1x, 16);
		double d1_b1b4_l1y_p[16], *d1_b1b4_l1y = d1_b1b4_l1y_p;
		int    d1_b1b4_l1y_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b4y_len, d1_b1b4y, l1y_len, l1y, &d1_b1b4_l1y, 16);
		double d1_b1b4_l1z_p[16], *d1_b1b4_l1z = d1_b1b4_l1z_p;
		int    d1_b1b4_l1z_len = o.Gen_Sum_With_PreAlloc(
      d1_b1b4z_len, d1_b1b4z, l1z_len, l1z, &d1_b1b4_l1z, 16);
		double d2_b2b4_l2x_p[16], *d2_b2b4_l2x = d2_b2b4_l2x_p;
		int    d2_b2b4_l2x_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b4x_len, d2_b2b4x, l2x_len, l2x, &d2_b2b4_l2x, 16);
		double d2_b2b4_l2y_p[16], *d2_b2b4_l2y = d2_b2b4_l2y_p;
		int    d2_b2b4_l2y_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b4y_len, d2_b2b4y, l2y_len, l2y, &d2_b2b4_l2y, 16);
		double d2_b2b4_l2z_p[16], *d2_b2b4_l2z = d2_b2b4_l2z_p;
		int    d2_b2b4_l2z_len = o.Gen_Sum_With_PreAlloc(
      d2_b2b4z_len, d2_b2b4z, l2z_len, l2z, &d2_b2b4_l2z, 16);
		double d3_b3b4_l3x_p[16], *d3_b3b4_l3x = d3_b3b4_l3x_p;
		int    d3_b3b4_l3x_len = o.Gen_Sum_With_PreAlloc(
      d3_b3b4x_len, d3_b3b4x, l3x_len, l3x, &d3_b3b4_l3x, 16);
		double d3_b3b4_l3y_p[16], *d3_b3b4_l3y = d3_b3b4_l3y_p;
		int    d3_b3b4_l3y_len = o.Gen_Sum_With_PreAlloc(
      d3_b3b4y_len, d3_b3b4y, l3y_len, l3y, &d3_b3b4_l3y, 16);
		double d3_b3b4_l3z_p[16], *d3_b3b4_l3z = d3_b3b4_l3z_p;
		int    d3_b3b4_l3z_len = o.Gen_Sum_With_PreAlloc(
      d3_b3b4z_len, d3_b3b4z, l3z_len, l3z, &d3_b3b4_l3z, 16);
		double d4d1_b1b4_l1x_p[16], *d4d1_b1b4_l1x = d4d1_b1b4_l1x_p;
		int    d4d1_b1b4_l1x_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d1_b1b4_l1x_len, d1_b1b4_l1x, &d4d1_b1b4_l1x, 16);
		double d4d1_b1b4_l1y_p[16], *d4d1_b1b4_l1y = d4d1_b1b4_l1y_p;
		int    d4d1_b1b4_l1y_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d1_b1b4_l1y_len, d1_b1b4_l1y, &d4d1_b1b4_l1y, 16);
		double d4d1_b1b4_l1z_p[16], *d4d1_b1b4_l1z = d4d1_b1b4_l1z_p;
		int    d4d1_b1b4_l1z_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d1_b1b4_l1z_len, d1_b1b4_l1z, &d4d1_b1b4_l1z, 16);
		double d4d2_b2b4_l2x_p[16], *d4d2_b2b4_l2x = d4d2_b2b4_l2x_p;
		int    d4d2_b2b4_l2x_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d2_b2b4_l2x_len, d2_b2b4_l2x, &d4d2_b2b4_l2x, 16);
		double d4d2_b2b4_l2y_p[16], *d4d2_b2b4_l2y = d4d2_b2b4_l2y_p;
		int    d4d2_b2b4_l2y_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d2_b2b4_l2y_len, d2_b2b4_l2y, &d4d2_b2b4_l2y, 16);
		double d4d2_b2b4_l2z_p[16], *d4d2_b2b4_l2z = d4d2_b2b4_l2z_p;
		int    d4d2_b2b4_l2z_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d2_b2b4_l2z_len, d2_b2b4_l2z, &d4d2_b2b4_l2z, 16);
		double d4d3_b3b4_l3x_p[16], *d4d3_b3b4_l3x = d4d3_b3b4_l3x_p;
		int    d4d3_b3b4_l3x_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d3_b3b4_l3x_len, d3_b3b4_l3x, &d4d3_b3b4_l3x, 16);
		double d4d3_b3b4_l3y_p[16], *d4d3_b3b4_l3y = d4d3_b3b4_l3y_p;
		int    d4d3_b3b4_l3y_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d3_b3b4_l3y_len, d3_b3b4_l3y, &d4d3_b3b4_l3y, 16);
		double d4d3_b3b4_l3z_p[16], *d4d3_b3b4_l3z = d4d3_b3b4_l3z_p;
		int    d4d3_b3b4_l3z_len = o.Gen_Product_With_PreAlloc(
      d4_len, d4, d3_b3b4_l3z_len, d3_b3b4_l3z, &d4d3_b3b4_l3z, 16);
		double p1p4x_p[16], *p1p4x = p1p4x_p;
		int p1p4x_len = o.Gen_Diff_With_PreAlloc(d4d1_b1b4_l1x_len, d4d1_b1b4_l1x,
		                                         d1p4x_len, d1p4x, &p1p4x, 16);
		double p1p4y_p[16], *p1p4y = p1p4y_p;
		int p1p4y_len = o.Gen_Diff_With_PreAlloc(d4d1_b1b4_l1y_len, d4d1_b1b4_l1y,
		                                         d1p4y_len, d1p4y, &p1p4y, 16);
		double p1p4z_p[16], *p1p4z = p1p4z_p;
		int p1p4z_len = o.Gen_Diff_With_PreAlloc(d4d1_b1b4_l1z_len, d4d1_b1b4_l1z,
		                                         d1p4z_len, d1p4z, &p1p4z, 16);
		double p2p4x_p[16], *p2p4x = p2p4x_p;
		int p2p4x_len = o.Gen_Diff_With_PreAlloc(d4d2_b2b4_l2x_len, d4d2_b2b4_l2x,
		                                         d2p4x_len, d2p4x, &p2p4x, 16);
		double p2p4y_p[16], *p2p4y = p2p4y_p;
		int p2p4y_len = o.Gen_Diff_With_PreAlloc(d4d2_b2b4_l2y_len, d4d2_b2b4_l2y,
		                                         d2p4y_len, d2p4y, &p2p4y, 16);
		double p2p4z_p[16], *p2p4z = p2p4z_p;
		int p2p4z_len = o.Gen_Diff_With_PreAlloc(d4d2_b2b4_l2z_len, d4d2_b2b4_l2z,
		                                         d2p4z_len, d2p4z, &p2p4z, 16);
		double p3p4x_p[16], *p3p4x = p3p4x_p;
		int p3p4x_len = o.Gen_Diff_With_PreAlloc(d4d3_b3b4_l3x_len, d4d3_b3b4_l3x,
		                                         d3p4x_len, d3p4x, &p3p4x, 16);
		double p3p4y_p[16], *p3p4y = p3p4y_p;
		int p3p4y_len = o.Gen_Diff_With_PreAlloc(d4d3_b3b4_l3y_len, d4d3_b3b4_l3y,
		                                         d3p4y_len, d3p4y, &p3p4y, 16);
		double p3p4z_p[16], *p3p4z = p3p4z_p;
		int p3p4z_len = o.Gen_Diff_With_PreAlloc(d4d3_b3b4_l3z_len, d4d3_b3b4_l3z,
		                                         d3p4z_len, d3p4z, &p3p4z, 16);
		double tmc_a_p[16], *tmc_a = tmc_a_p;
		int    tmc_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4y_len,
		                                               p2p4y, &tmc_a, 16);
		double tmc_b_p[16], *tmc_b = tmc_b_p;
		int    tmc_b_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4x_len,
		                                               p2p4x, &tmc_b, 16);
		double m01_p[16], *m01 = m01_p;
		int    m01_len =
		  o.Gen_Diff_With_PreAlloc(tmc_a_len, tmc_a, tmc_b_len, tmc_b, &m01, 16);
		double tmi_a_p[16], *tmi_a = tmi_a_p;
		int    tmi_a_len = o.Gen_Product_With_PreAlloc(p1p4x_len, p1p4x, p2p4z_len,
		                                               p2p4z, &tmi_a, 16);
		double tmi_b_p[16], *tmi_b = tmi_b_p;
		int    tmi_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4x_len,
		                                               p2p4x, &tmi_b, 16);
		double m02_p[16], *m02 = m02_p;
		int    m02_len =
		  o.Gen_Diff_With_PreAlloc(tmi_a_len, tmi_a, tmi_b_len, tmi_b, &m02, 16);
		double tma_a_p[16], *tma_a = tma_a_p;
		int    tma_a_len = o.Gen_Product_With_PreAlloc(p1p4y_len, p1p4y, p2p4z_len,
		                                               p2p4z, &tma_a, 16);
		double tma_b_p[16], *tma_b = tma_b_p;
		int    tma_b_len = o.Gen_Product_With_PreAlloc(p1p4z_len, p1p4z, p2p4y_len,
		                                               p2p4y, &tma_b, 16);
		double m12_p[16], *m12 = m12_p;
		int    m12_len =
		  o.Gen_Diff_With_PreAlloc(tma_a_len, tma_a, tma_b_len, tma_b, &m12, 16);
		double mt1_p[16], *mt1 = mt1_p;
		int    mt1_len =
		  o.Gen_Product_With_PreAlloc(m01_len, m01, p3p4z_len, p3p4z, &mt1, 16);
		double mt2_p[16], *mt2 = mt2_p;
		int    mt2_len =
		  o.Gen_Product_With_PreAlloc(m02_len, m02, p3p4y_len, p3p4y, &mt2, 16);
		double mt3_p[16], *mt3 = mt3_p;
		int    mt3_len =
		  o.Gen_Product_With_PreAlloc(m12_len, m12, p3p4x_len, p3p4x, &mt3, 16);
		double mtt_p[16], *mtt = mtt_p;
		int    mtt_len =
		  o.Gen_Diff_With_PreAlloc(mt2_len, mt2, mt1_len, mt1, &mtt, 16);
		double m012_p[16], *m012 = m012_p;
		int    m012_len =
		  o.Gen_Diff_With_PreAlloc(mtt_len, mtt, mt3_len, mt3, &m012, 16);

		return_value = m012[m012_len - 1];
		if (m012_p != m012)
			FreeDoubles(m012);
		if (mtt_p != mtt)
			FreeDoubles(mtt);
		if (mt3_p != mt3)
			FreeDoubles(mt3);
		if (mt2_p != mt2)
			FreeDoubles(mt2);
		if (mt1_p != mt1)
			FreeDoubles(mt1);
		if (m12_p != m12)
			FreeDoubles(m12);
		if (tma_b_p != tma_b)
			FreeDoubles(tma_b);
		if (tma_a_p != tma_a)
			FreeDoubles(tma_a);
		if (m02_p != m02)
			FreeDoubles(m02);
		if (tmi_b_p != tmi_b)
			FreeDoubles(tmi_b);
		if (tmi_a_p != tmi_a)
			FreeDoubles(tmi_a);
		if (m01_p != m01)
			FreeDoubles(m01);
		if (tmc_b_p != tmc_b)
			FreeDoubles(tmc_b);
		if (tmc_a_p != tmc_a)
			FreeDoubles(tmc_a);
		if (p3p4z_p != p3p4z)
			FreeDoubles(p3p4z);
		if (p3p4y_p != p3p4y)
			FreeDoubles(p3p4y);
		if (p3p4x_p != p3p4x)
			FreeDoubles(p3p4x);
		if (p2p4z_p != p2p4z)
			FreeDoubles(p2p4z);
		if (p2p4y_p != p2p4y)
			FreeDoubles(p2p4y);
		if (p2p4x_p != p2p4x)
			FreeDoubles(p2p4x);
		if (p1p4z_p != p1p4z)
			FreeDoubles(p1p4z);
		if (p1p4y_p != p1p4y)
			FreeDoubles(p1p4y);
		if (p1p4x_p != p1p4x)
			FreeDoubles(p1p4x);
		if (d4d3_b3b4_l3z_p != d4d3_b3b4_l3z)
			FreeDoubles(d4d3_b3b4_l3z);
		if (d4d3_b3b4_l3y_p != d4d3_b3b4_l3y)
			FreeDoubles(d4d3_b3b4_l3y);
		if (d4d3_b3b4_l3x_p != d4d3_b3b4_l3x)
			FreeDoubles(d4d3_b3b4_l3x);
		if (d4d2_b2b4_l2z_p != d4d2_b2b4_l2z)
			FreeDoubles(d4d2_b2b4_l2z);
		if (d4d2_b2b4_l2y_p != d4d2_b2b4_l2y)
			FreeDoubles(d4d2_b2b4_l2y);
		if (d4d2_b2b4_l2x_p != d4d2_b2b4_l2x)
			FreeDoubles(d4d2_b2b4_l2x);
		if (d4d1_b1b4_l1z_p != d4d1_b1b4_l1z)
			FreeDoubles(d4d1_b1b4_l1z);
		if (d4d1_b1b4_l1y_p != d4d1_b1b4_l1y)
			FreeDoubles(d4d1_b1b4_l1y);
		if (d4d1_b1b4_l1x_p != d4d1_b1b4_l1x)
			FreeDoubles(d4d1_b1b4_l1x);
		if (d3_b3b4_l3z_p != d3_b3b4_l3z)
			FreeDoubles(d3_b3b4_l3z);
		if (d3_b3b4_l3y_p != d3_b3b4_l3y)
			FreeDoubles(d3_b3b4_l3y);
		if (d3_b3b4_l3x_p != d3_b3b4_l3x)
			FreeDoubles(d3_b3b4_l3x);
		if (d2_b2b4_l2z_p != d2_b2b4_l2z)
			FreeDoubles(d2_b2b4_l2z);
		if (d2_b2b4_l2y_p != d2_b2b4_l2y)
			FreeDoubles(d2_b2b4_l2y);
		if (d2_b2b4_l2x_p != d2_b2b4_l2x)
			FreeDoubles(d2_b2b4_l2x);
		if (d1_b1b4_l1z_p != d1_b1b4_l1z)
			FreeDoubles(d1_b1b4_l1z);
		if (d1_b1b4_l1y_p != d1_b1b4_l1y)
			FreeDoubles(d1_b1b4_l1y);
		if (d1_b1b4_l1x_p != d1_b1b4_l1x)
			FreeDoubles(d1_b1b4_l1x);
		if (d3_b3b4z_p != d3_b3b4z)
			FreeDoubles(d3_b3b4z);
		if (d3_b3b4y_p != d3_b3b4y)
			FreeDoubles(d3_b3b4y);
		if (d3_b3b4x_p != d3_b3b4x)
			FreeDoubles(d3_b3b4x);
		if (d2_b2b4z_p != d2_b2b4z)
			FreeDoubles(d2_b2b4z);
		if (d2_b2b4y_p != d2_b2b4y)
			FreeDoubles(d2_b2b4y);
		if (d2_b2b4x_p != d2_b2b4x)
			FreeDoubles(d2_b2b4x);
		if (d1_b1b4z_p != d1_b1b4z)
			FreeDoubles(d1_b1b4z);
		if (d1_b1b4y_p != d1_b1b4y)
			FreeDoubles(d1_b1b4y);
		if (d1_b1b4x_p != d1_b1b4x)
			FreeDoubles(d1_b1b4x);
		if (d3p4z_p != d3p4z)
			FreeDoubles(d3p4z);
		if (d3p4y_p != d3p4y)
			FreeDoubles(d3p4y);
		if (d3p4x_p != d3p4x)
			FreeDoubles(d3p4x);
		if (d2p4z_p != d2p4z)
			FreeDoubles(d2p4z);
		if (d2p4y_p != d2p4y)
			FreeDoubles(d2p4y);
		if (d2p4x_p != d2p4x)
			FreeDoubles(d2p4x);
		if (d1p4z_p != d1p4z)
			FreeDoubles(d1p4z);
		if (d1p4y_p != d1p4y)
			FreeDoubles(d1p4y);
		if (d1p4x_p != d1p4x)
			FreeDoubles(d1p4x);
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
		if (l4x_p != l4x)
			FreeDoubles(l4x);
		if (l4y_p != l4y)
			FreeDoubles(l4y);
		if (l4z_p != l4z)
			FreeDoubles(l4z);
		if (d4_p != d4)
			FreeDoubles(d4);
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orient3D_IIII_exact<IT, ET>(p1, p2, p3, p4);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIII(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr)
{
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orient3D_IIII_filtered<IT, ET>(p1, p2, p3, p4, arr);
		if (is_sign_reliable(ret))
			return ret;
	}
	ret = orient3D_IIII_interval<IT, ET>(p1, p2, p3, p4);
	if (is_sign_reliable(ret))
		return ret;
	return orient3D_IIII_expansion<IT, ET>(p1, p2, p3, p4);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2x,
                               double p2y, double p3x, double p3y, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return Sign::UNCERTAIN;

	double b1p3x    = b1x - p3x;
	double b1p3y    = b1y - p3y;
	double d1_b1p3x = d1 * b1p3x;
	double d1_b1p3y = d1 * b1p3y;
	double ix       = d1_b1p3x + l1x;
	double iy       = d1_b1p3y + l1y;
	double p2p3x    = p2x - p3x;
	double p2p3y    = p2y - p3y;
	double t0       = ix * p2p3y;
	double t1       = iy * p2p3x;
	double det      = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2p3y)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 8.881784197001255e-15;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.197724203485299e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.0409451078885486e-12;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2x, IT p2y,
                               IT p3x, IT p3y)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p3x    = b1x - p3x;
	IT b1p3y    = b1y - p3y;
	IT d1_b1p3x = d1 * b1p3x;
	IT d1_b1p3y = d1 * b1p3y;
	IT ix       = d1_b1p3x + l1x;
	IT iy       = d1_b1p3y + l1y;
	IT p2p3x    = p2x - p3x;
	IT p2p3y    = p2y - p3y;
	IT t0       = ix * p2p3y;
	IT t1       = iy * p2p3x;
	IT det      = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2x, ET p2y,
                            ET p3x, ET p3y)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	ET b1p3x    = b1x - p3x;
	ET b1p3y    = b1y - p3y;
	ET d1_b1p3x = d1 * b1p3x;
	ET d1_b1p3y = d1 * b1p3y;
	ET ix       = d1_b1p3x + l1x;
	ET iy       = d1_b1p3y + l1y;
	ET p2p3x    = p2x - p3x;
	ET p2p3y    = p2y - p3y;
	ET t0       = ix * p2p3y;
	ET t1       = iy * p2p3x;
	ET det      = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2x,
                                double p2y, double p3x, double p3y)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          b1p3x[2];
		o.two_Diff(b1x, p3x, b1p3x);
		double b1p3y[2];
		o.two_Diff(b1y, p3y, b1p3y);
		double d1_b1p3x_p[128], *d1_b1p3x = d1_b1p3x_p;
		int    d1_b1p3x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3x, &d1_b1p3x, 128);
		double d1_b1p3y_p[128], *d1_b1p3y = d1_b1p3y_p;
		int    d1_b1p3y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3y, &d1_b1p3y, 128);
		double ix_p[128], *ix = ix_p;
		int    ix_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3x_len, d1_b1p3x, l1x_len, l1x, &ix, 128);
		double iy_p[128], *iy = iy_p;
		int    iy_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3y_len, d1_b1p3y, l1y_len, l1y, &iy, 128);
		double p2p3x[2];
		o.two_Diff(p2x, p3x, p2p3x);
		double p2p3y[2];
		o.two_Diff(p2y, p3y, p2p3y);
		double t0_p[128], *t0 = t0_p;
		int    t0_len = o.Gen_Product_With_PreAlloc(ix_len, ix, 2, p2p3y, &t0, 128);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(iy_len, iy, 2, p2p3x, &t1, 128);
		double det_p[128], *det = det_p;
		int det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 128);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (iy_p != iy)
			FreeDoubles(iy);
		if (ix_p != ix)
			FreeDoubles(ix);
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
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dxy_IEE_exact<IT, ET>(p1, p2x, p2y, p3x, p3y);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IEE(const GenericPoint3T<IT, ET> &p1, double p2x, double p2y,
                      double p3x, double p3y, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dxy_IEE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dxy_IEE_filtered<IT, ET>(p1, p2x, p2y, p3x, p3y, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dxy_IEE, arr);
	}
	ret = orientOn2Dxy_IEE_interval<IT, ET>(p1, p2x, p2y, p3x, p3y);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dxy_IEE, arr);
	return orientOn2Dxy_IEE_expansion<IT, ET>(p1, p2x, p2y, p3x, p3y);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dxy_IEE<IT, ET, WithSSFilter>(p1, p2.x(), p2.y(), p3.x(),
	                                              p3.y(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double p3x,
                               double p3y, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return Sign::UNCERTAIN;

	double b1p3x    = b1x - p3x;
	double b1p3y    = b1y - p3y;
	double b2p3x    = b2x - p3x;
	double b2p3y    = b2y - p3y;
	double d1_b1p3x = d1 * b1p3x;
	double d1_b1p3y = d1 * b1p3y;
	double i1x      = d1_b1p3x + l1x;
	double i1y      = d1_b1p3y + l1y;
	double d2_b2p3x = d2 * b2p3x;
	double d2_b2p3y = d2 * b2p3y;
	double i2x      = d2_b2p3x + l2x;
	double i2y      = d2_b2p3y + l2y;
	double t0       = i1x * i2y;
	double t1       = i1y * i2x;
	double det      = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p3y)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.684341886080809e-14;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.918578143578212e-13;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.314859663485564e-12;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 6.750832531876594e-13;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7707333863081813e-11;
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
		epsilon *= max_var;
		epsilon *= 4.189644187135872e-10;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT p3x, IT p3y)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p3x    = b1x - p3x;
	IT b1p3y    = b1y - p3y;
	IT b2p3x    = b2x - p3x;
	IT b2p3y    = b2y - p3y;
	IT d1_b1p3x = d1 * b1p3x;
	IT d1_b1p3y = d1 * b1p3y;
	IT i1x      = d1_b1p3x + l1x;
	IT i1y      = d1_b1p3y + l1y;
	IT d2_b2p3x = d2 * b2p3x;
	IT d2_b2p3y = d2 * b2p3y;
	IT i2x      = d2_b2p3x + l2x;
	IT i2y      = d2_b2p3y + l2y;
	IT t0       = i1x * i2y;
	IT t1       = i1y * i2x;
	IT det      = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET p3x, ET p3y)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	ET b1p3x    = b1x - p3x;
	ET b1p3y    = b1y - p3y;
	ET b2p3x    = b2x - p3x;
	ET b2p3y    = b2y - p3y;
	ET d1_b1p3x = d1 * b1p3x;
	ET d1_b1p3y = d1 * b1p3y;
	ET i1x      = d1_b1p3x + l1x;
	ET i1y      = d1_b1p3y + l1y;
	ET d2_b2p3x = d2 * b2p3x;
	ET d2_b2p3y = d2 * b2p3y;
	ET i2x      = d2_b2p3x + l2x;
	ET i2y      = d2_b2p3y + l2y;
	ET t0       = i1x * i2y;
	ET t1       = i1y * i2x;
	ET det      = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double p3x, double p3y,
                      PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dxy_IIE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dxy_IIE_filtered<IT, ET>(p1, p2, p3x, p3y, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dxy_IIE, arr);
	}
	ret = orientOn2Dxy_IIE_interval<IT, ET>(p1, p2, p3x, p3y);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dxy_IIE, arr);
	return orientOn2Dxy_IIE_expansion<IT, ET>(p1, p2, p3x, p3y);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dxy_IIE<IT, ET, WithSSFilter>(p1, p2, p3.x(), p3.y(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var))
		return Sign::UNCERTAIN;

	double b1b3x         = b1x - b3x;
	double b1b3y         = b1y - b3y;
	double b2b3x         = b2x - b3x;
	double b2b3y         = b2y - b3y;
	double d1_b1b3x      = d1 * b1b3x;
	double d1_b1b3y      = d1 * b1b3y;
	double d1_b1b3_l1x   = d1_b1b3x + l1x;
	double d1_b1b3_l1y   = d1_b1b3y + l1y;
	double d3d1_b1b3_l1x = d1_b1b3_l1x * d3;
	double d3d1_b1b3_l1y = d1_b1b3_l1y * d3;
	double d2_b2b3x      = d2 * b2b3x;
	double d2_b2b3y      = d2 * b2b3y;
	double d2_b2b3_l2x   = d2_b2b3x + l2x;
	double d2_b2b3_l2y   = d2_b2b3y + l2y;
	double d3d2_b2b3_l2x = d2_b2b3_l2x * d3;
	double d3d2_b2b3_l2y = d2_b2b3_l2y * d3;
	double l3d1x         = l3x * d1;
	double l3d1y         = l3y * d1;
	double l3d2x         = l3x * d2;
	double l3d2y         = l3y * d2;
	double i1x           = d3d1_b1b3_l1x - l3d1x;
	double i1y           = d3d1_b1b3_l1y - l3d1y;
	double i2x           = d3d2_b2b3_l2x - l3d2x;
	double i2y           = d3d2_b2b3_l2y - l3d2y;
	double t0            = i1x * i2y;
	double t1            = i1y * i2x;
	double det           = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b3y)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 8.455458555545215e-13;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 9.344448825832107e-12;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.816239768639228e-09;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.000987305878955e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1864709104081397e-08;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0999050320824756e-07;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 9.823293567468065e-11;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.742417331587001e-08;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 6.576550303227838e-07;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1310143236187382e-05;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1b3x         = b1x - b3x;
	IT b1b3y         = b1y - b3y;
	IT b2b3x         = b2x - b3x;
	IT b2b3y         = b2y - b3y;
	IT d1_b1b3x      = d1 * b1b3x;
	IT d1_b1b3y      = d1 * b1b3y;
	IT d1_b1b3_l1x   = d1_b1b3x + l1x;
	IT d1_b1b3_l1y   = d1_b1b3y + l1y;
	IT d3d1_b1b3_l1x = d1_b1b3_l1x * d3;
	IT d3d1_b1b3_l1y = d1_b1b3_l1y * d3;
	IT d2_b2b3x      = d2 * b2b3x;
	IT d2_b2b3y      = d2 * b2b3y;
	IT d2_b2b3_l2x   = d2_b2b3x + l2x;
	IT d2_b2b3_l2y   = d2_b2b3y + l2y;
	IT d3d2_b2b3_l2x = d2_b2b3_l2x * d3;
	IT d3d2_b2b3_l2y = d2_b2b3_l2y * d3;
	IT l3d1x         = l3x * d1;
	IT l3d1y         = l3y * d1;
	IT l3d2x         = l3x * d2;
	IT l3d2y         = l3y * d2;
	IT i1x           = d3d1_b1b3_l1x - l3d1x;
	IT i1y           = d3d1_b1b3_l1y - l3d1y;
	IT i2x           = d3d2_b2b3_l2x - l3d2x;
	IT i2y           = d3d2_b2b3_l2y - l3d2y;
	IT t0            = i1x * i2y;
	IT t1            = i1y * i2x;
	IT det           = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dxy_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	p3.getExactLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z);
	ET b1b3x         = b1x - b3x;
	ET b1b3y         = b1y - b3y;
	ET b2b3x         = b2x - b3x;
	ET b2b3y         = b2y - b3y;
	ET d1_b1b3x      = d1 * b1b3x;
	ET d1_b1b3y      = d1 * b1b3y;
	ET d1_b1b3_l1x   = d1_b1b3x + l1x;
	ET d1_b1b3_l1y   = d1_b1b3y + l1y;
	ET d3d1_b1b3_l1x = d1_b1b3_l1x * d3;
	ET d3d1_b1b3_l1y = d1_b1b3_l1y * d3;
	ET d2_b2b3x      = d2 * b2b3x;
	ET d2_b2b3y      = d2 * b2b3y;
	ET d2_b2b3_l2x   = d2_b2b3x + l2x;
	ET d2_b2b3_l2y   = d2_b2b3y + l2y;
	ET d3d2_b2b3_l2x = d2_b2b3_l2x * d3;
	ET d3d2_b2b3_l2y = d2_b2b3_l2y * d3;
	ET l3d1x         = l3x * d1;
	ET l3d1y         = l3y * d1;
	ET l3d2x         = l3x * d2;
	ET l3d2y         = l3y * d2;
	ET i1x           = d3d1_b1b3_l1x - l3d1x;
	ET i1y           = d3d1_b1b3_l1y - l3d1y;
	ET i2x           = d3d2_b2b3_l2x - l3d2x;
	ET i2y           = d3d2_b2b3_l2y - l3d2y;
	ET t0            = i1x * i2y;
	ET t1            = i1y * i2x;
	ET det           = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dxy_III, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dxy_III_filtered<IT, ET>(p1, p2, p3, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dxy_III, arr);
	}
	ret = orientOn2Dxy_III_interval<IT, ET>(p1, p2, p3);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dxy_III, arr);
	return orientOn2Dxy_III_expansion<IT, ET>(p1, p2, p3);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2y,
                               double p2z, double p3y, double p3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return Sign::UNCERTAIN;

	double b1p3y    = b1y - p3y;
	double b1p3z    = b1z - p3z;
	double d1_b1p3y = d1 * b1p3y;
	double d1_b1p3z = d1 * b1p3z;
	double iy       = d1_b1p3y + l1y;
	double iz       = d1_b1p3z + l1z;
	double p2p3y    = p2y - p3y;
	double p2p3z    = p2z - p3z;
	double t0       = iy * p2p3z;
	double t1       = iz * p2p3y;
	double det      = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2p3z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 8.881784197001255e-15;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.197724203485299e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.0409451078885486e-12;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2y, IT p2z,
                               IT p3y, IT p3z)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p3y    = b1y - p3y;
	IT b1p3z    = b1z - p3z;
	IT d1_b1p3y = d1 * b1p3y;
	IT d1_b1p3z = d1 * b1p3z;
	IT iy       = d1_b1p3y + l1y;
	IT iz       = d1_b1p3z + l1z;
	IT p2p3y    = p2y - p3y;
	IT p2p3z    = p2z - p3z;
	IT t0       = iy * p2p3z;
	IT t1       = iz * p2p3y;
	IT det      = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2y, ET p2z,
                            ET p3y, ET p3z)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	ET b1p3y    = b1y - p3y;
	ET b1p3z    = b1z - p3z;
	ET d1_b1p3y = d1 * b1p3y;
	ET d1_b1p3z = d1 * b1p3z;
	ET iy       = d1_b1p3y + l1y;
	ET iz       = d1_b1p3z + l1z;
	ET p2p3y    = p2y - p3y;
	ET p2p3z    = p2z - p3z;
	ET t0       = iy * p2p3z;
	ET t1       = iz * p2p3y;
	ET det      = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2y,
                                double p2z, double p3y, double p3z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          b1p3y[2];
		o.two_Diff(b1y, p3y, b1p3y);
		double b1p3z[2];
		o.two_Diff(b1z, p3z, b1p3z);
		double d1_b1p3y_p[128], *d1_b1p3y = d1_b1p3y_p;
		int    d1_b1p3y_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3y, &d1_b1p3y, 128);
		double d1_b1p3z_p[128], *d1_b1p3z = d1_b1p3z_p;
		int    d1_b1p3z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3z, &d1_b1p3z, 128);
		double iy_p[128], *iy = iy_p;
		int    iy_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3y_len, d1_b1p3y, l1y_len, l1y, &iy, 128);
		double iz_p[128], *iz = iz_p;
		int    iz_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3z_len, d1_b1p3z, l1z_len, l1z, &iz, 128);
		double p2p3y[2];
		o.two_Diff(p2y, p3y, p2p3y);
		double p2p3z[2];
		o.two_Diff(p2z, p3z, p2p3z);
		double t0_p[128], *t0 = t0_p;
		int    t0_len = o.Gen_Product_With_PreAlloc(iy_len, iy, 2, p2p3z, &t0, 128);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(iz_len, iz, 2, p2p3y, &t1, 128);
		double det_p[128], *det = det_p;
		int det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 128);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (iz_p != iz)
			FreeDoubles(iz);
		if (iy_p != iy)
			FreeDoubles(iy);
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
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dyz_IEE_exact<IT, ET>(p1, p2y, p2z, p3y, p3z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IEE(const GenericPoint3T<IT, ET> &p1, double p2y, double p2z,
                      double p3y, double p3z, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dyz_IEE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dyz_IEE_filtered<IT, ET>(p1, p2y, p2z, p3y, p3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dyz_IEE, arr);
	}
	ret = orientOn2Dyz_IEE_interval<IT, ET>(p1, p2y, p2z, p3y, p3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dyz_IEE, arr);
	return orientOn2Dyz_IEE_expansion<IT, ET>(p1, p2y, p2z, p3y, p3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dyz_IEE<IT, ET, WithSSFilter>(p1, p2.y(), p2.z(), p3.y(),
	                                              p3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double p3y,
                               double p3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return Sign::UNCERTAIN;

	double b1p3y    = b1y - p3y;
	double b1p3z    = b1z - p3z;
	double b2p3y    = b2y - p3y;
	double b2p3z    = b2z - p3z;
	double d1_b1p3y = d1 * b1p3y;
	double d1_b1p3z = d1 * b1p3z;
	double i1y      = d1_b1p3y + l1y;
	double i1z      = d1_b1p3z + l1z;
	double d2_b2p3y = d2 * b2p3y;
	double d2_b2p3z = d2 * b2p3z;
	double i2y      = d2_b2p3y + l2y;
	double i2z      = d2_b2p3z + l2z;
	double t0       = i1y * i2z;
	double t1       = i1z * i2y;
	double det      = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p3z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.684341886080809e-14;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.918578143578212e-13;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.314859663485564e-12;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 6.750832531876594e-13;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7707333863081813e-11;
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
		epsilon *= max_var;
		epsilon *= 4.189644187135872e-10;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT p3y, IT p3z)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p3y    = b1y - p3y;
	IT b1p3z    = b1z - p3z;
	IT b2p3y    = b2y - p3y;
	IT b2p3z    = b2z - p3z;
	IT d1_b1p3y = d1 * b1p3y;
	IT d1_b1p3z = d1 * b1p3z;
	IT i1y      = d1_b1p3y + l1y;
	IT i1z      = d1_b1p3z + l1z;
	IT d2_b2p3y = d2 * b2p3y;
	IT d2_b2p3z = d2 * b2p3z;
	IT i2y      = d2_b2p3y + l2y;
	IT i2z      = d2_b2p3z + l2z;
	IT t0       = i1y * i2z;
	IT t1       = i1z * i2y;
	IT det      = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET p3y, ET p3z)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	ET b1p3y    = b1y - p3y;
	ET b1p3z    = b1z - p3z;
	ET b2p3y    = b2y - p3y;
	ET b2p3z    = b2z - p3z;
	ET d1_b1p3y = d1 * b1p3y;
	ET d1_b1p3z = d1 * b1p3z;
	ET i1y      = d1_b1p3y + l1y;
	ET i1z      = d1_b1p3z + l1z;
	ET d2_b2p3y = d2 * b2p3y;
	ET d2_b2p3z = d2 * b2p3z;
	ET i2y      = d2_b2p3y + l2y;
	ET i2z      = d2_b2p3z + l2z;
	ET t0       = i1y * i2z;
	ET t1       = i1z * i2y;
	ET det      = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double p3y, double p3z,
                      PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dyz_IIE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dyz_IIE_filtered<IT, ET>(p1, p2, p3y, p3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dyz_IIE, arr);
	}
	ret = orientOn2Dyz_IIE_interval<IT, ET>(p1, p2, p3y, p3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dyz_IIE, arr);
	return orientOn2Dyz_IIE_expansion<IT, ET>(p1, p2, p3y, p3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dyz_IIE<IT, ET, WithSSFilter>(p1, p2, p3.y(), p3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var))
		return Sign::UNCERTAIN;

	double b1b3y         = b1y - b3y;
	double b1b3z         = b1z - b3z;
	double b2b3y         = b2y - b3y;
	double b2b3z         = b2z - b3z;
	double d1_b1b3y      = d1 * b1b3y;
	double d1_b1b3z      = d1 * b1b3z;
	double d1_b1b3_l1y   = d1_b1b3y + l1y;
	double d1_b1b3_l1z   = d1_b1b3z + l1z;
	double d3d1_b1b3_l1y = d1_b1b3_l1y * d3;
	double d3d1_b1b3_l1z = d1_b1b3_l1z * d3;
	double d2_b2b3y      = d2 * b2b3y;
	double d2_b2b3z      = d2 * b2b3z;
	double d2_b2b3_l2y   = d2_b2b3y + l2y;
	double d2_b2b3_l2z   = d2_b2b3z + l2z;
	double d3d2_b2b3_l2y = d2_b2b3_l2y * d3;
	double d3d2_b2b3_l2z = d2_b2b3_l2z * d3;
	double l3d1y         = l3y * d1;
	double l3d1z         = l3z * d1;
	double l3d2y         = l3y * d2;
	double l3d2z         = l3z * d2;
	double i1y           = d3d1_b1b3_l1y - l3d1y;
	double i1z           = d3d1_b1b3_l1z - l3d1z;
	double i2y           = d3d2_b2b3_l2y - l3d2y;
	double i2z           = d3d2_b2b3_l2z - l3d2z;
	double t0            = i1y * i2z;
	double t1            = i1z * i2y;
	double det           = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b3y)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b3z)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 8.455458555545215e-13;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 9.344448825832107e-12;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.816239768639228e-09;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.000987305878955e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1864709104081397e-08;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0999050320824756e-07;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 9.823293567468065e-11;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.742417331587001e-08;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 6.576550303227838e-07;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1310143236187382e-05;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1b3y         = b1y - b3y;
	IT b1b3z         = b1z - b3z;
	IT b2b3y         = b2y - b3y;
	IT b2b3z         = b2z - b3z;
	IT d1_b1b3y      = d1 * b1b3y;
	IT d1_b1b3z      = d1 * b1b3z;
	IT d1_b1b3_l1y   = d1_b1b3y + l1y;
	IT d1_b1b3_l1z   = d1_b1b3z + l1z;
	IT d3d1_b1b3_l1y = d1_b1b3_l1y * d3;
	IT d3d1_b1b3_l1z = d1_b1b3_l1z * d3;
	IT d2_b2b3y      = d2 * b2b3y;
	IT d2_b2b3z      = d2 * b2b3z;
	IT d2_b2b3_l2y   = d2_b2b3y + l2y;
	IT d2_b2b3_l2z   = d2_b2b3z + l2z;
	IT d3d2_b2b3_l2y = d2_b2b3_l2y * d3;
	IT d3d2_b2b3_l2z = d2_b2b3_l2z * d3;
	IT l3d1y         = l3y * d1;
	IT l3d1z         = l3z * d1;
	IT l3d2y         = l3y * d2;
	IT l3d2z         = l3z * d2;
	IT i1y           = d3d1_b1b3_l1y - l3d1y;
	IT i1z           = d3d1_b1b3_l1z - l3d1z;
	IT i2y           = d3d2_b2b3_l2y - l3d2y;
	IT i2z           = d3d2_b2b3_l2z - l3d2z;
	IT t0            = i1y * i2z;
	IT t1            = i1z * i2y;
	IT det           = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dyz_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	p3.getExactLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z);
	ET b1b3y         = b1y - b3y;
	ET b1b3z         = b1z - b3z;
	ET b2b3y         = b2y - b3y;
	ET b2b3z         = b2z - b3z;
	ET d1_b1b3y      = d1 * b1b3y;
	ET d1_b1b3z      = d1 * b1b3z;
	ET d1_b1b3_l1y   = d1_b1b3y + l1y;
	ET d1_b1b3_l1z   = d1_b1b3z + l1z;
	ET d3d1_b1b3_l1y = d1_b1b3_l1y * d3;
	ET d3d1_b1b3_l1z = d1_b1b3_l1z * d3;
	ET d2_b2b3y      = d2 * b2b3y;
	ET d2_b2b3z      = d2 * b2b3z;
	ET d2_b2b3_l2y   = d2_b2b3y + l2y;
	ET d2_b2b3_l2z   = d2_b2b3z + l2z;
	ET d3d2_b2b3_l2y = d2_b2b3_l2y * d3;
	ET d3d2_b2b3_l2z = d2_b2b3_l2z * d3;
	ET l3d1y         = l3y * d1;
	ET l3d1z         = l3z * d1;
	ET l3d2y         = l3y * d2;
	ET l3d2z         = l3z * d2;
	ET i1y           = d3d1_b1b3_l1y - l3d1y;
	ET i1z           = d3d1_b1b3_l1z - l3d1z;
	ET i2y           = d3d2_b2b3_l2y - l3d2y;
	ET i2z           = d3d2_b2b3_l2z - l3d2z;
	ET t0            = i1y * i2z;
	ET t1            = i1z * i2y;
	ET det           = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dyz_III, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dyz_III_filtered<IT, ET>(p1, p2, p3, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dyz_III, arr);
	}
	ret = orientOn2Dyz_III_interval<IT, ET>(p1, p2, p3);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dyz_III, arr);
	return orientOn2Dyz_III_expansion<IT, ET>(p1, p2, p3);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2x,
                               double p2z, double p3x, double p3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var))
		return Sign::UNCERTAIN;

	double b1p3z    = b1z - p3z;
	double b1p3x    = b1x - p3x;
	double d1_b1p3z = d1 * b1p3z;
	double d1_b1p3x = d1 * b1p3x;
	double iz       = d1_b1p3z + l1z;
	double ix       = d1_b1p3x + l1x;
	double p2p3z    = p2z - p3z;
	double p2p3x    = p2x - p3x;
	double t0       = iz * p2p3x;
	double t1       = ix * p2p3z;
	double det      = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(p2p3x)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::S:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 8.881784197001255e-15;
	}
	break;
	case PntArr3::L:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= 3.197724203485299e-14;
	}
	break;
	case PntArr3::T:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 1.0409451078885486e-12;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2x, IT p2z,
                               IT p3x, IT p3z)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p3z    = b1z - p3z;
	IT b1p3x    = b1x - p3x;
	IT d1_b1p3z = d1 * b1p3z;
	IT d1_b1p3x = d1 * b1p3x;
	IT iz       = d1_b1p3z + l1z;
	IT ix       = d1_b1p3x + l1x;
	IT p2p3z    = p2z - p3z;
	IT p2p3x    = p2x - p3x;
	IT t0       = iz * p2p3x;
	IT t1       = ix * p2p3z;
	IT det      = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2x, ET p2z,
                            ET p3x, ET p3z)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	ET b1p3z    = b1z - p3z;
	ET b1p3x    = b1x - p3x;
	ET d1_b1p3z = d1 * b1p3z;
	ET d1_b1p3x = d1 * b1p3x;
	ET iz       = d1_b1p3z + l1z;
	ET ix       = d1_b1p3x + l1x;
	ET p2p3z    = p2z - p3z;
	ET p2p3x    = p2x - p3x;
	ET t0       = iz * p2p3x;
	ET t1       = ix * p2p3z;
	ET det      = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2x,
                                double p2z, double p3x, double p3z)
{
	double return_value = NAN;
	#ifdef CHECK_FOR_XYZERFLOWS
	feclearexcept(FE_ALL_EXCEPT);
	#endif
	double l1x_p[128], *l1x = l1x_p, l1y_p[128], *l1y = l1y_p, l1z_p[128],
	                   *l1z = l1z_p, d1_p[128], *d1 = d1_p, b1x, b1y, b1z;
	int l1x_len = 128, l1y_len = 128, l1z_len = 128, d1_len = 128;
	p1.getExpansionLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1,
	                      d1_len, b1x, b1y, b1z);
	if ((d1[d1_len - 1] != 0))
	{
		expansionObject o;
		double          b1p3z[2];
		o.two_Diff(b1z, p3z, b1p3z);
		double b1p3x[2];
		o.two_Diff(b1x, p3x, b1p3x);
		double d1_b1p3z_p[128], *d1_b1p3z = d1_b1p3z_p;
		int    d1_b1p3z_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3z, &d1_b1p3z, 128);
		double d1_b1p3x_p[128], *d1_b1p3x = d1_b1p3x_p;
		int    d1_b1p3x_len =
		  o.Gen_Product_With_PreAlloc(d1_len, d1, 2, b1p3x, &d1_b1p3x, 128);
		double iz_p[128], *iz = iz_p;
		int    iz_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3z_len, d1_b1p3z, l1z_len, l1z, &iz, 128);
		double ix_p[128], *ix = ix_p;
		int    ix_len =
		  o.Gen_Sum_With_PreAlloc(d1_b1p3x_len, d1_b1p3x, l1x_len, l1x, &ix, 128);
		double p2p3z[2];
		o.two_Diff(p2z, p3z, p2p3z);
		double p2p3x[2];
		o.two_Diff(p2x, p3x, p2p3x);
		double t0_p[128], *t0 = t0_p;
		int    t0_len = o.Gen_Product_With_PreAlloc(iz_len, iz, 2, p2p3x, &t0, 128);
		double t1_p[128], *t1 = t1_p;
		int    t1_len = o.Gen_Product_With_PreAlloc(ix_len, ix, 2, p2p3z, &t1, 128);
		double det_p[128], *det = det_p;
		int det_len = o.Gen_Diff_With_PreAlloc(t0_len, t0, t1_len, t1, &det, 128);

		return_value = det[det_len - 1];
		if (det_p != det)
			FreeDoubles(det);
		if (t1_p != t1)
			FreeDoubles(t1);
		if (t0_p != t0)
			FreeDoubles(t0);
		if (ix_p != ix)
			FreeDoubles(ix);
		if (iz_p != iz)
			FreeDoubles(iz);
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
	}

	#ifdef CHECK_FOR_XYZERFLOWS
	if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
		return orientOn2Dzx_IEE_exact<IT, ET>(p1, p2x, p2z, p3x, p3z);
	#endif

	if (return_value > 0)
		return Sign::POSITIVE;
	if (return_value < 0)
		return Sign::NEGATIVE;
	if (return_value == 0)
		return Sign::ZERO;
	OMC_EXIT("Should not happen.");
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IEE(const GenericPoint3T<IT, ET> &p1, double p2x, double p2z,
                      double p3x, double p3z, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dzx_IEE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dzx_IEE_filtered<IT, ET>(p1, p2x, p2z, p3x, p3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dzx_IEE, arr);
	}
	ret = orientOn2Dzx_IEE_interval<IT, ET>(p1, p2x, p2z, p3x, p3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dzx_IEE, arr);
	return orientOn2Dzx_IEE_expansion<IT, ET>(p1, p2x, p2z, p3x, p3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dzx_IEE<IT, ET, WithSSFilter>(p1, p2.x(), p2.z(), p3.x(),
	                                              p3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double p3x,
                               double p3z, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var))
		return Sign::UNCERTAIN;

	double b1p3z    = b1z - p3z;
	double b1p3x    = b1x - p3x;
	double b2p3z    = b2z - p3z;
	double b2p3x    = b2x - p3x;
	double d1_b1p3z = d1 * b1p3z;
	double d1_b1p3x = d1 * b1p3x;
	double i1z      = d1_b1p3z + l1z;
	double i1x      = d1_b1p3x + l1x;
	double d2_b2p3z = d2 * b2p3z;
	double d2_b2p3x = d2 * b2p3x;
	double i2z      = d2_b2p3z + l2z;
	double i2x      = d2_b2p3x + l2x;
	double t0       = i1z * i2x;
	double t1       = i1x * i2z;
	double det      = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1p3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2p3x)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.684341886080809e-14;
	}
	break;
	case PntArr3::SL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.918578143578212e-13;
	}
	break;
	case PntArr3::ST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 5.314859663485564e-12;
	}
	break;
	case PntArr3::LL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= 6.750832531876594e-13;
	}
	break;
	case PntArr3::LT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.7707333863081813e-11;
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
		epsilon *= max_var;
		epsilon *= 4.189644187135872e-10;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT p3x, IT p3z)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1p3z    = b1z - p3z;
	IT b1p3x    = b1x - p3x;
	IT b2p3z    = b2z - p3z;
	IT b2p3x    = b2x - p3x;
	IT d1_b1p3z = d1 * b1p3z;
	IT d1_b1p3x = d1 * b1p3x;
	IT i1z      = d1_b1p3z + l1z;
	IT i1x      = d1_b1p3x + l1x;
	IT d2_b2p3z = d2 * b2p3z;
	IT d2_b2p3x = d2 * b2p3x;
	IT i2z      = d2_b2p3z + l2z;
	IT i2x      = d2_b2p3x + l2x;
	IT t0       = i1z * i2x;
	IT t1       = i1x * i2z;
	IT det      = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET p3x, ET p3z)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	ET b1p3z    = b1z - p3z;
	ET b1p3x    = b1x - p3x;
	ET b2p3z    = b2z - p3z;
	ET b2p3x    = b2x - p3x;
	ET d1_b1p3z = d1 * b1p3z;
	ET d1_b1p3x = d1 * b1p3x;
	ET i1z      = d1_b1p3z + l1z;
	ET i1x      = d1_b1p3x + l1x;
	ET d2_b2p3z = d2 * b2p3z;
	ET d2_b2p3x = d2 * b2p3x;
	ET i2z      = d2_b2p3z + l2z;
	ET i2x      = d2_b2p3x + l2x;
	ET t0       = i1z * i2x;
	ET t1       = i1x * i2z;
	ET det      = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double p3x, double p3z,
                      PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dzx_IIE, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dzx_IIE_filtered<IT, ET>(p1, p2, p3x, p3z, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dzx_IIE, arr);
	}
	ret = orientOn2Dzx_IIE_interval<IT, ET>(p1, p2, p3x, p3z);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dzx_IIE, arr);
	return orientOn2Dzx_IIE_expansion<IT, ET>(p1, p2, p3x, p3z);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	return orientOn2Dzx_IIE<IT, ET, WithSSFilter>(p1, p2, p3.x(), p3.z(), arr);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	double l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z,
	  l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var = 0;
	if (!p1.getFilteredLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z, max_var) ||
	    !p2.getFilteredLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z, max_var) ||
	    !p3.getFilteredLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z, max_var))
		return Sign::UNCERTAIN;

	double b1b3z         = b1z - b3z;
	double b1b3x         = b1x - b3x;
	double b2b3z         = b2z - b3z;
	double b2b3x         = b2x - b3x;
	double d1_b1b3z      = d1 * b1b3z;
	double d1_b1b3x      = d1 * b1b3x;
	double d1_b1b3_l1z   = d1_b1b3z + l1z;
	double d1_b1b3_l1x   = d1_b1b3x + l1x;
	double d3d1_b1b3_l1z = d1_b1b3_l1z * d3;
	double d3d1_b1b3_l1x = d1_b1b3_l1x * d3;
	double d2_b2b3z      = d2 * b2b3z;
	double d2_b2b3x      = d2 * b2b3x;
	double d2_b2b3_l2z   = d2_b2b3z + l2z;
	double d2_b2b3_l2x   = d2_b2b3x + l2x;
	double d3d2_b2b3_l2z = d2_b2b3_l2z * d3;
	double d3d2_b2b3_l2x = d2_b2b3_l2x * d3;
	double l3d1z         = l3z * d1;
	double l3d1x         = l3x * d1;
	double l3d2z         = l3z * d2;
	double l3d2x         = l3x * d2;
	double i1z           = d3d1_b1b3_l1z - l3d1z;
	double i1x           = d3d1_b1b3_l1x - l3d1x;
	double i2z           = d3d2_b2b3_l2z - l3d2z;
	double i2x           = d3d2_b2b3_l2x - l3d2x;
	double t0            = i1z * i2x;
	double t1            = i1x * i2z;
	double det           = t0 - t1;

	double _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b1b3x)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b3z)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(b2b3x)) > max_var)
		max_var = _tmp_fabs;
	double epsilon = max_var;
	switch (arr)
	{
	case PntArr3::SSS:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 8.455458555545215e-13;
	}
	break;
	case PntArr3::SSL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 9.344448825832107e-12;
	}
	break;
	case PntArr3::SST:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.816239768639228e-09;
	}
	break;
	case PntArr3::SLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.000987305878955e-11;
	}
	break;
	case PntArr3::SLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1864709104081397e-08;
	}
	break;
	case PntArr3::STT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 2.0999050320824756e-07;
	}
	break;
	case PntArr3::LLL:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 9.823293567468065e-11;
	}
	break;
	case PntArr3::LLT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 3.742417331587001e-08;
	}
	break;
	case PntArr3::LTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 6.576550303227838e-07;
	}
	break;
	case PntArr3::TTT:
	{
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= epsilon;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= max_var;
		epsilon *= 1.1310143236187382e-05;
	}
	break;
	default:
		OMC_EXIT("Unsopported points arrangement.");
	}

	return filter_sign(det, epsilon);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3)
{
	IT l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	if (!p1.getIntervalLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z) ||
	    !p2.getIntervalLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z) ||
	    !p3.getIntervalLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z))
		return Sign::UNCERTAIN;

	typename IT::Protector P;

	IT b1b3z         = b1z - b3z;
	IT b1b3x         = b1x - b3x;
	IT b2b3z         = b2z - b3z;
	IT b2b3x         = b2x - b3x;
	IT d1_b1b3z      = d1 * b1b3z;
	IT d1_b1b3x      = d1 * b1b3x;
	IT d1_b1b3_l1z   = d1_b1b3z + l1z;
	IT d1_b1b3_l1x   = d1_b1b3x + l1x;
	IT d3d1_b1b3_l1z = d1_b1b3_l1z * d3;
	IT d3d1_b1b3_l1x = d1_b1b3_l1x * d3;
	IT d2_b2b3z      = d2 * b2b3z;
	IT d2_b2b3x      = d2 * b2b3x;
	IT d2_b2b3_l2z   = d2_b2b3z + l2z;
	IT d2_b2b3_l2x   = d2_b2b3x + l2x;
	IT d3d2_b2b3_l2z = d2_b2b3_l2z * d3;
	IT d3d2_b2b3_l2x = d2_b2b3_l2x * d3;
	IT l3d1z         = l3z * d1;
	IT l3d1x         = l3x * d1;
	IT l3d2z         = l3z * d2;
	IT l3d2x         = l3x * d2;
	IT i1z           = d3d1_b1b3_l1z - l3d1z;
	IT i1x           = d3d1_b1b3_l1x - l3d1x;
	IT i2z           = d3d2_b2b3_l2z - l3d2z;
	IT i2x           = d3d2_b2b3_l2x - l3d2x;
	IT t0            = i1z * i2x;
	IT t1            = i1x * i2z;
	IT det           = t0 - t1;
	if (!det.is_sign_reliable())
		return Sign::UNCERTAIN;
	return OMC::sign(det);
}

template <typename IT, typename ET>
Sign orientOn2Dzx_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3)
{
	ET l1x, l1y, l1z, d1, b1x, b1y, b1z, l2x, l2y, l2z, d2, b2x, b2y, b2z, l3x,
	  l3y, l3z, d3, b3x, b3y, b3z;
	p1.getExactLambda(l1x, l1y, l1z, d1, b1x, b1y, b1z);
	p2.getExactLambda(l2x, l2y, l2z, d2, b2x, b2y, b2z);
	p3.getExactLambda(l3x, l3y, l3z, d3, b3x, b3y, b3z);
	ET b1b3z         = b1z - b3z;
	ET b1b3x         = b1x - b3x;
	ET b2b3z         = b2z - b3z;
	ET b2b3x         = b2x - b3x;
	ET d1_b1b3z      = d1 * b1b3z;
	ET d1_b1b3x      = d1 * b1b3x;
	ET d1_b1b3_l1z   = d1_b1b3z + l1z;
	ET d1_b1b3_l1x   = d1_b1b3x + l1x;
	ET d3d1_b1b3_l1z = d1_b1b3_l1z * d3;
	ET d3d1_b1b3_l1x = d1_b1b3_l1x * d3;
	ET d2_b2b3z      = d2 * b2b3z;
	ET d2_b2b3x      = d2 * b2b3x;
	ET d2_b2b3_l2z   = d2_b2b3z + l2z;
	ET d2_b2b3_l2x   = d2_b2b3x + l2x;
	ET d3d2_b2b3_l2z = d2_b2b3_l2z * d3;
	ET d3d2_b2b3_l2x = d2_b2b3_l2x * d3;
	ET l3d1z         = l3z * d1;
	ET l3d1x         = l3x * d1;
	ET l3d2z         = l3z * d2;
	ET l3d2x         = l3x * d2;
	ET i1z           = d3d1_b1b3_l1z - l3d1z;
	ET i1x           = d3d1_b1b3_l1x - l3d1x;
	ET i2z           = d3d2_b2b3_l2z - l3d2z;
	ET i2x           = d3d2_b2b3_l2x - l3d2x;
	ET t0            = i1z * i2x;
	ET t1            = i1x * i2z;
	ET det           = t0 - t1;
	return OMC::sign(det);
}

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr)
{
	OMC_PRED_PROFILE_INC_TOTAL(PredicateNames::_orientOn2Dzx_III, arr);
	Sign ret;
	if constexpr (WithSSFilter)
	{
		ret = orientOn2Dzx_III_filtered<IT, ET>(p1, p2, p3, arr);
		if (is_sign_reliable(ret))
			return ret;
		OMC_PRED_PROFILE_INC_SSFAIL(PredicateNames::_orientOn2Dzx_III, arr);
	}
	ret = orientOn2Dzx_III_interval<IT, ET>(p1, p2, p3);
	if (is_sign_reliable(ret))
		return ret;
	OMC_PRED_PROFILE_INC_DFAIL(PredicateNames::_orientOn2Dzx_III, arr);
	return orientOn2Dzx_III_expansion<IT, ET>(p1, p2, p3);
}

#endif

} // namespace OMC