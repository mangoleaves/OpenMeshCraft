#pragma once

#include "ImplicitPoint3T.h"

#include "OpenMeshCraft/Geometry/Predicates/ImplicitPointPredicates.h"

namespace OMC {

template <typename IT, typename ET>
ImplicitPoint3T_TPI<IT, ET>::ImplicitPoint3T_TPI() noexcept
  : GP(PointType::TPI)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_TPI<IT, ET>::ImplicitPoint3T_TPI(const EP &_v1, const EP &_v2,
                                                 const EP &_v3, const EP &_w1,
                                                 const EP &_w2, const EP &_w3,
                                                 const EP &_u1, const EP &_u2,
                                                 const EP &_u3) noexcept
  : GP(PointType::TPI)
  , iv1(&_v1)
  , iv2(&_v2)
  , iv3(&_v3)
  , iw1(&_w1)
  , iw2(&_w2)
  , iw3(&_w3)
  , iu1(&_u1)
  , iu2(&_u2)
  , iu3(&_u3)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_TPI<IT, ET>::~ImplicitPoint3T_TPI() noexcept
{
}

template <typename IT, typename ET>
ImplicitPoint3T_TPI<IT, ET>::ImplicitPoint3T_TPI(const IP &rhs) noexcept
  : GP(static_cast<const GP &>(rhs))
  , iv1(rhs.iv1)
  , iv2(rhs.iv2)
  , iv3(rhs.iv3)
  , iw1(rhs.iw1)
  , iw2(rhs.iw2)
  , iw3(rhs.iw3)
  , iu1(rhs.iu1)
  , iu2(rhs.iu2)
  , iu3(rhs.iu3)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_TPI<IT, ET>::ImplicitPoint3T_TPI(IP &&rhs) noexcept
  : GP(static_cast<GP &&>(rhs))
  , iv1(rhs.iv1)
  , iv2(rhs.iv2)
  , iv3(rhs.iv3)
  , iw1(rhs.iw1)
  , iw2(rhs.iw2)
  , iw3(rhs.iw3)
  , iu1(rhs.iu1)
  , iu2(rhs.iu2)
  , iu3(rhs.iu3)
{
}

template <typename IT, typename ET>
auto ImplicitPoint3T_TPI<IT, ET>::operator=(const IP &rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<const GP &>(rhs));
	iv1                      = rhs.iv1;
	iv2                      = rhs.iv2;
	iv3                      = rhs.iv3;
	iw1                      = rhs.iw1;
	iw2                      = rhs.iw2;
	iw3                      = rhs.iw3;
	iu1                      = rhs.iu1;
	iu2                      = rhs.iu2;
	iu3                      = rhs.iu3;
	return *this;
}

template <typename IT, typename ET>
auto ImplicitPoint3T_TPI<IT, ET>::operator=(IP &&rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<GP &&>(rhs));
	iv1                      = rhs.iv1;
	iv2                      = rhs.iv2;
	iv3                      = rhs.iv3;
	iw1                      = rhs.iw1;
	iw2                      = rhs.iw2;
	iw3                      = rhs.iw3;
	iu1                      = rhs.iu1;
	iu2                      = rhs.iu2;
	iu3                      = rhs.iu3;
	return *this;
}

template <typename IT, typename ET>
void ImplicitPoint3T_TPI<IT, ET>::get_Explicit(EP &e) const
{
	FT lambda_x, lambda_y, lambda_z, lambda_d, max_var = 0;
	if (!getFilteredLambda(lambda_x, lambda_y, lambda_z, lambda_d, max_var))
	{
		IT ilx, ily, ilz, id;
		if (!getIntervalLambda(ilx, ily, ilz, id))
		{
			ET elx, ely, elz, eld;
			getExactLambda(elx, ely, elz, eld);
			lambda_x = OMC::to_double(elx);
			lambda_y = OMC::to_double(ely);
			lambda_z = OMC::to_double(elz);
			lambda_d = OMC::to_double(eld);
		}
		else
		{
			lambda_x = ilx.sup() + ilx.inf();
			lambda_y = ily.sup() + ily.inf();
			lambda_z = ilz.sup() + ilz.inf();
			lambda_d = id.sup() + id.inf();
		}
	}
	e = EP(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
}

template <typename IT, typename ET>
auto ImplicitPoint3T_TPI<IT, ET>::to_Explicit() const -> EP
{
	EP e;
	get_Explicit(e);
	return e;
}

template <typename IT, typename ET>
bool ImplicitPoint3T_TPI<IT, ET>::getFilteredLambda(FT &lx, FT &ly, FT &lz,
                                                    FT &d, FT &mv) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.ssfilter_cached)
		{
			cv.ssfilter_cached  = true;
			cv.ssfilter_max_val = 0;
			if (!lambda3d_TPI_filtered(
			      V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(),
			      V3().x(), V3().y(), V3().z(), W1().x(), W1().y(), W1().z(),
			      W2().x(), W2().y(), W2().z(), W3().x(), W3().y(), W3().z(),
			      U1().x(), U1().y(), U1().z(), U2().x(), U2().y(), U2().z(),
			      U3().x(), U3().y(), U3().z(), cv.ssfilter_lambda_x,
			      cv.ssfilter_lambda_y, cv.ssfilter_lambda_z, cv.ssfilter_denominator,
			      cv.ssfilter_max_val))
				cv.ssfilter_denominator = 0;

			if (cv.ssfilter_denominator < 0)
			{
				cv.ssfilter_lambda_x    = -cv.ssfilter_lambda_x;
				cv.ssfilter_lambda_y    = -cv.ssfilter_lambda_y;
				cv.ssfilter_lambda_z    = -cv.ssfilter_lambda_z;
				cv.ssfilter_denominator = -cv.ssfilter_denominator;
			}
		}

		lx = cv.ssfilter_lambda_x;
		ly = cv.ssfilter_lambda_y;
		lz = cv.ssfilter_lambda_z;
		d  = cv.ssfilter_denominator;
		if (cv.ssfilter_denominator != 0 && cv.ssfilter_max_val > mv)
			mv = cv.ssfilter_max_val;
		return (cv.ssfilter_denominator != 0);
	}
	else
	{
		FT mv_ = 0;
		if (!lambda3d_TPI_filtered(V1().x(), V1().y(), V1().z(), V2().x(), V2().y(),
		                           V2().z(), V3().x(), V3().y(), V3().z(), W1().x(),
		                           W1().y(), W1().z(), W2().x(), W2().y(), W2().z(),
		                           W3().x(), W3().y(), W3().z(), U1().x(), U1().y(),
		                           U1().z(), U2().x(), U2().y(), U2().z(), U3().x(),
		                           U3().y(), U3().z(), lx, ly, lz, d, mv_))
			d = 0;

		if (d < 0)
		{
			lx = -lx;
			ly = -ly;
			lz = -lz;
			d  = -d;
		}
		if (d != 0 && mv_ > mv)
			mv = mv_;
		return d != 0;
	}
}

template <typename IT, typename ET>
bool ImplicitPoint3T_TPI<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &lz,
                                                    IT &d) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.dfilter_cached)
		{
			cv.dfilter_cached = true;
			lambda3d_TPI_interval<IT>(
			  V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(), V3().x(),
			  V3().y(), V3().z(), W1().x(), W1().y(), W1().z(), W2().x(), W2().y(),
			  W2().z(), W3().x(), W3().y(), W3().z(), U1().x(), U1().y(), U1().z(),
			  U2().x(), U2().y(), U2().z(), U3().x(), U3().y(), U3().z(),
			  cv.dfilter_lambda_x, cv.dfilter_lambda_y, cv.dfilter_lambda_z,
			  cv.dfilter_denominator);
			if (cv.dfilter_denominator.is_negative())
			{
				cv.dfilter_lambda_x.invert();
				cv.dfilter_lambda_y.invert();
				cv.dfilter_lambda_z.invert();
				cv.dfilter_denominator.invert();
			}
		}

		lx = cv.dfilter_lambda_x;
		ly = cv.dfilter_lambda_y;
		lz = cv.dfilter_lambda_z;
		d  = cv.dfilter_denominator;
		return (cv.dfilter_denominator.is_sign_reliable());
	}
	else
	{
		lambda3d_TPI_interval<IT>(V1().x(), V1().y(), V1().z(), V2().x(), V2().y(),
		                          V2().z(), V3().x(), V3().y(), V3().z(), W1().x(),
		                          W1().y(), W1().z(), W2().x(), W2().y(), W2().z(),
		                          W3().x(), W3().y(), W3().z(), U1().x(), U1().y(),
		                          U1().z(), U2().x(), U2().y(), U2().z(), U3().x(),
		                          U3().y(), U3().z(), lx, ly, lz, d);
		if (d.is_negative())
		{
			lx.invert();
			ly.invert();
			lz.invert();
			d.invert();
		}
		return d.is_sign_reliable();
	}
}

template <typename IT, typename ET>
void ImplicitPoint3T_TPI<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &lz,
                                                 ET &d) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.exact_cached)
		{
			cv.exact_cached      = true;
			cv.exact_denominator = new ET();
			cv.exact_lambda_x    = new ET();
			cv.exact_lambda_y    = new ET();
			cv.exact_lambda_z    = new ET();
			lambda3d_TPI_exact<ET>(
			  ET(V1().x()), ET(V1().y()), ET(V1().z()), ET(V2().x()), ET(V2().y()),
			  ET(V2().z()), ET(V3().x()), ET(V3().y()), ET(V3().z()), ET(W1().x()),
			  ET(W1().y()), ET(W1().z()), ET(W2().x()), ET(W2().y()), ET(W2().z()),
			  ET(W3().x()), ET(W3().y()), ET(W3().z()), ET(U1().x()), ET(U1().y()),
			  ET(U1().z()), ET(U2().x()), ET(U2().y()), ET(U2().z()), ET(U3().x()),
			  ET(U3().y()), ET(U3().z()), *cv.exact_lambda_x, *cv.exact_lambda_y,
			  *cv.exact_lambda_z, *cv.exact_denominator);
			if (OMC::sign(*cv.exact_denominator) == Sign::NEGATIVE)
			{
				*cv.exact_lambda_x    = -*cv.exact_lambda_x;
				*cv.exact_lambda_y    = -*cv.exact_lambda_y;
				*cv.exact_lambda_z    = -*cv.exact_lambda_z;
				*cv.exact_denominator = -*cv.exact_denominator;
			}
		}
		lx = *cv.exact_lambda_x;
		ly = *cv.exact_lambda_y;
		lz = *cv.exact_lambda_z;
		d  = *cv.exact_denominator;
	}
	else
	{
		lambda3d_TPI_exact<ET>(
		  ET(V1().x()), ET(V1().y()), ET(V1().z()), ET(V2().x()), ET(V2().y()),
		  ET(V2().z()), ET(V3().x()), ET(V3().y()), ET(V3().z()), ET(W1().x()),
		  ET(W1().y()), ET(W1().z()), ET(W2().x()), ET(W2().y()), ET(W2().z()),
		  ET(W3().x()), ET(W3().y()), ET(W3().z()), ET(U1().x()), ET(U1().y()),
		  ET(U1().z()), ET(U2().x()), ET(U2().y()), ET(U2().z()), ET(U3().x()),
		  ET(U3().y()), ET(U3().z()), lx, ly, lz, d);
		if (OMC::sign(d) == Sign::NEGATIVE)
		{
			lx = -lx;
			ly = -ly;
			lz = -lz;
			d  = -d;
		}
	}
}

template <typename IT, typename ET>
void ImplicitPoint3T_TPI<IT, ET>::getExpansionLambda(FT **lx, int &lx_len,
                                                     FT **ly, int &ly_len,
                                                     FT **lz, int &lz_len,
                                                     FT **d, int &d_len) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);

		if (!cv.expansion_cached)
		{
			cv.expansion_cached = true;
			lambda3d_TPI_expansion(
			  V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(), V3().x(),
			  V3().y(), V3().z(), W1().x(), W1().y(), W1().z(), W2().x(), W2().y(),
			  W2().z(), W3().x(), W3().y(), W3().z(), U1().x(), U1().y(), U1().z(),
			  U2().x(), U2().y(), U2().z(), U3().x(), U3().y(), U3().z(),
			  &cv.expansion_lambda_x, cv.expansion_lambda_x_len,
			  &cv.expansion_lambda_y, cv.expansion_lambda_y_len,
			  &cv.expansion_lambda_z, cv.expansion_lambda_z_len,
			  &cv.expansion_denominator, cv.expansion_d_len);
			if (cv.expansion_denominator[cv.expansion_d_len - 1] < 0)
			{
				expansionObject o;
				o.Gen_Invert(cv.expansion_lambda_x_len, cv.expansion_lambda_x);
				o.Gen_Invert(cv.expansion_lambda_y_len, cv.expansion_lambda_y);
				o.Gen_Invert(cv.expansion_lambda_z_len, cv.expansion_lambda_z);
				o.Gen_Invert(cv.expansion_d_len, cv.expansion_denominator);
			}
			normalizeLambda3D(cv.expansion_lambda_x, cv.expansion_lambda_x_len,
			                  cv.expansion_lambda_y, cv.expansion_lambda_y_len,
			                  cv.expansion_lambda_z, cv.expansion_lambda_z_len,
			                  cv.expansion_denominator, cv.expansion_d_len);
		}
		*lx    = cv.expansion_lambda_x;
		*ly    = cv.expansion_lambda_y;
		*lz    = cv.expansion_lambda_z;
		*d     = cv.expansion_denominator;
		lx_len = cv.expansion_lambda_x_len;
		ly_len = cv.expansion_lambda_y_len;
		lz_len = cv.expansion_lambda_z_len;
		d_len  = cv.expansion_d_len;
	}
	else
	{
		lambda3d_TPI_expansion(
		  V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(), V3().x(),
		  V3().y(), V3().z(), W1().x(), W1().y(), W1().z(), W2().x(), W2().y(),
		  W2().z(), W3().x(), W3().y(), W3().z(), U1().x(), U1().y(), U1().z(),
		  U2().x(), U2().y(), U2().z(), U3().x(), U3().y(), U3().z(), lx, lx_len,
		  ly, ly_len, lz, lz_len, d, d_len);
		if ((*d)[d_len - 1] < 0)
		{
			expansionObject o;
			o.Gen_Invert(lx_len, *lx);
			o.Gen_Invert(ly_len, *ly);
			o.Gen_Invert(lz_len, *lz);
			o.Gen_Invert(d_len, *d);
		}
		normalizeLambda3D(*lx, lx_len, *ly, ly_len, *lz, lz_len, *d, d_len);
	}
}

} // namespace OMC