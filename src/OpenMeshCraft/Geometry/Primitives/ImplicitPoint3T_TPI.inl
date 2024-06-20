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
#ifdef OMC_CACHE_SSF
	m_maxvar = 0;
	#if defined(OMC_OFFSET_PRED)
	FT bx, by, bz;
	if (!lambda3d_TPI_filtered(
	      V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(), V3().x(),
	      V3().y(), V3().z(), W1().x(), W1().y(), W1().z(), W2().x(), W2().y(),
	      W2().z(), W3().x(), W3().y(), W3().z(), U1().x(), U1().y(), U1().z(),
	      U2().x(), U2().y(), U2().z(), U3().x(), U3().y(), U3().z(), m_ld, m_lx,
	      m_ly, m_lz, bx, by, bz, m_maxvar))
		m_ld = 0;
	#elif defined(OMC_INDIRECT_PRED)
	if (!lambda3d_TPI_filtered(
	      V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(), V3().x(),
	      V3().y(), V3().z(), W1().x(), W1().y(), W1().z(), W2().x(), W2().y(),
	      W2().z(), W3().x(), W3().y(), W3().z(), U1().x(), U1().y(), U1().z(),
	      U2().x(), U2().y(), U2().z(), U3().x(), U3().y(), U3().z(), m_lx, m_ly,
	      m_lz, m_ld, m_maxvar))
		m_ld = 0;
	#endif

	if (m_ld < 0)
	{
		m_lx = -m_lx;
		m_ly = -m_ly;
		m_lz = -m_lz;
		m_ld = -m_ld;
	}
#endif
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
#ifdef OMC_CACHE_SSF
	m_lx     = rhs.m_lx;
	m_ly     = rhs.m_ly;
	m_lz     = rhs.m_lz;
	m_ld     = rhs.m_ld;
	m_maxvar = rhs.m_maxvar;
#endif
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
#ifdef OMC_CACHE_SSF
	m_lx     = rhs.m_lx;
	m_ly     = rhs.m_ly;
	m_lz     = rhs.m_lz;
	m_ld     = rhs.m_ld;
	m_maxvar = rhs.m_maxvar;
#endif
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
#ifdef OMC_CACHE_SSF
	m_lx     = rhs.m_lx;
	m_ly     = rhs.m_ly;
	m_lz     = rhs.m_lz;
	m_ld     = rhs.m_ld;
	m_maxvar = rhs.m_maxvar;
#endif
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
#ifdef OMC_CACHE_SSF
	m_lx     = rhs.m_lx;
	m_ly     = rhs.m_ly;
	m_lz     = rhs.m_lz;
	m_ld     = rhs.m_ld;
	m_maxvar = rhs.m_maxvar;
#endif
	return *this;
}

#if defined(OMC_INDIRECT_PRED)

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
	#ifndef OMC_CACHE_SSF
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
	#else
	lx = m_lx;
	ly = m_ly;
	lz = m_lz;
	d  = m_ld;
	if (d != 0 && m_maxvar > mv)
		mv = m_maxvar;
	return d != 0;
	#endif
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

#elif defined(OMC_OFFSET_PRED)

template <typename IT, typename ET>
void ImplicitPoint3T_TPI<IT, ET>::get_Explicit(EP &e) const
{
	FT lambda_x, lambda_y, lambda_z, lambda_d, beta_x, beta_y, beta_z,
	  max_var = 0;
	if (!getFilteredLambda(lambda_x, lambda_y, lambda_z, lambda_d, beta_x, beta_y,
	                       beta_z, max_var))
	{
		IT ilx, ily, ilz, id, ibx, iby, ibz;
		if (!getIntervalLambda(ilx, ily, ilz, id, ibx, iby, ibz))
		{
			ET elx, ely, elz, eld, ebx, eby, ebz;
			getExactLambda(elx, ely, elz, eld, ebx, eby, ebz);
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
	e = EP(beta_x + lambda_x / lambda_d, beta_y + lambda_y / lambda_d,
	       beta_z + lambda_z / lambda_d);
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
                                                    FT &d, FT &bx, FT &by,
                                                    FT &bz, FT &mv) const
{
	#ifndef OMC_CACHE_SSF
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
			      U3().x(), U3().y(), U3().z(), cv.ssfilter_denominator,
			      cv.ssfilter_lambda_x, cv.ssfilter_lambda_y, cv.ssfilter_lambda_z,
			      cv.ssfilter_beta_x, cv.ssfilter_beta_y, cv.ssfilter_beta_z,
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
		bx = cv.ssfilter_beta_x;
		by = cv.ssfilter_beta_y;
		bz = cv.ssfilter_beta_z;
		if (cv.ssfilter_denominator != 0 && cv.ssfilter_max_val > mv)
			mv = cv.ssfilter_max_val;
		return (cv.ssfilter_denominator != 0);
	}
	else
	{
		FT mv_ = 0;
		if (!lambda3d_TPI_filtered(
		      V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(), V3().x(),
		      V3().y(), V3().z(), W1().x(), W1().y(), W1().z(), W2().x(), W2().y(),
		      W2().z(), W3().x(), W3().y(), W3().z(), U1().x(), U1().y(), U1().z(),
		      U2().x(), U2().y(), U2().z(), U3().x(), U3().y(), U3().z(), d, lx, ly,
		      lz, bx, by, bz, mv_))
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
	#else
	lx = m_lx;
	ly = m_ly;
	lz = m_lz;
	d  = m_ld;
	bx = V1().x();
	by = V1().y();
	bz = V1().z();
	if (d != 0 && m_maxvar > mv)
		mv = m_maxvar;
	return d != 0;
	#endif
}

template <typename IT, typename ET>
bool ImplicitPoint3T_TPI<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &lz,
                                                    IT &d, IT &bx, IT &by,
                                                    IT &bz) const
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
			  cv.dfilter_denominator, cv.dfilter_lambda_x, cv.dfilter_lambda_y,
			  cv.dfilter_lambda_z, cv.dfilter_beta_x, cv.dfilter_beta_y,
			  cv.dfilter_beta_z);
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
		// beta xyz must be "double" floating point number
		bx = cv.dfilter_beta_x;
		by = cv.dfilter_beta_y;
		bz = cv.dfilter_beta_z;
		return (cv.dfilter_denominator.is_sign_reliable());
	}
	else
	{
		lambda3d_TPI_interval<IT>(V1().x(), V1().y(), V1().z(), V2().x(), V2().y(),
		                          V2().z(), V3().x(), V3().y(), V3().z(), W1().x(),
		                          W1().y(), W1().z(), W2().x(), W2().y(), W2().z(),
		                          W3().x(), W3().y(), W3().z(), U1().x(), U1().y(),
		                          U1().z(), U2().x(), U2().y(), U2().z(), U3().x(),
		                          U3().y(), U3().z(), d, lx, ly, lz, bx, by, bz);
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
void ImplicitPoint3T_TPI<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &lz, ET &d,
                                                 ET &bx, ET &by, ET &bz) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.exact_cached)
		{
			cv.alloc_ET();
			cv.exact_cached = true;
			lambda3d_TPI_exact<ET>(
			  ET(V1().x()), ET(V1().y()), ET(V1().z()), ET(V2().x()), ET(V2().y()),
			  ET(V2().z()), ET(V3().x()), ET(V3().y()), ET(V3().z()), ET(W1().x()),
			  ET(W1().y()), ET(W1().z()), ET(W2().x()), ET(W2().y()), ET(W2().z()),
			  ET(W3().x()), ET(W3().y()), ET(W3().z()), ET(U1().x()), ET(U1().y()),
			  ET(U1().z()), ET(U2().x()), ET(U2().y()), ET(U2().z()), ET(U3().x()),
			  ET(U3().y()), ET(U3().z()), *cv.exact_denominator, *cv.exact_lambda_x,
			  *cv.exact_lambda_y, *cv.exact_lambda_z, *cv.exact_beta_x,
			  *cv.exact_beta_y, *cv.exact_beta_z);
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
		// beta xyz must be "double" floating point number
		bx = *cv.exact_beta_x;
		by = *cv.exact_beta_y;
		bz = *cv.exact_beta_z;
	}
	else
	{
		lambda3d_TPI_exact<ET>(
		  ET(V1().x()), ET(V1().y()), ET(V1().z()), ET(V2().x()), ET(V2().y()),
		  ET(V2().z()), ET(V3().x()), ET(V3().y()), ET(V3().z()), ET(W1().x()),
		  ET(W1().y()), ET(W1().z()), ET(W2().x()), ET(W2().y()), ET(W2().z()),
		  ET(W3().x()), ET(W3().y()), ET(W3().z()), ET(U1().x()), ET(U1().y()),
		  ET(U1().z()), ET(U2().x()), ET(U2().y()), ET(U2().z()), ET(U3().x()),
		  ET(U3().y()), ET(U3().z()), d, lx, ly, lz, bx, by, bz);
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
                                                     FT **d, int &d_len, FT &bx,
                                                     FT &by, FT &bz) const
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
			  &cv.expansion_denominator, cv.expansion_d_len, &cv.expansion_lambda_x,
			  cv.expansion_lambda_x_len, &cv.expansion_lambda_y,
			  cv.expansion_lambda_y_len, &cv.expansion_lambda_z,
			  cv.expansion_lambda_z_len, cv.ssfilter_beta_x, cv.ssfilter_beta_y,
			  cv.ssfilter_beta_z);
			expansionObject o;
	#ifdef OMC_COMPRESS_EXPANSION
			o.CompressIf(cv.expansion_lambda_x_len, cv.expansion_lambda_x);
			o.CompressIf(cv.expansion_lambda_y_len, cv.expansion_lambda_y);
			o.CompressIf(cv.expansion_lambda_z_len, cv.expansion_lambda_z);
			o.CompressIf(cv.expansion_d_len, cv.expansion_denominator);
	#endif
			if (cv.expansion_denominator[cv.expansion_d_len - 1] < 0)
			{
				o.Gen_Invert(cv.expansion_lambda_x_len, cv.expansion_lambda_x);
				o.Gen_Invert(cv.expansion_lambda_y_len, cv.expansion_lambda_y);
				o.Gen_Invert(cv.expansion_lambda_z_len, cv.expansion_lambda_z);
				o.Gen_Invert(cv.expansion_d_len, cv.expansion_denominator);
			}
			normalizeLambda3D(cv.expansion_lambda_x, cv.expansion_lambda_x_len,
			                  cv.expansion_lambda_y, cv.expansion_lambda_y_len,
			                  cv.expansion_lambda_z, cv.expansion_lambda_z_len,
			                  cv.expansion_denominator, cv.expansion_d_len);
	#ifdef OMC_UPDATE_INTERVAL_BY_EXPANSION
			std::pair<double, double> it_; // -inf, sup
			it_ = o.To_Interval(cv.expansion_lambda_x_len, cv.expansion_lambda_x);
			cv.dfilter_lambda_x = IT(it_.first, it_.second);
			it_ = o.To_Interval(cv.expansion_lambda_y_len, cv.expansion_lambda_y);
			cv.dfilter_lambda_y = IT(it_.first, it_.second);
			it_ = o.To_Interval(cv.expansion_lambda_z_len, cv.expansion_lambda_z);
			cv.dfilter_lambda_z = IT(it_.first, it_.second);
			it_ = o.To_Interval(cv.expansion_d_len, cv.expansion_denominator);
			cv.dfilter_denominator = IT(it_.first, it_.second);
	#endif
		}
		*lx    = cv.expansion_lambda_x;
		*ly    = cv.expansion_lambda_y;
		*lz    = cv.expansion_lambda_z;
		*d     = cv.expansion_denominator;
		lx_len = cv.expansion_lambda_x_len;
		ly_len = cv.expansion_lambda_y_len;
		lz_len = cv.expansion_lambda_z_len;
		d_len  = cv.expansion_d_len;
		// beta xyz must be "double" floating point number
		bx     = cv.ssfilter_beta_x;
		by     = cv.ssfilter_beta_y;
		bz     = cv.ssfilter_beta_z;
	}
	else
	{
		lambda3d_TPI_expansion(
		  V1().x(), V1().y(), V1().z(), V2().x(), V2().y(), V2().z(), V3().x(),
		  V3().y(), V3().z(), W1().x(), W1().y(), W1().z(), W2().x(), W2().y(),
		  W2().z(), W3().x(), W3().y(), W3().z(), U1().x(), U1().y(), U1().z(),
		  U2().x(), U2().y(), U2().z(), U3().x(), U3().y(), U3().z(), d, d_len, lx,
		  lx_len, ly, ly_len, lz, lz_len, bx, by, bz);
		expansionObject o;
	#ifdef OMC_COMPRESS_EXPANSION
		o.CompressIf(lx_len, *lx);
		o.CompressIf(ly_len, *ly);
		o.CompressIf(lz_len, *lz);
		o.CompressIf(d_len, *d);
	#endif
		if ((*d)[d_len - 1] < 0)
		{
			o.Gen_Invert(lx_len, *lx);
			o.Gen_Invert(ly_len, *ly);
			o.Gen_Invert(lz_len, *lz);
			o.Gen_Invert(d_len, *d);
		}
		normalizeLambda3D(*lx, lx_len, *ly, ly_len, *lz, lz_len, *d, d_len);
	}
}

#endif

/// @brief This is just a function for profiling. It outputs maxvar under
/// indirect predicates to compare with other predicates
template <typename IT, typename ET>
auto ImplicitPoint3T_TPI<IT, ET>::getIndirectMaxVar() const -> FT
{
	FT ov1x = V1().x();
	FT ov1y = V1().y();
	FT ov1z = V1().z();
	FT ov2x = V2().x();
	FT ov2y = V2().y();
	FT ov2z = V2().z();
	FT ov3x = V3().x();
	FT ov3y = V3().y();
	FT ov3z = V3().z();
	FT ow1x = W1().x();
	FT ow1y = W1().y();
	FT ow1z = W1().z();
	FT ow2x = W2().x();
	FT ow2y = W2().y();
	FT ow2z = W2().z();
	FT ow3x = W3().x();
	FT ow3y = W3().y();
	FT ow3z = W3().z();
	FT ou1x = U1().x();
	FT ou1y = U1().y();
	FT ou1z = U1().z();
	FT ou2x = U2().x();
	FT ou2y = U2().y();
	FT ou2z = U2().z();
	FT ou3x = U3().x();
	FT ou3y = U3().y();
	FT ou3z = U3().z();

	FT v3x = ov3x - ov2x;
	FT v3y = ov3y - ov2y;
	FT v3z = ov3z - ov2z;
	FT v2x = ov2x - ov1x;
	FT v2y = ov2y - ov1y;
	FT v2z = ov2z - ov1z;
	FT w3x = ow3x - ow2x;
	FT w3y = ow3y - ow2y;
	FT w3z = ow3z - ow2z;
	FT w2x = ow2x - ow1x;
	FT w2y = ow2y - ow1y;
	FT w2z = ow2z - ow1z;
	FT u3x = ou3x - ou2x;
	FT u3y = ou3y - ou2y;
	FT u3z = ou3z - ou2z;
	FT u2x = ou2x - ou1x;
	FT u2y = ou2y - ou1y;
	FT u2z = ou2z - ou1z;

	FT _tmp_fabs, max_var = 0.;
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

	return max_var;
}

/// @brief This is just a function for profiling. It outputs maxvar under offset
/// predicates to compare with other predicates
template <typename IT, typename ET>
auto ImplicitPoint3T_TPI<IT, ET>::getOffsetMaxVar() const -> FT
{
	FT xa = V1().x();
	FT ya = V1().y();
	FT za = V1().z();
	FT xb = V2().x();
	FT yb = V2().y();
	FT zb = V2().z();
	FT xc = V3().x();
	FT yc = V3().y();
	FT zc = V3().z();
	FT xo = W1().x();
	FT yo = W1().y();
	FT zo = W1().z();
	FT xp = W2().x();
	FT yp = W2().y();
	FT zp = W2().z();
	FT xq = W3().x();
	FT yq = W3().y();
	FT zq = W3().z();
	FT xr = U1().x();
	FT yr = U1().y();
	FT zr = U1().z();
	FT xs = U2().x();
	FT ys = U2().y();
	FT zs = U2().z();
	FT xt = U3().x();
	FT yt = U3().y();
	FT zt = U3().z();

	FT xpo = xp - xo;
	FT ypo = yp - yo;
	FT zpo = zp - zo;
	FT xqo = xq - xo;
	FT yqo = yq - yo;
	FT zqo = zq - zo;
	FT xsr = xs - xr;
	FT ysr = ys - yr;
	FT zsr = zs - zr;
	FT xtr = xt - xr;
	FT ytr = yt - yr;
	FT ztr = zt - zr;
	FT xoa = xo - xa;
	FT yoa = yo - ya;
	FT zoa = zo - za;
	FT xra = xr - xa;
	FT yra = yr - ya;
	FT zra = zr - za;
	FT xba = xb - xa;
	FT yba = yb - ya;
	FT zba = zb - za;
	FT xca = xc - xa;
	FT yca = yc - ya;
	FT zca = zc - za;

	FT _tmp_fabs, max_var = 0.;
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

	return max_var;
}
} // namespace OMC