#pragma once

#include "ImplicitPoint3T.h"

#include "OpenMeshCraft/Geometry/Predicates/ImplicitPointPredicates.h"

namespace OMC {

template <typename IT, typename ET>
ImplicitPoint3T_LPI<IT, ET>::ImplicitPoint3T_LPI() noexcept
  : GP(PointType::LPI)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LPI<IT, ET>::ImplicitPoint3T_LPI(const EP &_p, const EP &_q,
                                                 const EP &_r, const EP &_s,
                                                 const EP &_t) noexcept
  : GP(PointType::LPI)
  , ip(&_p)
  , iq(&_q)
  , ir(&_r)
  , is(&_s)
  , it(&_t)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LPI<IT, ET>::~ImplicitPoint3T_LPI() noexcept
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LPI<IT, ET>::ImplicitPoint3T_LPI(const IP &rhs) noexcept
  : GP(static_cast<const GP &>(rhs))
  , ip(rhs.ip)
  , iq(rhs.iq)
  , ir(rhs.ir)
  , is(rhs.is)
  , it(rhs.it)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LPI<IT, ET>::ImplicitPoint3T_LPI(IP &&rhs) noexcept
  : GP(static_cast<GP &&>(rhs))
  , ip(rhs.ip)
  , iq(rhs.iq)
  , ir(rhs.ir)
  , is(rhs.is)
  , it(rhs.it)
{
}

template <typename IT, typename ET>
auto ImplicitPoint3T_LPI<IT, ET>::operator=(const IP &rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<const GP &>(rhs));
	ip                       = rhs.ip;
	iq                       = rhs.iq;
	ir                       = rhs.ir;
	is                       = rhs.is;
	it                       = rhs.it;
	return *this;
}

template <typename IT, typename ET>
auto ImplicitPoint3T_LPI<IT, ET>::operator=(IP &&rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<GP &&>(rhs));
	ip                       = rhs.ip;
	iq                       = rhs.iq;
	ir                       = rhs.ir;
	is                       = rhs.is;
	it                       = rhs.it;
	return *this;
}

#if defined(INDIRECT_PREDICATES)

template <typename IT, typename ET>
void ImplicitPoint3T_LPI<IT, ET>::get_Explicit(EP &e) const
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
auto ImplicitPoint3T_LPI<IT, ET>::to_Explicit() const -> EP
{
	EP e;
	get_Explicit(e);
	return e;
}

template <typename IT, typename ET>
bool ImplicitPoint3T_LPI<IT, ET>::getFilteredLambda(FT &lx, FT &ly, FT &lz,
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
			if (!lambda3d_LPI_filtered(
			      P().x(), P().y(), P().z(), Q().x(), Q().y(), Q().z(), R().x(),
			      R().y(), R().z(), S().x(), S().y(), S().z(), T().x(), T().y(),
			      T().z(), cv.ssfilter_denominator, cv.ssfilter_lambda_x,
			      cv.ssfilter_lambda_y, cv.ssfilter_lambda_z, cv.ssfilter_max_val))
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
		return (d != 0);
	}
	else
	{
		FT mv_ = 0;
		if (!lambda3d_LPI_filtered(P().x(), P().y(), P().z(), Q().x(), Q().y(),
		                           Q().z(), R().x(), R().y(), R().z(), S().x(),
		                           S().y(), S().z(), T().x(), T().y(), T().z(), d,
		                           lx, ly, lz, mv_))
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
bool ImplicitPoint3T_LPI<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &lz,
                                                    IT &d) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.dfilter_cached)
		{
			cv.dfilter_cached = true;
			lambda3d_LPI_interval<IT>(P().x(), P().y(), P().z(), Q().x(), Q().y(),
			                          Q().z(), R().x(), R().y(), R().z(), S().x(),
			                          S().y(), S().z(), T().x(), T().y(), T().z(),
			                          cv.dfilter_denominator, cv.dfilter_lambda_x,
			                          cv.dfilter_lambda_y, cv.dfilter_lambda_z);
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
		return (d.is_sign_reliable());
	}
	else
	{
		lambda3d_LPI_interval<IT>(P().x(), P().y(), P().z(), Q().x(), Q().y(),
		                          Q().z(), R().x(), R().y(), R().z(), S().x(),
		                          S().y(), S().z(), T().x(), T().y(), T().z(), d,
		                          lx, ly, lz);
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
void ImplicitPoint3T_LPI<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &lz,
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
			lambda3d_LPI_exact<ET>(ET(P().x()), ET(P().y()), ET(P().z()), ET(Q().x()),
			                       ET(Q().y()), ET(Q().z()), ET(R().x()), ET(R().y()),
			                       ET(R().z()), ET(S().x()), ET(S().y()), ET(S().z()),
			                       ET(T().x()), ET(T().y()), ET(T().z()),
			                       *cv.exact_denominator, *cv.exact_lambda_x,
			                       *cv.exact_lambda_y, *cv.exact_lambda_z);
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
		lambda3d_LPI_exact<ET>(ET(P().x()), ET(P().y()), ET(P().z()), ET(Q().x()),
		                       ET(Q().y()), ET(Q().z()), ET(R().x()), ET(R().y()),
		                       ET(R().z()), ET(S().x()), ET(S().y()), ET(S().z()),
		                       ET(T().x()), ET(T().y()), ET(T().z()), d, lx, ly,
		                       lz);
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
void ImplicitPoint3T_LPI<IT, ET>::getExpansionLambda(FT **lx, int &lx_len,
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
			lambda3d_LPI_expansion(P().x(), P().y(), P().z(), Q().x(), Q().y(),
			                       Q().z(), R().x(), R().y(), R().z(), S().x(),
			                       S().y(), S().z(), T().x(), T().y(), T().z(),
			                       &cv.expansion_denominator, cv.expansion_d_len,
			                       &cv.expansion_lambda_x, cv.expansion_lambda_x_len,
			                       &cv.expansion_lambda_y, cv.expansion_lambda_y_len,
			                       &cv.expansion_lambda_z, cv.expansion_lambda_z_len);
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
		lambda3d_LPI_expansion(P().x(), P().y(), P().z(), Q().x(), Q().y(), Q().z(),
		                       R().x(), R().y(), R().z(), S().x(), S().y(), S().z(),
		                       T().x(), T().y(), T().z(), d, d_len, lx, lx_len, ly,
		                       ly_len, lz, lz_len);
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

#elif defined(OFFSET_PREDICATES)

template <typename IT, typename ET>
void ImplicitPoint3T_LPI<IT, ET>::get_Explicit(EP &e) const
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
auto ImplicitPoint3T_LPI<IT, ET>::to_Explicit() const -> EP
{
	EP e;
	get_Explicit(e);
	return e;
}

template <typename IT, typename ET>
bool ImplicitPoint3T_LPI<IT, ET>::getFilteredLambda(FT &lx, FT &ly, FT &lz,
                                                    FT &d, FT &bx, FT &by,
                                                    FT &bz, FT &mv) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.ssfilter_cached)
		{
			cv.ssfilter_cached  = true;
			cv.ssfilter_max_val = 0;
			if (!lambda3d_LPI_filtered(
			      P().x(), P().y(), P().z(), Q().x(), Q().y(), Q().z(), R().x(),
			      R().y(), R().z(), S().x(), S().y(), S().z(), T().x(), T().y(),
			      T().z(), cv.ssfilter_denominator, cv.ssfilter_lambda_x,
			      cv.ssfilter_lambda_y, cv.ssfilter_lambda_z, cv.ssfilter_beta_x,
			      cv.ssfilter_beta_y, cv.ssfilter_beta_z, cv.ssfilter_max_val))
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
		return (d != 0);
	}
	else
	{
		FT mv_ = 0;
		if (!lambda3d_LPI_filtered(P().x(), P().y(), P().z(), Q().x(), Q().y(),
		                           Q().z(), R().x(), R().y(), R().z(), S().x(),
		                           S().y(), S().z(), T().x(), T().y(), T().z(), d,
		                           lx, ly, lz, bx, by, bz, mv_))
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
bool ImplicitPoint3T_LPI<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &lz,
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
			lambda3d_LPI_interval<IT>(
			  P().x(), P().y(), P().z(), Q().x(), Q().y(), Q().z(), R().x(), R().y(),
			  R().z(), S().x(), S().y(), S().z(), T().x(), T().y(), T().z(),
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
		return (d.is_sign_reliable());
	}
	else
	{
		lambda3d_LPI_interval<IT>(P().x(), P().y(), P().z(), Q().x(), Q().y(),
		                          Q().z(), R().x(), R().y(), R().z(), S().x(),
		                          S().y(), S().z(), T().x(), T().y(), T().z(), d,
		                          lx, ly, lz, bx, by, bz);
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
void ImplicitPoint3T_LPI<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &lz, ET &d,
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
			lambda3d_LPI_exact<ET>(
			  ET(P().x()), ET(P().y()), ET(P().z()), ET(Q().x()), ET(Q().y()),
			  ET(Q().z()), ET(R().x()), ET(R().y()), ET(R().z()), ET(S().x()),
			  ET(S().y()), ET(S().z()), ET(T().x()), ET(T().y()), ET(T().z()),
			  *cv.exact_denominator, *cv.exact_lambda_x, *cv.exact_lambda_y,
			  *cv.exact_lambda_z, *cv.exact_beta_x, *cv.exact_beta_y,
			  *cv.exact_beta_z);
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
		lambda3d_LPI_exact<ET>(ET(P().x()), ET(P().y()), ET(P().z()), ET(Q().x()),
		                       ET(Q().y()), ET(Q().z()), ET(R().x()), ET(R().y()),
		                       ET(R().z()), ET(S().x()), ET(S().y()), ET(S().z()),
		                       ET(T().x()), ET(T().y()), ET(T().z()), d, lx, ly, lz,
		                       bx, by, bz);
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
void ImplicitPoint3T_LPI<IT, ET>::getExpansionLambda(FT **lx, int &lx_len,
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
			lambda3d_LPI_expansion(
			  P().x(), P().y(), P().z(), Q().x(), Q().y(), Q().z(), R().x(), R().y(),
			  R().z(), S().x(), S().y(), S().z(), T().x(), T().y(), T().z(),
			  &cv.expansion_denominator, cv.expansion_d_len, &cv.expansion_lambda_x,
			  cv.expansion_lambda_x_len, &cv.expansion_lambda_y,
			  cv.expansion_lambda_y_len, &cv.expansion_lambda_z,
			  cv.expansion_lambda_z_len, cv.ssfilter_beta_x, cv.ssfilter_beta_y,
			  cv.ssfilter_beta_z);
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
		// beta xyz must be "double" floating point number
		bx     = cv.ssfilter_beta_x;
		by     = cv.ssfilter_beta_y;
		bz     = cv.ssfilter_beta_z;
	}
	else
	{
		lambda3d_LPI_expansion(P().x(), P().y(), P().z(), Q().x(), Q().y(), Q().z(),
		                       R().x(), R().y(), R().z(), S().x(), S().y(), S().z(),
		                       T().x(), T().y(), T().z(), d, d_len, lx, lx_len, ly,
		                       ly_len, lz, lz_len, bx, by, bz);
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

#endif

/// @brief This is just a function for profiling. It outputs maxvar under
/// indirect predicates to compare with other predicates
template <typename IT, typename ET>
auto ImplicitPoint3T_LPI<IT, ET>::getIndirectMaxVar() const -> FT
{
	FT a11   = P().x() - Q().x();
	FT a12   = P().y() - Q().y();
	FT a13   = P().z() - Q().z();
	FT a21   = S().x() - R().x();
	FT a22   = S().y() - R().y();
	FT a23   = S().z() - R().z();
	FT a31   = T().x() - R().x();
	FT a32   = T().y() - R().y();
	FT a33   = T().z() - R().z();
	FT px_rx = P().x() - R().x();
	FT py_ry = P().y() - R().y();
	FT pz_rz = P().z() - R().z();

	FT _tmp_fabs, max_var = 0.;
	if ((_tmp_fabs = fabs(P().x())) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(P().y())) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(P().z())) > max_var)
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

	return max_var;
}

/// @brief This is just a function for profiling. It outputs maxvar under offset
/// predicates to compare with other predicates
template <typename IT, typename ET>
auto ImplicitPoint3T_LPI<IT, ET>::getOffsetMaxVar() const -> FT
{
	FT xqp = Q().x() - P().x();
	FT yqp = Q().y() - P().y();
	FT zqp = Q().z() - P().z();
	FT xsr = S().x() - R().x();
	FT ysr = S().y() - R().y();
	FT zsr = S().z() - R().z();
	FT xtr = T().x() - R().x();
	FT ytr = T().y() - R().y();
	FT ztr = T().z() - R().z();
	FT xrq = R().x() - Q().x();
	FT yrq = R().y() - Q().y();
	FT zrq = R().z() - Q().z();

	FT _tmp_fabs, max_var = 0.;
	if ((_tmp_fabs = fabs(xqp)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yqp)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zqp)) > max_var)
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
	if ((_tmp_fabs = fabs(xrq)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yrq)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zrq)) > max_var)
		max_var = _tmp_fabs;

	return max_var;
}
} // namespace OMC