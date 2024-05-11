#pragma once

#include "ImplicitPoint3T.h"

#include "OpenMeshCraft/Geometry/Predicates/ImplicitPointPredicates.h"

namespace OMC {

#if defined(OMC_INDIRECT_PRED)

template <typename IT, typename ET>
ImplicitPoint3T_LNC<IT, ET>::ImplicitPoint3T_LNC() noexcept
  : GP(PointType::LNC)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LNC<IT, ET>::ImplicitPoint3T_LNC(const EP &_p, const EP &_q,
                                                 const double _t) noexcept
  : GP(PointType::LNC)
  , ip(&_p)
  , iq(&_q)
  , t(_t)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LNC<IT, ET>::~ImplicitPoint3T_LNC() noexcept
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LNC<IT, ET>::ImplicitPoint3T_LNC(const IP &rhs) noexcept
  : GP(static_cast<const GP &>(rhs))
  , ip(rhs.ip)
  , iq(rhs.iq)
  , t(rhs.t)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_LNC<IT, ET>::ImplicitPoint3T_LNC(IP &&rhs) noexcept
  : GP(static_cast<GP &&>(rhs))
  , ip(rhs.ip)
  , iq(rhs.iq)
  , t(rhs.t)
{
}

template <typename IT, typename ET>
auto ImplicitPoint3T_LNC<IT, ET>::operator=(const IP &rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<const GP &>(rhs));
	ip                       = rhs.ip;
	iq                       = rhs.iq;
	t                        = rhs.t;
	return *this;
}

template <typename IT, typename ET>
auto ImplicitPoint3T_LNC<IT, ET>::operator=(IP &&rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<GP &&>(rhs));
	ip                       = rhs.ip;
	iq                       = rhs.iq;
	t                        = rhs.t;
	return *this;
}

template <typename IT, typename ET>
void ImplicitPoint3T_LNC<IT, ET>::get_Explicit(EP &e) const
{
	FT lambda_x, lambda_y, lambda_z, lambda_d;
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
	e = EP(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
}

template <typename IT, typename ET>
auto ImplicitPoint3T_LNC<IT, ET>::to_Explicit() const -> EP
{
	EP e;
	get_Explicit(e);
	return e;
}

template <typename IT, typename ET>
bool ImplicitPoint3T_LNC<IT, ET>::getFilteredLambda(OMC_UNUSED FT &lx,
                                                    OMC_UNUSED FT &ly,
                                                    OMC_UNUSED FT &lz,
                                                    OMC_UNUSED FT &d,
                                                    OMC_UNUSED FT &mv) const
{
	OMC_EXIT("LNC point does not has filtered lambda.");
	// return false;
}

template <typename IT, typename ET>
bool ImplicitPoint3T_LNC<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &lz,
                                                    IT &d) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.dfilter_cached)
		{
			cv.dfilter_cached = true;
			lambda3d_LNC_interval<IT>(P().x(), P().y(), P().z(), Q().x(), Q().y(),
			                          Q().z(), T(), cv.dfilter_lambda_x,
			                          cv.dfilter_lambda_y, cv.dfilter_lambda_z,
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
		lambda3d_LNC_interval<IT>(P().x(), P().y(), P().z(), Q().x(), Q().y(),
		                          Q().z(), T(), lx, ly, lz, d);
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
void ImplicitPoint3T_LNC<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &lz,
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
			lambda3d_LNC_exact<ET>(ET(P().x()), ET(P().y()), ET(P().z()), ET(Q().x()),
			                       ET(Q().y()), ET(Q().z()), ET(T()),
			                       *cv.exact_lambda_x, *cv.exact_lambda_y,
			                       *cv.exact_lambda_z, *cv.exact_denominator);
		}
		lx = *cv.exact_lambda_x;
		ly = *cv.exact_lambda_y;
		lz = *cv.exact_lambda_z;
		d  = *cv.exact_denominator;
	}
	else
	{
		lambda3d_LNC_exact<ET>(ET(P().x()), ET(P().y()), ET(P().z()), ET(Q().x()),
		                       ET(Q().y()), ET(Q().z()), ET(T()), lx, ly, lz, d);
	}
}

template <typename IT, typename ET>
void ImplicitPoint3T_LNC<IT, ET>::getExpansionLambda(FT **lx, int &lx_len,
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

			lambda3d_LNC_expansion(P().x(), P().y(), P().z(), Q().x(), Q().y(),
			                       Q().z(), T(), &cv.expansion_lambda_x,
			                       cv.expansion_lambda_x_len, &cv.expansion_lambda_y,
			                       cv.expansion_lambda_y_len, &cv.expansion_lambda_z,
			                       cv.expansion_lambda_z_len,
			                       &cv.expansion_denominator, cv.expansion_d_len);
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
		lambda3d_LNC_expansion(P().x(), P().y(), P().z(), Q().x(), Q().y(), Q().z(),
		                       T(), lx, lx_len, ly, ly_len, lz, lz_len, d, d_len);
		normalizeLambda3D(*lx, lx_len, *ly, ly_len, *lz, lz_len, *d, d_len);
	}
}

#endif

} // namespace OMC