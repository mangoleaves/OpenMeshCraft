#pragma once

#include "ImplicitPoint2T.h"

#include "OpenMeshCraft/Geometry/Predicates/ImplicitPointPredicates.h"

namespace OMC {

template <typename IT, typename ET>
ImplicitPoint2T_SSI<IT, ET>::ImplicitPoint2T_SSI() noexcept
  : GP(PointType::SSI)
{
}

template <typename IT, typename ET>
ImplicitPoint2T_SSI<IT, ET>::ImplicitPoint2T_SSI(const EP &l11, const EP &l12,
                                                 const EP &l21,
                                                 const EP &l22) noexcept
  : GP(PointType::SSI)
  , l1_1(&l11)
  , l1_2(&l12)
  , l2_1(&l21)
  , l2_2(&l22)
{
}

template <typename IT, typename ET>
ImplicitPoint2T_SSI<IT, ET>::ImplicitPoint2T_SSI(const IP &rhs) noexcept
  : GP(static_cast<const GP &>(rhs))
  , l1_1(rhs.l1_1)
  , l1_2(rhs.l1_2)
  , l2_1(rhs.l2_1)
  , l2_2(rhs.l2_2)
{
}

template <typename IT, typename ET>
ImplicitPoint2T_SSI<IT, ET>::ImplicitPoint2T_SSI(IP &&rhs) noexcept
  : GP(static_cast<GP &&>(rhs))
  , l1_1(rhs.l1_1)
  , l1_2(rhs.l1_2)
  , l2_1(rhs.l2_1)
  , l2_2(rhs.l2_2)
{
}

template <typename IT, typename ET>
auto ImplicitPoint2T_SSI<IT, ET>::operator=(const IP &rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<const GP &>(rhs));
	l1_1                     = rhs.l1_1;
	l1_2                     = rhs.l1_2;
	l2_1                     = rhs.l2_1;
	l2_2                     = rhs.l2_2;
	return *this;
}

template <typename IT, typename ET>
auto ImplicitPoint2T_SSI<IT, ET>::operator=(IP &&rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<GP &&>(rhs));
	l1_1                     = rhs.l1_1;
	l1_2                     = rhs.l1_2;
	l2_1                     = rhs.l2_1;
	l2_2                     = rhs.l2_2;
	return *this;
}

template <typename IT, typename ET>
ImplicitPoint2T_SSI<IT, ET>::~ImplicitPoint2T_SSI() noexcept
{
}

template <typename IT, typename ET>
void ImplicitPoint2T_SSI<IT, ET>::get_Explicit(EP &e) const
{
	FT lambda_x, lambda_y, lambda_d, max_var = 0;
	if (!getFilteredLambda(lambda_x, lambda_y, lambda_d, max_var))
	{
		IT ilx, ily, id;
		if (!getIntervalLambda(ilx, ily, id))
		{
			ET elx, ely, eld;
			getExactLambda(elx, ely, eld);
			lambda_x = OMC::to_double(elx);
			lambda_y = OMC::to_double(ely);
			lambda_d = OMC::to_double(eld);
		}
		else
		{
			lambda_x = ilx.sup() + ilx.inf();
			lambda_y = ily.sup() + ily.inf();
			lambda_d = id.sup() + id.inf();
		}
	}
	e = EP(lambda_x / lambda_d, lambda_y / lambda_d);
}

template <typename IT, typename ET>
auto ImplicitPoint2T_SSI<IT, ET>::to_Explicit() const -> EP
{
	EP e;
	get_Explicit(e);
	return e;
}

template <typename IT, typename ET>
bool ImplicitPoint2T_SSI<IT, ET>::getFilteredLambda(FT &lx, FT &ly, FT &d,
                                                    FT &mv) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.ssfilter_cached)
		{
			cv.ssfilter_cached  = true;
			cv.ssfilter_max_val = 0;
			if (!lambda2d_SSI_filtered(L1_1().x(), L1_1().y(), L1_2().x(), L1_2().y(),
			                           L2_1().x(), L2_1().y(), L2_2().x(), L2_2().y(),
			                           cv.ssfilter_lambda_x, cv.ssfilter_lambda_y,
			                           cv.ssfilter_denominator, cv.ssfilter_max_val))
				cv.ssfilter_denominator = 0;

			if (cv.ssfilter_denominator < 0)
			{
				cv.ssfilter_lambda_x    = -cv.ssfilter_lambda_x;
				cv.ssfilter_lambda_y    = -cv.ssfilter_lambda_y;
				cv.ssfilter_denominator = -cv.ssfilter_denominator;
			}
		}

		lx = cv.ssfilter_lambda_x;
		ly = cv.ssfilter_lambda_y;
		d  = cv.ssfilter_denominator;
		if (cv.ssfilter_denominator != 0 && cv.ssfilter_max_val > mv)
			mv = cv.ssfilter_max_val;
		return (cv.ssfilter_denominator != 0);
	}
	else
	{
		FT mv_ = 0;
		if (!lambda2d_SSI_filtered(L1_1().x(), L1_1().y(), L1_2().x(), L1_2().y(),
		                           L2_1().x(), L2_1().y(), L2_2().x(), L2_2().y(),
		                           lx, ly, d, mv_))
			d = 0;

		if (d < 0)
		{
			lx = -lx;
			ly = -ly;
			d  = -d;
		}
		if (d != 0 && mv_ > mv)
			mv = mv_;
		return d != 0;
	}
}

template <typename IT, typename ET>
bool ImplicitPoint2T_SSI<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &d) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.dfilter_cached)
		{
			cv.dfilter_cached = true;
			lambda2d_SSI_interval<IT>(L1_1().x(), L1_1().y(), L1_2().x(), L1_2().y(),
			                          L2_1().x(), L2_1().y(), L2_2().x(), L2_2().y(),
			                          cv.dfilter_lambda_x, cv.dfilter_lambda_y,
			                          cv.dfilter_denominator);
			if (cv.dfilter_denominator.is_negative())
			{
				cv.dfilter_lambda_x.invert();
				cv.dfilter_lambda_y.invert();
				cv.dfilter_denominator.invert();
			}
		}

		lx = cv.dfilter_lambda_x;
		ly = cv.dfilter_lambda_y;
		d  = cv.dfilter_denominator;
		return (cv.dfilter_denominator.is_sign_reliable());
	}
	else
	{
		lambda2d_SSI_interval<IT>(L1_1().x(), L1_1().y(), L1_2().x(), L1_2().y(),
		                          L2_1().x(), L2_1().y(), L2_2().x(), L2_2().y(),
		                          lx, ly, d);
		if (d.is_negative())
		{
			lx.invert();
			ly.invert();
			d.invert();
		}
		return d.is_sign_reliable();
	}
}

template <typename IT, typename ET>
void ImplicitPoint2T_SSI<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &d) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.exact_cached)
		{
			cv.exact_cached      = true;
			cv.exact_lambda_x    = new ET();
			cv.exact_lambda_y    = new ET();
			cv.exact_denominator = new ET();
			lambda2d_SSI_exact<ET>(ET(L1_1().x()), ET(L1_1().y()), ET(L1_2().x()),
			                       ET(L1_2().y()), ET(L2_1().x()), ET(L2_1().y()),
			                       ET(L2_2().x()), ET(L2_2().y()), *cv.exact_lambda_x,
			                       *cv.exact_lambda_y, *cv.exact_denominator);
			if (OMC::sign(*cv.exact_denominator) == Sign::NEGATIVE)
			{
				*cv.exact_lambda_x    = -*cv.exact_lambda_x;
				*cv.exact_lambda_y    = -*cv.exact_lambda_y;
				*cv.exact_denominator = -*cv.exact_denominator;
			}
		}
		lx = *cv.exact_lambda_x;
		ly = *cv.exact_lambda_y;
		d  = *cv.exact_denominator;
	}
	else
	{
		lambda2d_SSI_exact<ET>(ET(L1_1().x()), ET(L1_1().y()), ET(L1_2().x()),
		                       ET(L1_2().y()), ET(L2_1().x()), ET(L2_1().y()),
		                       ET(L2_2().x()), ET(L2_2().y()), lx, ly, d);
		if (OMC::sign(d) == Sign::NEGATIVE)
		{
			lx = -lx;
			ly = -ly;
			d  = -d;
		}
	}
}

template <typename IT, typename ET>
void ImplicitPoint2T_SSI<IT, ET>::getExpansionLambda(FT **lx, int &lx_len,
                                                     FT **ly, int &ly_len,
                                                     FT **d, int &d_len) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);

		if (!cv.expansion_cached)
		{
			cv.expansion_cached = true;
			lambda2d_SSI_expansion(L1_1().x(), L1_1().y(), L1_2().x(), L1_2().y(),
			                       L2_1().x(), L2_1().y(), L2_2().x(), L2_2().y(),
			                       &cv.expansion_lambda_x, cv.expansion_lambda_x_len,
			                       &cv.expansion_lambda_y, cv.expansion_lambda_y_len,
			                       &cv.expansion_denominator, cv.expansion_d_len);
			if (cv.expansion_denominator[cv.expansion_d_len - 1] < 0)
			{
				expansionObject o;
				o.Gen_Invert(cv.expansion_lambda_x_len, cv.expansion_lambda_x);
				o.Gen_Invert(cv.expansion_lambda_y_len, cv.expansion_lambda_y);
				o.Gen_Invert(cv.expansion_d_len, cv.expansion_denominator);
			}
		}
		*lx    = cv.expansion_lambda_x;
		*ly    = cv.expansion_lambda_y;
		*d     = cv.expansion_denominator;
		lx_len = cv.expansion_lambda_x_len;
		ly_len = cv.expansion_lambda_y_len;
		d_len  = cv.expansion_d_len;
	}
	else
	{
		lambda2d_SSI_expansion(L1_1().x(), L1_1().y(), L1_2().x(), L1_2().y(),
		                       L2_1().x(), L2_1().y(), L2_2().x(), L2_2().y(), lx,
		                       lx_len, ly, ly_len, d, d_len);
		if ((*d)[d_len - 1] < 0)
		{
			expansionObject o;
			o.Gen_Invert(lx_len, *lx);
			o.Gen_Invert(ly_len, *ly);
			o.Gen_Invert(d_len, *d);
		}
	}
}

} // namespace OMC