#pragma once

#include "ImplicitPoint3T.h"

#include "OpenMeshCraft/Geometry/Predicates/ImplicitPointPredicates.h"

namespace OMC {

template <typename IT, typename ET>
ImplicitPoint3T_SSI<IT, ET>::ImplicitPoint3T_SSI() noexcept
  : GP(PointType::SSI)
{
}

template <typename IT, typename ET>
ImplicitPoint3T_SSI<IT, ET>::ImplicitPoint3T_SSI(const EP &_a, const EP &_b,
                                                 const EP &_p, const EP &_q,
                                                 int _plane) noexcept
  : GP(PointType::SSI)
  , ia(&_a)
  , ib(&_b)
  , ip(&_p)
  , iq(&_q)
  , plane(_plane)
{
#ifdef OMC_CACHE_SSF
	m_maxvar = 0;
	#if defined(OMC_OFFSET_PRED)
	if (plane == 0) // YZ
	{
		// first yz, then x
		if (!lambda3d_SSI_filtered(A().y(), A().z(), A().x(), B().y(), B().z(),
		                           B().x(), P().y(), P().z(), Q().y(), Q().z(),
		                           m_ld, m_ly, m_lz, m_lx, m_maxvar))
			m_ld = 0;
	}
	else if (plane == 1) // ZX
	{
		// first zx, then y
		if (!lambda3d_SSI_filtered(A().z(), A().x(), A().y(), B().z(), B().x(),
		                           B().y(), P().z(), P().x(), Q().z(), Q().x(),
		                           m_ld, m_lz, m_lx, m_ly, m_maxvar))
			m_ld = 0;
	}
	else if (plane == 2) // XY
	{
		// first xy, then z
		if (!lambda3d_SSI_filtered(A().x(), A().y(), A().z(), B().x(), B().y(),
		                           B().z(), P().x(), P().y(), Q().x(), Q().y(),
		                           m_ld, m_lx, m_ly, m_lz, m_maxvar))
			m_ld = 0;
	}
	else
	{
		OMC_EXIT("wrong plane in SSI.");
	}
	#elif defined(OMC_INDIRECT_PRED)
	if (plane == 0) // YZ
	{
		// first yz, then x
		if (!lambda3d_SSI_filtered(A().y(), A().z(), A().x(), B().y(), B().z(),
		                           B().x(), P().y(), P().z(), Q().y(), Q().z(),
		                           m_ly, m_lz, m_lx, m_ld, m_maxvar))
			m_ld = 0;
	}
	else if (plane == 1) // ZX
	{
		// first zx, then y
		if (!lambda3d_SSI_filtered(A().z(), A().x(), A().y(), B().z(), B().x(),
		                           B().y(), P().z(), P().x(), Q().z(), Q().x(),
		                           m_lz, m_lx, m_ly, m_ld, m_maxvar))
			m_ld = 0;
	}
	else if (plane == 2) // XY
	{
		// first xy, then z
		if (!lambda3d_SSI_filtered(A().x(), A().y(), A().z(), B().x(), B().y(),
		                           B().z(), P().x(), P().y(), Q().x(), Q().y(),
		                           m_lx, m_ly, m_lz, m_ld, m_maxvar))
			m_ld = 0;
	}
	else
	{
		OMC_EXIT("wrong plane in SSI.");
	}
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
ImplicitPoint3T_SSI<IT, ET>::~ImplicitPoint3T_SSI() noexcept
{
}

template <typename IT, typename ET>
ImplicitPoint3T_SSI<IT, ET>::ImplicitPoint3T_SSI(const IP &rhs) noexcept
  : GP(static_cast<const GP &>(rhs))
  , ia(rhs.ia)
  , ib(rhs.ib)
  , ip(rhs.ip)
  , iq(rhs.iq)
  , plane(rhs.plane)
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
ImplicitPoint3T_SSI<IT, ET>::ImplicitPoint3T_SSI(IP &&rhs) noexcept
  : GP(static_cast<GP &&>(rhs))
  , ia(rhs.ia)
  , ib(rhs.ib)
  , ip(rhs.ip)
  , iq(rhs.iq)
  , plane(rhs.plane)
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
auto ImplicitPoint3T_SSI<IT, ET>::operator=(const IP &rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<const GP &>(rhs));
	ia                       = rhs.ia;
	ib                       = rhs.ib;
	ip                       = rhs.ip;
	iq                       = rhs.iq;
	plane                    = rhs.plane;
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
auto ImplicitPoint3T_SSI<IT, ET>::operator=(IP &&rhs) -> IP &
{
	*static_cast<GP *>(this) = (static_cast<GP &&>(rhs));
	ia                       = rhs.ia;
	ib                       = rhs.ib;
	ip                       = rhs.ip;
	iq                       = rhs.iq;
	plane                    = rhs.plane;
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
void ImplicitPoint3T_SSI<IT, ET>::get_Explicit(EP &e) const
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
auto ImplicitPoint3T_SSI<IT, ET>::to_Explicit() const -> EP
{
	EP e;
	get_Explicit(e);
	return e;
}

template <typename IT, typename ET>
bool ImplicitPoint3T_SSI<IT, ET>::getFilteredLambda(FT &lx, FT &ly, FT &lz,
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
			if (plane == 0) // YZ
			{
				// first yz, then x
				if (!lambda3d_SSI_filtered(
				      A().y(), A().z(), A().x(), B().y(), B().z(), B().x(), P().y(),
				      P().z(), Q().y(), Q().z(), cv.ssfilter_lambda_y,
				      cv.ssfilter_lambda_z, cv.ssfilter_lambda_x,
				      cv.ssfilter_denominator, cv.ssfilter_max_val))
					cv.ssfilter_denominator = 0;
			}
			else if (plane == 1) // ZX
			{
				// first zx, then y
				if (!lambda3d_SSI_filtered(
				      A().z(), A().x(), A().y(), B().z(), B().x(), B().y(), P().z(),
				      P().x(), Q().z(), Q().x(), cv.ssfilter_lambda_z,
				      cv.ssfilter_lambda_x, cv.ssfilter_lambda_y,
				      cv.ssfilter_denominator, cv.ssfilter_max_val))
					cv.ssfilter_denominator = 0;
			}
			else if (plane == 2) // XY
			{
				// first xy, then z
				if (!lambda3d_SSI_filtered(
				      A().x(), A().y(), A().z(), B().x(), B().y(), B().z(), P().x(),
				      P().y(), Q().x(), Q().y(), cv.ssfilter_lambda_x,
				      cv.ssfilter_lambda_y, cv.ssfilter_lambda_z,
				      cv.ssfilter_denominator, cv.ssfilter_max_val))
					cv.ssfilter_denominator = 0;
			}
			else
			{
				OMC_ASSERT(false, "plane is not initialized.");
			}

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

		if (plane == 0) // YZ
		{
			// first yz, then x
			if (!lambda3d_SSI_filtered(A().y(), A().z(), A().x(), B().y(), B().z(),
			                           B().x(), P().y(), P().z(), Q().y(), Q().z(),
			                           ly, lz, lx, d, mv_))
				d = 0;
		}
		else if (plane == 1) // ZX
		{
			// first zx, then y
			if (!lambda3d_SSI_filtered(A().z(), A().x(), A().y(), B().z(), B().x(),
			                           B().y(), P().z(), P().x(), Q().z(), Q().x(),
			                           lz, lx, ly, d, mv_))
				d = 0;
		}
		else if (plane == 2) // XY
		{
			// first xy, then z
			if (!lambda3d_SSI_filtered(A().x(), A().y(), A().z(), B().x(), B().y(),
			                           B().z(), P().x(), P().y(), Q().x(), Q().y(),
			                           lx, ly, lz, d, mv_))
				d = 0;
		}
		else
		{
			OMC_ASSERT(false, "plane is not initialized.");
		}

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
bool ImplicitPoint3T_SSI<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &lz,
                                                    IT &d) const
{
	if (global_cached_values.is_enabled())
	{
		typename GCV::OnePointCachedValues &cv =
		  global_cached_values.get((void *)this);
		if (!cv.dfilter_cached)
		{
			cv.dfilter_cached = true;
			if (plane == 0) // YZ
			{
				// first yz, then x
				lambda3d_SSI_interval<IT>(A().y(), A().z(), A().x(), B().y(), B().z(),
				                          B().x(), P().y(), P().z(), Q().y(), Q().z(),
				                          cv.dfilter_lambda_y, cv.dfilter_lambda_z,
				                          cv.dfilter_lambda_x, cv.dfilter_denominator);
			}
			else if (plane == 1) // ZX
			{
				// first zx, then y
				lambda3d_SSI_interval<IT>(A().z(), A().x(), A().y(), B().z(), B().x(),
				                          B().y(), P().z(), P().x(), Q().z(), Q().x(),
				                          cv.dfilter_lambda_z, cv.dfilter_lambda_x,
				                          cv.dfilter_lambda_y, cv.dfilter_denominator);
			}
			else if (plane == 2) // XY
			{
				// first xy, then z
				lambda3d_SSI_interval<IT>(A().x(), A().y(), A().z(), B().x(), B().y(),
				                          B().z(), P().x(), P().y(), Q().x(), Q().y(),
				                          cv.dfilter_lambda_x, cv.dfilter_lambda_y,
				                          cv.dfilter_lambda_z, cv.dfilter_denominator);
			}
			else
			{
				OMC_ASSERT(false, "plane is not initialized.");
			}
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
		if (plane == 0) // YZ
		{
			// first yz, then x
			lambda3d_SSI_interval<IT>(A().y(), A().z(), A().x(), B().y(), B().z(),
			                          B().x(), P().y(), P().z(), Q().y(), Q().z(), ly,
			                          lz, lx, d);
		}
		else if (plane == 1) // ZX
		{
			// first zx, then y
			lambda3d_SSI_interval<IT>(A().z(), A().x(), A().y(), B().z(), B().x(),
			                          B().y(), P().z(), P().x(), Q().z(), Q().x(), lz,
			                          lx, ly, d);
		}
		else if (plane == 2) // XY
		{
			// first xy, then z
			lambda3d_SSI_interval<IT>(A().x(), A().y(), A().z(), B().x(), B().y(),
			                          B().z(), P().x(), P().y(), Q().x(), Q().y(), lx,
			                          ly, lz, d);
		}
		else
		{
			OMC_ASSERT(false, "plane is not initialized.");
		}
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
void ImplicitPoint3T_SSI<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &lz,
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
			if (plane == 0) // yz
			{
				// first yz, then x
				lambda3d_SSI_exact<ET>(
				  ET(A().y()), ET(A().z()), ET(A().x()), ET(B().y()), ET(B().z()),
				  ET(B().x()), ET(P().y()), ET(P().z()), ET(Q().y()), ET(Q().z()),
				  *cv.exact_lambda_y, *cv.exact_lambda_z, *cv.exact_lambda_x,
				  *cv.exact_denominator);
			}
			else if (plane == 1) // zx
			{
				// first zx, then y
				lambda3d_SSI_exact<ET>(
				  ET(A().z()), ET(A().x()), ET(A().y()), ET(B().z()), ET(B().x()),
				  ET(B().y()), ET(P().z()), ET(P().x()), ET(Q().z()), ET(Q().x()),
				  *cv.exact_lambda_z, *cv.exact_lambda_x, *cv.exact_lambda_y,
				  *cv.exact_denominator);
			}
			else if (plane == 2) // xy
			{
				// first xy, then z
				lambda3d_SSI_exact<ET>(
				  ET(A().x()), ET(A().y()), ET(A().z()), ET(B().x()), ET(B().y()),
				  ET(B().z()), ET(P().x()), ET(P().y()), ET(Q().x()), ET(Q().y()),
				  *cv.exact_lambda_x, *cv.exact_lambda_y, *cv.exact_lambda_z,
				  *cv.exact_denominator);
			}
			else
			{
				OMC_ASSERT(false, "plane not initialized.");
			}

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
		if (plane == 0) // yz
		{
			// first yz, then x
			lambda3d_SSI_exact<ET>(ET(A().y()), ET(A().z()), ET(A().x()), ET(B().y()),
			                       ET(B().z()), ET(B().x()), ET(P().y()), ET(P().z()),
			                       ET(Q().y()), ET(Q().z()), ly, lz, lx, d);
		}
		else if (plane == 1) // zx
		{
			// first zx, then y
			lambda3d_SSI_exact<ET>(ET(A().z()), ET(A().x()), ET(A().y()), ET(B().z()),
			                       ET(B().x()), ET(B().y()), ET(P().z()), ET(P().x()),
			                       ET(Q().z()), ET(Q().x()), lz, lx, ly, d);
		}
		else if (plane == 2) // xy
		{
			// first xy, then z
			lambda3d_SSI_exact<ET>(ET(A().x()), ET(A().y()), ET(A().z()), ET(B().x()),
			                       ET(B().y()), ET(B().z()), ET(P().x()), ET(P().y()),
			                       ET(Q().x()), ET(Q().y()), lx, ly, lz, d);
		}
		else
		{
			OMC_ASSERT(false, "plane not initialized.");
		}
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
void ImplicitPoint3T_SSI<IT, ET>::getExpansionLambda(FT **lx, int &lx_len,
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
			if (plane == 0) // YZ
			{
				// first yz, then x
				lambda3d_SSI_expansion(
				  A().y(), A().z(), A().x(), B().y(), B().z(), B().x(), P().y(),
				  P().z(), Q().y(), Q().z(), &cv.expansion_lambda_y,
				  cv.expansion_lambda_y_len, &cv.expansion_lambda_z,
				  cv.expansion_lambda_z_len, &cv.expansion_lambda_x,
				  cv.expansion_lambda_x_len, &cv.expansion_denominator,
				  cv.expansion_d_len);
			}
			else if (plane == 1) // ZX
			{
				// first zx, then y
				lambda3d_SSI_expansion(
				  A().z(), A().x(), A().y(), B().z(), B().x(), B().y(), P().z(),
				  P().x(), Q().z(), Q().x(), &cv.expansion_lambda_z,
				  cv.expansion_lambda_z_len, &cv.expansion_lambda_x,
				  cv.expansion_lambda_x_len, &cv.expansion_lambda_y,
				  cv.expansion_lambda_y_len, &cv.expansion_denominator,
				  cv.expansion_d_len);
			}
			else if (plane == 2) // XY
			{
				// first xy, then z
				lambda3d_SSI_expansion(
				  A().x(), A().y(), A().z(), B().x(), B().y(), B().z(), P().x(),
				  P().y(), Q().x(), Q().y(), &cv.expansion_lambda_x,
				  cv.expansion_lambda_x_len, &cv.expansion_lambda_y,
				  cv.expansion_lambda_y_len, &cv.expansion_lambda_z,
				  cv.expansion_lambda_z_len, &cv.expansion_denominator,
				  cv.expansion_d_len);
			}

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
		if (plane == 0) // YZ
		{
			// first yz, then x
			lambda3d_SSI_expansion(A().y(), A().z(), A().x(), B().y(), B().z(),
			                       B().x(), P().y(), P().z(), Q().y(), Q().z(), ly,
			                       ly_len, lz, lz_len, lx, lx_len, d, d_len);
		}
		else if (plane == 1) // ZX
		{
			// first zx, then y
			lambda3d_SSI_expansion(A().z(), A().x(), A().y(), B().z(), B().x(),
			                       B().y(), P().z(), P().x(), Q().z(), Q().x(), lz,
			                       lz_len, lx, lx_len, ly, ly_len, d, d_len);
		}
		else if (plane == 2) // XY
		{
			// first xy, then z
			lambda3d_SSI_expansion(A().x(), A().y(), A().z(), B().x(), B().y(),
			                       B().z(), P().x(), P().y(), Q().x(), Q().y(), lx,
			                       lx_len, ly, ly_len, lz, lz_len, d, d_len);
		}
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
void ImplicitPoint3T_SSI<IT, ET>::get_Explicit(EP &e) const
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
auto ImplicitPoint3T_SSI<IT, ET>::to_Explicit() const -> EP
{
	EP e;
	get_Explicit(e);
	return e;
}

template <typename IT, typename ET>
bool ImplicitPoint3T_SSI<IT, ET>::getFilteredLambda(FT &lx, FT &ly, FT &lz,
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
			if (plane == 0) // YZ
			{
				// first yz, then x
				if (!lambda3d_SSI_filtered(A().y(), A().z(), A().x(), B().y(), B().z(),
				                           B().x(), P().y(), P().z(), Q().y(), Q().z(),
				                           cv.ssfilter_denominator,
				                           cv.ssfilter_lambda_y, cv.ssfilter_lambda_z,
				                           cv.ssfilter_lambda_x, cv.ssfilter_max_val))
					cv.ssfilter_denominator = 0;
			}
			else if (plane == 1) // ZX
			{
				// first zx, then y
				if (!lambda3d_SSI_filtered(A().z(), A().x(), A().y(), B().z(), B().x(),
				                           B().y(), P().z(), P().x(), Q().z(), Q().x(),
				                           cv.ssfilter_denominator,
				                           cv.ssfilter_lambda_z, cv.ssfilter_lambda_x,
				                           cv.ssfilter_lambda_y, cv.ssfilter_max_val))
					cv.ssfilter_denominator = 0;
			}
			else if (plane == 2) // XY
			{
				// first xy, then z
				if (!lambda3d_SSI_filtered(A().x(), A().y(), A().z(), B().x(), B().y(),
				                           B().z(), P().x(), P().y(), Q().x(), Q().y(),
				                           cv.ssfilter_denominator,
				                           cv.ssfilter_lambda_x, cv.ssfilter_lambda_y,
				                           cv.ssfilter_lambda_z, cv.ssfilter_max_val))
					cv.ssfilter_denominator = 0;
			}
			else
			{
				OMC_ASSERT(false, "plane is not initialized.");
			}

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
		bx = A().x();
		by = A().y();
		bz = A().z();
		if (cv.ssfilter_denominator != 0 && cv.ssfilter_max_val > mv)
			mv = cv.ssfilter_max_val;
		return (d != 0);
	}
	else
	{
		FT mv_ = 0;

		if (plane == 0) // YZ
		{
			// first yz, then x
			if (!lambda3d_SSI_filtered(A().y(), A().z(), A().x(), B().y(), B().z(),
			                           B().x(), P().y(), P().z(), Q().y(), Q().z(), d,
			                           ly, lz, lx, mv_))
				d = 0;
		}
		else if (plane == 1) // ZX
		{
			// first zx, then y
			if (!lambda3d_SSI_filtered(A().z(), A().x(), A().y(), B().z(), B().x(),
			                           B().y(), P().z(), P().x(), Q().z(), Q().x(), d,
			                           lz, lx, ly, mv_))
				d = 0;
		}
		else if (plane == 2) // XY
		{
			// first xy, then z
			if (!lambda3d_SSI_filtered(A().x(), A().y(), A().z(), B().x(), B().y(),
			                           B().z(), P().x(), P().y(), Q().x(), Q().y(), d,
			                           lx, ly, lz, mv_))
				d = 0;
		}
		else
		{
			OMC_ASSERT(false, "plane is not initialized.");
		}

		if (d < 0)
		{
			lx = -lx;
			ly = -ly;
			lz = -lz;
			d  = -d;
		}
		bx = A().x();
		by = A().y();
		bz = A().z();
		if (d != 0 && mv_ > mv)
			mv = mv_;
		return d != 0;
	}
	#else
	lx = m_lx;
	ly = m_ly;
	lz = m_lz;
	d  = m_ld;
	bx = A().x();
	by = A().y();
	bz = A().z();
	if (d != 0 && m_maxvar > mv)
		mv = m_maxvar;
	return d != 0;
	#endif
}

template <typename IT, typename ET>
bool ImplicitPoint3T_SSI<IT, ET>::getIntervalLambda(IT &lx, IT &ly, IT &lz,
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
			if (plane == 0) // YZ
			{
				// first yz, then x
				lambda3d_SSI_interval<IT>(A().y(), A().z(), A().x(), B().y(), B().z(),
				                          B().x(), P().y(), P().z(), Q().y(), Q().z(),
				                          cv.dfilter_denominator, cv.dfilter_lambda_y,
				                          cv.dfilter_lambda_z, cv.dfilter_lambda_x);
			}
			else if (plane == 1) // ZX
			{
				// first zx, then y
				lambda3d_SSI_interval<IT>(A().z(), A().x(), A().y(), B().z(), B().x(),
				                          B().y(), P().z(), P().x(), Q().z(), Q().x(),
				                          cv.dfilter_denominator, cv.dfilter_lambda_z,
				                          cv.dfilter_lambda_x, cv.dfilter_lambda_y);
			}
			else if (plane == 2) // XY
			{
				// first xy, then z
				lambda3d_SSI_interval<IT>(A().x(), A().y(), A().z(), B().x(), B().y(),
				                          B().z(), P().x(), P().y(), Q().x(), Q().y(),
				                          cv.dfilter_denominator, cv.dfilter_lambda_x,
				                          cv.dfilter_lambda_y, cv.dfilter_lambda_z);
			}
			else
			{
				OMC_ASSERT(false, "plane is not initialized.");
			}
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
		bx = A().x();
		by = A().y();
		bz = A().z();
		return (d.is_sign_reliable());
	}
	else
	{
		if (plane == 0) // YZ
		{
			// first yz, then x
			lambda3d_SSI_interval<IT>(A().y(), A().z(), A().x(), B().y(), B().z(),
			                          B().x(), P().y(), P().z(), Q().y(), Q().z(), d,
			                          ly, lz, lx);
		}
		else if (plane == 1) // ZX
		{
			// first zx, then y
			lambda3d_SSI_interval<IT>(A().z(), A().x(), A().y(), B().z(), B().x(),
			                          B().y(), P().z(), P().x(), Q().z(), Q().x(), d,
			                          lz, lx, ly);
		}
		else if (plane == 2) // XY
		{
			// first xy, then z
			lambda3d_SSI_interval<IT>(A().x(), A().y(), A().z(), B().x(), B().y(),
			                          B().z(), P().x(), P().y(), Q().x(), Q().y(), d,
			                          lx, ly, lz);
		}
		else
		{
			OMC_ASSERT(false, "plane is not initialized.");
		}
		if (d.is_negative())
		{
			lx.invert();
			ly.invert();
			lz.invert();
			d.invert();
		}
		bx = A().x();
		by = A().y();
		bz = A().z();
		return d.is_sign_reliable();
	}
}

template <typename IT, typename ET>
void ImplicitPoint3T_SSI<IT, ET>::getExactLambda(ET &lx, ET &ly, ET &lz, ET &d,
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
			if (plane == 0) // yz
			{
				// first yz, then x
				lambda3d_SSI_exact<ET>(
				  ET(A().y()), ET(A().z()), ET(A().x()), ET(B().y()), ET(B().z()),
				  ET(B().x()), ET(P().y()), ET(P().z()), ET(Q().y()), ET(Q().z()),
				  *cv.exact_denominator, *cv.exact_lambda_y, *cv.exact_lambda_z,
				  *cv.exact_lambda_x);
			}
			else if (plane == 1) // zx
			{
				// first zx, then y
				lambda3d_SSI_exact<ET>(
				  ET(A().z()), ET(A().x()), ET(A().y()), ET(B().z()), ET(B().x()),
				  ET(B().y()), ET(P().z()), ET(P().x()), ET(Q().z()), ET(Q().x()),
				  *cv.exact_denominator, *cv.exact_lambda_z, *cv.exact_lambda_x,
				  *cv.exact_lambda_y);
			}
			else if (plane == 2) // xy
			{
				// first xy, then z
				lambda3d_SSI_exact<ET>(
				  ET(A().x()), ET(A().y()), ET(A().z()), ET(B().x()), ET(B().y()),
				  ET(B().z()), ET(P().x()), ET(P().y()), ET(Q().x()), ET(Q().y()),
				  *cv.exact_denominator, *cv.exact_lambda_x, *cv.exact_lambda_y,
				  *cv.exact_lambda_z);
			}
			else
			{
				OMC_ASSERT(false, "plane not initialized.");
			}

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
		bx = A().x();
		by = A().y();
		bz = A().z();
	}
	else
	{
		if (plane == 0) // yz
		{
			// first yz, then x
			lambda3d_SSI_exact<ET>(ET(A().y()), ET(A().z()), ET(A().x()), ET(B().y()),
			                       ET(B().z()), ET(B().x()), ET(P().y()), ET(P().z()),
			                       ET(Q().y()), ET(Q().z()), d, ly, lz, lx);
		}
		else if (plane == 1) // zx
		{
			// first zx, then y
			lambda3d_SSI_exact<ET>(ET(A().z()), ET(A().x()), ET(A().y()), ET(B().z()),
			                       ET(B().x()), ET(B().y()), ET(P().z()), ET(P().x()),
			                       ET(Q().z()), ET(Q().x()), d, lz, lx, ly);
		}
		else if (plane == 2) // xy
		{
			// first xy, then z
			lambda3d_SSI_exact<ET>(ET(A().x()), ET(A().y()), ET(A().z()), ET(B().x()),
			                       ET(B().y()), ET(B().z()), ET(P().x()), ET(P().y()),
			                       ET(Q().x()), ET(Q().y()), d, lx, ly, lz);
		}
		else
		{
			OMC_ASSERT(false, "plane not initialized.");
		}
		if (OMC::sign(d) == Sign::NEGATIVE)
		{
			lx = -lx;
			ly = -ly;
			lz = -lz;
			d  = -d;
		}
		bx = A().x();
		by = A().y();
		bz = A().z();
	}
}

template <typename IT, typename ET>
void ImplicitPoint3T_SSI<IT, ET>::getExpansionLambda(FT **lx, int &lx_len,
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
			if (plane == 0) // YZ
			{
				// first yz, then x
				lambda3d_SSI_expansion(
				  A().y(), A().z(), A().x(), B().y(), B().z(), B().x(), P().y(),
				  P().z(), Q().y(), Q().z(), &cv.expansion_denominator,
				  cv.expansion_d_len, &cv.expansion_lambda_y, cv.expansion_lambda_y_len,
				  &cv.expansion_lambda_z, cv.expansion_lambda_z_len,
				  &cv.expansion_lambda_x, cv.expansion_lambda_x_len);
			}
			else if (plane == 1) // ZX
			{
				// first zx, then y
				lambda3d_SSI_expansion(
				  A().z(), A().x(), A().y(), B().z(), B().x(), B().y(), P().z(),
				  P().x(), Q().z(), Q().x(), &cv.expansion_denominator,
				  cv.expansion_d_len, &cv.expansion_lambda_z, cv.expansion_lambda_z_len,
				  &cv.expansion_lambda_x, cv.expansion_lambda_x_len,
				  &cv.expansion_lambda_y, cv.expansion_lambda_y_len);
			}
			else if (plane == 2) // XY
			{
				// first xy, then z
				lambda3d_SSI_expansion(
				  A().x(), A().y(), A().z(), B().x(), B().y(), B().z(), P().x(),
				  P().y(), Q().x(), Q().y(), &cv.expansion_denominator,
				  cv.expansion_d_len, &cv.expansion_lambda_x, cv.expansion_lambda_x_len,
				  &cv.expansion_lambda_y, cv.expansion_lambda_y_len,
				  &cv.expansion_lambda_z, cv.expansion_lambda_z_len);
			}
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
		bx     = A().x();
		by     = A().y();
		bz     = A().z();
	}
	else
	{
		if (plane == 0) // YZ
		{
			// first yz, then x
			lambda3d_SSI_expansion(A().y(), A().z(), A().x(), B().y(), B().z(),
			                       B().x(), P().y(), P().z(), Q().y(), Q().z(), d,
			                       d_len, ly, ly_len, lz, lz_len, lx, lx_len);
		}
		else if (plane == 1) // ZX
		{
			// first zx, then y
			lambda3d_SSI_expansion(A().z(), A().x(), A().y(), B().z(), B().x(),
			                       B().y(), P().z(), P().x(), Q().z(), Q().x(), d,
			                       d_len, lz, lz_len, lx, lx_len, ly, ly_len);
		}
		else if (plane == 2) // XY
		{
			// first xy, then z
			lambda3d_SSI_expansion(A().x(), A().y(), A().z(), B().x(), B().y(),
			                       B().z(), P().x(), P().y(), Q().x(), Q().y(), d,
			                       d_len, lx, lx_len, ly, ly_len, lz, lz_len);
		}
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
		bx = A().x();
		by = A().y();
		bz = A().z();
	}
}

#endif

/// @brief This is just a function for profiling. It outputs maxvar under
/// indirect predicates to compare with other predicates
template <typename IT, typename ET>
auto ImplicitPoint3T_SSI<IT, ET>::getIndirectMaxVar() const -> FT
{
	FT xba = B().x() - A().x();
	FT yba = B().y() - A().y();
	FT zba = B().z() - A().z();
	FT xap = A().x() - P().x();
	FT yap = A().y() - P().y();
	FT yqp = Q().y() - P().y();
	FT xqp = Q().x() - P().x();

	FT _tmp_fabs, max_var = 0.;
	if ((_tmp_fabs = fabs(A().x())) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(A().y())) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(A().z())) > max_var)
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

	return max_var;
}

/// @brief This is just a function for profiling. It outputs maxvar under offset
/// predicates to compare with other predicates
template <typename IT, typename ET>
auto ImplicitPoint3T_SSI<IT, ET>::getOffsetMaxVar() const -> FT
{
	FT xap = A().x() - P().x();
	FT yap = A().y() - P().y();
	FT yqp = Q().y() - P().y();
	FT xqp = Q().x() - P().x();
	FT xba = B().x() - A().x();
	FT yba = B().y() - A().y();
	FT zba = B().z() - A().z();

	FT _tmp_fabs, max_var = 0.;
	if ((_tmp_fabs = fabs(xap)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yap)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yqp)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xqp)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(xba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(yba)) > max_var)
		max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(zba)) > max_var)
		max_var = _tmp_fabs;

	return max_var;
}
} // namespace OMC