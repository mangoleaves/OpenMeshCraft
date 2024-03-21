#pragma once

#include "IndirectPredicateDetails.h"
#include "IndirectPredicateDetailsHand.h"

#include "OpenMeshCraft/Utils/Exception.h"

#include <bitset>

namespace OMC {

#define NF IT, ET        // No semi-static Filter
#define SSF IT, ET, true // with Semi-Static Filter
#define DF IT, ET, false // skip semi-static filter, use Dynamic Filter
#define PntType(p) static_cast<uint32_t>(p.point_type())

template <typename FT, typename IT, typename ET>
Sign DotProductSign2D_Indirect<FT, IT, ET>::operator()(const PointT &p,
                                                       const PointT &r,
                                                       const PointT &q)
{
	PntArr2 arr = get_pnts_arr2<true>(PntType(p), PntType(r), PntType(q));

	// cases are simple, expand them by hand.
	switch (arr)
	{
	case PntArr2::EEE:
		return dotProductSign2D<SSF>(p, r, q);
	case PntArr2::EEI:
		return dotProductSign2D_EEI<NF>(p, r, q);
	case PntArr2::EIE:
		return dotProductSign2D_IEE<NF>(r, p, q);
	case PntArr2::EII:
		return dotProductSign2D_IEI<NF>(r, p, q);
	case PntArr2::IEE:
		return dotProductSign2D_IEE<NF>(p, r, q);
	case PntArr2::IEI:
		return dotProductSign2D_IEI<NF>(p, r, q);
	case PntArr2::IIE:
		return dotProductSign2D_IIE<NF>(p, r, q);
	default: // PntArr2::III
		return dotProductSign2D_III<NF>(p, r, q);
	}
}

template <typename FT, typename IT, typename ET>
Sign DotProductSign2D_Indirect<FT, IT, ET>::operator()(const PointT &p,
                                                       const PointT &r,
                                                       const PointT &q,
                                                       const PointT &s)
{
	PntArr2 arr =
	  get_pnts_arr2<true>(PntType(p), PntType(r), PntType(q), PntType(s));
	// cases are simple, expand them by hand.
	switch (arr)
	{
	case PntArr2::EEEE:
		return dotProductSign2D4P<SSF>(p, r, q, s);
	default:
		OMC_EXIT("DotProductSign2D4P - should not happen");
		return Sign::ZERO; // warning killer
	}
}

template <typename FT, typename IT, typename ET>
Sign DotProductSign3D_Indirect<FT, IT, ET>::operator()(const PointT &p,
                                                       const PointT &r,
                                                       const PointT &q)
{
	PntArr3 arr = get_pnts_arr3<true>(PntType(p), PntType(r), PntType(q));
	// cases are simple, expand them by hand.
	switch (arr)
	{
	case PntArr3::EEE:
		return dotProductSign3D<SSF>(p, r, q);
	case PntArr3::EEI:
		return dotProductSign3D_EEI<NF>(p, r, q);
	case PntArr3::EIE:
		return dotProductSign3D_IEE<NF>(r, p, q);
	case PntArr3::EII:
		return dotProductSign3D_IEI<NF>(r, p, q);
	case PntArr3::IEE:
		return dotProductSign3D_IEE<NF>(p, r, q);
	case PntArr3::IEI:
		return dotProductSign3D_IEI<NF>(p, r, q);
	case PntArr3::IIE:
		return dotProductSign3D_IIE<NF>(p, r, q);
	default: // PntArr3::III
		return dotProductSign3D_III<NF>(p, r, q);
	}
}

template <typename FT, typename IT, typename ET>
Sign DotProductSign3D_Indirect<FT, IT, ET>::operator()(const PointT &p,
                                                       const PointT &r,
                                                       const PointT &q,
                                                       const PointT &s)
{
	PntArr3 arr =
	  get_pnts_arr3<true>(PntType(p), PntType(r), PntType(q), PntType(s));
	// cases are simple, expand them by hand.
	switch (arr)
	{
	case PntArr3::EEEE:
		return dotProductSign3D4P<SSF>(p, r, q, s);
	default:
		OMC_EXIT("DotProductSign3D4P - should not happen");
		return Sign::ZERO; // warning killer
	}
}

template <typename FT, typename IT, typename ET>
Sign DotProductSignOn2D_Indirect<FT, IT, ET>::on_xy(const PointT &p,
                                                    const PointT &r,
                                                    const PointT &q)
{
	OMC_EXIT("DotProductSignOn2Dxy - should not happen");
	return Sign::ZERO; // warning killer
}

template <typename FT, typename IT, typename ET>
Sign DotProductSignOn2D_Indirect<FT, IT, ET>::on_xy(const PointT &p,
                                                    const PointT &r,
                                                    const PointT &q,
                                                    const PointT &s)
{
	PntArr3 arr =
	  get_pnts_arr3<true>(PntType(p), PntType(r), PntType(q), PntType(s));
	// cases are simple, expand them by hand.
	switch (arr)
	{
	case PntArr3::EEEE:
		return dotProductSignOn2Dxy4P<SSF>(p, r, q, s);
	default:
		OMC_EXIT("DotProductSignOn2Dxy4P - should not happen");
		return Sign::ZERO; // warning killer
	}
}

template <typename FT, typename IT, typename ET>
Sign DotProductSignOn2D_Indirect<FT, IT, ET>::on_yz(const PointT &p,
                                                    const PointT &r,
                                                    const PointT &q)
{
	OMC_EXIT("DotProductSignOn2Dyz - should not happen");
	return Sign::ZERO; // warning killer
}

template <typename FT, typename IT, typename ET>
Sign DotProductSignOn2D_Indirect<FT, IT, ET>::on_yz(const PointT &p,
                                                    const PointT &r,
                                                    const PointT &q,
                                                    const PointT &s)
{
	PntArr3 arr =
	  get_pnts_arr3<true>(PntType(p), PntType(r), PntType(q), PntType(s));
	// cases are simple, expand them by hand.
	switch (arr)
	{
	case PntArr3::EEEE:
		return dotProductSignOn2Dyz4P<SSF>(p, r, q, s);
	default:
		OMC_EXIT("DotProductSignOn2Dyz4P - should not happen");
		return Sign::ZERO; // warning killer
	}
}

template <typename FT, typename IT, typename ET>
Sign DotProductSignOn2D_Indirect<FT, IT, ET>::on_zx(const PointT &p,
                                                    const PointT &r,
                                                    const PointT &q)
{
	OMC_EXIT("DotProductSignOn2Dzx - should not happen");
	return Sign::ZERO; // warning killer
}

template <typename FT, typename IT, typename ET>
Sign DotProductSignOn2D_Indirect<FT, IT, ET>::on_zx(const PointT &p,
                                                    const PointT &r,
                                                    const PointT &q,
                                                    const PointT &s)
{
	PntArr3 arr =
	  get_pnts_arr3<true>(PntType(p), PntType(r), PntType(q), PntType(s));
	// cases are simple, expand them by hand.
	switch (arr)
	{
	case PntArr3::EEEE:
		return dotProductSignOn2Dzx4P<SSF>(p, r, q, s);
	default:
		OMC_EXIT("DotProductSignOn2Dzx4P - should not happen");
		return Sign::ZERO; // warning killer
	}
}

template <typename FT, typename IT, typename ET>
Sign Orient2D_Indirect<FT, IT, ET>::operator()(const PointT &p, const PointT &q,
                                               const PointT &query)
{
	PntArr2 arr = get_pnts_arr2<false>(PntType(p), PntType(q), PntType(query));
	// Orient2D is simple now, expand it by hand.
	switch (arr)
	{
	case PntArr2::EEE:
		return orient2d<IT, ET>(p, q, query);
	case PntArr2::EES:
		return orient2D_IEE<SSF>(query, p, q, PntArr2::S);
	case PntArr2::ESE:
		return orient2D_IEE<SSF>(q, query, p, PntArr2::S);
	case PntArr2::SEE:
		return orient2D_IEE<SSF>(p, q, query, PntArr2::S);
	case PntArr2::ESS:
		return orient2D_IIE<SSF>(q, query, p, PntArr2::SS);
	case PntArr2::SES:
		return orient2D_IIE<SSF>(query, p, q, PntArr2::SS);
	case PntArr2::SSE:
		return orient2D_IIE<SSF>(p, q, query, PntArr2::SS);
	case PntArr2::SSS:
		return orient2D_III<SSF>(p, q, query, PntArr2::SSS);
	default:
		OMC_EXIT("Orient2D - should not happen");
	}
}

template <typename FT, typename IT, typename ET>
Sign Orient2D_Indirect<FT, IT, ET>::operator()(const FT *p, const FT *q,
                                               const FT *query)
{
	return orient2d(p, q, query);
}

template <typename FT, typename IT, typename ET>
Sign SquareDistance2D_Indirect<FT, IT, ET>::operator()(const PointT &p,
                                                       const PointT &q,
                                                       FT            sqr_dis)
{
	PntArr2 arr = get_pnts_arr2<true>(PntType(p), PntType(q));
	// SquareDistance is simple now, expand it by hand.
	switch (arr)
	{
	case PntArr2::EE:
		return squareDistance2D<SSF>(p, q, sqr_dis);
	case PntArr2::EI:
		return squareDistance2D_IE<NF>(q, p, sqr_dis);
	case PntArr2::IE:
		return squareDistance2D_IE<NF>(p, q, sqr_dis);
	case PntArr2::II:
		return squareDistance2D_II<NF>(q, p, sqr_dis);
	default:
		OMC_ASSERT(false, "SquaredDistance2D - should not happen");
		return Sign::ZERO; // warning killer
	}
}

template <typename FT, typename IT, typename ET>
Sign SquareDistance3D_Indirect<FT, IT, ET>::operator()(const PointT &p,
                                                       const PointT &q,
                                                       FT            sqr_dis)
{
	PntArr3 arr = get_pnts_arr3<true>(PntType(p), PntType(q));
	// SquareDistance is simple now, expand it by hand.
	switch (arr)
	{
	case PntArr3::EE:
		return squareDistance3D<SSF>(p, q, sqr_dis);
	case PntArr3::EI:
		return squareDistance3D_IE<NF>(q, p, sqr_dis);
	case PntArr3::IE:
		return squareDistance3D_IE<NF>(p, q, sqr_dis);
	case PntArr3::II:
		return squareDistance3D_II<NF>(q, p, sqr_dis);
	default:
		OMC_ASSERT(false, "SquaredDistance3D - should not happen");
		return Sign::ZERO; // warning killer
	}
}

template <typename FT, typename IT, typename ET>
Sign InCircle_Indirect<FT, IT, ET>::operator()(const PointT &a, const PointT &b,
                                               const PointT &c, const PointT &d)
{
	int i = a.is_Explicit() + b.is_Explicit() + c.is_Explicit() + d.is_Explicit();
	// InCircle is simple now, expand it by hand.
	if (i == 4)
		return inCircle<SSF>(a, b, c, d);
	if (i == 3)
	{
		if (!a.is_Explicit())
			return inCircle_IEEE<NF>(a, b, c, d);
		if (!b.is_Explicit())
			return inCircle_IEEE<NF>(b, c, a, d);
		if (!c.is_Explicit())
			return inCircle_IEEE<NF>(c, d, a, b);
		return inCircle_IEEE<NF>(d, a, c, b);
	}
	if (i == 2)
	{
		if (c.is_Explicit() && d.is_Explicit())
			return inCircle_IIEE<NF>(a, b, c, d);
		if (b.is_Explicit() && d.is_Explicit())
			return inCircle_IIEE<NF>(a, c, d, b);
		if (a.is_Explicit() && d.is_Explicit())
			return inCircle_IIEE<NF>(b, c, a, d);
		if (b.is_Explicit() && c.is_Explicit())
			return inCircle_IIEE<NF>(d, a, c, b);
		if (a.is_Explicit() && c.is_Explicit())
			return inCircle_IIEE<NF>(d, b, a, c);
		return inCircle_IIEE<NF>(c, d, a, b);
	}
	if (i == 1)
	{
		if (d.is_Explicit())
			return inCircle_IIIE<NF>(a, b, c, d);
		if (c.is_Explicit())
			return inCircle_IIIE<NF>(d, b, a, c);
		if (b.is_Explicit())
			return inCircle_IIIE<NF>(a, c, d, b);
		return inCircle_IIIE<NF>(b, d, c, a);
	}
	return inCircle_IIII<NF>(a, b, c, d);
}

template <typename FT, typename IT, typename ET>
Sign InCircle_Indirect<FT, IT, ET>::operator()(const FT *a, const FT *b,
                                               const FT *c, const FT *d)
{
	// TODO optimize by shewchuk predicate
	return inCircle<SSF>(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]);
}

template <typename FT, typename IT, typename ET>
Sign Orient3D_Indirect<FT, IT, ET>::operator()(const PointT &a, const PointT &b,
                                               const PointT &c, const PointT &d)
{
	// Here we implicitly assume that points are 3D. Do not check.
	// clang-format off
	int i = a.is_Explicit() + b.is_Explicit() + c.is_Explicit() + d.is_Explicit();

	if (i == 4) return orient3d<IT, ET>(a, b, c, d);

	bool has_ssf = a.has_ssf() && b.has_ssf() && c.has_ssf() && d.has_ssf();

	if (has_ssf)
	{
		std::array<uint32_t, 4> pos{0,1,2,3}, types{PntType(a), PntType(b), PntType(c), PntType(d)};
		uint32_t swap_cnt;
		PntArr3 arr = sort_pnts_arr3(types, pos, swap_cnt);
		const PointT *pnts[4] = {&a, &b, &c, &d};
		Sign          sign    = Sign::ZERO;
		if (i == 3) sign = orient3D_IEEE<SSF>(*pnts[pos[0]], *pnts[pos[1]], *pnts[pos[2]], *pnts[pos[3]], arr);
		else if (i == 2) sign = orient3D_IIEE<SSF>(*pnts[pos[0]], *pnts[pos[1]], *pnts[pos[2]], *pnts[pos[3]], arr);
		else if (i == 1) sign = orient3D_IIIE<SSF>(*pnts[pos[0]], *pnts[pos[1]], *pnts[pos[2]], *pnts[pos[3]], arr);
		else sign = orient3D_IIII<SSF>(*pnts[pos[0]], *pnts[pos[1]], *pnts[pos[2]], *pnts[pos[3]], arr);

		return swap_cnt % 2 == 1 ? reverse_sign(sign) : sign;
	}
	else
	{
		if (i == 3)
		{
			if (!a.is_Explicit()) return orient3D_IEEE<DF>(a, b, c, d, PntArr3::I);
			if (!b.is_Explicit()) return orient3D_IEEE<DF>(b, c, a, d, PntArr3::I);
			if (!c.is_Explicit()) return orient3D_IEEE<DF>(c, d, a, b, PntArr3::I);
			/*if (!d.is_Explicit())*/ return orient3D_IEEE<DF>(d, a, c, b, PntArr3::I);
		}
		else if (i == 2)
		{
			if (c.is_Explicit() && d.is_Explicit()) return orient3D_IIEE<DF>(a, b, c, d, PntArr3::II);
			if (b.is_Explicit() && d.is_Explicit()) return orient3D_IIEE<DF>(a, c, d, b, PntArr3::II);
			if (a.is_Explicit() && d.is_Explicit()) return orient3D_IIEE<DF>(b, c, a, d, PntArr3::II);
			if (b.is_Explicit() && c.is_Explicit()) return orient3D_IIEE<DF>(d, a, c, b, PntArr3::II);
			if (a.is_Explicit() && c.is_Explicit()) return orient3D_IIEE<DF>(d, b, a, c, PntArr3::II);
			/*if (a.is_Explicit() && b.is_Explicit())*/ return orient3D_IIEE<DF>(c, d, a, b, PntArr3::II);
		}
		else if (i == 1)
		{
			if (d.is_Explicit()) return orient3D_IIIE<DF>(a, b, c, d, PntArr3::III);
			if (c.is_Explicit()) return orient3D_IIIE<DF>(d, b, a, c, PntArr3::III);
			if (b.is_Explicit()) return orient3D_IIIE<DF>(a, c, d, b, PntArr3::III);
			/*if (a.is_Explicit())*/ return orient3D_IIIE<DF>(b, d, c, a, PntArr3::III);
		}
		else
		{
			return orient3D_IIII<DF>(a, b, c, d, PntArr3::IIII);
		}
	}
	// clang-format on
}

template <typename FT, typename IT, typename ET>
Sign Orient3D_Indirect<FT, IT, ET>::operator()(const FT *a, const FT *b,
                                               const FT *c, const FT *d)
{
	return orient3d(a, b, c, d);
}

template <typename FT, typename IT, typename ET>
void Orient3D_Indirect<FT, IT, ET>::get_minors(const FT *a, const FT *b,
                                               const FT *c, FT *minor, FT *perm)
{
	orient3d_get_minors(a, b, c, minor, perm);
}

template <typename FT, typename IT, typename ET>
Sign Orient3D_Indirect<FT, IT, ET>::with_cached_minors(
  const FT *pa, const FT *pb, const FT *pc, const FT *pd, const FT *minor,
  const FT *perm)
{
	return orient3d_with_cached_minors(pa, pb, pc, pd, minor, perm);
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::operator()(const PointT &a,
                                                 const PointT &b,
                                                 const PointT &c, int n_max)
{
	if (n_max == 0)
		return on_yz(a, b, c);
	else if (n_max == 1)
		return on_zx(a, b, c);
	else
		return on_xy(a, b, c);
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::on_xy(const PointT &a, const PointT &b,
                                            const PointT &c)
{
	if (a.is_Explicit() && b.is_Explicit() && c.is_Explicit())
		return orient2dxy<IT, ET>(a, b, c);

	if (a.has_ssf() && b.has_ssf() && c.has_ssf())
	{
		int i = a.is_Explicit() + b.is_Explicit() + c.is_Explicit();
		if (i == 2)
		{ // cases are simple, expand them by hand.
			// clang-format off
			if (a.is_Implicit()) return orientOn2Dxy_IEE<SSF>(a, b, c, static_cast<PntArr3>(PntType(a)));
			if (b.is_Implicit()) return orientOn2Dxy_IEE<SSF>(b, c, a, static_cast<PntArr3>(PntType(b)));
			if (c.is_Implicit()) return orientOn2Dxy_IEE<SSF>(c, a, b, static_cast<PntArr3>(PntType(c)));
			// clang-format on
		}
		if (i == 1)
		{ // cases are simple, expand them by hand.
			// clang-format off
			if (a.is_Explicit())
			{
				if (b.point_type() <= c.point_type()) return orientOn2Dxy_IIE<SSF>(b, c, a, get_pnts_arr3<false>(PntType(b), PntType(c)));
				else return reverse_sign(orientOn2Dxy_IIE<SSF>( c, b, a, get_pnts_arr3<false>(PntType(c), PntType(b))));
			}
			if (b.is_Explicit())
			{
				if (a.point_type() <= c.point_type()) return reverse_sign(orientOn2Dxy_IIE<SSF>( a, c, b, get_pnts_arr3<false>(PntType(a), PntType(c))));
				else return orientOn2Dxy_IIE<SSF>( c, a, b, get_pnts_arr3<false>(PntType(c), PntType(a)));
			}
			if (c.is_Explicit())
			{
				if (a.point_type() <= b.point_type()) return orientOn2Dxy_IIE<SSF>( a, b, c, get_pnts_arr3<false>(PntType(a), PntType(b)));
				else return reverse_sign(orientOn2Dxy_IIE<SSF>( b, a, c, get_pnts_arr3<false>(PntType(b), PntType(a))));
			}
			// clang-format on
		}

		// cases are complex (S/L/T 27 cases), use arrangement.
		// clang-format off
		std::array<uint32_t, 3> pos{0, 1, 2}, types{PntType(a), PntType(b), PntType(c)};
		uint32_t      swap_cnt;
		PntArr3       arr     = sort_pnts_arr3(types, pos, swap_cnt);
		const PointT *pnts[3] = {&a, &b, &c};
		Sign          sign    = orientOn2Dxy_III<SSF>(*pnts[pos[0]], *pnts[pos[1]], *pnts[pos[2]], arr);
		// clang-format on
		return swap_cnt % 2 == 0 ? sign : reverse_sign(sign);
	}
	else
	{
		PntArr3 arr = get_pnts_arr3<true>(PntType(a), PntType(b), PntType(c));
		// cases are simple, expand them by hand.
		switch (arr)
		{
		case PntArr3::IEE:
			return orientOn2Dxy_IEE<DF>(a, b, c, PntArr3::I);
		case PntArr3::EIE:
			return orientOn2Dxy_IEE<DF>(b, c, a, PntArr3::I);
		case PntArr3::EEI:
			return orientOn2Dxy_IEE<DF>(c, a, b, PntArr3::I);
		case PntArr3::IIE:
			return orientOn2Dxy_IIE<DF>(a, b, c, PntArr3::II);
		case PntArr3::IEI:
			return orientOn2Dxy_IIE<DF>(c, a, b, PntArr3::II);
		case PntArr3::EII:
			return orientOn2Dxy_IIE<DF>(b, c, a, PntArr3::II);
		default: // PntArr3::III
			return orientOn2Dxy_III<DF>(a, b, c, PntArr3::III);
		}
	}
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::on_yz(const PointT &a, const PointT &b,
                                            const PointT &c)
{
	if (a.is_Explicit() && b.is_Explicit() && c.is_Explicit())
		return orient2dyz<IT, ET>(a, b, c);

	if (a.has_ssf() && b.has_ssf() && c.has_ssf())
	{
		int i = a.is_Explicit() + b.is_Explicit() + c.is_Explicit();
		if (i == 2)
		{ // cases are simple, expand them by hand.
			// clang-format off
			if (a.is_Implicit()) return orientOn2Dyz_IEE<SSF>(a, b, c, static_cast<PntArr3>(PntType(a)));
			if (b.is_Implicit()) return orientOn2Dyz_IEE<SSF>(b, c, a, static_cast<PntArr3>(PntType(b)));
			if (c.is_Implicit()) return orientOn2Dyz_IEE<SSF>(c, a, b, static_cast<PntArr3>(PntType(c)));
			// clang-format on
		}
		if (i == 1)
		{ // cases are simple, expand them by hand.
			// clang-format off
			if (a.is_Explicit())
			{
				if (b.point_type() <= c.point_type()) return orientOn2Dyz_IIE<SSF>( b, c, a, get_pnts_arr3<false>(PntType(b), PntType(c)));
				else return reverse_sign(orientOn2Dyz_IIE<SSF>( c, b, a, get_pnts_arr3<false>(PntType(c), PntType(b))));
			}
			if (b.is_Explicit())
			{
				if (a.point_type() <= c.point_type()) return reverse_sign(orientOn2Dyz_IIE<SSF>( a, c, b, get_pnts_arr3<false>(PntType(a), PntType(c))));
				else return orientOn2Dyz_IIE<SSF>( c, a, b, get_pnts_arr3<false>(PntType(c), PntType(a)));
			}
			if (c.is_Explicit())
			{
				if (a.point_type() <= b.point_type()) return orientOn2Dyz_IIE<SSF>( a, b, c, get_pnts_arr3<false>(PntType(a), PntType(b)));
				else return reverse_sign(orientOn2Dyz_IIE<SSF>( b, a, c, get_pnts_arr3<false>(PntType(b), PntType(a))));
			}
			// clang-format on
		}

		// cases are complex (S/L/T 27 cases), use arrangement.
		// clang-format off
		std::array<uint32_t, 3> pos{0, 1, 2}, types{PntType(a), PntType(b), PntType(c)};
		uint32_t      swap_cnt;
		PntArr3       arr     = sort_pnts_arr3(types, pos, swap_cnt);
		const PointT *pnts[3] = {&a, &b, &c};
		Sign          sign    = orientOn2Dyz_III<SSF>(*pnts[pos[0]], *pnts[pos[1]], *pnts[pos[2]], arr);
		// clang-format on
		return swap_cnt % 2 == 0 ? sign : reverse_sign(sign);
	}
	else
	{
		PntArr3 arr = get_pnts_arr3<true>(PntType(a), PntType(b), PntType(c));
		// cases are simple, expand them by hand.
		switch (arr)
		{
		case PntArr3::IEE:
			return orientOn2Dyz_IEE<DF>(a, b, c, PntArr3::I);
		case PntArr3::EIE:
			return orientOn2Dyz_IEE<DF>(b, c, a, PntArr3::I);
		case PntArr3::EEI:
			return orientOn2Dyz_IEE<DF>(c, a, b, PntArr3::I);
		case PntArr3::IIE:
			return orientOn2Dyz_IIE<DF>(a, b, c, PntArr3::II);
		case PntArr3::IEI:
			return orientOn2Dyz_IIE<DF>(c, a, b, PntArr3::II);
		case PntArr3::EII:
			return orientOn2Dyz_IIE<DF>(b, c, a, PntArr3::II);
		default: // PntArr3::III
			return orientOn2Dyz_III<DF>(a, b, c, PntArr3::III);
		}
	}
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::on_zx(const PointT &a, const PointT &b,
                                            const PointT &c)
{
	if (a.is_Explicit() && b.is_Explicit() && c.is_Explicit())
		return orient2dzx<IT, ET>(a, b, c);

	if (a.has_ssf() && b.has_ssf() && c.has_ssf())
	{
		int i = a.is_Explicit() + b.is_Explicit() + c.is_Explicit();
		if (i == 2)
		{ // cases are simple, expand them by hand.
			// clang-format off
			if (a.is_Implicit()) return orientOn2Dzx_IEE<SSF>(a, b, c, static_cast<PntArr3>(PntType(a)));
			if (b.is_Implicit()) return orientOn2Dzx_IEE<SSF>(b, c, a, static_cast<PntArr3>(PntType(b)));
			if (c.is_Implicit()) return orientOn2Dzx_IEE<SSF>(c, a, b, static_cast<PntArr3>(PntType(c)));
			// clang-format on
		}
		if (i == 1)
		{ // cases are simple, expand them by hand.
			// clang-format off
			if (a.is_Explicit())
			{
				if (b.point_type() <= c.point_type()) return orientOn2Dzx_IIE<SSF>( b, c, a, get_pnts_arr3<false>(PntType(b), PntType(c)));
				else return reverse_sign(orientOn2Dzx_IIE<SSF>( c, b, a, get_pnts_arr3<false>(PntType(c), PntType(b))));
			}
			if (b.is_Explicit())
			{
				if (a.point_type() <= c.point_type()) return reverse_sign(orientOn2Dzx_IIE<SSF>( a, c, b, get_pnts_arr3<false>(PntType(a), PntType(c))));
				else return orientOn2Dzx_IIE<SSF>( c, a, b, get_pnts_arr3<false>(PntType(c), PntType(a)));
			}
			if (c.is_Explicit())
			{
				if (a.point_type() <= b.point_type()) return orientOn2Dzx_IIE<SSF>( a, b, c, get_pnts_arr3<false>(PntType(a), PntType(b)));
				else return reverse_sign(orientOn2Dzx_IIE<SSF>( b, a, c, get_pnts_arr3<false>(PntType(b), PntType(a))));
			}
			// clang-format on
		}

		// cases are complex (S/L/T 27 cases), use arrangement.
		// clang-format off
		std::array<uint32_t, 3> pos{0, 1, 2}, types{PntType(a), PntType(b), PntType(c)};
		uint32_t      swap_cnt;
		PntArr3       arr     = sort_pnts_arr3(types, pos, swap_cnt);
		const PointT *pnts[3] = {&a, &b, &c};
		Sign sign             = orientOn2Dzx_III<SSF>(*pnts[pos[0]], *pnts[pos[1]], *pnts[pos[2]], arr);
		return swap_cnt % 2 == 0 ? sign : reverse_sign(sign);
		// clang-format on
	}
	else
	{
		PntArr3 arr = get_pnts_arr3<true>(PntType(a), PntType(b), PntType(c));
		// cases are simple, expand them by hand.
		switch (arr)
		{
		case PntArr3::IEE:
			return orientOn2Dzx_IEE<DF>(a, b, c, PntArr3::I);
		case PntArr3::EIE:
			return orientOn2Dzx_IEE<DF>(b, c, a, PntArr3::I);
		case PntArr3::EEI:
			return orientOn2Dzx_IEE<DF>(c, a, b, PntArr3::I);
		case PntArr3::IIE:
			return orientOn2Dzx_IIE<DF>(a, b, c, PntArr3::II);
		case PntArr3::IEI:
			return orientOn2Dzx_IIE<DF>(c, a, b, PntArr3::II);
		case PntArr3::EII:
			return orientOn2Dzx_IIE<DF>(b, c, a, PntArr3::II);
		default: // PntArr3::III
			return orientOn2Dzx_III<DF>(a, b, c, PntArr3::III);
		}
	}
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::operator()(const FT *a, const FT *b,
                                                 const FT *c, int n_max)
{
	if (n_max == 0)
		return on_yz(a, b, c);
	else if (n_max == 1)
		return on_zx(a, b, c);
	else
		return on_xy(a, b, c);
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::on_xy(const FT *a, const FT *b,
                                            const FT *c)
{
	return orient2dxy(a, b, c);
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::on_yz(const FT *a, const FT *b,
                                            const FT *c)
{
	return orient2dyz(a, b, c);
}

template <typename FT, typename IT, typename ET>
Sign OrientOn2D_Indirect<FT, IT, ET>::on_zx(const FT *a, const FT *b,
                                            const FT *c)
{
	return orient2dzx(a, b, c);
}

template <typename FT, typename IT, typename ET>
Sign LessThan3D_Indirect<FT, IT, ET>::on_x(const PointT &a, const PointT &b)
{
	if (a.is_Explicit() && b.is_Explicit())
		return static_cast<Sign>(((a.x() > b.x()) - (a.x() < b.x())));

	if (a.has_ssf() && b.has_ssf())
	{
		// cases are simple, expand them by hand.
		if (a.is_SSI() && b.is_Explicit())
			return lessThanOnX_IE<SSF>(a, b, PntArr3::S);
		if (a.is_LPI() && b.is_Explicit())
			return lessThanOnX_IE<SSF>(a, b, PntArr3::L);
		if (a.is_TPI() && b.is_Explicit())
			return lessThanOnX_IE<SSF>(a, b, PntArr3::T);
		if (a.is_Explicit() && b.is_SSI())
			return reverse_sign(lessThanOnX_IE<SSF>(b, a, PntArr3::S));
		if (a.is_Explicit() && b.is_LPI())
			return reverse_sign(lessThanOnX_IE<SSF>(b, a, PntArr3::L));
		if (a.is_Explicit() && b.is_TPI())
			return reverse_sign(lessThanOnX_IE<SSF>(b, a, PntArr3::T));
		if (a.is_SSI() && b.is_SSI())
			return lessThanOnX_II<SSF>(a, b, PntArr3::SS);
		if (a.is_SSI() && b.is_LPI())
			return lessThanOnX_II<SSF>(a, b, PntArr3::SL);
		if (a.is_SSI() && b.is_TPI())
			return lessThanOnX_II<SSF>(a, b, PntArr3::ST);
		if (a.is_LPI() && b.is_SSI())
			return reverse_sign(lessThanOnX_II<SSF>(b, a, PntArr3::SL));
		if (a.is_LPI() && b.is_LPI())
			return lessThanOnX_II<SSF>(a, b, PntArr3::LL);
		if (a.is_LPI() && b.is_TPI())
			return lessThanOnX_II<SSF>(a, b, PntArr3::LT);
		if (a.is_TPI() && b.is_SSI())
			return reverse_sign(lessThanOnX_II<SSF>(b, a, PntArr3::ST));
		if (a.is_TPI() && b.is_LPI())
			return reverse_sign(lessThanOnX_II<SSF>(b, a, PntArr3::LT));
		// if (a.is_TPI() && b.is_TPI())
		return lessThanOnX_II<SSF>(a, b, PntArr3::TT);
	}
	else
	{
		// cases are simple, expand them by hand.
		if (!a.is_Explicit() && b.is_Explicit())
			return lessThanOnX_IE<DF>(a, b, PntArr3::I);
		if (a.is_Explicit() && !b.is_Explicit())
			return reverse_sign(lessThanOnX_IE<DF>(b, a, PntArr3::I));
		return lessThanOnX_II<DF>(a, b, PntArr3::II);
	}
}

template <typename FT, typename IT, typename ET>
Sign LessThan3D_Indirect<FT, IT, ET>::on_y(const PointT &a, const PointT &b)
{
	if (a.is_Explicit() && b.is_Explicit())
		return static_cast<Sign>(((a.y() > b.y()) - (a.y() < b.y())));

	if (a.has_ssf() && b.has_ssf())
	{
		// cases are simple, expand them by hand.
		if (a.is_SSI() && b.is_Explicit())
			return lessThanOnY_IE<SSF>(a, b, PntArr3::S);
		if (a.is_LPI() && b.is_Explicit())
			return lessThanOnY_IE<SSF>(a, b, PntArr3::L);
		if (a.is_TPI() && b.is_Explicit())
			return lessThanOnY_IE<SSF>(a, b, PntArr3::T);
		if (a.is_Explicit() && b.is_SSI())
			return reverse_sign(lessThanOnY_IE<SSF>(b, a, PntArr3::S));
		if (a.is_Explicit() && b.is_LPI())
			return reverse_sign(lessThanOnY_IE<SSF>(b, a, PntArr3::L));
		if (a.is_Explicit() && b.is_TPI())
			return reverse_sign(lessThanOnY_IE<SSF>(b, a, PntArr3::T));
		if (a.is_SSI() && b.is_SSI())
			return lessThanOnY_II<SSF>(a, b, PntArr3::SS);
		if (a.is_SSI() && b.is_LPI())
			return lessThanOnY_II<SSF>(a, b, PntArr3::SL);
		if (a.is_SSI() && b.is_TPI())
			return lessThanOnY_II<SSF>(a, b, PntArr3::ST);
		if (a.is_LPI() && b.is_SSI())
			return reverse_sign(lessThanOnY_II<SSF>(b, a, PntArr3::SL));
		if (a.is_LPI() && b.is_LPI())
			return lessThanOnY_II<SSF>(a, b, PntArr3::LL);
		if (a.is_LPI() && b.is_TPI())
			return lessThanOnY_II<SSF>(a, b, PntArr3::LT);
		if (a.is_TPI() && b.is_SSI())
			return reverse_sign(lessThanOnY_II<SSF>(b, a, PntArr3::ST));
		if (a.is_TPI() && b.is_LPI())
			return reverse_sign(lessThanOnY_II<SSF>(b, a, PntArr3::LT));
		// if (a.is_TPI() && b.is_TPI())
		return lessThanOnY_II<SSF>(a, b, PntArr3::TT);
	}
	else
	{
		// cases are simple, expand them by hand.
		if (!a.is_Explicit() && b.is_Explicit())
			return lessThanOnY_IE<DF>(a, b, PntArr3::I);
		if (a.is_Explicit() && !b.is_Explicit())
			return reverse_sign(lessThanOnY_IE<DF>(b, a, PntArr3::I));
		return lessThanOnY_II<DF>(a, b, PntArr3::II);
	}
}

template <typename FT, typename IT, typename ET>
Sign LessThan3D_Indirect<FT, IT, ET>::on_z(const PointT &a, const PointT &b)
{
	if (a.is_Explicit() && b.is_Explicit())
		return static_cast<Sign>(((a.z() > b.z()) - (a.z() < b.z())));

	if (a.has_ssf() && b.has_ssf())
	{
		// cases are simple, expand them by hand.
		if (a.is_SSI() && b.is_Explicit())
			return lessThanOnZ_IE<SSF>(a, b, PntArr3::S);
		if (a.is_LPI() && b.is_Explicit())
			return lessThanOnZ_IE<SSF>(a, b, PntArr3::L);
		if (a.is_TPI() && b.is_Explicit())
			return lessThanOnZ_IE<SSF>(a, b, PntArr3::T);
		if (a.is_Explicit() && b.is_SSI())
			return reverse_sign(lessThanOnZ_IE<SSF>(b, a, PntArr3::S));
		if (a.is_Explicit() && b.is_LPI())
			return reverse_sign(lessThanOnZ_IE<SSF>(b, a, PntArr3::L));
		if (a.is_Explicit() && b.is_TPI())
			return reverse_sign(lessThanOnZ_IE<SSF>(b, a, PntArr3::T));
		if (a.is_SSI() && b.is_SSI())
			return lessThanOnZ_II<SSF>(a, b, PntArr3::SS);
		if (a.is_SSI() && b.is_LPI())
			return lessThanOnZ_II<SSF>(a, b, PntArr3::SL);
		if (a.is_SSI() && b.is_TPI())
			return lessThanOnZ_II<SSF>(a, b, PntArr3::ST);
		if (a.is_LPI() && b.is_SSI())
			return reverse_sign(lessThanOnZ_II<SSF>(b, a, PntArr3::SL));
		if (a.is_LPI() && b.is_LPI())
			return lessThanOnZ_II<SSF>(a, b, PntArr3::LL);
		if (a.is_LPI() && b.is_TPI())
			return lessThanOnZ_II<SSF>(a, b, PntArr3::LT);
		if (a.is_TPI() && b.is_SSI())
			return reverse_sign(lessThanOnZ_II<SSF>(b, a, PntArr3::ST));
		if (a.is_TPI() && b.is_LPI())
			return reverse_sign(lessThanOnZ_II<SSF>(b, a, PntArr3::LT));
		// if (a.is_TPI() && b.is_TPI())
		return lessThanOnZ_II<SSF>(a, b, PntArr3::TT);
	}
	else
	{
		if (!a.is_Explicit() && b.is_Explicit())
			return lessThanOnZ_IE<DF>(a, b, PntArr3::I);
		if (a.is_Explicit() && !b.is_Explicit())
			return reverse_sign(lessThanOnZ_IE<DF>(b, a, PntArr3::I));
		return lessThanOnZ_II<DF>(a, b, PntArr3::II);
	}
}

template <typename FT, typename IT, typename ET>
std::array<Sign, 3> LessThan3D_Indirect<FT, IT, ET>::on_all(const PointT &a,
                                                            const PointT &b)
{
	if (a.is_Explicit() && b.is_Explicit())
		return lessThanOnAll_EE(a.x(), a.y(), a.z(), b.x(), b.y(), b.z());

	if (a.has_ssf() && b.has_ssf())
	{
		// cases are simple, expand them by hand.
		if (a.is_SSI() && b.is_Explicit())
			return lessThanOnAll_IE<SSF>(a, b, PntArr3::S);
		if (a.is_LPI() && b.is_Explicit())
			return lessThanOnAll_IE<SSF>(a, b, PntArr3::L);
		if (a.is_TPI() && b.is_Explicit())
			return lessThanOnAll_IE<SSF>(a, b, PntArr3::T);
		if (a.is_Explicit() && b.is_SSI())
			return reverse_signs(lessThanOnAll_IE<SSF>(b, a, PntArr3::S));
		if (a.is_Explicit() && b.is_LPI())
			return reverse_signs(lessThanOnAll_IE<SSF>(b, a, PntArr3::L));
		if (a.is_Explicit() && b.is_TPI())
			return reverse_signs(lessThanOnAll_IE<SSF>(b, a, PntArr3::T));
		if (a.is_SSI() && b.is_SSI())
			return lessThanOnAll_II<SSF>(a, b, PntArr3::SS);
		if (a.is_SSI() && b.is_LPI())
			return lessThanOnAll_II<SSF>(a, b, PntArr3::SL);
		if (a.is_SSI() && b.is_TPI())
			return lessThanOnAll_II<SSF>(a, b, PntArr3::ST);
		if (a.is_LPI() && b.is_SSI())
			return reverse_signs(lessThanOnAll_II<SSF>(b, a, PntArr3::SL));
		if (a.is_LPI() && b.is_LPI())
			return lessThanOnAll_II<SSF>(a, b, PntArr3::LL);
		if (a.is_LPI() && b.is_TPI())
			return lessThanOnAll_II<SSF>(a, b, PntArr3::LT);
		if (a.is_TPI() && b.is_SSI())
			return reverse_signs(lessThanOnAll_II<SSF>(b, a, PntArr3::ST));
		if (a.is_TPI() && b.is_LPI())
			return reverse_signs(lessThanOnAll_II<SSF>(b, a, PntArr3::LT));
		if (a.is_TPI() && b.is_TPI())
			return lessThanOnAll_II<SSF>(a, b, PntArr3::TT);
	}
	else
	{
		if (!a.is_Explicit() && b.is_Explicit())
			return lessThanOnAll_IE<DF>(a, b, PntArr3::I);
		if (a.is_Explicit() && !b.is_Explicit())
			return reverse_signs(lessThanOnAll_IE<DF>(b, a, PntArr3::I));
		return lessThanOnAll_II<DF>(a, b, PntArr3::II);
	}

	OMC_EXIT("LessThan3D - should not happen");
	// return std::array<Sign, 3>{Sign::ZERO, Sign::ZERO, Sign::ZERO};
	// warning killer
}

template <typename FT, typename IT, typename ET>
Sign LessThan3D_Indirect<FT, IT, ET>::on_x(const PointT &a, const FT *b)
{
	if (a.is_Explicit())
		return static_cast<Sign>(((a.x() > b.x()) - (a.x() < b.x())));
	if (a.has_ssf())
		return lessThanOnX_IE<SSF>(a, b[0], static_cast<PntArr3>(PntType(a)));
	else
		return lessThanOnX_IE<DF>(a, b[0], PntArr3::I);

	OMC_EXIT("LessThan3D - should not happen");
	return Sign::ZERO; // warning killer
}

template <typename FT, typename IT, typename ET>
Sign LessThan3D_Indirect<FT, IT, ET>::on_y(const PointT &a, const FT *b)
{
	if (a.is_Explicit())
		return static_cast<Sign>(((a.y() > b.y()) - (a.y() < b.y())));
	if (a.has_ssf())
		return lessThanOnY_IE<SSF>(a, b[1], static_cast<PntArr3>(PntType(a)));
	else
		return lessThanOnY_IE<DF>(a, b[1], PntArr3::I);

	OMC_EXIT("LessThan3D - should not happen");
	return Sign::ZERO; // warning killer
}

template <typename FT, typename IT, typename ET>
Sign LessThan3D_Indirect<FT, IT, ET>::on_z(const PointT &a, const FT *b)
{
	if (a.is_Explicit())
		return static_cast<Sign>(((a.z() > b.z()) - (a.z() < b.z())));
	if (a.has_ssf())
		return lessThanOnZ_IE<SSF>(a, b[2], static_cast<PntArr3>(PntType(a)));
	else
		return lessThanOnZ_IE<DF>(a, b[2], PntArr3::I);

	OMC_EXIT("LessThan3D - should not happen");
	return Sign::ZERO; // warning killer
}

template <typename FT, typename IT, typename ET>
std::array<Sign, 3> LessThan3D_Indirect<FT, IT, ET>::on_all(const PointT &a,
                                                            const FT     *b)
{
	if (a.is_Explicit())
		return lessThanOnAll_EE(a.x(), a.y(), a.z(), b[0], b[1], b[2]);
	if (a.has_ssf())
		return lessThanOnAll_IE<SSF>(a, b[0], b[1], b[2],
		                             static_cast<PntArr3>(PntType(a)));
	else
		return lessThanOnAll_IE<DF>(a, b[0], b[1], b[2]);

	OMC_EXIT("LessThan3D - should not happen");
	return std::array<Sign, 3>{Sign::ZERO, Sign::ZERO,
	                           Sign::ZERO}; // warning killer
}

template <typename FT, typename IT, typename ET>
Sign LessThan3D_Indirect<FT, IT, ET>::operator()(const PointT &a,
                                                 const PointT &b)
{
	if (a.is_Explicit() && b.is_Explicit())
		return lessThan_EE<IT, ET>(a, b);

	if (a.has_ssf() && b.has_ssf())
	{
		// cases are simple, expand them by hand.
		if (a.is_SSI() && b.is_Explicit())
			return lessThan_IE<SSF>(a, b, PntArr3::S);
		if (a.is_LPI() && b.is_Explicit())
			return lessThan_IE<SSF>(a, b, PntArr3::L);
		if (a.is_TPI() && b.is_Explicit())
			return lessThan_IE<SSF>(a, b, PntArr3::T);
		if (a.is_Explicit() && b.is_SSI())
			return reverse_sign(lessThan_IE<SSF>(b, a, PntArr3::S));
		if (a.is_Explicit() && b.is_LPI())
			return reverse_sign(lessThan_IE<SSF>(b, a, PntArr3::L));
		if (a.is_Explicit() && b.is_TPI())
			return reverse_sign(lessThan_IE<SSF>(b, a, PntArr3::T));
		if (a.is_SSI() && b.is_SSI())
			return lessThan_II<SSF>(a, b, PntArr3::SS);
		if (a.is_SSI() && b.is_LPI())
			return lessThan_II<SSF>(a, b, PntArr3::SL);
		if (a.is_SSI() && b.is_TPI())
			return lessThan_II<SSF>(a, b, PntArr3::ST);
		if (a.is_LPI() && b.is_SSI())
			return reverse_sign(lessThan_II<SSF>(b, a, PntArr3::SL));
		if (a.is_LPI() && b.is_LPI())
			return lessThan_II<SSF>(a, b, PntArr3::LL);
		if (a.is_LPI() && b.is_TPI())
			return lessThan_II<SSF>(a, b, PntArr3::LT);
		if (a.is_TPI() && b.is_SSI())
			return reverse_sign(lessThan_II<SSF>(b, a, PntArr3::ST));
		if (a.is_TPI() && b.is_LPI())
			return reverse_sign(lessThan_II<SSF>(b, a, PntArr3::LT));
		// if (a.is_TPI() && b.is_TPI())
		return lessThan_II<SSF>(a, b, PntArr3::TT);
	}
	else
	{
		// cases are simple, expand them by hand.
		if (!a.is_Explicit() && b.is_Explicit())
			return lessThan_IE<DF>(a, b, PntArr3::I);
		if (a.is_Explicit() && !b.is_Explicit())
			return reverse_sign(lessThan_IE<DF>(b, a, PntArr3::I));
		return lessThan_II<DF>(a, b, PntArr3::II);
	}
}

template <typename FT, typename IT, typename ET>
bool CollinearPoints3D_Indirect<FT, IT, ET>::misaligned(const PointT &A,
                                                        const PointT &B,
                                                        const PointT &C)
{
	return (is_sign_posneg(OrientOn2D().on_xy(A, B, C)) ||
	        is_sign_posneg(OrientOn2D().on_yz(A, B, C)) ||
	        is_sign_posneg(OrientOn2D().on_zx(A, B, C)));
}

template <typename FT, typename IT, typename ET>
bool CollinearPoints3D_Indirect<FT, IT, ET>::misaligned(const FT *A,
                                                        const FT *B,
                                                        const FT *C)
{
	return (is_sign_posneg(orient2dxy(A, B, C)) ||
	        is_sign_posneg(orient2dyz(A, B, C)) ||
	        is_sign_posneg(orient2dzx(A, B, C)));
}

template <typename FT, typename IT, typename ET>
bool CollinearPoints3D_Indirect<FT, IT, ET>::misaligned(const PointT &A,
                                                        const PointT &B,
                                                        const PointT &C,
                                                        int           n_max)
{
	return ((n_max == 2 && is_sign_posneg(OrientOn2D().on_xy(A, B, C))) ||
	        (n_max == 0 && is_sign_posneg(OrientOn2D().on_yz(A, B, C))) ||
	        (n_max == 1 && is_sign_posneg(OrientOn2D().on_zx(A, B, C))));
}

template <typename FT, typename IT, typename ET>
auto CollinearSort3D_Indirect<FT, IT, ET>::operator()(const PointT &p,
                                                      const PointT &q,
                                                      const PointT &r)
  -> std::tuple<const PointT &, const PointT &, const PointT &>
{
	using CR = const PointT &;
	using CP = const PointT *;

#define SORT_ON_AXIS(axis)                              \
	if (LessThan3D().on_##axis(*a, *b) == Sign::POSITIVE) \
		std::swap(a, b);                                    \
	if (LessThan3D().on_##axis(*b, *c) == Sign::POSITIVE) \
		std::swap(b, c);                                    \
	if (LessThan3D().on_##axis(*a, *c) == Sign::POSITIVE) \
		std::swap(a, c);

	CP a = &p, b = &q, c = &r;
	if (LessThan3D().on_x(p, q) != Sign::ZERO ||
	    LessThan3D().on_x(p, r) != Sign::ZERO ||
	    LessThan3D().on_x(q, r) != Sign::ZERO)
	{
		SORT_ON_AXIS(x);
	}
	else if (LessThan3D().on_y(p, q) != Sign::ZERO ||
	         LessThan3D().on_y(p, r) != Sign::ZERO ||
	         LessThan3D().on_y(q, r) != Sign::ZERO)
	{
		SORT_ON_AXIS(y);
	}
	else if (LessThan3D().on_z(p, q) != Sign::ZERO ||
	         LessThan3D().on_z(p, r) != Sign::ZERO ||
	         LessThan3D().on_z(q, r) != Sign::ZERO)
	{
		SORT_ON_AXIS(z);
	}
#undef SORT_ON_AXIS

	return std::make_tuple<CR, CR, CR>(*a, *b, *c);
}

template <typename FT, typename IT, typename ET>
int MaxComponentInTriangleNormal<FT, IT, ET>::operator()(FT ov1x, FT ov1y,
                                                         FT ov1z, FT ov2x,
                                                         FT ov2y, FT ov2z,
                                                         FT ov3x, FT ov3y,
                                                         FT ov3z)
{
	return maxComponentInTriangleNormal(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x,
	                                    ov3y, ov3z);
}

template <typename FT, typename IT, typename ET>
Sign InSphere_Indirect<FT, IT, ET>::operator()(const PointT &a, const PointT &b,
                                               const PointT &c, const PointT &d,
                                               const PointT &e)
{
	const int num_explicit = a.is_Explicit() + b.is_Explicit() + c.is_Explicit() +
	                         d.is_Explicit() + e.is_Explicit();
	if (num_explicit == 5)
		return inSphere<SSF>(a, b, c, d, e);

	// clang-format off
	std::array<uint32_t, 5> pos{0, 1, 2, 3, 4},
	  types{static_cast<uint32_t>(a.is_Explicit() ? PointT::PointType::Explicit : PointT::PointType::Implicit),
	        static_cast<uint32_t>(b.is_Explicit() ? PointT::PointType::Explicit : PointT::PointType::Implicit),
	        static_cast<uint32_t>(c.is_Explicit() ? PointT::PointType::Explicit : PointT::PointType::Implicit),
	        static_cast<uint32_t>(d.is_Explicit() ? PointT::PointType::Explicit : PointT::PointType::Implicit),
	        static_cast<uint32_t>(e.is_Explicit() ? PointT::PointType::Explicit : PointT::PointType::Implicit)};
	uint32_t swap_cnt;
	sort_pnts_arr3(types, pos, swap_cnt);
	const PointT *A[5] = {&a, &b, &c, &d, &e};
	Sign sign;
	if (num_explicit == 4) sign = inSphere_IEEEE<NF>(*A[pos[0]], *A[pos[1]], *A[pos[2]], *A[pos[3]], *A[pos[4]]);
	else if (num_explicit == 3) sign = inSphere_IIEEE<NF>(*A[pos[0]], *A[pos[1]], *A[pos[2]], *A[pos[3]], *A[pos[4]]);
	else if (num_explicit == 2) sign = inSphere_IIIEE<NF>(*A[pos[0]], *A[pos[1]], *A[pos[2]], *A[pos[3]], *A[pos[4]]);
	else if (num_explicit == 1) sign = inSphere_IIIIE<NF>(*A[pos[0]], *A[pos[1]], *A[pos[2]], *A[pos[3]], *A[pos[4]]);
	else sign = inSphere_IIIII<NF>(*A[pos[0]], *A[pos[1]], *A[pos[2]], *A[pos[3]], *A[pos[4]]);
	// clang-format on
	return swap_cnt % 2 == 1 ? reverse_sign(sign) : sign;
}

template <typename FT, typename IT, typename ET>
Sign InSphere_Indirect<FT, IT, ET>::operator()(const FT *a, const FT *b,
                                               const FT *c, const FT *d,
                                               const FT *e)
{
	// TODO optimize by shewchuk predicate
	return inSphere<IT, ET>(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2],
	                        d[0], d[1], d[2], e[0], e[1], e[2]);
}

#undef NF
#undef SSF
#undef DF
#undef PntType

} // namespace OMC