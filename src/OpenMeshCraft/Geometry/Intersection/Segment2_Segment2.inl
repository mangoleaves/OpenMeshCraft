#pragma once

#include "Segment2_Segment2.h"

namespace OMC {

template <typename Kernel>
bool Segment2_Segment2_Do_Intersect<Kernel>::operator()(
  const SegmentT &seg0, const SegmentT &seg1) const
{
	return intersection_type(seg0, seg1) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
SimplexIntersectionType
Segment2_Segment2_Do_Intersect<Kernel>::intersection_type(
  const SegmentT &seg0, const SegmentT &seg1) const
{
	return intersection_type(seg0.start(), seg0.end(), seg1.start(), seg1.end());
}

template <typename Kernel>
SimplexIntersectionType
Segment2_Segment2_Do_Intersect<Kernel>::intersection_type(
  const GPointT &s00, const GPointT &s01, const GPointT &s10,
  const GPointT &s11) const
{
	Sign s00_wrt_s1 = Orient2D()(s10, s11, s00);
	Sign s01_wrt_s1 = Orient2D()(s10, s11, s01);
	Sign s10_wrt_s0 = Orient2D()(s00, s01, s10);
	Sign s11_wrt_s0 = Orient2D()(s00, s01, s11);

	// segments intersect at a single point
	if (s00_wrt_s1 != s01_wrt_s1 && s10_wrt_s0 != s11_wrt_s0)
	{
		// edges share an endpoint
		if (LessThan2D().coincident(s00, s10) ||
		    LessThan2D().coincident(s00, s11) ||
		    LessThan2D().coincident(s01, s10) || LessThan2D().coincident(s01, s11))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		// at least one segment endpoint is involved in the intersection
		return SimplexIntersectionType::INTERSECT;
	}

	// degenerate case: collinear segments
	if (s00_wrt_s1 == Sign::ZERO && s01_wrt_s1 == Sign::ZERO &&
	    s10_wrt_s0 == Sign::ZERO && s11_wrt_s0 == Sign::ZERO)
	{
		// coincident segments
		if ((LessThan2D().coincident(s00, s10) &&
		     LessThan2D().coincident(s01, s11)) ||
		    (LessThan2D().coincident(s00, s11) &&
		     LessThan2D().coincident(s01, s10)))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		Sign Xmin_s0 = LessThan2D().on_x(s00, s01);
		Sign Ymin_s0 = LessThan2D().on_y(s00, s01);

		Sign Xmin_s1 = LessThan2D().on_x(s10, s11);
		Sign Ymin_s1 = LessThan2D().on_y(s10, s11);

		// clang-format off
		// test s00 against s1 range
		if (Xmin_s1 != Sign::ZERO &&
		    (Xmin_s1 == Sign::NEGATIVE && LessThan2D().on_x(s10, s00) == Sign::NEGATIVE && LessThan2D().on_x(s00, s11) == Sign::NEGATIVE ||
		     Xmin_s1 == Sign::POSITIVE && LessThan2D().on_x(s11, s00) == Sign::NEGATIVE && LessThan2D().on_x(s00, s10) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		if (Ymin_s1 != Sign::ZERO &&
		    (Ymin_s1 == Sign::NEGATIVE && LessThan2D().on_y(s10, s00) == Sign::NEGATIVE && LessThan2D().on_y(s00, s11) == Sign::NEGATIVE ||
		     Ymin_s1 == Sign::POSITIVE && LessThan2D().on_y(s11, s00) == Sign::NEGATIVE && LessThan2D().on_y(s00, s10) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		// test s01 against s1 range
		if (Xmin_s1 != Sign::ZERO &&
		    (Xmin_s1 == Sign::NEGATIVE && LessThan2D().on_x(s10, s01) == Sign::NEGATIVE && LessThan2D().on_x(s01, s11) == Sign::NEGATIVE ||
		     Xmin_s1 == Sign::POSITIVE && LessThan2D().on_x(s11, s01) == Sign::NEGATIVE && LessThan2D().on_x(s01, s10) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		if (Ymin_s1 != Sign::ZERO &&
		    (Ymin_s1 == Sign::NEGATIVE && LessThan2D().on_y(s10, s01) == Sign::NEGATIVE && LessThan2D().on_y(s01, s11) == Sign::NEGATIVE ||
		     Ymin_s1 == Sign::POSITIVE && LessThan2D().on_y(s11, s01) == Sign::NEGATIVE && LessThan2D().on_y(s01, s10) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		// test s10 against s0 range
		if (Xmin_s0 != Sign::ZERO &&
		    (Xmin_s0 == Sign::NEGATIVE && LessThan2D().on_x(s00, s10) == Sign::NEGATIVE && LessThan2D().on_x(s10, s01) == Sign::NEGATIVE ||
		     Xmin_s0 == Sign::POSITIVE && LessThan2D().on_x(s01, s10) == Sign::NEGATIVE && LessThan2D().on_x(s10, s00) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		if (Ymin_s0 != Sign::ZERO &&
		    (Ymin_s0 == Sign::NEGATIVE && LessThan2D().on_y(s00, s10) == Sign::NEGATIVE && LessThan2D().on_y(s10, s01) == Sign::NEGATIVE ||
		     Ymin_s0 == Sign::POSITIVE && LessThan2D().on_y(s01, s10) == Sign::NEGATIVE && LessThan2D().on_y(s10, s00) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		// test s11 against s0 range
		if (Xmin_s0 != Sign::ZERO &&
		    (Xmin_s0 == Sign::NEGATIVE && LessThan2D().on_x(s00, s11) == Sign::NEGATIVE && LessThan2D().on_x(s11, s01) == Sign::NEGATIVE ||
		     Xmin_s0 == Sign::POSITIVE && LessThan2D().on_x(s01, s11) == Sign::NEGATIVE && LessThan2D().on_x(s11, s00) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		if (Ymin_s0 != Sign::ZERO &&
		    (Ymin_s0 == Sign::NEGATIVE && LessThan2D().on_y(s00, s11) == Sign::NEGATIVE && LessThan2D().on_y(s11, s01) == Sign::NEGATIVE ||
		     Ymin_s0 == Sign::POSITIVE && LessThan2D().on_y(s01, s11) == Sign::NEGATIVE && LessThan2D().on_y(s11, s00) == Sign::NEGATIVE))
			return SimplexIntersectionType::OVERLAP;
		// clang-format on
	}
	return SimplexIntersectionType::DO_NOT_INTERSECT;
}

template <typename Kernel>
SimplexIntersectionType
Segment2_Segment2_Do_Intersect<Kernel>::intersection_type(const NT *s00,
                                                          const NT *s01,
                                                          const NT *s10,
                                                          const NT *s11) const
{
	Sign s00_wrt_s1 = Orient2D()(s10, s11, s00);
	Sign s01_wrt_s1 = Orient2D()(s10, s11, s01);
	Sign s10_wrt_s0 = Orient2D()(s00, s01, s10);
	Sign s11_wrt_s0 = Orient2D()(s00, s01, s11);

	// segments intersect at a single point
	if (s00_wrt_s1 != s01_wrt_s1 && s10_wrt_s0 != s11_wrt_s0)
	{
		// edges share an endpoint
		if (vec_equals_2d(s00, s10) || vec_equals_2d(s00, s11) ||
		    vec_equals_2d(s01, s10) || vec_equals_2d(s01, s11))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		// at least one segment endpoint is involved in the intersection
		return SimplexIntersectionType::INTERSECT;
	}

	// degenerate case: collinear segments
	if (s00_wrt_s1 == Sign::ZERO && s01_wrt_s1 == Sign::ZERO &&
	    s10_wrt_s0 == Sign::ZERO && s11_wrt_s0 == Sign::ZERO)
	{
		// coincident segments
		if ((vec_equals_2d(s00, s10) && vec_equals_2d(s01, s11)) ||
		    (vec_equals_2d(s00, s11) && vec_equals_2d(s01, s10)))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		NT Xmin_s1 = s10[0] < s11[0] ? s10[0] : s11[0];
		NT Xmax_s1 = s10[0] < s11[0] ? s11[0] : s10[0];
		NT Ymin_s1 = s10[1] < s11[1] ? s10[1] : s11[1];
		NT Ymax_s1 = s10[1] < s11[1] ? s11[1] : s10[1];
		NT Xmin_s0 = s00[0] < s01[0] ? s00[0] : s01[0];
		NT Xmax_s0 = s00[0] < s01[0] ? s01[0] : s00[0];
		NT Ymin_s0 = s00[1] < s01[1] ? s00[1] : s01[1];
		NT Ymax_s0 = s00[1] < s01[1] ? s01[1] : s00[1];

		if ( // test s0 endpoints against s1 range
		  (s00[0] > Xmin_s1 && s00[0] < Xmax_s1) ||
		  (s00[1] > Ymin_s1 && s00[1] < Ymax_s1) ||
		  (s01[0] > Xmin_s1 && s01[0] < Xmax_s1) ||
		  (s01[1] > Ymin_s1 && s01[1] < Ymax_s1) ||
		  // test s1 endpoints against s0 range
		  (s10[0] > Xmin_s0 && s10[0] < Xmax_s0) ||
		  (s10[1] > Ymin_s0 && s10[1] < Ymax_s0) ||
		  (s11[0] > Xmin_s0 && s11[0] < Xmax_s0) ||
		  (s11[1] > Ymin_s0 && s11[1] < Ymax_s0))
		{
			return SimplexIntersectionType::OVERLAP;
		}
	}
	return SimplexIntersectionType::DO_NOT_INTERSECT;
}

} // namespace OMC