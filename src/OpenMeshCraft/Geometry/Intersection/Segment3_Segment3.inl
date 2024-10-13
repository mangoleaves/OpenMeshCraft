#pragma once

#include "Segment3_Segment3.h"

namespace OMC {

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::operator()(
  const SegmentT &seg0, const SegmentT &seg1) const
{
	if (Orient3D()(seg0.start(), seg0.end(), seg1.start(), seg1.end()) !=
	    Sign::ZERO)
		return false;
	return intersection_type(seg0, seg1) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::operator()(const GPointT &A,
                                                        const GPointT &B,
                                                        const GPointT &P,
                                                        const GPointT &Q) const
{
	if (Orient3D()(A, B, P, Q) != Sign::ZERO)
		return false;
	return intersection_type(A, B, P, Q) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::cross(const SegmentT &seg0,
                                                   const SegmentT &seg1) const
{
	return cross(seg0.start(), seg0.end(), seg1.start(), seg1.end());
}

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::cross_inner(
  const SegmentT &seg0, const SegmentT &seg1) const
{
	return cross_inner(seg0.start(), seg0.end(), seg1.start(), seg1.end());
}

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::cross(const GPointT &A,
                                                   const GPointT &B,
                                                   const GPointT &P,
                                                   const GPointT &Q) const
{
	Sign o11, o12, o21, o22;

	o11 = OrientOn2D().on_xy(P, A, B);
	o12 = OrientOn2D().on_xy(Q, B, A);
	o21 = OrientOn2D().on_xy(A, P, Q);
	o22 = OrientOn2D().on_xy(B, Q, P);
	if ((o11 != Sign::ZERO || o12 != Sign::ZERO) &&
	    (static_cast<int>(o11) * static_cast<int>(o12) >= 0) &&
	    (o21 != Sign::ZERO || o22 != Sign::ZERO) &&
	    (static_cast<int>(o21) * static_cast<int>(o22) >= 0))
		return true;

	o11 = OrientOn2D().on_yz(P, A, B);
	o12 = OrientOn2D().on_yz(Q, B, A);
	o21 = OrientOn2D().on_yz(A, P, Q);
	o22 = OrientOn2D().on_yz(B, Q, P);
	if ((o11 != Sign::ZERO || o12 != Sign::ZERO) &&
	    (static_cast<int>(o11) * static_cast<int>(o12) >= 0) &&
	    (o21 != Sign::ZERO || o22 != Sign::ZERO) &&
	    (static_cast<int>(o21) * static_cast<int>(o22) >= 0))
		return true;

	o11 = OrientOn2D().on_zx(P, A, B);
	o12 = OrientOn2D().on_zx(Q, B, A);
	o21 = OrientOn2D().on_zx(A, P, Q);
	o22 = OrientOn2D().on_zx(B, Q, P);
	if ((o11 != Sign::ZERO || o12 != Sign::ZERO) &&
	    (static_cast<int>(o11) * static_cast<int>(o12) >= 0) &&
	    (o21 != Sign::ZERO || o22 != Sign::ZERO) &&
	    (static_cast<int>(o21) * static_cast<int>(o22) >= 0))
		return true;

	return false;
}

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::cross_inner(const GPointT &A,
                                                         const GPointT &B,
                                                         const GPointT &P,
                                                         const GPointT &Q) const
{
	Sign o11, o12, o21, o22;

	o11 = OrientOn2D().on_xy(P, A, B);
	o12 = OrientOn2D().on_xy(Q, B, A);
	o21 = OrientOn2D().on_xy(A, P, Q);
	o22 = OrientOn2D().on_xy(B, Q, P);
	if (o11 != Sign::ZERO || o21 != Sign::ZERO || o12 != Sign::ZERO ||
	    o22 != Sign::ZERO)
		return (o11 == o12 && o21 == o22);

	o11 = OrientOn2D().on_yz(P, A, B);
	o12 = OrientOn2D().on_yz(Q, B, A);
	o21 = OrientOn2D().on_yz(A, P, Q);
	o22 = OrientOn2D().on_yz(B, Q, P);
	if (o11 != Sign::ZERO || o21 != Sign::ZERO || o12 != Sign::ZERO ||
	    o22 != Sign::ZERO)
		return (o11 == o12 && o21 == o22);

	o11 = OrientOn2D().on_zx(P, A, B);
	o12 = OrientOn2D().on_zx(Q, B, A);
	o21 = OrientOn2D().on_zx(A, P, Q);
	o22 = OrientOn2D().on_zx(B, Q, P);
	if (o11 != Sign::ZERO || o21 != Sign::ZERO || o12 != Sign::ZERO ||
	    o22 != Sign::ZERO)
		return (o11 == o12 && o21 == o22);

	return false;
}

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::cross(const GPointT &A,
                                                   const GPointT &B,
                                                   const GPointT &P,
                                                   const GPointT &Q,
                                                   int            n_max) const
{
	Sign o11, o12, o21, o22;

	if (n_max == 2)
	{
		o11 = OrientOn2D().on_xy(P, A, B);
		o12 = OrientOn2D().on_xy(Q, B, A);
		o21 = OrientOn2D().on_xy(A, P, Q);
		o22 = OrientOn2D().on_xy(B, Q, P);
	}
	else if (n_max == 0)
	{
		o11 = OrientOn2D().on_yz(P, A, B);
		o12 = OrientOn2D().on_yz(Q, B, A);
		o21 = OrientOn2D().on_yz(A, P, Q);
		o22 = OrientOn2D().on_yz(B, Q, P);
	}
	else
	{
		o11 = OrientOn2D().on_zx(P, A, B);
		o12 = OrientOn2D().on_zx(Q, B, A);
		o21 = OrientOn2D().on_zx(A, P, Q);
		o22 = OrientOn2D().on_zx(B, Q, P);
	}

	return ((o11 != Sign::ZERO || o12 != Sign::ZERO) &&
	        (static_cast<int>(o11) * static_cast<int>(o12) >= 0) &&
	        (o21 != Sign::ZERO || o22 != Sign::ZERO) &&
	        (static_cast<int>(o21) * static_cast<int>(o22) >= 0));
}

template <typename Kernel>
bool Segment3_Segment3_Do_Intersect<Kernel>::cross_inner(const GPointT &A,
                                                         const GPointT &B,
                                                         const GPointT &P,
                                                         const GPointT &Q,
                                                         int n_max) const
{
	Sign o11, o12, o21, o22;

	if (n_max == 2)
	{
		o11 = OrientOn2D().on_xy(P, A, B);
		o12 = OrientOn2D().on_xy(Q, B, A);
		o21 = OrientOn2D().on_xy(A, P, Q);
		o22 = OrientOn2D().on_xy(B, Q, P);
	}
	else if (n_max == 0)
	{
		o11 = OrientOn2D().on_yz(P, A, B);
		o12 = OrientOn2D().on_yz(Q, B, A);
		o21 = OrientOn2D().on_yz(A, P, Q);
		o22 = OrientOn2D().on_yz(B, Q, P);
	}
	else
	{
		o11 = OrientOn2D().on_zx(P, A, B);
		o12 = OrientOn2D().on_zx(Q, B, A);
		o21 = OrientOn2D().on_zx(A, P, Q);
		o22 = OrientOn2D().on_zx(B, Q, P);
	}

	return ((o11 != Sign::ZERO || o21 != Sign::ZERO || o12 != Sign::ZERO ||
	         o22 != Sign::ZERO) &&
	        o11 == o12 && o21 == o22);
}

template <typename Kernel>
SimplexIntersectionType
Segment3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const SegmentT &seg0, const SegmentT &seg1) const
{
	return intersection_type(seg0.start(), seg0.end(), seg1.start(), seg0.end());
}

template <typename Kernel>
SimplexIntersectionType
Segment3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &s00, const GPointT &s01, const GPointT &s10,
  const GPointT &s11) const
{
	bool s00_s10_same = LessThan3D().coincident(s00, s10);
	bool s00_s11_same = LessThan3D().coincident(s00, s11);
	bool s01_s10_same = LessThan3D().coincident(s01, s10);
	bool s01_s11_same = LessThan3D().coincident(s01, s11);

	// check for coincident segments
	bool s00_is_shared = (s00_s10_same || s00_s11_same);
	bool s01_is_shared = (s01_s10_same || s01_s11_same);
	bool s10_is_shared = (s00_s10_same || s01_s10_same);
	bool s11_is_shared = (s00_s11_same || s01_s11_same);

	// s0 and s1 are coincident or one edge is degenerate and coincides with one
	// vertex of the other
	// (sufficient but not neccessary condition)
	if (s00_is_shared && s01_is_shared && s10_is_shared && s11_is_shared)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	// check 2D projections of the segments
	// OPT: find a projection plane where both segments do not degenerate,
	// check 2D projections of the segments only once.

	// Do not intersect on one projection plane <=> Do not intersect in 3D
	// Intersect on one projection plane <=> Intersect in 3D
	// (sufficient and necessary condition)
	SimplexIntersectionType x_res = intersection_type_on2d(s00, s01, s10, s11, 0);
	if (x_res == SimplexIntersectionType::DO_NOT_INTERSECT ||
	    x_res == SimplexIntersectionType::INTERSECT)
		return x_res;

	SimplexIntersectionType y_res = intersection_type_on2d(s00, s01, s10, s11, 1);
	if (y_res == SimplexIntersectionType::DO_NOT_INTERSECT ||
	    y_res == SimplexIntersectionType::INTERSECT)
		return y_res;

	SimplexIntersectionType z_res = intersection_type_on2d(s00, s01, s10, s11, 2);
	if (z_res == SimplexIntersectionType::DO_NOT_INTERSECT ||
	    z_res == SimplexIntersectionType::INTERSECT)
		return z_res;

	// Overlap in three projection planes <=> Overlap in 3D
	// (sufficient and necessary condition)
	if (x_res == SimplexIntersectionType::OVERLAP &&
	    y_res == SimplexIntersectionType::OVERLAP &&
	    z_res == SimplexIntersectionType::OVERLAP)
		return SimplexIntersectionType::OVERLAP;

	if (x_res == SimplexIntersectionType::SIMPLICIAL_COMPLEX ||
	    y_res == SimplexIntersectionType::SIMPLICIAL_COMPLEX ||
	    z_res == SimplexIntersectionType::SIMPLICIAL_COMPLEX)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	OMC_ASSERT(false, "Impossible case.");
	return SimplexIntersectionType::INTERSECT; // just kill compiler warning
}

template <typename Kernel>
SimplexIntersectionType
Segment3_Segment3_Do_Intersect<Kernel>::intersection_type(const NT *s00,
                                                          const NT *s01,
                                                          const NT *s10,
                                                          const NT *s11) const
{
	bool s00_s10_same = vec_equals_3d(s00, s10);
	bool s00_s11_same = vec_equals_3d(s00, s11);
	bool s01_s10_same = vec_equals_3d(s01, s10);
	bool s01_s11_same = vec_equals_3d(s01, s11);

	// check for coincident segments
	bool s00_is_shared = (s00_s10_same || s00_s11_same);
	bool s01_is_shared = (s01_s10_same || s01_s11_same);
	bool s10_is_shared = (s00_s10_same || s01_s10_same);
	bool s11_is_shared = (s00_s11_same || s01_s11_same);

	// s0 and s1 are coincident or one edge is degenerate and coincides with one
	// vertex of the other
	// (sufficient but not neccessary condition)
	if (s00_is_shared && s01_is_shared && s10_is_shared && s11_is_shared)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	// check 2D projections of the segments
	// OPT: find a projection plane where both segments do not degenerate,
	// check 2D projections of the segments only once.

	// Do not intersect on one projection plane <=> Do not intersect in 3D
	// Intersect on one projection plane <=> Intersect in 3D
	// (sufficient and necessary condition)
	NT s00_dropX[2] = {s00[1], s00[2]}, s01_dropX[2] = {s01[1], s01[2]};
	NT s10_dropX[2] = {s10[1], s10[2]}, s11_dropX[2] = {s11[1], s11[2]};
	SimplexIntersectionType x_res = Segment2_Segment2_DoInter().intersection_type(
	  s00_dropX, s01_dropX, s10_dropX, s11_dropX);
	if (x_res == SimplexIntersectionType::DO_NOT_INTERSECT ||
	    x_res == SimplexIntersectionType::INTERSECT)
		return x_res;

	NT s00_dropY[2] = {s00[0], s00[2]}, s01_dropY[2] = {s01[0], s01[2]};
	NT s10_dropY[2] = {s10[0], s10[2]}, s11_dropY[2] = {s11[0], s11[2]};
	SimplexIntersectionType y_res = Segment2_Segment2_DoInter().intersection_type(
	  s00_dropY, s01_dropY, s10_dropY, s11_dropY);
	if (y_res == SimplexIntersectionType::DO_NOT_INTERSECT ||
	    y_res == SimplexIntersectionType::INTERSECT)
		return y_res;

	NT s00_dropZ[2] = {s00[0], s00[1]}, s01_dropZ[2] = {s01[0], s01[1]};
	NT s10_dropZ[2] = {s10[0], s10[1]}, s11_dropZ[2] = {s11[0], s11[1]};
	SimplexIntersectionType z_res = Segment2_Segment2_DoInter().intersection_type(
	  s00_dropZ, s01_dropZ, s10_dropZ, s11_dropZ);
	if (z_res == SimplexIntersectionType::DO_NOT_INTERSECT ||
	    z_res == SimplexIntersectionType::INTERSECT)
		return z_res;

	// Overlap in three projection planes <=> Overlap in 3D
	// (sufficient and necessary condition)
	if (x_res == SimplexIntersectionType::OVERLAP &&
	    y_res == SimplexIntersectionType::OVERLAP &&
	    z_res == SimplexIntersectionType::OVERLAP)
		return SimplexIntersectionType::OVERLAP;

	if (x_res == SimplexIntersectionType::SIMPLICIAL_COMPLEX ||
	    y_res == SimplexIntersectionType::SIMPLICIAL_COMPLEX ||
	    z_res == SimplexIntersectionType::SIMPLICIAL_COMPLEX)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	OMC_ASSERT(false, "Impossible case.");
	return SimplexIntersectionType::INTERSECT; // just kill compiler warning
}

template <typename Kernel>
SimplexIntersectionType
Segment3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const NT *s00, const NT *s01, const NT *s10, const NT *s11, int n1_max) const
{
	bool s00_s10_same = vec_equals_3d(s00, s10);
	bool s00_s11_same = vec_equals_3d(s00, s11);
	bool s01_s10_same = vec_equals_3d(s01, s10);
	bool s01_s11_same = vec_equals_3d(s01, s11);

	// check for coincident segments
	bool s00_is_shared = (s00_s10_same || s00_s11_same);
	bool s01_is_shared = (s01_s10_same || s01_s11_same);
	bool s10_is_shared = (s00_s10_same || s01_s10_same);
	bool s11_is_shared = (s00_s11_same || s01_s11_same);

	// s0 and s1 are coincident or one edge is degenerate and coincides with one
	// vertex of the other
	// (sufficient but not neccessary condition)
	if (s00_is_shared && s01_is_shared && s10_is_shared && s11_is_shared)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	// check 2D projections of the segments

	if (n1_max == 0)
	{
		NT s00_dropX[2] = {s00[1], s00[2]}, s01_dropX[2] = {s01[1], s01[2]};
		NT s10_dropX[2] = {s10[1], s10[2]}, s11_dropX[2] = {s11[1], s11[2]};
		return Segment2_Segment2_DoInter().intersection_type(s00_dropX, s01_dropX,
		                                                     s10_dropX, s11_dropX);
	}
	else if (n1_max == 1)
	{
		NT s00_dropY[2] = {s00[0], s00[2]}, s01_dropY[2] = {s01[0], s01[2]};
		NT s10_dropY[2] = {s10[0], s10[2]}, s11_dropY[2] = {s11[0], s11[2]};
		return Segment2_Segment2_DoInter().intersection_type(s00_dropY, s01_dropY,
		                                                     s10_dropY, s11_dropY);
	}
	else if (n1_max == 2)
	{
		NT s00_dropZ[2] = {s00[0], s00[1]}, s01_dropZ[2] = {s01[0], s01[1]};
		NT s10_dropZ[2] = {s10[0], s10[1]}, s11_dropZ[2] = {s11[0], s11[1]};
		return Segment2_Segment2_DoInter().intersection_type(s00_dropZ, s01_dropZ,
		                                                     s10_dropZ, s11_dropZ);
	}

	OMC_ASSERT(false, "Impossible case.");
	return SimplexIntersectionType::INTERSECT; // just kill compiler warning
}

template <typename Kernel>
SimplexIntersectionType
Segment3_Segment3_Do_Intersect<Kernel>::intersection_type_on2d(
  const GPointT &s00, const GPointT &s01, const GPointT &s10,
  const GPointT &s11, int on2d) const
{
#define Same_On(p, q, dim) (LessThan3D().on_##dim(p, q) == Sign::ZERO)
#define Less_On(p, q, dim) (LessThan3D().on_##dim(p, q) == Sign::NEGATIVE)

	// Almostly same as Segment2_Segment2_Do_Intersect.
	// So I just copy it here and modify it to macro to test intersection on
	// different 2d planes.
#define On_Plane(xy, x, y)                                                     \
	Sign s00_wrt_s1 = OrientOn2D().on_##xy(s10, s11, s00);                       \
	Sign s01_wrt_s1 = OrientOn2D().on_##xy(s10, s11, s01);                       \
	Sign s10_wrt_s0 = OrientOn2D().on_##xy(s00, s01, s10);                       \
	Sign s11_wrt_s0 = OrientOn2D().on_##xy(s00, s01, s11);                       \
	if ((s00_wrt_s1 == s01_wrt_s1 && s00_wrt_s1 != Sign::ZERO) ||                \
	    (s10_wrt_s0 == s11_wrt_s0 && s10_wrt_s0 != Sign::ZERO))                  \
	{                                                                            \
		return SimplexIntersectionType::DO_NOT_INTERSECT;                          \
	}                                                                            \
	if (s00_wrt_s1 != s01_wrt_s1 && s10_wrt_s0 != s11_wrt_s0)                    \
	{                                                                            \
		if ((Same_On(s00, s10, x) && Same_On(s00, s10, y)) ||                      \
		    (Same_On(s00, s11, x) && Same_On(s00, s11, y)) ||                      \
		    (Same_On(s01, s10, x) && Same_On(s01, s10, y)) ||                      \
		    (Same_On(s01, s11, x) && Same_On(s01, s11, y)))                        \
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;                      \
                                                                               \
		return SimplexIntersectionType::INTERSECT;                                 \
	}                                                                            \
	if ((Same_On(s00, s10, x) && Same_On(s00, s10, y) && Same_On(s01, s11, x) && \
	     Same_On(s01, s11, y)) ||                                                \
	    (Same_On(s00, s11, x) && Same_On(s00, s11, y) && Same_On(s01, s10, x) && \
	     Same_On(s01, s10, y)))                                                  \
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;                        \
	Sign Xmin_s0 = LessThan3D().on_##x(s00, s01);                                \
	Sign Ymin_s0 = LessThan3D().on_##y(s00, s01);                                \
	Sign Xmin_s1 = LessThan3D().on_##x(s10, s11);                                \
	Sign Ymin_s1 = LessThan3D().on_##y(s10, s11);                                \
	if (Xmin_s1 != Sign::ZERO &&                                                 \
	    ((Xmin_s1 == Sign::NEGATIVE && Less_On(s10, s00, x) &&                   \
	      Less_On(s00, s11, x)) ||                                               \
	     (Xmin_s1 == Sign::POSITIVE && Less_On(s11, s00, x) &&                   \
	      Less_On(s00, s10, x))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	if (Ymin_s1 != Sign::ZERO &&                                                 \
	    ((Ymin_s1 == Sign::NEGATIVE && Less_On(s10, s00, y) &&                   \
	      Less_On(s00, s11, y)) ||                                               \
	     (Ymin_s1 == Sign::POSITIVE && Less_On(s11, s00, y) &&                   \
	      Less_On(s00, s10, y))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	if (Xmin_s1 != Sign::ZERO &&                                                 \
	    ((Xmin_s1 == Sign::NEGATIVE && Less_On(s10, s01, x) &&                   \
	      Less_On(s01, s11, x)) ||                                               \
	     (Xmin_s1 == Sign::POSITIVE && Less_On(s11, s01, x) &&                   \
	      Less_On(s01, s10, x))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	if (Ymin_s1 != Sign::ZERO &&                                                 \
	    ((Ymin_s1 == Sign::NEGATIVE && Less_On(s10, s01, y) &&                   \
	      Less_On(s01, s11, y)) ||                                               \
	     (Ymin_s1 == Sign::POSITIVE && Less_On(s11, s01, y) &&                   \
	      Less_On(s01, s10, y))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	if (Xmin_s0 != Sign::ZERO &&                                                 \
	    ((Xmin_s0 == Sign::NEGATIVE && Less_On(s00, s10, x) &&                   \
	      Less_On(s10, s01, x)) ||                                               \
	     (Xmin_s0 == Sign::POSITIVE && Less_On(s01, s10, x) &&                   \
	      Less_On(s10, s00, x))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	if (Ymin_s0 != Sign::ZERO &&                                                 \
	    ((Ymin_s0 == Sign::NEGATIVE && Less_On(s00, s10, y) &&                   \
	      Less_On(s10, s01, y)) ||                                               \
	     (Ymin_s0 == Sign::POSITIVE && Less_On(s01, s10, y) &&                   \
	      Less_On(s10, s00, y))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	if (Xmin_s0 != Sign::ZERO &&                                                 \
	    ((Xmin_s0 == Sign::NEGATIVE && Less_On(s00, s11, x) &&                   \
	      Less_On(s11, s01, x)) ||                                               \
	     (Xmin_s0 == Sign::POSITIVE && Less_On(s01, s11, x) &&                   \
	      Less_On(s11, s00, x))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	if (Ymin_s0 != Sign::ZERO &&                                                 \
	    ((Ymin_s0 == Sign::NEGATIVE && Less_On(s00, s11, y) &&                   \
	      Less_On(s11, s01, y)) ||                                               \
	     (Ymin_s0 == Sign::POSITIVE && Less_On(s01, s11, y) &&                   \
	      Less_On(s11, s00, y))))                                                \
		return SimplexIntersectionType::OVERLAP;                                   \
	return SimplexIntersectionType::DO_NOT_INTERSECT;

	if (on2d == 0)
	{
		On_Plane(yz, y, z);
	}
	else if (on2d == 1)
	{
		On_Plane(zx, z, x);
	}
	else // if (on2d == 2)
	{
		On_Plane(xy, x, y);
	}
}

} // namespace OMC