#pragma once

#include "Segment3_Point3.h"

namespace OMC {

template <typename Kernel>
bool Segment3_Point3_Do_Intersect<Kernel>::operator()(const SegmentT &seg,
                                                      const GPointT  &pnt) const
{
	const GPointT &s0 = seg.start();
	const GPointT &s1 = seg.end();
	return operator()(s0, s1, pnt);
}

template <typename Kernel>
bool Segment3_Point3_Do_Intersect<Kernel>::operator()(const GPointT &s0,
                                                      const GPointT &s1,
                                                      const GPointT &pnt) const
{
	if (CollinearPoints3D().misaligned(s0, s1, pnt))
		return false; // strictly outside

	Sign seq_x = LessThan3D().on_x(s0, s1);
	// clang-format off
	if (seq_x != Sign::ZERO)
	{
	  if ((seq_x == Sign::NEGATIVE && LessThan3D().on_x(s0, pnt) <= Sign::ZERO && LessThan3D().on_x(pnt, s1) <= Sign::ZERO) ||
	      (seq_x == Sign::POSITIVE && LessThan3D().on_x(s1, pnt) <= Sign::ZERO && LessThan3D().on_x(pnt, s0) <= Sign::ZERO))
			return true;	// inside
		else
			return false;	// strictly outside
	}

	Sign seq_y = LessThan3D().on_y(s0, s1);
	if (seq_y != Sign::ZERO)
	{
	  if ((seq_y == Sign::NEGATIVE && LessThan3D().on_y(s0, pnt) <= Sign::ZERO && LessThan3D().on_y(pnt, s1) <= Sign::ZERO) ||
	      (seq_y == Sign::POSITIVE && LessThan3D().on_y(s1, pnt) <= Sign::ZERO && LessThan3D().on_y(pnt, s0) <= Sign::ZERO))
			return true;	// inside
		else
			return false;	// strictly outside
	}

	Sign seq_z = LessThan3D().on_z(s0, s1);
	if (seq_z != Sign::ZERO)
	{
	  if ((seq_z == Sign::NEGATIVE && LessThan3D().on_z(s0, pnt) <= Sign::ZERO && LessThan3D().on_z(pnt, s1) <= Sign::ZERO) ||
	      (seq_z == Sign::POSITIVE && LessThan3D().on_z(s1, pnt) <= Sign::ZERO && LessThan3D().on_z(pnt, s0) <= Sign::ZERO))
			return true;	// inside
		else
			return false;	// strictly outside
	}
	// clang-format on

	return false; // strictly outside;
}

template <typename Kernel>
PointInType
Segment3_Point3_Do_Intersect<Kernel>::in_segment(const SegmentT &seg,
                                                 const GPointT  &pnt) const
{
	return in_segment(seg.start(), seg.end(), pnt);
}

template <typename Kernel>
PointInType Segment3_Point3_Do_Intersect<Kernel>::in_segment(
  const GPointT &s0, const GPointT &s1, const GPointT &pnt) const
{
	PointInSimplexType type = intersection_type(s0, s1, pnt);
	if (type == PointInSimplexType::STRICTLY_INSIDE)
		return PointInType::STRICTLY_INSIDE;
	else if (type == PointInSimplexType::STRICTLY_OUTSIDE)
		return PointInType::STRICTLY_OUTSIDE;
	else
		return PointInType::ON_BOUNDARY;
}

template <typename Kernel>
PointInType Segment3_Point3_Do_Intersect<Kernel>::in_segment_collinear(
  const SegmentT &seg, const GPointT &pnt) const
{
	return in_segment_collinear(seg.start(), seg.end(), pnt);
}

template <typename Kernel>
PointInType Segment3_Point3_Do_Intersect<Kernel>::in_segment_collinear(
  const GPointT &s0, const GPointT &s1, const GPointT &pnt) const
{
	if (LessThan3D().coincident(pnt, s0))
		return PointInType::ON_BOUNDARY;
	if (LessThan3D().coincident(pnt, s1))
		return PointInType::ON_BOUNDARY;

	Sign seq_x = LessThan3D().on_x(s0, s1);
	// clang-format off
	if (seq_x != Sign::ZERO)
	{
	  if ((seq_x == Sign::NEGATIVE && LessThan3D().on_x(s0, pnt) == Sign::NEGATIVE && LessThan3D().on_x(pnt, s1) == Sign::NEGATIVE) ||
	      (seq_x == Sign::POSITIVE && LessThan3D().on_x(s1, pnt) == Sign::NEGATIVE && LessThan3D().on_x(pnt, s0) == Sign::NEGATIVE))
			return PointInType::STRICTLY_INSIDE;
		else
			return PointInType::STRICTLY_OUTSIDE;
	}

	Sign seq_y = LessThan3D().on_y(s0, s1);
	if (seq_y != Sign::ZERO)
	{
	  if ((seq_y == Sign::NEGATIVE && LessThan3D().on_y(s0, pnt) == Sign::NEGATIVE && LessThan3D().on_y(pnt, s1) == Sign::NEGATIVE) ||
	      (seq_y == Sign::POSITIVE && LessThan3D().on_y(s1, pnt) == Sign::NEGATIVE && LessThan3D().on_y(pnt, s0) == Sign::NEGATIVE))
			return PointInType::STRICTLY_INSIDE;
		else
			return PointInType::STRICTLY_OUTSIDE;
	}

	Sign seq_z = LessThan3D().on_z(s0, s1);
	if (seq_z != Sign::ZERO)
	{
	  if ((seq_z == Sign::NEGATIVE && LessThan3D().on_z(s0, pnt) == Sign::NEGATIVE && LessThan3D().on_z(pnt, s1) == Sign::NEGATIVE) ||
	      (seq_z == Sign::POSITIVE && LessThan3D().on_z(s1, pnt) == Sign::NEGATIVE && LessThan3D().on_z(pnt, s0) == Sign::NEGATIVE))
		  return PointInType::STRICTLY_INSIDE;
		else
			return PointInType::STRICTLY_OUTSIDE;
	}
	// clang-format on

	return PointInType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
PointInSimplexType Segment3_Point3_Do_Intersect<Kernel>::intersection_type(
  const SegmentT &seg, const GPointT &pnt) const
{
	return intersection_type(seg.start(), seg.end(), pnt);
}

template <typename Kernel>
PointInSimplexType Segment3_Point3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &s0, const GPointT &s1, const GPointT &p) const
{
	if (LessThan3D().coincident(p, s0))
		return PointInSimplexType::ON_VERT0;
	if (LessThan3D().coincident(p, s1))
		return PointInSimplexType::ON_VERT1;

	if (CollinearPoints3D().misaligned(s0, s1, p))
		return PointInSimplexType::STRICTLY_OUTSIDE;

	// clang-format off
	Sign seq_x = LessThan3D().on_x(s0, s1);

	if (seq_x != Sign::ZERO)
	{
	  if ((seq_x == Sign::NEGATIVE && LessThan3D().on_x(s0, p) == Sign::NEGATIVE && LessThan3D().on_x(p, s1) == Sign::NEGATIVE) ||
	      (seq_x == Sign::POSITIVE && LessThan3D().on_x(s1, p) == Sign::NEGATIVE && LessThan3D().on_x(p, s0) == Sign::NEGATIVE))
			return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}

	Sign seq_y = LessThan3D().on_y(s0, s1);
	if (seq_y != Sign::ZERO)
	{
	  if ((seq_y == Sign::NEGATIVE && LessThan3D().on_y(s0, p) == Sign::NEGATIVE && LessThan3D().on_y(p, s1) == Sign::NEGATIVE) ||
	      (seq_y == Sign::POSITIVE && LessThan3D().on_y(s1, p) == Sign::NEGATIVE && LessThan3D().on_y(p, s0) == Sign::NEGATIVE))
			return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}

	Sign seq_z = LessThan3D().on_z(s0, s1);
	if (seq_z != Sign::ZERO)
	{
	  if ((seq_z == Sign::NEGATIVE && LessThan3D().on_z(s0, p) == Sign::NEGATIVE && LessThan3D().on_z(p, s1) == Sign::NEGATIVE) ||
	      (seq_z == Sign::POSITIVE && LessThan3D().on_z(s1, p) == Sign::NEGATIVE && LessThan3D().on_z(p, s0) == Sign::NEGATIVE))
		  return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}
	// clang-format on

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
PointInSimplexType Segment3_Point3_Do_Intersect<Kernel>::intersection_type(
  const NT *s0, const NT *s1, const NT *p) const
{
	if (vec_equals_3d(p, s0))
		return PointInSimplexType::ON_VERT0;
	if (vec_equals_3d(p, s1))
		return PointInSimplexType::ON_VERT1;

	if (CollinearPoints3D().misaligned(s0, s1, p))
		return PointInSimplexType::STRICTLY_OUTSIDE;

	if ((p[0] > std::min(s0[0], s1[0]) && p[0] < std::max(s0[0], s1[0])) ||
	    (p[1] > std::min(s0[1], s1[1]) && p[1] < std::max(s0[1], s1[1])) ||
	    (p[2] > std::min(s0[2], s1[2]) && p[2] < std::max(s0[2], s1[2])))
	{
		return PointInSimplexType::STRICTLY_INSIDE;
	}

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
bool Segment3_Point3_Do_Intersect<Kernel>::operator()(const SegmentT &seg,
                                                      const GPointT  &pnt,
                                                      int n_max) const
{
	const GPointT &s0 = seg.start();
	const GPointT &s1 = seg.end();

	if (OrientOn2D().on_xy(s0, s1, pnt, n_max) != Sign::ZERO)
		return false; // strictly outside;

	Sign seq_x = LessThan3D().on_x(s0, s1);
	// clang-format off
	if (seq_x != Sign::ZERO)
	{
	  if ((seq_x == Sign::NEGATIVE && LessThan3D().on_x(s0, pnt) <= Sign::ZERO && LessThan3D().on_x(pnt, s1) <= Sign::ZERO) ||
	      (seq_x == Sign::POSITIVE && LessThan3D().on_x(s1, pnt) <= Sign::ZERO && LessThan3D().on_x(pnt, s0) <= Sign::ZERO))
			return true;	// inside
		else
			return false;	// strictly outside
	}

	Sign seq_y = LessThan3D().on_y(s0, s1);
	if (seq_y != Sign::ZERO)
	{
	  if ((seq_y == Sign::NEGATIVE && LessThan3D().on_y(s0, pnt) <= Sign::ZERO && LessThan3D().on_y(pnt, s1) <= Sign::ZERO) ||
	      (seq_y == Sign::POSITIVE && LessThan3D().on_y(s1, pnt) <= Sign::ZERO && LessThan3D().on_y(pnt, s0) <= Sign::ZERO))
			return true;	// inside
		else
			return false;	// strictly outside
	}

	Sign seq_z = LessThan3D().on_z(s0, s1);
	if (seq_z != Sign::ZERO)
	{
	  if ((seq_z == Sign::NEGATIVE && LessThan3D().on_z(s0, pnt) <= Sign::ZERO && LessThan3D().on_z(pnt, s1) <= Sign::ZERO) ||
	      (seq_z == Sign::POSITIVE && LessThan3D().on_z(s1, pnt) <= Sign::ZERO && LessThan3D().on_z(pnt, s0) <= Sign::ZERO))
			return true;	// inside
		else
			return false;	// strictly outside
	}
	// clang-format on

	return false; // strictly outside;
}

template <typename Kernel>
PointInType Segment3_Point3_Do_Intersect<Kernel>::in_segment(
  const SegmentT &seg, const GPointT &pnt, int n_max) const
{
	return in_segment(seg.start(), seg.end(), pnt, n_max);
}

template <typename Kernel>
PointInType Segment3_Point3_Do_Intersect<Kernel>::in_segment(const GPointT &s0,
                                                             const GPointT &s1,
                                                             const GPointT &pnt,
                                                             int n_max) const
{
	PointInSimplexType type = intersection_type(s0, s1, pnt, n_max);
	if (type == PointInSimplexType::STRICTLY_INSIDE)
		return PointInType::STRICTLY_INSIDE;
	else if (type == PointInSimplexType::STRICTLY_OUTSIDE)
		return PointInType::STRICTLY_OUTSIDE;
	else
		return PointInType::ON_BOUNDARY;
}

template <typename Kernel>
PointInSimplexType Segment3_Point3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &s0, const GPointT &s1, const GPointT &p, int n_max) const
{
	if (LessThan3D().coincident(p, s0))
		return PointInSimplexType::ON_VERT0;
	if (LessThan3D().coincident(p, s1))
		return PointInSimplexType::ON_VERT1;

	if (OrientOn2D()(s0, s1, p, n_max) != Sign::ZERO)
		return PointInSimplexType::STRICTLY_OUTSIDE;

	// clang-format off
	Sign seq_x = LessThan3D().on_x(s0, s1);
	if (seq_x != Sign::ZERO)
	{
	  if ((seq_x == Sign::NEGATIVE && LessThan3D().on_x(s0, p) == Sign::NEGATIVE && LessThan3D().on_x(p, s1) == Sign::NEGATIVE) ||
	      (seq_x == Sign::POSITIVE && LessThan3D().on_x(s1, p) == Sign::NEGATIVE && LessThan3D().on_x(p, s0) == Sign::NEGATIVE))
			return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}

	Sign seq_y = LessThan3D().on_y(s0, s1);
	if (seq_y != Sign::ZERO)
	{
	  if ((seq_y == Sign::NEGATIVE && LessThan3D().on_y(s0, p) == Sign::NEGATIVE && LessThan3D().on_y(p, s1) == Sign::NEGATIVE) ||
	      (seq_y == Sign::POSITIVE && LessThan3D().on_y(s1, p) == Sign::NEGATIVE && LessThan3D().on_y(p, s0) == Sign::NEGATIVE))
			return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}

	Sign seq_z = LessThan3D().on_z(s0, s1);
	if (seq_z != Sign::ZERO)
	{
	  if ((seq_z == Sign::NEGATIVE && LessThan3D().on_z(s0, p) == Sign::NEGATIVE && LessThan3D().on_z(p, s1) == Sign::NEGATIVE) ||
	      (seq_z == Sign::POSITIVE && LessThan3D().on_z(s1, p) == Sign::NEGATIVE && LessThan3D().on_z(p, s0) == Sign::NEGATIVE))
		  return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}
	// clang-format on

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
PointInSimplexType Segment3_Point3_Do_Intersect<Kernel>::intersection_type(
  const NT *s0, const NT *s1, const NT *p, int n_max) const
{
	if (vec_equals_3d(p, s0))
		return PointInSimplexType::ON_VERT0;
	if (vec_equals_3d(p, s1))
		return PointInSimplexType::ON_VERT1;

	if (OrientOn2D()(s0, s1, p, n_max) != Sign::ZERO)
		return PointInSimplexType::STRICTLY_OUTSIDE;

	if ((p[0] > std::min(s0[0], s1[0]) && p[0] < std::max(s0[0], s1[0])) ||
	    (p[1] > std::min(s0[1], s1[1]) && p[1] < std::max(s0[1], s1[1])) ||
	    (p[2] > std::min(s0[2], s1[2]) && p[2] < std::max(s0[2], s1[2])))
	{
		return PointInSimplexType::STRICTLY_INSIDE;
	}

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

} // namespace OMC