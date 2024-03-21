#pragma once

#include "Triangle3_Segment3.h"

namespace OMC {

template <typename Kernel>
bool Triangle3_Segment3_Do_Intersect<Kernel>::operator()(
  const TriangleT &tri, const SegmentT &seg) const
{
	return intersection_type(tri, seg) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
bool Triangle3_Segment3_Do_Intersect<Kernel>::operator()(
  const GPointT &v1, const GPointT &v2, const GPointT &v3, const GPointT &s1,
  const GPointT &s2) const
{
	return intersection_type(v1, v2, v3, s1, s2) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
bool Triangle3_Segment3_Do_Intersect<Kernel>::cross(const TriangleT &tri,
                                                    const SegmentT  &seg) const
{
	return cross(tri.v0(), tri.v1(), tri.v2(), seg.start(), seg.end());
}

template <typename Kernel>
bool Triangle3_Segment3_Do_Intersect<Kernel>::cross_inner(
  const TriangleT &tri, const SegmentT &seg) const
{
	return cross_inner(tri.v0(), tri.v1(), tri.v2(), seg.start(), seg.end());
}

template <typename Kernel>
bool Triangle3_Segment3_Do_Intersect<Kernel>::cross(const GPointT &v1,
                                                    const GPointT &v2,
                                                    const GPointT &v3,
                                                    const GPointT &s1,
                                                    const GPointT &s2) const
{
	Sign o1 = Orient3D()(v1, v2, v3, s1);
	Sign o2 = Orient3D()(v1, v2, v3, s2);

	if (o1 == Sign::ZERO && o2 == Sign::ZERO) // coplanar
		return false;
	if ((o1 > Sign::ZERO && o2 > Sign::ZERO) ||
	    (o1 < Sign::ZERO && o2 < Sign::ZERO)) // above or below
		return false;

	// s intersects t (borders included), if the signs of the three tetrahedra
	// obtained combining s with the three edges of t are all equal
	o1 = Orient3D()(s1, s2, v1, v2);
	o2 = Orient3D()(s1, s2, v2, v3);
	if ((o1 > Sign::ZERO && o2 < Sign::ZERO) ||
	    (o1 < Sign::ZERO && o2 > Sign::ZERO))
		return false;
	Sign o3 = Orient3D()(s1, s2, v3, v1);
	if ((o1 > Sign::ZERO && o3 < Sign::ZERO) ||
	    (o1 < Sign::ZERO && o3 > Sign::ZERO))
		return false;
	if ((o2 > Sign::ZERO && o3 < Sign::ZERO) ||
	    (o2 < Sign::ZERO && o3 > Sign::ZERO))
		return false;

	return true;
}

template <typename Kernel>
bool Triangle3_Segment3_Do_Intersect<Kernel>::cross_inner(
  const GPointT &v1, const GPointT &v2, const GPointT &v3, const GPointT &s1,
  const GPointT &s2) const
{
	Sign o1 = Orient3D()(v1, v2, v3, s1);
	if (o1 == Sign::ZERO)
		return false;
	Sign o2 = Orient3D()(v1, v2, v3, s2);
	if (o2 == Sign::ZERO)
		return false;

	if ((o1 > Sign::ZERO && o2 > Sign::ZERO) ||
	    (o1 < Sign::ZERO && o2 < Sign::ZERO))
		return false;
	o1 = Orient3D()(s1, s2, v1, v2);
	o2 = Orient3D()(s1, s2, v2, v3);
	if ((o1 >= Sign::ZERO && o2 <= Sign::ZERO) ||
	    (o1 <= Sign::ZERO && o2 >= Sign::ZERO))
		return false;
	Sign o3 = Orient3D()(s1, s2, v3, v1);
	if ((o1 >= Sign::ZERO && o3 <= Sign::ZERO) ||
	    (o1 <= Sign::ZERO && o3 >= Sign::ZERO))
		return false;
	if ((o2 >= Sign::ZERO && o3 <= Sign::ZERO) ||
	    (o2 <= Sign::ZERO && o3 >= Sign::ZERO))
		return false;
	return true;
}

template <typename Kernel>
SimplexIntersectionType
Triangle3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const TriangleT &tri, const SegmentT &seg) const
{
	return intersection_type(tri.v0(), tri.v1(), tri.v2(), seg.start(),
	                         seg.end());
}

template <typename Kernel>
SimplexIntersectionType
Triangle3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &t0, const GPointT &t1, const GPointT &t2, const GPointT &s0,
  const GPointT &s1) const
{
	bool s0_in_vertices =
	  (LessThan3D().coincident(s0, t0) || LessThan3D().coincident(s0, t1) ||
	   LessThan3D().coincident(s0, t2));
	bool s1_in_vertices =
	  (LessThan3D().coincident(s1, t0) || LessThan3D().coincident(s1, t1) ||
	   LessThan3D().coincident(s1, t2));
	if (s0_in_vertices && s1_in_vertices)
	{
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;
	}

	Sign vol_s0_t, vol_s1_t;
	vol_s0_t = Orient3D()(t0, t1, t2, s0);
	vol_s1_t = Orient3D()(t0, t1, t2, s1);

	if (vol_s0_t > Sign::ZERO && vol_s1_t > Sign::ZERO) // s is above t
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	if (vol_s0_t < Sign::ZERO && vol_s1_t < Sign::ZERO) // s is below t
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	if (vol_s0_t == Sign::ZERO && vol_s1_t == Sign::ZERO) // s and t are coplanar
	{
		// same code as the 2D version, I just copied it here....
		if (Triangle3_Point3_DoInter().intersection_type(t0, t1, t2, s0) !=
		      PointInSimplexType::STRICTLY_OUTSIDE ||
		    Triangle3_Point3_DoInter().intersection_type(t0, t1, t2, s1) !=
		      PointInSimplexType::STRICTLY_OUTSIDE)
			return SimplexIntersectionType::INTERSECT;

		switch (Segment3_Segment3_DoInter().intersection_type(s0, s1, t0, t1))
		{
		case SimplexIntersectionType::INTERSECT:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::OVERLAP:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::SIMPLICIAL_COMPLEX:
			break;
		case SimplexIntersectionType::DO_NOT_INTERSECT:
			break;
		}

		switch (Segment3_Segment3_DoInter().intersection_type(s0, s1, t1, t2))
		{
		case SimplexIntersectionType::INTERSECT:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::OVERLAP:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::SIMPLICIAL_COMPLEX:
			break;
		case SimplexIntersectionType::DO_NOT_INTERSECT:
			break;
		}

		switch (Segment3_Segment3_DoInter().intersection_type(s0, s1, t2, t0))
		{
		case SimplexIntersectionType::INTERSECT:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::OVERLAP:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::SIMPLICIAL_COMPLEX:
			break;
		case SimplexIntersectionType::DO_NOT_INTERSECT:
			break;
		}

		// if it is a simplicial complex from any view, then it really is...
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	}

	// s intersects t (borders included), if the signs of the three tetrahedra
	// obtained combining s with the three edges of t are all equal

	// if one point is coplanar and coincides with a triangle vertex, then they
	// form a valid complex
	if (s0_in_vertices || s1_in_vertices)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	Sign vol_s_t01, vol_s_t12;
	vol_s_t01 = Orient3D()(s0, s1, t0, t1);
	vol_s_t12 = Orient3D()(s0, s1, t1, t2);

	if ((vol_s_t01 > Sign::ZERO && vol_s_t12 < Sign::ZERO) ||
	    (vol_s_t01 < Sign::ZERO && vol_s_t12 > Sign::ZERO))
		return SimplexIntersectionType::DO_NOT_INTERSECT;

	Sign vol_s_t20 = Orient3D()(s0, s1, t2, t0);

	if ((vol_s_t12 > Sign::ZERO && vol_s_t20 < Sign::ZERO) ||
	    (vol_s_t12 < Sign::ZERO && vol_s_t20 > Sign::ZERO))
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	if ((vol_s_t20 > Sign::ZERO && vol_s_t01 < Sign::ZERO) ||
	    (vol_s_t20 < Sign::ZERO && vol_s_t01 > Sign::ZERO))
		return SimplexIntersectionType::DO_NOT_INTERSECT;

	return SimplexIntersectionType::INTERSECT;
}

template <typename Kernel>
SimplexIntersectionType
Triangle3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const NT *t0, const NT *t1, const NT *t2, const NT *s0, const NT *s1,
  int &n_max, const NT *t_min, const NT *t_perm) const
{
	OMC_INTER_PROFILE_INC_TOTAL(IntersectionNames::T3S3);

	bool s0_in_vertices =
	  (vec_equals_3d(s0, t0) || vec_equals_3d(s0, t1) || vec_equals_3d(s0, t2));
	bool s1_in_vertices =
	  (vec_equals_3d(s1, t0) || vec_equals_3d(s1, t1) || vec_equals_3d(s1, t2));
	if (s0_in_vertices && s1_in_vertices)
	{
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;
	}

	Sign vol_s0_t, vol_s1_t;
	if (t_min != nullptr)
	{
		vol_s0_t = Orient3D().with_cached_minors(t0, t1, t2, s0, t_min, t_perm);
		vol_s1_t = Orient3D().with_cached_minors(t0, t1, t2, s1, t_min, t_perm);
	}
	else
	{
		vol_s0_t = Orient3D()(t0, t1, t2, s0);
		vol_s1_t = Orient3D()(t0, t1, t2, s1);
	}

	if (vol_s0_t > Sign::ZERO && vol_s1_t > Sign::ZERO) // s is above t
	{
		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 0);
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	}
	if (vol_s0_t < Sign::ZERO && vol_s1_t < Sign::ZERO) // s is below t
	{
		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 1);
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	}
	if (vol_s0_t == Sign::ZERO && vol_s1_t == Sign::ZERO) // s and t are coplanar
	{
		if (n_max == -1)
			n_max = MaxCompInTriNormal()(t0, t1, t2);

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 2);

		// same code as the 2D version, I just copied it here....
		if ((Triangle3_Point3_DoInter().in_triangle(t0, t1, t2, s0, n_max) !=
		       PointInType::STRICTLY_OUTSIDE ||
		     Triangle3_Point3_DoInter().in_triangle(t0, t1, t2, s1, n_max) !=
		       PointInType::STRICTLY_OUTSIDE))
			return SimplexIntersectionType::INTERSECT;

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 3);

		switch (
		  Segment3_Segment3_DoInter().intersection_type(s0, s1, t0, t1, n_max))
		{
		case SimplexIntersectionType::INTERSECT:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::OVERLAP:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::SIMPLICIAL_COMPLEX:
			break;
		case SimplexIntersectionType::DO_NOT_INTERSECT:
			break;
		}

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 4);

		switch (
		  Segment3_Segment3_DoInter().intersection_type(s0, s1, t1, t2, n_max))
		{
		case SimplexIntersectionType::INTERSECT:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::OVERLAP:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::SIMPLICIAL_COMPLEX:
			break;
		case SimplexIntersectionType::DO_NOT_INTERSECT:
			break;
		}

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 5);

		switch (
		  Segment3_Segment3_DoInter().intersection_type(s0, s1, t2, t0, n_max))
		{
		case SimplexIntersectionType::INTERSECT:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::OVERLAP:
			return SimplexIntersectionType::INTERSECT;
		case SimplexIntersectionType::SIMPLICIAL_COMPLEX:
			break;
		case SimplexIntersectionType::DO_NOT_INTERSECT:
			break;
		}

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 6);

		return SimplexIntersectionType::DO_NOT_INTERSECT;
	}

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 7);

	// s intersects t (borders included), if the signs of the three tetrahedra
	// obtained combining s with the three edges of t are all equal

	// if one point is coplanar and coincides with a triangle vertex, then they
	// form a valid complex
	if (s0_in_vertices || s1_in_vertices)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	Sign vol_s_t01, vol_s_t12;
	vol_s_t01 = Orient3D()(s0, s1, t0, t1);
	vol_s_t12 = Orient3D()(s0, s1, t1, t2);

	if ((vol_s_t01 > Sign::ZERO && vol_s_t12 < Sign::ZERO) ||
	    (vol_s_t01 < Sign::ZERO && vol_s_t12 > Sign::ZERO))
		return SimplexIntersectionType::DO_NOT_INTERSECT;

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 8);

	Sign vol_s_t20 = Orient3D()(s0, s1, t2, t0);

	if ((vol_s_t12 > Sign::ZERO && vol_s_t20 < Sign::ZERO) ||
	    (vol_s_t12 < Sign::ZERO && vol_s_t20 > Sign::ZERO))
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	if ((vol_s_t20 > Sign::ZERO && vol_s_t01 < Sign::ZERO) ||
	    (vol_s_t20 < Sign::ZERO && vol_s_t01 > Sign::ZERO))
		return SimplexIntersectionType::DO_NOT_INTERSECT;

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3S3, 9);

	return SimplexIntersectionType::INTERSECT;
}

} // namespace OMC