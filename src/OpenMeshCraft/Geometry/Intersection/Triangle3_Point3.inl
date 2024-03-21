#pragma once

#include "Triangle3_Point3.h"

namespace OMC {

template <typename Kernel>
bool Triangle3_Point3_Do_Intersect<Kernel>::operator()(const TriangleT &tri,
                                                       const GPointT &pnt) const
{
	// I copy intersection_type here, just ignore checking point in segment.
	const GPointT &t0 = tri.v0(), &t1 = tri.v1(), &t2 = tri.v2();

	return operator()(t0, t1, t2, pnt);
}

template <typename Kernel>
bool Triangle3_Point3_Do_Intersect<Kernel>::operator()(const GPointT &t0,
                                                       const GPointT &t1,
                                                       const GPointT &t2,
                                                       const GPointT &pnt) const
{
	// I copy intersection_type here, just ignore checking point in segment.

	if (Orient3D()(t0, t1, t2, pnt) != Sign::ZERO)
		return false;

	// test for point in vert
	if (LessThan3D().coincident(pnt, t0) || LessThan3D().coincident(pnt, t1) ||
	    LessThan3D().coincident(pnt, t2))
		return true;

	Sign o1, o2, o3;

	o1 = OrientOn2D().on_yz(t0, t1, pnt);
	o2 = OrientOn2D().on_yz(t1, t2, pnt);
	o3 = OrientOn2D().on_yz(t2, t0, pnt);

	if (o1 != Sign::ZERO && o2 != Sign::ZERO && o3 != Sign::ZERO)
	{
		return (o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
		       (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO);
	}

	o1 = OrientOn2D().on_zx(t0, t1, pnt);
	o2 = OrientOn2D().on_zx(t1, t2, pnt);
	o3 = OrientOn2D().on_zx(t2, t0, pnt);

	if (o1 != Sign::ZERO && o2 != Sign::ZERO && o3 != Sign::ZERO)
	{
		return (o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
		       (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO);
	}

	o1 = OrientOn2D().on_xy(t0, t1, pnt);
	o2 = OrientOn2D().on_xy(t1, t2, pnt);
	o3 = OrientOn2D().on_xy(t2, t0, pnt);

	return (o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
	       (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO);
}

template <typename Kernel>
PointInType
Triangle3_Point3_Do_Intersect<Kernel>::in_triangle(const TriangleT &tri,
                                                   const GPointT   &pnt) const
{
	return in_triangle(tri.v0(), tri.v1(), tri.v2(), pnt);
}

template <typename Kernel>
PointInType Triangle3_Point3_Do_Intersect<Kernel>::in_triangle(
  const GPointT &t0, const GPointT &t1, const GPointT &t2,
  const GPointT &p) const
{
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, p) == Sign::ZERO,
	                     "point is not on triangle.");
	// test for point in vert
	if (LessThan3D().coincident(p, t0) || LessThan3D().coincident(p, t1) ||
	    LessThan3D().coincident(p, t2))
		return PointInType::ON_BOUNDARY;

	Sign o1, o2, o3;

	o1 = OrientOn2D().on_yz(t0, t1, p);
	o2 = OrientOn2D().on_yz(t1, t2, p);
	o3 = OrientOn2D().on_yz(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
		return check_inout(o1, o2, o3);

	o1 = OrientOn2D().on_zx(t0, t1, p);
	o2 = OrientOn2D().on_zx(t1, t2, p);
	o3 = OrientOn2D().on_zx(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
		return check_inout(o1, o2, o3);

	o1 = OrientOn2D().on_xy(t0, t1, p);
	o2 = OrientOn2D().on_xy(t1, t2, p);
	o3 = OrientOn2D().on_xy(t2, t0, p);

	return check_inout(o1, o2, o3);
}

template <typename Kernel>
PointInType Triangle3_Point3_Do_Intersect<Kernel>::in_triangle(
  const NT *t0, const NT *t1, const NT *t2, const NT *p) const
{
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, p) == Sign::ZERO,
	                     "point is not on triangle.");
	// test for point in vert
	if (vec_equals_3d(p, t0) || vec_equals_3d(p, t1) || vec_equals_3d(p, t2))
		return PointInType::ON_BOUNDARY;

	Sign o1, o2, o3;

	o1 = OrientOn2D().on_yz(t0, t1, p);
	o2 = OrientOn2D().on_yz(t1, t2, p);
	o3 = OrientOn2D().on_yz(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
		return check_inout(o1, o2, o3);

	o1 = OrientOn2D().on_zx(t0, t1, p);
	o2 = OrientOn2D().on_zx(t1, t2, p);
	o3 = OrientOn2D().on_zx(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
		return check_inout(o1, o2, o3);

	o1 = OrientOn2D().on_xy(t0, t1, p);
	o2 = OrientOn2D().on_xy(t1, t2, p);
	o3 = OrientOn2D().on_xy(t2, t0, p);

	return check_inout(o1, o2, o3);
}

template <typename Kernel>
PointInSimplexType Triangle3_Point3_Do_Intersect<Kernel>::intersection_type(
  const TriangleT &tri, const GPointT &pnt) const
{
	return intersection_type(tri.v0(), tri.v1(), tri.v2(), pnt);
}

template <typename Kernel>
PointInSimplexType Triangle3_Point3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &t0, const GPointT &t1, const GPointT &t2,
  const GPointT &p) const
{
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, p) == Sign::ZERO,
	                     "point is not on triangle.");

	// test for point in vert
	if (LessThan3D().coincident(p, t0))
		return PointInSimplexType::ON_VERT0;
	if (LessThan3D().coincident(p, t1))
		return PointInSimplexType::ON_VERT1;
	if (LessThan3D().coincident(p, t2))
		return PointInSimplexType::ON_VERT2;

	// test for point in edge in 3D
	if (Segment3_Point3_DoInter().intersection_type(t0, t1, p) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE0;
	if (Segment3_Point3_DoInter().intersection_type(t1, t2, p) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE1;
	if (Segment3_Point3_DoInter().intersection_type(t2, t0, p) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE2;

	// test for the interior: project t on XYZ and, if the check is never false in
	// any of the projections, then p must be inside t.
	// (I copy partial of point_in_triangle_2d here.)
	Sign o1, o2, o3;

	o1 = OrientOn2D().on_yz(t0, t1, p);
	o2 = OrientOn2D().on_yz(t1, t2, p);
	o3 = OrientOn2D().on_yz(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
	{
		if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
		    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
			return PointInSimplexType::STRICTLY_INSIDE;
	}

	o1 = OrientOn2D().on_zx(t0, t1, p);
	o2 = OrientOn2D().on_zx(t1, t2, p);
	o3 = OrientOn2D().on_zx(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
	{
		if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
		    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
			return PointInSimplexType::STRICTLY_INSIDE;
	}

	o1 = OrientOn2D().on_xy(t0, t1, p);
	o2 = OrientOn2D().on_xy(t1, t2, p);
	o3 = OrientOn2D().on_xy(t2, t0, p);

	if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
	    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
		return PointInSimplexType::STRICTLY_INSIDE;

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
PointInSimplexType Triangle3_Point3_Do_Intersect<Kernel>::intersection_type(
  const NT *t0, const NT *t1, const NT *t2, const NT *p) const
{
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, p) == Sign::ZERO,
	                     "point is not on triangle.");

	// test for point in vert
	if (vec_equals_3d(p, t0))
		return PointInSimplexType::ON_VERT0;
	if (vec_equals_3d(p, t1))
		return PointInSimplexType::ON_VERT1;
	if (vec_equals_3d(p, t2))
		return PointInSimplexType::ON_VERT2;

	// test for point in edge in 3D
	if (Segment3_Point3_DoInter().intersection_type(t0, t1, p) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE0;
	if (Segment3_Point3_DoInter().intersection_type(t1, t2, p) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE1;
	if (Segment3_Point3_DoInter().intersection_type(t2, t0, p) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE2;

	// test for the interior: project t on XYZ and, if the check is never false in
	// any of the projections, then p must be inside t.
	// (I copy partial of point_in_triangle_2d here.)
	Sign o1, o2, o3;

	o1 = OrientOn2D().on_yz(t0, t1, p);
	o2 = OrientOn2D().on_yz(t1, t2, p);
	o3 = OrientOn2D().on_yz(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
	{
		if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
		    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
			return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}

	o1 = OrientOn2D().on_zx(t0, t1, p);
	o2 = OrientOn2D().on_zx(t1, t2, p);
	o3 = OrientOn2D().on_zx(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
	{
		if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
		    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
			return PointInSimplexType::STRICTLY_INSIDE;
		else
			return PointInSimplexType::STRICTLY_OUTSIDE;
	}

	o1 = OrientOn2D().on_xy(t0, t1, p);
	o2 = OrientOn2D().on_xy(t1, t2, p);
	o3 = OrientOn2D().on_xy(t2, t0, p);

	if (o1 != Sign::ZERO || o2 != Sign::ZERO || o3 != Sign::ZERO)
	{
		if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
		    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
			return PointInSimplexType::STRICTLY_INSIDE;
	}

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
bool Triangle3_Point3_Do_Intersect<Kernel>::operator()(const TriangleT &tri,
                                                       const GPointT   &pnt,
                                                       int n_max) const
{
	// I copy intersection_type here, just ignore checking point in segment.
	const GPointT &t0 = tri.v0(), &t1 = tri.v1(), &t2 = tri.v2();
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, pnt) == Sign::ZERO,
	                     "point is not on triangle.");

	if (Orient3D()(t0, t1, t2, pnt) != Sign::ZERO)
		return false;

	// test for point in vert
	if (LessThan3D().coincident(pnt, t0) || LessThan3D().coincident(pnt, t1) ||
	    LessThan3D().coincident(pnt, t2))
		return true;

	Sign o1 = OrientOn2D()(t0, t1, pnt, n_max);
	Sign o2 = OrientOn2D()(t1, t2, pnt, n_max);
	Sign o3 = OrientOn2D()(t2, t0, pnt, n_max);

	return (o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
	       (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO);
}

template <typename Kernel>
PointInType Triangle3_Point3_Do_Intersect<Kernel>::in_triangle(
  const TriangleT &tri, const GPointT &pnt, int n_max) const
{
	return in_triangle(tri.v0(), tri.v1(), tri.v2(), pnt, n_max);
}

template <typename Kernel>
PointInType Triangle3_Point3_Do_Intersect<Kernel>::in_triangle(
  const GPointT &t0, const GPointT &t1, const GPointT &t2, const GPointT &p,
  int n_max) const
{
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, p) == Sign::ZERO,
	                     "point is not on triangle.");
	// I copy intersection_type here, just ignore checking point in segment.
	Sign o1 = OrientOn2D()(t0, t1, p, n_max);
	Sign o2 = OrientOn2D()(t1, t2, p, n_max);
	Sign o3 = OrientOn2D()(t2, t0, p, n_max);

	return check_inout(o1, o2, o3);
}

template <typename Kernel>
PointInType Triangle3_Point3_Do_Intersect<Kernel>::in_triangle(
  const NT *t0, const NT *t1, const NT *t2, const NT *p, int n_max) const
{
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, p) == Sign::ZERO,
	                     "point is not on triangle.");
	// test for point in vert
	if (vec_equals_3d(p, t0) || vec_equals_3d(p, t1) || vec_equals_3d(p, t2))
		return PointInType::ON_BOUNDARY;

	Sign o1 = OrientOn2D()(t0, t1, p, n_max);
	Sign o2 = OrientOn2D()(t1, t2, p, n_max);
	Sign o3 = OrientOn2D()(t2, t0, p, n_max);

	return check_inout(o1, o2, o3);
}

template <typename Kernel>
PointInSimplexType Triangle3_Point3_Do_Intersect<Kernel>::intersection_type(
  const NT *t0, const NT *t1, const NT *t2, const NT *p, int n_max) const
{
	OMC_EXPENSIVE_ASSERT(Orient3D()(t0, t1, t2, p) == Sign::ZERO,
	                     "point is not on triangle.");
	// test for point in vert
	if (vec_equals_3d(p, t0))
		return PointInSimplexType::ON_VERT0;
	if (vec_equals_3d(p, t1))
		return PointInSimplexType::ON_VERT1;
	if (vec_equals_3d(p, t2))
		return PointInSimplexType::ON_VERT2;

	// test for point in edge in 3D
	if (Segment3_Point3_DoInter().intersection_type(t0, t1, p, n_max) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE0;
	if (Segment3_Point3_DoInter().intersection_type(t1, t2, p, n_max) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE1;
	if (Segment3_Point3_DoInter().intersection_type(t2, t0, p, n_max) ==
	    PointInSimplexType::STRICTLY_INSIDE)
		return PointInSimplexType::ON_EDGE2;

	// test for the interior: project t on XYZ and, if the check is never false in
	// any of the projections, then p must be inside t

	// test for the interior: project t on XYZ and, if the check is never false in
	// any of the projections, then p must be inside t.
	// (I copy partial of point_in_triangle_2d here.)
	Sign o1 = OrientOn2D()(t0, t1, p, n_max);
	Sign o2 = OrientOn2D()(t1, t2, p, n_max);
	Sign o3 = OrientOn2D()(t2, t0, p, n_max);

	if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
	    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
		return PointInSimplexType::STRICTLY_INSIDE;
	return PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
PointInType Triangle3_Point3_Do_Intersect<Kernel>::check_inout(Sign o1, Sign o2,
                                                               Sign o3) const
{
	if ((o1 > Sign::ZERO && o2 > Sign::ZERO && o3 > Sign::ZERO) ||
	    (o1 < Sign::ZERO && o2 < Sign::ZERO && o3 < Sign::ZERO))
		return PointInType::STRICTLY_INSIDE;

	if ((o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
	    (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO))
		return PointInType::ON_BOUNDARY;

	return PointInType::STRICTLY_OUTSIDE;
}

} // namespace OMC