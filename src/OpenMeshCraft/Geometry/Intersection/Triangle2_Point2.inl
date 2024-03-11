#pragma once

#include "Triangle2_Point2.h"

namespace OMC {

template <typename Kernel>
bool Triangle2_Point2_Do_Intersect<Kernel>::operator()(const TriangleT &tri,
                                                       const GPointT &pnt) const
{
	return intersection_type(tri, pnt) >= PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
PointInType
Triangle2_Point2_Do_Intersect<Kernel>::in_triangle(const TriangleT &tri,
                                                   const GPointT   &pnt) const
{
	PointInSimplexType type = intersection_type(tri, pnt);
	if (type == PointInSimplexType::STRICTLY_INSIDE)
		return PointInType::STRICTLY_INSIDE;
	else if (type == PointInSimplexType::STRICTLY_OUTSIDE)
		return PointInType::STRICTLY_OUTSIDE;
	else
		return PointInType::ON_BOUNDARY;
}

template <typename Kernel>
PointInSimplexType Triangle2_Point2_Do_Intersect<Kernel>::intersection_type(
  const TriangleT &tri, const GPointT &pnt) const
{
	return intersection_type(tri.v0(), tri.v1(), tri.v2(), pnt);
}

template <typename Kernel>
PointInSimplexType Triangle2_Point2_Do_Intersect<Kernel>::intersection_type(
  const GPointT &t0, const GPointT &t1, const GPointT &t2,
  const GPointT &p) const
{
	if (LessThan2D().coincident(p, t0))
		return PointInSimplexType::ON_VERT0;
	if (LessThan2D().coincident(p, t1))
		return PointInSimplexType::ON_VERT1;
	if (LessThan2D().coincident(p, t2))
		return PointInSimplexType::ON_VERT2;

	Sign o1 = Orient2D()(t0, t1, p);
	Sign o2 = Orient2D()(t1, t2, p);
	Sign o3 = Orient2D()(t2, t0, p);

	bool hit = (o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
	           (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO);

	if (hit)
	{
		if (o1 == Sign::ZERO)
			return PointInSimplexType::ON_EDGE0;
		if (o2 == Sign::ZERO)
			return PointInSimplexType::ON_EDGE1;
		if (o3 == Sign::ZERO)
			return PointInSimplexType::ON_EDGE2;

		return PointInSimplexType::STRICTLY_INSIDE;
	}

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

template <typename Kernel>
PointInSimplexType Triangle2_Point2_Do_Intersect<Kernel>::intersection_type(
  const NT *t0, const NT *t1, const NT *t2, const NT *p) const
{
	if (vec_equals_2d(p, t0))
		return PointInSimplexType::ON_VERT0;
	if (vec_equals_2d(p, t1))
		return PointInSimplexType::ON_VERT1;
	if (vec_equals_2d(p, t2))
		return PointInSimplexType::ON_VERT2;

	Sign o1 = Orient2D()(t0, t1, p);
	Sign o2 = Orient2D()(t1, t2, p);
	Sign o3 = Orient2D()(t2, t0, p);

	bool hit = (o1 >= Sign::ZERO && o2 >= Sign::ZERO && o3 >= Sign::ZERO) ||
	           (o1 <= Sign::ZERO && o2 <= Sign::ZERO && o3 <= Sign::ZERO);

	if (hit)
	{
		if (o1 == Sign::ZERO)
			return PointInSimplexType::ON_EDGE0;
		if (o2 == Sign::ZERO)
			return PointInSimplexType::ON_EDGE1;
		if (o3 == Sign::ZERO)
			return PointInSimplexType::ON_EDGE2;

		return PointInSimplexType::STRICTLY_INSIDE;
	}

	return PointInSimplexType::STRICTLY_OUTSIDE;
}

} // namespace OMC