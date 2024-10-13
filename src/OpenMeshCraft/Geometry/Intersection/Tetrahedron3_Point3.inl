#pragma once

#include "Tetrahedron3_Point3.h"

namespace OMC {

template <typename Kernel>
bool Tetrahedron3_Point3_Do_Intersect<Kernel>::operator()(
  const TetrahedronT &tet, const GPointT &pnt) const
{
	const GPointT &t0 = tet.v0(), &t1 = tet.v1(), &t2 = tet.v2(), &t3 = tet.v3();
	return operator()(t0, t1, t2, t3, pnt);
}

template <typename Kernel>
bool Tetrahedron3_Point3_Do_Intersect<Kernel>::operator()(
  const GPointT &v0, const GPointT &v1, const GPointT &v2, const GPointT &v3,
  const GPointT &p) const
{
	return intersection_type(v0, v1, v2, v3, p) >=
	       PointInSimplexType::STRICTLY_INSIDE;
}

template <typename Kernel>
PointInSimplexType Tetrahedron3_Point3_Do_Intersect<Kernel>::intersection_type(
  const TetrahedronT &tet, const GPointT &pnt) const
{
	const GPointT &t0 = tet.v0(), &t1 = tet.v1(), &t2 = tet.v2(), &t3 = tet.v3();
	return intersection_type(t0, t1, t2, t3, pnt);
}

template <typename Kernel>
PointInSimplexType Tetrahedron3_Point3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &t0, const GPointT &t1, const GPointT &t2, const GPointT &t3,
  const GPointT &p) const
{
	Sign p_wrt_f0 = Orient3D()(t1, t3, t2, p); // face opposite to t0
	Sign p_wrt_f1 = Orient3D()(t0, t2, t3, p); // face opposite to t1
	Sign p_wrt_f2 = Orient3D()(t0, t3, t1, p); // face opposite to t2
	Sign p_wrt_f3 = Orient3D()(t0, t1, t2, p); // face opposite to t3

	OMC_ASSERT(p_wrt_f0 != Sign::ZERO || p_wrt_f1 != Sign::ZERO ||
	             p_wrt_f2 != Sign::ZERO || p_wrt_f3 != Sign::ZERO,
	           "Degenerate tetrahedron.");

	// four non-zero
	if ((p_wrt_f0 > Sign::ZERO && p_wrt_f1 > Sign::ZERO &&
	     p_wrt_f2 > Sign::ZERO && p_wrt_f3 > Sign::ZERO) ||
	    (p_wrt_f0 < Sign::ZERO && p_wrt_f1 < Sign::ZERO &&
	     p_wrt_f2 < Sign::ZERO && p_wrt_f3 < Sign::ZERO))
	{
		return PointInSimplexType::STRICTLY_INSIDE;
	}

	// clang-format off

  // three non-zero: on one vertex
	if (p_wrt_f1 == Sign::ZERO && p_wrt_f2 == Sign::ZERO && p_wrt_f3 == Sign::ZERO)
		return PointInSimplexType::ON_VERT0;
	if (p_wrt_f0 == Sign::ZERO && p_wrt_f2 == Sign::ZERO && p_wrt_f3 == Sign::ZERO)
		return PointInSimplexType::ON_VERT1;
	if (p_wrt_f0 == Sign::ZERO && p_wrt_f1 == Sign::ZERO && p_wrt_f3 == Sign::ZERO)
		return PointInSimplexType::ON_VERT2;
	if (p_wrt_f0 == Sign::ZERO && p_wrt_f1 == Sign::ZERO && p_wrt_f2 == Sign::ZERO)
		return PointInSimplexType::ON_VERT3;

  // two non-zero: colinear with one edge
	if (p_wrt_f2 == Sign::ZERO && p_wrt_f3 == Sign::ZERO && p_wrt_f0 == p_wrt_f1)
		return PointInSimplexType::ON_EDGE0;// edge: v0-v1
	if (p_wrt_f1 == Sign::ZERO && p_wrt_f3 == Sign::ZERO && p_wrt_f0 == p_wrt_f2)
		return PointInSimplexType::ON_EDGE1;// edge: v0-v2
	if (p_wrt_f1 == Sign::ZERO && p_wrt_f2 == Sign::ZERO && p_wrt_f0 == p_wrt_f3)
		return PointInSimplexType::ON_EDGE2;// edge: v0-v3
	if (p_wrt_f0 == Sign::ZERO && p_wrt_f3 == Sign::ZERO && p_wrt_f1 == p_wrt_f2)
		return PointInSimplexType::ON_EDGE3;// edge: v1-v2
	if (p_wrt_f0 == Sign::ZERO && p_wrt_f2 == Sign::ZERO && p_wrt_f1 == p_wrt_f3)
		return PointInSimplexType::ON_EDGE4;// edge: v1-v3
	if (p_wrt_f0 == Sign::ZERO && p_wrt_f1 == Sign::ZERO && p_wrt_f2 == p_wrt_f3)
		return PointInSimplexType::ON_EDGE5;// edge: v2-v3
  
  // one non-zero: coplanar with one face
  if (p_wrt_f0 == Sign::ZERO && p_wrt_f1 == p_wrt_f2 && p_wrt_f2 == p_wrt_f3)
    return PointInSimplexType::ON_FACE0;// face: v1-v2-v3
  if (p_wrt_f1 == Sign::ZERO && p_wrt_f0 == p_wrt_f2 && p_wrt_f2 == p_wrt_f3)
    return PointInSimplexType::ON_FACE1;// face: v0-v2-v3
  if (p_wrt_f2 == Sign::ZERO && p_wrt_f0 == p_wrt_f1 && p_wrt_f1 == p_wrt_f3)
    return PointInSimplexType::ON_FACE2;// face: v0-v1-v3
  if (p_wrt_f3 == Sign::ZERO && p_wrt_f0 == p_wrt_f1 && p_wrt_f1 == p_wrt_f2)
    return PointInSimplexType::ON_FACE3;// face: v0-v1-v2

	// clang-format on
	return PointInSimplexType::STRICTLY_OUTSIDE;
}

} // namespace OMC