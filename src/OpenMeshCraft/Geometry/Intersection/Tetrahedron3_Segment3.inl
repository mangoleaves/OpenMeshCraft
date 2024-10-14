#pragma once

#include "Tetrahedron3_Segment3.h"

namespace OMC {

template <typename Kernel>
bool Tetrahedron3_Segment3_Do_Intersect<Kernel>::operator()(
  const TetrahedronT &tet, const SegmentT &seg) const
{
	const GPointT &t0 = tet.v0(), &t1 = tet.v1(), &t2 = tet.v2(), &t3 = tet.v3();
	const GPointT &s0 = seg.start(), &s1 = seg.end();
	return operator()(t0, t1, t2, t3, s0, s1);
}

template <typename Kernel>
bool Tetrahedron3_Segment3_Do_Intersect<Kernel>::operator()(
  const GPointT &t0, const GPointT &t1, const GPointT &t2, const GPointT &t3,
  const GPointT &s0, const GPointT &s1) const
{
	return intersection_type(t0, t1, t2, t3, s0, s1) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
SimplexIntersectionType
Tetrahedron3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const TetrahedronT &tet, const SegmentT &seg) const
{
	const GPointT &t0 = tet.v0(), &t1 = tet.v1(), &t2 = tet.v2(), &t3 = tet.v3();
	const GPointT &s0 = seg.start(), &s1 = seg.end();
	return intersection_type(t0, t1, t2, t3, s0, s1);
}

template <typename Kernel>
SimplexIntersectionType
Tetrahedron3_Segment3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &t0, const GPointT &t1, const GPointT &t2, const GPointT &t3,
  const GPointT &s0, const GPointT &s1) const
{
	// the tetrahedron is abbreviated as t, and the segment as s
	PointInSimplexType s0_in_tet =
	  Tetrahedron3_Point3_DoInter().intersection_type(t0, t1, t2, t3, s0);
	PointInSimplexType s1_in_tet =
	  Tetrahedron3_Point3_DoInter().intersection_type(t0, t1, t2, t3, s1);

	bool s0_in_verts = s0_in_tet >= PointInSimplexType::ON_VERT0 &&
	                   s0_in_tet <= PointInSimplexType::ON_VERT3;
	bool s1_in_verts = s1_in_tet >= PointInSimplexType::ON_VERT0 &&
	                   s1_in_tet <= PointInSimplexType::ON_VERT3;

	// s is an edge of t, form a valid complex
	if (s0_in_verts && s1_in_verts)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	// s intersects t at an endpoint
	if (s0_in_tet == PointInSimplexType::STRICTLY_INSIDE ||
	    s0_in_tet >= PointInSimplexType::ON_EDGE0 ||
	    s1_in_tet == PointInSimplexType::STRICTLY_INSIDE ||
	    s1_in_tet >= PointInSimplexType::ON_EDGE0)
		return SimplexIntersectionType::INTERSECT;

	// Now, the two endpoints of the segment are in the following cases:
	// - Both are outside the tetrahedron
	//   - Intersect
	//   - Do not intersect
	// - One is on a vertex, and the other is outside the tetrahedron.
	//   - Simplicial complex
	//   - Do not intersect

	SimplexIntersectionType seg_f0 =
	  Triangle3_Segment3_DoInter().intersection_type(t1, t3, t2, s0, s1);
	if (seg_f0 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	SimplexIntersectionType seg_f1 =
	  Triangle3_Segment3_DoInter().intersection_type(t0, t2, t3, s0, s1);
	if (seg_f1 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	SimplexIntersectionType seg_f2 =
	  Triangle3_Segment3_DoInter().intersection_type(t0, t3, t1, s0, s1);
	if (seg_f2 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	SimplexIntersectionType seg_f3 =
	  Triangle3_Segment3_DoInter().intersection_type(t0, t1, t2, s0, s1);
	if (seg_f3 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

#ifdef OMC_ENABLE_EXPENSIVE_ASSERT
	uint32_t sc_count = // simplicial complex count
	  (uint32_t)(seg_f0 == SimplexIntersectionType::SIMPLICIAL_COMPLEX) +
	  (uint32_t)(seg_f1 == SimplexIntersectionType::SIMPLICIAL_COMPLEX) +
	  (uint32_t)(seg_f2 == SimplexIntersectionType::SIMPLICIAL_COMPLEX) +
	  (uint32_t)(seg_f3 == SimplexIntersectionType::SIMPLICIAL_COMPLEX);
	uint32_t dni_count = // do not intersect count
	  (uint32_t)(seg_f0 == SimplexIntersectionType::DO_NOT_INTERSECT) +
	  (uint32_t)(seg_f1 == SimplexIntersectionType::DO_NOT_INTERSECT) +
	  (uint32_t)(seg_f2 == SimplexIntersectionType::DO_NOT_INTERSECT) +
	  (uint32_t)(seg_f3 == SimplexIntersectionType::DO_NOT_INTERSECT);
#endif

	// If one endpoint is on a vertex, then they form a simplicial complex
	if (s0_in_verts || s1_in_verts)
	{
		OMC_EXPENSIVE_ASSERT(
		  sc_count == 3,
		  "segment should form simplicial complex with exactly three faces.");
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;
	}
	else
	{
		OMC_EXPENSIVE_ASSERT(dni_count == 4,
		                     "segment should not intersect with any face.");
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	}
}

} // namespace OMC