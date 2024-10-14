#pragma once

#include "Tetrahedron3_Triangle3.h"

namespace OMC {

template <typename Kernel>
bool Tetrahedron3_Triangle3_Do_Intersect<Kernel>::operator()(
  const TetrahedronT &tet, const TriangleT &tri) const
{
	const GPointT &t0 = tet.v0(), &t1 = tet.v1(), &t2 = tet.v2(), &t3 = tet.v3();
	const GPointT &v0 = tri.v0(), &v1 = tri.v1(), &v2 = tri.v2();
	return operator()(t0, t1, t2, t3, v0, v1, v2);
}

template <typename Kernel>
bool Tetrahedron3_Triangle3_Do_Intersect<Kernel>::operator()(
  const GPointT &t0, const GPointT &t1, const GPointT &t2, const GPointT &t3,
  const GPointT &v0, const GPointT &v1, const GPointT &v2) const
{
	return intersection_type(t0, t1, t2, t3, v0, v1, v2) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
SimplexIntersectionType
Tetrahedron3_Triangle3_Do_Intersect<Kernel>::intersection_type(
  const TetrahedronT &tet, const TriangleT &tri) const
{
	const GPointT &t0 = tet.v0(), &t1 = tet.v1(), &t2 = tet.v2(), &t3 = tet.v3();
	const GPointT &v0 = tri.v0(), &v1 = tri.v1(), &v2 = tri.v2();
	return intersection_type(t0, t1, t2, t3, v0, v1, v2);
}

template <typename Kernel>
SimplexIntersectionType
Tetrahedron3_Triangle3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &t0, const GPointT &t1, const GPointT &t2, const GPointT &t3,
  const GPointT &v0, const GPointT &v1, const GPointT &v2) const
{
	// the tetrahedron is abbreviated as tet, and the triangle as tri.
	// the vertices of tet is called ti, and the vertices of tri is called vi.
	PointInSimplexType v0_in_tet =
	  Tetrahedron3_Point3_DoInter().intersection_type(t0, t1, t2, t3, v0);
	PointInSimplexType v1_in_tet =
	  Tetrahedron3_Point3_DoInter().intersection_type(t0, t1, t2, t3, v1);
	PointInSimplexType v2_in_tet =
	  Tetrahedron3_Point3_DoInter().intersection_type(t0, t1, t2, t3, v2);

	bool v0_in_verts = v0_in_tet >= PointInSimplexType::ON_VERT0 &&
	                   v0_in_tet <= PointInSimplexType::ON_VERT3;
	bool v1_in_verts = v1_in_tet >= PointInSimplexType::ON_VERT0 &&
	                   v1_in_tet <= PointInSimplexType::ON_VERT3;
	bool v2_in_verts = v2_in_tet >= PointInSimplexType::ON_VERT0 &&
	                   v2_in_tet <= PointInSimplexType::ON_VERT3;

	// tri is a face of tet, form a valid complex.
	if (v0_in_verts && v1_in_verts && v2_in_verts)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	// tri intersects tet at a triangle vertex
	if (v0_in_tet == PointInSimplexType::STRICTLY_INSIDE ||
	    v0_in_tet >= PointInSimplexType::ON_EDGE0 ||
	    v1_in_tet == PointInSimplexType::STRICTLY_INSIDE ||
	    v1_in_tet >= PointInSimplexType::ON_EDGE0 ||
	    v2_in_tet == PointInSimplexType::STRICTLY_INSIDE ||
	    v2_in_tet >= PointInSimplexType::ON_EDGE0)
		return SimplexIntersectionType::INTERSECT;

	// Now, the three vertices the triangle are in the following cases:
	// - All are outside the tetrahedron
	//   - Intersect
	//   - Do not intersect
	// - One or two are on tet's vertices, and the others are outside the
	//   tetrahedron.
	//   - Simplicial complex
	//   - Do not intersect

	SimplexIntersectionType tri_f0 =
	  Triangle3_Triangle3_DoInter().intersection_type(t1, t3, t2, v0, v1, v2);
	if (tri_f0 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	SimplexIntersectionType tri_f1 =
	  Triangle3_Triangle3_DoInter().intersection_type(t0, t2, t3, v0, v1, v2);
	if (tri_f1 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	SimplexIntersectionType tri_f2 =
	  Triangle3_Triangle3_DoInter().intersection_type(t0, t3, t1, v0, v1, v2);
	if (tri_f2 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	SimplexIntersectionType tri_f3 =
	  Triangle3_Triangle3_DoInter().intersection_type(t0, t1, t2, v0, v1, v2);
	if (tri_f3 == SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

#ifdef OMC_ENABLE_EXPENSIVE_ASSERT
	uint32_t sc_count = // simplicial complex count
	  (uint32_t)(tri_f0 == SimplexIntersectionType::SIMPLICIAL_COMPLEX) +
	  (uint32_t)(tri_f1 == SimplexIntersectionType::SIMPLICIAL_COMPLEX) +
	  (uint32_t)(tri_f2 == SimplexIntersectionType::SIMPLICIAL_COMPLEX) +
	  (uint32_t)(tri_f3 == SimplexIntersectionType::SIMPLICIAL_COMPLEX);
	uint32_t dni_count = // do not intersect count
	  (uint32_t)(tri_f0 == SimplexIntersectionType::DO_NOT_INTERSECT) +
	  (uint32_t)(tri_f1 == SimplexIntersectionType::DO_NOT_INTERSECT) +
	  (uint32_t)(tri_f2 == SimplexIntersectionType::DO_NOT_INTERSECT) +
	  (uint32_t)(tri_f3 == SimplexIntersectionType::DO_NOT_INTERSECT);
#endif

	// If one tri's vertex is on a tet's vertex, then they form a simplicial
	// complex
	if (v0_in_verts || v1_in_verts || v2_in_verts)
	{
		OMC_EXPENSIVE_ASSERT(
		  sc_count >= 3,
		  "triangle should form simplicial complex with three or four faces.");
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;
	}
	else
	{
		OMC_EXPENSIVE_ASSERT(dni_count == 4,
		                     "triangle should not intersect with any face.");
		return SimplexIntersectionType::DO_NOT_INTERSECT;
	}
}

} // namespace OMC