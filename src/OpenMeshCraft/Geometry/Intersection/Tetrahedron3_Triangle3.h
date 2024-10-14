#pragma once

#include "IntersectionUtils.h"
#include "Tetrahedron3_Segment3.h"
#include "Triangle3_Triangle3.h"

namespace OMC {

/**
 * @brief Check if Tetrahedron3 and Triangle3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Tetrahedron3_Triangle3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT         = typename K::Vec3;
	using GPointT      = typename K::GPoint3;
	using SegmentT     = typename K::Segment3;
	using TriangleT    = typename K::Triangle3;
	using TetrahedronT = typename K::Tetrahedron3;

	using LessThan3D        = typename K::LessThan3D;
	using OrientOn2D        = typename K::OrientOn2D;
	using Orient3D          = typename K::Orient3D;
	using CollinearPoints3D = typename K::CollinearPoints3D;

	using Triangle3_Triangle3_DoInter = Triangle3_Triangle3_Do_Intersect<Kernel>;
	using Tetrahedron3_Point3_DoInter = Tetrahedron3_Point3_Do_Intersect<Kernel>;
	using Tetrahedron3_Segment3_DoInter =
	  Tetrahedron3_Segment3_Do_Intersect<Kernel>;

public:
	/**
	 * @brief Check if a tetrahedron and a triangle intersect.
	 * @note Assume that no tetrahedron or triangle is degenerate.
	 */
	bool operator()(const TetrahedronT &tet, const TriangleT &tri) const;

	/**
	 * @brief Check if a tetrahedron and a triangle intersect.
	 * @param t0_t1_t2_t3 four vertices of the tetrahedron.
	 * @param v0_v1_v2 three vertices of the triangle.
	 * @note Assume that no tetrahedron or triangle is degenerate.
	 */
	bool operator()(const GPointT &t0, const GPointT &t1, const GPointT &t2,
	                const GPointT &t3, const GPointT &v0, const GPointT &v1,
	                const GPointT &v2) const;

	/**
	 * @brief Get the intersection type between a tetrahedron and a triangle.
	 * @note Assume that no tetrahedron or triangle is degenerate.
	 */
	SimplexIntersectionType intersection_type(const TetrahedronT &tet,
	                                          const TriangleT    &tri) const;

	/**
	 * @brief Get the intersection type between a tetrahedron and a triangle.
	 * @param t0_t1_t2_t3 four vertices of the tetrahedron.
	 * @param v0_v1_v2 three vertices of the triangle.
	 * @note Assume that no tetrahedron or triangle is degenerate.
	 */
	SimplexIntersectionType
	intersection_type(const GPointT &t0, const GPointT &t1, const GPointT &t2,
	                  const GPointT &t3, const GPointT &v0, const GPointT &v1,
	                  const GPointT &v2) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Tetrahedron3_Triangle3.inl"
#endif