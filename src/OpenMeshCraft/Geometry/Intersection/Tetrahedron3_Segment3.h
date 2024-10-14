#pragma once

#include "IntersectionUtils.h"
#include "Segment3_Segment3.h"
#include "Tetrahedron3_Point3.h"
#include "Triangle3_Segment3.h"

namespace OMC {

/**
 * @brief Check if Tetrahedron3 and Segment3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Tetrahedron3_Segment3_Do_Intersect
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

	using Tetrahedron3_Point3_DoInter = Tetrahedron3_Point3_Do_Intersect<Kernel>;
	using Triangle3_Segment3_DoInter  = Triangle3_Segment3_Do_Intersect<Kernel>;
	using Segment3_Segment3_DoInter   = Segment3_Segment3_Do_Intersect<Kernel>;

public:
	/**
	 * @brief Check if tetrahedron and segment intersect.
	 * @note Assume no one is degenerate.
	 */
	bool operator()(const TetrahedronT &tet, const SegmentT &seg) const;

	/**
	 * @brief Check if tetrahedron and segment intersect.
	 * @param v0_v1_v2_v3 four vertices of tetrahedron
	 * @param s0_s1 two vertices of segment
	 * @note Assume no one is degenerate.
	 */
	bool operator()(const GPointT &t0, const GPointT &t1, const GPointT &t2,
	                const GPointT &t3, const GPointT &s0,
	                const GPointT &s1) const;

	/**
	 * @brief Get the intersection type between tetrahedron and segment.
	 * @note Assume no one is degenerate.
	 */
	SimplexIntersectionType intersection_type(const TetrahedronT &tet,
	                                          const SegmentT     &seg) const;

	/**
	 * @brief Get the intersection type between tetrahedron and segment.
	 * @param v0_v1_v2_v3 four vertices of tetrahedron
	 * @param s0_s1 two vertices of segment
	 * @note Assume no one is degenerate.
	 */
	SimplexIntersectionType
	intersection_type(const GPointT &t0, const GPointT &t1, const GPointT &t2,
	                  const GPointT &t3, const GPointT &s0,
	                  const GPointT &s1) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Tetrahedron3_Segment3.inl"
#endif