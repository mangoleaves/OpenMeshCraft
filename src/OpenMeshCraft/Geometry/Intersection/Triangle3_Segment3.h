#pragma once

#include "IntersectionUtils.h"
#include "Segment3_Segment3.h"
#include "Triangle3_Point3.h"

namespace OMC {

/**
 * @brief Check if Triangle3 and Segment3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Triangle3_Segment3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT      = typename K::Vec3;
	using GPointT   = typename K::GPoint3;
	using SegmentT  = typename K::Segment3;
	using TriangleT = typename K::Triangle3;

	using LessThan3D         = typename K::LessThan3D;
	using OrientOn2D         = typename K::OrientOn2D;
	using Orient3D           = typename K::Orient3D;
	using MaxCompInTriNormal = typename K::MaxCompInTriNormal;

	using Triangle3_Point3_DoInter  = Triangle3_Point3_Do_Intersect<Kernel>;
	using Segment3_Segment3_DoInter = Segment3_Segment3_Do_Intersect<Kernel>;

public:
	/**
	 * @brief Check if triangle and segment intersect.
	 * @note Assume no one is degenerate.
	 */
	bool operator()(const TriangleT &tri, const SegmentT &seg) const;

	/**
	 * @brief Check if triangle and segment intersect.
	 * @param v1_v2_v3 three vertices of triangle
	 * @param s1_s2 two vertices of segment
	 * @note Assume no one is degenerate.
	 */
	bool operator()(const GPointT &v1, const GPointT &v2, const GPointT &v3,
	                const GPointT &s1, const GPointT &s2) const;

	/**
	 * @brief Check if the closure of triangle and segment cross at a single point
	 * (overlap/coplanar is excluded).
	 * @return TRUE if they intersect at a point.
	 * @note Assume no one is degenerate.
	 */
	bool cross(const TriangleT &tri, const SegmentT &seg) const;

	/**
	 * @brief Check if the interior of triangle and segment cross at a single
	 * point (overlap/coplanar is excluded).
	 * @return TRUE if they intersect at a point.
	 * @note Assume no one is degenerate.
	 */
	bool cross_inner(const TriangleT &tri, const SegmentT &seg) const;

	/**
	 * @brief Check if the closure of triangle and segment cross at a single point
	 * (overlap/coplanar is excluded).
	 * @param v1_v2_v3 three vertices of triangle
	 * @param s1_s2 two vertices of segment
	 * @return TRUE if they intersect at a point.
	 * @note Assume no one is degenerate.
	 */
	bool cross(const GPointT &v1, const GPointT &v2, const GPointT &v3,
	           const GPointT &s1, const GPointT &s2) const;

	/**
	 * @brief Check if the interior of triangle and segment cross at a single
	 * point (overlap/coplanar is excluded).
	 * @param v1_v2_v3 three vertices of triangle
	 * @param s1_s2 two vertices of segment
	 * @return TRUE if they intersect at a point.
	 * @note Assume no one is degenerate.
	 */
	bool cross_inner(const GPointT &v1, const GPointT &v2, const GPointT &v3,
	                 const GPointT &s1, const GPointT &s2) const;

	/**
	 * @brief Get the intersection type between triangle and segment.
	 * @note Assume no one is degenerate.
	 */
	SimplexIntersectionType intersection_type(const TriangleT &tri,
	                                          const SegmentT  &seg) const;

	/**
	 * @brief Get the intersection type between triangle and segment.
	 * @param v1_v2_v3 three vertices of triangle
	 * @param s1_s2 two vertices of segment
	 * @note Assume no one is degenerate.
	 */
	SimplexIntersectionType
	intersection_type(const GPointT &v1, const GPointT &v2, const GPointT &v3,
	                  const GPointT &s1, const GPointT &s2) const;

	/**
	 * @brief Get the intersection type between triangle and segment.
	 * @param [in] v1_v2_v3 three vertices of triangle
	 * @param [in] s1_s2 two vertices of segment
	 * @param [in] t_min cached minor for triangle
	 * @param [in] t_perm cached permutation for triangle
	 * @param [inout] n_max max component of triangle normal. If input is -1, it
	 * will be updated to 0/1/2 if it is calculated inernally. If input is 0/1/2,
	 * it will be used.
	 * @note Assume no one is degenerate.
	 */
	SimplexIntersectionType intersection_type(const NT *v1, const NT *v2,
	                                          const NT *v3, const NT *s1,
	                                          const NT *s2, int &n_max,
	                                          const NT *t_min  = nullptr,
	                                          const NT *t_perm = nullptr) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Triangle3_Segment3.inl"
#endif