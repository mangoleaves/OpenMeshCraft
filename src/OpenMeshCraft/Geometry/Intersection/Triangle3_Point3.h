#pragma once

#include "IntersectionUtils.h"
#include "Segment3_Point3.h"

namespace OMC {

/**
 * @brief Check if Triangle3 and Point3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Triangle3_Point3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT      = typename K::Vec3;
	using GPointT   = typename K::GPoint3;
	using TriangleT = typename K::Triangle3;

	using LessThan3D        = typename K::LessThan3D;
	using OrientOn2D        = typename K::OrientOn2D;
	using Orient3D          = typename K::Orient3D;
	using CollinearPoints3D = typename K::CollinearPoints3D;

	using Segment3_Point3_DoInter = Segment3_Point3_Do_Intersect<Kernel>;

public:
	/**
	 * @brief Check if triangle and point intersect.
	 * @note assume that triangle is not degenerate.
	 */
	bool operator()(const TriangleT &tri, const GPointT &pnt) const;

	/**
	 * @brief Check if triangle and point intersect.
	 * @param v0_v1_v2 three vertices of triangle
	 * @param p point
	 * @note assume that triangle is not degenerate.
	 */
	bool operator()(const GPointT &v0, const GPointT &v1, const GPointT &v2,
	                const GPointT &p) const;

	// Equivalent but faster. They assume that the
	// dominant normal component is n_max (see MaxComponentInTriangleNormal).

	bool operator()(const TriangleT &tri, const GPointT &pnt, int n_max) const;

	/*************************************************************************/
	/* NOTE below functions all assume that point is COPLANAR to triangle!!! */
	/*************************************************************************/

	/**
	 * @brief Get the relative position of point to triangle.
	 * @return Sign indicating relative position of point to triangle.
	 * @note assume that point is coplanar to triangle.
	 * @note assume that triangle is not degenerate.
	 */
	PointInType in_triangle(const TriangleT &tri, const GPointT &pnt) const;

	/**
	 * @brief Get the relative position of point to triangle.
	 * @param v0_v1_v2 three vertices of triangle
	 * @param p point
	 * @return Sign indicating relative position of point to triangle.
	 * @note assume that point is coplanar to triangle.
	 * @note assume that triangle is not degenerate.
	 */
	PointInType in_triangle(const GPointT &v0, const GPointT &v1,
	                        const GPointT &v2, const GPointT &p) const;

	/**
	 * @brief Get the relative position of point to triangle.
	 * @param v0_v1_v2 three vertices of triangle
	 * @param p point
	 * @return Sign indicating relative position of point to triangle.
	 * @note assume that point is coplanar to triangle.
	 * @note assume that triangle is not degenerate.
	 */
	PointInType in_triangle(const NT *v0, const NT *v1, const NT *v2,
	                        const NT *p) const;

	/**
	 * @brief Get the intersection type between triangle and point.
	 * @note assume that point is coplanar to triangle.
	 * @note assume that triangle is not degenerate.
	 */
	PointInSimplexType intersection_type(const TriangleT &tri,
	                                     const GPointT   &pnt) const;

	/**
	 * @brief Get the intersection type between triangle and point.
	 * @param v0_v1_v2 three vertices of triangle
	 * @param p point
	 * @note assume that point is coplanar to triangle.
	 * @note assume that triangle is not degenerate.
	 */
	PointInSimplexType intersection_type(const GPointT &v0, const GPointT &v1,
	                                     const GPointT &v2,
	                                     const GPointT &p) const;

	/**
	 * @brief Get the intersection type between triangle and point.
	 * @param v0_v1_v2 three vertices of triangle
	 * @param p point
	 * @note assume that point is coplanar to triangle.
	 * @note assume that triangle is not degenerate.
	 */
	PointInSimplexType intersection_type(const NT *v0, const NT *v1, const NT *v2,
	                                     const NT *p) const;

	// Equivalent but faster. They assume that the
	// dominant normal component is n_max (see MaxComponentInTriangleNormal).

	PointInType in_triangle(const TriangleT &tri, const GPointT &pnt,
	                        int n_max) const;

	PointInType in_triangle(const GPointT &v0, const GPointT &v1,
	                        const GPointT &v2, const GPointT &p, int n_max) const;

	PointInType in_triangle(const NT *v0, const NT *v1, const NT *v2, const NT *p,
	                        int n_max) const;

	PointInSimplexType intersection_type(const NT *v0, const NT *v1, const NT *v2,
	                                     const NT *p, int n_max) const;

protected:
	PointInType check_inout(Sign o1, Sign o2, Sign o3) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Triangle3_Point3.inl"
#endif