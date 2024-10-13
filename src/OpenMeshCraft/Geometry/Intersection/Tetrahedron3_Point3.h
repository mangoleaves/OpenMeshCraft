#pragma once

#include "IntersectionUtils.h"
#include "Triangle3_Point3.h"

namespace OMC {

/**
 * @brief Check if Tetrahedron3 and Point3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Tetrahedron3_Point3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT         = typename K::Vec3;
	using GPointT      = typename K::GPoint3;
	using TetrahedronT = typename K::Tetrahedron3;

	using LessThan3D        = typename K::LessThan3D;
	using OrientOn2D        = typename K::OrientOn2D;
	using Orient3D          = typename K::Orient3D;
	using CollinearPoints3D = typename K::CollinearPoints3D;

public:
	/**
	 * @brief Check if tetrahedron and point intersect.
	 * @note assume that tetrahedron is not degenerate.
	 */
	bool operator()(const TetrahedronT &tet, const GPointT &pnt) const;

	/**
	 * @brief Check if tetrahedron and point intersect.
	 * @param v0_v1_v2_v3 four vertices of tetrahedron
	 * @param p point
	 * @note assume that tetrahedron is not degenerate.
	 */
	bool operator()(const GPointT &v0, const GPointT &v1, const GPointT &v2,
	                const GPointT &v3, const GPointT &p) const;

	/**
	 * @brief Get the intersection type between tetrahedron and point.
	 * @note assume that tetrahedron is not degenerate.
	 */
	PointInSimplexType intersection_type(const TetrahedronT &tet,
	                                     const GPointT      &pnt) const;

	/**
	 * @brief Get the intersection type between tetrahedron and point.
	 * @param v0_v1_v2_v3 four vertices of tetrahedron
	 * @param p point
	 * @note assume that tetrahedron is not degenerate.
	 */
	PointInSimplexType intersection_type(const GPointT &v0, const GPointT &v1,
	                                     const GPointT &v2, const GPointT &v3,
	                                     const GPointT &p) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Tetrahedron3_Point3.inl"
#endif