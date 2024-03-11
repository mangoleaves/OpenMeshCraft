#pragma once

#include "IntersectionUtils.h"
#include "Triangle3_Segment3.h"

namespace OMC {

/**
 * @brief Check if Triangle3 and Triangle3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Triangle3_Triangle3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT      = typename K::Vec3;
	using GPointT   = typename K::GPoint3;
	using SegmentT  = typename K::Segment3;
	using TriangleT = typename K::Triangle3;

	using LessThan3D        = typename K::LessThan3D;
	using OrientOn2D        = typename K::OrientOn2D;
	using Orient3D          = typename K::Orient3D;
	using CollinearPoints3D = typename K::CollinearPoints3D;

	using Triangle3_Segment3_DoInter = Triangle3_Segment3_Do_Intersect<Kernel>;

public:
	/**
	 * @brief Check if two triangles intersect.
	 * @note Assume that no triangle is degenerate.
	 */
	bool operator()(const TriangleT &tri0, const TriangleT &tri1) const;

	/**
	 * @brief Check if two triangles intersect.
	 * @param t00_t01_t02 three vertices of the first triangle.
	 * @param t10_t11_t12 three vertices of the second triangle.
	 * @note Assume that no triangle is degenerate.
	 */
	bool operator()(const GPointT &t00, const GPointT &t01, const GPointT &t02,
	                const GPointT &t10, const GPointT &t11,
	                const GPointT &t12) const;

	/**
	 * @brief Get the intersection type between two triangles.
	 * @note Assume that no triangle is degenerate.
	 */
	SimplexIntersectionType intersection_type(const TriangleT &tri0,
	                                          const TriangleT &tri1) const;

	/**
	 * @brief Get the intersection type between two triangles.
	 * @param t00_t01_t02 three vertices of the first triangle.
	 * @param t10_t11_t12 three vertices of the second triangle.
	 * @note Assume that no triangle is degenerate.
	 */
	SimplexIntersectionType
	intersection_type(const GPointT &t00, const GPointT &t01, const GPointT &t02,
	                  const GPointT &t10, const GPointT &t11,
	                  const GPointT &t12) const;

	/**
	 * @brief Get the intersection type between two triangles.
	 * @param t00_t01_t02 three vertices of the first triangle.
	 * @param t10_t11_t12 three vertices of the second triangle.
	 * @param t0_minor_perm minor and permutation for triangle t0.
	 * @param t1_minor_perm minor and permutation for triangle t1.
	 * @note Assume that no triangle is degenerate.
	 */
	SimplexIntersectionType
	intersection_type(const NT *t00, const NT *t01, const NT *t02, const NT *t10,
	                  const NT *t11, const NT *t12, const NT *t0_min = nullptr,
	                  const NT *t0_perm = nullptr, const NT *t1_min = nullptr,
	                  const NT *t1_perm = nullptr) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Triangle3_Triangle3.inl"
#endif