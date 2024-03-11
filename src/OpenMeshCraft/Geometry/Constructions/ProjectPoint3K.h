#pragma once

#include "OpenMeshCraft/Utils/Exception.h"

#include <variant>

namespace OMC {

/**
 * @brief ProjectPoint3K implement several project algorithms that project a
 * point to different primitives, such as segment and triangle.
 * @tparam Kernel.
 */
template <typename Kernel>
class ProjectPoint3K
{
public:
	using K = Kernel;

	using NT           = typename K::NT;
	using VecT         = typename K::Vec3;
	using EPointT      = typename K::EPoint3;
	using GPointT      = typename K::GPoint3;
	using BoundingBoxT = typename K::BoundingBox3;
	using SegmentT     = typename K::Segment3;
	using TriangleT    = typename K::Triangle3;

	using ToEP = typename K::ToEP;

	using Orient3D   = typename K::Orient3D;
	using LessThan3D = typename K::LessThan3D;

	using CheckDegenerate3 = typename K::CheckDegenerate3;

	using ConstructNormal3 = typename K::ConstructNormal3;

public:
	/**
	 * @brief Project a point to segment.
	 * @return the projected point.
	 */
	EPointT operator()(const SegmentT &segment, const GPointT &point) const;

	/**
	 * @brief Project a point to triangle.
	 * @return the projected point.
	 */
	EPointT operator()(const TriangleT &triangle, const GPointT &point) const;

	/**
	 * @brief Project a point to box.
	 * @return the projected point.
	 */
	EPointT operator()(const BoundingBoxT &bbox, const GPointT &point) const;

private:
	/**
	 * @brief When segment is not degenerate, project point to it.
	 * @return GPointT the projected point.
	 */
	EPointT proj_to_segment(const SegmentT &segment, const GPointT &point) const;

	/**
	 * @brief When triangle is not degenerate, project point to it.
	 * @return GPointT the projected point.
	 */
	EPointT proj_to_triangle(const TriangleT &triangle,
	                         const GPointT   &point) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ProjectPoint3K.inl"
#endif