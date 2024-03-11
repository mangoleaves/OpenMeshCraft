#pragma once

#include "CheckDegenerate3K.h"

namespace OMC {

/**
 * @brief Check if a segment is degenerate. If it is degenerate, return the
 * degeneration type and geometry.
 * @return std::variant Void if no degeneration, otherwise it contains
 * primitives whose types are defined in this class \ref CheckDegenerate3K.
 * @retval PointT If the segment degenerate to point, return the point.
 */
template <typename Kernel>
auto CheckDegenerate3K<Kernel>::operator()(const SegmentT &segment) -> DgnType
{
	if (LessThan3D().coincident(segment.start(), segment.end()))
		return DgnType{segment.start()}; // degenerate to point
	else
		return DgnType(); // contains NoDgn.
}

/**
 * @brief Check if a triangle is degenerate. If it is degenerate, return the
 * degeneration type and geometry.
 * @return std::variant Void if no degeneration, otherwise it contains
 * primitives whose types are defined in this class \ref CheckDegenerate3K.
 * @retval PointT If the triangle degenerates to point.
 * @retval SegmentT If the trianngle degenerates to segment.
 */
template <typename Kernel>
auto CheckDegenerate3K<Kernel>::operator()(const TriangleT &triangle) -> DgnType
{
	auto &t = triangle;

	if (LessThan3D().coincident(t.v0(), t.v1()) &&
	    LessThan3D().coincident(t.v1(), t.v2()))
	{
		return DgnType{t.v0()}; // degenerate to point
	}
	else if (CollinearPoints3D()(t.v0(), t.v1(), t.v2()))
	{
		auto [a, b, c] = CollinearSort3D()(t.v0(), t.v1(), t.v2());
		return DgnType{SegmentT(a, c)}; // degenerate to segment
	}
	else
		return DgnType(); // NoDgn
}

} // namespace OMC