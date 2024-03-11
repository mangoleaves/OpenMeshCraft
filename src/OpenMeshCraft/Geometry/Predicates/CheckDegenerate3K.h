#pragma once

#include <variant>

namespace OMC {

/**
 * @brief CheckDegenerate3K implements several algorithms that check if a
 * primitive in 3D is degenerate.
 * @tparam Kernel
 */
template <typename Kernel>
class CheckDegenerate3K
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using GPointT   = typename K::GPoint3;
	using SegmentT  = typename K::Segment3;
	using TriangleT = typename K::Triangle3;

	using LessThan3D       = typename K::LessThan3D;
	using CollinearPoints3D = typename K::CollinearPoints3D;
	using CollinearSort3D   = typename K::CollinearSort3D;

	struct NoDgn
	{
	};
	using DgnType = std::variant<NoDgn, GPointT, SegmentT>;

public:
	/**
	 * @brief Check if a segment is degenerate. If it is degenerate, return the
	 * degeneration type and geometry.
	 * @return std::variant Void if no degeneration, otherwise it contains
	 * primitives whose types are defined in this class \ref CheckDegenerate3K.
	 * @retval PointT If the segment degenerate to point, return the point.
	 */
	DgnType operator()(const SegmentT &segment);

	/**
	 * @brief Check if a triangle is degenerate. If it is degenerate, return the
	 * degeneration type and geometry.
	 * @return std::variant Void if no degeneration, otherwise it contains
	 * primitives whose types are defined in this class \ref CheckDegenerate3K.
	 * @retval PointT If the triangle degenerates to point.
	 * @retval SegmentT If the trianngle degenerates to segment.
	 */
	DgnType operator()(const TriangleT &triangle);
};

} // namespace OMC

#include "CheckDegenerate3K.inl"