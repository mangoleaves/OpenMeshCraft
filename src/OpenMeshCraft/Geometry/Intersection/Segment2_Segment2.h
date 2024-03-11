#pragma once

#include "IntersectionUtils.h"

namespace OMC {

/**
 * @brief Check if Triangle2 and Point2 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Segment2_Segment2_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT     = typename K::Vec2;
	using GPointT  = typename K::GPoint2;
	using SegmentT = typename K::Segment2;

	using LessThan2D = typename K::LessThan2D;
	using Orient2D   = typename K::Orient2D;

public:
	/**
	 * @brief Check if segment and segment intersect.
	 * @note assume thant no segment is degenerate.
	 */
	bool operator()(const SegmentT &seg0, const SegmentT &seg1) const;

	/**
	 * @brief Get the intersection type between segment and segment.
	 * @note assume thant no segment is degenerate.
	 */
	SimplexIntersectionType intersection_type(const SegmentT &seg0,
	                                          const SegmentT &seg1) const;

	/**
	 * @brief Get the intersection type between segment and segment.
	 * @param s00_s01 two vertices of the segment
	 * @param s10_s11 two vertices of the segment
	 * @note assume thant no segment is degenerate.
	 */
	SimplexIntersectionType intersection_type(const GPointT &s00,
	                                          const GPointT &s01,
	                                          const GPointT &s10,
	                                          const GPointT &s11) const;

	/**
	 * @brief Get the intersection type between segment and segment.
	 * @param s00_s01 two vertices of the segment
	 * @param s10_s11 two vertices of the segment
	 * @note assume thant no segment is degenerate.
	 */
	SimplexIntersectionType intersection_type(const NT *s00, const NT *s01,
	                                          const NT *s10, const NT *s11) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Segment2_Segment2.inl"
#endif