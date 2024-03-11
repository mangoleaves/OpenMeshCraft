#pragma once

#include "IntersectionUtils.h"

namespace OMC {

/**
 * @brief Check if Segment3 and Point3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Segment3_Point3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT     = typename K::Vec3;
	using GPointT  = typename K::GPoint3;
	using SegmentT = typename K::Segment3;

	using OrientOn2D       = typename K::OrientOn2D;
	using Orient3D         = typename K::Orient3D;
	using CollinearPoints3D = typename K::CollinearPoints3D;
	using LessThan3D       = typename K::LessThan3D;

public:
	/**
	 * @brief Check if segment and point intersect.
	 */
	bool operator()(const SegmentT &seg, const GPointT &pnt) const;

	/**
	 * @brief Check if segment and point intersect.
	 */
	bool operator()(const GPointT &s0, const GPointT &s1, const GPointT &pnt) const;

	/**
	 * @brief Get the relative position of point to segment.
	 * @return relative position of point to segment.
	 */
	PointInType in_segment(const SegmentT &seg, const GPointT &pnt) const;

	/**
	 * @brief Get the relative position of point to segment.
	 * @param s0_s1 two vertices of the segment.
	 * @param pnt point.
	 * @return relative position of point to segment.
	 */
	PointInType in_segment(const GPointT &s0, const GPointT &s1,
	                       const GPointT &pnt) const;

	/**
	 * @brief Get the relative position of point to segment.
	 * @return relative position of point to segment.
	 * @note this function assume that segment and point are collinear.
	 */
	PointInType in_segment_collinear(const SegmentT &seg,
	                                const GPointT  &pnt) const;

	/**
	 * @brief Get the relative position of point to segment.
	 * @param s0_s1 two vertices of the segment.
	 * @param pnt point.
	 * @return relative position of point to segment.
	 * @note this function assume that three points are collinear.
	 */
	PointInType in_segment_collinear(const GPointT &s0, const GPointT &s1,
	                                const GPointT &pnt) const;

	/**
	 * @brief Get the intersection type between segment and point.
	 */
	PointInSimplexType intersection_type(const SegmentT &seg,
	                                     const GPointT  &pnt) const;

	/**
	 * @brief Get the intersection type between segment and point.
	 * @param s0_s1 two vertices of segment
	 * @param p point
	 */
	PointInSimplexType intersection_type(const GPointT &s0, const GPointT &s1,
	                                     const GPointT &p) const;

	/**
	 * @brief Get the intersection type between segment and point.
	 * @param s0_s1 two vertices of segment
	 * @param p point
	 */
	PointInSimplexType intersection_type(const NT *s0, const NT *s1,
	                                     const NT *p) const;

	// Equivalent but faster. They assume that the
	// dominant normal component is n_max (see MaxComponentInTriangleNormal).

	bool operator()(const SegmentT &seg, const GPointT &pnt, int n_max) const;

	PointInType in_segment(const SegmentT &seg, const GPointT &pnt,
	                       int n_max) const;

	PointInType in_segment(const GPointT &s0, const GPointT &s1,
	                       const GPointT &pnt, int n_max) const;

	PointInSimplexType intersection_type(const GPointT &s0, const GPointT &s1,
	                                     const GPointT &p, int n_max) const;

	PointInSimplexType intersection_type(const NT *s0, const NT *s1, const NT *p,
	                                     int n_max) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Segment3_Point3.inl"
#endif