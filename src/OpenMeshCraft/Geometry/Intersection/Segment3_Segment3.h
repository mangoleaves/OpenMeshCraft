#pragma once

#include "IntersectionUtils.h"
#include "Segment2_Segment2.h"

namespace OMC {

/**
 * @brief Check if two Segment3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Segment3_Segment3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using VecT     = typename K::Vec3;
	using GPointT  = typename K::GPoint3;
	using SegmentT = typename K::Segment3;

	using LessThan3D = typename K::LessThan3D;
	using OrientOn2D = typename K::OrientOn2D;
	using Orient3D   = typename K::Orient3D;

	using Segment2_Segment2_DoInter = Segment2_Segment2_Do_Intersect<Kernel>;

public:
	/**
	 * @brief Check if two segments intersect. Two segments
	 * are assumed to be coplanar, otherwise the result is undetermined.
	 * @note Assume that no segment is degenerate.
	 */
	bool operator()(const SegmentT &seg0, const SegmentT &seg1) const;

	/**
	 * @brief Check if two segments intersect. Two segments
	 * are assumed to be coplanar, otherwise the result is undetermined.
	 * @param AB Point A and B form the first segment
	 * @param PQ Point P and Q form the second segment
	 * @note Assume that no segment is degenerate.
	 */
	bool operator()(const GPointT &A, const GPointT &B, const GPointT &P,
	                const GPointT &Q) const;

	/*************************************************************************/
	/******* Below functions all assume that two segments are coplanar *******/
	/*************************************************************************/

	/**
	 * @brief Check if the closure of two segments intersect at a single point
	 * (overlap is excluded). Two segments are assumed to be coplanar, otherwise
	 * the result is undetermined.
	 * @return TRUE if the closure of two segments intersect at a single point.
	 * @note Assume that no segment is degenerate and two segments are coplanar.
	 */
	bool cross(const SegmentT &seg0, const SegmentT &seg1) const;

	/**
	 * @brief Check if the interior of two segments in 3D intersect at a single
	 * point (overlap is excluded). Two segments are assumed to be coplanar,
	 * otherwise the result is undetermined.
	 * @return TRUE if the interior of two segments intersect at a single point.
	 * @note Assume that no segment is degenerate and two segments are coplanar.
	 */
	bool cross_inner(const SegmentT &seg0, const SegmentT &seg1) const;

	/**
	 * @brief Check if the closure of two segments intersect at a single point
	 * (overlap is excluded). Two segments are assumed to be coplanar, otherwise
	 * the result is undetermined.
	 * @param AB Point A and B form the first segment
	 * @param PQ Point P and Q form the second segment
	 * @return TRUE if the closure of two segments intersect at a single point.
	 * @note Assume that no segment is degenerate and two segments are coplanar.
	 */
	bool cross(const GPointT &A, const GPointT &B, const GPointT &P,
	           const GPointT &Q) const;

	/**
	 * @brief Check if the interior of two segments in 3D intersect at a single
	 * point (overlap is excluded). Two segments are assumed to be coplanar,
	 * otherwise the result is undetermined.
	 * @param AB Point A and B form the first segment
	 * @param PQ Point P and Q form the second segment
	 * @return TRUE if the interior of two segments intersect at a single point.
	 * @note Assume that no segment is degenerate and two segments are coplanar.
	 */
	bool cross_inner(const GPointT &A, const GPointT &B, const GPointT &P,
	                 const GPointT &Q) const;

	// Equivalent but faster. They assume that points are coplanar and the
	// dominant normal component is n_max (see MaxComponentInTriangleNormal).

	bool cross(const GPointT &A, const GPointT &B, const GPointT &P,
	           const GPointT &Q, int n_max) const;

	bool cross_inner(const GPointT &A, const GPointT &B, const GPointT &P,
	                 const GPointT &Q, int n_max) const;

	/**
	 * @brief Get the intersection type between two segments.
	 * @note Assume that two segments are coplanar.
	 */
	SimplexIntersectionType intersection_type(const SegmentT &seg0,
	                                          const SegmentT &seg1) const;

	/**
	 * @brief Get the intersection type between two segments.
	 * @param AB two points form the first segment.
	 * @param PQ two points form the second segment.
	 * @note Assume that two segments are coplanar.
	 */
	SimplexIntersectionType intersection_type(const GPointT &A, const GPointT &B,
	                                          const GPointT &P,
	                                          const GPointT &Q) const;

	/**
	 * @brief Get the intersection type between two segments.
	 * @param AB two points form the first segment.
	 * @param PQ two points form the second segment.
	 * @note Assume that two segments are coplanar.
	 */
	SimplexIntersectionType intersection_type(const NT *A, const NT *B,
	                                          const NT *P, const NT *Q) const;

	// Equivalent but faster. They assume that the
	// dominant normal component is n_max (see MaxComponentInTriangleNormal).

	/**
	 * @brief Get the intersection type between two segments.
	 * @param AB two points form the first segment.
	 * @param PQ two points form the second segment.
	 * @param n1_max max normal component of segment PQ.
	 * @note Assume that two segments are coplanar.
	 */
	SimplexIntersectionType intersection_type(const NT *A, const NT *B,
	                                          const NT *P, const NT *Q,
	                                          int n1_max) const;

protected:
	/**
	 * @brief Get the intersection type between two segments on 2D plane.
	 * @param AB two points form the first segment.
	 * @param PQ two points form the second segment.
	 * @param on2d 0: yz plane; 1: zx plane; 2: xy plane.
	 * @note almostly same as Segment2_Segment2_Do_Intersect.
	 */
	SimplexIntersectionType
	intersection_type_on2d(const GPointT &A, const GPointT &B, const GPointT &P,
	                       const GPointT &Q, int on2d) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Segment3_Segment3.inl"
#endif