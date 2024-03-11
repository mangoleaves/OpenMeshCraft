#pragma once

#include "ProjectPoint3K.h"

#include <variant>

namespace OMC {

template <typename Kernel>
auto ProjectPoint3K<Kernel>::operator()(const SegmentT &segment,
                                        const GPointT  &point) const -> EPointT
{
	typename CheckDegenerate3::DgnType degeneration = CheckDegenerate3()(segment);

	if (std::holds_alternative<typename CheckDegenerate3::NoDgn>(degeneration))
	{
		return proj_to_segment(segment, point);
	}
	else if (std::holds_alternative<GPointT>(degeneration))
	{
		// segment is degenerate to point, point is projected to point.
		return std::get<GPointT>(std::move(degeneration));
	}

#if OMC_ENABLE_ASSERT
	OMC_ASSERT(false, "Unexpected degenenration type");
#else
	return EPointT();
#endif
}

template <typename Kernel>
auto ProjectPoint3K<Kernel>::operator()(const TriangleT &triangle,
                                        const GPointT   &point) const -> EPointT
{

	typename CheckDegenerate3::DgnType degeneration =
	  CheckDegenerate3()(triangle);

	if (std::holds_alternative<typename CheckDegenerate3::NoDgn>(degeneration))
	{
		return proj_to_triangle(triangle, point);
	}
	else if (std::holds_alternative<SegmentT>(degeneration))
	{
		// segment is degenerate to point, point is projected to point.
		return proj_to_segment(std::get<SegmentT>(degeneration), point);
	}
	else if (std::holds_alternative<GPointT>(degeneration))
	{
		// segment is degenerate to point, point is projected to point.
		return std::get<GPointT>(std::move(degeneration));
	}
#if OMC_ENABLE_ASSERT
	OMC_ASSERT(false, "Unexpected degenenration type");
#else
	return EPointT();
#endif
}

template <typename Kernel>
auto ProjectPoint3K<Kernel>::operator()(const BoundingBoxT &bbox,
                                        const GPointT &point) const -> EPointT
{
	EPointT result = ToEP()(point);
#define COMPARE_AXIS(axis)                                                \
	if (LessThan3D().on_##axis(point, bbox.min_bound()) <= Sign::ZERO)      \
		result.axis() = bbox.min_bound().axis();                              \
	else if (LessThan3D().on_##axis(point, bbox.max_bound()) >= Sign::ZERO) \
		result.axis() = bbox.max_bound().axis();

	COMPARE_AXIS(x);
	COMPARE_AXIS(y);
	COMPARE_AXIS(z);
#undef COMPARE_AXIS
	return result;
}

template <typename Kernel>
auto ProjectPoint3K<Kernel>::proj_to_segment(const SegmentT &segment,
                                             const GPointT  &point) const
  -> EPointT
{
	// No degeneration, project point to the segment.
	VecT segment_vec = segment.end() - segment.start();
	NT   numerator   = segment_vec.dot(point - segment.start());

	if (numerator /*segment_to_vector.dot(query - *segment.first)*/ <= NT(0.0))
	{
		return segment.start();
	}
	else if (segment_vec.dot(point - segment.end()) >= NT(0.0))
	{
		return segment.end();
	}
	else
	{
		NT dominator = segment_vec.sqrnorm();
		return segment.start() + (numerator / dominator) * segment_vec;
	}
}

template <typename Kernel>
auto ProjectPoint3K<Kernel>::proj_to_triangle(const TriangleT &triangle,
                                              const GPointT   &point) const
  -> EPointT
{
	const EPointT &a = triangle.v0(), &b = triangle.v1(), &c = triangle.v2();
	EPointT        ep = ToEP()(point);

	// 1. Project point to the support plane of the triangle.
	VecT    normal      = ConstructNormal3()(triangle);
	NT      numerator   = (ep - a).dot(normal);
	NT      denominator = normal.sqrnorm();
	EPointT proj_point  = ep - (numerator / denominator) * normal;
	// proj_point is the project point on the support plane.

	// 2. Check if the projected point is inside the triangle.
	// * If it is inside the triangle, then return it as result.
	// * If it is outside the triangle, further project proj_point to one segment
	// depending on its position.

	// If the proj_point is outside the triangle, we project it to triangle edge
	// lines and check if it is on the triangle edge segments. We need to try no
	// more than three edges.

	bool    outside = false;
	EPointT result;

	auto check_seg = [&ep, &normal, &proj_point, &outside,
	                  &result](const EPointT &start,
	                           const EPointT &end) -> bool // found result
	{
		// Offset points of triangle by the normal of triangle to construct a plane
		// which is orthogonal to the triangle.
		EPointT offset_start = start + normal;

		// The orientation of point with respect to orthogonal plane defined on ij,
		// where ij means vertex vi and vj.
		// These orientations are used to determine whether the proj_point is
		// inside/outside the triangle.
		Sign ori_point = Orient3D()(start, end, offset_start, ep);
		if (ori_point != Sign::NEGATIVE)
		{
			outside = true;

			// 3. The proj_point is outside the triangle.
			// Project proj_point to Line, check if it is on the Segment.
			VecT seg_vec        = end - start;
			VecT start_to_point = proj_point - start;
			VecT end_to_point   = proj_point - end;

			NT ori_start_to_end = start_to_point.dot(seg_vec);
			NT ori_end_to_start = end_to_point.dot(seg_vec);
			if (ori_start_to_end >= NT(0.0) && ori_end_to_start <= NT(0.0))
			{
				// 4. The proj_proj_point is on the Segment
				NT seg_sqrnorm = seg_vec.sqrnorm();
				result         = start + (ori_start_to_end / seg_sqrnorm) * seg_vec;
				return true;
			}
			// 5. The proj_proj_point is not on the Segment.
			// It is possibly on another segment or vertices.
		}
		return false;
	};

	if (check_seg(a, b) || check_seg(b, c) || check_seg(c, a))
		return result;

	if (outside)
	{
		NT dis_to_a = (ep - a).sqrnorm();
		NT dis_to_b = (ep - b).sqrnorm();
		NT dis_to_c = (ep - c).sqrnorm();
		if (dis_to_a <= dis_to_b)
		{
			if (dis_to_a <= dis_to_c)
				return a;
			else
				return c;
		}
		else
		{
			if (dis_to_b <= dis_to_c)
				return b;
			else
				return c;
		}
	}
	else // proj_point is inside the triangle, return it.
		return proj_point;
}

} // namespace OMC