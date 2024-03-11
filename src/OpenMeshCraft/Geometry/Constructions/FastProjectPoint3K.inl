#pragma once

#include "FastProjectPoint3K.h"

#include <type_traits>

namespace OMC {

template <typename Kernel>
FastProjectPoint3K<Kernel>::AuxTriangle::AuxTriangle(const EPointT &p1,
                                                     const EPointT &p2,
                                                     const EPointT &p3)
{
	update_aux_triangle(p1, p2, p3);
}

template <typename Kernel>
void FastProjectPoint3K<Kernel>::AuxTriangle::update_aux_triangle(
  const EPointT &p1, const EPointT &p2, const EPointT &p3)
{
	VecT p12 = p2 - p1, p23 = p3 - p2, p31 = p1 - p3;
	// find obtuse angle, move obtuse angle point to ap_v0
	if (p12.dot(p31) > 0)
	{
		ap_has_obtuse_angle = true;
		ap_v0 = p1, ap_v1 = p2, ap_v2 = p3;
	}
	else if (p12.dot(p23) > 0)
	{
		ap_has_obtuse_angle = true;
		ap_v0 = p2, ap_v1 = p3, ap_v2 = p1;
	}
	else if (p23.dot(p31) > 0)
	{
		ap_has_obtuse_angle = true;
		ap_v0 = p3, ap_v1 = p1, ap_v2 = p2;
	}
	else
	{
		ap_has_obtuse_angle = false;
		ap_v0 = p1, ap_v1 = p2, ap_v2 = p3;
	}
	// calculate edge vector and squared length
	ap_edge_vec01    = ap_v1 - ap_v0;
	ap_edge_vec02    = ap_v2 - ap_v0;
	ap_edge_vec12    = ap_v2 - ap_v1;
	ap_edge_sqrlen01 = ap_edge_vec01.sqrnorm();
	ap_edge_sqrlen02 = ap_edge_vec02.sqrnorm();
	ap_edge_sqrlen12 = ap_edge_vec12.sqrnorm();

	// calculate face normal
	ap_face_normal         = ap_edge_vec01.cross(ap_edge_vec02);
	ap_is_plane_degenerate = std::abs(ap_face_normal.x()) < 1e-10 &&
	                         std::abs(ap_face_normal.y()) < 1e-10 &&
	                         std::abs(ap_face_normal.z()) < 1e-10;
	ap_face_normal.normalize();

	// calculate edge normal
	if (!ap_is_plane_degenerate)
	{
		ap_edge_normal01 = (ap_edge_vec01).cross(ap_face_normal).normalized();
		ap_edge_normal02 = (-ap_edge_vec02).cross(ap_face_normal).normalized();
		ap_edge_normal12 = (ap_edge_vec12).cross(ap_face_normal).normalized();
	}
}

template <typename Kernel>
template <typename PrimT>
auto FastProjectPoint3K<Kernel>::operator()(const PrimT   &prim,
                                            const GPointT &query) const
  -> EPointT
{
	if constexpr (std::is_base_of_v<AuxTriangle, PrimT>)
	{
		return project_to_aux_triangle(prim, query);
	}
	// else if constexpr (std::is_base_of_v<OtherAuxType, PrimT>)
	// {
	// 		project to other aux type.
	// }
	else
	{
		return EPointT();
	}
}

template <typename Kernel>
template <typename TriT>
auto FastProjectPoint3K<Kernel>::project_to_aux_triangle(
  const TriT &triangle, const GPointT &query) const -> EPointT
{
	EPointT     q = ToEP()(query);
	const TriT &t = triangle;
	// If the plane is degenerate, then the triangle is degenerate, and
	// one tries to find to which segment it is equivalent.
	// Otherwise project to the non-degenerate triangle.
	if (triangle.ap_is_plane_degenerate)
	{
		SegmentT longest_segment;
		if (t.ap_edge_sqrlen01 > t.ap_edge_sqrlen12)
		{
			if (t.ap_edge_sqrlen01 > t.ap_edge_sqrlen02)
				longest_segment.start() = t.ap_v0, longest_segment.end() = t.ap_v1;
			else
				longest_segment.start() = t.ap_v2, longest_segment.end() = t.ap_v0;
		}
		else
		{
			if (t.ap_edge_sqrlen12 > t.ap_edge_sqrlen02)
				longest_segment.start() = t.ap_v1, longest_segment.end() = t.ap_v2;
			else
				longest_segment.start() = t.ap_v2, longest_segment.end() = t.ap_v0;
		}
		ProjectPoint3 project;
		return project(longest_segment, query);
	}
	else
	{
		if (!t.ap_has_obtuse_angle)
		{
			// if there is no obtuse angle, when the proj_point is outside the
			// triangle, the proj_point can be simplify projected to edge segment.
			VecT v0_to_query = q - t.ap_v0;
			if (v0_to_query.dot(t.ap_edge_normal01) > 0.0)
			{
				// if it is outside edge 01, project it to edge 01
				NT numerator = v0_to_query.dot(t.ap_edge_vec01);
				if (numerator < 0.0)
					return t.ap_v0;
				else if (numerator > t.ap_edge_sqrlen01)
					return t.ap_v1;
				else
					return t.ap_v0 + (numerator / t.ap_edge_sqrlen01) * t.ap_edge_vec01;
			}
			else if (v0_to_query.dot(t.ap_edge_normal02) > 0.0)
			{
				// if it is outside edge 02, project it to edge 02
				NT numerator = v0_to_query.dot(t.ap_edge_vec02);
				if (numerator < 0.0)
					return t.ap_v0;
				else if (numerator > t.ap_edge_sqrlen02)
					return t.ap_v2;
				else
					return t.ap_v0 + (numerator / t.ap_edge_sqrlen02) * t.ap_edge_vec02;
			}
		}
		else
		{
			// if there is obtuse angle, when the proj_point is outside the triangle,
			// the proj_point can **NOT** be simplify projected to edge segment.
			VecT v0_to_query = q - t.ap_v0;
			if (v0_to_query.dot(t.ap_edge_normal01) > 0.0)
			{
				// if it is outside edge 01, try to project it to edge 01
				NT numerator01 = v0_to_query.dot(t.ap_edge_vec01);
				if (numerator01 > t.ap_edge_sqrlen01)
					return t.ap_v1;
				else if (numerator01 < 0.0)
				{
					// must also outside edge 02
					NT numerator02 = v0_to_query.dot(t.ap_edge_vec02);
					if (numerator02 > t.ap_edge_sqrlen02)
						return t.ap_v2;
					else if (numerator02 < 0.0)
						return t.ap_v0;
					else
						return t.ap_v0 +
						       (numerator02 / t.ap_edge_sqrlen02) * t.ap_edge_vec02;
				}
				else
					return t.ap_v0 + (numerator01 / t.ap_edge_sqrlen01) * t.ap_edge_vec01;
			}
			else if (v0_to_query.dot(t.ap_edge_normal02) > 0.0)
			{
				// if it is outside edge 02, project it to edge 02
				NT numerator = v0_to_query.dot(t.ap_edge_vec02);
				if (numerator > t.ap_edge_sqrlen02)
					return t.ap_v2;
				else
					return t.ap_v0 + (numerator / t.ap_edge_sqrlen02) * t.ap_edge_vec02;
			}
		}
		VecT v1_to_query = q - t.ap_v1;
		if (v1_to_query.dot(t.ap_edge_normal12) > 0.0) // is outside edge 12 ?
		{
			NT numerator = v1_to_query.dot(t.ap_edge_vec12);
			if (numerator < 0.0)
				return t.ap_v1;
			else if (numerator > t.ap_edge_sqrlen12)
				return t.ap_v2;
			else
				return t.ap_v1 + (numerator / t.ap_edge_sqrlen12) * t.ap_edge_vec12;
		}
		NT propotion = (v1_to_query).dot(t.ap_face_normal); // is inside triangle!
		return q - propotion * t.ap_face_normal;
	}
}

} // namespace OMC