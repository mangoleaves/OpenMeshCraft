#pragma once

#include "OpenMeshCraft/Geometry/Primitives/Triangle3T.h"
#include "OpenMeshCraft/Utils/Exception.h"

#include <type_traits>

namespace OMC {

/**
 * @brief The triangle with additional geometry information for fast projection.
 * @tparam Kernel
 */
template <typename Kernel>
class FastProjectPoint3K
{
	static_assert(std::is_floating_point_v<typename Kernel::NT>,
	              "Only support floating point number.");

public:
	using K = Kernel;

	using NT           = typename K::NT;
	using VecT         = typename K::Vec3;
	using GPointT      = typename K::GPoint3;
	using EPointT      = typename K::EPoint3;
	using BoundingBoxT = typename K::BoundingBox3;
	using SegmentT     = typename K::Segment3;
	using TriangleT    = typename K::Triangle3;

	using ToEP = typename K::ToEP;

	using ProjectPoint3 = typename K::ProjectPoint3;

	/* Auxiliary data structures for fast projection */

	/**
	 * @brief An auxiliary data structure of base triangle for fast projecting a
	 * point to a triangle.
	 * This class can be inherited by a user-defined triangle.
	 * For example, a user may inherit the base triangle, this aux tri and other
	 * aux tris to get a triangle they want:
	 * `class Tri : public Triangle3T, AuxTriangle0, AuxTriangle1...`
	 */
	class AuxTriangle
	{
	public:
		/// @brief Initialize the aux triangle with three points from base triangle.
		AuxTriangle(const EPointT &p1, const EPointT &p2, const EPointT &p3);

		/// @brief Update the aux triangle with three points from base triangle.
		void update_aux_triangle(const EPointT &p1, const EPointT &p2,
		                         const EPointT &p3);

	private:
		// used for accelerating projecting a point to this triangle.
		EPointT ap_v0, ap_v1, ap_v2;
		VecT    ap_face_normal;
		VecT    ap_edge_vec01, ap_edge_vec02, ap_edge_vec12;
		VecT    ap_edge_normal01, ap_edge_normal02, ap_edge_normal12;
		NT      ap_edge_sqrlen01, ap_edge_sqrlen02, ap_edge_sqrlen12;
		bool    ap_has_obtuse_angle;
		bool    ap_is_plane_degenerate;

		friend class FastProjectPoint3K<Kernel>;
	};

public:
	template <typename PrimT>
	EPointT operator()(const PrimT &prim, const GPointT &point) const;

private:
	/* methods for AuxTriangle */
	template <typename TriT>
	EPointT project_to_aux_triangle(const TriT    &triangle,
	                                const GPointT &query) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "FastProjectPoint3K.inl"
#endif