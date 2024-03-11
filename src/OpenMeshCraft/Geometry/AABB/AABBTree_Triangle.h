#pragma once

#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

namespace OMC {

/***********************************************/
/* 1. Define the reference point of the triangle
/***********************************************/

/**
 * @brief The type of reference point.
 *
 */
enum class AABB_Triangle_ReferencePointType
{
	First,   /// The first point of the primitive, e.g., v0 of a triangle.
	Centroid /// The centroid of the primitive.
};

/**
 * @brief Get the reference point of triangle.
 * The reference point will be used to split triangles along an axis.
 * @tparam TriT
 * @tparam RefPntType
 */
template <typename TriT, AABB_Triangle_ReferencePointType RefPntType>
class AABB_Triangle_ReferencePoint
{
	static_assert(
	  RefPntType == AABB_Triangle_ReferencePointType::First ||
	    RefPntType == AABB_Triangle_ReferencePointType::Centroid,
	  "AABB_Triangle_ReferencePoint havn't support other reference point type.");

public:
	using NT = typename TriT::NT;

public:
	auto operator()(const TriT &tri)
	  -> remove_cvref_t<decltype(std::declval<TriT>().v0())>
	{
		if constexpr (RefPntType == AABB_Triangle_ReferencePointType::First)
		{
			return tri.v0();
		}
		else if constexpr (RefPntType == AABB_Triangle_ReferencePointType::Centroid)
		{
			return (tri.v0() + tri.v1() + tri.v2()) / NT(3.);
		}
	}
};

} // namespace OMC