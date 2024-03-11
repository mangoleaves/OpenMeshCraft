#pragma once

#include "AABBTree_Triangle.h"

#include "AABBTraits.h"
#include "AABBTraversalTraits.h"
#include "AABBTree.h"

#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/Primitives/PrimitiveWithAttribute.h"

namespace OMC {

/***********************************************/
/* 2. Define minimum set of AABB traits
/***********************************************/
template <typename Kernel>
class AABBMinimumTraits_Triangle_Intersection
{
public:
	/* Belows are used by AABBTree */

	// Triangle type
	// using TriT = typename APAC::Triangle3;
	// Attribute type
	// NOTE: If you don't want attach attribute, comment it and replace PrimT with
	// PrimitiveType.
	using K    = Kernel;
	using TriT = typename K::Triangle3;

	using PrimAttrT = size_t;

	// Primitive type
	using PrimT = PrimitiveWithAttribute<TriT, PrimAttrT>;
	// Primitive reference point
	using PrimReferencePoint =
	  AABB_Triangle_ReferencePoint<PrimT,
	                               AABB_Triangle_ReferencePointType::First>;
	// Calculate bounding box
	using CalcBbox = typename K::CalcBoundingBox3;
};

/****************************************************/
/* 3. Automatically deduce the complete AABB traits
/****************************************************/
template <typename Kernel>
using AABBTraits_Triangle_Intersection =
  AABBAutoDeduceTraits<AABBMinimumTraits_Triangle_Intersection<Kernel>>;

/****************************************************/
/* 4. Define the AABB Tree used for intersection.
/****************************************************/

template <typename Kernel>
class AABBTree_Triangle_Intersection
  : public AABBTree<AABBTraits_Triangle_Intersection<Kernel>>
{
public:
	using K = Kernel;

	using BaseT  = AABBTree<AABBTraits_Triangle_Intersection<K>>;
	using ThisT  = AABBTree_Triangle_Intersection<K>;
	using Traits = AABBTraits_Triangle_Intersection<K>;

	using TriT      = typename Traits::TriT;
	using PrimT     = typename Traits::PrimT;
	using PrimAttrT = typename Traits::PrimAttrT;
	// using Point_Projection = typename Traits::Point_Projection;
	//  Other types will be deduced by AABBTree

	using Indices = std::vector<size_t>;

public:
	template <typename QPrimT>
	void all_intersections(const QPrimT &query, Indices &results) const;

	const typename BaseT::BboxT &get_Bbox(const size_t i) const
	{
		return this->m_nodes[i].bbox();
	}

protected:
	template <typename _QPrimT>
	class BoxInterTraits
	{
	public:
		using QPrimT = _QPrimT;

		using NT     = typename Traits::NT;
		using PointT = typename Traits::PointT;
		using PrimT  = typename Traits::PrimT;
		using BboxT  = typename Traits::BboxT;

		using DoIntersect = typename K::DoIntersect;
		using CalcBbox    = typename Traits::CalcBbox;
	};

	template <typename QPrimT>
	using BoxTrav = AABB_BoxInterTraversal<BoxInterTraits<QPrimT>>;
};

template <typename Kernel>
template <typename QPrimT>
inline void AABBTree_Triangle_Intersection<Kernel>::all_intersections(
  const QPrimT &query, Indices &results) const
{
	BoxTrav<QPrimT> box_trav(query);
	this->traversal(box_trav);

	auto &prim_ptrs = box_trav.result();
	results.reserve(prim_ptrs.size());
	for (auto pp : prim_ptrs)
		results.push_back(pp->attribute());
}

} // namespace OMC
