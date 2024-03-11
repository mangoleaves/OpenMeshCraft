#pragma once

#include "OpenMeshCraft/Geometry/Primitives/PrimitiveWithAttribute.h"

#include "OpenMeshCraft/Utils/SFINAE.h"

#include <vector>

namespace OMC {



template <typename KdTraits>
class KdAutoDeduceTraits
{
public:
	/* Below types must be provided */

	// point type
	using PointT       = typename KdTraits::PointT;
	// bounding box type
	using BboxT        = typename KdTraits::BboxT;
	// project point to bounding box to calculate the distance between them.
	using ProjectPoint = typename KdTraits::ProjectPoint;

	/* Try to deduce below types*/
	// number type
	GET_TYPE_OTHERWISE_VOID(PointT, NT, DeducedNT);
	GET_TYPE_OTHERWISE_DEFAULT(KdTraits, NT, DeducedNT, NT);
	static_assert(!std::is_void_v<NT>,
	              "NT is not provided and can't be deduced.");
	// point attribute type. if it is not provided, we use an arbitrary type.
	GET_TYPE_OTHERWISE_DEFAULT(KdTraits, PointAttrT, int, PointAttrT);

	/* Below types are automatically deduced. */
	using KdPoint         = PrimitiveWithAttribute<PointT, PointAttrT>;
	using KdPoints        = std::vector<KdPoint>;
	using KdPointsIter    = typename KdPoints::iterator;
	using KdPointPtr      = KdPoint *;
	using KdPointPtrs     = std::vector<KdPointPtr>;
	using KdPointPtrsIter = typename KdPointPtrs::iterator;

	class KdBox : public BboxT
	{
	public:
		KdBox() = default;
		KdBox(const KdPoint &minB, const KdPoint &maxB)
		  : BboxT(minB.primitive(), maxB.primitive())
		{
		}

		void update_from_point_pointers(KdPointPtrsIter begin, KdPointPtrsIter end)
		{
			if (begin == end)
				return;

			this->min_bound() = (**begin).primitive();
			this->max_bound() = (**begin).primitive();
			for (begin++; begin != end; begin++)
				(*this) += (**begin).primitive();
		}
	};
};



} // namespace OMC