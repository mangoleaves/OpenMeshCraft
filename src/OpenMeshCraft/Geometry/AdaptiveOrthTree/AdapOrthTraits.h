#pragma once

#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"
#include "OpenMeshCraft/Utils/IndexDef.h"
#include "OpenMeshCraft/Utils/SFINAE.h"

namespace OMC {

template <typename AdapOrthTraits>
class AdapOrthAutoDeduceTraits
{
public:
	/* Below values and types must be provided */
	static constexpr size_t Dimension = AdapOrthTraits::Dimension;

	// bounding box type
	using BboxT           = typename AdapOrthTraits::BboxT;
	// split predicate
	using SplitPred       = typename AdapOrthTraits::SplitPred;
	// split predicate
	using ShapeRefinePred = typename AdapOrthTraits::ShapeRefinePred;
	// do intersect
	using DoIntersect     = typename AdapOrthTraits::DoIntersect;
	// calculate bounding box
	using CalcBbox        = typename AdapOrthTraits::CalcBbox;

	/* Try to deduce below values */
	GET_VALUE_OTHERWISE_DEFAULT(AdapOrthTraits, size_t, MaxDepth, 32, MaxDepth);
	GET_VALUE_OTHERWISE_DEFAULT(AdapOrthTraits, bool, EnableVertices, false,
	                            EnableVertices);
	GET_VALUE_OTHERWISE_DEFAULT(AdapOrthTraits, bool, StoreBoxesInInternalNodes,
	                            false, StoreBoxesInInternalNodes);

	/* Try to deduce below types */
	// number type
	GET_TYPE_OTHERWISE_DEFAULT(AdapOrthTraits, NT, double, NT);
	static_assert(!std::is_void_v<NT>,
	              "NT is not provided and can't be deduced.");

	// node attribute
	GET_TYPE_OTHERWISE_DEFAULT(AdapOrthTraits, NodeAttrT, double, NodeAttrT);

	/* Below types are automatically deduced. */
	class OrBboxT : public BboxT
	{
	public:
		const BboxT &bbox() const { return *static_cast<const BboxT *>(this); }
		BboxT       &bbox() { return *static_cast<BboxT *>(this); }

		index_t  id() const { return m_id; }
		index_t &id() { return m_id; }

	protected:
		index_t m_id;
	};

	class calc_box_from_boxes
	{
	public:
		template <typename BoxIter>
		BboxT operator()(BoxIter begin, BoxIter end)
		{
			BboxT box;
			if (begin == end)
				return box;

			box = *begin;
			for (++begin; begin != end; ++begin)
				box += *begin;

			return box;
		}
	};

	class calc_box_from_box_indices
	{
	public:
		template <typename Tree, typename BoxIdxIter>
		BboxT operator()(const Tree &tree, BoxIdxIter begin, BoxIdxIter end)
		{
			BboxT box;
			if (begin == end)
				return box;

			box = tree.box(*begin);
			for (++begin; begin != end; ++begin)
				box += tree.box(*begin);

			return box;
		}
	};
};

#define AdapOrthTreeAbbreviate(Type) \
	using Type##Ptr  = Type *;         \
	using Type##CPtr = const Type *;   \
	using Type##Ref  = Type &;         \
	using Type##CRef = const Type &

} // namespace OMC