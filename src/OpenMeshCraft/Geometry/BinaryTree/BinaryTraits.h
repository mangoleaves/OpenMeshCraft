#pragma once

#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"
#include "OpenMeshCraft/Utils/IndexDef.h"
#include "OpenMeshCraft/Utils/SFINAE.h"

namespace OMC {

template <typename BinaryTraits>
class BinaryAutoDeduceTraits
{
public:
	/* Below values and types must be provided */
	static constexpr size_t Dimension = BinaryTraits::Dimension;

	// bounding box type
	using BboxT           = typename BinaryTraits::BboxT;
	// split predicate: when to split and when to stop.
	using SplitPred       = typename BinaryTraits::SplitPred;
	// split manner: where to split.
	using SplitManner     = typename BinaryTraits::SplitManner;
	// shape refine predicate: when to refine and when to stop.
	using ShapeRefinePred = typename BinaryTraits::ShapeRefinePred;
	// do intersect
	using DoIntersect     = typename BinaryTraits::DoIntersect;
	// calculate bounding box
	using CalcBbox        = typename BinaryTraits::CalcBbox;

	/* Try to deduce below values */
	GET_VALUE_OTHERWISE_DEFAULT(BinaryTraits, size_t, MaxDepth, 32, MaxDepth);

	/* Try to deduce below types */
	// number type
	GET_TYPE_OTHERWISE_DEFAULT(BinaryTraits, NT, double, NT);
	static_assert(!std::is_void_v<NT>,
	              "NT is not provided and can't be deduced.");

	// node attribute
	GET_TYPE_OTHERWISE_DEFAULT(BinaryTraits, NodeAttrT, double, NodeAttrT);

	static_assert(std::is_trivially_copyable<BboxT>::value);

	/* Below types are automatically deduced. */
	class TreeBboxT : public BboxT
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

#define BinaryTreeAbbreviate(Type) \
	using Type##Ptr  = Type *;       \
	using Type##CPtr = const Type *; \
	using Type##Ref  = Type &;       \
	using Type##CRef = const Type &

} // namespace OMC