#pragma once

#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"
#include "OpenMeshCraft/Utils/IndexDef.h"
#include "OpenMeshCraft/Utils/SFINAE.h"

namespace OMC {

template <typename OrthTraits>
class OrthAutoDeduceTraits
{
public:
	/* Below values and types must be provided */
	static constexpr size_t Dimension = OrthTraits::Dimension;

	// bounding box type
	using BboxT       = typename OrthTraits::BboxT;
	// split predicate
	using SplitPred   = typename OrthTraits::SplitPred;
	// do intersect
	using DoIntersect = typename OrthTraits::DoIntersect;
	// calculate bounding box
	using CalcBbox    = typename OrthTraits::CalcBbox;

	/* Try to deduce below values */
	GET_VALUE_OTHERWISE_DEFAULT(OrthTraits, size_t, MaxDepth, 32, MaxDepth);
	GET_VALUE_OTHERWISE_DEFAULT(OrthTraits, bool, EnableVertices, false,
	                            EnableVertices);
	GET_VALUE_OTHERWISE_DEFAULT(OrthTraits, bool, StoreBoxesInInternalNodes,
	                            false, StoreBoxesInInternalNodes);

	/* Try to deduce below types */
	// number type
	GET_TYPE_OTHERWISE_DEFAULT(OrthTraits, NT, double, NT);
	static_assert(!std::is_void_v<NT>,
	              "NT is not provided and can't be deduced.");

	// node attribute
	GET_TYPE_OTHERWISE_DEFAULT(OrthTraits, NodeAttrT, double, NodeAttrT);

	// vertex attribute
	GET_TYPE_OTHERWISE_DEFAULT(OrthTraits, VertexAttrT, double, VertexAttrT);

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

			// TODO parallel by std::minmax_element(std::execution::unseq_par)?
			box = *begin;
			for (++begin; begin != end; ++begin)
				box += *begin;

			return box;
		}
	};

	class calc_box_from_box_pointers
	{
	public:
		template <typename BoxPtrIter>
		BboxT operator()(BoxPtrIter begin, BoxPtrIter end)
		{
			BboxT box;
			if (begin == end)
				return box;

			// TODO parallel by std::minmax_element(std::execution::unseq_par)?
			box = (*begin)->box();
			for (++begin; begin != end; ++begin)
				box += (*begin)->box();

			return box;
		}
	};
};

#define OrthTreeAbbreviate(Type)   \
	using Type##Ptr  = Type *;       \
	using Type##CPtr = const Type *; \
	using Type##Ref  = Type &;       \
	using Type##CRef = const Type &

} // namespace OMC