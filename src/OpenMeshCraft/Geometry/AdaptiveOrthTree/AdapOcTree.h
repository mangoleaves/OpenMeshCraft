#pragma once

#include "AdapOrthTree.h"

#include <algorithm>

namespace OMC {

/**
 * @brief Adaptive octree, i.e., orthogonal tree in 3D.
 * @details Details about octree or orthogonal tree are in OrthogonalTree.
 */
template <typename _Traits>
class AdapOcTree : public AdapOrthTree<_Traits>
{
public:
	using Traits = _Traits;

	static constexpr size_t MaxDepth  = Traits::MaxDepth;
	static constexpr size_t Dimension = Traits::Dimension;
	static constexpr size_t Degree    = (1u << Dimension);

	using NT = typename Traits::NT;

	using Bbox = typename Traits::BboxT;

	using OrBbox = typename Traits::OrBboxT;
	AdapOrthTreeAbbreviate(OrBbox);

	using OrPoint = remove_cvref_t<decltype(std::declval<OrBbox>().min_bound())>;
	AdapOrthTreeAbbreviate(OrPoint);

	using Node = AdapOrthNode<Traits>;
	AdapOrthTreeAbbreviate(Node);

protected:
	virtual void
	calc_box_for_children(NodeRef nd, OrPointCRef center,
	                      std::array<Bbox, Degree> &child_boxes) final;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "AdapOcTree.inl"
#endif