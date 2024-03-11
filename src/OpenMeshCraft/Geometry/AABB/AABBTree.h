#pragma once

#include "AABBNode.h"

#include <vector>

namespace OMC {

/**
 * @brief Axis-Aligned Bounding Box Tree.
 * @tparam Traits The traits defines tree's behaviors.
 */
template <typename Traits>
class AABBTree
{
public:
	/* The tratis should provide followings: */

	using NT        = typename Traits::NT;
	using PointT    = typename Traits::PointT;
	using PrimT     = typename Traits::PrimT;
	using Prims     = typename Traits::Prims;
	using PrimsIter = typename Traits::PrimsIter;
	using BboxT     = typename Traits::BboxT;

	using CalcBbox           = typename Traits::CalcBbox;
	using SplitPrims         = typename Traits::SplitPrims;
	using PrimReferencePoint = typename Traits::PrimReferencePoint;

	/* Types used by AABB tree */

	// AABB Node
	using NodeT    = AABBNode<PrimT, BboxT>;
	using NodePtr  = NodeT *;
	using NodeCPtr = const NodeT *;
	// container to store nodes
	using Nodes    = std::vector<NodeT>;

public:
	AABBTree() {}

	AABBTree(Prims &&primitives);

	AABBTree(PrimsIter first, PrimsIter beyond);

	~AABBTree();

	/// @brief Insert primitives to the tree. Won't call build.
	void insert(const Prims &primitives);
	/// @brief Insert primitives to the tree. Won't call build.
	void insert(Prims &&primitives);
	/// @brief Insert primitives to the tree. Won't call build.
	void insert(PrimsIter first, PrimsIter beyond);
	/// @brief Build the AABB tree with given primitives (by insert).
	void build();

	void clear();

	inline size_t size() const { return m_primitives.size(); }
	inline bool   empty() const { return m_primitives.empty(); }

	template <typename TraversalTrait>
	void traversal(TraversalTrait &traits) const;

private:
	// build functions
	void expand(NodePtr node, PrimsIter first, PrimsIter beyond);

	template <typename TraversalTrait>
	bool traversal_node(NodeCPtr node, TraversalTrait &traits,
	                    const size_t nb_primitives) const;

protected:
	Prims m_primitives;
	Nodes m_nodes;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "AABBTree.inl"
#endif