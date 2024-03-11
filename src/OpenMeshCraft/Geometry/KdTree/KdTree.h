#pragma once

#include "KdNode.h"
#include "KdPointContainer.h"
#include "KdTraversal.h"

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/ExtendedTypeTraits.h"

#include <algorithm>
#include <deque>
#include <vector>

namespace OMC {

/**
 * @brief KdTree. It now provide below functions:
 * 1. orthogonal nearest search.
 * @tparam KdTraits The traits provided to KdTree. See \ref KdAutoDeduceTraits
 * for how to define it.
 */
template <typename KdTraits>
class KdTree
{
public:
	using NT         = typename KdTraits::NT;
	using PointT     = typename KdTraits::PointT;
	using PointAttrT = typename KdTraits::PointAttrT;

	using KdPoint         = typename KdTraits::KdPoint;
	using KdPoints        = typename KdTraits::KdPoints;
	using KdPointsIter    = typename KdTraits::KdPointsIter;
	using KdPointPtr      = typename KdTraits::KdPointPtr;
	using KdPointPtrs     = typename KdTraits::KdPointPtrs;
	using KdPointPtrsIter = typename KdTraits::KdPointPtrsIter;
	using KdBox           = typename KdTraits::KdBox;

	using Node         = KdNode<KdTraits>;
	using NodePtr      = Node *;
	using InternalNode = KdInternalNode<KdTraits>;
	using InternalPtr  = InternalNode *;
	using LeafNode     = KdLeafNode<KdTraits>;
	using LeafPtr      = LeafNode *;

	using PointContainer = KdPointContainer<KdTraits>;

public:
	KdTree() = default;

	template <typename Points, typename Attrs>
	KdTree(const Points &points, const Attrs &attrs)
	{
		insert(points, attrs);
		build();
	}

	/**
	 * @brief Insert points and attrs to KdTree. Won't build KdTree.
	 * @param points points to initialize KdTree.
	 * @param attrs attributes that attached to points, e.g., the index of each
	 * point.
	 */
	template <typename Points, typename Attrs>
	void insert(const Points &points, const Attrs &attrs);

	/**
	 * @brief Build the KdTree.
	 */
	void build();

	inline bool empty() { return pts.empty(); }

	/**
	 * @brief Search the neatest point of the query in KdTree.
	 * @param query The query point.
	 * @return
	 */
	std::pair<PointT, PointAttrT> search_nearest_point(const PointT &query);

private:
	NodePtr create_leaf_node(PointContainer &container);

	NodePtr new_internal_node();

	/**
	 * @brief Recursively create the internal node of KdTree.
	 * @param n The internal node.
	 * @param c The point container of this node.
	 */
	void create_internal_node(NodePtr n, PointContainer &container);

	/**
	 * @brief Find the best separate dimension and separate value and then
	 * split points to two groups.
	 * @param[out] sep_dim The best separate dimension
	 * @param[out] sep_val The best separate value
	 * @param[inout] c_origin The points to split. It also stores the part of
	 * split result.
	 * @param[out] c_low It stores the other part of split result.
	 */
	void split(size_t &sep_dim, NT &sep_val, PointContainer &container_origin,
	           PointContainer &container_low);

	/**
	 * @brief Update the lower/higher low/high value of a node depending on its
	 * children's box.
	 */
	void handle_extended_node(InternalPtr nh, PointContainer &container,
	                          PointContainer &container_low);

private:
	std::deque<InternalNode> internal_nodes;
	std::deque<LeafNode>     leaf_nodes;

	NodePtr     tree_root;
	KdBox       bbox;
	KdPoints    pts;
	KdPointPtrs data;

	size_t bucket_size = 10;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "KdTree.inl"
#endif