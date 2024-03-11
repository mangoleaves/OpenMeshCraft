#pragma once

#include "KdNode.h"

namespace OMC {

/**
 * @brief Orthogonal nearest search. If find the nearest point in KdTree.
 * @todo Support the case when find more than one neatest point.
 */
template <typename KdTraits>
class OrthogonalNearestSeach
{
public:
	using NT           = typename KdTraits::NT;
	using PointT       = typename KdTraits::PointT;
	using PointAttrT   = typename KdTraits::PointAttrT;
	using ProjectPoint = typename KdTraits::ProjectPoint;

	using KdPoint      = typename KdTraits::KdPoint;
	using KdPointsIter = typename KdTraits::KdPointsIter;
	using KdBox        = typename KdTraits::KdBox;

	using Node         = KdNode<KdTraits>;
	using NodePtr      = Node *;
	using InternalNode = KdInternalNode<KdTraits>;
	using InternalPtr  = InternalNode *;
	using LeafNode     = KdLeafNode<KdTraits>;
	using LeafPtr      = LeafNode *;

public:
	/**
	 * @brief Construct the object and find the nearest point at the same time.
	 * @param tree_root The pointer to root node of KdTree.
	 * @param tree_bbox The bounding box of the whole KdTree.
	 * @param query The query point.
	 */
	OrthogonalNearestSeach(NodePtr tree_root, const KdBox &tree_bbox,
	                       const PointT &query);

	inline KdPointsIter closest_point_iter() { return m_result; }

private:
	/**
	 * @brief Recursive function to traversal KdTree to find the nearest point.
	 * @param node The pointer to node of KdTree.
	 * @param sqrdis The distance to current closest point.
	 */
	void compute_nearest_neightbors_orthogonally(NodePtr node, NT sqrdis);

	/**
	 * @brief Search nearest point in leaf node.
	 * @param node A leaf node.
	 */
	void search_nearest_in_leaf(LeafPtr node);

	inline NT new_square_distance(NT dist, NT old_off, NT new_off)
	{
		return dist + (new_off * new_off - old_off * old_off);
	}

private:
	PointT m_dists;
	PointT m_query;

	KdPointsIter m_result;
	NT           m_square_distance;

	ProjectPoint m_project_point;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "KdTraversal.inl"
#endif