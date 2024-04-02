#pragma once

#include "Utils.h"

namespace OMC {

template <typename Traits>
class DetectBBI
{
public:
	using NT     = typename Traits::NT;
	using EPoint = typename Traits::EPoint;
	using GPoint = typename Traits::GPoint;

	// used to check Box-Box intersection.
	using DoIntersect = typename Traits::DoIntersect;

	using Tree = Arr_Tree_Intersection<Traits>;

	DetectBBI(const std::vector<GPoint *> &_verts,
	          const std::vector<index_t>  &_tris,
	          const std::vector<Label> &_labels, const size_t _num_labels,
	          const Tree &_tree, std::vector<UIPair> &_BBI_pairs,
	          bool _ignore_same_label, MeshArrangements_Stats *_stats = nullptr,
	          bool _verbose = false);

private:
	/**
	 * @brief Leaf nodes of OcTree are divided into two groups.
	 *
	 * If number of boxes in a leaf node is less than a specific number (e.g.,
	 * 1000), the node is assigned to \p small_leaf_nodes, otherwise it is
	 * assigned to \p large_leaf_nodes.
	 *
	 * Then, for \p small_leaf_nodes, we parallelize intersection detection on
	 * these nodes. For \p large_leaf_nodes, we collect pairs of boxes in one node
	 * and parallelize detection on these pairs.
	 *
	 * Why do we do that?
	 * We want the scale of parallelism is relatively large, hence there won't be
	 * significant overhead on parallelism.
	 * @param [in] leaf_nodes All leaf nodes of OcTree.
	 * @param [out] small_leaf_nodes	See explanation in brief.
	 * @param [out] large_leaf_nodes See explanation in brief.
	 */
	void partitionNodes(const std::vector<index_t> &leaf_nodes,
	                    std::vector<index_t>       &small_leaf_nodes,
	                    std::vector<index_t>       &large_leaf_nodes);

	/**
	 * @brief Parallelize intersection detection on \p nodes.
	 * @param leaf_nodes The small_leaf_nodes from partitionNodes.
	 */
	void parallelOnSmallNodes(const std::vector<index_t> &nodes);

	/**
	 * @brief Parallelize intersection detection on pairs of boxes in each node in
	 * \p nodes.
	 * @param leaf_nodes The large_leaf_nodes from partitionNodes.
	 */
	void parallelOnLargeNodes(const std::vector<index_t> &nodes);

	void cacheBoxesInNode(const typename Tree::Node             &node,
	                      CStyleVector<typename Tree::TreeBbox> &cached_boxes,
	                      bool                                   cache_labels,
	                      CStyleVector<Label>                   &cached_labels);

protected:
	/* Input data */
	const std::vector<GPoint *> &verts;
	const std::vector<index_t>  &tris;
	const std::vector<Label>    &labels;
	const size_t                 num_labels;
	const Tree                  &tree;
	/* Output data */
	std::vector<UIPair>         &BBI_pairs;
	size_t                       num_BBI_pairs;

protected:
	/* ignore intersection between triangles with same label */
	bool                    ignore_same_label;
	/* statistics */
	MeshArrangements_Stats *stats;
	/* Behavior control flags */
	bool                    verbose;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "DetectBBI.inl"
#endif