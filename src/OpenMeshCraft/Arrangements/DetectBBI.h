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
	 * 1000), the node is assigned to \p uniq_leaf_nodes, otherwise it is assigned
	 * to \p dupl_leaf_nodes.
	 *
	 * Then, for \p uniq_leaf_nodes, we parallelize intersection detection on
	 * these nodes. For \p dupl_leaf_nodes, we collect pairs of boxes in one node
	 * and parallelize detection on these pairs.
	 *
	 * Why do we do that?
	 * We want the scale of parallelism is relatively large, hence there won't be
	 * significant overhead on parallelism.
	 * @param [in] leaf_nodes All leaf nodes of OcTree.
	 * @param [out] uniq_leaf_nodes	See explanation in brief.
	 * @param [out] dupl_leaf_nodes See explanation in brief.
	 */
	void partitionLeafNodes(const std::vector<index_t> &leaf_nodes,
	                        std::vector<index_t>       &uniq_leaf_nodes,
	                        std::vector<index_t>       &dupl_leaf_nodes);

	/**
	 * @brief Parallelize intersection detection on \p leaf_nodes.
	 * @param leaf_nodes The uniq_leaf_nodes from partitionLeafNodes.
	 */
	void parallelOnUniqNodes(const std::vector<index_t> &leaf_nodes);

	/**
	 * @brief Parallelize intersection detection on pairs of boxes in each node in
	 * \p leaf_nodes.
	 * @param leaf_nodes The dupl_leaf_nodes from partitionLeafNodes.
	 */
	void parallelOnDuplNodes(const std::vector<index_t> &leaf_nodes);

protected:
	/* Input data */
	const std::vector<GPoint *> &verts;
	const std::vector<index_t>  &tris;
	const std::vector<Label>    &labels;
	const size_t                 num_labels;
	const Tree                  &tree;
	/* Output data */
	std::vector<UIPair>         &BBI_pairs;

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