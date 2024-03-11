#pragma once

#include "Utils.h"

namespace OMC {

template <typename Traits>
class DetectIntersections
{
public:
	using NT     = typename Traits::NT;
	using EPoint = typename Traits::EPoint;
	using GPoint = typename Traits::GPoint;

	using AsEP = typename Traits::AsEP;

	using Orient3D    = typename Traits::Orient3D;
	using DoIntersect = typename Traits::DoIntersect;
	using Triangle3_Triangle3_DoIntersect =
	  typename Traits::Triangle3_Triangle3_DoIntersect;

	using Tree = Arr_OcTree_Intersection<Traits>;

	DetectIntersections(
	  const std::vector<GPoint *> &_verts, const std::vector<index_t> &_tris,
	  const std::vector<Label> &_labels, const size_t _num_labels, Tree &_tree,
	  std::vector<UIPair> &_intersection_list, std::vector<ShewchukCache> &_cache,
	  bool                    _ignore_inter_with_same_label,
	  MeshArrangements_Stats *_stats = nullptr, bool _verbose = false);

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
	 * We hope for a large-scale parallelism, hence there won't be significant
	 * overhead on parallelism.
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
	void parallelOnLeafNodes(const std::vector<index_t> &leaf_nodes);

	/**
	 * @brief Parallelize intersection detection on pairs of boxes in each node in
	 * \p leaf_nodes.
	 * @param leaf_nodes The dupl_leaf_nodes from partitionLeafNodes.
	 */
	void parallelOnUniqPairs(const std::vector<index_t> &leaf_nodes);

	/**
	 * @brief A simple wrap of triangle-triangle intersection test.
	 */
	bool intersectsTriangle(const EPoint &t1_v0, const EPoint &t1_v1,
	                        const EPoint &t1_v2, const EPoint &t2_v0,
	                        const EPoint &t2_v1, const EPoint &t2_v2,
	                        NT *t1_min = nullptr, NT *t1_perm = nullptr,
	                        NT *t2_min = nullptr, NT *t2_perm = nullptr);

protected:
	/* Input data */
	const std::vector<GPoint *> &verts;
	const std::vector<index_t>  &tris;
	const std::vector<Label>    &labels;
	const size_t                 num_labels;
	Tree                        &tree;
	/* Output data */
	std::vector<UIPair>         &intersection_list;

public:
	/* Middle auxiliary data */
	/* cached data, will be used by ClassifyIntersection later. */
	std::vector<ShewchukCache> &cache;
	std::vector<std::once_flag> cached;

protected:
	/* ignore intersection between triangles with same label */
	bool                    ignore_inter_with_same_label;
	/* statistics */
	MeshArrangements_Stats *stats;
	/* Behavior control flags */
	bool                    verbose;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "DetectIntersections.inl"
#endif