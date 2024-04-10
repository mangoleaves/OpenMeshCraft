#pragma once

#include "TriangleSoup.h"

namespace OMC {

template <typename Traits>
class DetectClassifyTTIs
{
public:
	using NT         = typename Traits::NT;
	using EPoint     = typename Traits::EPoint;
	using GPoint     = typename Traits::GPoint;
	using IPoint_SSI = typename Traits::IPoint_SSI;
	using IPoint_LPI = typename Traits::IPoint_LPI;
	using IPoint_TPI = typename Traits::IPoint_TPI;
	using AsGP       = typename Traits::AsGP;
	using AsEP       = typename Traits::AsEP;
	using ToEP       = typename Traits::ToEP;
	using CreateSSI  = typename Traits::CreateSSI;
	using CreateLPI  = typename Traits::CreateLPI;
	using CreateTPI  = typename Traits::CreateTPI;

	using DoIntersect        = typename Traits::DoIntersect;
	using Orient3D           = typename Traits::Orient3D;
	using CollinearPoints3D  = typename Traits::CollinearPoints3D;
	using OrientOn2D         = typename Traits::OrientOn2D;
	using LessThan3D         = typename Traits::LessThan3D;
	using MaxCompInTriNormal = typename Traits::MaxCompInTriNormal;

	using PntArena = PointArena<Traits>;
	using TriSoup  = TriangleSoup<Traits>;

	using Tree = Arr_Tree_Intersection<Traits>;

	DetectClassifyTTIs(TriSoup &_ts, const Tree &_tree, bool _ignore_same_label,
	                   MeshArrangements_Stats *_stats, bool _verbose);

protected:
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
	                      bool cache_labels, CStyleVector<Label> &cached_labels);

	void propagateCoplanarTrianglesIntersections();

	bool pointInsideTriangle(index_t p_id, index_t t_id);

protected:
	TriSoup               &ts;
	std::vector<PntArena> &pnt_arenas;
	std::vector<IdxArena> &idx_arenas;

	const std::vector<Label> &labels;
	const Tree               &tree;

	/* ignore intersection between triangles with same label */
	bool                    ignore_same_label;
	/* statistics */
	MeshArrangements_Stats *stats;
	/* Behavior control flags */
	bool                    verbose;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "DetectClassifyTTIs.inl"
#endif