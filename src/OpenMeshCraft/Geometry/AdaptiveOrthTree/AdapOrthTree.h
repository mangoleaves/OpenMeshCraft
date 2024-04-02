#pragma once

#include "AdapOrthNode.h"

#include "OpenMeshCraft/Utils/CStyleVector.h"

#include "tbb/tbb.h"

#include <deque>
#include <queue>

namespace OMC {

/**
 * @brief Adaptive orthogonal tree. In 2D it is mainly QuadTree, in 3D it is
 * mainly OcTree. But it will divide node to less children depending on whether
 * primitives are partitionable in all dimensions.
 * @tparam _Traits _Traits WON'T be automatically deduced in the class.
 */
template <typename _Traits>
class AdapOrthTree
{
public: /* Types and Declarations *******************************************/
	using Traits = _Traits;

	/// The maximal depth of orthogonal tree.
	/// Root node is at depth 0 and is counted into depth.
	/// For example, when MaxDepth is 2, the tree is allowed to have a root node
	/// (at depth 0) and root node's children (at depth 1). Deeper nodes are not
	/// allowd to exist in the tree.
	static constexpr size_t MaxDepth = Traits::MaxDepth;

	/// Dimension, typically 2 or 3, or higher n
	static constexpr size_t Dimension = Traits::Dimension;

	/// size of children, typically 4 or 8, or higher 2^n.
	static constexpr size_t Degree = (1u << Dimension);

	/// Store pointers to boxes in internal node.
	/// (Of course, pointers to boxes will always be stored in leaf nodes.)
	static constexpr bool StoreBoxesInInternalNodes =
	  Traits::StoreBoxesInInternalNodes;

	using NT = typename Traits::NT;

	using Bbox = typename Traits::BboxT;

	using TreeBbox = typename Traits::TreeBboxT;
	AdapOrthTreeAbbreviate(TreeBbox);
	using TreeBboxes = CStyleVector<TreeBbox>;
	using TreePoint =
	  remove_cvref_t<decltype(std::declval<TreeBbox>().min_bound())>;
	AdapOrthTreeAbbreviate(TreePoint);

	using CalcBbox    = typename Traits::CalcBbox;
	using DoIntersect = typename Traits::DoIntersect;

	using SplitPred = typename Traits::SplitPred;
	// using CollapsePred = typename Traits::CollapsePred;

	using ShapeRefinePred = typename Traits::ShapeRefinePred;

	using calc_box_from_boxes       = typename Traits::calc_box_from_boxes;
	using calc_box_from_box_indices = typename Traits::calc_box_from_box_indices;

	using Node = AdapOrthNode<Traits>;
	AdapOrthTreeAbbreviate(Node);

	using Nodes     = tbb::concurrent_vector<Node>;
	using NodesIter = typename Nodes::iterator;

public: /* Constructors and Destructor *************************************/
	AdapOrthTree() = default;

	AdapOrthTree(const AdapOrthTree &rhs) { shallow_copy(rhs); }

	AdapOrthTree &operator=(const AdapOrthTree &rhs)
	{
		shallow_copy(rhs);
		return *this;
	}

	/**
	 * @brief Only copy topology and shape, won't copy data and attributes.
	 */
	void shallow_copy(const AdapOrthTree &rhs);

	/**
	 * @brief Insert primitives and indices to tree. Will clear tree before
	 * inserting, Won't build tree until calling construction.
	 * @param primitives primitives to initialize tree. Won't store primitives
	 * internally, only store their boxes.
	 * @param indices indices of primitives.
	 * @note Once the primitives are inserted, they can't be changed, except
	 * clear, reinsert and rebuild.
	 */
	template <typename Primitives, typename Indices>
	void insert_primitives(const Primitives &primitives, const Indices &indices);

	/**
	 * @brief Insert bounding boxes and indices to tree. Will clear tree before
	 * inserting, Won't build tree until calling construction.
	 * @param bboxes bounding boxes to initialize tree. If your primitives are
	 * composed by hybrid types, construct bboxes for them and just insert bboxes.
	 * @param indices indices of bounding boxes.
	 * @note Once the bboxes are inserted, they can't be changed, except
	 * clear, reinsert and rebuild.
	 */
	template <typename Bboxes, typename Indices>
	void insert_boxes(const Bboxes &bboxes, const Indices &indices);

	/**
	 * @brief Construct the orthogonal tree after inserting points.
	 * @param compact_box if set to false, the box will be an equal length cube,
	 * otherwise be a unequal length cuboid to bound input as compact as possible.
	 * @param enlarge_ratio enlarge the bounding box of the inserting points.
	 * @param adaptive_thres control the adaptive strategy when splitting nodes.
	 */
	void construct(bool compact_box = false, NT enlarge_ratio = 1.2,
	               NT adaptive_thres = 0.1);

	/**
	 * @brief Refine shape of tree (e.g., split a node if it is too large).
	 * Won't reassign data(boxes).
	 */
	void shape_refine();

	/** @brief Clear the tree.  */
	void clear();

protected: /* Modifiers ******************************************************/
	/**
	 * @brief Allocate \p n_children nodes as new children.
	 * @return The index to the first child, other children can be
	 * accessed by offset the index, we guarantee that these children are
	 * sequential.
	 */
	index_t new_children(size_t n_children);

	/**
	 * @brief Split this node to children node.
	 * @param node_idx The node to split.
	 * @return true if this node is splitted.
	 */
	bool split(index_t node_idx);

	/**
	 * @brief Collaspe children node to its parent node.
	 * @param nd The node to collapse.
	 * @param collect_boxes If set to true, collect boxes in children.
	 * @todo How to re-cycle garbage?
	 */
	void collapse(index_t node_idx);

	/** @brief Assign boxes when splitting a node to children.  */
	void assign_boxes(NodeRef nd, TreePointCRef center,
	                  std::array<std::vector<index_t>, Degree> &assign_res,
	                  std::array<size_t, Dimension>            &lower,
	                  std::array<size_t, Dimension>            &higher);

	/**
	 * @brief Compare a box with a center point.
	 * @param box any box.
	 * @param center any center point.
	 * @return
	 * for the first bitset in pair: for any dimension i, if the box overlaps the
	 * lower half space of the center point, set bitset[i] to true.
	 * for the second bitset in pair: for any dimension i, if the box overlaps the
	 * higher half space of the center point, set bitset[i] to true.
	 */
	static std::pair<std::bitset<Dimension>, std::bitset<Dimension>>
	compare_box_with_center(TreeBboxCRef box, TreePointCRef center);

	/**
	 * @brief If boxes are not partitionable on dimension i, collapse nodes on
	 * dimension i. Boxes in collasped nodes will be moved to existed nodes.
	 * This function calculate destination for collapsed nodes.
	 * @param partitionable whether boxes are partitionable on all dimensions.
	 * @param destination boxes on i-th child node are moved to destination[i]-th
	 * child node.
	 * @return true if moving happens.
	 */
	void collapse_destination(const std::array<bool, Dimension> &partitionable,
	                          std::array<index_t, Degree>       &destination);

	virtual void calc_box_for_children(NodeRef nd, TreePointCRef center,
	                                   std::array<Bbox, Degree> &child_boxes) = 0;

	void calc_tbox_for_children(NodeRef nd);

public: /* Queries */
	/// @brief get the root box
	const Bbox &box() const { return m_bbox; }

	TreeBboxCRef box(index_t id) const { return m_boxes[id]; }
	TreeBboxRef  box(index_t id) { return m_boxes[id]; }

	index_t root_node_idx() const { return m_root_idx; }

	NodeRef  root_node() { return node(m_root_idx); }
	NodeCRef root_node() const { return node(m_root_idx); }

	NodeRef  node(index_t idx) { return m_nodes[idx]; }
	NodeCRef node(index_t idx) const { return m_nodes[idx]; }

	size_t size() const { return m_boxes.size(); }

	TreePoint node_center(NodeCRef nd) const;

	std::vector<index_t> all_leaf_nodes() const;

public: /* Traversal *********************************************************/
	template <typename TraversalTrait>
	void traversal(TraversalTrait &traits) const;

protected:
	template <typename TraversalTrait>
	bool traversal_node(NodeCRef nd, TraversalTrait &traits) const;

	// clang-format off
	template <typename Iter, typename OutIter, typename LessPred, typename EqualPred>
	void merge_unique(Iter b1, Iter e1, Iter b2, Iter e2, OutIter o, LessPred lp, EqualPred ep);
	// clang-format off

protected:
	/***** tree's kernel data *****/

	/// nodes of the orthogonal tree. root node is stored in m_nodes[0].
	Nodes m_nodes;

	constexpr static index_t m_root_idx = 0;

	tbb::spin_mutex m_new_children_mutex;

	/// Bounding box of all input primitives
	Bbox m_bbox;

	/// Boxes (contain boxes and indices).
	TreeBboxes m_boxes;

	/***** Predicates *****/

	/// Split predicate, check whether a node need to split.
	SplitPred m_split_pred;

	/// shape refine predicate, check node's shape and decide whether to split it.
	ShapeRefinePred m_shape_refine_pred;

	/// Do intersect
	DoIntersect m_do_intersect;

	/// Calculate boungding box
	CalcBbox m_calc_bbox;

	/***** Behavior control flags and data *****/

	/// bounding box enlarge ratio
	NT m_enlarge_ratio = 1.5;

	NT m_adaptive_thres = 0.1;
};



} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "AdapOrthTree.inl"
#endif