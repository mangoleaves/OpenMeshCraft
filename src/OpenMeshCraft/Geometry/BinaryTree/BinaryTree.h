#pragma once

#include "BinaryNode.h"

#include "OpenMeshCraft/Utils/CStyleVector.h"

#include "tbb/tbb.h"

#include <deque>
#include <queue>
#include <stack>

namespace OMC {

/**
 * @brief Binary tree.
 * @tparam _Traits _Traits WON'T be automatically deduced in the class.
 */
template <typename _Traits>
class BinaryTree
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

	/// size of children, constantly being 2 in binary tree.
	static constexpr size_t Degree = 2;

	using NT = typename Traits::NT;

	using Bbox = typename Traits::BboxT;

	using TreeBbox = typename Traits::TreeBboxT;
	BinaryTreeAbbreviate(TreeBbox);

	using TreeBboxes = CStyleVector<TreeBbox>;

	using Node = BinaryNode<Traits>;
	BinaryTreeAbbreviate(Node);
	using Nodes     = tbb::concurrent_vector<Node>;
	using NodesIter = typename Nodes::iterator;

	using TreeBboxesContainer = Node::BboxesContainer;

	using TreePoint =
	  remove_cvref_t<decltype(std::declval<TreeBbox>().min_bound())>;
	BinaryTreeAbbreviate(TreePoint);

	using CalcBbox    = typename Traits::CalcBbox;
	using DoIntersect = typename Traits::DoIntersect;

	using SplitPred   = typename Traits::SplitPred;
	using SplitManner = typename Traits::SplitManner;

	using ShapeRefinePred = typename Traits::ShapeRefinePred;

	using calc_box_from_boxes       = typename Traits::calc_box_from_boxes;
	using calc_box_from_box_indices = typename Traits::calc_box_from_box_indices;

public: /* Constructors and Destructor *************************************/
	BinaryTree() = default;

	BinaryTree(const BinaryTree &rhs) { shallow_copy(rhs); }

	BinaryTree &operator=(const BinaryTree &rhs)
	{
		shallow_copy(rhs);
		return *this;
	}

	/**
	 * @brief Only copy topology and shape, won't copy data and attributes.
	 */
	void shallow_copy(const BinaryTree &rhs);

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
	 */
	void construct(bool compact_box = false, NT enlarge_ratio = 1.2);

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

	size_t depth() const;

	std::vector<std::pair<index_t, size_t>> all_nodes_with_height() const;

	std::vector<index_t> all_nodes() const;

	std::vector<index_t> all_leaf_nodes() const;

public: /* Traversal *********************************************************/
	template <typename TraversalTrait>
	void traversal(TraversalTrait &traits) const;

	template <typename TraversalTrait>
	bool traversal_node(NodeCRef nd, TraversalTrait &traits) const;

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

	/// Split predicate, test whether a node need to split.
	SplitPred m_split_pred;

	/// Split manner, decide how to split a node.
	SplitManner m_split_manner;

	/// Shape refine predicate, test node's shape and decide whether to split it.
	ShapeRefinePred m_shape_refine_pred;

	/// Do intersect
	DoIntersect m_do_intersect;

	/// Calculate boungding box
	CalcBbox m_calc_bbox;

	/***** Behavior control flags and data *****/

	/// bounding box enlarge ratio
	NT m_enlarge_ratio = 1.5;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "BinaryTree.inl"
#endif