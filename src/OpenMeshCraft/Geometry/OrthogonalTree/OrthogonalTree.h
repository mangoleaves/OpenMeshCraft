#pragma once

#include "OrthogonalNode.h"

#include <deque>
#include <queue>

namespace OMC {

/**
 * @brief Orthogonal tree. In 2D it is QuadTree, in 3D it is OcTree.
 * @tparam _Traits _Traits WON'T be automatically deduced in the class.
 */
template <typename _Traits>
class OrthogonalTree
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

	/// Enable vertices stored in tree and nodes.
	static constexpr bool EnableVertices = Traits::EnableVertices;

	/// Store pointers to boxes in internal node.
	/// (Of course, pointers to boxes will always be stored in leaf nodes.)
	static constexpr bool StoreBoxesInInternalNodes =
	  Traits::StoreBoxesInInternalNodes;

	using NT = typename Traits::NT;

	using Bbox = typename Traits::BboxT;

	using OrBbox = typename Traits::OrBboxT;
	OrthTreeAbbreviate(OrBbox);

	using OrBboxes       = std::vector<OrBbox>;
	using OrBboxesIter   = typename OrBboxes::iterator;
	using OrBboxPtrs     = std::vector<OrBboxPtr>;
	using OrBboxPtrsIter = typename OrBboxPtrs ::iterator;

	using OrPoint = remove_cvref_t<decltype(std::declval<OrBbox>().min_bound())>;
	OrthTreeAbbreviate(OrPoint);

	using CalcBbox    = typename Traits::CalcBbox;
	using DoIntersect = typename Traits::DoIntersect;

	using SplitPred = typename Traits::SplitPred;
	// using CollapsePred = typename Traits::CollapsePred;

	using calc_box_from_boxes = typename Traits::calc_box_from_boxes;
	using calc_box_from_box_pointers =
	  typename Traits::calc_box_from_box_pointers;

	using Node = OrthogonalNode<Traits>;
	OrthTreeAbbreviate(Node);

	using LocalCoordinates  = typename Node::LocalCoordinates;
	using GlobalCoordinates = typename Node::GlobalCoordinates;

	using Vertex = OrthogonalVertex<Traits>;
	OrthTreeAbbreviate(Vertex);

public: /* Constructors and Destructor *************************************/
	OrthogonalTree() = default;

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
	 * @param dupl_thres if a node's duplication degree is larger than dupl_thres,
	 * stop splitting this node.
	 * @param depth_delta the difference between adjacent nodes is less than
	 * depth_delta.	if set to MaxDepth or a larger number, gradation won't be
	 * operated.
	 */
	void construct(bool compact_box = false, NT enlarge_ratio = 1.2,
	               NT dupl_thres = 2., index_t depth_delta = 2);

	/**
	 * @brief Clear the tree.
	 */
	void clear();

protected: /* Modifiers ******************************************************/
	/**
	 * @brief Allocate `Degree` nodes as new children.
	 * @return The index to the first child, other children can be
	 * accessed by offset the index, we guarantee that these children are
	 * sequential.
	 */
	index_t new_children();

	/**
	 * @brief Create a new vertex and return the index of it.
	 * @return index_t the index of the new vertex.
	 */
	index_t new_vertex();

	/**
	 * @brief refines the orthtree such that the difference of depth between two
	 * immediate neighbor leaves is never more than given depth delta.
	 */
	void grade(NT dupl_thres, index_t depth_delta);

	/**
	 * @brief A child node inherits neccessary information from its parent.
	 * @param child the child node.
	 * @param parent the parent node.
	 * @param local_coordinates the local coordinates of this child node with
	 * respect to its parent node.
	 */
	void child_inherit_parent(index_t child_idx, index_t parent_idx,
	                          const LocalCoordinates &local_coord);

	/**
	 * @brief Split this node to children node.
	 * @param node_idx The node to split.
	 */
	void split(index_t node_idx);

	/**
	 * @brief Collaspe children node to its parent node.
	 * @param node_idx The node to collapse.
	 */
	void collapse(index_t node_idx);

	/** @brief Assign boxes when splitting a node to children.  */
	void assign_boxes(NodeRef nd, OrPointCRef center);

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
	compare_box_with_center(OrBboxCRef box, OrPointCRef center);

	virtual void build_vertices() = 0;

	virtual void calc_box_for_children(NodeRef nd, OrPointCRef center) = 0;

public: /* Queries */
	index_t depth() const { return m_side_length_per_depth.size() - 1; }

	index_t root_node_idx() const { return m_root_idx; }

	NodeRef  root_node() { return node(m_root_idx); }
	NodeCRef root_node() const { return node(m_root_idx); }

	NodeRef  node(index_t idx) { return m_nodes[idx]; }
	NodeCRef node(index_t idx) const { return m_nodes[idx]; }

	size_t size() const { return m_boxes.size(); }

	/**
	 * @brief Calculate the center of a given node by its global coordinates and
	 * side length.
	 */
	OrPoint node_center(NodeCRef nd) const;

	/** @brief Get the side length of a node. */
	OrPoint node_side_length(NodeCRef nd) const;

	std::deque<Vertex>       &vertices() { return m_vertices; }
	const std::deque<Vertex> &vertices() const { return m_vertices; }

	VertexRef  vertex(index_t idx) { return m_vertices[idx]; }
	VertexCRef vertex(index_t idx) const { return m_vertices[idx]; }

	std::vector<index_t> all_leaf_nodes() const;

	SplitPred       &split_pred() { return m_split_pred; }
	const SplitPred &split_pred() const { return m_split_pred; }

	NT &duplication_threshold() { return m_dupl_thres; }
	NT  duplication_threshold() const { return m_dupl_thres; }

public: /* Traversal and adjacency *******************************************/
	/**
	 * @brief Find the adjacent node along the direction in specific dimension.
	 * @details Adjacent nodes are found according to several properties:
	 *  - adjacent nodes may be shallower than the seek node, but never deeper.
	 *  - a node has at most `2 * Dimension` different adjacent nodes
	 *    (e.g., in 3D: left, right, up, down, front, back).
	 *  - adjacent nodes are not required to be leaf nodes.
	 *
	 * Here's a diagram demonstrating the concept for a Quadtree
	 * (kindly provided by CGAL):
	 *
	 *  ```
	 *  +---------------+---------------+
	 *  |               |               |
	 *  |               |               |
	 *  |               |               |
	 *  |       A       |               |
	 *  |               |               |
	 *  |               |               |
	 *  |               |               |
	 *  +-------+-------+---+---+-------+
	 *  |       |       |   |   |       |
	 *  |   A   |  (S)  +---A---+       |
	 *  |       |       |   |   |       |
	 *  +---+---+-------+---+---+-------+
	 *  |   |   |       |       |       |
	 *  +---+---+   A   |       |       |
	 *  |   |   |       |       |       |
	 *  +---+---+-------+-------+-------+
	 *  ```
	 *
	 *  - (S) : Seek node
	 *  - A  : Adjacent node
	 *
	 * Note how the top adjacent node is larger than the seek node.
	 * The right adjacent node is the same size, even though it contains
	 * further subdivisions.
	 *
	 * This implementation returns the adjacent node if it's found.  If
	 * there is no adjacent node in that direction, it returns a null
	 * node.
	 * @param dim dimension
	 * @param dir dirtion: 0 (negative along the dimension), 1 (positive along the
	 * dimension)
	 * @return NodePtr
	 */
	index_t adjacent_node(NodeCRef nd, index_t dim, bool dir) const;

	/**
	 * @brief Check if two Node is topologically same.
	 */
	bool is_topologically_same(NodeCRef lhs, NodeCRef rhs) const;

	/**
	 * @brief find the box where the point locates.
	 * @param point query point
	 * @return index_t index to the box.
	 */
	template <typename PointT>
	index_t locate(const PointT &point) const;

	template <typename TraversalTrait>
	void traversal(TraversalTrait &traits) const;

public: /* Auxiliary or debug utils ******************************************/
	const std::vector<size_t> &assign_count() const { return m_assign_cnt; }

protected:
	template <typename TraversalTrait>
	bool traversal_node(NodeCRef nd, TraversalTrait &traits) const;

protected:
	/***** tree's kernel data *****/

	/// nodes of the orthogonal tree. root node is stored in m_nodes[0].
	std::deque<Node> m_nodes;

	constexpr static index_t m_root_idx = 0;

	/// Bounding box of all input points
	Bbox m_bbox;

	/// Side length for nodes at different depth
	std::vector<OrPoint> m_side_length_per_depth;

	/// Boxes (contain boxes and indices).
	OrBboxes m_boxes;

	/// (unique) vertices of boxes
	std::deque<Vertex> m_vertices;

	/***** Predicates *****/

	/// Split predicate, check whether a node need to split.
	SplitPred m_split_pred;

	/// Do intersect
	DoIntersect m_do_intersect;

	/// Calculate boungding box
	CalcBbox m_calc_bbox;

	/***** Auxiliary data *****/

	/// counts how many times a box is assigned to different leaf nodes
	std::vector<size_t> m_assign_cnt;

	/***** Behavior control flags and data *****/

	/// bounding box enlarge ratio
	NT m_enlarge_ratio = 1.5;

	/// duplication threshold
	NT m_dupl_thres = 8.0;

	/// maximal depth delta between adjacent leaf nodes
	size_t m_depth_delta = 2;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "OrthogonalTree.inl"
#endif