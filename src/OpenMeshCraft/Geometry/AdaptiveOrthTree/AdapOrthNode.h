#pragma once

#include "AdapOrthTraits.h"

#include "OpenMeshCraft/Utils/Exception.h"

#include <array>
#include <bitset>
#include <memory>
#include <vector>

namespace OMC {

/**
 * @brief Node in adaptive orthogonal tree.
 * Node maintains a local topology of the tree (parent and children).
 * Node also contains geometry information, but never modify it.
 * @tparam Traits
 */
template <typename Traits>
class AdapOrthNode
{
public:
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

	static_assert(MaxDepth <= 64);

	/// Geometry information
	using Bbox     = typename Traits::BboxT;
	using TreeBbox = typename Traits::TreeBboxT;
	AdapOrthTreeAbbreviate(TreeBbox);

	/// Attribute type defined on node.
	using NodeAttrT = typename Traits::NodeAttrT;

	/// Abbreviation of AdapOrthNode.
	using Node = AdapOrthNode<Traits>;
	AdapOrthTreeAbbreviate(Node);

	/** @brief Array storing a bundle of children.
	 * @details The orthtree subdivides the space in 2 on each dimension
	 * available, so a child node can be accessed by selecting a Boolean
	 * value for each dimension. The `index` parameter is thus
	 * interpreted as a bitmap, where each bit matches a dimension
	 * (starting by the least significant bit for coordinate X).
	 *
	 * For example, in the case of an octree (dimension 3):
	 *
	 * - index 0 (000 in binary) is the children on the "minimum corner" (xmin,
	 * ymin, zmin)
	 * - index 1 (001 in binary) is on (xmax, ymin, zmin)
	 * - index 2 (010 in binary) is on (xmin, ymax, zmin)
	 * - index 6 (101 in binary) is on (xmax, ymin, zmax)
	 * @note When splitting a node, allocate children in bundle, not one by one.
	 * @note We only use a pointer to represent children, we guarantee that
	 * children are sequantial in memory. It is almost equivalent to
	 *      using Children = std::array<NodePtr, Degree>;
	 */
	using Children = index_t;

public: /* Constructors (Copy, Move, Assign) and Destructor */
	AdapOrthNode();

	AdapOrthNode(const Node &src) { *this = src; }

	AdapOrthNode(Node &&src) { *this = std::move(src); }

	NodeRef operator=(const Node &src);

	NodeRef operator=(Node &&src);

	void shallow_copy(NodeCRef rhs);

	~AdapOrthNode() {}

public: /* Queries */
	bool is_root() const { return m_depth == 0; }

	bool is_internal() const { return is_valid_idx(m_children); }

	bool is_leaf() const { return !is_valid_idx(m_children); }

public: /* Data access */
	/// @brief Access parent.
	index_t       &parent() { return m_parent; }
	/// @brief Access parent.
	const index_t &parent() const { return m_parent; }

	/// @brief Access one child by local index (0 ~ Degree-1).
	index_t child(index_t index) const;

	/// @brief Access children (start index of children)
	index_t       &children() { return m_children; }
	/// @brief Access children (start index of children)
	const index_t &children() const { return m_children; }

	/// @brief Access children size
	size_t       &children_size() { return m_children_size; }
	/// @brief Access children size
	const size_t &children_size() const { return m_children_size; }

	index_t       &depth() { return m_depth; }
	const index_t &depth() const { return m_depth; }

	/// Access loose bounding box of this node
	Bbox       &box() { return m_box; }
	const Bbox &box() const { return m_box; }

	/// Access tight bounding box of this node
	Bbox       &tbox() { return m_tbox; }
	const Bbox &tbox() const { return m_tbox; }

	/// Access geometry information
	std::vector<index_t>       &boxes() { return m_boxes; }
	const std::vector<index_t> &boxes() const { return m_boxes; }

	/// Access size
	size_t       &size() { return m_size; }
	const size_t &size() const { return m_size; }

	/// children map
	const std::array<index_t, Degree> &child_map() const { return m_child_map; }
	std::array<index_t, Degree>       &child_map() { return m_child_map; }

	/// Access node attribute
	NodeAttrT       &attribute() { return m_attribute; }
	const NodeAttrT &attribute() const { return m_attribute; }

protected:
	/// index of parent, InvalidIndex if parent doesn't exist.
	index_t                     m_parent;
	/// children of this node.
	Children                    m_children;
	size_t                      m_children_size;
	/// start from 0 (Root), end with MaxDepth - 1.
	index_t                     m_depth;
	/// loose box of this node
	Bbox                        m_box;
	/// tight box of this node
	Bbox                        m_tbox;
	/// boxes interseted with this node.
	std::vector<index_t>        m_boxes;
	/// size of all boxes stored in this node (including children).
	size_t                      m_size;
	/// Map children to all corners
	std::array<index_t, Degree> m_child_map;
	/// attribute
	NodeAttrT                   m_attribute;
};

template <typename Traits>
AdapOrthNode<Traits>::AdapOrthNode()
{
	m_parent        = InvalidIndex;
	m_children      = InvalidIndex;
	m_children_size = 0;
	m_depth         = 0;
	m_size          = 0;
}

/**
 * @brief Shallow copy.
 * @param src source node.
 */
template <typename Traits>
auto AdapOrthNode<Traits>::operator=(const Node &src) -> NodeRef
{
	m_parent        = src.m_parent;
	m_children      = src.m_children;
	m_children_size = src.m_children_size;
	m_depth         = src.m_depth;
	m_box           = src.m_box;
	m_tbox          = src.m_tbox;
	m_boxes         = src.m_boxes;
	m_size          = src.m_size;
	m_child_map     = src.m_child_map;
	m_attribute     = src.m_attribute;
	return *this;
}

/**
 * @brief Move (Shallow copy).
 * @param src source node.
 */
template <typename Traits>
auto AdapOrthNode<Traits>::operator=(Node &&src) -> NodeRef
{
	m_parent        = std::move(src.m_parent);
	m_children      = std::move(src.m_children);
	m_children_size = std::move(src.m_children_size);
	m_depth         = std::move(src.m_depth);
	m_box           = std::move(src.m_box);
	m_tbox          = std::move(src.m_tbox);
	m_boxes         = std::move(src.m_boxes);
	m_size          = std::move(src.m_size);
	m_child_map     = std::move(src.m_child_map);
	m_attribute     = std::move(src.m_attribute);
	return *this;
}

/**
 * @brief Only copy topology and box shape.
 * Won't copy contained data and attributes.
 * @return Node The copied node.
 */
template <typename Traits>
void AdapOrthNode<Traits>::shallow_copy(NodeCRef rhs)
{
	m_parent        = rhs.m_parent;
	m_children      = rhs.m_children;
	m_children_size = rhs.m_children_size;
	m_depth         = rhs.m_depth;
	m_box           = rhs.m_box;
	m_tbox          = rhs.m_tbox;
	m_size          = 0;
	m_child_map     = rhs.m_child_map;
}

/// @brief Access one child by local index.
template <typename Traits>
auto AdapOrthNode<Traits>::child(index_t index) const -> index_t
{
	OMC_EXPENSIVE_ASSERT(is_valid_idx(m_children), "leaf node has no child.");
	OMC_EXPENSIVE_ASSERT(index < children_size(), "index {} out of range.",
	                     index);
	return m_children + index;
}

} // namespace OMC