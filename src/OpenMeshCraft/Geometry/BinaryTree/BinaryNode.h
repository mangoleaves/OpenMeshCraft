#pragma once

#include "BinaryTraits.h"

#include "OpenMeshCraft/Utils/Exception.h"

#include <array>
#include <bitset>
#include <memory>
#include <ranges>
#include <vector>

namespace OMC {

/**
 * @brief Node in binary tree.
 * Node maintains a local topology of the tree (parent and children).
 * Node also contains geometry information, but never modify it.
 * @tparam Traits
 */
template <typename Traits>
class BinaryNode
{
public:
	/// The maximal depth of binary tree.
	/// Root node is at depth 0 and is counted into depth.
	/// For example, when MaxDepth is 2, the tree is allowed to have a root node
	/// (at depth 0) and root node's children (at depth 1). Deeper nodes are not
	/// allowd to exist in the tree.
	static constexpr size_t MaxDepth = Traits::MaxDepth;

	/// Dimension, typically 2 or 3, or higher n
	static constexpr size_t Dimension = Traits::Dimension;

	/// size of children, constantly being 2 in binary tree.
	static constexpr size_t Degree = 2;

	static_assert(MaxDepth <= 64);

	/// Geometry information
	using NT       = typename Traits::NT;
	using Bbox     = typename Traits::BboxT;
	using TreeBbox = typename Traits::TreeBboxT;
	BinaryTreeAbbreviate(TreeBbox);

	/// Attribute type defined on node.
	using NodeAttrT = typename Traits::NodeAttrT;

	/// Abbreviation of BinaryNode.
	using Node = BinaryNode<Traits>;
	BinaryTreeAbbreviate(Node);

	/// Two children of this node in binary tree
	/// lower: idx. higher: idx+1.
	using Children = index_t;

public:
#if 0
	using BboxesContainer = IotaView<index_t>;
#else
	using BboxesContainer = std::ranges::iota_view<index_t, index_t>;
#endif

public: /* Constructors (Copy, Move, Assign) and Destructor */
	BinaryNode();

	BinaryNode(const Node &src) { *this = src; }

	BinaryNode(Node &&src) { *this = std::move(src); }

	NodeRef operator=(const Node &src);

	NodeRef operator=(Node &&src);

	void shallow_copy(NodeCRef rhs);

	~BinaryNode() {}

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
	const size_t children_size() const { return is_internal() ? 2 : 0; }

	index_t       &depth() { return m_depth; }
	const index_t &depth() const { return m_depth; }

	/// Access loose bounding box of this node
	Bbox       &box() { return m_box; }
	const Bbox &box() const { return m_box; }

	/// Access tight bounding box of this node
	Bbox       &tbox() { return m_tbox; }
	const Bbox &tbox() const { return m_tbox; }

	/// Access geometry information
	BboxesContainer       &boxes() { return m_boxes; }
	const BboxesContainer &boxes() const { return m_boxes; }

	/// Access size
	size_t       &size() { return m_size; }
	const size_t &size() const { return m_size; }

	index_t       &split_dim() { return m_split_dim; }
	const index_t &split_dim() const { return m_split_dim; }

	NT       &split_coord() { return m_split_coord; }
	const NT &split_coord() const { return m_split_coord; }

	/// Access node attribute
	NodeAttrT       &attribute() { return m_attribute; }
	const NodeAttrT &attribute() const { return m_attribute; }

protected:
	/// index of parent, InvalidIndex if parent doesn't exist.
	index_t         m_parent;
	/// children of this node.
	Children        m_children;
	/// start from 0 (Root), end with MaxDepth - 1.
	index_t         m_depth;
	/// loose box of this node
	Bbox            m_box;
	/// tight box of this node
	Bbox            m_tbox;
	/// boxes interseted with this node.
	BboxesContainer m_boxes;
	/// size of all boxes stored in this node (including children).
	size_t          m_size;
	/// axis/dimension where this node is split
	index_t         m_split_dim;
	/// coordinate where this node is split
	NT              m_split_coord;
	/// attribute
	NodeAttrT       m_attribute;
};

template <typename Traits>
BinaryNode<Traits>::BinaryNode()
{
	m_parent   = InvalidIndex;
	m_children = InvalidIndex;
	m_depth    = 0;
	m_size     = 0;
}

/**
 * @brief Shallow copy.
 * @param src source node.
 */
template <typename Traits>
auto BinaryNode<Traits>::operator=(const Node &src) -> NodeRef
{
	m_parent    = src.m_parent;
	m_children  = src.m_children;
	m_depth     = src.m_depth;
	m_box       = src.m_box;
	m_tbox      = src.m_tbox;
	m_boxes     = src.m_boxes;
	m_size      = src.m_size;
	m_attribute = src.m_attribute;
	return *this;
}

/**
 * @brief Move (Shallow copy).
 * @param src source node.
 */
template <typename Traits>
auto BinaryNode<Traits>::operator=(Node &&src) -> NodeRef
{
	m_parent    = std::move(src.m_parent);
	m_children  = std::move(src.m_children);
	m_depth     = std::move(src.m_depth);
	m_box       = std::move(src.m_box);
	m_tbox      = std::move(src.m_tbox);
	m_boxes     = std::move(src.m_boxes);
	m_size      = std::move(src.m_size);
	m_attribute = std::move(src.m_attribute);
	return *this;
}

/**
 * @brief Only copy topology and box shape.
 * Won't copy contained data and attributes.
 * @return Node The copied node.
 */
template <typename Traits>
void BinaryNode<Traits>::shallow_copy(NodeCRef rhs)
{
	m_parent   = rhs.m_parent;
	m_children = rhs.m_children;
	m_depth    = rhs.m_depth;
	m_box      = rhs.m_box;
	m_tbox     = rhs.m_tbox;
	m_size     = 0;
}

/// @brief Access one child by local index.
template <typename Traits>
auto BinaryNode<Traits>::child(index_t index) const -> index_t
{
	OMC_EXPENSIVE_ASSERT(is_valid_idx(m_children), "leaf node has no child.");
	OMC_EXPENSIVE_ASSERT(index < children_size(), "index {} out of range.",
	                     index);
	return m_children + index;
}

} // namespace OMC