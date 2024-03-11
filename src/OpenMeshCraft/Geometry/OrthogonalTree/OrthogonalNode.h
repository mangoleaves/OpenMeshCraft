#pragma once

#include "OrthogonalVertex.h"

#include <array>
#include <bitset>
#include <memory>
#include <vector>

namespace OMC {

/**
 * @brief Node in orthogonal tree.
 * Node maintains a local topology of the tree (parent and children),
 * and provides some functions to modify and traverse the tree's topology.
 * Node also contains geometry information, but never modify it.
 * @tparam Traits
 */
template <typename Traits>
class OrthogonalNode
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

	/// enable vertices (corner vertices and face vertices).
	static constexpr bool EnableVertices = Traits::EnableVertices;

	static_assert(MaxDepth <= 64);

	/// Geometry information
	using Bbox   = typename Traits::BboxT;
	using OrBbox = typename Traits::OrBboxT;
	OrthTreeAbbreviate(OrBbox);

	/// Attribute type defined on node.
	using NodeAttrT = typename Traits::NodeAttrT;

	/// Abbreviation of OrthogonalNode.
	using Node = OrthogonalNode<Traits>;
	OrthTreeAbbreviate(Node);

	/// Abbreviation of OrthogonalVertex.
	using Vertex = OrthogonalVertex<Traits>;
	OrthTreeAbbreviate(Vertex);

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

	/**
	 * @brief Array storing corner vertices of this node.
	 * A node has and only has Degree corner vertices.
	 * @note A node hold a vertex, but the vertex may not be adjacent to the node.
	 * For example in 2D, the big cell holds 4 corner vertices (marked in `c`),
	 * but all 4 corner vertices are adjacent to small cells.
	 *          c-------f-------c
	 *          |       |       |
	 *          |       |       |
	 *          |       |       |
	 *          f-------+-------f
	 *          |       |       |
	 *          |       |       |
	 *          |       |       |
	 *          c-------f-------c
	 */
	using Vertices = std::array<index_t, Degree>;

	/**
	 * @brief unique pointer, allocate memory if need to store vertices.
	 */
	using VerticesUPtr = std::unique_ptr<Vertices>;

	/**
	 * @brief Represent this node's relationship to its parent.
	 * LocalCoord[i] indicates whether the i-th dimention of this node is
	 * lower (bit=0) or higher (bit=1).
	 * @note The local coordinates for root node is all zero.
	 */
	using LocalCoordinates = std::bitset<Dimension>;

	/// type used in global coordinates, we expect it to be as small as possible.
	using gc_t = std::conditional_t<MaxDepth <= 32, uint32_t, uint64_t>;
	/**
	 * @brief Represent this node's relatioship to the root.
	 * The global coordinates store all the local coordinates in the path from
	 * root to this node.
	 * If n is the depth of this node, LocalCoord_d is the local coordinates
	 * at depth d, 0 <= d <= n. The we have:
	 *       GlobalCoord[i][n-d] = LocalCoord_d[i].
	 * GlobalCoord is well defined when d>=1
	 */
	using GlobalCoordinates = std::array<gc_t, Dimension>;

public: /* Constructors (Copy, Move, Assign) and Destructor */
	OrthogonalNode();

	OrthogonalNode(const Node &src) { *this = src; }

	OrthogonalNode(Node &&src) { *this = std::move(src); }

	NodeRef operator=(const Node &src);

	NodeRef operator=(Node &&src);

	~OrthogonalNode() {}

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

	/// @brief Access one vertex by local index	(0 ~ Degree-1).
	index_t       &vertex(index_t index);
	/// @brief Access one vertex by local index	(0 ~ Degree-1).
	const index_t &vertex(index_t index) const;

	/// @brief Access vertices.
	Vertices       &vertices();
	/// @brief Access vertices.
	const Vertices &vertices() const;

	index_t       &depth() { return m_depth; }
	const index_t &depth() const { return m_depth; }

	bool             local_coordinates(size_t dim) const;
	void             local_coordinates(LocalCoordinates &lc) const;
	LocalCoordinates local_coordinates() const;

	/// @brief Access global coordinates.
	const GlobalCoordinates &global_coordinates() const
	{
		return m_global_coordinates;
	}
	/// @brief Access global coordinates.
	GlobalCoordinates &global_coordinates() { return m_global_coordinates; }

	/// Access bounding box of this node
	Bbox       &box() { return m_box; }
	const Bbox &box() const { return m_box; }

	/// Access geometry information
	std::vector<OrBboxCPtr>       &boxes() { return m_boxes; }
	const std::vector<OrBboxCPtr> &boxes() const { return m_boxes; }

	/// Access size
	size_t       &size() { return m_size; }
	const size_t &size() const { return m_size; }

	/// Access duplication degree
	float       &dupl_degree() { return m_dupl_degree; }
	const float &dupl_degree() const { return m_dupl_degree; }

	/// Access node attribute
	NodeAttrT       &attribute() { return m_attribute; }
	const NodeAttrT &attribute() const { return m_attribute; }

protected:
	/// index of parent, InvalidIndex if parent doesn't exist.
	index_t                 m_parent;
	/// children of this node.
	Children                m_children;
	/// vertices of this node.
	VerticesUPtr            m_vertices;
	/// global coordinates of this node.
	GlobalCoordinates       m_global_coordinates;
	/// start from 0 (Root), end with MaxDepth - 1.
	index_t                 m_depth;
	/// box of this node
	Bbox                    m_box;
	/// boxes interseted with this node.
	std::vector<OrBboxCPtr> m_boxes;
	/// size of all boxes stored in this node (including children).
	size_t                  m_size;
	/// duplication degree = sum of m_size of children / m_size of itself.
	float                   m_dupl_degree;
	/// attribute
	NodeAttrT               m_attribute;
};

template <typename Traits>
OrthogonalNode<Traits>::OrthogonalNode()
{
	m_parent   = InvalidIndex;
	m_children = InvalidIndex;
	m_depth    = 0;
	m_size     = 0;

	std::fill(m_global_coordinates.begin(), m_global_coordinates.end(), 0);

	if constexpr (EnableVertices)
		m_vertices = std::make_unique<Vertices>();
}

/**
 * @brief Shallow copy.
 * @param src source node.
 */
template <typename Traits>
auto OrthogonalNode<Traits>::operator=(const Node &src) -> NodeRef
{
	m_boxes              = src.m_boxes;
	m_parent             = src.m_parent;
	m_children           = src.m_children;
	m_global_coordinates = src.m_global_coordinates;
	m_depth              = src.m_depth;
	m_box                = src.m_box;
	m_boxes              = src.m_boxes;
	m_size               = src.m_size;
	m_attribute          = src.m_attribute;
	if (src.m_vertices)
		m_vertices = std::make_unique<Vertices>(*src.m_vertices);
	return *this;
}

/**
 * @brief Move (Shallow copy).
 * @param src source node.
 */
template <typename Traits>
auto OrthogonalNode<Traits>::operator=(Node &&src) -> NodeRef
{
	m_boxes              = std::move(src.m_boxes);
	m_parent             = std::move(src.m_parent);
	m_children           = std::move(src.m_children);
	m_vertices           = std::move(src.m_vertices);
	m_global_coordinates = std::move(src.m_global_coordinates);
	m_depth              = std::move(src.m_depth);
	m_box                = std::move(src.m_box);
	m_boxes              = std::move(src.m_boxes);
	m_size               = std::move(src.m_size);
	m_attribute          = std::move(src.m_attribute);
	return *this;
}

/// @brief Access one child by local index.
template <typename Traits>
auto OrthogonalNode<Traits>::child(index_t index) const -> index_t
{
	OMC_EXPENSIVE_ASSERT(is_valid_idx(m_children), "leaf node has no child.");
	OMC_EXPENSIVE_ASSERT(index < Degree, "index {} out of range.", index);
	return m_children + index;
}

/// @brief Access vertices.
/// There are Degree vertices in Degree subspaces.
/// They can be accessed by local coordinates.
template <typename Traits>
auto OrthogonalNode<Traits>::vertex(index_t index) -> index_t &
{
	OMC_EXPENSIVE_ASSERT(m_vertices, "Null vertices.");
	OMC_EXPENSIVE_ASSERT(index < Degree, "index {} out of range.", index);
	return (*m_vertices)[index];
}

/// @brief Access vertices.
/// There are Degree vertices in Degree subspaces.
/// They can be accessed by local coordinates.
template <typename Traits>
auto OrthogonalNode<Traits>::vertex(index_t index) const -> const index_t &
{
	OMC_EXPENSIVE_ASSERT(m_vertices, "Null vertices.");
	OMC_EXPENSIVE_ASSERT(index < Degree, "index {} out of range.", index);
	return (*m_vertices)[index];
}

/// @brief Access vertices.
template <typename Traits>
auto OrthogonalNode<Traits>::vertices() -> Vertices &
{
	OMC_EXPENSIVE_ASSERT(m_vertices, "Null vertices.");
	return *m_vertices;
}
/// @brief Access vertices.
template <typename Traits>
auto OrthogonalNode<Traits>::vertices() const -> const Vertices &
{
	OMC_EXPENSIVE_ASSERT(m_vertices, "Null vertices.");
	return *m_vertices;
}

/**
 * @brief Query the local coordinate at dimension \p dim.
 */
template <typename Traits>
bool OrthogonalNode<Traits>::local_coordinates(size_t dim) const
{
	OMC_EXPENSIVE_ASSERT(dim < Dimension, "invalid dimension {}", dim);
	return m_global_coordinates[dim] & gc_t(1);
}

/**
 * @brief Get the local coordinate and store in \p lc.
 */
template <typename Traits>
void OrthogonalNode<Traits>::local_coordinates(LocalCoordinates &lc) const
{
	// OPT: cost is expensive.
	if constexpr (Dimension >= 1)
		lc[0] = m_global_coordinates[0] & gc_t(1);
	if constexpr (Dimension >= 2)
		lc[1] = m_global_coordinates[1] & gc_t(1);
	if constexpr (Dimension >= 3)
		lc[2] = m_global_coordinates[2] & gc_t(1);
	if constexpr (Dimension > 3)
		for (size_t i = 3; i < Dimension; i++)
			lc[i] = m_global_coordinates[i] & gc_t(1);
}

/**
 * @brief Get the local coordinate.
 */
template <typename Traits>
auto OrthogonalNode<Traits>::local_coordinates() const -> LocalCoordinates
{
	LocalCoordinates result;
	local_coordinates(result);
	return result;
}

} // namespace OMC
