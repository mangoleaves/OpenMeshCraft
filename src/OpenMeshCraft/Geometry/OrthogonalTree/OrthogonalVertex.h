#pragma once

#include "OrthogonalTraits.h"

#include "OpenMeshCraft/Utils/Exception.h"

#include <array>
#include <bitset>
#include <memory>

namespace OMC {

/**
 * @brief Vertex in orthogonal tree.
 * Vertex is the corner of nodes in orthogonal tree.
 *  - A node have 2^n vertices, where n is dimension of the tree.
 *  - A vertex is shared by at most 2^n nodes.
 *    An example to illustrate the size of adjacent nodes.
 *    ```
 *          1---------------2---------------1        ymax
 *          |               |               |        /|\
 *          |               |               |         |
 *          |               |               |         |
 *          |    [2][3]     |               |         |
 *          |               |               |         |
 *          |               |               |         |
 *          |               |               |         |
 *          2-------3-------4---3---3-------2         |
 *          |       |       |   |   |       |         |
 *          |  [0]  |  [1]  3---4---3       |         |
 *          |       |       |   |   |       |         |
 *          2---3---4-------4---3---4-------2         |
 *          |   |   |       |       |       |         |
 *          2---4---3       |       |       |         |
 *          |   |   |       |       |       |        \|/
 *          1---2---2-------2-------2-------1       ymin
 *
 *    xmin <---------------------------------> xmax
 *    ```
 *  - A vertex will store all adjacent nodes in an arrays of size 2^n.
 *    An adjacent node can be accessed by selecting a Boolean value for each
 *    dimension. The `index` parameter is thus interpreted as a bitmap, where
 *    each bit matches a dimension (starting by the least significant bit for
 *    coordinate X).
 *
 *    For example, in the case of an octree (dimension 2):
 *    - index [0] (00 in binary) is on (xmin, ymin)
 *    - index [1] (01 in binary) is on (xmax, ymin)
 *    - index [2] (10 in binary) is on (xmin, ymax)
 *    - index [3] (11 in binary) is on (xmax, ymax)
 *    An example is illustrate in the above figure.
 *
 *    NOTE Different indices may refer to the same node because the node is at
 *    smaller depth and hasn't been splitted.
 * @tparam NodeTraits
 */
template <typename Traits>
class OrthogonalVertex
{
public:
	/// Dimension, typically 2 or 3, or higher n
	static constexpr size_t Dimension = Traits::Dimension;

	/// A node has Degree vertices, a vertex is shared by at most Degree nodes.
	/// Degree is typically 4 or 8, or higher 2^n.
	static constexpr size_t Degree = (1ul << Dimension);

	using OrBbox  = typename Traits::OrBboxT;
	using OrPoint = remove_cvref_t<decltype(std::declval<OrBbox>().min_bound())>;

	using VertexAttrT = typename Traits::VertexAttrT;

public: /* Constructors and Destructor */
	OrthogonalVertex() = default;

public: /* Data access */
	OrPoint       &position() { return m_position; }
	const OrPoint &position() const { return m_position; }

	VertexAttrT       &attribute() { return m_attribute; }
	const VertexAttrT &attribute() const { return m_attribute; }

protected:
	OrPoint     m_position;
	VertexAttrT m_attribute;
};

} // namespace OMC
