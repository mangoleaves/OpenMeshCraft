#pragma once

#include "AuxStructure.h"
#include "Utils.h"

#include "OpenMeshCraft/Geometry/AdaptiveOrthTree/AdapOcTree.h"
#include "OpenMeshCraft/Geometry/AdaptiveOrthTree/AdapOrthTraversalTraits.h"

#include <fstream>
#include <vector>

namespace OMC {

/********************************************************************/
/* Split tree, used in triangulating single triangle                */
/* by splitting faces or edges.                                     */
/********************************************************************/

class SplitTree
{
public:
	struct Node
	{
		Node()
		  : v{InvalidIndex, InvalidIndex, InvalidIndex}
		  , c{InvalidIndex, InvalidIndex, InvalidIndex}
		{
		}

		Node(const index_t &_v0, const index_t &_v1, const index_t &_v2)
		  : v{_v0, _v1, _v2}
		  , c{InvalidIndex, InvalidIndex, InvalidIndex}
		{
		}

		index_t v[3]; // triangle vertices
		index_t c[3]; // child nodes, c2 may be invalid
		index_t ps;   // splitting point/vertex
	};

public:
	inline SplitTree() {}

	inline SplitTree(const size_t &size) { nodes.reserve(size); }

	inline index_t addNode(const index_t &v0, const index_t &v1,
	                       const index_t &v2)
	{
		nodes.emplace_back(v0, v1, v2);
		return nodes.size() - 1;
	}

	inline const Node &getNode(const index_t &node_id) const
	{
		OMC_EXPENSIVE_ASSERT(node_id < nodes.size(), "out fo range node id");
		return nodes[node_id];
	}

	inline void addChildren(const index_t &node_id, const index_t &c0,
	                        const index_t &c1)
	{
		OMC_EXPENSIVE_ASSERT(node_id < nodes.size(), "out fo range node id");
		OMC_EXPENSIVE_ASSERT(!is_valid_idx(nodes[node_id].c[0]),
		                     "assigning no empty children list");

		nodes[node_id].c[0] = c0;
		nodes[node_id].c[1] = c1;
	}

	inline void addChildren(const index_t &node_id, const index_t &c0,
	                        const index_t &c1, const index_t &c2)
	{
		OMC_EXPENSIVE_ASSERT(node_id < nodes.size(), "out fo range node id");
		OMC_EXPENSIVE_ASSERT(!is_valid_idx(nodes[node_id].c[0]),
		                     "assigning no empty children list");

		nodes[node_id].c[0] = c0;
		nodes[node_id].c[1] = c1;
		nodes[node_id].c[2] = c2;
	}

	inline void setSplitPoint(const index_t &node_id, const index_t &vid)
	{
		OMC_EXPENSIVE_ASSERT(node_id < nodes.size(), "out fo range node id");
		nodes[node_id].ps = vid;
	}

	inline index_t splitPoint(const index_t &node_id)
	{
		OMC_EXPENSIVE_ASSERT(node_id < nodes.size(), "out fo range node id");
		return nodes[node_id].ps;
	}

	inline void offsetVertesInNode(const index_t &node_id, size_t off)
	{
		OMC_EXPENSIVE_ASSERT(node_id < nodes.size(), "out fo range node id");
		Node   &nd = nodes[node_id];
		index_t v0 = nd.v[0], v1 = nd.v[1], v2 = nd.v[2];
		nd.v[(0 + off) % 3] = v0;
		nd.v[(1 + off) % 3] = v1;
		nd.v[(2 + off) % 3] = v2;
	}

private:
	std::vector<Node> nodes;
};

/********************************************************************/
/* Octree, used in detecting intersections                          */
/********************************************************************/

class Arr_TreeSplitPred
{
public:
	Arr_TreeSplitPred()
	  : threshold(1000)
	{
	}

	Arr_TreeSplitPred(size_t _threshold)
	  : threshold(_threshold)
	{
	}

	template <typename Tree, typename OrthNode>
	bool operator()(OMC_UNUSED const Tree &tree, const OrthNode &node)
	{
		return node.size() > threshold;
	}

private:
	size_t threshold;
};

class Arr_TreeShapeRefinePred
{
public:
	Arr_TreeShapeRefinePred()
	  : scale(1.)
	{
	}

	Arr_TreeShapeRefinePred(double _scale)
	  : scale(_scale)
	{
	}

	template <typename Tree, typename Node>
	bool operator()(const Tree &tree, const Node &node,
	                std::array<bool, 3> &partitionable);

private:
	double scale;
};

template <typename AppTraits>
class Arr_Tree_MinimalTraits
{
public:
	static constexpr size_t Dimension = 3;
	static constexpr size_t MaxDepth  = 16;

	using NT    = typename AppTraits::NT;
	using BboxT = typename AppTraits::BoundingBox;

	using SplitPred       = Arr_TreeSplitPred;
	using ShapeRefinePred = Arr_TreeShapeRefinePred;
	using DoIntersect     = typename AppTraits::DoIntersect;
	using CalcBbox        = typename AppTraits::CalcBbox;

	struct NodeAttrT
	{
		tbb::spin_mutex                      mutex;
		AuxPointMap_ConcurrentMap<AppTraits> point_map;

		NodeAttrT() = default;
		NodeAttrT(const NodeAttrT &rhs) { operator=(rhs); }
		NodeAttrT(NodeAttrT &&rhs) { operator=(std::move(rhs)); }
		NodeAttrT &operator=(const NodeAttrT &rhs)
		{
			point_map = rhs.point_map;
			return *this;
		}
		NodeAttrT &operator=(NodeAttrT &&rhs)
		{
			point_map = std::move(rhs.point_map);
			return *this;
		}
	};
};

template <typename AppTraits>
using Arr_TreeTraits =
  AdapOrthAutoDeduceTraits<Arr_Tree_MinimalTraits<AppTraits>>;

template <typename AppTraits>
class Arr_Tree_Intersection : public AdapOcTree<Arr_TreeTraits<AppTraits>>
{
public:
	using GPoint     = typename AppTraits::GPoint;
	using AsEP       = typename AppTraits::AsEP;
	using LessThan3D = typename AppTraits::LessThan3D;

	using TreeTraits = Arr_TreeTraits<AppTraits>;
	using BaseT      = AdapOcTree<TreeTraits>;
	using ThisT      = Arr_Tree_Intersection<AppTraits>;

	using NT    = typename TreeTraits::NT;
	using BboxT = typename TreeTraits::BboxT;

	using CalcBbox    = typename TreeTraits::CalcBbox;
	using DoIntersect = typename TreeTraits::DoIntersect;

	/// used to store indices of intersected primitives
	using Indices = std::vector<index_t>;

public: /* Constructors ***************************************************/
	Arr_Tree_Intersection(MeshArrangements_Stats *_stats)
	  : BaseT()
	  , stats(_stats)
	{
	}

	/// @brief it is shallow copy (see details in AdapOrthTree)
	Arr_Tree_Intersection(const Arr_Tree_Intersection &rhs)
	  : BaseT(rhs)
	  , stats(rhs.stats)
	{
	}

public: /* build and refine ***********************************************/
	void init_from_triangle_soup(const std::vector<GPoint *>   &verts,
	                             const std::vector<index_t>    &tris,
	                             const MeshArrangements_Config &config);

	void shape_refine(size_t num_intersection_pairs);

public: /* Query and update Interfaces ************************************/
	template <typename QPrimT>
	void all_intersections(const QPrimT &query, Indices &results) const;

	/// @brief insert new primitives in the tree. only insert it to existed nodes.
	template <typename PrimT>
	void insert_primitive(const PrimT &ins_prim, index_t ins_id)
	{
		insert_box(CalcBbox()(ins_prim), ins_id);
	}

	void insert_triangle(const GPoint *v0, const GPoint *v1, const GPoint *v2,
	                     index_t ins_id);

	/// @brief insert new box in the tree. only insert it to existed nodes.
	void insert_box(const BboxT &ins_box, index_t ins_id);

	typename BaseT::NodeRef locate_point(const GPoint *pp);

	std::pair<std::atomic<index_t> *, bool>
	insert_point(const GPoint *pp, std::atomic<index_t> *ip);

	/**
	 * @brief insert a point into tree.
	 * @param pp point's pointer
	 * @param ip index's pointer
	 * @param get_idx A function object that inserts point and get index.
	 * @return std::pair<index_t, bool> index of (existed) point and boolean
	 * indicating that whether the point is inserted.
	 */
	template <typename GetIndex>
	std::pair<index_t, bool>
	insert_point_F(const GPoint *pp, std::atomic<index_t> *ip, GetIndex get_idx);

	std::pair<index_t, bool> find(const GPoint *pp);

	std::atomic<index_t> &at(const GPoint *pp);

	void clear_points();

protected: /* Internal types, functions and data. **************************/
	         /* used in traversal */

	using BoxTrav = AdapOrth_BoxInterTraversal<TreeTraits>;

	/* used in statistics */

	MeshArrangements_Stats *stats = nullptr;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Tree.inl"
#endif