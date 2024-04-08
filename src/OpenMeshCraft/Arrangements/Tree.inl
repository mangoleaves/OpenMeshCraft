#pragma once

#include "Tree.h"

namespace OMC {

/********************************************************************/
/* Octree, used in detecting intersections                          */
/********************************************************************/

template <typename AppTraits>
void Arr_Tree_Intersection<AppTraits>::init_from_triangle_soup(
  const std::vector<GPoint *> &verts, const std::vector<index_t> &tris,
  const size_t reserve_tris, const MeshArrangements_Config &config)
{
	this->clear();

	size_t num_tris = tris.size() / 3;
	this->m_boxes.reserve(reserve_tris);
	this->m_boxes.resize(num_tris);
	tbb::parallel_for((size_t)0, num_tris,
	                  [this, &verts, &tris](size_t i)
	                  {
		                  BboxT box(AsEP()(*verts[tris[i * 3]]));
		                  box += AsEP()(*verts[tris[i * 3 + 1]]);
		                  box += AsEP()(*verts[tris[i * 3 + 2]]);
		                  this->m_boxes[i].bbox() = box;
		                  this->m_boxes[i].id()   = i;
	                  });
	// set parameters
	this->m_split_pred = Arr_TreeSplitPred(config.tree_split_size_thres);
	// construct the tree
	this->construct(/*compact*/ true, config.tree_enlarge_ratio,
	                config.tree_adaptive_thres);

#if 0 && defined(OMC_ARR_PROFILE)

	std::fstream fout;
	fout.open("./data/test_output/Arrangements/octree_nodes_size.txt",
	          std::ios::out | std::ios::app);
	if (fout.is_open())
	{
		std::vector<index_t> leaf_nodes = all_leaf_nodes();
		for (index_t ni : leaf_nodes)
		{
			if (node(ni).size() != 0)
				fout << node(ni).size() << std::endl;
		}
	}
	fout.close();

#endif
}

template <typename Tree, typename Node>
bool Arr_TreeShapeRefinePred::operator()(const Tree &tree, const Node &node,
                                         std::array<bool, 3> &partitionable)
{
	auto diag_length = (tree.box().max_bound() - tree.box().min_bound()).length();
	auto node_length = node.box().max_bound() - node.box().min_bound();
	partitionable[0] = node_length[0] > scale * diag_length;
	partitionable[1] = node_length[1] > scale * diag_length;
	partitionable[2] = node_length[2] > scale * diag_length;
	return partitionable[0] || partitionable[1] || partitionable[2];
}

template <typename AppTraits>
void Arr_Tree_Intersection<AppTraits>::shape_refine()
{
	this->m_shape_refine_pred = typename BaseT::ShapeRefinePred(0.1);

	BaseT::shape_refine();
}

template <typename AppTraits>
template <typename QPrimT>
void Arr_Tree_Intersection<AppTraits>::all_intersections(const QPrimT &query,
                                                         Indices &results) const
{
	BoxTrav box_trav(query);
	this->traversal(box_trav);
	results = box_trav.result();
}

template <typename AppTraits>
void Arr_Tree_Intersection<AppTraits>::insert_triangle(const GPoint *v0,
                                                       const GPoint *v1,
                                                       const GPoint *v2,
                                                       index_t       ins_id)
{
	BboxT box(AsEP()(*v0));
	box += AsEP()(*v1);
	box += AsEP()(*v2);
	insert_box(box, ins_id);
}

template <typename AppTraits>
void Arr_Tree_Intersection<AppTraits>::insert_box(const BboxT &ins_box,
                                                  index_t      ins_id)
{
	OMC_EXPENSIVE_ASSERT(this->m_nodes.size() != 0, "empty tree.");

	this->m_boxes.resize(this->m_boxes.size() + 1);
	this->m_boxes.back().bbox() = ins_box;
	this->m_boxes.back().id()   = ins_id;

	std::queue<index_t> nodes_to_insert;
	nodes_to_insert.push(this->root_node_idx());

	while (!nodes_to_insert.empty())
	{
		index_t                 nd_idx = nodes_to_insert.front();
		typename BaseT::NodeRef nd     = this->node(nd_idx);
		nodes_to_insert.pop();

		if (nd.is_leaf())
		{
			nd.boxes().push_back(this->m_boxes.size() - 1);
			nd.tbox() += ins_box;
			nd.tbox().min_bound().maximize(nd.box().min_bound());
			nd.tbox().max_bound().minimize(nd.box().max_bound());
			nd.size() += 1;
		}
		else // nd.is_internal()
		{
			for (index_t ch_i = 0; ch_i < nd.children_size(); ch_i++)
			{
				typename BaseT::NodeCRef ch_nd = this->node(nd.child(ch_i));
				if (DoIntersect()(ch_nd.box(), ins_box))
					nodes_to_insert.push(nd.child(ch_i));
			}
		}
	}
}

template <typename AppTraits>
auto Arr_Tree_Intersection<AppTraits>::locate_point(const GPoint *pp) ->
  typename BaseT::NodeRef
{ // in arrangements context, inserted point is always in one leaf node.
	OMC_EXPENSIVE_ASSERT(this->m_nodes.size() != 0, "empty tree.");

	typename BaseT::NodePtr nd_ptr = &this->root_node();
	while (nd_ptr->is_internal())
	{
		typename BaseT::TreePoint center = this->node_center(*nd_ptr);

		std::array<Sign, 3> cmp_sign = LessThan3D().on_all(*pp, center);

		size_t cmp_res = (static_cast<size_t>(cmp_sign[0] >= Sign::ZERO)) |
		                 (static_cast<size_t>(cmp_sign[1] >= Sign::ZERO) << 1) |
		                 (static_cast<size_t>(cmp_sign[2] >= Sign::ZERO) << 2);

		size_t ch_idx = nd_ptr->child_map()[cmp_res];

		nd_ptr = &this->node(nd_ptr->child(ch_idx));
	}
	return *nd_ptr;
}

template <typename AppTraits>
std::pair<std::atomic<index_t> *, bool>
Arr_Tree_Intersection<AppTraits>::insert_point(const GPoint         *pp,
                                               std::atomic<index_t> *ip)
{
	OMC_EXPENSIVE_ASSERT(this->m_nodes.size() != 0, "empty tree.");

	typename BaseT::NodeRef leaf_nd = this->locate_point(pp);

	// reach the target leaf node, insert point to concurrent map.
	return leaf_nd.attribute().point_map.insert(pp, ip);
}

template <typename AppTraits>
template <typename GetIndex>
std::pair<index_t, bool> Arr_Tree_Intersection<AppTraits>::insert_point_F(
  const GPoint *pp, std::atomic<index_t> *ip, GetIndex get_idx)
{
	OMC_EXPENSIVE_ASSERT(this->m_nodes.size() != 0, "empty tree.");

	typename BaseT::NodeRef leaf_nd = this->locate_point(pp);

	// reach the target leaf node, insert point to concurrent map.
	std::pair<std::atomic<index_t> *, bool> ins =
	  leaf_nd.attribute().point_map.insert(pp, ip);

	if (ins.second) // succeed to insert.
	{
		get_idx(pp, ip); // assign a valid index
		return std::pair<index_t, bool>(ip->load(), true);
	}
	else // the point already exists, fail to insert.
	{
		int spin_count = 1;
		// we may get an invalid index in case that the point is just inserted.
		// spin until the index is valid, the waiting time depends on get_idx.
		while (ins.first->load(std::memory_order_relaxed) == InvalidIndex)
		{
			tbb::detail::machine_pause(spin_count);
			spin_count = std::min(spin_count * 2, 16);
		}
		return std::pair<index_t, bool>(ins.first->load(std::memory_order_relaxed),
		                                false);
	}
}

template <typename AppTraits>
std::pair<index_t, bool>
Arr_Tree_Intersection<AppTraits>::find(const GPoint *pp)
{
	typename BaseT::NodeRef leaf_nd = this->locate_point(pp);
	return leaf_nd.attribute().point_map.find(pp);
}

template <typename AppTraits>
std::atomic<index_t> &Arr_Tree_Intersection<AppTraits>::at(const GPoint *pp)
{
	typename BaseT::NodeRef leaf_nd = this->locate_point(pp);
	return leaf_nd.attribute().point_map.at(pp);
}

template <typename AppTraits>
void Arr_Tree_Intersection<AppTraits>::clear_points()
{
	for (typename BaseT::NodeRef nd : this->m_nodes)
	{
		if (nd.is_leaf())
		{
			nd.attribute().point_map = AuxPointMap_ConcurrentMap<AppTraits>();
		}
	}
}

} // namespace OMC
