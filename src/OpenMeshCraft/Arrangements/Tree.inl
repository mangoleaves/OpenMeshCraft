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

#if defined(OMC_ARR_PROFILE) && defined(OMC_ARR_PROF_TREE_NODE)

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

} // namespace OMC
