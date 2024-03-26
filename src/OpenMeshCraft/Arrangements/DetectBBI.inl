#pragma once

#include "DetectBBI.h"

namespace OMC {

template <typename Traits>
DetectBBI<Traits>::DetectBBI(const std::vector<GPoint *> &_verts,
                             const std::vector<index_t>  &_tris,
                             const std::vector<Label>    &_labels,
                             const size_t _num_labels, const Tree &_tree,
                             std::vector<UIPair>    &_BBI_pairs,
                             bool                    _ignore_same_label,
                             MeshArrangements_Stats *_stats, bool _verbose)
  : verts(_verts)
  , tris(_tris)
  , labels(_labels)
  , num_labels(_num_labels)
  , tree(_tree)
  , BBI_pairs(_BBI_pairs)
  , ignore_same_label(_ignore_same_label)
  , stats(_stats)
  , verbose(_verbose)
{
	BBI_pairs.clear();
	BBI_pairs.reserve(_tris.size() / 3);

	std::vector<index_t> leaf_nodes = tree.all_leaf_nodes();

	std::vector<index_t> uniq_leaf_nodes, dupl_leaf_nodes;
	partitionLeafNodes(leaf_nodes, uniq_leaf_nodes, dupl_leaf_nodes);
	parallelOnUniqNodes(uniq_leaf_nodes);
	parallelOnDuplNodes(dupl_leaf_nodes);

	parallel_remove_duplicates(BBI_pairs);

	OMC_ARR_PROFILE_INC_REACH_CNT(ArrFuncNames::D_BBI, 0, BBI_pairs.size());

	if (verbose)
		Logger::info(std::format("[OpenMeshCraft Arrangements] Detect unique pairs "
		                         "of box-box intersections {}.",
		                         BBI_pairs.size()));
}

template <typename Traits>
void DetectBBI<Traits>::partitionLeafNodes(
  const std::vector<index_t> &leaf_nodes, std::vector<index_t> &uniq_leaf_nodes,
  std::vector<index_t> &dupl_leaf_nodes)
{
	uniq_leaf_nodes.reserve(leaf_nodes.size());
	dupl_leaf_nodes.reserve(leaf_nodes.size());

	for (index_t ni : leaf_nodes)
	{
		if (tree.node(ni).size() == 0)
			continue;
		if (tree.node(ni).size() > 400)
			dupl_leaf_nodes.push_back(ni);
		else
			uniq_leaf_nodes.push_back(ni);
	}

#ifdef OMC_ARR_PROFILE
	#if 0
	std::fstream fout;
	fout.open("./data/test_output/Arrangements/nodes_size.txt",
	          std::ios::out | std::ios::app);
	if (fout.is_open())
	{
		for (index_t ni : leaf_nodes)
		{
			if (tree.node(ni).size() != 0)
				fout << tree.node(ni).size() << std::endl;
		}
	}
	fout.close();
	#endif
#endif
}

template <typename Traits>
void DetectBBI<Traits>::parallelOnUniqNodes(
  const std::vector<index_t> &leaf_nodes)
{
	std::vector<std::vector<Label>> cache_labels_mt(
	  tbb::this_task_arena::max_concurrency());
	std::vector<std::vector<typename Tree::OrBbox>> cache_boxes_mt(
	  tbb::this_task_arena::max_concurrency());
	std::vector<std::vector<UIPair>> cache_check_pairs_mt(
	  tbb::this_task_arena::max_concurrency());
	std::vector<std::vector<UIPair>> check_pairs_mt(
	  tbb::this_task_arena::max_concurrency());

	auto evaluate_possible_BBI_pairs = [this](size_t s, index_t node_idx)
	{
		size_t num_boxes = tree.node(node_idx).size();
		size_t new_s     = num_boxes == 0
		                     ? size_t(0)
		                     : size_t(((num_boxes - 1) * num_boxes / 2) * 0.1);
		return s + new_s;
	};

	size_t possible_pairs =
	  std::accumulate(leaf_nodes.begin(), leaf_nodes.end(), size_t(0),
	                  evaluate_possible_BBI_pairs);
	for (std::vector<UIPair> &check_pairs : check_pairs_mt)
		check_pairs.reserve(std::max(
		  size_t(possible_pairs / tbb::this_task_arena::max_concurrency() * 2),
		  size_t(1000)));

	// Run intersection detection in a leaf node.
	auto on_leaf_node = [this, &cache_labels_mt, &cache_boxes_mt,
	                     &cache_check_pairs_mt, &check_pairs_mt](size_t node_idx)
	{
		int thread_id = tbb::this_task_arena::current_thread_index();

		std::vector<typename Tree::OrBbox> &cache_boxes = cache_boxes_mt[thread_id];
		std::vector<Label>  &cache_labels = cache_labels_mt[thread_id];
		std::vector<UIPair> &cache_pairs  = cache_check_pairs_mt[thread_id];
		std::vector<UIPair> &check_pairs  = check_pairs_mt[thread_id];

		const size_t num_boxes = tree.node(node_idx).size();
		const auto  &boxes     = tree.node(node_idx).boxes();

		// for better memory cache (it really saves a lot of time :) )
		cache_boxes.resize(num_boxes);
		for (size_t i = 0; i < num_boxes; i++)
			cache_boxes[i] = tree.box(boxes[i]);

		if (ignore_same_label)
		{
			cache_labels.resize(num_boxes);
			for (size_t i = 0; i < num_boxes; i++)
				cache_labels[i] = labels[tree.box(boxes[i]).id()];
		}

		cache_pairs.clear();
		cache_pairs.reserve(
		  std::max(size_t(num_boxes * (num_boxes - 1) / 2 * 0.2), size_t(1000)));

		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
		                              num_boxes * (num_boxes - 1) / 2);
		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI_UNIQ,
		                              num_boxes * (num_boxes - 1) / 2);

		for (index_t bi0 = 0; bi0 < num_boxes; bi0++)
		{
			const typename Tree::OrBbox &b0 = cache_boxes[bi0];

			for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
			{
				if ((ignore_same_label &&
				     (cache_labels[bi0] & cache_labels[bi1]).any()))
					continue;

				const typename Tree::OrBbox &b1 = cache_boxes[bi1];

				if (!DoIntersect()(b0.bbox(), b1.bbox()))
					continue; // early reject

				OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
				cache_pairs.push_back(uniquePair(b0.id(), b1.id()));
			}
		}

		check_pairs.insert(check_pairs.end(), cache_pairs.begin(),
		                   cache_pairs.end());
	};

#if 0
	tbb::parallel_for_each(leaf_nodes.begin(), leaf_nodes.end(), on_leaf_node);
#else
	std::for_each(leaf_nodes.begin(), leaf_nodes.end(), on_leaf_node);
#endif

	// Collect unique pairs
	size_t new_size = std::accumulate(
	  check_pairs_mt.begin(), check_pairs_mt.end(), size_t(0),
	  [](size_t s, const std::vector<UIPair> &cp) { return s + cp.size(); });

	OMC_ARR_PROFILE_INC_REACH_CNT(ArrFuncNames::D_BBI_UNIQ, 0, new_size);

	BBI_pairs.reserve(BBI_pairs.size() + new_size);
	for (std::vector<UIPair> &cp : check_pairs_mt)
		BBI_pairs.insert(BBI_pairs.end(), cp.begin(), cp.end());
}

template <typename Traits>
void DetectBBI<Traits>::parallelOnDuplNodes(
  const std::vector<index_t> &leaf_nodes)
{
	std::vector<Label>                 cache_labels;
	std::vector<typename Tree::OrBbox> cache_boxes;
	std::vector<std::vector<UIPair>>   check_pairs_mt(
    tbb::this_task_arena::max_concurrency());
	for (std::vector<UIPair> &check_pairs : check_pairs_mt)
		check_pairs.reserve(10000);

	for (index_t node_idx : leaf_nodes)
	{
		// sort boxes for faster traversal
		const auto  &boxes     = tree.node(node_idx).boxes();
		const size_t num_boxes = boxes.size();

		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
		                              num_boxes * (num_boxes - 1) / 2);
		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI_DUPL,
		                              num_boxes * (num_boxes - 1) / 2);

		{ // for better memory cache :>
			cache_labels.resize(num_boxes);
			cache_boxes.resize(num_boxes);
			for (size_t i = 0; i < num_boxes; i++)
			{
				const typename Tree::OrBbox &b = tree.box(boxes[i]);

				cache_boxes[i]  = b;
				cache_labels[i] = labels[b.id()];
			}
			for (std::vector<UIPair> &check_pairs : check_pairs_mt)
				check_pairs.clear();
		}

		// collect (possibly duplicate) pairs
		auto collect_pairs = [this, num_boxes, &cache_labels, &cache_boxes,
		                      &check_pairs_mt](index_t bi0)
		{
			index_t thread_id = tbb::this_task_arena::current_thread_index();
			std::vector<UIPair> &check_pairs = check_pairs_mt[thread_id];

			const typename Tree::OrBbox &b0 = cache_boxes[bi0];
			for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
			{
				if ((ignore_same_label &&
				     (cache_labels[bi0] & cache_labels[bi1]).any()))
					continue;

				const typename Tree::OrBbox &b1 = cache_boxes[bi1];

				if (!DoIntersect()(b0.bbox(), b1.bbox()))
					continue; // early reject

				OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
				check_pairs.push_back(uniquePair(b0.id(), b1.id()));
			}
		};
#if 1
		tbb::parallel_for(index_t(0), num_boxes, collect_pairs);
#else
		for (index_t i = 0; i < num_boxes; i++)
			collect_pairs(i);
#endif

		// Collect unique pairs
		size_t new_size = std::accumulate(
		  check_pairs_mt.begin(), check_pairs_mt.end(), size_t(0),
		  [](size_t s, const std::vector<UIPair> &cp) { return s + cp.size(); });

		OMC_ARR_PROFILE_INC_REACH_CNT(ArrFuncNames::D_BBI_DUPL, 0, new_size);

		for (std::vector<UIPair> &cp : check_pairs_mt)
			BBI_pairs.insert(BBI_pairs.end(), cp.begin(), cp.end());
	}
}

} // namespace OMC