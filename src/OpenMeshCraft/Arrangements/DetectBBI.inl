#pragma once

#include "DetectBBI.h"

#define OMC_ARR_DETECT_BBI_PARA
#define OMC_ARR_COLLECT_BBI_PAIRS

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

	std::vector<index_t> small_nodes, large_nodes;
	partitionNodes(leaf_nodes, small_nodes, large_nodes);
	parallelOnSmallNodes(small_nodes);
	parallelOnLargeNodes(large_nodes);

#ifdef OMC_ARR_COLLECT_BBI_PAIRS
	parallel_remove_duplicates(BBI_pairs);
#endif

	OMC_ARR_PROFILE_INC_REACH_CNT(ArrFuncNames::D_BBI, 0, BBI_pairs.size());

	if (verbose)
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
		Logger::info(std::format("[OpenMeshCraft Arrangements] Detect unique pairs "
		                         "of box-box intersections {}.",
		                         BBI_pairs.size()));
#else
		Logger::info(std::format("[OpenMeshCraft Arrangements] Detect unique pairs "
		                         "of box-box intersections {}.",
		                         num_BBI_pairs));
#endif
}

template <typename Traits>
void DetectBBI<Traits>::partitionNodes(const std::vector<index_t> &nodes,
                                       std::vector<index_t>       &small_nodes,
                                       std::vector<index_t>       &large_nodes)
{
	small_nodes.reserve(nodes.size());
	large_nodes.reserve(nodes.size());

	for (index_t ni : nodes)
	{
		if (tree.node(ni).size() == 0)
			continue;
		if (tree.node(ni).size() > 400)
			large_nodes.push_back(ni);
		else
			small_nodes.push_back(ni);
	}
}

template <typename Traits>
void DetectBBI<Traits>::parallelOnSmallNodes(const std::vector<index_t> &nodes)
{
	// Run intersection detection in a node.
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
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

	size_t possible_pairs = std::accumulate(nodes.begin(), nodes.end(), size_t(0),
	                                        evaluate_possible_BBI_pairs);
	for (std::vector<UIPair> &check_pairs : check_pairs_mt)
		check_pairs.reserve(std::max(
		  size_t(possible_pairs / tbb::this_task_arena::max_concurrency() * 2),
		  size_t(1000)));
#endif

	// Run intersection detection in a leaf node.
	auto on_small_node = [&](size_t node_idx)
	{
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
		int thread_id = tbb::this_task_arena::current_thread_index();
		std::vector<UIPair> &cache_pairs = cache_check_pairs_mt[thread_id];
		std::vector<UIPair> &check_pairs = check_pairs_mt[thread_id];
#endif

		const typename Tree::Node &node = tree.node(node_idx);

		CStyleVector<typename Tree::TreeBbox> cache_boxes;
		CStyleVector<Label>                   cache_labels;

		cacheBoxesInNode(node, cache_boxes, ignore_same_label, cache_labels);

		const size_t num_boxes = tree.node(node_idx).size();

#ifdef OMC_ARR_COLLECT_BBI_PAIRS
		cache_pairs.clear();
		cache_pairs.reserve(
		  std::max(size_t(num_boxes * (num_boxes - 1) / 2 * 0.5), size_t(1000)));
#else
		size_t local_num_BBI_pairs = 0;
#endif

		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
		                              num_boxes * (num_boxes - 1) / 2);
		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI_UNIQ,
		                              num_boxes * (num_boxes - 1) / 2);

		if (ignore_same_label)
		{
			for (index_t bi0 = 0; bi0 < num_boxes; bi0++)
			{
				const typename Tree::TreeBbox &b0 = cache_boxes[bi0];

				for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
				{
					if ((cache_labels[bi0] & cache_labels[bi1]).any())
						continue;

					const typename Tree::TreeBbox &b1 = cache_boxes[bi1];

					if (!DoIntersect()(b0.bbox(), b1.bbox()))
						continue; // early reject

					OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
					cache_pairs.push_back(uniquePair(b0.id(), b1.id()));
#else
					local_num_BBI_pairs++;
#endif
				}
			}
		}
		else
		{
			for (index_t bi0 = 0; bi0 < num_boxes; bi0++)
			{
				const typename Tree::TreeBbox &b0 = cache_boxes[bi0];

				for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
				{
					const typename Tree::TreeBbox &b1 = cache_boxes[bi1];

					if (!DoIntersect()(b0.bbox(), b1.bbox()))
						continue; // early reject

					OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
					cache_pairs.push_back(uniquePair(b0.id(), b1.id()));
#else
					local_num_BBI_pairs++;
#endif
				}
			}
		}
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
		check_pairs.insert(check_pairs.end(), cache_pairs.begin(),
		                   cache_pairs.end());
#else
		num_BBI_pairs += local_num_BBI_pairs;
#endif
	};

#ifdef OMC_ARR_DETECT_BBI_PARA
	tbb::parallel_for_each(nodes.begin(), nodes.end(), on_small_node);
#else
	std::for_each(nodes.begin(), nodes.end(), on_small_node);
#endif

#ifdef OMC_ARR_COLLECT_BBI_PAIRS
	// Collect unique pairs
	OMC_UNUSED size_t new_size = std::accumulate(
	  check_pairs_mt.begin(), check_pairs_mt.end(), size_t(0),
	  [](size_t s, const std::vector<UIPair> &cp) { return s + cp.size(); });

	OMC_ARR_PROFILE_INC_REACH_CNT(ArrFuncNames::D_BBI_UNIQ, 0, new_size);

	BBI_pairs.reserve(BBI_pairs.size() + new_size);
	for (std::vector<UIPair> &cp : check_pairs_mt)
		BBI_pairs.insert(BBI_pairs.end(), cp.begin(), cp.end());
#endif
}

template <typename Traits>
void DetectBBI<Traits>::parallelOnLargeNodes(const std::vector<index_t> &nodes)
{
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
	std::vector<std::vector<UIPair>> check_pairs_mt(
	  tbb::this_task_arena::max_concurrency());
	for (std::vector<UIPair> &check_pairs : check_pairs_mt)
		check_pairs.reserve(10000);
#endif

	size_t num_threads = tbb::this_task_arena::max_concurrency();

	for (index_t node_idx : nodes)
	{
		const typename Tree::Node &node = tree.node(node_idx);

		size_t num_boxes = node.size();

		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
		                              num_boxes * (num_boxes - 1) / 2);
		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI_DUPL,
		                              num_boxes * (num_boxes - 1) / 2);

		CStyleVector<typename Tree::TreeBbox> cache_boxes;
		CStyleVector<Label>                   cache_labels;
		cacheBoxesInNode(node, cache_boxes, ignore_same_label, cache_labels);

		// check Box-Box-Intersection in current node.
		auto on_large_node = [&](index_t thread_id)
		{
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
			std::vector<UIPair> &check_pairs = check_pairs_mt[thread_id];
			check_pairs.clear();
#else
			size_t local_num_BBI_pairs = 0;
#endif

			if (ignore_same_label)
			{
				for (index_t bi0 = thread_id; bi0 < num_boxes; bi0 += num_threads)
				{
					const typename Tree::TreeBbox &b0 = cache_boxes[bi0];
					for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
					{
						if ((cache_labels[bi0] & cache_labels[bi1]).any())
							continue;

						const typename Tree::TreeBbox &b1 = cache_boxes[bi1];
						if (!DoIntersect()(b0, b1))
							continue; // early reject.

						OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
						check_pairs.push_back(uniquePair(b0.id(), b1.id()));
#else
						local_num_BBI_pairs++;
#endif
					}
				}
			}
			else
			{
				for (index_t bi0 = thread_id; bi0 < num_boxes; bi0 += num_threads)
				{
					const typename Tree::TreeBbox &b0 = cache_boxes[bi0];
					for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
					{
						const typename Tree::TreeBbox &b1 = cache_boxes[bi1];
						if (!DoIntersect()(b0, b1))
							continue; // early reject.

						OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
#ifdef OMC_ARR_COLLECT_BBI_PAIRS
						check_pairs.push_back(uniquePair(b0.id(), b1.id()));
#else
						local_num_BBI_pairs++;
#endif
					}
				}
			}
#ifndef OMC_ARR_COLLECT_BBI_PAIRS
			num_BBI_pairs += local_num_BBI_pairs;
#endif
		};
#ifdef OMC_ARR_DETECT_BBI_PARA
		tbb::parallel_for(size_t(0), num_threads, on_large_node);
#else
		std::ranges::for_each(std::ranges::iota_view(size_t(0), num_threads),
		                      on_large_node);
#endif

#ifdef OMC_ARR_COLLECT_BBI_PAIRS
		// Collect unique pairs
		OMC_UNUSED size_t new_size = std::accumulate(
		  check_pairs_mt.begin(), check_pairs_mt.end(), size_t(0),
		  [](size_t s, const std::vector<UIPair> &cp) { return s + cp.size(); });

		OMC_ARR_PROFILE_INC_REACH_CNT(ArrFuncNames::D_BBI_DUPL, 0, new_size);

		for (std::vector<UIPair> &cp : check_pairs_mt)
			BBI_pairs.insert(BBI_pairs.end(), cp.begin(), cp.end());
#endif
	}
}

template <typename Traits>
void DetectBBI<Traits>::cacheBoxesInNode(
  const typename Tree::Node             &node,
  CStyleVector<typename Tree::TreeBbox> &cached_boxes, bool cache_labels,
  CStyleVector<Label> &cached_labels)
{
	const auto  &boxes     = node.boxes();
	const size_t num_boxes = node.size();

	cached_boxes.resize(num_boxes, /*keep_data*/ false);
	for (size_t i = 0; i < num_boxes; i++)
		cached_boxes[i] = tree.box(boxes[i]);

	if (cache_labels)
	{
		cached_labels.resize(num_boxes, /*keep_data*/ false);
		for (size_t i = 0; i < num_boxes; i++)
			cached_labels[i] = labels[tree.box(boxes[i]).id()];
	}
}

} // namespace OMC