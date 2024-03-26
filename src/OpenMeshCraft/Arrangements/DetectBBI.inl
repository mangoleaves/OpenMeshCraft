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
	parallelOnLeafNodes(uniq_leaf_nodes);
	parallelOnUniqPairs(dupl_leaf_nodes);

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
		if (tree.node(ni).size() > 1000)
			dupl_leaf_nodes.push_back(ni);
		else
			uniq_leaf_nodes.push_back(ni);
	}
}

template <typename Traits>
void DetectBBI<Traits>::parallelOnLeafNodes(
  const std::vector<index_t> &leaf_nodes)
{
	std::vector<std::vector<Label>> cache_labels_mt(
	  tbb::this_task_arena::max_concurrency());
	std::vector<std::vector<typename Tree::OrBbox>> cache_boxes_mt(
	  tbb::this_task_arena::max_concurrency());
	std::vector<std::vector<UIPair>> check_pairs_mt(
	  tbb::this_task_arena::max_concurrency());

	// Run intersection detection in a leaf node.
	auto on_leaf_node =
	  [this, &cache_labels_mt, &cache_boxes_mt, &check_pairs_mt](size_t node_idx)
	{
		int thread_id = tbb::this_task_arena::current_thread_index();

		std::vector<Label> &cache_labels = cache_labels_mt[thread_id];
		std::vector<typename Tree::OrBbox> &cache_boxes = cache_boxes_mt[thread_id];
		std::vector<UIPair>                &check_pairs = check_pairs_mt[thread_id];

		const size_t num_boxes = tree.node(node_idx).size();
		const auto  &boxes     = tree.node(node_idx).boxes();

		// for better memory cache (it really saves a lot of time :) )
		cache_labels.resize(num_boxes);
		cache_boxes.resize(num_boxes);
		for (size_t i = 0; i < num_boxes; i++)
		{
			const typename Tree::OrBbox &b = tree.box(boxes[i]);

			cache_boxes[i]  = b;
			cache_labels[i] = labels[b.id()];
		}

		if (ignore_same_label)
		{
			// partition boxes by labels.
			// NOTE assume that labels are compact.
			std::vector<std::vector<index_t>> boxes_with_labels;
			{
				boxes_with_labels.resize(num_labels);
				for (std::vector<index_t> &bwl : boxes_with_labels)
					bwl.reserve(num_boxes / num_labels);
				for (index_t bi = 0; bi < num_boxes; bi++)
					for (index_t li = 0; li < num_labels; li++)
						if (cache_labels[bi][li])
							boxes_with_labels[li].push_back(bi);
			}

			for (index_t bwl_0 = 0; bwl_0 < boxes_with_labels.size(); bwl_0++)
			{
				for (index_t bwl_1 = bwl_0 + 1; bwl_1 < boxes_with_labels.size();
				     bwl_1++)
				{
					for (index_t bi0 : boxes_with_labels[bwl_0])
					{
						const typename Tree::OrBbox &b0 = cache_boxes[bi0];

						for (index_t bi1 : boxes_with_labels[bwl_1])
						{
							if (bi0 >= bi1)
								continue;

							const typename Tree::OrBbox &b1 = cache_boxes[bi1];

							if (!DoIntersect()(b0.bbox(), b1.bbox()))
								continue; // early reject

							OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
							check_pairs.push_back(uniquePair(b0.id(), b1.id()));
						}
					}
				}
			}
		}
		else
		{
			OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
			                              num_boxes * (num_boxes - 1) / 2);
			for (index_t bi0 = 0; bi0 < num_boxes; bi0++)
			{
				const typename Tree::OrBbox &b0 = cache_boxes[bi0];

				for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
				{
					const typename Tree::OrBbox &b1 = cache_boxes[bi1];

					if (!DoIntersect()(b0.bbox(), b1.bbox()))
						continue; // early reject

					OMC_EXPENSIVE_ASSERT(b0.id() != b1.id(), "duplicate triangles.");
					check_pairs.push_back(uniquePair(b0.id(), b1.id()));
				}
			}
		}
	};

	tbb::parallel_for_each(leaf_nodes.begin(), leaf_nodes.end(), on_leaf_node);

	// Collect unique pairs
	BBI_pairs.reserve(BBI_pairs.size() +
	                  std::accumulate(check_pairs_mt.begin(),
	                                  check_pairs_mt.end(), size_t(0),
	                                  [](size_t s, const std::vector<UIPair> &cp)
	                                  { return s + cp.size(); }));

	for (std::vector<UIPair> &cp : check_pairs_mt)
		BBI_pairs.insert(BBI_pairs.end(), cp.begin(), cp.end());
}

template <typename Traits>
void DetectBBI<Traits>::parallelOnUniqPairs(
  const std::vector<index_t> &leaf_nodes)
{
	std::vector<Label>                 cache_labels;
	std::vector<typename Tree::OrBbox> cache_boxes;

	for (index_t node_idx : leaf_nodes)
	{
		// sort boxes for faster traversal
		const auto  &boxes     = tree.node(node_idx).boxes();
		const size_t num_boxes = boxes.size();

		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
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
		}

		std::vector<std::vector<UIPair>> check_pairs_mt(
		  tbb::this_task_arena::max_concurrency());

		// collect (possibly duplicate) pairs
		auto collect_pairs = [this, num_boxes, &cache_labels, &cache_boxes,
		                      &check_pairs_mt](index_t bi0)
		{
			index_t thread_id = tbb::this_task_arena::current_thread_index();
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
				check_pairs_mt[thread_id].push_back(uniquePair(b0.id(), b1.id()));
			}
		};
		tbb::parallel_for(index_t(0), num_boxes, collect_pairs);

		// Collect unique pairs
		size_t new_size =
		  BBI_pairs.size() +
		  std::accumulate(check_pairs_mt.begin(), check_pairs_mt.end(), size_t(0),
		                  [](size_t s, const std::vector<UIPair> &cp)
		                  { return s + cp.size(); });

		if (new_size > BBI_pairs.capacity())
			BBI_pairs.reserve(new_size * 2);

		for (std::vector<UIPair> &cp : check_pairs_mt)
			BBI_pairs.insert(BBI_pairs.end(), cp.begin(), cp.end());
	}
}

} // namespace OMC