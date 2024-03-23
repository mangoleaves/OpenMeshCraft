#pragma once

#include "DetectIntersections.h"

namespace OMC {

template <typename Traits>
DetectIntersections<Traits>::DetectIntersections(
  const std::vector<GPoint *> &_verts, const std::vector<index_t> &_tris,
  const std::vector<Label> &_labels, const size_t _num_labels, Tree &_tree,
  std::vector<UIPair> &_intersection_list, std::vector<ShewchukCache> &_cache,
  bool _ignore_inter_with_same_label, MeshArrangements_Stats *_stats,
  bool _verbose)
  : verts(_verts)
  , tris(_tris)
  , labels(_labels)
  , num_labels(_num_labels)
  , tree(_tree)
  , intersection_list(_intersection_list)
  , cache(_cache)
  , ignore_inter_with_same_label(_ignore_inter_with_same_label)
  , stats(_stats)
  , verbose(_verbose)
{
	intersection_list.clear();
	intersection_list.reserve(_tris.size() / 3);

#ifndef OMC_ARR_DI_NO_CACHE
	cache  = std::vector<ShewchukCache>(tree.size());
	cached = std::vector<std::once_flag>(tree.size());
#endif

	std::vector<index_t> leaf_nodes = tree.all_leaf_nodes();

	std::vector<index_t> uniq_leaf_nodes, dupl_leaf_nodes;
	partitionLeafNodes(leaf_nodes, uniq_leaf_nodes, dupl_leaf_nodes);
	parallelOnLeafNodes(uniq_leaf_nodes);
	parallelOnUniqPairs(dupl_leaf_nodes);

	parallel_remove_duplicates(intersection_list);

	/* #region(collapsed) log */
	if (verbose)
		Logger::info(
		  std::format("[OpenMeshCraft Arrangements] Detect intersections "
		              "(only unique pairs): {} pairs.",
		              intersection_list.size()));
	/* #endregion */
}

template <typename Traits>
void DetectIntersections<Traits>::partitionLeafNodes(
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
void DetectIntersections<Traits>::parallelOnLeafNodes(
  const std::vector<index_t> &leaf_nodes)
{
	// Detect intersections
	tbb::spin_mutex mutex;

	// Test if two triangles intersect and cache intermediate results.
	auto test_tri_tri = [this, &mutex](index_t t0_id, const EPoint *T0_v[],
	                                   index_t t1_id, const EPoint *T1_v[])
	{
#ifndef OMC_ARR_DI_NO_CACHE
		std::call_once(cached[t0_id],
		               [this, &T0_v, &t0_id]()
		               {
			               Orient3D::get_minors(T0_v[0]->data(), T0_v[1]->data(),
			                                    T0_v[2]->data(), cache[t0_id].minor,
			                                    cache[t0_id].perm);
		               });
		std::call_once(cached[t1_id],
		               [this, &T1_v, &t1_id]()
		               {
			               Orient3D::get_minors(T1_v[0]->data(), T1_v[1]->data(),
			                                    T1_v[2]->data(), cache[t1_id].minor,
			                                    cache[t1_id].perm);
		               });
#endif
		if (intersectsTriangle(*T0_v[0], *T0_v[1], *T0_v[2], *T1_v[0], *T1_v[1],
		                       *T1_v[2]
#ifndef OMC_ARR_DI_NO_CACHE
		                       ,
		                       cache[t0_id].minor, cache[t0_id].perm,
		                       cache[t1_id].minor, cache[t1_id].perm
#endif
		                       ))
		{
			std::lock_guard<tbb::spin_mutex> guard(mutex);
			intersection_list.push_back(uniquePair(t0_id, t1_id));
		}
	};

	// Run intersection detection in a leaf node.
	auto on_leaf_node = [this, &test_tri_tri](size_t node_idx)
	{
		const size_t num_boxes = tree.node(node_idx).size();
		const auto  &boxes     = tree.node(node_idx).boxes();

		// for better memory cache (it really saves a lot of time :) )
		// OPT cache labels and boxes for threads.
		std::vector<Label>                 cache_labels;
		std::vector<typename Tree::OrBbox> cache_boxes;
		{
			cache_labels.reserve(num_boxes);
			cache_boxes.reserve(num_boxes);
			for (size_t i = 0; i < num_boxes; i++)
			{
				const typename Tree::OrBbox &b = tree.box(boxes[i]);
				cache_boxes.push_back(b);
				cache_labels.push_back(labels[b.id()]);
			}
		}

		if (ignore_inter_with_same_label)
		{
			// partition boxes by labels.
			// NOTE assume that labels are compact.
			std::vector<std::vector<index_t>> boxes_with_labels;
			{
				boxes_with_labels.resize(num_labels);
				for (std::vector<index_t> &bwl : boxes_with_labels)
					bwl.reserve(num_boxes / num_labels);
				// OPT avoid double for
				for (index_t bi = 0; bi < num_boxes; bi++)
					for (index_t li = 0; li < num_labels; li++)
						if (cache_labels[bi][li])
							boxes_with_labels[li].push_back(bi);
			}

			// OPT so much for
			for (index_t bwl_0 = 0; bwl_0 < boxes_with_labels.size(); bwl_0++)
			{
				for (index_t bwl_1 = bwl_0 + 1; bwl_1 < boxes_with_labels.size();
				     bwl_1++)
				{
					for (index_t bi0 : boxes_with_labels[bwl_0])
					{
						const typename Tree::OrBbox &b0 = cache_boxes[bi0];

						const EPoint *T0_v[3] = {&AsEP()(*verts[tris[b0.id() * 3 + 0]]),
						                         &AsEP()(*verts[tris[b0.id() * 3 + 1]]),
						                         &AsEP()(*verts[tris[b0.id() * 3 + 2]])};
						for (index_t bi1 : boxes_with_labels[bwl_1])
						{
							if (bi0 >= bi1)
								continue;

							const typename Tree::OrBbox &b1 = cache_boxes[bi1];

							if (!DoIntersect()(b0.bbox(), b1.bbox()))
								continue; // early reject

							const EPoint *T1_v[3] = {&AsEP()(*verts[tris[b1.id() * 3 + 0]]),
							                         &AsEP()(*verts[tris[b1.id() * 3 + 1]]),
							                         &AsEP()(*verts[tris[b1.id() * 3 + 2]])};

							test_tri_tri(b0.id(), T0_v, b1.id(), T1_v);
						}
					}
				}
			}
		}
		else
		{
			for (index_t bi0 = 0; bi0 < num_boxes; bi0++)
			{
				const typename Tree::OrBbox &b0 = cache_boxes[bi0];

				const EPoint *T0_v[3] = {&AsEP()(*verts[tris[b0.id() * 3 + 0]]),
				                         &AsEP()(*verts[tris[b0.id() * 3 + 1]]),
				                         &AsEP()(*verts[tris[b0.id() * 3 + 2]])};
				for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
				{
					const typename Tree::OrBbox &b1 = cache_boxes[bi1];

					if (!DoIntersect()(b0.bbox(), b1.bbox()))
						continue; // early reject

					const EPoint *T1_v[3] = {&AsEP()(*verts[tris[b1.id() * 3 + 0]]),
					                         &AsEP()(*verts[tris[b1.id() * 3 + 1]]),
					                         &AsEP()(*verts[tris[b1.id() * 3 + 2]])};

					test_tri_tri(b0.id(), T0_v, b1.id(), T1_v);
				}
			}
		}
	};

	tbb::parallel_for_each(leaf_nodes.begin(), leaf_nodes.end(), on_leaf_node);
}

template <typename Traits>
void DetectIntersections<Traits>::parallelOnUniqPairs(
  const std::vector<index_t> &leaf_nodes)
{
	for (index_t node_idx : leaf_nodes)
	{
		// sort boxes for faster traversal
		const auto  &boxes     = tree.node(node_idx).boxes();
		const size_t num_boxes = boxes.size();

		std::vector<Label>                 cache_labels;
		std::vector<typename Tree::OrBbox> cache_boxes;
		{ // for better memory cache :>
			cache_labels.reserve(num_boxes);
			cache_boxes.reserve(num_boxes);
			for (size_t i = 0; i < num_boxes; i++)
			{
				const typename Tree::OrBbox &b = tree.box(boxes[i]);
				cache_boxes.push_back(b);
				cache_labels.push_back(labels[b.id()]);
			}
		}

		std::vector<std::vector<UIPair>> concurrency_check_pairs(
		  tbb::this_task_arena::max_concurrency());

		// collect (possibly duplicate) pairs
		auto collect_pairs = [this, num_boxes, &cache_labels, &cache_boxes,
		                      &concurrency_check_pairs](index_t bi0)
		{
			index_t thread_id = tbb::this_task_arena::current_thread_index();
			const typename Tree::OrBbox &b0 = cache_boxes[bi0];
			for (index_t bi1 = bi0 + 1; bi1 < num_boxes; bi1++)
			{
				if ((ignore_inter_with_same_label &&
				     (cache_labels[bi0] & cache_labels[bi1]).any()))
					continue;

				const typename Tree::OrBbox &b1 = cache_boxes[bi1];

				if (!DoIntersect()(b0.bbox(), b1.bbox()))
					continue; // early reject

				concurrency_check_pairs[thread_id].push_back(
				  uniquePair(b0.id(), b1.id()));
			}
		};
		tbb::parallel_for(index_t(0), num_boxes, collect_pairs);

		// Collect unique pairs
		std::vector<UIPair> check_pairs;
		check_pairs.reserve(std::accumulate(
		  concurrency_check_pairs.begin(), concurrency_check_pairs.end(), size_t(0),
		  [](size_t s, const std::vector<UIPair> &cp) { return s + cp.size(); }));

		for (std::vector<UIPair> &cp : concurrency_check_pairs)
			check_pairs.insert(check_pairs.end(), cp.begin(), cp.end());

		// Detect intersections
		tbb::spin_mutex mutex;

		auto test_tri_tri = [this, &mutex](index_t t0_id, const EPoint *T0_v[],
		                                   index_t t1_id, const EPoint *T1_v[])
		{
#ifndef OMC_ARR_DI_NO_CACHE
			std::call_once(cached[t0_id],
			               [this, &T0_v, &t0_id]()
			               {
				               Orient3D::get_minors(T0_v[0]->data(), T0_v[1]->data(),
				                                    T0_v[2]->data(), cache[t0_id].minor,
				                                    cache[t0_id].perm);
			               });
			std::call_once(cached[t1_id],
			               [this, &T1_v, &t1_id]()
			               {
				               Orient3D::get_minors(T1_v[0]->data(), T1_v[1]->data(),
				                                    T1_v[2]->data(), cache[t1_id].minor,
				                                    cache[t1_id].perm);
			               });
#endif
			if (intersectsTriangle(*T0_v[0], *T0_v[1], *T0_v[2], *T1_v[0], *T1_v[1],
			                       *T1_v[2]
#ifndef OMC_ARR_DI_NO_CACHE
			                       ,
			                       cache[t0_id].minor, cache[t0_id].perm,
			                       cache[t1_id].minor, cache[t1_id].perm
#endif
			                       ))
			{
				std::lock_guard<tbb::spin_mutex> guard(mutex);
				intersection_list.push_back(uniquePair(t0_id, t1_id));
			}
		};

		auto on_check_pair = [this, &test_tri_tri](const UIPair &check_pair)
		{
			index_t t0 = check_pair.first, t1 = check_pair.second;

			const EPoint *T0_v[3] = {&AsEP()(*verts[tris[t0 * 3 + 0]]),
			                         &AsEP()(*verts[tris[t0 * 3 + 1]]),
			                         &AsEP()(*verts[tris[t0 * 3 + 2]])};
			const EPoint *T1_v[3] = {&AsEP()(*verts[tris[t1 * 3 + 0]]),
			                         &AsEP()(*verts[tris[t1 * 3 + 1]]),
			                         &AsEP()(*verts[tris[t1 * 3 + 2]])};

			test_tri_tri(t0, T0_v, t1, T1_v);
		};

		tbb::parallel_for_each(check_pairs.begin(), check_pairs.end(),
		                       on_check_pair);
	}
}

template <typename Traits>
bool DetectIntersections<Traits>::intersectsTriangle(
  const EPoint &t1_v0, const EPoint &t1_v1, const EPoint &t1_v2,
  const EPoint &t2_v0, const EPoint &t2_v1, const EPoint &t2_v2, NT *t1_min,
  NT *t1_perm, NT *t2_min, NT *t2_perm)
{
	auto res = Triangle3_Triangle3_DoIntersect().intersection_type(
	  t1_v0.data(), t1_v1.data(), t1_v2.data(), t2_v0.data(), t2_v1.data(),
	  t2_v2.data(), t1_min, t1_perm, t2_min, t2_perm);

	return (res > SimplexIntersectionType::SIMPLICIAL_COMPLEX);
}

} // namespace OMC