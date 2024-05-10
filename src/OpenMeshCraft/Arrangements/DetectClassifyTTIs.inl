#pragma once

#include "DetectClassifyTTI.h"
#include "DetectClassifyTTIs.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"
#include "OpenMeshCraft/Utils/Exception.h"


#if defined(OMC_ARR_PROFILE) && 0
	#define COLLECT_INTERSECTING_TRIANGLE                                     \
		if (Tri3_Tri3_DoInter().intersection_type(                              \
		      ts.triVertPtr(b0.id(), 0), ts.triVertPtr(b0.id(), 1),             \
		      ts.triVertPtr(b0.id(), 2), ts.triVertPtr(b1.id(), 0),             \
		      ts.triVertPtr(b1.id(), 1),                                        \
		      ts.triVertPtr(b1.id(), 2)) >= SimplexIntersectionType::INTERSECT) \
		{                                                                       \
			intersecting_triangle_pairs.push_back(uniquePair(b0.id(), b1.id()));  \
		}

	#define REPORT_INTERSECTING_TRIANGLE                                \
		remove_duplicates(intersecting_triangle_pairs);                   \
		OMC::Logger::info(std::format(                                    \
		  "[OpenMeshCraft Arrangements] intersecting triangle pairs: {}", \
		  intersecting_triangle_pairs.size()));
#else
	#define COLLECT_INTERSECTING_TRIANGLE
	#define REPORT_INTERSECTING_TRIANGLE
#endif

namespace OMC {

template <typename Traits>
DetectClassifyTTIs<Traits>::DetectClassifyTTIs(TriSoup &_ts, const Tree &_tree,
                                               bool _ignore_same_label,
                                               MeshArrangements_Stats *_stats,
                                               bool                    _verbose)
  : ts(_ts)
  , pnt_arenas(*ts.pnt_arenas)
  , labels(ts.tri_labels)
  , tree(_tree)
  , ignore_same_label(_ignore_same_label)
  , stats(_stats)
  , verbose(_verbose)
{
	/* Detect and classify Triangle-Triangle-Intersection */

	// calculate orthogonal planes
	ts.calcOrthogonalPlane();

	// check Triangle-Triangle-Intersection
	std::vector<index_t> leaf_nodes = tree.all_leaf_nodes();
	std::vector<index_t> small_nodes, large_nodes;
	partitionNodes(leaf_nodes, small_nodes, large_nodes);

	GPoint::enable_global_cached_values(tbb::this_task_arena::max_concurrency());
	parallelOnSmallNodes(small_nodes);
	parallelOnLargeNodes(large_nodes);
	GPoint::disable_global_cached_values();

	REPORT_INTERSECTING_TRIANGLE;
	/* Post fix of Triangle-Triangle-Intersection */

	// add end points to all ts.edge2pts
	ts.addEndPointsToE2P();
	// remove duplicates during parallel execution
	ts.removeDuplicatesBeforeFix();
#ifndef OMC_ARR_GLOBAL_POINT_SET
	// post fix of indices
	ts.fixColinearEdgesIntersections();
#endif
	ts.fixAllIndices();
	// calculate orthogonal planes and orientations
	ts.calcTriangleOrient();

	/* propagate intersection points from edges to coplanar triangles */

	propagateCoplanarIntersections();
}

template <typename Traits>
void DetectClassifyTTIs<Traits>::partitionNodes(
  const std::vector<index_t> &nodes, std::vector<index_t> &small_nodes,
  std::vector<index_t> &large_nodes)
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
void DetectClassifyTTIs<Traits>::parallelOnSmallNodes(
  const std::vector<index_t> &nodes)
{
	// Run intersection detection in a leaf node.
	auto on_small_node = [&](size_t node_idx)
	{
		const typename Tree::Node &node = tree.node(node_idx);

		index_t thread_id = tbb::this_task_arena::current_thread_index();

		CStyleVector<typename Tree::TreeBbox> cache_boxes;
		CStyleVector<Label>                   cache_labels;

		cacheBoxesInNode(node, cache_boxes, ignore_same_label, cache_labels);

		const size_t num_boxes = tree.node(node_idx).size();

		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
		                              num_boxes * (num_boxes - 1) / 2);
		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI_SMALL,
		                              num_boxes * (num_boxes - 1) / 2);

		GPoint::clear_global_cached_values();
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

					DetectClassifyTTI<Traits> dc(ts, pnt_arenas[thread_id]);
					dc.check_TTI(b0.id(), b1.id());

					COLLECT_INTERSECTING_TRIANGLE;
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

					DetectClassifyTTI<Traits> dc(ts, pnt_arenas[thread_id]);
					dc.check_TTI(b0.id(), b1.id());

					COLLECT_INTERSECTING_TRIANGLE;
				}
			}
		}
	};

#ifdef OMC_ARR_DC_TTI_PARA
	tbb::parallel_for_each(nodes.begin(), nodes.end(), on_small_node);
#else
	std::for_each(nodes.begin(), nodes.end(), on_small_node);
#endif
}

template <typename Traits>
void DetectClassifyTTIs<Traits>::parallelOnLargeNodes(
  const std::vector<index_t> &nodes)
{
	size_t num_threads = tbb::this_task_arena::max_concurrency();

	for (index_t node_idx : nodes)
	{
		const typename Tree::Node &node = tree.node(node_idx);

		size_t num_boxes = node.size();

		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI,
		                              num_boxes * (num_boxes - 1) / 2);
		OMC_ARR_PROFILE_INC_TOTAL_CNT(ArrFuncNames::D_BBI_LARGE,
		                              num_boxes * (num_boxes - 1) / 2);

		CStyleVector<typename Tree::TreeBbox> cache_boxes;
		CStyleVector<Label>                   cache_labels;
		cacheBoxesInNode(node, cache_boxes, ignore_same_label, cache_labels);

		// check Box-Box-Intersection in current node.
		auto on_large_node = [&](index_t thread_id)
		{
			GPoint::clear_global_cached_values();
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
						DetectClassifyTTI<Traits> dc(ts, pnt_arenas[thread_id]);
						dc.check_TTI(b0.id(), b1.id());

						COLLECT_INTERSECTING_TRIANGLE;
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
						DetectClassifyTTI<Traits> dc(ts, pnt_arenas[thread_id]);
						dc.check_TTI(b0.id(), b1.id());

						COLLECT_INTERSECTING_TRIANGLE;
					}
				}
			}
		};
#ifdef OMC_ARR_DC_TTI_PARA
		tbb::parallel_for(size_t(0), num_threads, on_large_node);
#else
		std::ranges::for_each(std::ranges::iota_view(size_t(0), num_threads),
		                      on_large_node);
#endif
	}
}

template <typename Traits>
void DetectClassifyTTIs<Traits>::cacheBoxesInNode(
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

template <typename Traits>
void DetectClassifyTTIs<Traits>::propagateCoplanarIntersections()
{
	// mutexes on tri
	std::vector<tbb::spin_mutex> tri_mutexes(ts.numTris());

	auto propagate_copl_edge = [this, &tri_mutexes](index_t t_id)
	{
		const concurrent_vector<typename TriSoup::CCrEdgeInfo> &ccr_edge_infos =
		  ts.coplanarEdges(t_id);

		for (const typename TriSoup::CCrEdgeInfo &edge_info : ccr_edge_infos)
		{
			index_t copl_e_id = edge_info.e_id;
			index_t v0_id     = edge_info.v0_id;
			index_t v1_id     = edge_info.v1_id;

			// intersection points inside triangle
			bool inside = false;

			const concurrent_vector<index_t>     &t2p = ts.trianglePointsList(t_id);
			const typename TriSoup::Edge2PntsSet &e2p = ts.edgePointsList(copl_e_id);

			for (index_t p_id : e2p)
			{
				if (!inside)
				{
					if (p_id == v0_id || p_id == v1_id) // step from outside to inside
						inside = true;
				}
				else // already inside
				{
					if (p_id == v0_id || p_id == v1_id) // step from inside to outside
						break;                            // end traversal current edge.
					// else, current p_id is inside triangle
					// check if point is already in t2p
					if (std::find(t2p.begin(), t2p.end(), p_id) == t2p.end())
						// if not, add it
						ts.addVertexInTriangle(t_id, p_id);
				}
			}
		}
	};

#ifdef OMC_ARR_DC_TTI_PARA
	tbb::parallel_for(size_t(0), ts.numTris(), propagate_copl_edge);
#else
	for (size_t t_id = 0; t_id < ts.numTris(); t_id++)
		propagate_copl_edge(t_id);
#endif
}

} // namespace OMC

#undef COLLECT_INTERSECTING_TRIANGLE
#undef REPORT_INTERSECTING_TRIANGLE