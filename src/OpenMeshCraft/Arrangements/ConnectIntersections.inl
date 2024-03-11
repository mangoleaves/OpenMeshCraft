#pragma once

#include "ConnectIntersections.h"

namespace OMC {

template <typename Traits>
ConnectIntersections<Traits>::ConnectIntersections(
  const std::vector<index_t> &_tris, std::vector<UIPair> &_intersection_list,
  std::vector<std::vector<UIPair>> &_connected_intersection_lists,
  MeshArrangements_Stats *_stats, bool _verbose)
  : tris(_tris)
  , intersection_list(_intersection_list)
  , connected_intersection_lists(_connected_intersection_lists)
  , stats(_stats)
  , verbose(_verbose)
{
	connected_intersection_lists.clear();

	if (intersection_list.empty())
		return;

	if (intersection_list.size() < 1000)
	{
		/* #region(collapsed) log */
		if (verbose)
			Logger::info("[OpenMeshCraft Arrangements] Intersection pairs are too "
			             "little, so merge them all.");
		/* #endregion */
		connected_intersection_lists.push_back(intersection_list);
		return;
	}

	constructConnectedComponent();

	balenceConnectedComponent();
}

template <typename Traits>
void ConnectIntersections<Traits>::constructConnectedComponent()
{
	// Only build graph on triangles that have intersections.
	// So, first we collect triangles with intersections.
	size_t               num_tris = tris.size() / 3;
	std::vector<uint8_t> tri_has_intersection(num_tris, false);

	for (const UIPair &uip : intersection_list)
	{
		tri_has_intersection[uip.first]  = true;
		tri_has_intersection[uip.second] = true;
	}

	std::vector<index_t> map_idx(num_tris, InvalidIndex);
	index_t              map_t_id = 0;

	for (index_t t_id = 0; t_id < num_tris; t_id++)
		if (tri_has_intersection[t_id])
			map_idx[t_id] = map_t_id++;

	Graph graph;
	// Then, collect connected edges in graph.
	// A pair of intersected triangles are seen as connected.
	for (const UIPair &uip : intersection_list)
		boost::add_edge(map_idx[uip.first], map_idx[uip.second], graph);

	// construct connected components
	std::vector<size_t> component(boost::num_vertices(graph));
	index_t num_components = boost::connected_components(graph, component.data());

	// collect connected intersection list
	connected_intersection_lists.resize(num_components);
	for (const UIPair &uip : intersection_list)
	{
		OMC_EXPENSIVE_ASSERT(
		  component[map_idx[uip.first]] == component[map_idx[uip.second]],
		  "intersection pair does not appear in one connected component");

		connected_intersection_lists[component[map_idx[uip.first]]].push_back(uip);
	}
}

template <typename Traits>
void ConnectIntersections<Traits>::balenceConnectedComponent()
{
	if (connected_intersection_lists.empty())
		return;

	size_t total_pairs = intersection_list.size();

	tbb::parallel_sort(
	  connected_intersection_lists.begin(), connected_intersection_lists.end(),
	  [](const std::vector<UIPair> &lhs, const std::vector<UIPair> &rhs)
	  { return lhs.size() > rhs.size(); });

	if (connected_intersection_lists.front().size() > total_pairs / 2)
	{
		// first connected component is too large,
		// connected components are imbalence, merge all of them.
		connected_intersection_lists.clear();
		connected_intersection_lists.push_back(intersection_list);
		connected_intersection_lists.resize(1);
		/* #region(collapsed) log */
		if (verbose)
			Logger::info("[OpenMeshCraft Arrangements] Connected component are "
			             "imbalence, so merge them all.");
		/* #endregion */
		return;
	}

	// merge too little connected component.
	// too much connected components, merge and balence them.
	expect_num_pairs   = std::max(size_t(1000), total_pairs / 32);
	size_t dst_list_id = 0;
	while (dst_list_id < connected_intersection_lists.size() - 1)
	{
		while (dst_list_id < connected_intersection_lists.size() - 1 &&
		       connected_intersection_lists[dst_list_id].size() < expect_num_pairs)
		{
			connected_intersection_lists[dst_list_id].insert(
			  connected_intersection_lists[dst_list_id].end(),
			  connected_intersection_lists.back().begin(),
			  connected_intersection_lists.back().end());
			connected_intersection_lists.pop_back();
		}
		++dst_list_id;
	}
	/* #region(collapsed) log */
	if (verbose)
	{
		std::string msg = "";
		for (const std::vector<UIPair> &inter_list : connected_intersection_lists)
			msg += std::to_string(inter_list.size()) + ", ";
		Logger::info(std::format(
		  "[OpenMeshCraft Arrangements] Balence connected componnets: {}", msg));
	}
	/* #endregion */
}

} // namespace OMC