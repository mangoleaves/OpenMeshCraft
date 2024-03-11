#pragma once

#include "Utils.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace OMC {

template <typename Traits>
class ConnectIntersections
{
public:
	ConnectIntersections(
	  const std::vector<index_t> &_tris, std::vector<UIPair> &_intersection_list,
	  std::vector<std::vector<UIPair>> &_connected_intersection_lists,
	  MeshArrangements_Stats *_stats = nullptr, bool _verbose = false);

private:
	void constructConnectedComponent();

	void balenceConnectedComponent();

private:
	using Graph =
	  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;

	size_t expect_num_pairs = 1000;

	/* Input data */
	const std::vector<index_t>       &tris;
	std::vector<UIPair>              &intersection_list;
	/* Output data */
	std::vector<std::vector<UIPair>> &connected_intersection_lists;
	/* statistics */
	MeshArrangements_Stats           *stats;
	/* Behavior control flags */
	bool                              verbose;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ConnectIntersections.inl"
#endif