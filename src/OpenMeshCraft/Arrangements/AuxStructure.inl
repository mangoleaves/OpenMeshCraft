#pragma once

#include "AuxStructure.h"

namespace OMC {

template <typename Traits>
std::pair<std::atomic<index_t> *, bool>
AuxPointMap_BTree<Traits>::insert(const GPoint *pp, std::atomic<index_t> *idx)
{
	auto ins = map.insert({AuxPoint<Traits>(pp), idx});
	return std::make_pair(
	  // the position of v (pos if first time, or the existed position otherwise)
	  ins.first->second,
	  // the result of the insert operation (true or false)
	  ins.second);
}

template <typename Traits>
std::pair<std::atomic<index_t> *, bool>
AuxPointMap_ConcurrentMap<Traits>::insert(const GPoint         *pp,
                                          std::atomic<index_t> *idx)
{
	auto ins = map.insert({AuxPoint<Traits>(pp), idx});
	return std::make_pair(
	  // the position of v (pos if first time, or the existed position otherwise)
	  ins.first->second,
	  // the result of the insert operation (true or false)
	  ins.second);
}

template <typename Traits>
std::pair<std::atomic<index_t> *, bool>
AuxPointMap_Tree<Traits>::insert(const GPoint *pp, std::atomic<index_t> *idx)
{
	return tree->insert_point(pp, idx);
}

template <typename Traits>
template <typename GetIndex>
std::pair<index_t, bool>
AuxPointMap_Tree<Traits>::insert_F(const GPoint *pp, std::atomic<index_t> *idx,
                                   GetIndex get_idx)
{
	return tree->insert_point_F(pp, idx, get_idx);
}

} // namespace OMC