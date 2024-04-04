#pragma once

#include "Utils.h"

// clang-format off
#include "OpenMeshCraft/Utils/DisableWarnings.h"
#include "parallel_hashmap/btree.h"
#include "OpenMeshCraft/Utils/EnableWarnings.h"
// clang-format on

#include <mutex>
#include <set>

namespace OMC {

/// @brief Wrap of GPoint*
template <typename Traits>
class AuxPoint
{
public:
	using GPoint     = typename Traits::GPoint;
	using LessThan3D = typename Traits::LessThan3D;

	const GPoint *pt;
	AuxPoint(const GPoint *pt_)
	  : pt(pt_)
	{
	}

	bool operator<(const AuxPoint &rhs) const
	{
		return LessThan3D()(*pt, *(rhs.pt)) == Sign::NEGATIVE;
	}
	bool operator==(const AuxPoint &rhs) const
	{
		return LessThan3D().coincident(*pt, *(rhs.pt));
	}
};

/// @brief Map GPoint to a unique index
template <typename Traits>
class AuxPointMap
{
public:
	using GPoint = typename Traits::GPoint;

public:
	AuxPointMap() = default;
	virtual ~AuxPointMap() {}

public:
	/// @brief Insert point into map.
	/// @param item pair of geometry implementation and input index
	/// @return pair of index (input index if the point is inserted in map first
	/// time, otherwise the previous saved index) and the result of the insert
	/// operation (true if succeed to insert otherwise false).
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) = 0;
};

template <typename Traits>
class AuxPointMap_BTree : public AuxPointMap<Traits>
{
public:
	using GPoint = typename Traits::GPoint;

public:
	AuxPointMap_BTree() = default;
	virtual ~AuxPointMap_BTree() {}

public:
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) override final;

protected:
	phmap::btree_map<AuxPoint<Traits>, std::atomic<index_t> *> map;
};

template <typename Traits>
class AuxPointMap_ConcurrentMap : public AuxPointMap<Traits>
{
public:
	using GPoint = typename Traits::GPoint;

public:
	AuxPointMap_ConcurrentMap() = default;
	virtual ~AuxPointMap_ConcurrentMap() {}

public:
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) override final;

protected:
	tbb::concurrent_map<AuxPoint<Traits>, std::atomic<index_t> *> map;
};

// Forward declaration
template <typename Traits>
class Arr_Tree_Intersection;

template <typename Traits>
class AuxPointMap_Tree : public AuxPointMap<Traits>
{
public:
	using GPoint = typename Traits::GPoint;
	using Tree   = Arr_Tree_Intersection<Traits>;

public:
	AuxPointMap_Tree(Tree *_tree)
	  : tree(_tree)
	{
	}
	virtual ~AuxPointMap_Tree() {}

public:
	virtual std::pair<std::atomic<index_t> *, bool>
	insert(const GPoint *pp, std::atomic<index_t> *idx) override final;

	template <typename GetIndex>
	std::pair<index_t, bool> insert_F(const GPoint *pp, std::atomic<index_t> *idx,
	                                  GetIndex get_idx);

private:
	Tree *tree;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "AuxStructure.inl"
#endif