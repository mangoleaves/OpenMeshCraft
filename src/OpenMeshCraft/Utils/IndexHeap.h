#pragma once

#include "IndexDef.h"

#include "Exception.h"

#include <vector>

namespace OMC {

/**
 * @brief Each element is a pair of key and value (k, v).
 * Key is supposed to be an index with type index_t.
 * Value is usually something to be sorted, e.g., priority.
 * IndexDenseHeap use a heap to sort values by comparing value.
 * It also can modify value by key and resort the values.
 * @note We suppose that key is densely distributed, otherwise
 * use unordered_map to implement m_pos_in_heap.
 */
template <typename value_type, typename comparator = std::less<value_type>>
class IndexDenseHeap
{
public:
	using key_type = index_t;
	using kv_pair  = std::pair<key_type, value_type>;

public:
	/// @brief Default constructor
	IndexDenseHeap() { m_heap.resize(min_size()); }

	constexpr size_t min_size() { return 1; }
	size_t           heap_size() const { return m_heap.size() - 1; }

	void reserve(size_t s)
	{
		m_heap.reserve(s + 1);
		m_pos_in_heap.reserve(s);
	}

	bool empty() const { return heap_size() == 0; }

	/**
	 * @brief push a (key, value) pair into heap.
	 */
	template <bool AllowUpdate = true>
	void push(key_type key, const value_type &value)
	{
		if constexpr (AllowUpdate)
		{
			if (exist(key))
			{
				update(key, value);
				return;
			}
		}
		else
		{
			OMC_EXPENSIVE_ASSERT(!exist(key), "Add duplicate key into heap.");
		}

		m_heap.emplace_back(key, value);
		if (key < m_pos_in_heap.size())
		{
			m_pos_in_heap[key] = heap_size(); // start from 1, not 0.
		}
		else
		{
			m_pos_in_heap.resize(key + 1, InvalidIndex);
			m_pos_in_heap[key] = heap_size();
		}
		up_move(heap_size());
	}

	/**
	 * @brief get the (key, value) pair on the top of heap.
	 */
	const kv_pair &top() const
	{
		OMC_EXPENSIVE_ASSERT(!empty(), "Heap is empty.");
		return m_heap[first_pos()];
	}

	/**
	 * @brief pop the (key, value) pair on the top of heap.
	 */
	void pop()
	{
		if (empty())
			return;
		// get element in the first position, move the last element to the first
		// position.
		m_pos_in_heap[m_heap[first_pos()].first] = InvalidIndex;
		if (heap_size() > 1)
		{
			m_heap[first_pos()] = std::move(m_heap[last_pos()]);
			m_heap.pop_back();
			// down move the the element in first position.
			down_move(first_pos());
		}
		else
		{
			m_heap.pop_back();
		}
	}

	bool exist(key_type key) const
	{
		return is_valid_idx(key) && key < m_pos_in_heap.size() &&
		       is_valid_idx(m_pos_in_heap[key]);
	}

	/**
	 * @brief update a value specified by key.
	 * @param value the new value.
	 */
	void update(key_type key, const value_type &value)
	{
		OMC_EXPENSIVE_ASSERT(exist(key),
		                     "Key does not exist in heap, fail to update.");
		// find the pos in heap, only update value, but won't update key
		index_t           pos       = m_pos_in_heap[key];
		const value_type &old_value = m_heap[pos].second;

		bool move_up_or_down = comparator()(value, old_value);
		m_heap[pos].second   = value;

		if (move_up_or_down)
			up_move(pos);
		else
			down_move(pos);
	}

	/**
	 * @brief remove a pair of key and value specified by key.
	 */
	void remove(key_type key)
	{
		// similar to update, find pos in heap, but update key and value at the same
		// time. the new key and value are from last position.
		index_t pos        = m_pos_in_heap[key];
		m_pos_in_heap[key] = InvalidIndex;

		if (pos == last_pos())
		{
			m_heap.pop_back();
			return;
		}

		const value_type &old_value  = m_heap[pos].second;
		const value_type &last_value = m_heap[last_pos()].second;

		bool move_up_or_down             = comparator()(last_value, old_value);
		m_heap[pos]                      = std::move(m_heap[last_pos()]);
		m_pos_in_heap[m_heap[pos].first] = pos;
		m_heap.pop_back();

		// move kv pair on the last position to the position to remove.
		if (move_up_or_down)
			up_move(pos);
		else
			down_move(pos);
	}

protected:
	constexpr size_t empty_pos() const { return 0; }
	constexpr size_t first_pos() const { return 1; }
	size_t           last_pos() const { return m_heap.size() - 1; }

	index_t leftc(index_t i) const { return 2 * i; }
	index_t rightc(index_t i) const { return 2 * i + 1; }
	index_t parent(index_t i) const { return i >> 1; }

	/**
	 * @brief Compare two kv_pair by value in heap indexed by i0 and i1.
	 */
	bool compare(index_t i0, index_t i1) const
	{
		return comparator()(m_heap[i0].second, m_heap[i1].second);
	}

	/**
	 * @brief move element in heap[idx] up.
	 * @param idx the position in heap.
	 */
	void up_move(index_t idx)
	{
		m_heap[empty_pos()] = std::move(m_heap[idx]);
		index_t p           = parent(idx);

		while (p >= first_pos() && compare(empty_pos(), p))
		{
			m_heap[idx]                      = std::move(m_heap[p]);
			m_pos_in_heap[m_heap[idx].first] = idx;

			idx = p;
			p   = parent(idx);
		}

		m_heap[idx]                      = std::move(m_heap[empty_pos()]);
		m_pos_in_heap[m_heap[idx].first] = idx;
	}

	/**
	 * @brief move element in heap[idx] down.
	 * @param idx the position in heap.
	 */
	void down_move(index_t idx)
	{
		m_heap[empty_pos()] = std::move(m_heap[idx]);

		while (idx <= (heap_size() >> 1))
		{
			index_t lc = leftc(idx);  // left child must exist
			index_t rc = rightc(idx); // right child may not exist

			index_t min_idx = idx; // index in (idx, lc, rc) with minimal value

			if (compare(lc, empty_pos()))
			{
				min_idx = rc <= heap_size() && compare(rc, lc) ? rc : lc;
			}
			else if (rc <= heap_size() && compare(rc, empty_pos()))
			{
				min_idx = rc;
			}

			if (min_idx != idx)
			{
				m_heap[idx]                      = std::move(m_heap[min_idx]);
				m_pos_in_heap[m_heap[idx].first] = idx;
				idx                              = min_idx;
			}
			else
				break;
		}
		m_heap[idx]                      = std::move(m_heap[empty_pos()]);
		m_pos_in_heap[m_heap[idx].first] = idx;
	}

protected:
	std::vector<key_type> m_pos_in_heap; /// < map key to position in heap
	std::vector<kv_pair>  m_heap;        /// < sort (key, value).
};

} // namespace OMC