#pragma once

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/Macros.h"

#include "parallel_hashmap/phmap.h"
#include "tbb/tbb.h"

#include <atomic>
#include <deque>
#include <memory>
#include <thread>
#include <unordered_map>

namespace OMC {

template <typename IT, typename ET>
class OnePointCachedValues2
{
public:
	using FT = double;

public:
	// floating-point filter
	// NaN demoninator is used to check if it is cached.
	FT ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_denominator,
	  ssfilter_max_val;
	// dynamic filter (interval number)
	// NaN demoninator is used to check if it is cached.
	IT  dfilter_lambda_x, dfilter_lambda_y, dfilter_denominator;
	// exact number
	// exact cached is used to check if it is cached.
	ET *exact_lambda_x = nullptr, *exact_lambda_y = nullptr,
	   *exact_denominator  = nullptr;
	// expansion number
	// nullptr denominator is used to check if it is cached.
	FT *expansion_lambda_x = nullptr, *expansion_lambda_y = nullptr,
	   *expansion_denominator  = nullptr;
	int expansion_lambda_x_len = 0, expansion_lambda_y_len = 0,
	    expansion_d_len = 0;

	bool ssfilter_cached  = false;
	bool dfilter_cached   = false;
	bool exact_cached     = false;
	bool expansion_cached = false;

	OnePointCachedValues2() noexcept
	{
		expansion_lambda_x     = nullptr;
		expansion_lambda_y     = nullptr;
		expansion_denominator  = nullptr;
		expansion_lambda_x_len = 0;
		expansion_lambda_y_len = 0;
		expansion_d_len        = 0;

		ssfilter_cached  = false;
		dfilter_cached   = false;
		exact_cached     = false;
		expansion_cached = false;
	}

	OnePointCachedValues2(const OnePointCachedValues2 &rhs) = delete;
	OnePointCachedValues2(OnePointCachedValues2 &&rhs)      = delete;

	void operator=(const OnePointCachedValues2 &rhs) = delete;
	void operator=(OnePointCachedValues2 &&rhs)      = delete;

	~OnePointCachedValues2()
	{
		if (expansion_cached)
		{
			free(expansion_lambda_x);
			free(expansion_lambda_y);
			free(expansion_denominator);
		}
		if (exact_cached)
		{
			delete exact_lambda_x;
			delete exact_lambda_y;
			delete exact_denominator;
		}
	}
};

#if defined(OMC_INDIRECT_PRED)

template <typename IT, typename ET>
class OnePointCachedValues3
{
public:
	using FT = double;

public:
	// floating-point filter
	// NaN demoninator is used to check if it is cached.
	#ifdef OMC_CACHE_SSF
	FT ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_lambda_z,
	  ssfilter_denominator, ssfilter_max_val;
	#endif
	// dynamic filter (interval number)
	// NaN demoninator is used to check if it is cached.
	IT  dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator;
	// exact number
	// exact cached is used to check if it is cached.
	ET *exact_lambda_x = nullptr, *exact_lambda_y = nullptr,
	   *exact_lambda_z = nullptr, *exact_denominator = nullptr;
	// expansion number
	// nullptr denominator is used to check if it is cached.
	FT *expansion_lambda_x = nullptr, *expansion_lambda_y = nullptr,
	   *expansion_lambda_z = nullptr, *expansion_denominator = nullptr;
	int expansion_lambda_x_len = 0, expansion_lambda_y_len = 0,
	    expansion_lambda_z_len = 0, expansion_d_len = 0;

	bool ssfilter_cached  = false;
	bool dfilter_cached   = false;
	bool exact_cached     = false;
	bool expansion_cached = false;

	OnePointCachedValues3() noexcept
	{
		exact_lambda_x    = nullptr;
		exact_lambda_y    = nullptr;
		exact_lambda_z    = nullptr;
		exact_denominator = nullptr;

		expansion_lambda_x     = nullptr;
		expansion_lambda_y     = nullptr;
		expansion_lambda_z     = nullptr;
		expansion_denominator  = nullptr;
		expansion_lambda_x_len = 0;
		expansion_lambda_y_len = 0;
		expansion_lambda_z_len = 0;
		expansion_d_len        = 0;

		ssfilter_cached  = false;
		dfilter_cached   = false;
		exact_cached     = false;
		expansion_cached = false;
	}

	OnePointCachedValues3(const OnePointCachedValues3 &rhs) = delete;
	OnePointCachedValues3(OnePointCachedValues3 &&rhs)      = delete;

	void operator=(const OnePointCachedValues3 &rhs) = delete;
	void operator=(OnePointCachedValues3 &&rhs)      = delete;

	~OnePointCachedValues3() noexcept
	{
		if (expansion_cached)
		{
			free(expansion_lambda_x);
			free(expansion_lambda_y);
			free(expansion_lambda_z);
			free(expansion_denominator);
		}
		if (exact_cached)
		{
			delete exact_lambda_x;
			delete exact_lambda_y;
			delete exact_lambda_z;
			delete exact_denominator;
		}
	}
};

#elif defined(OMC_OFFSET_PRED)

template <typename IT, typename ET>
class OnePointCachedValues3
{
public:
	using FT = double;

public:
	// floating-point filter
	#ifdef OMC_CACHE_SSF
	FT ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_lambda_z,
	  ssfilter_denominator, ssfilter_beta_x, ssfilter_beta_y, ssfilter_beta_z,
	  ssfilter_max_val;
	#endif
	// dynamic filter (interval number)
	IT dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator,
	  dfilter_beta_x, dfilter_beta_y, dfilter_beta_z;
	// exact number
	ET *exact_lambda_x = nullptr, *exact_lambda_y = nullptr,
	   *exact_lambda_z = nullptr, *exact_denominator = nullptr,
	   *exact_beta_x = nullptr, *exact_beta_y = nullptr, *exact_beta_z = nullptr;
	// expansion number
	FT *expansion_lambda_x = nullptr, *expansion_lambda_y = nullptr,
	   *expansion_lambda_z = nullptr, *expansion_denominator = nullptr;
	int expansion_lambda_x_len = 0, expansion_lambda_y_len = 0,
	    expansion_lambda_z_len = 0, expansion_d_len = 0;
	// use beta in ssfilter

	bool ssfilter_cached  = false;
	bool dfilter_cached   = false;
	bool exact_cached     = false;
	bool expansion_cached = false;

	OnePointCachedValues3() noexcept
	{
		exact_lambda_x    = nullptr;
		exact_lambda_y    = nullptr;
		exact_lambda_z    = nullptr;
		exact_denominator = nullptr;
		exact_beta_x      = nullptr;
		exact_beta_y      = nullptr;
		exact_beta_z      = nullptr;

		expansion_lambda_x     = nullptr;
		expansion_lambda_y     = nullptr;
		expansion_lambda_z     = nullptr;
		expansion_denominator  = nullptr;
		expansion_lambda_x_len = 0;
		expansion_lambda_y_len = 0;
		expansion_lambda_z_len = 0;
		expansion_d_len        = 0;

		ssfilter_cached  = false;
		dfilter_cached   = false;
		exact_cached     = false;
		expansion_cached = false;
	}

	OnePointCachedValues3(const OnePointCachedValues3 &rhs) = delete;
	OnePointCachedValues3(OnePointCachedValues3 &&rhs)      = delete;

	void operator=(const OnePointCachedValues3 &rhs) = delete;
	void operator=(OnePointCachedValues3 &&rhs)      = delete;

	~OnePointCachedValues3() noexcept
	{
		if (expansion_cached)
		{
			free(expansion_lambda_x);
			free(expansion_lambda_y);
			free(expansion_lambda_z);
			free(expansion_denominator);
		}
		if (exact_cached || exact_lambda_x)
		{
			delete[] exact_lambda_x;
		}
	}

	void alloc_ET()
	{
		ET *new_et        = new ET[7];
		exact_lambda_x    = new_et + 0;
		exact_lambda_y    = new_et + 1;
		exact_lambda_z    = new_et + 2;
		exact_denominator = new_et + 3;
		exact_beta_x      = new_et + 4;
		exact_beta_y      = new_et + 5;
		exact_beta_z      = new_et + 6;
	}
};

#endif

template <typename IT_, typename ET_, typename OnePointCachedValues_>
class GlobalCachedValues
{
public:
	using FT = double;
	using IT = IT_;
	using ET = ET_;

	using OnePointCachedValues = OnePointCachedValues_;
	using PointsCachedValuesMap =
	  phmap::flat_hash_map<void *, OnePointCachedValues *>;
	using PointsCachedValuesArena = std::deque<OnePointCachedValues>;

public:
	/// @brief resize to \p thread_num maps, \p thread_num is often the max
	/// thread number.
	void resize(size_t thread_num)
	{
		// spin lock
		while (spin_lock.test_and_set(std::memory_order_acquire))
			;
		m_maps.clear();
		m_maps.resize(thread_num);
		m_arenas.clear();
		m_arenas.resize(thread_num);
		spin_lock.clear(std::memory_order_release);
	}

	OnePointCachedValues &get(void *point_ptr)
	{
		int thread_id = tbb::this_task_arena::current_thread_index();
		OMC_EXPENSIVE_ASSERT((size_t)thread_id < m_maps.size(),
		                     "thread id {} excceed maps size {}", thread_id,
		                     m_maps.size());
		PointsCachedValuesMap &map  = m_maps[thread_id];
		auto                   iter = map.find(point_ptr);
		if (iter != map.end())
		{
			return *iter->second;
		}
		else
		{
			PointsCachedValuesArena &arena = m_arenas[thread_id];
			arena.emplace_back();
			OnePointCachedValues *cv = &arena.back();
			map.insert(std::pair<void *, OnePointCachedValues *>(point_ptr, cv));
			return *cv;
		}
	}

	void remove(void *point_ptr)
	{
		int thread_id = tbb::this_task_arena::current_thread_index();
		OMC_EXPENSIVE_ASSERT((size_t)thread_id < m_maps.size(),
		                     "thread id {} excceed maps size {}", thread_id,
		                     m_maps.size());
		m_maps[thread_id].erase(point_ptr);
	}

	void clear_cached_values()
	{
		int thread_id = tbb::this_task_arena::current_thread_index();
		OMC_EXPENSIVE_ASSERT((size_t)thread_id < m_maps.size(),
		                     "thread id {} excceed maps size {}", thread_id,
		                     m_maps.size());
		m_maps[thread_id].clear();
		m_arenas[thread_id].clear();
	}

	void enable()
	{
		// spin lock
		while (spin_lock.test_and_set(std::memory_order_acquire))
			;
		m_enabled = true;
		spin_lock.clear(std::memory_order_release);
	}

	void disable()
	{
		// spin lock
		while (spin_lock.test_and_set(std::memory_order_acquire))
			;
		m_enabled = false;
		m_maps.clear();
		spin_lock.clear(std::memory_order_release);
	}

	bool is_enabled() const { return m_enabled; }

protected:
	/// one map per thread, map generic point's pointer to cached value
	std::deque<PointsCachedValuesMap>   m_maps;
	std::deque<PointsCachedValuesArena> m_arenas;

	bool m_enabled = false;

	/// mutex for multi-thread safety
	std::atomic_flag spin_lock = ATOMIC_FLAG_INIT;
};

} // namespace OMC