#pragma once

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/Macros.h"

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

#if defined(INDIRECT_PREDICATES)

template <typename IT, typename ET>
class OnePointCachedValues3
{
public:
	using FT = double;

public:
	// floating-point filter
	// NaN demoninator is used to check if it is cached.
	FT ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_lambda_z,
	  ssfilter_denominator, ssfilter_max_val;
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

#elif defined(OFFSET_PREDICATES)

template <typename IT, typename ET>
class OnePointCachedValues3
{
public:
	using FT = double;

public:
	// floating-point filter
	// NaN demoninator is used to check if it is cached.
	FT ssfilter_lambda_x, ssfilter_lambda_y, ssfilter_lambda_z,
	  ssfilter_denominator, ssfilter_beta_x, ssfilter_beta_y, ssfilter_beta_z,
	  ssfilter_max_val;
	// dynamic filter (interval number)
	// NaN demoninator is used to check if it is cached.
	IT dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator,
	  dfilter_beta_x, dfilter_beta_y, dfilter_beta_z;
	// exact number
	// exact cached is used to check if it is cached.
	ET *exact_lambda_x = nullptr, *exact_lambda_y = nullptr,
	   *exact_lambda_z = nullptr, *exact_denominator = nullptr,
	   *exact_beta_x = nullptr, *exact_beta_y = nullptr, *exact_beta_z = nullptr;
	// expansion number
	// nullptr denominator is used to check if it is cached.
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
		exact_beta_x      = nullptr;
		exact_beta_y      = nullptr;
		exact_beta_z      = nullptr;
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
	  std::unordered_map<void *, OnePointCachedValues>;

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
		spin_lock.clear(std::memory_order_release);
	}

	OnePointCachedValues &get(void *point_ptr)
	{
		int thread_id = tbb::this_task_arena::current_thread_index();
		OMC_EXPENSIVE_ASSERT((size_t)thread_id < m_maps.size(),
		                     "thread id {} excceed maps size {}", thread_id,
		                     m_maps.size());
		OnePointCachedValues &value = m_maps[thread_id][point_ptr];
		return value;
	}

	void remove(void *point_ptr)
	{
		if (m_enabled)
		{
			int thread_id = tbb::this_task_arena::current_thread_index();
			OMC_EXPENSIVE_ASSERT((size_t)thread_id < m_maps.size(),
			                     "thread id {} excceed maps size {}", thread_id,
			                     m_maps.size());
			m_maps[thread_id].erase(point_ptr);
		}
	}

	void clear_cached_values()
	{
		if (m_enabled)
		{
			int thread_id = tbb::this_task_arena::current_thread_index();
			OMC_EXPENSIVE_ASSERT((size_t)thread_id < m_maps.size(),
			                     "thread id {} excceed maps size {}", thread_id,
			                     m_maps.size());
			m_maps[thread_id].clear();
		}
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
	std::deque<PointsCachedValuesMap> m_maps;

	bool m_enabled = false;

	/// mutex for multi-thread safety
	std::atomic_flag spin_lock = ATOMIC_FLAG_INIT;
};

} // namespace OMC