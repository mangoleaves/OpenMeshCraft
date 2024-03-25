#pragma once

// clang-format off
#include "OpenMeshCraft/Utils/DisableWarnings.h"
#include "tbb/tbb.h"
#include "parallel_hashmap/phmap.h"
#include "OpenMeshCraft/Utils/EnableWarnings.h"
// clang-format on

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/IndexDef.h"

#include "OpenMeshCraft/Geometry/Intersection/IntersectionUtils.h"
#include "OpenMeshCraft/NumberTypes/NumberUtils.h"

#include <algorithm>
#include <bitset>
#include <deque>
#include <execution>
#include <iterator>
#include <memory>
#include <vector>

namespace OMC {

constexpr int NBIT = 32;

/// * Label for each triangle in arrangements (boolean, and other applications).
///   Label indicates where this triangle locates on.
///   (For example, if Label[0] and Label[1] is true, this triangle locates on
///   triangle soup with index 0 and triangle soup with index 1 at the same
///   time.)
/// * Users can use this label to transfer attributes from input triangle soups
///   to the output triangle soup.
/// * Consider efficiency, we limit maximal label count to be less than
///   NBIT(32).
using Label = std::bitset<NBIT>;

struct MeshArrangements_Stats
{
	/* Preprocessing *************************************************/

	double pp_elapsed   = 0.; // timings of preprocessing
	double mdv_elapsed  = 0.; // timings of merging duplicate vertices
	double rdt_elapsed  = 0.; // timings of removing degenerate and
	                          // duplicate triangles
	double tree_elapsed = 0.; // timings of building tree

	/* Detect intersections ******************************************/

	double di_elapsed = 0.; // timings of detecting intersection

	/* connect intersetion *******************************************/

	double cn_elapsed = 0.; // timings of connecting intersection

	/* classify intersection *****************************************/

	double ci_elapsed = 0.; // timings of classifying intersection

	/* triangulation *************************************************/

	double tr_elapsed = 0.; // timings of triangulation
};

/// @brief Each output triangle has a surface label and inside label.
struct Labels
{
	/// surface (mesh) that a triangle belongs to.
	std::vector<Label> surface;
	/// surface (mesh) that a triangle is inside.
	std::vector<Label> inside;
	/// how many surfaces (meshes) in arrangements.
	size_t             num;
};

using UIPair = std::pair<index_t, index_t>;

enum Plane
{
	YZ = 0,
	ZX = 1,
	XY = 2
};

inline size_t LabelToIdx(const Label &b)
{
	OMC_EXPENSIVE_ASSERT(b.count() == 1, "more than 1 bit set to 1");
	uint32_t x = b.to_ulong();
#ifdef __GNUC__
	return __builtin_ffs(x) - 1; // GCC
#endif
#ifdef _MSC_VER
	return 32 - __lzcnt(x) - 1; // Visual studio
#endif
}

inline UIPair uniquePair(index_t i, index_t j)
{
	return (i < j) ? std::make_pair(i, j) : std::make_pair(j, i);
}

/// "<" in case of unique, "==" in case of two invalid indices.
inline bool isUnique(const UIPair &p) { return p.first <= p.second; }

struct ShewchukCache
{
	double minor[3];
	double perm[3];
};

using ShewchukCachePtr = const ShewchukCache *;

inline Plane intToPlane(const int &norm) { return static_cast<Plane>(norm); }
inline int   planeToInt(const Plane &p) { return static_cast<int>(p); }

template <typename Traits>
struct PointArena
{
public:
	using EPoint     = typename Traits::EPoint;
	using GPoint     = typename Traits::GPoint;
	using IPoint_SSI = typename Traits::IPoint_SSI;
	using IPoint_LPI = typename Traits::IPoint_LPI;
	using IPoint_TPI = typename Traits::IPoint_TPI;

public:
	std::vector<EPoint>    init; // explicit points
	std::deque<IPoint_SSI> ssi;  // SSI points
	std::deque<IPoint_LPI> lpi;  // LPI points
	std::deque<IPoint_TPI> tpi;  // TPI points

public:
	void recycle(IPoint_SSI *ssi_ptr) { recycled_ssi.push(ssi_ptr); }
	void recycle(IPoint_LPI *lpi_ptr) { recycled_lpi.push(lpi_ptr); }
	void recycle(IPoint_TPI *tpi_ptr) { recycled_tpi.push(tpi_ptr); }

	IPoint_SSI *emplace(IPoint_SSI &&p)
	{
		IPoint_SSI *ssi_ptr = nullptr;
		if (recycled_ssi.empty())
		{
			ssi.emplace_back(std::move(p));
			ssi_ptr = &ssi.back();
		}
		else
		{
			ssi_ptr = recycled_ssi.front();
			recycled_ssi.pop();
			*ssi_ptr = std::move(p);
		}
		return ssi_ptr;
	}

	IPoint_LPI *emplace(IPoint_LPI &&p)
	{
		IPoint_LPI *lpi_ptr = nullptr;
		if (recycled_lpi.empty())
		{
			lpi.emplace_back(std::move(p));
			lpi_ptr = &lpi.back();
		}
		else
		{
			lpi_ptr = recycled_lpi.front();
			recycled_lpi.pop();
			*lpi_ptr = std::move(p);
		}
		return lpi_ptr;
	}

	IPoint_TPI *emplace(IPoint_TPI &&p)
	{
		IPoint_TPI *tpi_ptr = nullptr;
		if (recycled_tpi.empty())
		{
			tpi.emplace_back(std::move(p));
			tpi_ptr = &tpi.back();
		}
		else
		{
			tpi_ptr = recycled_tpi.front();
			recycled_tpi.pop();
			*tpi_ptr = std::move(p);
		}
		return tpi_ptr;
	}

	void reserve_ssi(size_t new_n)
	{
		for (size_t i = 0; i < new_n; i++)
			recycle(&ssi.emplace_back());
	}

	void reserve_lpi(size_t new_n)
	{
		for (size_t i = 0; i < new_n; i++)
			recycle(&lpi.emplace_back());
	}

	void reserve_tpi(size_t new_n)
	{
		for (size_t i = 0; i < new_n; i++)
			recycle(&tpi.emplace_back());
	}

private:
	std::queue<IPoint_SSI *> recycled_ssi;
	std::queue<IPoint_LPI *> recycled_lpi;
	std::queue<IPoint_TPI *> recycled_tpi;
};

class IdxArena
{
public:
	std::deque<std::atomic<index_t>> indices; // points' indices
	/// NOTE indices is not global unique, each local patch may have
	/// different local index of a global unique point.

public:
	void recycle(std::atomic<index_t> *idx_ptr) { recycled_idx.push(idx_ptr); }

	std::atomic<index_t> *emplace(index_t idx)
	{
		std::atomic<index_t> *idx_ptr = nullptr;
		if (recycled_idx.empty())
		{
			indices.emplace_back(idx);
			idx_ptr = &indices.back();
		}
		else
		{
			idx_ptr = recycled_idx.front();
			recycled_idx.pop();
			idx_ptr->store(idx);
		}
		return idx_ptr;
	}

	/// @brief reserve new number of indices
	void reserve(size_t new_n)
	{
		for (size_t i = 0; i < new_n; i++)
			recycle(&indices.emplace_back());
	}

private:
	using RecycleQueue = std::queue<std::atomic<index_t> *>;

	RecycleQueue recycled_idx;
};

/// @note Only used this inlined vector to store simple types!!!
template <typename T, size_t N>
class InlinedVector
{
	// N is not zero and N is power of 2.
	static_assert(N != 0 && (N & (N - 1)) == 0);

public:
	using value_type = T;

public: /* Constructors *****************************************************/
	InlinedVector() noexcept
	  : cur_elem_size(0)
	{
	}

	InlinedVector(const InlinedVector &rhs) noexcept { operator=(rhs); }

	InlinedVector &operator=(const InlinedVector &rhs) noexcept
	{
		if (rhs.cur_elem_size <= N)
		{
			std::copy(rhs.inlined_data.begin(),
			          rhs.inlined_data.begin() + rhs.cur_elem_size,
			          inlined_data.begin());
			heap_data.clear();
		}
		else
		{
			inlined_data = rhs.inlined_data;
			heap_data    = rhs.heap_data;
		}
		cur_elem_size = rhs.cur_elem_size;
		return *this;
	}

	InlinedVector(InlinedVector &&rhs) noexcept { operator=(std::move(rhs)); }

	InlinedVector &operator=(InlinedVector &&rhs) noexcept
	{
		if (rhs.cur_elem_size <= N)
		{
			std::copy(rhs.inlined_data.begin(),
			          rhs.inlined_data.begin() + rhs.cur_elem_size,
			          inlined_data.begin());
			heap_data.clear();
		}
		else
		{
			inlined_data = rhs.inlined_data;
			heap_data    = std::move(rhs.heap_data);
		}
		cur_elem_size = rhs.cur_elem_size;
		return *this;
	}

public: /* Modifier **********************************************************/
	void push_back(const T &elem)
	{
		if (cur_elem_size < N)
			inlined_data[cur_elem_size] = elem;
		else
			heap_data.push_back(elem);
		cur_elem_size++;
	}
	void push_back(T &&elem)
	{
		if (cur_elem_size < N)
			inlined_data[cur_elem_size] = std::move(elem);
		else
			heap_data.push_back(elem);
		cur_elem_size++;
	}

	void pop_back()
	{
		if (cur_elem_size > N)
			heap_data.pop_back();
		cur_elem_size--;
	}

	void clear()
	{
		cur_elem_size = 0;
		heap_data.clear();
	}

	void resize(size_t s)
	{ /* do nothing */
	}

	void reserve(size_t s)
	{ /* do nothing */
	}

public: /* Access ************************************************************/
	const T &operator[](size_t idx) const
	{
		if (idx < N)
			return inlined_data[idx];
		else
			return heap_data[idx - N];
	}

	T &operator[](size_t idx)
	{
		if (idx < N)
			return inlined_data[idx];
		else
			return heap_data[idx - N];
	}

	T       &back() { return operator[](cur_elem_size - 1); }
	const T &back() const { return operator[](cur_elem_size - 1); }

	size_t size() const { return (size_t)cur_elem_size; }

	bool empty() const { return cur_elem_size == 0; }

public: /* Iterator **********************************************************/
	      // clang-format off
	class ForwardNoSkipIterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = T*;
		using reference         = T&;
	public:
		ForwardNoSkipIterator(InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		ForwardNoSkipIterator &operator++() { cur_idx++; return *this; }
		ForwardNoSkipIterator &operator--() { cur_idx--; return *this; }
		ForwardNoSkipIterator operator++(int) { ForwardNoSkipIterator cpy(*this); cur_idx++; return cpy; }
		ForwardNoSkipIterator operator--(int) { ForwardNoSkipIterator cpy(*this); cur_idx--; return cpy; }
		T &operator*() const { return vec->operator[](cur_idx); }
		T &operator*() { return vec->operator[](cur_idx); }
		T *operator->() const { return &vec->operator[](cur_idx); }
		T *operator->() { return &vec->operator[](cur_idx); }
		bool operator==(const ForwardNoSkipIterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const ForwardNoSkipIterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		InlinedVector   *vec;
		difference_type  cur_idx;
	};

	ForwardNoSkipIterator begin() {return ForwardNoSkipIterator(this, 0);}
	ForwardNoSkipIterator end() {return ForwardNoSkipIterator(this, (typename ForwardNoSkipIterator::difference_type)cur_elem_size);}

	class ForwardNoSkipConstIterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = const T*;
		using reference         = const T&;
	public:
		ForwardNoSkipConstIterator(const InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		ForwardNoSkipConstIterator &operator++() { cur_idx++; return *this; }
		ForwardNoSkipConstIterator &operator--() { cur_idx--; return *this; }
		ForwardNoSkipConstIterator operator++(int) { ForwardNoSkipConstIterator cpy(*this); cur_idx++; return cpy; }
		ForwardNoSkipConstIterator operator--(int) { ForwardNoSkipConstIterator cpy(*this); cur_idx--; return cpy; }
		const T &operator*() const { return vec->operator[](cur_idx); }
		const T &operator*() { return vec->operator[](cur_idx); }
		const T *operator->() const { return &vec->operator[](cur_idx); }
		const T *operator->() { return &vec->operator[](cur_idx); }
		bool operator==(const ForwardNoSkipConstIterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const ForwardNoSkipConstIterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		const InlinedVector   *vec;
		difference_type       cur_idx;
	};
	ForwardNoSkipConstIterator begin() const {return ForwardNoSkipConstIterator(this, 0);}
	ForwardNoSkipConstIterator end() const {return ForwardNoSkipConstIterator(this, (typename ForwardNoSkipIterator::difference_type)cur_elem_size);}
	      // clang-format on

public:
	std::array<T, N> inlined_data;
	std::vector<T>   heap_data;

	size_t cur_elem_size;
};

struct DuplTriInfo
{
	index_t t_id; // triangle id
	index_t l_id; // label id
	bool    w;    // winding (CW or CCW)
};

template <typename Points, typename Triangles, typename NT>
void load(const Points &points, const Triangles &triangles, const size_t label,
          std::vector<NT> &coords, std::vector<index_t> &flat_tris,
          std::vector<size_t> &labels)
{
	size_t p_off = coords.size() / 3; // prev num verts
	coords.resize(coords.size() + points.size() * 3);
	tbb::parallel_for(size_t(0), points.size(),
	                  [&coords, &points, &p_off](size_t p_id)
	                  {
		                  coords[(p_off + p_id) * 3]     = points[p_id][0];
		                  coords[(p_off + p_id) * 3 + 1] = points[p_id][1];
		                  coords[(p_off + p_id) * 3 + 2] = points[p_id][2];
	                  });

	size_t t_off = flat_tris.size() / 3; // prev num tris
	flat_tris.resize(flat_tris.size() + triangles.size() * 3);
	tbb::parallel_for(
	  size_t(0), triangles.size(),
	  [&flat_tris, &triangles, &t_off, &p_off](size_t t_id)
	  {
		  flat_tris[(t_off + t_id) * 3]     = p_off + triangles[t_id][0];
		  flat_tris[(t_off + t_id) * 3 + 1] = p_off + triangles[t_id][1];
		  flat_tris[(t_off + t_id) * 3 + 2] = p_off + triangles[t_id][2];
	  });

	size_t l_off = labels.size();
	labels.resize(labels.size() + triangles.size());
	std::fill(std::execution::par_unseq, labels.begin() + l_off, labels.end(),
	          label);
}

// t0 -> vertex ids of triangle t0, t1 -> vertex ids of triangle t1
inline bool consistentWinding(const index_t *t0, const index_t *t1)
{
	int j = 0;
	while (j < 3 && t0[0] != t1[j])
		j++;
	OMC_EXPENSIVE_ASSERT(j < 3, "not same triangle");
	return t0[1] == t1[(j + 1) % 3] && t0[2] == t1[(j + 2) % 3];
}

inline bool triContainsVert(index_t t_id, index_t v_id,
                            const std::vector<index_t> &in_tris)
{
	return in_tris[3 * t_id] == v_id || in_tris[3 * t_id + 1] == v_id ||
	       in_tris[3 * t_id + 2] == v_id;
}

template <typename T>
inline void remove_duplicates(std::vector<T> &values)
{
	std::sort(values.begin(), values.end());
	values.erase(std::unique(values.begin(), values.end()), values.end());
}

template <typename T, size_t N>
inline size_t remove_duplicates(std::array<T, N> &values)
{
	std::sort(values.begin(), values.end());
	return (size_t)(std::unique(values.begin(), values.end()) - values.begin());
}

template <typename T>
inline void parallel_remove_duplicates(std::vector<T> &values)
{
	tbb::parallel_sort(values.begin(), values.end());
	values.erase(std::unique(values.begin(), values.end()), values.end());
}

template <typename T, typename Alloc,
          template <typename T_, typename Alloc_> class Container>
inline bool contains(const Container<T, Alloc> &values, T value)
{
	return std::find(values.begin(), values.end(), value) != values.end();
}

template <typename T, typename Alloc,
          template <typename T_, typename Alloc_> class Container>
inline void reserve(Container<T, Alloc> &values, size_t s)
{
	for (T &v : values)
		v.reserve(s);
}

} // namespace OMC

namespace std {

template <typename T, size_t N>
struct hash<std::array<T, N>>
{
	std::hash<T>  hasher;
	inline size_t operator()(const std::array<T, N> &p) const
	{
		size_t seed = 0;
		for (size_t i = 0; i < N; i++)
			seed ^= hasher(p[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

template <typename T>
struct hash<std::pair<T, T>>
{
	std::hash<T>  hasher;
	inline size_t operator()(const std::pair<T, T> &p) const
	{
		size_t seed = 0;
		seed ^= hasher(p.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		seed ^= hasher(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

template <typename T>
struct hash<std::vector<T>>
{
	std::hash<T>  hasher;
	inline size_t operator()(const std::vector<T> &p) const
	{
		size_t seed = 0;
		for (size_t i = 0; i < p.size(); i++)
			seed ^= hasher(p[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

} // namespace std

namespace OMC {
#ifdef OMC_ARR_PROFILE

enum class ArrFuncNames : size_t
{
	D_BBI = 0,
	DC_TTI,
	CNT
};

struct ArrProfile
{
	static const uint32_t BRANCH_CNT = 64;

	static size_t total_count[(size_t)ArrFuncNames::CNT];
	static size_t reach_count[(size_t)ArrFuncNames::CNT][BRANCH_CNT];
	static size_t reach_line[(size_t)ArrFuncNames::CNT][BRANCH_CNT];

	static inline void initialize();
	static inline void inc_total(ArrFuncNames name);
	static inline void inc_reach(ArrFuncNames name, uint32_t branch_flag,
	                             uint32_t branch_line);

	static inline void inc_total(ArrFuncNames name, size_t count);
	static inline void inc_reach(ArrFuncNames name, uint32_t branch_flag,
	                             uint32_t branch_line, size_t count);

	static inline void print();
};

inline void ArrProfile::initialize()
{
	for (size_t i = 0; i < (size_t)ArrFuncNames::CNT; i++)
	{
		total_count[i] = 0;
		for (size_t j = 0; j < BRANCH_CNT; j++)
		{
			reach_count[i][j] = 0;
		}
	}
}

inline void ArrProfile::inc_total(ArrFuncNames name)
{
	total_count[(size_t)name] += 1;
}

inline void ArrProfile::inc_reach(ArrFuncNames name, uint32_t branch_flag,
                                  uint32_t branch_line)
{
	reach_count[(size_t)name][branch_flag] += 1;
	reach_line[(size_t)name][branch_flag] = branch_line;
}

inline void ArrProfile::inc_total(ArrFuncNames name, size_t count)
{
	total_count[(size_t)name] += count;
}

inline void ArrProfile::inc_reach(ArrFuncNames name, uint32_t branch_flag,
                                  uint32_t branch_line, size_t count)
{
	reach_count[(size_t)name][branch_flag] += count;
	reach_line[(size_t)name][branch_flag] = branch_line;
}

inline void ArrProfile::print()
{
	// clang-format off
  std::vector<std::string> func_names = {
		"Detect BBI",
	  "Detect & Classify TTI"
  };

	// clang-format on

	for (size_t i = 0; i < (size_t)ArrFuncNames::CNT; i++)
	{
		std::cout << std::format("{}:\n", func_names[i]);
		int last_branch_flag = 0;
		// find the last non-zero branch flag
		for (int j = (int)BRANCH_CNT - 1; j >= 0; j--)
		{
			if (reach_count[i][j] != 0)
			{
				last_branch_flag = j;
				break;
			}
		}

		for (int j = 0; j <= last_branch_flag; j++)
		{
			double reach_raio = (double)reach_count[i][j] / (double)total_count[i];
			std::cout << std::format("  line {}: {:.2f}%\n", reach_line[i][j],
			                         reach_raio * 100.);
		}
	}
}

	#define OMC_ARR_PROFILE_INIT OMC::ArrProfile::initialize()
	#define OMC_ARR_PROFILE_PRINT OMC::ArrProfile::print()

	#define OMC_ARR_PROFILE_INC_TOTAL(func) OMC::ArrProfile::inc_total(func)

	#define OMC_ARR_PROFILE_INC_REACH(func, branch_flag) \
		OMC::ArrProfile::inc_reach(func, branch_flag, __LINE__)

	#define OMC_ARR_PROFILE_INC_TOTAL_CNT(func, count) OMC::ArrProfile::inc_total(func, count)

	#define OMC_ARR_PROFILE_INC_REACH_CNT(func, branch_flag, count) \
		OMC::ArrProfile::inc_reach(func, branch_flag, __LINE__, count)

#else

	#define OMC_ARR_PROFILE_INIT
	#define OMC_ARR_PROFILE_PRINT

	#define OMC_ARR_PROFILE_INC_TOTAL(func)
	#define OMC_ARR_PROFILE_INC_REACH(func, branch_flag)

	#define OMC_ARR_PROFILE_INC_TOTAL_CNT(func, count)
	#define OMC_ARR_PROFILE_INC_REACH_CNT(func, branch_flag, count)

#endif
} // namespace OMC

#define OMC_ARR_START_ELAPSE(name) auto name = OMC::Logger::elapse_reset();

#define OMC_ARR_SAVE_ELAPSED(name, dst_name, dscrpt)                     \
	if (stats != nullptr)                                                  \
		stats->dst_name = OMC::Logger::elapsed(name).count();                \
	if (verbose)                                                           \
	{                                                                      \
		OMC::Logger::info(std::format("[OpenMeshCraft Arrangements] " dscrpt \
		                              " time : {} s",                        \
		                              OMC::Logger::elapsed(name).count()));  \
	}
