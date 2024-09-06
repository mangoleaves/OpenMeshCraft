#pragma once

#include <array>
#include <cstddef>

namespace OMC {

/* Index definition ********************************************************/

/// @brief Index unsed in mesh
using index_t = size_t;

/// @brief A value indicates invalid index. When index is unsigned type, `-1`
/// will be cast to the maximum unsigned value.
/// @todo Use type_traits for a safe initialization.
constexpr index_t InvalidIndex = index_t(-1);

/// @brief Make an index invalid.
inline void invalidate_idx(index_t &idx) { idx = InvalidIndex; }

/// @brief Check if an index is valid.
inline bool is_valid_idx(index_t idx) { return idx != InvalidIndex; }

/* Other index definition ***************************************************/

/// index pair
using IPair = std::pair<index_t, index_t>;
/// unique index pair (expect first <= second)
using UIPair = std::pair<index_t, index_t>;

inline UIPair unique_pair(index_t i, index_t j)
{
	return (i < j) ? std::make_pair(i, j) : std::make_pair(j, i);
}

/// "<" in case of unique, "==" in case of two invalid indices.
inline bool is_unique(const UIPair &p) { return p.first <= p.second; }

} // namespace OMC