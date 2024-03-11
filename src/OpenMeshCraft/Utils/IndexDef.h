#pragma once

#include <cstddef>

namespace OMC {

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

} // namespace OMC