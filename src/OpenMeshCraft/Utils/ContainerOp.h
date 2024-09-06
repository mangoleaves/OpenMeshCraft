#pragma once

// clang-format off
#include "OpenMeshCraft/Utils/DisableWarnings.h"
#include "tbb/tbb.h"
#include "OpenMeshCraft/Utils/EnableWarnings.h"
// clang-format on

#include <array>
#include <vector>

namespace OMC {
// This function removes duplicate elements from a vector of type T.
// It first sorts the vector, then resizes it to remove duplicates.
template <typename T>
inline void remove_duplicates(std::vector<T> &values)
{
	std::sort(values.begin(), values.end());
	values.resize(std::unique(values.begin(), values.end()) - values.begin());
}

// This function removes duplicate elements from a fixed-size array of type T.
// It sorts the array and returns the new size after removing duplicates.
template <typename T, size_t N>
inline size_t remove_duplicates(std::array<T, N> &values)
{
	std::sort(values.begin(), values.end());
	return (size_t)(std::unique(values.begin(), values.end()) - values.begin());
}

// This function removes duplicate elements from a concurrent vector of type T.
// It sorts the vector and then resizes it to eliminate duplicates.
template <typename T>
inline void remove_duplicates(tbb::concurrent_vector<T> &values)
{
	std::sort(values.begin(), values.end());
	values.resize(std::unique(values.begin(), values.end()) - values.begin());
}

// This function removes duplicate elements from a vector of type T in parallel.
// It uses parallel sorting to improve performance before resizing to remove duplicates.
template <typename T>
inline void parallel_remove_duplicates(std::vector<T> &values)
{
	tbb::parallel_sort(values.begin(), values.end());
	values.resize(std::unique(values.begin(), values.end()) - values.begin());
}

// This function checks if a given value exists in a container of type T.
// It returns true if the value is found, otherwise false.
template <typename T, typename Alloc,
          template <typename T_, typename Alloc_> class Container>
inline bool contains(const Container<T, Alloc> &values, T value)
{
	return std::find(values.begin(), values.end(), value) != values.end();
}

// This function reserves space in each element of a container of type T.
// It iterates through the values and reserves the specified size.
template <typename T, typename Alloc,
          template <typename T_, typename Alloc_> class Container>
inline void reserve(Container<T, Alloc> &values, size_t s)
{
	for (T &v : values)
		v.reserve(s);
}
} // namespace OMC