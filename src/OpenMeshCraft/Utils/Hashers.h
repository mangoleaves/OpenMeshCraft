#pragma once

#include <array>
#include <vector>

namespace OMC {

// This struct defines a hash function for std::array types.
// It uses a custom hash algorithm to combine the hashes of the elements in the
// array.
template <class Key>
struct hash;

template <typename T, size_t N>
struct hash<std::array<T, N>>
{
	std::hash<T> hasher; // Hash function for the element type T.

	// This operator overload calculates the hash for a given std::array.
	inline size_t operator()(const std::array<T, N> &p) const
	{
		size_t seed = 0; // Initialize the seed for the hash.
		for (size_t i = 0; i < N; i++)
			// Combine the hash of each element with the seed using a custom
			// algorithm.
			seed ^= hasher(p[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed; // Return the final computed hash.
	}
};

// This struct defines a hash function for std::pair types.
// It combines the hashes of the two elements in the pair.
template <typename T>
struct hash<std::pair<T, T>>
{
	std::hash<T> hasher; // Hash function for the element type T.

	// This operator overload calculates the hash for a given std::pair.
	inline size_t operator()(const std::pair<T, T> &p) const
	{
		size_t seed = 0; // Initialize the seed for the hash.
		// Combine the hash of the first element with the seed.
		seed ^= hasher(p.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		// Combine the hash of the second element with the seed.
		seed ^= hasher(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed; // Return the final computed hash.
	}
};

// This struct defines a hash function for std::vector types.
// It combines the hashes of the elements in the vector.
template <typename T>
struct hash<std::vector<T>>
{
	std::hash<T> hasher; // Hash function for the element type T.

	// This operator overload calculates the hash for a given std::vector.
	inline size_t operator()(const std::vector<T> &p) const
	{
		size_t seed = 0; // Initialize the seed for the hash.
		for (size_t i = 0; i < p.size(); i++)
			// Combine the hash of each element with the seed using a custom
			// algorithm.
			seed ^= hasher(p[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed; // Return the final computed hash.
	}
};

} // namespace OMC