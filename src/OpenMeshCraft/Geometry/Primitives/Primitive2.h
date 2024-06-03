#pragma once

#include <cstddef>

namespace OMC {

/**
 * @brief Primitive2 is the base class of all primitives in 2D.
 * @tparam _NT number type.
 */
template <typename _NT>
class Primitive2
{
public:
	using NT = _NT;

	static constexpr size_t dim = 2;
};

} // namespace OMC