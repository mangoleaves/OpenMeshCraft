#pragma once

namespace OMC {

/**
 * @brief Primitive3 is the base class of all primitives in 3D.
 * @tparam _NT number type.
 */
template <typename _NT>
class Primitive3
{
public:
	using NT = _NT;

	static constexpr size_t dim = 3;
};

} // namespace OMC