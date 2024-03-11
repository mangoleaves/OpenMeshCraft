#pragma once

#include "OpenMeshCraft/Utils/Macros.h"

#include <cassert>
#include <iostream>
#include <limits>

#undef FALSE
#undef TRUE

namespace OMC {

enum class Sign
{
	POSITIVE  = 1,
	NEGATIVE  = -1,
	ZERO      = 0,
	UNCERTAIN = 2
};

/// @brief Check if sign is POSITIVE or NEGATIVE.
inline bool is_sign_posneg(Sign s)
{
	return s == Sign::POSITIVE || s == Sign::NEGATIVE;
}

/// @brief Sign is reliable if it is not UNCERTAIN.
inline bool is_sign_reliable(Sign s) { return s != Sign::UNCERTAIN; }

/// @brief Used in numeric filter.
inline Sign filter_sign(double num, double eps)
{
	if (num > eps)
		return Sign::POSITIVE;
	else if (num < -eps)
		return Sign::NEGATIVE;
	else
		return Sign::UNCERTAIN;
}

/// @brief when sign is NEG, ZERO or POS, reverse the sign to POS, ZERO or NEG.
/// reverse_sign is undefined when s is uncertain.
inline Sign reverse_sign(Sign s)
{
	return static_cast<Sign>(-static_cast<int>(s));
}

/// sign stored in any container
template <typename Signs>
Signs reverse_signs(Signs signs)
{
	for (Sign &s : signs)
		s = reverse_sign(s);
	return signs;
}

enum class Certainty
{
	FALSE     = 0,
	TRUE      = 1,
	UNCERTAIN = 2
};

inline bool is_certain(Certainty c) { return c != Certainty::UNCERTAIN; }
inline bool is_uncertain(Certainty c) { return c == Certainty::UNCERTAIN; }
inline bool get_certain(Certainty c) { return c == Certainty::TRUE; }

/**
 * @brief A wrapper of std::numeric_limits.
 * @note When define new number types, extend NumericLimits to it.
 * @tparam NT number type.
 */
template <typename NT>
class NumericLimits
{
	// TODO: Add more numeric limits for LazeExactNumber
public:
	static NT min() { return std::numeric_limits<NT>::min(); }

	static NT max() { return std::numeric_limits<NT>::max(); }

	static NT lowest() { return std::numeric_limits<NT>::lowest(); }

	static NT epsilon() { return std::numeric_limits<NT>::epsilon(); }
};

/**
 * @brief A bundle of unary operators on different number types
 * @tparam NT number type
 */
template <typename NT>
class UnaryOperators
{
public:
	/// @brief Get the sign of a number.
	static Sign sign(const NT &n)
	{
		if (n > NT(0))
			return Sign::POSITIVE;
		else if (n < NT(0))
			return Sign::NEGATIVE;
		else
			return Sign::ZERO;
	}
	/// @brief absolute value
	static NT abs(const NT &n) { return std::abs(n); };
	/// @brief negate given number
	static NT negate(const NT &n) { return -n; }
};

/// @brief Get the sign of a number.
template <typename NT>
inline Sign sign(const NT &n)
{
	return UnaryOperators<NT>::sign(n);
}

/// @brief absolute value
template <typename NT>
inline NT abs(const NT &n)
{
	return UnaryOperators<NT>::abs(n);
}

template <typename NT>
class absolute
{
public:
	NT operator()(const NT &n) { return abs<NT>(n); }
};

/// @brief negate given number
template <typename NT>
inline NT negate(const NT &n)
{
	return UnaryOperators<NT>::negate(n);
}

template <typename NT>
inline double to_double(const NT &n)
{
	return static_cast<double>(n);
}

} // namespace OMC