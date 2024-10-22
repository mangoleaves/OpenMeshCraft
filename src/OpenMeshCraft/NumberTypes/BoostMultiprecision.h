#pragma once

#include "NumberUtils.h"

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/Macros.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#if defined(_WIN32) && defined(__GNUC__)
  // When I compile on windows with gcc14, I get errors about duplicate
  // definitions from STL and Boost. Adding below macro succeeds to bypass the
  // problem. But the macro should probably be removed if no such error occurs
  // in the future.
  // In addtion, such problem does not occur when compiling on windows with MSVC
  // or on linux-OS with gcc. Strange, right?
	#define BOOST_CSTDFLOAT_NO_LIBQUADMATH_SUPPORT
#endif
#include "boost/multiprecision/gmp.hpp"
#include "boost/operators.hpp"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

namespace OMC {

using BoostMpzInt   = boost::multiprecision::mpz_int;
using BoostRational = boost::multiprecision::mpq_rational;
using BoostFloat    = boost::multiprecision::mpf_float;

/***** Unary operators for BoostRational *****/

template <>
class UnaryOperators<BoostRational>
{
public:
	using NT = BoostRational;

public:
	static Sign sign(const NT &x)
	{
		if (x > 0)
			return Sign::POSITIVE;
		else if (x < 0)
			return Sign::NEGATIVE;
		else
			return Sign::ZERO;
	}

	static NT abs(const NT &x) { return boost::multiprecision::abs(x); }

	static NT negate(const NT &x) { return -x; }
};

template <>
class UnaryOperators<BoostFloat>
{
public:
	using NT = BoostFloat;

public:
	static Sign sign(const NT &x)
	{
		if (x > 0)
			return Sign::POSITIVE;
		else if (x < 0)
			return Sign::NEGATIVE;
		else
			return Sign::ZERO;
	}

	static NT abs(const NT &x) { return boost::multiprecision::abs(x); }

	static NT negate(const NT &x) { return -x; }
};

/***** Convertors for Boost numbers *****/

class BoostMpToInterval
{
public:
	std::pair<double, double> operator()(const BoostRational &x);

	std::pair<double, double> operator()(const BoostFloat &x);

	std::pair<double, double> operator()(const BoostMpzInt &x);
};

class BoostMpToDouble
{
public:
	double operator()(const BoostRational &x) { return x.convert_to<double>(); }

	double operator()(const BoostFloat &x) { return x.convert_to<double>(); }
};

template <>
inline double to_double(const BoostRational &n)
{
	return BoostMpToDouble()(n);
}

template <>
inline double to_double(const BoostFloat &n)
{
	return BoostMpToDouble()(n);
}

} // namespace OMC