#pragma once

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "boost/multiprecision/gmp.hpp"
#include "boost/operators.hpp"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

namespace OMC {

using BoostMpzInt   = boost::multiprecision::mpz_int;
using BoostRational = boost::multiprecision::mpq_rational;
using BoostFloat    = boost::multiprecision::mpf_float;

template <typename Protected>
class IntervalNumber;

template <typename HasApprox, typename Protected, typename ThreadSafe>
class LazyNumber;

} // namespace OMC