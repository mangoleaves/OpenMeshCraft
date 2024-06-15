#pragma once

#include "OpenMeshCraft/Utils/Exception.h"

#include <cfenv>
#include <type_traits>

namespace OMC {

// For MSVC:
// Make sure adding compile option "/fp:strict"

// FOR GCC:
// Semantics of Floating Point Math in GCC
//    -frounding-math -fsignaling-nans
// Support infinities, NaNs, gradual underflow, signed zeros, exception flags
// and traps, setting rounding modes.
// Compare with C99's
//    #pragma STDC FENV ACCESS ON. #pragma STDC FENV_ACCESS ON

using fpu_round_t = decltype(FE_UPWARD);

inline fpu_round_t FPU_get_round() { return fegetround(); }

inline void FPU_set_round(fpu_round_t round) { fesetround(round); }

inline fpu_round_t FPU_get_then_set_round(fpu_round_t round)
{
	fpu_round_t old = FPU_get_round();
	FPU_set_round(round);
	return old;
}

/// @brief A class whose constructor sets the FPU mode to +inf, saves a backup
/// of it, and whose destructor resets it back to the saved state.
/// @tparam Protected Whether the rounding mode is protected in the context.
template <typename Protected = std::false_type>
class FPU_RoundingProtector;

/// @brief Rounding mode is protected in the context, this protector do nothing.
template <>
class FPU_RoundingProtector<std::true_type>
{
public:
	FPU_RoundingProtector()
	{
		OMC_EXPENSIVE_ASSERT(FPU_get_round() == FE_UPWARD,
		                     "round mode is not protected.");
	}
	~FPU_RoundingProtector() {}
};

/// @brief Rounding mode is not protected, this protector will try to protect.
template <>
class FPU_RoundingProtector<std::false_type>
{
public:
	FPU_RoundingProtector()
	  : backup(FPU_get_then_set_round(FE_UPWARD))
	{
	}
	~FPU_RoundingProtector() { FPU_set_round(backup); }

private:
	fpu_round_t backup;
};

} // namespace OMC