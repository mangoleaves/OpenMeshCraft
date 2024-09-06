#pragma once

#include "Exception.h"

#include <bitset>

namespace OMC
{

///  To ensure efficiency, we limit the maximum label count to less than 
///   LABEL_NBIT (32).
constexpr int LABEL_NBIT = 32;

using Label = std::bitset<LABEL_NBIT>;

inline size_t LabelToIdx(const Label &b)
{
	OMC_EXPENSIVE_ASSERT(b.count() == 1, "more than 1 bit set to 1");
	uint32_t x = b.to_ulong();
#ifdef __GNUC__
	return __builtin_ffs(x) - 1; // GCC
#endif
#ifdef _MSC_VER
	return 32 - __lzcnt(x) - 1; // Visual studio
#endif
}
}