#include "MeshArrangements.h"

namespace OMC {

template class MeshArrangements<EIAC, TriSoupTraits>;

#ifdef OMC_ARR_PROFILE

std::atomic_size_t ArrProfile::total_count[(size_t)ArrFuncNames::CNT];
std::atomic_size_t ArrProfile::reach_count[(size_t)ArrFuncNames::CNT][BRANCH_CNT];
std::atomic_size_t ArrProfile::reach_line[(size_t)ArrFuncNames::CNT][BRANCH_CNT];

#endif

} // namespace OMC