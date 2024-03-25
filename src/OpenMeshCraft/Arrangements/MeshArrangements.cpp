#include "MeshArrangements.h"

namespace OMC {

template class MeshArrangements<EIAC, TriSoupTraits>;

#ifdef OMC_ARR_PROFILE

size_t ArrProfile::total_count[(size_t)ArrFuncNames::CNT];
size_t ArrProfile::reach_count[(size_t)ArrFuncNames::CNT][BRANCH_CNT];
size_t ArrProfile::reach_line[(size_t)ArrFuncNames::CNT][BRANCH_CNT];

#endif

} // namespace OMC