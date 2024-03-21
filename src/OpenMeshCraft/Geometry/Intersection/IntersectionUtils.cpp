#include "IntersectionUtils.h"

#include <format>
#include <iostream>
#include <string>
#include <vector>

namespace OMC {

#ifdef OMC_INTER_PROFILE

uint32_t IntersectionProfile::total_count[(size_t)IntersectionNames::CNT];
uint32_t IntersectionProfile::reach_count[(size_t)IntersectionNames::CNT]
                                         [BRANCH_CNT];
uint32_t IntersectionProfile::reach_line[(size_t)IntersectionNames::CNT]
                                        [BRANCH_CNT];

void IntersectionProfile::initialize()
{
	for (size_t i = 0; i < (size_t)IntersectionNames::CNT; i++)
	{
		total_count[i] = 0;
		for (size_t j = 0; j < BRANCH_CNT; j++)
		{
			reach_count[i][j] = 0;
		}
	}
}

void IntersectionProfile::inc_total(IntersectionNames name)
{
	total_count[(size_t)name] += 1;
}

void IntersectionProfile::inc_reach(IntersectionNames name,
                                    uint32_t branch_flag, uint32_t branch_line)
{
	reach_count[(size_t)name][branch_flag] += 1;
	reach_line[(size_t)name][branch_flag] = branch_line;
}

void IntersectionProfile::print()
{
	// clang-format off
  std::vector<std::string> func_names = {
	  "Triangle3_Triangle3",
	  "Triangle3_Segment3"
  };

	// clang-format on

	for (size_t i = 0; i < (size_t)IntersectionNames::CNT; i++)
	{
		std::cout << std::format("{}:\n", func_names[i]);
		int last_branch_flag = 0;
		// find the last non-zero branch flag
		for (int j = (int)BRANCH_CNT - 1; j >= 0; j--)
		{
			if (reach_count[i][j] != 0)
			{
				last_branch_flag = j;
				break;
			}
		}

		for (int j = 0; j < last_branch_flag; j++)
		{
			double reach_raio = (double)reach_count[i][j] / (double)total_count[i];
			std::cout << std::format("  line {}: {:.2f}%\n", reach_line[i][j],
			                         reach_raio * 100.);
		}
	}
}

#endif

} // namespace OMC