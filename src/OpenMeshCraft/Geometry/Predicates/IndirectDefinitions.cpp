#include "IndirectDefinitions.h"

#include <format>
#include <iostream>
#include <string>
#include <vector>

namespace OMC {

#ifdef OMC_PRED_PROFILE

uint32_t PredicatesProfile::total_count[(size_t)PredicateNames::CNT][ARR_CNT];
uint32_t PredicatesProfile::ss_fail_count[(size_t)PredicateNames::CNT][ARR_CNT];
uint32_t PredicatesProfile::d_fail_count[(size_t)PredicateNames::CNT][ARR_CNT];

void PredicatesProfile::initialize()
{
	for (size_t i = 0; i < (size_t)PredicateNames::CNT; i++)
	{
		for (size_t j = 0; j < ARR_CNT; j++)
		{
			total_count[i][j]   = 0;
			ss_fail_count[i][j] = 0;
			d_fail_count[i][j]  = 0;
		}
	}
}

void PredicatesProfile::inc_total(PredicateNames name, PntArr3 arr)
{
	// clang-format off
	switch (arr)
	{
	case PntArr3::S:   total_count[(size_t)name][0] += 1;break;
	case PntArr3::L:   total_count[(size_t)name][1] += 1;break;
	case PntArr3::T:   total_count[(size_t)name][2] += 1;break;
	case PntArr3::SS:  total_count[(size_t)name][3] += 1;break;
	case PntArr3::SL:  total_count[(size_t)name][4] += 1;break;
	case PntArr3::ST:  total_count[(size_t)name][5] += 1;break;
	case PntArr3::LL:  total_count[(size_t)name][6] += 1;break;
	case PntArr3::LT:  total_count[(size_t)name][7] += 1;break;
	case PntArr3::TT:  total_count[(size_t)name][8] += 1;break;
	case PntArr3::SSS: total_count[(size_t)name][9] += 1;break;
	case PntArr3::SSL: total_count[(size_t)name][10] += 1;break;
	case PntArr3::SST: total_count[(size_t)name][11] += 1;break;
	case PntArr3::SLL: total_count[(size_t)name][12] += 1;break;
	case PntArr3::SLT: total_count[(size_t)name][13] += 1;break;
	case PntArr3::STT: total_count[(size_t)name][14] += 1;break;
	case PntArr3::LLL: total_count[(size_t)name][15] += 1;break;
	case PntArr3::LLT: total_count[(size_t)name][16] += 1;break;
	case PntArr3::LTT: total_count[(size_t)name][17] += 1;break;
	case PntArr3::TTT: total_count[(size_t)name][18] += 1;break;
	}
	// clang-format on
}

void PredicatesProfile::inc_ss_fail(PredicateNames name, PntArr3 arr)
{
	// clang-format off
	switch (arr)
	{
	case PntArr3::S:   ss_fail_count[(size_t)name][0] += 1;break;
	case PntArr3::L:   ss_fail_count[(size_t)name][1] += 1;break;
	case PntArr3::T:   ss_fail_count[(size_t)name][2] += 1;break;
	case PntArr3::SS:  ss_fail_count[(size_t)name][3] += 1;break;
	case PntArr3::SL:  ss_fail_count[(size_t)name][4] += 1;break;
	case PntArr3::ST:  ss_fail_count[(size_t)name][5] += 1;break;
	case PntArr3::LL:  ss_fail_count[(size_t)name][6] += 1;break;
	case PntArr3::LT:  ss_fail_count[(size_t)name][7] += 1;break;
	case PntArr3::TT:  ss_fail_count[(size_t)name][8] += 1;break;
	case PntArr3::SSS: ss_fail_count[(size_t)name][9] += 1;break;
	case PntArr3::SSL: ss_fail_count[(size_t)name][10] += 1;break;
	case PntArr3::SST: ss_fail_count[(size_t)name][11] += 1;break;
	case PntArr3::SLL: ss_fail_count[(size_t)name][12] += 1;break;
	case PntArr3::SLT: ss_fail_count[(size_t)name][13] += 1;break;
	case PntArr3::STT: ss_fail_count[(size_t)name][14] += 1;break;
	case PntArr3::LLL: ss_fail_count[(size_t)name][15] += 1;break;
	case PntArr3::LLT: ss_fail_count[(size_t)name][16] += 1;break;
	case PntArr3::LTT: ss_fail_count[(size_t)name][17] += 1;break;
	case PntArr3::TTT: ss_fail_count[(size_t)name][18] += 1;break;
	}
	// clang-format on
}

void PredicatesProfile::inc_d_fail(PredicateNames name, PntArr3 arr)
{
	// clang-format off
	switch (arr)
	{
	case PntArr3::S:   d_fail_count[(size_t)name][0] += 1;break;
	case PntArr3::L:   d_fail_count[(size_t)name][1] += 1;break;
	case PntArr3::T:   d_fail_count[(size_t)name][2] += 1;break;
	case PntArr3::SS:  d_fail_count[(size_t)name][3] += 1;break;
	case PntArr3::SL:  d_fail_count[(size_t)name][4] += 1;break;
	case PntArr3::ST:  d_fail_count[(size_t)name][5] += 1;break;
	case PntArr3::LL:  d_fail_count[(size_t)name][6] += 1;break;
	case PntArr3::LT:  d_fail_count[(size_t)name][7] += 1;break;
	case PntArr3::TT:  d_fail_count[(size_t)name][8] += 1;break;
	case PntArr3::SSS: d_fail_count[(size_t)name][9] += 1;break;
	case PntArr3::SSL: d_fail_count[(size_t)name][10] += 1;break;
	case PntArr3::SST: d_fail_count[(size_t)name][11] += 1;break;
	case PntArr3::SLL: d_fail_count[(size_t)name][12] += 1;break;
	case PntArr3::SLT: d_fail_count[(size_t)name][13] += 1;break;
	case PntArr3::STT: d_fail_count[(size_t)name][14] += 1;break;
	case PntArr3::LLL: d_fail_count[(size_t)name][15] += 1;break;
	case PntArr3::LLT: d_fail_count[(size_t)name][16] += 1;break;
	case PntArr3::LTT: d_fail_count[(size_t)name][17] += 1;break;
	case PntArr3::TTT: d_fail_count[(size_t)name][18] += 1;break;
	}
	// clang-format on
}

void PredicatesProfile::print()
{
	// clang-format off
  std::vector<std::string> pred_names = {
	  "lessThanOnX_IE",
	  "lessThanOnX_II",
	  "lessThanOnY_IE",
	  "lessThanOnY_II",
	  "lessThanOnZ_IE",
	  "lessThanOnZ_II",
	  "orientOn2Dxy_IEE",
	  "orientOn2Dxy_IIE",
	  "orientOn2Dxy_III",
	  "orientOn2Dyz_IEE",
	  "orientOn2Dyz_IIE",
	  "orientOn2Dyz_III",
	  "orientOn2Dzx_IEE",
	  "orientOn2Dzx_IIE",
	  "orientOn2Dzx_III",
  };

	std::vector<std::string> arr_names = {
 	  "S"  ,
	  "L"  ,
	  "T"  ,
	  "SS" ,
	  "SL" ,
	  "ST" ,
	  "LL" ,
	  "LT" ,
	  "TT" ,
	  "SSS",
	  "SSL",
	  "SST",
	  "SLL",
	  "SLT",
	  "STT",
	  "LLL",
	  "LLT",
	  "LTT",
	  "TTT"
	};
	// clang-format on

	for (size_t i = 0; i < (size_t)PredicateNames::CNT; i++)
	{
		std::cout << std::format("{}:\n", pred_names[i]);
		for (size_t j = 0; j < ARR_CNT; j++)
		{
			if (total_count[i][j] != 0)
			{
				double ss_succeed =
				  1. - (double)ss_fail_count[i][j] / (double)total_count[i][j];
				double d_succeed =
				  (1. - (double)d_fail_count[i][j] / (double)total_count[i][j]) -
				  ss_succeed;
				double e_succeed =
				  (double)d_fail_count[i][j] / (double)total_count[i][j];
				std::cout << std::format("  {}: {:.2f}% {:.2f}% {:.2f}%\n",
				                         arr_names[j], ss_succeed * 100.,
				                         d_succeed * 100., e_succeed * 100.);
			}
		}
	}
}

#endif

} // namespace OMC