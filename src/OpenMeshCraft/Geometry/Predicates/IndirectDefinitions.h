#pragma once

#include "OpenMeshCraft/Geometry/Primitives/GenericPoint2T.h"
#include "OpenMeshCraft/Geometry/Primitives/GenericPoint3T.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"

#include <array>
#include <atomic>

namespace OMC {

// Here we implement a type converter, to convert multiple generic points' types
// to an enum class type. The enum can be used in branch to make a jump table.

#define E2 0 // Explicit
#define I2 1 // Implicit
#define S2 2 // SSI

#define E3 0  // Explicit
#define I3 1  // Implicit
#define S3 2  // SSI
#define L3 3  // LPI
#define T3 4  // TPI

#define BEP 4 // bits for each point
// clang-format off
#define Map2(p0, p1)   \
  ((p1 << (BEP * 1)) + \
   (p0 << (BEP * 0)))
#define Map3(p0, p1, p2) \
  ((p2 << (BEP * 2)) +   \
   (p1 << (BEP * 1)) +   \
   (p0 << (BEP * 0)))
#define Map4(p0, p1, p2, p3) \
  ((p3 << (BEP * 3)) +       \
   (p2 << (BEP * 2)) +       \
   (p1 << (BEP * 1)) +       \
   (p0 << (BEP * 0)))
#define Map5(p0, p1, p2, p3, p4) \
  ((p4 << (BEP * 4)) +           \
   (p3 << (BEP * 3)) +           \
   (p2 << (BEP * 2)) +           \
   (p1 << (BEP * 1)) +           \
   (p0 << (BEP * 0)))
// clang-format on

// [smallest, ..., largest]
enum class PntArr2 : size_t
{
	E    = E2,
	I    = I2,
	S    = S2,
	// #################################
	EE   = Map2(E2, E2),
	// E/I * 2  (full arrangement)
	IE   = Map2(I2, E2),
	EI   = Map2(E2, I2),
	II   = Map2(I2, I2),
	// E/S * 2  (full arrangement)
	SE   = Map2(S2, E2),
	ES   = Map2(E2, S2),
	SS   = Map2(S2, S2),
	// #################################
	EEE  = Map3(E2, E2, E2),
	// E/I * 3  (full arrangement)
	IEE  = Map3(I2, E2, E2),
	EIE  = Map3(E2, I2, E2),
	IIE  = Map3(I2, I2, E2),
	EEI  = Map3(E2, E2, I2),
	IEI  = Map3(I2, E2, I2),
	EII  = Map3(E2, I2, I2),
	III  = Map3(I2, I2, I2),
	// E/S * 3  (full arrangement)
	SEE  = Map3(S2, E2, E2),
	ESE  = Map3(E2, S2, E2),
	SSE  = Map3(S2, S2, E2),
	EES  = Map3(E2, E2, S2),
	SES  = Map3(S2, E2, S2),
	ESS  = Map3(E2, S2, S2),
	SSS  = Map3(S2, S2, S2),
	// #################################
	EEEE = Map4(E2, E2, E2, E2)

	// remains are not implemented, we need auto generator... :(
};

enum class PntArr3 : size_t
{
	E     = E3,
	I     = I3,
	S     = S3,
	L     = L3,
	T     = T3,
	// #################################
	EE    = Map2(E3, E3),
	// E/I * 2  (full arrangement)
	IE    = Map2(I3, E3),
	EI    = Map2(E3, I3),
	II    = Map2(I3, I3),
	// E/S/L/T * 2  (not full arrangement)
	SS    = Map2(S3, S3),
	LS    = Map2(L3, S3),
	TS    = Map2(T3, S3),
	LL    = Map2(L3, L3),
	TL    = Map2(T3, L3),
	TT    = Map2(T3, T3),
	// #################################
	EEE   = Map3(E3, E3, E3),
	// E/I * 3  (full arrangement)
	IEE   = Map3(I3, E3, E3),
	EIE   = Map3(E3, I3, E3),
	IIE   = Map3(I3, I3, E3),
	EEI   = Map3(E3, E3, I3),
	IEI   = Map3(I3, E3, I3),
	EII   = Map3(E3, I3, I3),
	III   = Map3(I3, I3, I3),
	// E/S/L/T * 3  (not full arrangement)
	SSS   = Map3(S3, S3, S3),
	LSS   = Map3(L3, S3, S3),
	TSS   = Map3(T3, S3, S3),
	LLS   = Map3(L3, L3, S3),
	TLS   = Map3(T3, L3, S3),
	TTS   = Map3(T3, T3, S3),
	LLL   = Map3(L3, L3, L3),
	TLL   = Map3(T3, L3, L3),
	TTL   = Map3(T3, T3, L3),
	TTT   = Map3(T3, T3, T3),
	// #################################
	// E * 4
	EEEE  = Map4(E3, E3, E3, E3),
	// E/I * 4  (full arrangement)
	IEEE  = Map4(I3, E3, E3, E3),
	EIEE  = Map4(E3, I3, E3, E3),
	IIEE  = Map4(I3, I3, E3, E3),
	EEIE  = Map4(E3, E3, I3, E3),
	IEIE  = Map4(I3, E3, I3, E3),
	EIIE  = Map4(E3, I3, I3, E3),
	IIIE  = Map4(I3, I3, I3, E3),
	EEEI  = Map4(E3, E3, E3, I3),
	IEEI  = Map4(I3, E3, E3, I3),
	EIEI  = Map4(E3, I3, E3, I3),
	IIEI  = Map4(I3, I3, E3, I3),
	EEII  = Map4(E3, E3, I3, I3),
	IEII  = Map4(I3, E3, I3, I3),
	EIII  = Map4(E3, I3, I3, I3),
	IIII  = Map4(I3, I3, I3, I3),
	// E/S/L/T * 4 (not full arrangement)
	// 4
	SSSS  = Map4(S3, S3, S3, S3),
	LSSS  = Map4(L3, S3, S3, S3),
	TSSS  = Map4(T3, S3, S3, S3),
	LLSS  = Map4(L3, L3, S3, S3),
	TLSS  = Map4(T3, L3, S3, S3),
	TTSS  = Map4(T3, T3, S3, S3),
	LLLS  = Map4(L3, L3, L3, S3),
	TLLS  = Map4(T3, L3, L3, S3),
	TTLS  = Map4(T3, T3, L3, S3),
	TTTS  = Map4(T3, T3, T3, S3),
	LLLL  = Map4(L3, L3, L3, L3),
	TLLL  = Map4(T3, L3, L3, L3),
	TTLL  = Map4(T3, T3, L3, L3),
	TTTL  = Map4(T3, T3, T3, L3),
	TTTT  = Map4(T3, T3, T3, T3),
	// #################################
	// E * 5
	EEEEE = Map5(E3, E3, E3, E3, E3),
	// E/I * 5 (not full arrangement)
	IEEEE = Map5(I3, E3, E3, E3, E3),
	IIEEE = Map5(I3, I3, E3, E3, E3),
	IIIEE = Map5(I3, I3, I3, E3, E3),
	IIIIE = Map5(I3, I3, I3, I3, E3),
	IIIII = Map5(I3, I3, I3, I3, I3)
};

template <int N>
void sort_types(std::array<uint32_t, N> &pos, std::array<uint32_t, N> &types,
                uint32_t &swap_cnt)
{
	// a simple bubble sort
	swap_cnt = 0;
	for (int i = 1; i < N; ++i)
	{
		uint32_t p = pos[i];
		uint32_t t = types[i];
		if (t == E3)
			continue;

		int j = i - 1;
		while (j >= 0 && t > types[j])
		{
			pos[j + 1]   = pos[j];
			types[j + 1] = types[j];
			++swap_cnt;
			--j;
		}
		pos[j + 1]   = p;
		types[j + 1] = t;
	}
}

template <bool ToEI>
PntArr2 get_pnts_arr2(uint32_t p0, uint32_t p1)
{
	if constexpr (ToEI)
	{
		p0 = p0 > I2 ? I2 : E2;
		p1 = p1 > I2 ? I2 : E2;
	}
	return static_cast<PntArr2>(Map2(p0, p1));
}

template <bool ToEI>
PntArr2 get_pnts_arr2(uint32_t p0, uint32_t p1, uint32_t p2)
{
	if constexpr (ToEI)
	{
		p0 = p0 > I2 ? I2 : E2;
		p1 = p1 > I2 ? I2 : E2;
		p2 = p2 > I2 ? I2 : E2;
	}
	return static_cast<PntArr2>(Map3(p0, p1, p2));
}

template <bool ToEI>
PntArr2 get_pnts_arr2(uint32_t p0, uint32_t p1, uint32_t p2, uint32_t p3)
{
	if constexpr (ToEI)
	{
		p0 = p0 > I2 ? I2 : E2;
		p1 = p1 > I2 ? I2 : E2;
		p2 = p2 > I2 ? I2 : E2;
		p3 = p3 > I2 ? I2 : E2;
	}
	return static_cast<PntArr2>(Map4(p0, p1, p2, p3));
}

template <bool ToEI>
PntArr3 get_pnts_arr3(uint32_t p0, uint32_t p1)
{
	if constexpr (ToEI)
	{
		p0 = p0 > I3 ? I3 : E3;
		p1 = p1 > I3 ? I3 : E3;
	}
	return static_cast<PntArr3>(Map2(p0, p1));
}

template <bool ToEI>
PntArr3 get_pnts_arr3(uint32_t p0, uint32_t p1, uint32_t p2)
{
	if constexpr (ToEI)
	{
		p0 = p0 > I3 ? I3 : E3;
		p1 = p1 > I3 ? I3 : E3;
		p2 = p2 > I3 ? I3 : E3;
	}
	return static_cast<PntArr3>(Map3(p0, p1, p2));
}

template <bool ToEI>
PntArr3 get_pnts_arr3(uint32_t p0, uint32_t p1, uint32_t p2, uint32_t p3)
{
	if constexpr (ToEI)
	{
		p0 = p0 > I3 ? I3 : E3;
		p1 = p1 > I3 ? I3 : E3;
		p2 = p2 > I3 ? I3 : E3;
		p3 = p3 > I3 ? I3 : E3;
	}
	return static_cast<PntArr3>(Map4(p0, p1, p2, p3));
}

template <bool ToEI>
PntArr3 get_pnts_arr3(uint32_t p0, uint32_t p1, uint32_t p2, uint32_t p3,
                      uint32_t p4)
{
	if constexpr (ToEI)
	{
		p0 = p0 > I3 ? I3 : E3;
		p1 = p1 > I3 ? I3 : E3;
		p2 = p2 > I3 ? I3 : E3;
		p3 = p3 > I3 ? I3 : E3;
		p4 = p4 > I3 ? I3 : E3;
	}
	return static_cast<PntArr3>(Map5(p0, p1, p2, p3, p4));
}

/// @brief get two points' arrangement in 3D
inline PntArr2 sort_pnts_arr2(std::array<uint32_t, 2> &types,
                              std::array<uint32_t, 2> &pos, uint32_t &swap_cnt)
{
	// sort points
	sort_types<2>(pos, types, swap_cnt);
	return static_cast<PntArr2>(Map2(types[0], types[1]));
}

/// @brief get three points' arrangement in 3D
inline PntArr2 sort_pnts_arr2(std::array<uint32_t, 3> &types,
                              std::array<uint32_t, 3> &pos, uint32_t &swap_cnt)
{
	// sort points
	sort_types<3>(pos, types, swap_cnt);
	return static_cast<PntArr2>(Map3(types[0], types[1], types[2]));
}

/// @brief get two points' arrangement in 3D
inline PntArr3 sort_pnts_arr3(std::array<uint32_t, 2> &types,
                              std::array<uint32_t, 2> &pos, uint32_t &swap_cnt)
{
	// sort points
	sort_types<2>(pos, types, swap_cnt);
	return static_cast<PntArr3>(Map2(types[0], types[1]));
}

/// @brief get three points' arrangement in 3D
inline PntArr3 sort_pnts_arr3(std::array<uint32_t, 3> &types,
                              std::array<uint32_t, 3> &pos, uint32_t &swap_cnt)
{
	// sort points
	sort_types<3>(pos, types, swap_cnt);
	return static_cast<PntArr3>(Map3(types[0], types[1], types[2]));
}

/// @brief get four points' arrangement in 3D
inline PntArr3 sort_pnts_arr3(std::array<uint32_t, 4> &types,
                              std::array<uint32_t, 4> &pos, uint32_t &swap_cnt)
{
	// sort points
	sort_types<4>(pos, types, swap_cnt);
	return static_cast<PntArr3>(Map4(types[0], types[1], types[2], types[3]));
}

/// @brief get five points' arrangement in 3D
inline PntArr3 sort_pnts_arr3(std::array<uint32_t, 5> &types,
                              std::array<uint32_t, 5> &pos, uint32_t &swap_cnt)
{
	// sort points
	sort_types<5>(pos, types, swap_cnt);
	return static_cast<PntArr3>(Map5(types[0], types[1], types[2], types[3], types[4]));
}

#undef E2
#undef I2
#undef S2
#undef E3
#undef I3
#undef S3
#undef L3
#undef T3
#undef Map2
#undef Map3
#undef Map4
#undef Map5
#undef BEP

#ifdef OMC_PRED_PROFILE

// === profile control
// #define OMC_PRED_PROFILE_FILTER
// #define OMC_PRED_PROFILE_IPOINT
// #define OMC_PRED_PROFILE_LENGTH
// #define OMC_PRED_PROFILE_COMPRS

enum class PredicateNames : size_t
{
	// below use filter,ss_fail,d_fail,real_zero count
	_lessThanOnX_IE = 0,
	_lessThanOnX_II,
	_lessThanOnY_IE,
	_lessThanOnY_II,
	_lessThanOnZ_IE,
	_lessThanOnZ_II,
	_orientOn2Dxy_IEE,
	_orientOn2Dxy_IIE,
	_orientOn2Dxy_III,
	_orientOn2Dyz_IEE,
	_orientOn2Dyz_IIE,
	_orientOn2Dyz_III,
	_orientOn2Dzx_IEE,
	_orientOn2Dzx_IIE,
	_orientOn2Dzx_III,
	_orient3D_IEEE,
	// below use total count
	_ssi_filter,
	_ssi_interval,
	_ssi_expansion,
	_lpi_filter,
	_lpi_interval,
	_lpi_expansion,
	_tpi_filter,
	_tpi_interval,
	_tpi_expansion,
	// below use total, branch count
	_orientOn2D_III,
	// below use total, branch count
	_expan_reducable_len,
	_expan_reduced_len,
	CNT
};

struct PredicatesProfile
{
	static constexpr size_t ARR_CNT    = 32;
	static constexpr size_t MAG_CNT    = 100;
	static constexpr size_t BRANCH_CNT = 256;

	// clang-format off
	static std::atomic_size_t filter_count[(size_t)PredicateNames::CNT][ARR_CNT];
	static std::atomic_size_t ss_fail_count[(size_t)PredicateNames::CNT][ARR_CNT];
	static std::atomic_size_t d_fail_count[(size_t)PredicateNames::CNT][ARR_CNT];
	static std::atomic_size_t real_zero_count[(size_t)PredicateNames::CNT][ARR_CNT];

	static std::atomic_size_t total_count[(size_t)PredicateNames::CNT];

	static std::atomic_size_t branch_count[(size_t)PredicateNames::CNT][BRANCH_CNT];
	// clang-format on

	static void initialize();
	static void inc_filter(PredicateNames name, PntArr3 arr);
	static void inc_ss_fail(PredicateNames name, PntArr3 arr);
	static void inc_d_fail(PredicateNames name, PntArr3 arr);
	static void inc_real_zero(PredicateNames name, PntArr3 arr);

	static void inc_total(PredicateNames name, size_t count);

	static void inc_branch(PredicateNames name, size_t branch);

	static void print();
};

	#define OMC_PRED_PROFILE_INIT OMC::PredicatesProfile::initialize()
	#define OMC_PRED_PROFILE_PRINT OMC::PredicatesProfile::print()

	#ifdef OMC_PRED_PROFILE_FILTER
		#define OMC_PRED_PROFILE_INC_FILTER(pred, arr) \
			OMC::PredicatesProfile::inc_filter(pred, arr)
		#define OMC_PRED_PROFILE_INC_SSFAIL(pred, arr) \
			OMC::PredicatesProfile::inc_ss_fail(pred, arr)
		#define OMC_PRED_PROFILE_INC_DFAIL(pred, arr) \
			OMC::PredicatesProfile::inc_d_fail(pred, arr)
		#define OMC_PRED_PROFILE_INC_REALZERO(ret, pred, arr) \
			if (ret == Sign::ZERO)                              \
			{                                                   \
				OMC::PredicatesProfile::inc_real_zero(pred, arr); \
			}
	#else
		#define OMC_PRED_PROFILE_INC_FILTER(pred, arr)
		#define OMC_PRED_PROFILE_INC_SSFAIL(pred, arr)
		#define OMC_PRED_PROFILE_INC_DFAIL(pred, arr)
		#define OMC_PRED_PROFILE_INC_REALZERO(ret, pred, arr)
	#endif

	#ifdef OMC_PRED_PROFILE_IPOINT
		#define OMC_PRED_PROFILE_INC_IP_TOTAL(name) \
			OMC::PredicatesProfile::inc_total(name, 1)
	#else
		#define OMC_PRED_PROFILE_INC_IP_TOTAL(name)
	#endif

	#ifdef OMC_PRED_PROFILE_LENGTH
		#define OMC_PRED_PROFILE_INC_LEN_TOTAL(name, count) \
			OMC::PredicatesProfile::inc_total(name, count)
		#define OMC_PRED_PROFILE_INC_LEN(name, branch) \
			OMC::PredicatesProfile::inc_branch(name, branch)
	#else
		#define OMC_PRED_PROFILE_INC_LEN_TOTAL(name, count)
		#define OMC_PRED_PROFILE_INC_LEN(name, branch)
	#endif

	#ifdef OMC_PRED_PROFILE_COMPRS
		#define OMC_PRED_PROFILE_SAVE_COMPRS(orig_len, new_len)                    \
			if ((new_len) < (orig_len))                                              \
			{                                                                        \
				OMC::PredicatesProfile::inc_total(                                     \
				  PredicateNames::_expan_reducable_len, 1);                            \
				OMC::PredicatesProfile::inc_total(PredicateNames::_expan_reduced_len,  \
				                                  1);                                  \
				OMC::PredicatesProfile::inc_branch(                                    \
				  PredicateNames::_expan_reducable_len, orig_len);                     \
				OMC::PredicatesProfile::inc_branch(PredicateNames::_expan_reduced_len, \
				                                   (orig_len) - (new_len));            \
			}
	#else
		#define OMC_PRED_PROFILE_SAVE_COMPRS(orig_len, new_len)
	#endif

#else

	#define OMC_PRED_PROFILE_INIT
	#define OMC_PRED_PROFILE_PRINT

	#define OMC_PRED_PROFILE_INC_FILTER(pred, arr)
	#define OMC_PRED_PROFILE_INC_SSFAIL(pred, arr)
	#define OMC_PRED_PROFILE_INC_DFAIL(pred, arr)
	#define OMC_PRED_PROFILE_INC_REALZERO(ret, pred, arr)

	#define OMC_PRED_PROFILE_INC_IP_TOTAL(name)

	#define OMC_PRED_PROFILE_INC_LEN_TOTAL(name, count)
	#define OMC_PRED_PROFILE_INC_LEN(name, branch)

	#define OMC_PRED_PROFILE_SAVE_COMPRS(orig_len, new_len)
#endif

} // namespace OMC