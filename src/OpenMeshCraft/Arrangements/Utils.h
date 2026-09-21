#pragma once

// clang-format off
#include "OpenMeshCraft/Utils/DisableWarnings.h"
#include "tbb/tbb.h"
#include "OpenMeshCraft/Utils/EnableWarnings.h"
// clang-format on

#include "OpenMeshCraft/Utils/FormatCompat.h"
#include "OpenMeshCraft/Utils/IndexDef.h"
#include "OpenMeshCraft/Utils/Label.h"
#include "OpenMeshCraft/Utils/Logger.h"

#include <algorithm>
#include <deque>
#include <queue>
#include <vector>

namespace OMC {

// === parallel control for each part
#define OMC_ARR_TS_PARA
#define OMC_ARR_DC_TTI_PARA
#define OMC_ARR_TR_PARA
// === accelerate control for components
#define OMC_ARR_DC_FILTER_O3D
#define OMC_ARR_AVOID_TPI
// === ablation control for features
// #define OMC_ARR_GLOBAL_POINT_SET
// #define OMC_ARR_AUX_LPI
// #define OMC_ARR_3D_PREDS

struct MeshArrangements_Config
{
  /*   behavior   */
  // report more information during running
  bool   verbose                = false;
  // If set to true, the algorithm will ignore intersections between triangles
  // in the same mesh. This feature is utilized by mesh boolean operations.
  bool   ignore_same_mesh       = false;
  // If set to true, the explicit results (points and triangles) will be saved
  // in the output mesh specified by setTriMeshAsOutput. All intermediate data
  // will be cleared.
  bool   output_explicit_result = false;
  /*  tree parameters  */
  double tree_enlarge_ratio     = 1.01;
  size_t tree_split_size_thres  = 400;
  /*  only for OcTree  */
  double tree_adaptive_thres    = 0.2;
};

struct MeshArrangements_Stats
{
  /* Preprocess and build tree *************************************/

  double pp_elapsed   = 0.; // timings of preprocessing
  double tree_elapsed = 0.; // timings of building tree

  /* detect and classify intersection ******************************/

  double ci_elapsed = 0.; // timings of detection and classification

  /* triangulation *************************************************/

  double tr_elapsed = 0.; // timings of triangulation
};

/// @brief Each output triangle has a surface label and inside label.
struct ArrLabels
{
  /// surface (mesh) that a triangle belongs to.
  std::vector<Label> surface;
  /// surface (mesh) that a triangle is inside.
  std::vector<Label> inside;
  /// how many surfaces (meshes) in arrangements.
  size_t             num;
};

template <typename Traits>
struct ArrPointArena
{
public:
  using EPoint3      = typename Traits::EPoint3;
  using GPoint3      = typename Traits::GPoint3;
  using IPoint3T_SSI = typename Traits::IPoint3T_SSI;
  using IPoint3T_LPI = typename Traits::IPoint3T_LPI;
  using IPoint3T_TPI = typename Traits::IPoint3T_TPI;

public:
  std::deque<IPoint3T_SSI> ssi; // SSI points
  std::deque<IPoint3T_LPI> lpi; // LPI points
  std::deque<IPoint3T_TPI> tpi; // TPI points

public:
  void recycle(IPoint3T_SSI *ssi_ptr) { recycled_ssi.push(ssi_ptr); }
  void recycle(IPoint3T_LPI *lpi_ptr) { recycled_lpi.push(lpi_ptr); }
  void recycle(IPoint3T_TPI *tpi_ptr) { recycled_tpi.push(tpi_ptr); }

  IPoint3T_SSI *emplace(IPoint3T_SSI &&p)
  {
    IPoint3T_SSI *ssi_ptr = nullptr;
    if (recycled_ssi.empty())
    {
      ssi.emplace_back(std::move(p));
      ssi_ptr = &ssi.back();
    }
    else
    {
      ssi_ptr = recycled_ssi.front();
      recycled_ssi.pop();
      *ssi_ptr = std::move(p);
    }
    return ssi_ptr;
  }

  IPoint3T_LPI *emplace(IPoint3T_LPI &&p)
  {
    IPoint3T_LPI *lpi_ptr = nullptr;
    if (recycled_lpi.empty())
    {
      lpi.emplace_back(std::move(p));
      lpi_ptr = &lpi.back();
    }
    else
    {
      lpi_ptr = recycled_lpi.front();
      recycled_lpi.pop();
      *lpi_ptr = std::move(p);
    }
    return lpi_ptr;
  }

  IPoint3T_TPI *emplace(IPoint3T_TPI &&p)
  {
    IPoint3T_TPI *tpi_ptr = nullptr;
    if (recycled_tpi.empty())
    {
      tpi.emplace_back(std::move(p));
      tpi_ptr = &tpi.back();
    }
    else
    {
      tpi_ptr = recycled_tpi.front();
      recycled_tpi.pop();
      *tpi_ptr = std::move(p);
    }
    return tpi_ptr;
  }

  void reserve_ssi(size_t new_n)
  {
    for (size_t i = 0; i < new_n; i++)
      recycle(&ssi.emplace_back());
  }

  void reserve_lpi(size_t new_n)
  {
    for (size_t i = 0; i < new_n; i++)
      recycle(&lpi.emplace_back());
  }

  void reserve_tpi(size_t new_n)
  {
    for (size_t i = 0; i < new_n; i++)
      recycle(&tpi.emplace_back());
  }

private:
  std::queue<IPoint3T_SSI *> recycled_ssi;
  std::queue<IPoint3T_LPI *> recycled_lpi;
  std::queue<IPoint3T_TPI *> recycled_tpi;
};

struct ArrDuplTriInfo
{
  index_t t_id; // triangle id
  index_t l_id; // label id
  bool    w;    // winding (CW or CCW)
};

} // namespace OMC

#ifdef OMC_ARR_PROFILE

// === profile control
// #define OMC_ARR_PROF_TTI
// #define OMC_ARR_PROF_MAXVAR
// #define OMC_ARR_PROF_TPI_LENGTH
// #define OMC_ARR_PROF_ES_PNTS
// #define OMC_ARR_PROF_TREE_NODE

  #include <fstream>

namespace OMC {

enum class ArrFuncNames : size_t
{
  D_BBI = 0,
  D_BBI_SMALL,
  D_BBI_LARGE,
  DC_TTI,
  IP_CNT,
  IP_MAXVAR_ORDER,
  SSI_EXP_LEN,
  LPI_EXP_LEN,
  TPI_EXP_LEN,
  CNT
};

struct ArrProfile
{
  static const uint32_t BRANCH_CNT = 256;

  static std::atomic_size_t total_count[(size_t)ArrFuncNames::CNT];
  static std::atomic_size_t reach_count[(size_t)ArrFuncNames::CNT][BRANCH_CNT];
  static std::atomic_size_t reach_line[(size_t)ArrFuncNames::CNT][BRANCH_CNT];

  static inline void initialize();
  static inline void inc_total(ArrFuncNames name);
  static inline void inc_reach(ArrFuncNames name, uint32_t branch_flag,
                               uint32_t branch_line);

  static inline void inc_total(ArrFuncNames name, size_t count);
  static inline void inc_reach(ArrFuncNames name, uint32_t branch_flag,
                               uint32_t branch_line, size_t count);

  static inline void print();
};

inline void ArrProfile::initialize()
{
  for (size_t i = 0; i < (size_t)ArrFuncNames::CNT; i++)
  {
    total_count[i] = 0;
    for (size_t j = 0; j < BRANCH_CNT; j++)
    {
      reach_count[i][j] = 0;
    }
  }
}

inline void ArrProfile::inc_total(ArrFuncNames name)
{
  total_count[(size_t)name] += 1;
}

inline void ArrProfile::inc_reach(ArrFuncNames name, uint32_t branch_flag,
                                  uint32_t branch_line)
{
  reach_count[(size_t)name][branch_flag] += 1;
  reach_line[(size_t)name][branch_flag] = branch_line;
}

inline void ArrProfile::inc_total(ArrFuncNames name, size_t count)
{
  total_count[(size_t)name] += count;
}

inline void ArrProfile::inc_reach(ArrFuncNames name, uint32_t branch_flag,
                                  uint32_t branch_line, size_t count)
{
  reach_count[(size_t)name][branch_flag] += count;
  reach_line[(size_t)name][branch_flag] = branch_line;
}

inline void ArrProfile::print()
{
  // clang-format off
  std::vector<std::string> func_names = {
    "Detect BBI",
    "Detect BBI small node",
    "Detect BBI large node",
    "Detect & Classify TTI",
    "Implicit Points Count",
    "Implicit Points MaxVar Order",
    "SSI Expansion Length",
    "LPI Expansion Length",
    "TPI Expansion Length"
  };

  // clang-format on

  for (size_t i = 0; i < (size_t)ArrFuncNames::CNT; i++)
  {
    std::cout << OMC::format("{}:\n", func_names[i]);
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

    for (int j = 0; j <= last_branch_flag; j++)
    {
      double reach_raio = (double)reach_count[i][j] / (double)total_count[i];
      std::cout << OMC::format("  line {}, branch {}: {:.2f}%, {}\n",
                               reach_line[i][j].load(), j, reach_raio * 100.,
                               reach_count[i][j].load());
    }
  }
}
} // namespace OMC

  #define OMC_ARR_PROFILE_INIT OMC::ArrProfile::initialize()
  #define OMC_ARR_PROFILE_PRINT OMC::ArrProfile::print()

  #define OMC_ARR_PROFILE_INC_TOTAL(func) OMC::ArrProfile::inc_total(func)

  #define OMC_ARR_PROFILE_INC_REACH(func, branch_flag) \
    OMC::ArrProfile::inc_reach(func, branch_flag, __LINE__)

  #define OMC_ARR_PROFILE_INC_TOTAL_CNT(func, count) \
    OMC::ArrProfile::inc_total(func, count)

  #define OMC_ARR_PROFILE_INC_REACH_CNT(func, branch_flag, count) \
    OMC::ArrProfile::inc_reach(func, branch_flag, __LINE__, count)

#else

  #define OMC_ARR_PROFILE_INIT
  #define OMC_ARR_PROFILE_PRINT

  #define OMC_ARR_PROFILE_INC_TOTAL(func)
  #define OMC_ARR_PROFILE_INC_REACH(func, branch_flag)

  #define OMC_ARR_PROFILE_INC_TOTAL_CNT(func, count)
  #define OMC_ARR_PROFILE_INC_REACH_CNT(func, branch_flag, count)

#endif

#define OMC_ARR_START_ELAPSE(name) auto name = OMC::Logger::elapse_reset();

#define OMC_ARR_SAVE_ELAPSED(name, dst_name, dscrpt)                     \
  if (stats != nullptr)                                                  \
    stats->dst_name = OMC::Logger::elapsed(name).count();                \
  if (config.verbose)                                                    \
  {                                                                      \
    OMC::Logger::info(OMC::format("[OpenMeshCraft Arrangements] " dscrpt \
                                  " time : {} s",                        \
                                  OMC::Logger::elapsed(name).count()));  \
  }
