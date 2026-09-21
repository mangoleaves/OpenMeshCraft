#pragma once

#include "OpenMeshCraft/Utils/FormatCompat.h"

#include <deque>
#include <queue>

namespace OMC {

// Enable shuffling missing segments in each loop of segment recovery.
// #define OMC_CDT_SHUFFLE_MISSING_SEGMENTS

// Enable exact inSphere predicate in finding encroaching point.
#define OMC_CDT_EXACT_ENCROACH_TEST

struct ConstrDelTet_Config
{
  /// output more log information
  bool verbose                = false;
  /// If set to true, the explicit result (points and tetrahedra) will be saved
  /// in output mesh set by setTriMeshAsOutput.
  bool output_explicit_result = false;

  /// The threshold of the number of Steiner points inserted in segment
  /// recovery. The threshold is a multiple of the number of vertices.
  /// Set to 0 to disable the threshold.
  size_t Steiner_point_thres = 0;
};

struct ConstrDelTet_Stats
{
  /* Preprocessing *************************************************/

  double pp_elapsed = 0.; // timings of preprocessing

  /* Delaunay tetrahedralization ***********************************/

  double dt_elapsed = 0.; // timings of Delaunay tetrahedralization

  /* Constraints Recovery ******************************************/

  double seg_elapsed  = 0.; // timings of segment recovery
  double face_elapsed = 0.; // timings of face recovery

  size_t seg_steiner = 0; // number of Steiner points
};

template <typename Traits>
struct CDTPointArena
{
public:
  using EPoint3      = typename Traits::EPoint3;
  using GPoint3      = typename Traits::GPoint3;
  using IPoint3T_LNC = typename Traits::IPoint3T_LNC;

public:
  std::deque<IPoint3T_LNC> lnc; // LNC points

public:
  void recycle(IPoint3T_LNC *lnc_ptr) { recycled_lnc.push(lnc_ptr); }

  IPoint3T_LNC *emplace(IPoint3T_LNC &&p)
  {
    IPoint3T_LNC *lnc_ptr = nullptr;
    if (recycled_lnc.empty())
    {
      lnc.emplace_back(std::move(p));
      lnc_ptr = &lnc.back();
    }
    else
    {
      lnc_ptr = recycled_lnc.front();
      recycled_lnc.pop();
      *lnc_ptr = std::move(p);
    }
    return lnc_ptr;
  }

  void reserve_lnc(size_t new_n)
  {
    for (size_t i = 0; i < new_n; i++)
      recycle(&lnc.emplace_back());
  }

private:
  std::queue<IPoint3T_LNC *> recycled_lnc;
};

#define OMC_CDT_START_ELAPSE(name) auto name = OMC::Logger::elapse_reset();

#define OMC_CDT_SAVE_ELAPSED(name, dst_name, dscrpt)                    \
  if (stats != nullptr)                                                 \
    stats->dst_name = OMC::Logger::elapsed(name).count();               \
  if (config.verbose)                                                   \
  {                                                                     \
    OMC::Logger::info(OMC::format("[OpenMeshCraft CDT] " dscrpt         \
                                  " time : {} s",                       \
                                  OMC::Logger::elapsed(name).count())); \
  }

} // namespace OMC
