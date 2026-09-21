#pragma once

#include "SegmentRecover.h"

#include "Utils.h"

#include "boost/container/flat_set.hpp"
#include "OpenMeshCraft/Utils/FormatCompat.h"

#include <iostream>

namespace OMC {

/**
 * @brief Initialize with the Delaunay tetrahedral mesh and the input
 * constrained triangles. Prepare for the constraits recovery.
 * @param _tet_mesh the Delaunay tetrahedral mesh
 * @param _plc the input constrained piecewise linear complex
 */
template <typename Traits>
SegmentRecover<Traits>::SegmentRecover(std::vector<GPoint3 *> &_verts,
                                       std::vector<PntArena>  &_pnt_arenas,
                                       TetMesh &_tet_mesh, PLC &_plc,
                                       ConstrDelTet_Config _config,
                                       ConstrDelTet_Stats *_stats)
  : verts(_verts)
  , pnt_arenas(_pnt_arenas)
  , tet_mesh(_tet_mesh)
  , plc(_plc)
  , config(_config)
  , stats(_stats)
{
  // initialize states
  protecting_sphere_initialized = false;
}

template <typename Traits>
void SegmentRecover<Traits>::segmentRecovery()
{
  segmentRecovery_SiHang();

#if defined(OMC_ENABLE_EXPENSIVE_ASSERT)
  // check if all segments are recovered
  for (index_t eid = 0; eid < plc.numEdges(); eid++)
  {
    const PLCEdge &e = plc.edge(eid);
    if (e.isConstraint())
    {
      OMC_ASSERT(tet_mesh.edgeExists(e.ep0(), e.ep1()), "Missing segment {}.",
                 eid);
    }
  }
#endif
}

/**
 * @brief Create a new vertex and update `tet_mesh` and `plc` at the same time.
 * @param new_pnt geometric implementation of the new vertex
 * @return the index to the new vertex.
 */
template <typename Traits>
template <typename PointType>
index_t SegmentRecover<Traits>::newVtx(PointType new_pnt)
{
  // create the new vertex
  index_t new_vid     = verts.size();
  // Put the new vertex into the point arena
  auto   *new_pnt_ptr = pnt_arenas[0].emplace(std::move(new_pnt));
  // Put the new vertex into the vertex list
  verts.push_back(static_cast<GPoint3 *>(new_pnt_ptr));

  // create auxiliary data in TetMesh & PLC
  tet_mesh.newVtx(new_vid);
  plc.newVtx(new_vid);

  return new_vid;
}

/**
 * @brief segment recovery using the Si-Hang method.
 * @param num_loop the number of loops to run the segment recovery.
 * If `num_loop` is 0, the segment recovery will run until all segments are
 * recovered.
 * @pre PLC are well constructed, and the tetrahedral mesh is Delaunay.
 */
template <typename Traits>
void SegmentRecover<Traits>::segmentRecovery_SiHang(size_t num_loop)
{
  // traverse all edges in the PLC to find missing segments
  std::vector<index_t> missing_segments;
  for (index_t eid = 0; eid < plc.numEdges(); eid++)
  {
    const PLCEdge &e = plc.edge(eid);
    if (e.isConstraint() && !tet_mesh.edgeExists(e.ep0(), e.ep1()))
    {
      missing_segments.push_back(eid);
    }
  }

  // Once a segment is split, new missing segments may appear near it.
  // So, mark adjacent vertices as `TO_CHECK`, and check them after a loop.
  for (index_t vid = 0; vid < tet_mesh.sizeVerts(); vid++)
  {
    tet_mesh.unmark(vid, VTX_MARK::TO_CHECK);
    tet_mesh.unmark(vid, VTX_MARK::ENCROACHED);
  }

  size_t orig_vn     = plc.input_nv;
  size_t split_count = 0;
  size_t loop_count  = 0;

  // Main loop: split missing segments to recover them.
  while (!missing_segments.empty())
  {
    // Split the existing constrained segment
    while (!missing_segments.empty())
    {
      // Get the constrained segment to split
      index_t eid = missing_segments.back();
      missing_segments.pop_back();
      const PLCEdge &e   = plc.edge(eid);
      index_t        ep0 = e.ep0(), ep1 = e.ep1();

      // Check if the segment is still missing
      if (tet_mesh.edgeExists(ep0, ep1))
        continue;

      // Split the segment
      index_t new_vid = splitMissingSegment(eid);

      // Visit the neighbor vertices
      // New missing edges will appear near visited vertices
      tet_mesh.mark(ep0, VTX_MARK::TO_CHECK);
      tet_mesh.mark(ep1, VTX_MARK::TO_CHECK);
      tet_mesh.mark(new_vid, VTX_MARK::TO_CHECK);
      InlinedVector64<index_t> local_vv;
      tet_mesh.VV(new_vid, local_vv);
      for (index_t vid : local_vv)
        tet_mesh.mark(vid, VTX_MARK::TO_CHECK);

      // log and output
      split_count++;
      if (config.verbose && split_count % 100 == 0)
      {
        std::cout << OMC::format(
          "\r[OpenMeshCraft CDT] {} segments are split. "
          "{} segments are missing.                     ",
          split_count, missing_segments.size());
      }
    }
    loop_count++;
    if (num_loop > 0 && loop_count == num_loop)
      break;
    // Find new missing edges around visited vertices
    for (index_t eid = 0; eid < plc.numEdges(); eid++)
    {
      const PLCEdge &e = plc.edge(eid);
      if (e.isConstraint() &&
          (tet_mesh.isMarked(e.ep0(), VTX_MARK::TO_CHECK) ||
           tet_mesh.isMarked(e.ep1(), VTX_MARK::TO_CHECK)) &&
          !tet_mesh.edgeExists(e.ep0(), e.ep1()))
      {
        missing_segments.push_back(eid);
      }
    }
    OMC_ASSERT(config.Steiner_point_thres == 0 ||
                 split_count < orig_vn * config.Steiner_point_thres,
               "Too many Steiner points are inserted.");
#ifdef OMC_CDT_SHUFFLE_MISSING_SEGMENTS
    std::shuffle(missing_segments.begin(), missing_segments.end(),
                 std::default_random_engine(std::random_device()()));
#endif
    // clear the `TO_CHECK` mark for all vertices
    for (index_t vid = 0; vid < tet_mesh.sizeVerts(); vid++)
      tet_mesh.unmark(vid, VTX_MARK::TO_CHECK);
  }
  if (config.verbose) // output a new line
    std::cout << std::endl;

  tet_mesh.removeDeletedTets();
}

template <typename Traits>
index_t SegmentRecover<Traits>::splitMissingSegment(index_t eid)
{
  IPoint3T_LNC new_pnt;
  index_t      curr_tet = InvalidIndex;

  const PLCEdge &edge = plc.edge(eid);

  if (edge.type == PLCEdgeType::BOTH_ACUTE_VERTEX)
  {
    // This split strategy is used to split segment with two acute vertices, and
    // generating two sub-edges with type `ONE_ACUTE_VERTEX`.
    curr_tet = TetMesh::toIdOff(tet_mesh.incTet(edge.ep0()));
    new_pnt  = middlePoint(edge);
  }
  else // ONE_ACUTE_VERTEX or NO_ACUTE_VERTEX
  {
    index_t ref_vid;
    findReferenceEncroachingPoint(eid, ref_vid, curr_tet);

#ifdef OMC_CDT_EXACT_ENCROACH_TEST
    OMC_EXPENSIVE_ASSERT(is_valid_idx(ref_vid),
                         "Could not find a valid reference encroaching point.");
#else
    if (!is_valid_idx(ref_vid))
    { //  can't find the reference encroaching point due to inexact predicate.
      // just (use BothAcute strategy to) split the segment at its middle point.
      curr_tet = TetMesh::toIdOff(tet_mesh.incTet(edge.ep0()));
      new_pnt  = middlePoint(edge);
    }
    else
#endif
    {
      if (edge.type == PLCEdgeType::NO_ACUTE_VERTEX)
        new_pnt = splitSegment_NoAcuteVertex(eid, ref_vid);
      else // ONE_ACUTE_VERTEX
        new_pnt = splitSegment_OneAcuteVertex(eid, ref_vid);
    }
  }

  // Add the new point
  index_t new_vid = newVtx(new_pnt);
  // Split edge by the new point
  plc.splitPLCEdge(eid, new_vid);
  // insert the splitting point into the Delaunay tetrahedral mesh
  DelTet DT(tet_mesh);
  DT.insertVertex(new_vid, curr_tet);

  OMC_EXPENSIVE_ASSERT(DT.localVerify(new_vid),
                       "Invalid Delaunay tetrahedralization.");
  return new_vid;
}

/**
 * @brief Find reference encroaching point for segment `eid`.
 *
 * The segment has two endpoints `v1` and `v2`.
 * The segment defines a diametral sphere `D`.
 * Encroaching points are points enclosed (or visited?) by `D`.
 * Reference encroaching point is the point `r` in encroaching points, such
 * that `v1`, `v2` and `r` define a circle with maximum radius.
 * @param eid index to the segment
 * @param ref_vid index to the reference encroaching point
 * @param ref_tidoff index_off to the tetrahedron containing the reference
 * encroaching point
 * @param enc_verts indices to the encroaching points
 * @see Section 3.3 Segment recovery, in [Robust CDT].
 */
template <typename Traits>
template <typename Container>
void SegmentRecover<Traits>::findReferenceEncroachingPoint(index_t  eid,
                                                           index_t &ref_vid,
                                                           index_t &ref_tidoff,
                                                           Container *enc_verts)
{
  InlinedVector64<index_t> encroach_tets;
  const PLCEdge           &edge = plc.edge(eid);

  // find tetrahedra adjacent to the first endpoint `ep0`
  tet_mesh.VT(edge.ep0(), encroach_tets);
  // and mark the tetrahedra as visited
  for (index_t tet_idoff : encroach_tets)
    tet_mesh.mark(tet_idoff, TET_MARK::VISITED);

  tet_mesh.mark(edge.ep0(), VTX_MARK::VISITED);
  tet_mesh.mark(edge.ep1(), VTX_MARK::VISITED);
  tet_mesh.mark(edge.ep0(), VTX_MARK::ENCROACHED);
  tet_mesh.mark(edge.ep1(), VTX_MARK::ENCROACHED);

  const GPoint3 &p0    = point(edge.ep0());
  const GPoint3 &p1    = point(edge.ep1());
  const GPoint3 *ref_p = nullptr;
  ref_vid              = InvalidIndex;
  ref_tidoff           = InvalidIndex;

  for (index_t i = 0; i < encroach_tets.size(); i++)
  {
    index_t tet_idoff = encroach_tets[i];

    // check if the tetrahedron has an encroaching point
    for (index_t j = 0; j < 4; j++)
    {
      index_t vid = tet_mesh.tetNode(tet_idoff + j);
      if (tet_mesh.isMarked(vid, VTX_MARK::VISITED) ||
          tet_mesh.isMarked(vid, VTX_MARK::ENCROACHED))
        continue;
      tet_mesh.mark(vid, VTX_MARK::VISITED);
      const GPoint3 &curr_p = point(vid);

      // check if the vertex is encroaching
      if (inSphere(p0, p1, curr_p))
      {
        tet_mesh.mark(vid, VTX_MARK::ENCROACHED);
        if (enc_verts) // store the encroaching vertices if required
          enc_verts->push_back(vid);
        // check if it is the reference encroaching point
        if (ref_vid == InvalidIndex || largerSphere(p0, p1, curr_p, *ref_p))
        {
          ref_vid    = vid;
          ref_tidoff = tet_idoff;
          ref_p      = &curr_p;
        }
      }
    }

    // clang-format off
    const int is_encroached[] = {
      tet_mesh.isMarked(tet_mesh.tetNode(tet_idoff), VTX_MARK::ENCROACHED),
      tet_mesh.isMarked(tet_mesh.tetNode(tet_idoff + 1), VTX_MARK::ENCROACHED),
      tet_mesh.isMarked(tet_mesh.tetNode(tet_idoff + 2), VTX_MARK::ENCROACHED),
      tet_mesh.isMarked(tet_mesh.tetNode(tet_idoff + 3), VTX_MARK::ENCROACHED),
    };
    const int total_encroached = is_encroached[0] + is_encroached[1] + is_encroached[2] + is_encroached[3];
    // clang-format on

    // Add neighboring tetrahedra to `encroach_tets` for subsequent encroachment
    // checks. The cases are divided based on the number of encroaching vertices
    // in the current tetrahedron:
    //
    // - If the current tetrahedron has no encroaching vertices:
    //     No neighboring tetrahedron is added.
    // - If the current tetrahedron has only one encroaching vertex:
    //    Neighboring tetrahedra at corners other than the encroaching corner
    //    are added.
    // - If the current tetrahedron has more than one encroaching vertex:
    //    All neighboring tetrahedra are added.

    for (index_t j = 0; j < 4; j++)
    {
      index_t neigh_idoff = tet_mesh.tetNeigh(tet_idoff + j);

      if (tet_mesh.isMarked(neigh_idoff, TET_MARK::VISITED) ||
          !tet_mesh.isFiniteTet(neigh_idoff))
        continue;

      if (total_encroached - is_encroached[j] > 0)
      {
        encroach_tets.push_back(TetMesh::clipId(neigh_idoff));
        tet_mesh.mark(neigh_idoff, TET_MARK::VISITED);
      }
    }
  }

  // clear all marks
  tet_mesh.unmark(edge.ep0(), VTX_MARK::VISITED);
  tet_mesh.unmark(edge.ep1(), VTX_MARK::VISITED);
  tet_mesh.unmark(edge.ep0(), VTX_MARK::ENCROACHED);
  tet_mesh.unmark(edge.ep1(), VTX_MARK::ENCROACHED);
  for (index_t idoff : encroach_tets)
  {
    tet_mesh.unmark(idoff, TET_MARK::VISITED);
    for (index_t j = 0; j < 4; j++)
    {
      index_t vid = tet_mesh.tetNode(idoff + j);
      tet_mesh.unmark(vid, VTX_MARK::VISITED);
      tet_mesh.unmark(vid, VTX_MARK::ENCROACHED);
    }
  }
}

/**
 * @brief Find the splitting point to split the constrained edge `eid` that has
 * no acute vertices, but do not really split the edge.
 * @param eid The constrained edge to split.
 * @param ref_vid The reference encroaching point.
 * @return the splitting point.
 * @see The strategy is described in Section 3.3 of [Robust CDT].
 */
template <typename Traits>
auto SegmentRecover<Traits>::splitSegment_NoAcuteVertex(index_t eid,
                                                        index_t ref_vid) const
  -> IPoint3T_LNC
{
  const PLCEdge &edge = plc.edge(eid);

  const GPoint3 &ep0_pnt = point(edge.ep0());
  const GPoint3 &ep1_pnt = point(edge.ep1());
  const GPoint3 &ref_pnt = point(ref_vid);

  if (isLessThanHalfDistance(ep0_pnt, ref_pnt, ep1_pnt))
  { // The `ref_pnt` is closer to the endpoint `ep0`, and the distance between
    // `ref_pnt` and `ep0` is less than half the distance between `ep0` and
    // `ep1`.
    return lineSphereIntersection_noAc(eid, false, ref_vid);
  }
  else if (isLessThanHalfDistance(ep1_pnt, ref_pnt, ep0_pnt))
  { // The `ref_pnt` is closer to the endpoint `ep1`, and the distance between
    // `ref_pnt` and `ep1` is less than half the distance between `ep0` and
    // `ep1`.
    return lineSphereIntersection_noAc(eid, true, ref_vid);
  }
  else
  { // The distances between both <`ref_pnt`, `ep0`> and <`ref_pnt`, `ep1`>
    // are larger than half the distance between `ep0` and `ep1`.
    // Split the segment at the middle point.
    return middlePoint(edge);
  }
}

/**
 * @brief Find the splitting point to split the constrained edge `eid` that has
 * only one acute vertex, but do not really split the segment.
 * @param eid The constrained edge to split.
 * @param ref_vid the reference encroaching point.
 * @return the new splitting point.
 * @see The strategy is described in Section 3.3 of [Robust CDT].
 */
template <typename Traits>
auto SegmentRecover<Traits>::splitSegment_OneAcuteVertex(index_t eid,
                                                         index_t ref_vid) const
  -> IPoint3T_LNC
{
  const PLCEdge &edge = plc.edge(eid);

  IPoint3T_LNC new_pnt =
    lineSphereIntersection_oneAc(eid, edge.acute_vid, ref_vid);

  if (isLessThanDistance(new_pnt, point(edge.ep1()), point(ref_vid)))
  { // The new point is closer to the non-acute endpoint `ep1` than the
    // reference encroaching point, we should switch to another split strategy
    // (described in Section 3.3 in [Robust CDT]), but it is not really
    // necessary. Just using midpoint provides better performances
    new_pnt = middlePoint(edge);
  }

  return new_pnt;
}

/**
 * @brief Two points `a` and `b` form a smallest (diametral) sphere, check if
 * the query point `c` lies inside or touches the sphere.
 * @return True if the point lies inside or touches the sphere, false otherwise.
 */
template <typename Traits>
bool SegmentRecover<Traits>::inSphere(const GPoint3 &a, const GPoint3 &b,
                                      const GPoint3 &c)
{
  // In [Robust CDT 2023], this predicate is implemented in floating-point
  // numbers, I reimplemnted it in exact predicates.

#ifdef OMC_CDT_EXACT_ENCROACH_TEST
  return InSphere()(a, b, c) >= Sign::ZERO;
#else
  EPoint3 a_ep = ToEP()(a), b_ep = ToEP()(b), c_ep = ToEP()(c);
  return (c_ep - a_ep).sqrnorm() + (c_ep - b_ep).sqrnorm() <=
         (a_ep - b_ep).sqrnorm();
#endif
}

/**
 * @brief Points `a`, `b` and `c` form a smallest sphere S(abc), and points
 * `a`, `b` and `d` form another smallest sphere S(abd). We want to know if
 * S(abc) is larger than S(abd).
 * @return true if S(abc) is larger than S(abd), false otherwise.
 */
template <typename Traits>
bool SegmentRecover<Traits>::largerSphere(const GPoint3 &a, const GPoint3 &b,
                                          const GPoint3 &c, const GPoint3 &d)
{
  // It is complex to generate a complete point arrangement to calculate an
  // exact predicate. This predicate tolerate floating-point errors, so I
  // directly use the inexact predicate.

#if 0
  return InSphere().largerSphere(a, b, c, d) == Sign::POSITIVE;
#else
  // calculate vectors between points
  EPoint3 a_ep = ToEP()(a), b_ep = ToEP()(b), c_ep = ToEP()(c),
          d_ep = ToEP()(d);
  Vec3 ac = (c_ep - a_ep), bc = (c_ep - b_ep);
  Vec3 ad = (d_ep - a_ep), bd = (d_ep - b_ep);

  // let the angle between ac and bc be <abc>, and the angle between ad and bd
  // be <abd>.
  // cos(<abc>) = (ac dot bc) / (|ac| * |bc|), so as cos(<abd>).
  // S(abc) is larger than S(abd) if and only if cos(<abc>) < cos(<abd>).

  NT len_c = ac.sqrnorm() * bc.sqrnorm();
  NT len_d = ad.sqrnorm() * bd.sqrnorm();
  NT dot_c = ac.dot(bc);
  NT dot_d = ad.dot(bd);

  // cos(<abc>) < cos(<abd>) is equivalent to the below inequality.
  return (dot_d * dot_d) * len_c < (dot_c * dot_c) * len_d;
#endif
}

/**
 * @brief Check if the distance between `a` and `b` is less than the distance
 * between `a` and `c`.
 * @return true if less, false otherwise.
 */
template <typename Traits>
bool SegmentRecover<Traits>::isLessThanDistance(const GPoint3 &a,
                                                const GPoint3 &b,
                                                const GPoint3 &c)
{
  // This predicate tolerate floating-point errors, so I directly use the
  // inexact predicate.

#if 0
  return SquaredDistance3D()(a, b, c) == Sign::NEGATIVE;
#else
  EPoint3 a_ep = ToEP()(a), b_ep = ToEP()(b), c_ep = ToEP()(c);
  return (b_ep - a_ep).sqrnorm() < (c_ep - a_ep).sqrnorm();
#endif
}

/**
 * @brief Check if the distance between `a` and `b` is less than half the
 * distance between `a` and `c`.
 * @return true if less, false otherwise.
 */
template <typename Traits>
bool SegmentRecover<Traits>::isLessThanHalfDistance(const GPoint3 &a,
                                                    const GPoint3 &b,
                                                    const GPoint3 &c)
{
  // This predicate tolerate floating-point errors, so I directly use the
  // inexact predicate.

#if 0
  // scale the squared distance between a and b by 4 to avoid square root.
  return SquaredDistance3D()(a, b, c, /*ab_scale*/ 4) == Sign::NEGATIVE;
#else
  EPoint3 a_ep = ToEP()(a), b_ep = ToEP()(b), c_ep = ToEP()(c);
  return (b_ep - a_ep).sqrnorm() * 4 < (c_ep - a_ep).sqrnorm();
#endif
}

/**
 * @brief Simply get the inerpolation parameters for the sub-edge.
 * ep0 = oep0 + t0 * (oep1 - oep0),
 * ep1 = oep0 + t1 * (oep1 - oep0).
 * We assume that `oep0.t <= ep0.t < ep1.t <= oep1.t`.
 *
 * @param oep0 one endpoint of the original PLC edge
 * @param oep1 the other endpoint of the original PLC edge
 * @param ep0  one endpoint of the sub-edge
 * @param ep1 the other endpoint of the sub-edge
 * @return A pair of interpolation parameters t0 and t1.
 */
template <typename Traits>
std::pair<double, double>
SegmentRecover<Traits>::getInterpolateT(index_t oep0, index_t oep1, index_t ep0,
                                        index_t ep1) const
{
  double t0 = -1.0, t1 = -1.0;

  // calculate the interpolation parameter for `ep0`
  if (ep0 == oep0)
    t0 = 0.0;
  else
  {
    OMC_EXPENSIVE_ASSERT(&(point(ep0).LNC().P()) == &point(oep0) &&
                           &(point(ep0).LNC().Q()) == &point(oep1),
                         "Wrong orders.");
    t0 = point(ep0).LNC().T();
  }

  // calculate the interpolation parameter for `ep1`
  if (ep1 == oep1)
    t1 = 1.0;
  else
  {
    OMC_EXPENSIVE_ASSERT(&(point(ep1).LNC().P()) == &point(oep0) &&
                           &(point(ep1).LNC().Q()) == &point(oep1),
                         "Wrong orders.");
    t1 = point(ep1).LNC().T();
  }

  OMC_EXPENSIVE_ASSERT(t0 < t1, "Wrong orders.");
  return std::pair<double, double>(t0, t1);
}

/**
 * @brief Get the middle point of a PLC edge represented by LNC implicit point.
 * @param e The given PLC edge.
 * @return IPoint3T_LNC The middle point of the edge.
 */
template <typename Traits>
auto SegmentRecover<Traits>::middlePoint(const PLCEdge &e) const -> IPoint3T_LNC
{
  index_t ep0 = e.ep0(), ep1 = e.ep1();

  if (e.hasAncestor()) // The edge is a sub-edge of a split edge.
  {
    const PLCEdge &oe   = plc.edge(e.ancestor_id);
    index_t        oep0 = oe.ep0(), oep1 = oe.ep1();

    auto [t0, t1] = getInterpolateT(oep0, oep1, ep0, ep1);
    OMC_EXPENSIVE_ASSERT((t0 + t1) * 0.5 != t0 && (t0 + t1) * 0.5 != t1,
                         "The edge can not be split at midpoint.");
    return CreateLNC()(point(oep0), point(oep1), (t0 + t1) * 0.5);
  }
  else // The edge is not split yet.
  {
    return CreateLNC()(point(ep0), point(ep1), /*interpolation T*/ 0.5);
  }
}

/**
 * @brief Get the intersection point between the line and the sphere,
 * represented as an LNC implicit point.
 *
 * - The line is defined by the constrained edge `eid`.
 *
 * - The sphere is centered at one endpoint (determined by `reverse`) of the
 * edge, with the radius being the distance between the endpoint and the
 * reference encroaching point `ref_vid`.
 *
 * @param [in] eid The index of the constrained edge.
 * @param [in] reverse TRUE if the reference encroaching point is closer to the
 * endpoint `ep1`, FALSE if it is closer to the endpoint `ep0`. The sphere
 * centers at the endpoint that is closer to the reference encroaching point.
 * @param [in] ref_vid The index of the reference encroaching point.
 * @note `noAc` means no acute vertex.
 * @return IPoint3T_LNC The intersection point represented in LNC.
 */
template <typename Traits>
auto SegmentRecover<Traits>::lineSphereIntersection_noAc(index_t eid,
                                                         bool    reverse,
                                                         index_t ref_vid) const
  -> IPoint3T_LNC
{
  const PLCEdge &e = plc.edge(eid);

  // Get the endpoints of the edge and its original edge.
  index_t ep0 = e.ep0(), ep1 = e.ep1();
  index_t oep0 = e.hasAncestor() ? plc.edge(e.ancestor_id).ep0() : ep0;
  index_t oep1 = e.hasAncestor() ? plc.edge(e.ancestor_id).ep1() : ep1;

  OMC_EXPENSIVE_ASSERT(point(oep0).is_explicit() && point(oep1).is_explicit(),
                       "Input points contain implicit points.");
  // Get the vectors of related points.
  Vec3 oe0_v    = AsEP()(point(oep0)).as_vec();
  Vec3 oe1_v    = AsEP()(point(oep1)).as_vec();
  Vec3 ref_v    = ToEP()(point(ref_vid)).as_vec();
  Vec3 end_v    = ToEP()(point(reverse ? ep1 : ep0)).as_vec();
  // Get the interpolation parameters
  auto [t0, t1] = getInterpolateT(oep0, oep1, ep0, ep1);
  // Parameterize the sphere radius to the original segment
  double radius =
    std::sqrt((ref_v - end_v).sqrnorm() / (oe1_v - oe0_v).sqrnorm());
  // Get the parameter of the intersection point
  double t = reverse ? t1 - radius : t0 + radius;
  // Check if the intersection point is inside the edge
  if (t <= t0 || t >= t1)
  { // if no (maybe caused by numerical error), return the middle point
    t = (t0 + t1) * 0.5;
  }
  OMC_EXPENSIVE_ASSERT((t0 < t && t < t1), "The point is outside the edge.");
  return CreateLNC()(point(oep0), point(oep1), t);
}

/**
 * @brief  Get the intersection point between the line and the sphere,
 * represented as an LNC implicit point.
 *
 * - The line is defined by the constrained edge `eid`.
 *
 * - The sphere is centered at the acute endpoint of the edge, with the radius
 * being the distance between the acute endpoint and the reference encroaching
 * point `ref_vid`.
 *
 * - Acute endpoint is always put in the first position.
 *
 * @param [in] eid The index of the constrained edge.
 * @param [in] acute_vid The index of the acute vertex.
 * @param [in] ref_vid The index of the reference encroaching point.
 * @note `oneAc` means one acute vertex.
 * @return IPoint3T_LNC The intersection point represented in LNC.
 */
template <typename Traits>
auto SegmentRecover<Traits>::lineSphereIntersection_oneAc(index_t eid,
                                                          index_t acute_vid,
                                                          index_t ref_vid) const
  -> IPoint3T_LNC
{
  const PLCEdge &e = plc.edge(eid);

  // Get the endpoints of the edge and its original edge.
  index_t ep0 = e.ep0(), ep1 = e.ep1();
  index_t oep0 = e.hasAncestor() ? plc.edge(e.ancestor_id).ep0() : ep0;
  index_t oep1 = e.hasAncestor() ? plc.edge(e.ancestor_id).ep1() : ep1;

  OMC_EXPENSIVE_ASSERT((acute_vid == oep0 || acute_vid == oep1),
                       "Wrong acute vertices.");
  OMC_EXPENSIVE_ASSERT(point(oep0).is_explicit() && point(oep1).is_explicit(),
                       "Input points contain implicit points.");

  // Get the vectors of related points.
  Vec3 oe0_v    = AsEP()(point(oep0)).as_vec();
  Vec3 oe1_v    = AsEP()(point(oep1)).as_vec();
  Vec3 ref_v    = ToEP()(point(ref_vid)).as_vec();
  // Get the interpolation parameters
  auto [t0, t1] = getInterpolateT(oep0, oep1, ep0, ep1);
  // Parameterize the sphere radius to the original segment
  double t =
    acute_vid == oep0
      ? std::sqrt((ref_v - oe0_v).sqrnorm() / (oe1_v - oe0_v).sqrnorm())
      : 1.0 - std::sqrt((ref_v - oe1_v).sqrnorm() / (oe1_v - oe0_v).sqrnorm());
  // Ensure that the intersection point is inside the edge,
  // and make the sub segments as long as possible.
  double eps = (t1 - t0) * 0.2;
  if (t <= (t0 + eps) || t >= (t1 - eps))
  { // Otherwise return the middle point
    t = (t0 + t1) * 0.5;
  }
  OMC_EXPENSIVE_ASSERT((t0 < t && t < t1), "The point is outside the edge.");
  return CreateLNC()(point(oep0), point(oep1), t);
}

/**
 * @brief Constructs protecting spheres for vertices in the PLC.
 *
 * If a vertex needs protection, it constructs an initial protecting sphere
 * based on the shortest edge lengths connected to the vertex.
 * It then checks for intersections with other segments and shrinks the sphere
 * if necessary to ensure the vertex is well-protected.
 * @see Shewchuk, J. R. Constrained Delaunay Tetrahedralizations and Provably
 * Good Boundary Recovery. In Proceedings of the International Meshing
 * Roundtable Conference (2002).
 */
template <typename Traits>
void SegmentRecover<Traits>::buildProtectingSphere()
{
  if (protecting_sphere_initialized)
    return;

  // Build segment tree for input constrained segments =====================

  // segment tree
  std::vector<IndexedSegment> indexed_segments;
  indexed_segments.reserve(plc.numEdges());
  // collect primitives
  for (index_t eid = 0; eid < plc.numEdges(); eid++)
  {
    const PLCEdge &e = plc.edge(eid);
    if (e.isConstraint())
    {
      // endpoints
      index_t        ev0 = e.ep0(), ev1 = e.ep1();
      const GPoint3 &gp0 = point(ev0), &gp1 = point(ev1);
      // segment tree
      indexed_segments.emplace_back(Segment3(AsEP()(gp0), AsEP()(gp1)), eid);
    }
  }
  // segment tree
  seg_tree.insert(indexed_segments.begin(), indexed_segments.end());
  seg_tree.build();

  // Calculate protecting spheres for vertices =============================

  protecting_sphere_squared_radius.resize(plc.numVertices(), -1.0);

  auto approxEdgeLength = [this](const PLCEdge &e)
  { return (point(e.ep0()) - point(e.ep1())).sqrnorm(); };

  auto buildSphere = [this, &approxEdgeLength](index_t vid)
  {
    // check if the vertex need to be protected,
    // get the shortest edge length at the same time.
    bool   need_protect            = false;
    double shortest_acute_edge     = DBL_MAX; // l_i in reference paper
    double shortest_non_acute_edge = DBL_MAX; // d_i in reference paper

    boost::container::flat_set<index_t> inc_edges;
    auto [inc_edges_begin, inc_edges_end] = plc.vertIncEdges(vid);
    for (auto iter = inc_edges_begin; iter != inc_edges_end; iter++)
    {
      index_t        eid = *iter;
      const PLCEdge &e   = plc.edge(eid);
      inc_edges.insert(eid);

      if (e.type == PLCEdgeType::NO_ACUTE_VERTEX)
        continue;
      else if (e.type == PLCEdgeType::ONE_ACUTE_VERTEX)
      {
        need_protect = true;
        if (e.acute_vid == vid)
        {
          shortest_acute_edge =
            std::min(shortest_acute_edge, approxEdgeLength(e));
        }
        else
        {
          shortest_non_acute_edge =
            std::min(shortest_non_acute_edge, approxEdgeLength(e));
        }
      }
      else // e.type == PLCEdgeType::BOTH_ACUTE_VERTEX
      {
        need_protect = true;
        shortest_acute_edge =
          std::min(shortest_acute_edge, approxEdgeLength(e));
      }
    }

    if (!need_protect)
      return;

    // construct a initial protecting sphere based on the above lengths.
    Sphere3 sphere(AsEP()(point(vid)),
                   std::min(shortest_acute_edge * (1.0 / 4.0),
                            shortest_non_acute_edge * (4.0 / 9.0)));

    // find segments that intersect the currect sphere
    InlinedVector64<index_t> intersected_edges;
    {
      // find possible intersected segments by AABB tree
      InlinedVector64<index_t> possible_intersected_edges;
      seg_tree.all_intersections(sphere, possible_intersected_edges);
      // find really intersected segments
      boost::container::flat_set<index_t> _intersected_edges;
      for (index_t eid : possible_intersected_edges)
      {
        const PLCEdge &e = plc.edge(eid);
        Segment3       seg(AsEP()(point(e.ep0())), AsEP()(point(e.ep1())));
        if (DoIntersect()(sphere, seg))
          _intersected_edges.insert(eid);
      }
      // remove incident edges in found intersected edges.
      std::set_difference(_intersected_edges.begin(), _intersected_edges.end(),
                          inc_edges.begin(), inc_edges.end(),
                          std::back_inserter(intersected_edges));
    }

    // if no intersected segments, the vertex is well-protected.
    if (intersected_edges.empty())
    {
      protecting_sphere_squared_radius[vid] = sphere.squared_radius();
      return;
    }

    // otherwise, we need to shrink the sphere.
    InlinedVector64<Segment3> intersected_segs;
    for (index_t eid : intersected_edges)
    {
      const PLCEdge &e = plc.edge(eid);
      intersected_segs.emplace_back(AsEP()(point(e.ep0())),
                                    AsEP()(point(e.ep1())));
    }

    // the shortest projection distance from the vertex to the intersected
    // segments is also called `local feature size` (lfs).
    // we need to shrink the sphere to be less than the local feature size.
    double min_proj_dis = DBL_MAX;
    for (const Segment3 &seg : intersected_segs)
    {
      EPoint3 proj_pnt = ProjectPoint3()(seg, sphere.center());
      double  proj_dis = (proj_pnt - sphere.center()).sqrnorm();
      min_proj_dis     = std::min(min_proj_dis, proj_dis);
    }

    // reduce the projection distance to avoid numerical error
    sphere.squared_radius() =
      std::min(min_proj_dis * 0.5, sphere.squared_radius());
    // check if the sphere is still intersected by segments
    while (!intersected_segs.empty())
    {
      if (DoIntersect()(sphere, intersected_segs.back()))
        sphere.squared_radius() *= 0.5;
      else
        intersected_segs.pop_back();
    }

    protecting_sphere_squared_radius[vid] = sphere.squared_radius();
  };

  tbb::parallel_for(index_t(0), plc.numVertices(), buildSphere);
  protecting_sphere_initialized = true;
}

template <typename Traits>
void SegmentRecover<Traits>::protectVertex(index_t       eid,
                                           IPoint3T_LNC &steiner_point) const
{
  // use a protecting sphere to protect edges around the acute vertex.
  const PLCEdge &anc_edge = plc.edge(plc.ancestorEdge(eid));
  index_t        oep0 = anc_edge.ep0(), oep1 = anc_edge.ep1();

  if (protecting_sphere_squared_radius[oep0] > 0.0 &&
      DoIntersect()(
        Sphere3(AsEP()(point(oep0)), protecting_sphere_squared_radius[oep0]),
        AsGP()(steiner_point)))
  {
    steiner_point = splitSegment_ProtectingSphere(eid, oep0);
  }
  else if (protecting_sphere_squared_radius[oep1] > 0.0 &&
           DoIntersect()(Sphere3(AsEP()(point(oep1)),
                                 protecting_sphere_squared_radius[oep1]),
                         AsGP()(steiner_point)))
  {
    steiner_point = splitSegment_ProtectingSphere(eid, oep1);
  }
}

/**
 * @brief Get the Steiner point that is the intersection of the segment and
 * its protecting sphere.
 * @param eid The edge to split.
 * @param center_vid The center vertex of the protecting sphere.
 * @return The Steiner point.
 */
template <typename Traits>
auto SegmentRecover<Traits>::splitSegment_ProtectingSphere(
  index_t eid, index_t center_vid) const -> IPoint3T_LNC
{
  const PLCEdge &e     = plc.edge(eid);
  const PLCEdge &anc_e = plc.edge(plc.ancestorEdge(eid));

  // Get the endpoints of the edge and its original edge.
  index_t ep0 = e.ep0(), ep1 = e.ep1();
  index_t oep0 = anc_e.ep0(), oep1 = anc_e.ep1();
  // Get the interpolation parameters
  auto [t0, t1] = getInterpolateT(oep0, oep1, ep0, ep1);

  OMC_EXPENSIVE_ASSERT(
    center_vid == oep0 || center_vid == oep1,
    "The center vertex is not an endpoint of the original edge.");
  OMC_EXPENSIVE_ASSERT(point(oep0).is_explicit() && point(oep1).is_explicit(),
                       "Input points contain implicit points.");

  // Get the vectors of related points.
  Vec3 oe0_v = AsEP()(point(oep0)).as_vec();
  Vec3 oe1_v = AsEP()(point(oep1)).as_vec();

  // Parameterize the sphere radius to the original segment
  double t = std::sqrt(protecting_sphere_squared_radius[center_vid] /
                       (oe1_v - oe0_v).sqrnorm());
  t        = center_vid == oep0 ? t : 1.0 - t;
  // WARN t still has numerical error
  OMC_EXPENSIVE_ASSERT((t0 < t && t < t1), "The point is outside the edge.");
  return CreateLNC()(point(oep0), point(oep1), t);
}

} // namespace OMC
