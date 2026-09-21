#pragma once

#include "DelaunayRefine.h"

#ifdef OMC_DELAUNAY_REFINE_PROFILE
  #include "OpenMeshCraft/Utils/FormatCompat.h"
  #include "OpenMeshCraft/Utils/Logger.h"

  #include <iostream>
#endif

namespace OMC {

template <typename T, typename D>
DelaunayRefine<T, D>::DelaunayRefine(const Domain   &_domain,
                                     const Criteria &_criteria, Points &_points,
                                     Weights &_weights, TetMesh &_tet_mesh)
  : domain(_domain)
  , criteria(_criteria)
  , points(_points)
  , weights(_weights)
  , tet_mesh(_tet_mesh)
{
}

template <typename T, typename D>
DelaunayRefine<T, D>::ProtectBall::ProtectBall(
  const Point3 &_center, NT _radius_sq,
  FeatureVertexIndex _feature_vertex_index, std::true_type /*is_vertex*/)
{
  center            = _center;
  radius_sq         = _radius_sq;
  feature_index     = _feature_vertex_index;
  is_feature_vertex = true;
}

template <typename T, typename D>
DelaunayRefine<T, D>::ProtectBall::ProtectBall(
  const Point3 &_center, NT _radius_sq, FeatureEdgeIndex _feature_edge_index,
  std::false_type /*is_not_vertex*/)
{
  center            = _center;
  radius_sq         = _radius_sq;
  feature_index     = _feature_edge_index;
  is_feature_vertex = false;
}

template <typename T, typename D>
void DelaunayRefine<T, D>::preserveFeatures()
{
  if (!domain.hasFeatures())
    return; // no features to preserve

  // construct protecting balls on feature vertices
  constructBallsOnVertices();

  // construct protecting balls on feature edges
  constructBallsOnEdges();

  // Finally, insert all protect balls as weighted vertex into tetrahedral mesh.
  insertAllBalls();
}

template <typename T, typename D>
void DelaunayRefine<T, D>::constructBallsOnVertices()
{
  std::vector<std::pair<FeatureVertexIndex, Point3>> feature_vertices;
  domain.getFeatureVertices(std::back_inserter(feature_vertices));
  feature_vertex_balls.resize(domain.numberFeatureVertices(), InvalidIndex);

  if (feature_vertices.size() == 0)
    return; // no feature vertices

  if (feature_vertices.size() == 1)
  {
    // only one feature vertex, no need to check the distance between them
    const auto &[idx, pnt] = feature_vertices[0];
    index_t ball_idx       = protect_balls.size();
    protect_balls.emplace_back(
      /*center=*/pnt, /*weight=*/0.0, /*index=*/idx,
      /*is_vertex=*/std::true_type());
    feature_vertex_balls[idx] = ball_idx;
    return;
  }

  for (const auto &[idx, pnt] : feature_vertices)
  {
    // The radius must be less than the distance to the nearest feature
    // vertex/edge except the current feature vertex and its adjacent feature
    // edge.
    Point3   nearest_feature_pnt = domain.nearestFeature(pnt, idx, /*dim=*/0);
    const NT nearest_sq_dist     = (pnt - nearest_feature_pnt).sqrnorm();
    NT       radius_sq           = nearest_sq_dist / 9.0;

    // Ensure the radius is less than the size suggested by sizing field.
    if (domain.hasSizingField())
    {
      // need an exact feature size? pass feature_vertex_index to domain?
      const NT size    = domain.pointSize(pnt);
      const NT size_sq = size * size;
      radius_sq        = std::min(radius_sq, size_sq);
    }

    // avoid too small balls
    if (criteria.hasProtectBallMinimalWeight())
    {
      radius_sq = std::max(radius_sq, criteria.protectBallMinimalWeight());
    }

    index_t ball_idx = protect_balls.size();
    protect_balls.emplace_back(/*center=*/pnt, /*weight=*/radius_sq,
                               /*index=*/idx, /*is_vertex=*/std::true_type());
    feature_vertex_balls[idx] = ball_idx;
  }
}

template <typename T, typename D>
void DelaunayRefine<T, D>::constructBallsOnEdges()
{
  std::vector<std::pair<FeatureEdgeIndex, FeatureEdge>> feature_edges;
  domain.getFeatureEdges(std::back_inserter(feature_edges));
  feature_edge_balls.resize(domain.numberFeatureEdges());

  if (feature_edges.size() == 0)
    return; // no feature edges

  for (const auto &[edge_idx, edge] : feature_edges)
  {
    const auto &[v0_idx, v1_idx]     = edge;
    std::vector<index_t> &edge_balls = feature_edge_balls[edge_idx];

    edge_balls = std::vector<index_t>(
      {feature_vertex_balls[v0_idx], feature_vertex_balls[v1_idx]});

    if (edgeDenseSampled(edge_idx, /*first=*/0, /*last=*/1))
      continue; // edge is already dense sampled, no ball is inserted.

    // edge is not dense sampled, populate more balls
    populateEdge(edge_idx, /*first=*/0, /*last=*/1);

    // Shrink balls to be far from nearest feature vertex/edge to satisfy
    // [Disjoint Non-Overlapping] requirements.
    // We ensure each ball satisfies condition:
    //    `Radius < distance to nearest feature vertex/edge / 3.0`
    size_t shrink_iterations = 0;
    do
    {
      // shrinking may violate [Dense Sampling] requirements.
      // find violating ranges
      InlinedVector64<std::pair<index_t, index_t>> sparse_ranges =
        findSparseSampledRange(edge_idx);
      // repopulate balls on violating ranges
      for (const auto &[first, last] : sparse_ranges)
      {
        // clang-format off
        index_t first_local = static_cast<index_t>(std::distance(edge_balls.begin(), std::find(edge_balls.begin(), edge_balls.end(), first)));
        index_t last_local  = static_cast<index_t>(std::distance(edge_balls.begin(), std::find(edge_balls.begin(), edge_balls.end(), last)));
        // clang-format on
        populateEdge(edge_idx, first_local, last_local);
      }
    } while ((shrink_iterations++) < shrink_max_iterations &&
             shrinkBalls(edge_idx));
  }
}

template <typename T, typename D>
bool DelaunayRefine<T, D>::shrinkBalls(FeatureEdgeIndex feature_edge_idx)
{
  // Shrink balls to be far from nearest feature vertex/edge.
  // Radius < distance to nearest feature vertex/edge / 3.0
  const std::vector<index_t> &edge_balls = feature_edge_balls[feature_edge_idx];
  bool                        shrunk     = false;

  for (size_t i = 1; i < edge_balls.size() - 1; i++)
  { // skip two endpoints
    index_t ball_idx = edge_balls[i];

    ProtectBall &ball = protect_balls[ball_idx];

    // Find the nearest feature vertex/edge to the ball center.
    Point3 nearest_feature_pnt =
      domain.nearestFeature(ball.center, ball.feature_index, /*dim=*/1);
    const NT nearest_sq_dist  = (ball.center - nearest_feature_pnt).sqrnorm();
    NT target_radius_sq = nearest_sq_dist / 9.0;

    // TODO adjust size by sizing field?

    if (criteria.hasProtectBallMinimalWeight())
    { // avoid too small balls
      target_radius_sq = std::max(target_radius_sq,
                                   criteria.protectBallMinimalWeight());
    }

    if (ball.radius_sq <= target_radius_sq)
      continue; // already satisfied

    // Shrink the ball
    ball.radius_sq = target_radius_sq;
    shrunk         = true;
  }

  return shrunk;
}

template <typename T, typename D>
bool DelaunayRefine<T, D>::edgeDenseSampled(FeatureEdgeIndex feature_edge_idx,
                                            index_t first_ball_local_idx,
                                            index_t last_ball_local_idx) const
{
  OMC_EXPENSIVE_ASSERT(first_ball_local_idx + 1 == last_ball_local_idx,
                       "Non-adjacent balls: {} - {}", first_ball_local_idx,
                       last_ball_local_idx);

  index_t first_ball_idx =
    feature_edge_balls.at(feature_edge_idx)[first_ball_local_idx];
  index_t last_ball_idx =
    feature_edge_balls.at(feature_edge_idx)[last_ball_local_idx];
  const ProtectBall &first_ball = protect_balls[first_ball_idx];
  const ProtectBall &last_ball  = protect_balls[last_ball_idx];

  NT distance     = (first_ball.center - last_ball.center).norm();
  NT first_radius = std::sqrt(first_ball.radius_sq);
  NT last_radius  = std::sqrt(last_ball.radius_sq);
  // swap to ensure first_radius <= last_radius
  if (first_radius > last_radius)
    std::swap(first_radius, last_radius);

  // ensure balls intersect significantly
  return distance < dense_sampled_factor * first_radius + last_radius;
}

template <typename T, typename D>
bool DelaunayRefine<T, D>::verifyDenseSampled(
  FeatureEdgeIndex feature_edge_idx) const
{
  for (index_t i = 0; i < feature_edge_balls[feature_edge_idx].size() - 1; ++i)
  {
    if (!edgeDenseSampled(feature_edge_idx, /*first=*/i, /*last=*/i + 1))
      return false;
  }
  return true;
}

template <typename T, typename D>
void DelaunayRefine<T, D>::populateEdge(FeatureEdgeIndex feature_edge_idx,
                                        index_t          first_ball_local_idx,
                                        index_t          last_ball_local_idx)
{
  std::vector<index_t> &edge_balls = feature_edge_balls.at(feature_edge_idx);

  OMC_EXPENSIVE_ASSERT(first_ball_local_idx < last_ball_local_idx,
                       "The first ball {} must be before the last ball {}.",
                       first_ball_local_idx, last_ball_local_idx);

  // if there are balls between the first and last balls, these balls are not
  // dense sampled, we need to delete them and populate new balls.
  if (last_ball_local_idx - first_ball_local_idx > 1)
  {
    balls_deleted.insert(balls_deleted.end(),
                         edge_balls.begin() + first_ball_local_idx + 1,
                         edge_balls.begin() + last_ball_local_idx);
    edge_balls.erase(edge_balls.begin() + first_ball_local_idx + 1,
                     edge_balls.begin() + last_ball_local_idx);
    last_ball_local_idx = first_ball_local_idx + 1;
  }

  // populate new balls between the first and last balls
  index_t first_ball_idx = edge_balls[first_ball_local_idx];
  index_t last_ball_idx  = edge_balls[last_ball_local_idx];

  const ProtectBall &first_ball = protect_balls[first_ball_idx];
  const ProtectBall &last_ball  = protect_balls[last_ball_idx];

  const NT first_radius = std::sqrt(first_ball.radius_sq);
  const NT last_radius  = std::sqrt(last_ball.radius_sq);

  const bool adjust_orientation = first_radius > last_radius;

  // Populate method from CGAL::Mesh_3::Protect_edges_sizing_field.
  //
  // Notations:
  // sp = first_radius , sq = last_radius,  d = first_last_geodesic_distance
  // n = nb_points,   r = delta_step_size
  //
  // Hypothesis:
  // sp <= sq
  //
  // Let's define
  // P0 = p, Pn+1 = q, d(Pi,Pi+1) = ai
  //
  // The following constraints should be verified:
  // a0 = sp + r, an = sq,
  // ai+1 = ai + r
  // d = Sum(ai)
  //
  // The following could be obtained:
  // r = (sq - sp) / (n+1)
  // n = 2(d-sq) / (sp+sq)
  //
  // =======================
  // Calculus details:
  // ai+1 = ai + r
  // ai+1 = a0 + r*(i+1)
  //   an = a0 + r*n
  //   sq = sp + r + r*n
  //    r = (sq-sp) / (n+1)
  //
  //   d = Sum(ai)
  //   d = Sum(sp + (i+1)*r)
  //   d = (n+1)*sp + (n+1)(n+2)/2 * r
  //   d = (n+1)*sp + (n+1)(n+2)/2 * (sq-sp) / (n+1)
  // 2*d = 2(n+1)*sp + (n+2)*sq - (n+2)*sp
  // 2*d = n*sp + (n+2)*sq
  //   n = 2(d-sq) / (sp+sq)
  // =======================

  auto populate_balls =
    [this](const Point3 &p, NT sp, const Point3 &q, const NT sq, const NT d)
  {
    int  n = static_cast<int>(std::ceil(NT(2.0) * (d - sq) / (sp + sq)));
    NT   r = (sq - sp) / NT(n + 1);
    Vec3 v = (q - p).normalized();

    // Since n is underestimated, we need to adjust size of steps
    // CD = covered distance
    NT CD             = sp * NT(n + 1) + NT((n + 1) * (n + 2)) / NT(2) * r;
    NT dleft_frac     = d / CD;
    // Initialize step sizes
    NT step_size      = sp + r;
    NT norm_step_size = dleft_frac * step_size;
    // Initial distance
    NT pt_dist        = norm_step_size;

    // If there is some place to insert one point, insert it
    if ((0 == n) && (d >= sp + sq))
    {
      n              = 1;
      step_size      = sp + (d - sp - sq) / NT(2);
      pt_dist        = step_size;
      norm_step_size = step_size;
    }

    std::vector<std::pair<Point3, NT>> balls;
    // Launch balls
    for (int i = 1; i <= n; ++i)
    {
      // New point position
      Point3 new_point = p + v * pt_dist;

      // New point weight
      // (min between norm_step_size and linear interpolation,
      //  linear interpolation size is used in repopulation process)
      NT current_size = std::min(norm_step_size, sp + pt_dist / d * (sq - sp));
      NT point_weight = current_size * current_size;

      // TODO adjust by sizing field?

      // Avoid too small balls
      if (criteria.hasProtectBallMinimalWeight())
      {
        point_weight =
          std::max(point_weight, criteria.protectBallMinimalWeight());
      }

      balls.push_back({new_point, point_weight});

      // Step size
      step_size += r;
      norm_step_size = dleft_frac * step_size;

      // Increment distance
      pt_dist += norm_step_size;
    }

    return balls;
  };

  NT edge_length = (last_ball.center - first_ball.center).norm();
  std::vector<std::pair<Point3, NT>> balls =
    adjust_orientation
      ? populate_balls(last_ball.center, last_radius, first_ball.center,
                       first_radius, edge_length)
      : populate_balls(first_ball.center, first_radius, last_ball.center,
                       last_radius, edge_length);

  auto new_ball = [this](const Point3 &pnt, NT radius_sq,
                         index_t feature_edge_idx) -> index_t
  {
    if (balls_deleted.empty())
    {
      protect_balls.emplace_back(pnt, radius_sq, feature_edge_idx,
                                 /*is_vertex=*/std::false_type());
      return protect_balls.size() - 1;
    }
    else
    {
      index_t ball_idx = balls_deleted.back();
      balls_deleted.pop_back();
      protect_balls[ball_idx] = ProtectBall(pnt, radius_sq, feature_edge_idx,
                                            /*is_vertex=*/std::false_type());
      return ball_idx;
    }
  };

  // Construct new balls
  InlinedVector64<index_t> new_balls;
  for (const auto &[pnt, radius_sq] : balls)
    new_balls.push_back(new_ball(pnt, radius_sq, feature_edge_idx));

  // Insert balls into the feature edge
  if (adjust_orientation)
    edge_balls.insert(edge_balls.begin() + last_ball_local_idx,
                      new_balls.rbegin(), new_balls.rend());
  else
    edge_balls.insert(edge_balls.begin() + last_ball_local_idx,
                      new_balls.begin(), new_balls.end());
  last_ball_local_idx = first_ball_local_idx + new_balls.size() + 1;
}

template <typename T, typename D>
InlinedVector64<std::pair<index_t, index_t>>
DelaunayRefine<T, D>::findSparseSampledRange(
  FeatureEdgeIndex feature_edge_idx) const
{
  // Identify ranges of balls on the same edge that violate the [Dense Sampling]
  // condition. Then split these ranges into subranges where the ball radius
  // decreases monotonically after repopulation.
  //
  // During repopulation, each new ball's radius is smaller than the linearly
  // interpolated radius at its center. Therefore, if a ball's current radius is
  // larger than the interpolated radius, the new radius will be smaller as
  // well. In this case, it is safe to repopulate balls between these points.

  const std::vector<index_t> &edge_balls =
    feature_edge_balls.at(feature_edge_idx);
  const index_t begin = 0, end = edge_balls.size() - 1;
  const Point3 &begin_center = protect_balls[edge_balls[begin]].center;

  InlinedVector64<NT> radius(edge_balls.size(), 0.0);
  InlinedVector64<NT> length(edge_balls.size(), 0.0);
  for (size_t i = 0; i < edge_balls.size(); i++)
  {
    index_t            ball_idx = edge_balls[i];
    const ProtectBall &ball     = protect_balls[ball_idx];

    radius[i] = std::sqrt(ball.radius_sq);
    length[i] = (ball.center - begin_center).norm();
  }

  // Verify if the linear interpolation value, determined by the points (rp, lp)
  // and (rn, ln), is smaller than the actual value `rn` at the position `ln`.
  // `p` = previous, `c` = current, `n` = next
  auto isLinearInterpolationDecrease =
    [](NT rp, NT lp, NT rc, NT lc, NT rn, NT ln)
  { return (rn - rp) * (lc - lp) <= (rc - rp) * (ln - lp); };

  // Iterate through the balls on the edge and find ranges of balls that are
  // not dense sampled and decrease monotonically after interpolation.
  InlinedVector64<std::pair<index_t, index_t>> sparse_ranges;
  index_t first_ball_idx = begin, last_ball_idx = begin;
  while (last_ball_idx < end)
  {
    if (edgeDenseSampled(feature_edge_idx, last_ball_idx, last_ball_idx + 1))
    { // [last, last + 1] is dense sampled
      if (last_ball_idx > first_ball_idx)
      { // if a sparse range is found, add it to the result
        sparse_ranges.push_back(
          {edge_balls[first_ball_idx], edge_balls[last_ball_idx]});
      }
      // go to check subsequent range [last + 1, ...)
      first_ball_idx = last_ball_idx + 1;
      last_ball_idx  = first_ball_idx;
      continue;
    }

    // now [last, last + 1] is not dense sampled

    if (last_ball_idx > first_ball_idx &&
        !isLinearInterpolationDecrease(
          /*rp=*/radius[first_ball_idx], /*lp=*/length[first_ball_idx],
          /*rc=*/radius[last_ball_idx], /*lc=*/length[last_ball_idx],
          /*rn=*/radius[last_ball_idx + 1], /*ln=*/length[last_ball_idx + 1]))
    { // at least three balls to check (first, last and last + 1),
      // and linear interpolation does not decrease if adding `last + 1`

      // stop checking the current range, add it to the result
      sparse_ranges.push_back(
        {edge_balls[first_ball_idx], edge_balls[last_ball_idx]});
      // go to check subsequent range [last + 1, ...)
      first_ball_idx = last_ball_idx;
      last_ball_idx  = first_ball_idx;
      continue;
    }

    // now [last, last + 1] is not dense sampled,
    // and linear interpolation on [first, last + 1] is decreasing.
    last_ball_idx++; // check next ball
  }
  if (last_ball_idx > first_ball_idx) // add the last sparse range if it exists
    sparse_ranges.push_back(
      {edge_balls[first_ball_idx], edge_balls[last_ball_idx]});

  return sparse_ranges;
}

template <typename T, typename D>
void DelaunayRefine<T, D>::insertAllBalls()
{
  // clang-format off
#ifdef OMC_DELAUNAY_REFINE_PROFILE
  size_t num_inserted_point    = 0;
  size_t num_coincident_vertex = 0;
  size_t num_hidden_point      = 0;
  #define COLLECT_BALL_INSERT_RESULT(cs)                                               \
    if (cs == BallConflictStatus::OK) { num_inserted_point++; }                        \
    else if (cs == BallConflictStatus::HIDDEN_POINT) { num_hidden_point++; }           \
    else if (cs == BallConflictStatus::COINCIDENT_VERTEX) { num_coincident_vertex++; }
#else
  #define COLLECT_BALL_INSERT_RESULT(cs)
#endif
  // clang-format on

  // Find a tetrahedron to start walking from
  index_t tet = 0;
  while (tet < tet_mesh.sizeTets() * 4 && tet_mesh.isTetDeleted(tet))
    tet += 4;
  OMC_ASSERT(tet < tet_mesh.sizeTets() * 4, "All tets are deleted.");

  // insert balls of all feature vertices
  for (index_t ball_idx : feature_vertex_balls)
  {
    const ProtectBall &ball = protect_balls[ball_idx];
    BallConflictStatus cs   = insertBall(ball, tet);
    COLLECT_BALL_INSERT_RESULT(cs);
  }

  // insert balls on feature edges
  for (const std::vector<index_t> &edge_balls : feature_edge_balls)
  {
    for (size_t i = 1; i < edge_balls.size() - 1; i++)
    { // skip endpoints
      index_t ball_idx = edge_balls[i];

      const ProtectBall &ball = protect_balls[ball_idx];
      BallConflictStatus cs   = insertBall(ball, tet);
      COLLECT_BALL_INSERT_RESULT(cs);
    }
  }

#ifdef OMC_DELAUNAY_REFINE_PROFILE
  Logger::info("\nDelaunayRefine FeaturePreserving profile:");
  Logger::info(OMC::format("  Inserted points: {}", num_inserted_point));
  Logger::info(OMC::format("  Hidden points: {}", num_hidden_point));
  Logger::info(OMC::format("  Coincident vertices: {}", num_coincident_vertex));
#endif
}

template <typename T, typename D>
auto DelaunayRefine<T, D>::insertBall(const ProtectBall &ball, index_t &tet)
  -> BallConflictStatus
{
  DelTet dt(tet_mesh);

  // Locate the point in Delaunay tetrahedralization =======================

  int dimension = -1;
  dt.walk(/*point=*/ball.center, tet, &dimension);

  if (dimension == 0)
  {
    // The point is coincident with a vertex of the tetrahedral mesh.
    return BallConflictStatus::COINCIDENT_VERTEX;
  }

  // Find the conflict zone of the point ==================================

  InlinedVector64<index_t> conflict_tets;
  InlinedVector64<index_t> conflict_corners;
  dt.conflict(/*point=*/ball.center, /*weight=*/ball.radius_sq, tet,
              conflict_tets, conflict_corners);

  if (conflict_tets.empty())
  {
    // The weighted point is hidden by an existing weighted vertex.
    return BallConflictStatus::HIDDEN_POINT;
  }

  // Insert the point into the mesh ========================================
  dt.removeConflicts(conflict_tets);
  index_t new_vtx_vid = newVertex(ball.center, /*weight=*/ball.radius_sq);
  dt.filling(new_vtx_vid, conflict_corners);

  tet = TetMesh::toIdOff(tet_mesh.incTet(new_vtx_vid));
  return BallConflictStatus::OK;
}

template <typename T, typename D>
void DelaunayRefine<T, D>::scanFaces()
{
  size_t num_faces = tet_mesh.sizeFaces();

  for (index_t face_idoff = 0; face_idoff < num_faces; face_idoff++)
  {
    checkNewFace(face_idoff);
  }
}

template <typename T, typename D>
bool DelaunayRefine<T, D>::isFaceRefineDone()
{
  // TODO: An external force stop control
  // TODO: Other internal stop conditions
  return face_queue.empty();
}

template <typename T, typename D>
auto DelaunayRefine<T, D>::processOneFace() -> FaceConflictStatus
{
  // Get the next face to process from the queue ========================

  auto [face_idoff, face_to_refine] = face_queue.top();
  face_queue.pop();

  OMC_EXPENSIVE_ASSERT(!tet_mesh.isTetDeleted(face_idoff),
                       "The face to process is deleted.");

  index_t mirror_idoff = tet_mesh.mirrorFace(face_idoff);
  index_t finite_idoff =
    tet_mesh.isFiniteTet(face_idoff) ? face_idoff : mirror_idoff;

  // Get the refinement point ===========================================

  Point3 refine_point = face_to_refine.intersection;

  // Locate the refinement point in Delaunay tetrahedralization =========

  DelTet del_tet(tet_mesh);

  index_t tet_id    = TetMesh::clipId(finite_idoff);
  int     dimension = -1;
  del_tet.walk(refine_point, tet_id, &dimension);

  if (dimension == 0)
  {
    // The dual point is coincident with a vertex of the tetrahedral mesh.
    return FaceConflictStatus::COINCIDENT_VERTEX;
  }

  // Find the conflict zone of the refinement point =====================

  InlinedVector64<index_t> conflict_tets;
  InlinedVector64<index_t> conflict_corners;
  del_tet.conflict(refine_point, /*weight*/ 0, tet_id, conflict_tets,
                   conflict_corners);

  if (conflict_tets.empty())
  {
    // The dual point is hidden by an existing weighted vertex.
    return FaceConflictStatus::HIDDEN_POINT;
  }

  // check if the face to refine is conflict with the dual point
  if (std::find(conflict_tets.begin(), conflict_tets.end(),
                TetMesh::clipId(face_idoff)) == conflict_tets.end() ||
      std::find(conflict_tets.begin(), conflict_tets.end(),
                TetMesh::clipId(mirror_idoff)) == conflict_tets.end())
  {
    // clear marks
    for (index_t tet_idoff : conflict_tets)
      tet_mesh.unmark(tet_idoff, TET_MARK::CONFLICT);
    for (index_t corner_idoff : conflict_corners)
      tet_mesh.unmark(corner_idoff, TET_MARK::VISITED);
    return FaceConflictStatus::FACE_NOT_CONFLICT;
  }

  // Remove faces in conflict zone from the queue ========================
  for (index_t cf_tet : conflict_tets)
  {
    for (index_t i = 0; i < 4; i++)
    {
      index_t cf_face        = cf_tet + i;
      index_t mirror_cf_face = tet_mesh.mirrorFace(cf_face);
      if (tet_mesh.isMarked(cf_face, FACE_MARK::RESTRICTED))
      {
        if (face_queue.exist(cf_face))
          face_queue.remove(cf_face);
        if (face_queue.exist(mirror_cf_face))
          face_queue.remove(mirror_cf_face);
      }
    }
  }

  // Insert the refinement point into the mesh ============================
  del_tet.removeConflicts(conflict_tets);
  index_t new_vtx_vid = newVertex(refine_point, /*weight*/ 0);
  del_tet.filling(new_vtx_vid, conflict_corners);

  // Check new faces in the conflict zone =================================

  for (index_t corner : conflict_corners) // for each corner `c`
  {
    tet_mesh.unmark(corner, FACE_MARK::VISITED);

    index_t oppo  = tet_mesh.tetNeigh(corner);
    index_t idoff = TetMesh::clipId(oppo);
    for (index_t i = 0; i < 4; i++)
      checkNewFace(idoff + i);
  }

  return FaceConflictStatus::OK;
}

template <typename T, typename D>
void DelaunayRefine<T, D>::refineFaces()
{
#ifdef OMC_DELAUNAY_REFINE_PROFILE
  size_t num_inserted_point    = 0;
  size_t num_coincident_vertex = 0;
  size_t num_hidden_point      = 0;
  size_t num_face_not_conflict = 0;
#endif

  scanFaces();

  while (!isFaceRefineDone())
  {
    FaceConflictStatus cs = processOneFace();

#ifdef OMC_DELAUNAY_REFINE_PROFILE
    switch (cs)
    {
    case FaceConflictStatus::COINCIDENT_VERTEX:
      num_coincident_vertex += 1;
      break;
    case FaceConflictStatus::HIDDEN_POINT:
      num_hidden_point += 1;
      break;
    case FaceConflictStatus::FACE_NOT_CONFLICT:
      num_face_not_conflict += 1;
      break;
    case FaceConflictStatus::OK:
      num_inserted_point += 1;
      if (num_inserted_point % 100 == 0)
        std::cout << OMC::format("\rInserted points {}         ",
                                 num_inserted_point);
      break;
    }
#endif
  }

#ifdef OMC_DELAUNAY_REFINE_PROFILE
  Logger::info("\nDelaunayRefine FaceRefine profile:");
  Logger::info(OMC::format("  Inserted points: {}", num_inserted_point));
  Logger::info(OMC::format("  Coincident vertex: {}", num_coincident_vertex));
  Logger::info(OMC::format("  Hidden point: {}", num_hidden_point));
  Logger::info(OMC::format("  Face not conflict: {}", num_face_not_conflict));
#endif

  tet_mesh.removeDeletedTets();
}

/**
 * @brief Check a new face to determine if it is restricted to a surface patch
 * and evaluate its quality. If the face is found to be restricted, it is added
 * to the refinement queue for further processing.
 *
 * @param face_idoff The idoff of the face in the tetrahedral mesh.
 *
 * @note While checking a new face, no modifications should be made to its
 * associated conflict zone to ensure consistency and correctness.
 */
template <typename T, typename D>
void DelaunayRefine<T, D>::checkNewFace(index_t face_idoff)
{
  index_t mirror_idoff = tet_mesh.mirrorFace(face_idoff);

  if (!tet_mesh.isFiniteTet(face_idoff) && !tet_mesh.isFiniteTet(mirror_idoff))
  {
    // face and its mirror face are both infinite.
    return;
  }

  index_t min_idoff = std::min(face_idoff, mirror_idoff);
  index_t max_idoff = std::max(face_idoff, mirror_idoff);

  if (!(tet_mesh.mark(min_idoff, FACE_MARK::VISITED) &&
        tet_mesh.mark(max_idoff, FACE_MARK::VISITED)))
  { // face or its mirror face have been visited.
    // this scenario can occur during parallel processing when multiple threads
    // attempt to process the same face or its mirror face simultaneously.
    // so, we always try to mark the smaller idoff first.
    return;
  }

  SurfacePatchIndex surface_patch;
  Point3            intersection;
  if (isFaceRestricted(face_idoff, surface_patch, intersection))
  {
    // face is restricted to a surface patch
    tet_mesh.mark(face_idoff, FACE_MARK::RESTRICTED);
    tet_mesh.mark(mirror_idoff, FACE_MARK::RESTRICTED);
    // evaluate the quality
    FaceQuality face_quality =
      criteria.faceQuality(tet_mesh, face_idoff, intersection);

    if (face_quality)
    { // face should be refined.
      // we put the min_idoff to the queue to ensure identifier unique.
      FaceToRefine face_to_refine(face_quality, intersection);
      face_queue.push(min_idoff, face_to_refine);
    }
  }
  else
  { // face is not restricted
    tet_mesh.unmark(face_idoff, FACE_MARK::RESTRICTED);
    tet_mesh.unmark(mirror_idoff, FACE_MARK::RESTRICTED);
  }
}

/**
 * @brief Checks if a face of the tetrahedral mesh is restricted to a
 * surface patch specified by the domain.
 *
 * A face is considered `restricted` if the dual of the face in the dual
 * Voronoi diagram (usually a segment or a ray) intersects with the surface
 * patch.
 *
 * @param face_idoff The idoff of the face in the tetrahedral mesh.
 * @param surface_patch The index of the surface patch.
 * @param intersection The intersection point between the dual of the face
 * and the surface patch.
 * @return true if the face is restricted to a surface patch, otherwise
 * false.
 */
template <typename T, typename D>
bool DelaunayRefine<T, D>::isFaceRestricted(index_t            face_idoff,
                                            SurfacePatchIndex &surface_patch,
                                            Point3 &intersection) const
{
  // Check if the face is adjacent to two finite tetrahedra
  index_t mirror_idoff = tet_mesh.mirrorFace(face_idoff);

  OMC_EXPENSIVE_ASSERT(tet_mesh.isFiniteTet(face_idoff) ||
                         tet_mesh.isFiniteTet(mirror_idoff),
                       "Invalid face idoff.");

  if (tet_mesh.isFiniteTet(face_idoff) && tet_mesh.isFiniteTet(mirror_idoff))
  {
    // Get the dual segment of the face
    Segment3 dual_seg = tet_mesh.faceDualSegment(face_idoff);

    if (dual_seg.start() == dual_seg.end())
      return false; // degenerate segment

    // Adjust the endpoints of the dual segment to get always the same
    // intersection
    if (canonicalCompare(dual_seg.start(), dual_seg.end()) == Sign::POSITIVE)
      dual_seg = Segment3(dual_seg.end(), dual_seg.start());

    // Check if the dual segment intersects with the surface patch
    if (domain.doIntersectSurface(dual_seg, surface_patch))
    {
      OMC_UNUSED int dimension;
      domain.surfaceIntersection(dual_seg, surface_patch, intersection,
                                 dimension);
      return true;
    }
  }
  else
  {
    // Get the dual ray of the face
    index_t finite_idoff =
      tet_mesh.isFiniteTet(face_idoff) ? face_idoff : mirror_idoff;
    Ray3 dual_ray = tet_mesh.faceDualRay(finite_idoff);

    if (dual_ray.start() == dual_ray.start() + dual_ray.direction())
      return false; // degenerate ray

    if (domain.doIntersectSurface(dual_ray, surface_patch))
    {
      OMC_UNUSED int dimension;
      domain.surfaceIntersection(dual_ray, surface_patch, intersection,
                                 dimension);
      return true;
    }
  }

  return false;
}

template <typename T, typename D>
void DelaunayRefine<T, D>::scanCells()
{
  size_t num_tets = tet_mesh.sizeTets();

  for (index_t tet_id = 0; tet_id < num_tets; tet_id++)
  {
    checkNewCell(TetMesh::toIdOff(tet_id));
  }
}

template <typename T, typename D>
bool DelaunayRefine<T, D>::isCellRefineDone()
{
  // TODO: An external force stop control
  // TODO: Other internal stop conditions
  return cell_queue.empty();
}

template <typename T, typename D>
auto DelaunayRefine<T, D>::processOneCell() -> CellConflictStatus
{
  // Get the next cell to process from the queue ========================

  auto [tet_idoff, cell_to_refine] = cell_queue.top();
  cell_queue.pop();

  OMC_EXPENSIVE_ASSERT(!tet_mesh.isTetDeleted(tet_idoff),
                       "The tetrahedron to process is deleted.");

  // Get the refinement point ===========================================

  Point3 refine_point = tet_mesh.tetDualPoint(tet_idoff);

  // Locate the refinement point in Delaunay tetrahedralization =========

  DelTet del_tet(tet_mesh);

  index_t tet_id    = TetMesh::clipId(tet_idoff);
  int     dimension = -1;
  del_tet.walk(refine_point, tet_id, &dimension);

  if (dimension == 0)
  {
    // The dual point is coincident with a vertex of the tetrahedral mesh.
    return CellConflictStatus::COINCIDENT_VERTEX;
  }

  // Find the conflict zone of the refinement point =====================

  InlinedVector64<index_t> conflict_tets;
  InlinedVector64<index_t> conflict_corners;
  del_tet.conflict(refine_point, /*weight*/ 0, tet_id, conflict_tets,
                   conflict_corners);

  if (conflict_tets.empty())
  {
    // The dual point is hidden by an existing weighted vertex.
    return CellConflictStatus::HIDDEN_POINT;
  }

  // Remove cells in conflict zone from the queue ========================
  for (index_t cf_tet : conflict_tets)
  {
    if (cell_queue.exist(cf_tet))
      cell_queue.remove(cf_tet);
  }

  // Insert the refinement point into the mesh ============================
  del_tet.removeConflicts(conflict_tets);
  index_t new_vtx_vid = newVertex(refine_point, /*weight*/ 0);
  del_tet.filling(new_vtx_vid, conflict_corners);

  // Check new cells in the conflict zone =================================

  for (index_t corner : conflict_corners) // for each corner `c`
  {
    checkNewCell(TetMesh::clipId(tet_mesh.tetNeigh(corner)));
  }

  return CellConflictStatus::OK;
}

template <typename T, typename D>
void DelaunayRefine<T, D>::refineCells()
{
  scanCells();

#ifdef OMC_DELAUNAY_REFINE_PROFILE
  size_t num_inserted_point    = 0;
  size_t num_coincident_vertex = 0;
  size_t num_hidden_point      = 0;
#endif

  while (!isCellRefineDone())
  {
    CellConflictStatus cs = processOneCell();

#ifdef OMC_DELAUNAY_REFINE_PROFILE
    switch (cs)
    {
    case CellConflictStatus::COINCIDENT_VERTEX:
      num_coincident_vertex += 1;
      break;
    case CellConflictStatus::HIDDEN_POINT:
      num_hidden_point += 1;
      break;
    case CellConflictStatus::OK:
      num_inserted_point += 1;
      if (num_inserted_point % 100 == 0)
        std::cout << OMC::format("\rInserted points {}         ",
                                 num_inserted_point);
      break;
    }
#endif
  }

#ifdef OMC_DELAUNAY_REFINE_PROFILE
  Logger::info("\nDelaunayRefine CellRefine profile:");
  Logger::info(OMC::format("  Inserted points: {}", num_inserted_point));
  Logger::info(OMC::format("  Coincident vertex: {}", num_coincident_vertex));
  Logger::info(OMC::format("  Hidden point: {}", num_hidden_point));
#endif

  tet_mesh.removeDeletedTets();
}

/**
 * @brief Check a new cell to determine if it need to be refined.
 *
 * @param tet_idoff The idoff of the tetrahedron in the tetrahedral mesh.
 *
 * @note While checking a new cell, no modifications should be made to its
 * associated conflict zone to ensure consistency and correctness.
 */
template <typename T, typename D>
void DelaunayRefine<T, D>::checkNewCell(index_t tet_idoff)
{
  // ensure unique idoff
  tet_idoff = TetMesh::clipId(tet_idoff);

  if (!tet_mesh.isFiniteTet(tet_idoff))
  { // the tetrahedron is infinite.
    return;
  }

  Point3 dual_point = tet_mesh.tetDualPoint(tet_idoff);
  if (!domain.isInDomain(dual_point))
  { // the dual point is outside the domain.
    return;
  }

  CellQuality cell_quality = criteria.cellQuality(tet_mesh, tet_idoff);

  if (cell_quality)
  { // the tetrahedron should be refined.
    cell_queue.push(tet_idoff, cell_quality);
  }
}

template <typename T, typename D>
index_t DelaunayRefine<T, D>::newVertex(const Point3 &point, NT weight)
{
  // TODO lock in multi-threaded env.

  index_t new_vtx_id = points.size();
  points.push_back(point);
  weights.push_back(weight);
  tet_mesh.newVtx(points.size() - 1);
  return new_vtx_id;
}

template <typename T, typename D>
Sign DelaunayRefine<T, D>::canonicalCompare(const Point3 &p1, const Point3 &p2)
{
  return p1.x() < p2.x()   ? Sign::NEGATIVE
         : p1.x() > p2.x() ? Sign::POSITIVE
         : p1.y() < p2.y() ? Sign::NEGATIVE
         : p1.y() > p2.y() ? Sign::POSITIVE
         : p1.z() < p2.z() ? Sign::NEGATIVE
         : p1.z() > p2.z() ? Sign::POSITIVE
                           : Sign::ZERO;
}

} // namespace OMC
