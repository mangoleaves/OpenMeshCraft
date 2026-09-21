#pragma once

#include "FaceRecover.h"

#include "OpenMeshCraft/Utils/FormatCompat.h"
#include "OpenMeshCraft/Utils/Hashers.h"

#include <iostream>

namespace OMC {

/**
 * @brief Initialize with the Delaunay tetrahedral mesh and the input
 * constrained triangles. Prepare for the constraits recovery.
 * @param _tet_mesh the Delaunay tetrahedral mesh
 * @param _plc the input constrained piecewise linear complex
 */
template <typename Traits>
FaceRecover<Traits>::FaceRecover(std::vector<GPoint3 *> &_verts,
                                 TetMesh &_tet_mesh, PLC &_plc,
                                 ConstrDelTet_Config _config,
                                 ConstrDelTet_Stats *_stats)
  : verts(_verts)
  , tet_mesh(_tet_mesh)
  , plc(_plc)
  , config(_config)
  , stats(_stats)
{
}

template <typename Traits>
void FaceRecover<Traits>::faceRecovery()
{
  // initialize auxiliary data for face recovery.
  {
    v_orient.clear();
    v_cached_orient.clear();
    v_count.clear();
    v_reindex.clear();

    v_orient.resize(verts.size(), Sign::UNCERTAIN);
    v_count.resize(verts.size(), 0);
    v_reindex.resize(verts.size(), InvalidIndex);

    for (index_t tid = 0; tid < tet_mesh.sizeTets(); tid++)
      tet_mesh.unmark(TetMesh::toIdOff(tid), TET_MARK::VISITED);
  }

  // traverse all faces in the PLC to recover the missing faces
  bool   need_recursion  = false;
  size_t recover_succeed = 0, recover_fail = 0;

  do
  {
    need_recursion = false;

    for (index_t i = 0; i < plc.numFaces(); i++)
    {
      std::vector<index_t> tets;
      getTetsIntersectingFace(i, tets);

      if (tets.empty())
        continue; // this face is already recovered.

      bool succeed = false, expanded = false;
      recoverFace_cavityExpanding(i, tets, succeed, expanded);

      // log and output
      if (succeed)
        recover_succeed++;
      else
        recover_fail++;
      if (config.verbose)
        std::cout << OMC::format(
          "\r[OpenMeshCraft CDT] {} faces are recovered. {} faces are missing.",
          recover_succeed, recover_fail);

      // A recovered face may be destroyed by the recovery of another face
      // when expansion is needed.
      // OPT: Record relation between tet face and PLC face to detect destroyed
      // faces more efficiently.
      if (expanded)
        need_recursion = true;
    }
  } while (need_recursion && recover_fail == 0);
  if (config.verbose) // output a new line
    std::cout << std::endl;
  OMC_ASSERT(recover_fail == 0, "Fail to recover {} faces.", recover_fail);
}

/**
 * @brief Finds the tetrahedra intersecting a given PLC face.
 *
 * This function identifies the tetrahedra that intersect with a specified PLC
 * face. It assumes that all PLC segments have already been recovered.
 *
 * @param fid The index of the PLC face to check for intersections.
 * @param tets A vector to store the indices of the intersecting tetrahedra.
 *
 * @pre All PLC segments are recovered.
 * @post The vector `tets` contains the indices of the tetrahedra intersecting
 * the PLC face.
 *
 * @details
 * The function performs the following main steps:
 * 1. Finds part of tetrahedra intersecting the face around flat vertices or
 * around a bounding edge.
 * 2. Expands the search to find all intersected tetrahedra by traversing
 * neighbors.
 * 3. Checks if found tetrahedra penetrate the face and outputs the results.
 */
template <typename Traits>
void FaceRecover<Traits>::getTetsIntersectingFace(index_t               fid,
                                                  std::vector<index_t> &tets)
{
  // =========================================================================
  // # Find the tetrahedra intersecting the PLC face (not all)
  // ## pre-condition: all PLC segments are recovered.
  // ## post-condition: find the tetrahedra intersecting the PLC face (not all).

  // ## Get and initialize data.

  // Get the PLC face
  const PLCFace &face = plc.face(fid);

  // Get the first bounding edge `e0` of the face
  // - `tri0` is the incident input triangle of `e0`.
  // - `rev0` is the orientation of `e0`.
  // - `orig_e0` is the original edge of `e0` and an edge of `tri0`
  index_t        tri0 = InvalidIndex;
  const PLCEdge &e0   = plc.boundingEdge(face, 0, &tri0);
  const PLCEdge &orig_e0 =
    is_valid_idx(e0.ancestor_id) ? plc.edge(e0.ancestor_id) : e0;

  // Get the endpoints of `e0` and `orig_e0`
  index_t e0p0 = e0.ep0(), e0p1 = e0.ep1();
  index_t oe0p0 = orig_e0.ep0(), oe0p1 = orig_e0.ep1();
  // Get the three vertices of `tri0`
  index_t tri_v[3] = {plc.triVtx(tri0, 0), plc.triVtx(tri0, 1),
                      plc.triVtx(tri0, 2)};
  int     n_max    = -1; // the maximal component of triangle normal.
  Sign    tri_ori  = Sign::UNCERTAIN; // triangle vertex orientation.
  OMC_EXPENSIVE_ASSERT(
    (oe0p0 == tri_v[0] || oe0p0 == tri_v[1] || oe0p0 == tri_v[2]) &&
      (oe0p1 == tri_v[0] || oe0p1 == tri_v[1] || oe0p1 == tri_v[2]),
    "The original edge is not an edge of the input triangle.");

#ifdef OMC_ENABLE_EXPENSIVE_ASSERT
  // check if all bounding edges are recovered.
  for (const BoundingEdge &be : face.bounding_edges)
  {
    for (index_t i = be.range.start; i < be.range.start + be.range.size; i++)
    {
      const PLCEdge &edge = plc.subEdge(i);
      OMC_ASSERT(tet_mesh.edgeExists(edge.ep0(), edge.ep1()),
                 "The edge is not recovered.");
    }
  }
#endif

  // Adjust the order of the vertices of `tri0` to make sure that `oe0p0` and
  // `oe0p1` are the first two vertices of `tri0`.
  tri_v[2] += tri_v[0] + tri_v[1];
  tri_v[0] = oe0p0, tri_v[1] = oe0p1, tri_v[2] = tri_v[2] - oe0p0 - oe0p1;

  // Initialize four vertices of adajcent tetrahedron.
  // First essemble the two endpoints of the first bounding edge.
  index_t                  tet_v[4] = {e0p0, e0p1, InvalidIndex, InvalidIndex};
  // Get the incident tets of the first bounding edge
  InlinedVector64<index_t> edge_incident_tets;
  tet_mesh.ET(tet_v[0], tet_v[1], edge_incident_tets);

  // ## If the face just has one triangle and is not splitted,
  //    check if the face is already recovered.
  if (face.triangles.size() == 1 && face.bounding_vertices.size() == 3 &&
      face.flat_vertices.empty())
  {
    // The opposite vertex to the edge `e0`
    index_t opp_v = (tri_v[0] != tet_v[0] && tri_v[0] != tet_v[1])   ? tri_v[0]
                    : (tri_v[1] != tet_v[0] && tri_v[1] != tet_v[1]) ? tri_v[1]
                                                                     : tri_v[2];
    // Check if the opposite vertex is in the incident tets
    for (index_t tet_idoff : edge_incident_tets)
    {
      if (tet_mesh.tetHasVertex(tet_idoff, opp_v))
      {         // The tetrahedron contains the opposite vertex
        return; // The face is already recovered, exit.
      }
    }
  }

  // ## initialize data for further check

  // initialize vertex orientation and increase count
  // v_orient is initialized to UNCERTAIN and v_count is initialized to zero.
  // so, remember to reset them before and after calling this function.
  for (index_t vid : face.bounding_vertices)
  {
    v_orient[vid] = Sign::ZERO;
    v_count[vid]++;
  }
  for (index_t vid : face.flat_vertices)
  {
    v_orient[vid] = Sign::ZERO;
  }

  auto isVtxBounding = [this](index_t vid) { return v_count[vid] > 0; };

  // The vector to store intersected tetrahedra
  InlinedVector64<index_t> B;

  // ## Find the intersected tetrahedra around flat vertices
  if (!face.flat_vertices.empty())
  {
    // find in VT
    for (index_t vid : face.flat_vertices)
    {
      InlinedVector64<index_t> tmp_B;
      tet_mesh.VT(vid, tmp_B);
      B.insert(B.end(), tmp_B.begin(), tmp_B.end());
    }
    // unique
    B.erase(std::unique(B.begin(), B.end()), B.end());
  }

  // ## Find a tet `t0` in ET(e) intersecting the face interior.
  if (!B.empty())               // if `B` is not empty.
    edge_incident_tets.clear(); // skip finding.
  for (index_t tet_idoff : edge_incident_tets)
  {
    // Get the edge opposite to `e0`
    // (`e0` is stored in tet_v[0,1], the oppo edge is stored in tet_v[3,4])
    tet_mesh.oppoEdge(tet_idoff, tet_v[0], tet_v[1], tet_v[2], tet_v[3]);

    // Calculate the orientation of two endpoints of the opposite edge with
    // respect to the plane defined by the PLC face.
    //                        ========== plane ===========  query point
    Sign ot2 = orient3dCached(tri_v[0], tri_v[1], tri_v[2], tet_v[2]);
    Sign ot3 = orient3dCached(tri_v[0], tri_v[1], tri_v[2], tet_v[3]);

    OMC_EXPENSIVE_ASSERT(ot2 != Sign::ZERO || ot3 != Sign::ZERO,
                         "Degenerate tetrahedron.");

    if (ot2 != Sign::ZERO && ot2 == ot3) // totally above or below the plane
      continue;                          // ==> no intersection, skip
    else if (ot2 == Sign::ZERO || ot3 == Sign::ZERO) // one endpoint is coplanar
    {
      index_t copl_vid = ot2 == Sign::ZERO ? tet_v[2] : tet_v[3];

      if (!isVtxBounding(copl_vid)) // this point is not a bounding vertex,
        continue;                   // so it is outside, skip it.

      const GPoint3 &tri_p0 = point(tri_v[0]), &tri_p1 = point(tri_v[1]),
                    &tri_p2 = point(tri_v[2]), &copl_p = point(copl_vid);

      if (n_max == -1)
      {
        n_max = MaxCompInTriNormal()(
          AsEP()(tri_p0).data(), AsEP()(tri_p1).data(), AsEP()(tri_p2).data());
        tri_ori = OrientOn2D()(tri_p0, tri_p1, tri_p2, n_max);
      }

      Sign copl_ori = OrientOn2D()(tri_p0, tri_p1, copl_p, n_max);
      OMC_EXPENSIVE_ASSERT(tri_ori != Sign::ZERO && copl_ori != Sign::ZERO,
                           "Degenerate triangle.");

      if (tri_ori == copl_ori) // coplanar point is inside the face
      {
        B.push_back(tet_idoff); // find an intersected tetrahedron
        break;                  // exit loop
      }
    }
    else if (ot2 != ot3 && segCrossesFace(tet_v[2], tet_v[3], face))
    {                         // the opposite edge crosses the face
      B.push_back(tet_idoff); // find an intersected tetrahedron
      break;                  // exit loop
    }
  }

  // mark found tetrahedra as visited
  for (index_t tet_idoff : B)
    tet_mesh.mark(tet_idoff, TET_MARK::VISITED);

  OMC_EXPENSIVE_ASSERT(!B.empty(), "No intersected tetrahedra found.");

  // =========================================================================
  // # Find all the intersected tetrahedra by expanding
  //   from the found tetrahedra.
  // ## pre-condition: part intersected tetrahedra are found
  // ## post-condition: all intersected tetrahedra are found

  auto collectTet = [this, &B](index_t tet_idoff)
  {
    B.push_back(tet_idoff);
    tet_mesh.mark(tet_idoff, TET_MARK::VISITED);
  };

  phmap::flat_hash_set<IPair, hash<IPair>> vtx_adj_set;
  plc.buildBoundingVtxAdjSet(face, vtx_adj_set);

  auto areVertsAdj = [&vtx_adj_set](index_t v0, index_t v1) -> bool
  { return vtx_adj_set.find(unique_pair(v0, v1)) != vtx_adj_set.end(); };

  for (index_t k = 0; k < B.size(); k++)
  {
    index_t tet_idoff = B[k];

    // Traverse the four neighbors of the tetrahedron
    for (index_t j = 0; j < 4; j++)
    {
      index_t nb_idoff = TetMesh::clipId(tet_mesh.tetNeigh(tet_idoff + j));
      // If the neighbor is not finite or visited, skip
      if (!tet_mesh.isFiniteTet(nb_idoff) ||
          tet_mesh.isMarked(nb_idoff, TET_MARK::VISITED))
        continue;
      // The vertices of the common triangle between two tets.
      index_t cv[3]  = {tet_mesh.tetNode(tet_idoff + ((j + 1) & 3)),
                        tet_mesh.tetNode(tet_idoff + ((j + 2) & 3)),
                        tet_mesh.tetNode(tet_idoff + ((j + 3) & 3))};
      // Is the common vertex a bounding vertex?
      bool    bv[3]  = {isVtxBounding(cv[0]), isVtxBounding(cv[1]),
                        isVtxBounding(cv[2])};
      // POSITIVE = 1, ZERO = 0, NEGATIVE = -1
      int     o3d[3] = {
        static_cast<int>(orient3dCached(tri_v[0], tri_v[1], tri_v[2], cv[0])),
        static_cast<int>(orient3dCached(tri_v[0], tri_v[1], tri_v[2], cv[1])),
        static_cast<int>(orient3dCached(tri_v[0], tri_v[1], tri_v[2], cv[2]))};

      if (bv[0] && bv[1] && bv[2])
      { // All bounding vertices, this common face is on the PLC face.
        // Thus, the opposite tetrahedron is intersected.
        collectTet(nb_idoff);
      }
      else if (bv[0] && bv[1])
      { // Only two bounding vertices.
        // The edge connecting the two bounding vertices is inside the PLC face
        // if the two bounding vertices are NOT adjacent, otherwise the edge is
        // a bounding edge.
        // The opposite tetrahedron is intersected in the former case.
        if (!areVertsAdj(cv[0], cv[1]))
          collectTet(nb_idoff);
      }
      else if (bv[1] && bv[2])
      { // Same as the above case
        if (!areVertsAdj(cv[1], cv[2]))
          collectTet(nb_idoff);
      }
      else if (bv[2] && bv[0])
      { // Same as the above case
        if (!areVertsAdj(cv[2], cv[0]))
          collectTet(nb_idoff);
      }
      else if (bv[0])
      { // Only one bounding vertex.
        // Condition: The edge opposite to the bounding vertex in the common
        // triangle crosses the PLC face
        // ==> The common triangle intersects the PLC face
        // ==> The opposite tetrahedron is intersected.
        if (o3d[1] * o3d[2] < 0)
          collectTet(nb_idoff);
      }
      else if (bv[1])
      { // Same as the above case
        if (o3d[2] * o3d[0] < 0)
          collectTet(nb_idoff);
      }
      else if (bv[2])
      { // Same as the above case
        if (o3d[0] * o3d[1] < 0)
          collectTet(nb_idoff);
      }
      else
      { // No bounding vertex.
        // Condition: Any edge of the common triangle crosses the PLC face
        // ==> The common triangle intersects the PLC face
        // ==> The opposite tetrahedron is intersected.
        if (o3d[1] * o3d[2] < 0 || o3d[2] * o3d[0] < 0 || o3d[0] * o3d[1] < 0)
          collectTet(nb_idoff);
      }
    }
  }

  // Check if found tetrahedra penetrate the face and output
  for (index_t tet_idoff : B)
  {
    tet_mesh.unmark(tet_idoff, TET_MARK::VISITED);
    index_t v[4] = {
      tet_mesh.tetNode(tet_idoff + 0), tet_mesh.tetNode(tet_idoff + 1),
      tet_mesh.tetNode(tet_idoff + 2), tet_mesh.tetNode(tet_idoff + 3)};
    // POSITIVE = 1, ZERO = 0, NEGATIVE = -1
    int ov[4] = {
      static_cast<int>(orient3dCached(tri_v[0], tri_v[1], tri_v[2], v[0])),
      static_cast<int>(orient3dCached(tri_v[0], tri_v[1], tri_v[2], v[1])),
      static_cast<int>(orient3dCached(tri_v[0], tri_v[1], tri_v[2], v[2])),
      static_cast<int>(orient3dCached(tri_v[0], tri_v[1], tri_v[2], v[3]))};

    if (!((ov[0] >= 0 && ov[1] >= 0 && ov[2] >= 0 && ov[3] >= 0) ||
          (ov[0] <= 0 && ov[1] <= 0 && ov[2] <= 0 && ov[3] <= 0)))
    { // the tet penetrate the face, output
      OMC_EXPENSIVE_ASSERT(tetIntersectsFace(tet_idoff, face), "Wrong tet.");
      tets.push_back(tet_idoff);
    }
    // otherwise, the tet is just visiting the face, and the face is partially
    // recovered.
  }

  // Clear auxiliary orientation and count
  for (index_t vid : v_cached_orient)
    v_orient[vid] = Sign::UNCERTAIN;
  v_cached_orient.clear();
  for (index_t vid : face.bounding_vertices)
  {
    v_count[vid]  = 0;
    v_orient[vid] = Sign::UNCERTAIN;
  }
  for (index_t vid : face.flat_vertices)
    v_orient[vid] = Sign::UNCERTAIN;

#ifdef OMC_ENABLE_EXPENSIVE_ASSERT
  // check if all bounding and flat vertices are traversed.
  for (index_t tet_idoff : B)
  {
    index_t v[4] = {
      tet_mesh.tetNode(tet_idoff + 0), tet_mesh.tetNode(tet_idoff + 1),
      tet_mesh.tetNode(tet_idoff + 2), tet_mesh.tetNode(tet_idoff + 3)};
    v_count[v[0]] = v_count[v[1]] = v_count[v[2]] = v_count[v[3]] = 1;
  }
  for (index_t vid : face.bounding_vertices)
  {
    OMC_ASSERT(v_count[vid], "missing bounding vertex in cavity");
  }
  for (index_t vid : face.flat_vertices)
  {
    OMC_ASSERT(v_count[vid], "missing flat vertex in cavity");
  }
  for (index_t tet_idoff : B)
  {
    index_t v[4] = {
      tet_mesh.tetNode(tet_idoff + 0), tet_mesh.tetNode(tet_idoff + 1),
      tet_mesh.tetNode(tet_idoff + 2), tet_mesh.tetNode(tet_idoff + 3)};
    v_count[v[0]] = v_count[v[1]] = v_count[v[2]] = v_count[v[3]] = 0;
  }
#endif
}

template <typename Traits>
void FaceRecover<Traits>::recoverFace_cavityExpanding(
  index_t fid, const std::vector<index_t> &tets, bool &succeed, bool &expanded)
{
  // Get the PLC face
  const PLCFace &face = plc.face(fid);

  index_t tid = face.triangles[0];
  index_t tv0 = plc.triVtx(tid, 0), tv1 = plc.triVtx(tid, 1),
          tv2 = plc.triVtx(tid, 2);

  { // Initialize vertex orientation and count
    for (index_t vid : face.bounding_vertices)
    {
      v_orient[vid] = Sign::ZERO;
      v_count[vid]  = 0;
    }
    for (index_t vid : face.flat_vertices)
    {
      v_orient[vid] = Sign::ZERO;
      v_count[vid]  = 0;
    }

    for (index_t tet_idoff : tets)
    {
      for (index_t j = 0; j < 4; j++)
      {
        index_t vid  = tet_mesh.tetNode(tet_idoff + j);
        v_count[vid] = 0;
        orient3dCached(tv0, tv1, tv2, vid);
      }
    }
  }

  // =========================================================================
  // # Divide cavity into top and bottom half cavity
  // ## pre-condition: all tetrahedra intersecting the PLC face are found.
  // ## post-condition: cavity vertices and boundary faces are divided.

  // visit all tetrahedra in `tets` for dividing cavity
  for (index_t tet_idoff : tets)
    tet_mesh.mark(tet_idoff, TET_MARK::VISITED);

  // divide cavity into top and bottom half cavity.

  // indices of the vertices of the two half cavities.
  // - vertices on the PLC face belong to both sides.
  InlinedVector64<index_t> top_vertices, bottom_vertices;

  // boundary faces of the two half cavities.
  // - each boundary face has a corresponding corner in an adjacent unvisited
  // tet.
  //   (see explanation of `corner` in TetMesh::tet_node.(3))
  // - store the indices of corners in below vectors.
  // - faces on the PLC face are ignored.
  InlinedVector64<index_t> top_faces, bottom_faces;

  for (index_t tet_idoff : tets)
  {
    // divide vertices into top and bottom cavity.
    for (index_t j = 0; j < 4; j++)
    {
      index_t vid = tet_mesh.tetNode(tet_idoff + j);
      if (v_count[vid]) // skip traversed vertices
        continue;

      v_count[vid] = 1;
      Sign ov      = v_orient[vid];
      OMC_EXPENSIVE_ASSERT(is_sign_reliable(ov), "uninitialized orientation.");
      if (ov >= Sign::ZERO)
        top_vertices.push_back(vid);
      if (ov <= Sign::ZERO)
        bottom_vertices.push_back(vid);
    }
    // divide boundary faces into top and bottom cavity.
    for (index_t j = 0; j < 4; j++)
    {
      // the neighboring corner
      index_t neigh_idoff = tet_mesh.tetNeigh(tet_idoff + j);
      if (tet_mesh.isMarked(neigh_idoff, TET_MARK::VISITED))
        continue;
      // Get the vertices of the common face, which is also the boundary face of
      // the cavity.
      index_t v0, v1, v2;
      tet_mesh.faceVertices(neigh_idoff, v0, v1, v2);

      bool top = v_orient[v0] >= Sign::ZERO && v_orient[v1] >= Sign::ZERO &&
                 v_orient[v2] >= Sign::ZERO;
      bool bottom = v_orient[v0] <= Sign::ZERO && v_orient[v1] <= Sign::ZERO &&
                    v_orient[v2] <= Sign::ZERO;
      OMC_EXPENSIVE_ASSERT(
        (top || bottom) && (!top || !bottom),
        "(1) Face belongs to neither top nor bottom half cavity."
        "(2) Face belongs to both top and bottom half cavity.");
      // check if the face belongs to the top or bottom half cavity.
      if (top)
        top_faces.push_back(neigh_idoff);
      else if (bottom)
        bottom_faces.push_back(neigh_idoff);
    }
  }

  // Sort vertices and faces so that we can build a sequential map from global
  // mesh to local meshed cavity. Ensure that the relative order of vertices are
  // same in both global and local meshes, because the simulation of simplicity
  // relies on this.
  std::sort(top_vertices.begin(), top_vertices.end());
  std::sort(bottom_vertices.begin(), bottom_vertices.end());
  std::sort(top_faces.begin(), top_faces.end());
  std::sort(bottom_faces.begin(), bottom_faces.end());

  // clear vertex count
  for (index_t vid : top_vertices)
    v_count[vid] = 0;
  for (index_t vid : bottom_vertices)
    v_count[vid] = 0;
  // clear visited marks
  for (index_t tet_idoff : tets)
    tet_mesh.unmark(tet_idoff, TET_MARK::VISITED);
  // v_orient is kept for further use.

  // =========================================================================
  // # Delaunay tetrahedralize the top and bottom half cavity.
  //   The cavity may be expanded to recover missing cavity boundary faces.
  // ## pre-condition: cavity is properly divided.
  // ## post-condition: two half cavities are properly tetrahedralized.

  // Vertices of the top and bottom half cavity, used for Delaunay
  // tetrahedralization
  std::vector<GPoint3 *> top_vps, bottom_vps;
  top_vps.reserve(2 * top_vertices.size());
  bottom_vps.reserve(2 * bottom_vertices.size());

  // TetMeshes of the top and bottom half cavity
  std::unique_ptr<TetMesh> top_mesh, bottom_mesh;
  // new tets after expanding, waiting to be removed.
  std::vector<index_t>     tets_to_remove;
  // missing boundary face during meshing and embedding the cavity
  index_t                  missing_face = InvalidIndex;

  // cavity is ok when it does not cross the PLC plane.
  bool cavity_ok = true;

  // tetrahedralize the half cavity.
  // if any cavity face is missing, expand the cavity to recover the missing
  // face, until no missing face.
  while (cavity_ok)
  {
    // load global vertices to local mesh, keep the order of vertices.
    top_vps.clear();
    for (index_t vid : top_vertices)
      top_vps.push_back(&point(vid));
    // Delaunay tetrahedralize the cavity
    top_mesh = std::make_unique<TetMesh>(top_vps);
    DelTet top_dt(*top_mesh);
    top_dt.tetrahedralize();
    // check if any cavity face is missing
    if (cavityHasMissingFace(*top_mesh, top_vertices, top_faces, missing_face))
    {
      // expand the cavity to recover the missing boundary face.
      index_t new_tet = InvalidIndex, new_vertex = InvalidIndex;
      expandCavity(top_vertices, top_faces, missing_face, new_tet, new_vertex);
      // mark the tetrahedron to remove
      tets_to_remove.push_back(new_tet);
      if (is_valid_idx(new_vertex))
      {
        // calculate orientation for new vertex
        v_orient[new_vertex] = orient3dCached(tv0, tv1, tv2, new_vertex);
        // check if the cavity is still ok after expanding
        if (v_orient[new_vertex] < Sign::ZERO)
          cavity_ok = false;
      }
      // mark the cavity is expanded
      expanded = true;
    }
    else // no missing face found, tethraheralization is done.
      break;
  }

  while (cavity_ok)
  {
    // load global vertices to local mesh, keep the order of vertices.
    bottom_vps.clear();
    for (index_t vid : bottom_vertices)
      bottom_vps.push_back(&point(vid));
    // Delaunay tetrahedralize the cavity
    bottom_mesh = std::make_unique<TetMesh>(bottom_vps);
    DelTet bottom_dt(*bottom_mesh);
    bottom_dt.tetrahedralize();
    // check if any cavity face is missing
    if (cavityHasMissingFace(*bottom_mesh, bottom_vertices, bottom_faces,
                             missing_face))
    {
      // expand the cavity to recover the missing boundary face.
      index_t new_tet = InvalidIndex, new_vertex = InvalidIndex;
      expandCavity(bottom_vertices, bottom_faces, missing_face, new_tet,
                   new_vertex);
      // mark the tetrahedron to remove
      tets_to_remove.push_back(new_tet);
      if (is_valid_idx(new_vertex))
      {
        // calculate orientation for new vertex
        v_orient[new_vertex] = orient3dCached(tv0, tv1, tv2, new_vertex);
        // add the vertex if it is not in the cavity
        if (v_orient[new_vertex] > Sign::ZERO)
          cavity_ok = false;
      }
      // mark the cavity is expanded
      expanded = true;
    }
    else // no missing face found, tethraheralization is done.
      break;
  }

  if (cavity_ok)
  { // Really modify the global mesh.
    size_t n_top_faces = top_faces.size(), n_bottom_faces = bottom_faces.size();
    // first, embed the tetrahedralization of the cavity to the global mesh
    embedMeshedCavity(*top_mesh, top_vertices, top_faces, bottom_faces);
    embedMeshedCavity(*bottom_mesh, bottom_vertices, bottom_faces, top_faces);

    OMC_ASSERT(n_top_faces == top_faces.size(),
               "The number of top faces is not consistent.");
    OMC_ASSERT(n_bottom_faces < bottom_faces.size(),
               "The number of bottom faces is not consistent.");

#ifdef OMC_ENABLE_EXPENSIVE_ASSERT
    { // check validity of the embedded part
      size_t  new_finite_tets = top_mesh->sizeTets() + bottom_mesh->sizeTets();
      DelTet  DT(tet_mesh);
      index_t tet_idoff_start = (tet_mesh.sizeTets() - new_finite_tets) * 4;
      index_t tet_idoff_end   = tet_mesh.sizeTets() * 4;
      for (index_t tet_idoff = tet_idoff_start; tet_idoff < tet_idoff_end;
           tet_idoff += 4)
      {
        OMC_ASSERT(DT.verifyVolume(tet_idoff), "Negative volume.");
        OMC_ASSERT(DT.verifyNeighbor(tet_idoff), "Connectivity error.");
      }
    }
#endif

    // then, remove the tetrahedra of cavity
    for (index_t idoff : tets) // original part of the cavity
      tet_mesh.markTetAsDeleted(idoff);
    for (index_t idoff : tets_to_remove) // expanded part of the cavity
      tet_mesh.markTetAsDeleted(idoff);
    tet_mesh.removeDeletedTets();
  }
  succeed = cavity_ok;

#if 0
  if (!cavity_ok)
  {
    std::fstream fout;
    fout.open("./data/test_output/top_cavity.obj", std::ios::out);
    for (index_t vid : top_vertices)
    {
      EPoint3 p = ToEP()(point(vid));
      fout << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
    fout.close();
    fout.open("./data/test_output/bottom_cavity.obj", std::ios::out);
    for (index_t vid : bottom_vertices)
    {
      EPoint3 p = ToEP()(point(vid));
      fout << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
    fout.close();
  }
#endif

  { // Clear cached vertex orientation
    for (index_t vid : v_cached_orient)
      v_orient[vid] = Sign::UNCERTAIN;
    v_cached_orient.clear();
    for (index_t vid : face.bounding_vertices)
      v_orient[vid] = Sign::UNCERTAIN;
    for (index_t vid : face.flat_vertices)
      v_orient[vid] = Sign::UNCERTAIN;
  }
}

/**
 * @brief Checks if the cavity has any missing faces.
 *
 * This function determines if any face in the given list of faces is missing
 * in the local mesh. It maps the global indices of vertices to local indices,
 * traverses all faces to check for their existence in the local mesh, and
 * identifies any missing face.
 *
 * @param local_mesh The local tetrahedral mesh to check against.
 * @param vertices A vector of vertex indices in the global mesh.
 * @param faces A vector of face indices in the global mesh.
 * @param missing_face Output parameter that will hold the index of the missing
 * face if any.
 * @return True if there is a missing face, false otherwise.
 * @note `v_reindex` is used, not thread safe.
 */
template <typename Traits>
bool FaceRecover<Traits>::cavityHasMissingFace(
  const TetMesh &local_mesh, const InlinedVector64<index_t> &vertices,
  const InlinedVector64<index_t> &faces, index_t &missing_face)
{
  missing_face = InvalidIndex;
  // map global index in `tet_mesh` to local index in `local_mesh`
  for (index_t i = 0; i < vertices.size(); i++)
    v_reindex[vertices[i]] = i;

  // traverse all faces to check if any one is missing
  for (index_t idoff : faces)
  {
    // map global face to local face
    index_t lv0, lv1, lv2;
    tet_mesh.faceVertices(idoff, lv0, lv1, lv2);
    lv0 = v_reindex[lv0], lv1 = v_reindex[lv1], lv2 = v_reindex[lv2];
    // check if local face exists
    if (!local_mesh.faceExists(lv0, lv1, lv2))
    {
      missing_face = idoff;
      break;
    }
  }

  // reset the map in `v_reindex`
  for (index_t i = 0; i < vertices.size(); i++)
    v_reindex[vertices[i]] = InvalidIndex;
  return is_valid_idx(missing_face);
}

/**
 * @brief Expand the cavity to recover the missing boundary face.
 *
 * The missing face is adjacent to two tetrahedra.
 * One tetrahedron belongs to the cavity while the other does not.
 * We expand the cavity by adding the new tetrahedron into the cavity.
 *
 * @param [in] vertices Vertices of the half cavity.
 * @param [in] faces Boundary faces of the half cavity.
 * @param [in] missing_face The missing boundary face when meshing the cavity.
 * @param [out] new_tet The new tetrahedron of the cavity after expanding.
 * @param [out] new_vertex The new vertex of the cavity after expanding.
 */
template <typename Traits>
void FaceRecover<Traits>::expandCavity(InlinedVector64<index_t> &vertices,
                                       InlinedVector64<index_t> &faces,
                                       index_t missing_face, index_t &new_tet,
                                       index_t &new_vertex)
{
  OMC_EXPENSIVE_ASSERT(tet_mesh.isFiniteTet(missing_face),
                       "The missing face belongs to an infinite tet.");

  // (1) Get the `new_tet`
  // The new tetrahedron containing `missing_face` will be added to the cavity.
  // `missing_face` corresponds to a corner in `new_tet`
  new_tet = TetMesh::clipId(missing_face);

  // (2) Remove `missing_face` from `faces`
  auto faces_iter = std::lower_bound(faces.begin(), faces.end(), missing_face);
  OMC_EXPENSIVE_ASSERT(faces_iter != faces.end() && *faces_iter == missing_face,
                       "The missing face is not in the boundary faces.");
  faces.erase(faces_iter);

  // (3) Add the `new_tet`.
  // For each corner in `new_tet` except the one corresponding to (corr. to)
  // `missing_face`, check if the face corr. to the corner needs to be
  // expanded.
  size_t added_corner_count = 0;
  for (index_t i = 0; i < 4; i++)
  {
    index_t corner = new_tet + i;
    if (corner == missing_face) // Skip the corner corr. to the missing face
      continue;
    // Check if the face corr. to this corner is already in the boundary faces
    faces_iter = std::lower_bound(faces.begin(), faces.end(), corner);
    if (faces_iter != faces.end() && *faces_iter == corner)
    { // If the corner is in the boundary faces, remove it due to expanding.
      faces.erase(faces_iter);
    }
    else
    { // Otherwise, add the face corr. to this corner of the cavity.
      // (we actually add the opposite corner.)
      index_t corner_to_add = tet_mesh.tetNeigh(corner);
      auto    it = std::lower_bound(faces.begin(), faces.end(), corner_to_add);
      faces.insert(it, corner_to_add);
      added_corner_count++;
    }
  }

  // (4) A new vertex is possibly added to the cavity.
  // If three corners are added, add their common vertex to the cavity. The
  // common vertex is the opposite vertex to the missing face in the `new_tet`.
  if (added_corner_count == 3)
  {
    new_vertex = tet_mesh.tetNode(missing_face);
    auto vit   = std::lower_bound(vertices.begin(), vertices.end(), new_vertex);
    if (vit == vertices.end() || *vit != new_vertex)
      vertices.insert(vit, new_vertex);
    else
      new_vertex = InvalidIndex;
  }
  else
    new_vertex = InvalidIndex;
}

/**
 * @brief Embeds a meshed cavity into the global tetrahedral mesh.
 *
 * This function takes a local tetrahedral mesh representing a cavity and embeds
 * it into a global tetrahedral mesh. It performs the following steps:
 *
 * 1. Maps global corners to local ones and identifies boundary corners.
 *
 * 2. Classifies local tetrahedra as inside or outside based on boundary
 * corners.
 *
 * 3. Removes outside tetrahedra from the local mesh.
 *
 * 4. Embeds the local mesh into the global mesh.
 *
 * 5. Collects base faces (corners) that are boundary faces of the cavity.
 *
 * @tparam Traits The traits class providing necessary types and constants.
 * @param local_mesh The local tetrahedral mesh representing the cavity.
 * @param vertices A vector of vertex indices in the global mesh.
 * @param faces A vector of face indices in the global mesh.
 * @param base A vector to store the base faces (corners) of the cavity.
 */
template <typename Traits>
void FaceRecover<Traits>::embedMeshedCavity(
  TetMesh &local_mesh, const InlinedVector64<index_t> &vertices,
  const InlinedVector64<index_t> &faces, InlinedVector64<index_t> &base)
{
  std::vector<uint8_t> corner_is_boundary(local_mesh.sizeTets() * 4, false);

  typedef struct BoundaryCornerPair
  {
    index_t c0;  // one corner in the local mesh.
    index_t c1;  // the other corner in the local mesh.
    index_t bnd; // boundary corner in the global mesh.
  } BCP;
  InlinedVector64<BCP> bcpairs;

  // (1) Map the global corners (with corresponding boundary face) to local ones

  // Build a map from global vertices to local vertices
  for (index_t i = 0; i < vertices.size(); i++)
    v_reindex[vertices[i]] = i;

  // Find the boundary corners in the local mesh
  for (index_t idoff : faces)
  {
    // map global face to local face
    index_t lv0, lv1, lv2;
    tet_mesh.faceVertices(idoff, lv0, lv1, lv2);
    lv0 = v_reindex[lv0], lv1 = v_reindex[lv1], lv2 = v_reindex[lv2];
    // Get the two corners corresponding to the local face
    index_t c0, c1;
    local_mesh.faceCorners(lv0, lv1, lv2, c0, c1);
    // Mark them as boundary corners and record the pair
    corner_is_boundary[c0] = true;
    corner_is_boundary[c1] = true;
    bcpairs.push_back({c0, c1, idoff});
  }

  // Reset the map
  for (index_t i = 0; i < vertices.size(); i++)
    v_reindex[vertices[i]] = InvalidIndex;

  // (2) Classify the local tetrahedra to inside/outside based on the boundary
  // corners. We will embed inner tetrahedra into the global mesh. Note that
  // not all tetrahedra or vertices will be embedded.

  // Find an infinite tet in the local mesh that has at least a vertex not on
  // the PLC face. Start from this tet to classify the tetrahedra.
  index_t start_tet = InvalidIndex;
  for (index_t tid = 0; tid < local_mesh.sizeTets(); tid++)
  {
    index_t tet_idoff = TetMesh::toIdOff(tid);
    if (!local_mesh.isFiniteTet(tet_idoff) &&
        (v_orient[vertices[local_mesh.tetNode(tet_idoff)]] != Sign::ZERO ||
         v_orient[vertices[local_mesh.tetNode(tet_idoff + 1)]] != Sign::ZERO ||
         v_orient[vertices[local_mesh.tetNode(tet_idoff + 2)]] != Sign::ZERO))
    {
      start_tet = tet_idoff;
      break;
    }
  }
  OMC_ASSERT(is_valid_idx(start_tet), "No start tet found.");

  local_mesh.classifyInOut(corner_is_boundary, start_tet);

  // (3) Remove the outside tetrahedra, since we only embed the inside
  // tetrahedra to the global mesh.
  std::vector<index_t> tet_reindex(local_mesh.sizeTets());
  std::iota(tet_reindex.begin(), tet_reindex.end(), 0);

  std::vector<uint8_t> remain_vtx(vertices.size(), false);

  // Depart the outside tetrahedra from the local mesh
  for (index_t tid = 0; tid < local_mesh.sizeTets(); tid++)
  {
    index_t idoff = TetMesh::toIdOff(tid);
    if (local_mesh.isMarked(idoff, TET_MARK::OUTSIDE))
    {
      // depart the outside tetrahedron from its inside neighbors.
      for (index_t j = 0; j < 4; j++)
      {
        index_t neigh_idoff = local_mesh.tetNeigh(idoff + j);
        if (local_mesh.isMarked(neigh_idoff, TET_MARK::INSIDE))
          local_mesh.tetNeigh(neigh_idoff) = InvalidIndex;
      }
    }
    else
    {
      OMC_EXPENSIVE_ASSERT(local_mesh.isMarked(idoff, TET_MARK::INSIDE),
                           "Wrong mark.");
      // remaining vertices only relate to the inside tetrahedra.
      for (index_t j = 0; j < 4; j++)
      {
        index_t vid            = local_mesh.tetNode(idoff + j);
        local_mesh.incTet(vid) = tid;
        remain_vtx[vid]        = true;
      }
    }
  }

  // To remove the outside tetrahedra, we move the inside tetrahedra from the
  // end to replace the outside tetrahedra at the beginning.
  index_t last_tid = local_mesh.sizeTets() - 1;
  for (index_t tid = 0; tid < last_tid; tid++)
  {
    index_t idoff = TetMesh::toIdOff(tid);
    if (!local_mesh.isMarked(idoff, TET_MARK::OUTSIDE))
      continue;

    // move to the last inside tetrahedron
    while (local_mesh.isMarked(TetMesh::toIdOff(last_tid), TET_MARK::OUTSIDE))
      last_tid--; // there is at least one inside tet.
    if (tid >= last_tid)
      break; // all outside tets are removed, break the loop.

    index_t last_idoff = TetMesh::toIdOff(last_tid);
    // Move the last inside tetrahedron to the current position.
    for (index_t j = 0; j < 4; j++)
    {
      TetMesh &LM   = local_mesh;
      // Update the node information.
      index_t &node = LM.tetNode(idoff + j);
      node          = LM.tetNode(last_idoff + j);

      // Update the neighbor information.
      index_t neigh          = LM.tetNeigh(last_idoff + j);
      LM.tetNeigh(idoff + j) = neigh;
      if (is_valid_idx(neigh))
        LM.tetNeigh(neigh) = idoff + j;

      // Update the incident tetrahedron information.
      LM.incTet(node) = tid;
    }
    // Update the mark.
    local_mesh.template getMark<TET_MARK>(idoff) =
      local_mesh.template getMark<TET_MARK>(last_idoff);
    // Update the reindex map.
    tet_reindex[last_tid] = tid;
    tet_reindex[tid]      = last_tid;
    // move the last tetrahedron.
    last_tid--;
  }
  while (local_mesh.isMarked(TetMesh::toIdOff(last_tid), TET_MARK::OUTSIDE))
    last_tid--;
  // Finally, all outside tetrahedra are removed by the above loop.
  local_mesh.resizeTets(last_tid + 1);

  // (4) Embed the local mesh to the global

  // record size before embedding
  size_t n_local_tets  = local_mesh.sizeTets();
  size_t n_local_verts = local_mesh.sizeVerts();
  size_t n_global_tets = tet_mesh.sizeTets();

  // update the node in the local mesh to the corresponding global node
  for (index_t tid = 0; tid < n_local_tets; tid++)
  {
    for (index_t idoff = TetMesh::toIdOff(tid), j = 0; j < 4; j++)
    {
      OMC_EXPENSIVE_ASSERT(local_mesh.tetNode(idoff + j) !=
                             TetMesh::INFINITE_VERTEX,
                           "Infinite vertex.");
      local_mesh.tetNode(idoff + j) = vertices[local_mesh.tetNode(idoff + j)];
    }
  }
  // attach the local nodes to the global mesh
  tet_mesh.tet_node.insert(tet_mesh.tet_node.end(), local_mesh.tet_node.begin(),
                           local_mesh.tet_node.end());

  // update the vertex-tetrahedron incident relation in the global mesh
  for (index_t vid = 0; vid < n_local_verts; vid++)
  {
    if (remain_vtx[vid])
    {
      OMC_EXPENSIVE_ASSERT(local_mesh.incTet(vid) < n_local_tets,
                           "Wrong incTet.");
      tet_mesh.incTet(vertices[vid]) = local_mesh.incTet(vid) + n_global_tets;
    }
  }
  // update the neighbor in the local mesh
  size_t _4n_global_tets = n_global_tets * 4;
  for (index_t tid = 0; tid < n_local_tets; tid++)
  {
    for (index_t idoff = TetMesh::toIdOff(tid), j = 0; j < 4; j++)
    {
      index_t &neigh = local_mesh.tetNeigh(idoff + j);
      if (is_valid_idx(neigh))
        neigh += _4n_global_tets;
    }
  }
  // attach the local neighbor to the global mesh
  tet_mesh.tet_neigh.insert(tet_mesh.tet_neigh.end(),
                            local_mesh.tet_neigh.begin(),
                            local_mesh.tet_neigh.end());

  OMC_EXPENSIVE_ASSERT(tet_mesh.tet_neigh.size() == 4 * tet_mesh.sizeTets(),
                       "Wrong tet_neigh size.");
  // connect the neighbors adjacent to the cavity boundary
  for (const BCP &bcp : bcpairs)
  {
    index_t t0 = tet_reindex[TetMesh::getId(bcp.c0)];
    index_t t1 = tet_reindex[TetMesh::getId(bcp.c1)];
    OMC_EXPENSIVE_ASSERT((t0 < n_local_tets || t1 < n_local_tets) &&
                           (t0 >= n_local_tets || t1 >= n_local_tets),
                         "One inside tet and one outside tet.");
    index_t inner_corner             = t0 < n_local_tets
                                         ? TetMesh::toIdOff(t0) + TetMesh::clipOff(bcp.c0)
                                         : TetMesh::toIdOff(t1) + TetMesh::clipOff(bcp.c1);
    index_t bnd                      = bcp.bnd;
    index_t global_corner            = _4n_global_tets + inner_corner;
    tet_mesh.tetNeigh(global_corner) = bnd;
    tet_mesh.tetNeigh(bnd)           = global_corner;
  }
  // update the tet marks in the global mesh
  tet_mesh.resizeTets(tet_mesh.sizeTets());

  // (5) Collect the base faces (corners). They are the boundary faces of the
  // cavity, located on the PLC face.
  for (index_t tid = n_global_tets; tid < tet_mesh.sizeTets(); tid++)
  {
    index_t idoff = TetMesh::toIdOff(tid);
    for (index_t j = 0; j < 4; j++)
    {
      if (is_valid_idx(tet_mesh.tetNeigh(idoff + j)))
        continue;
#ifdef OMC_ENABLE_EXPENSIVE_ASSERT
      index_t v0, v1, v2;
      tet_mesh.faceVertices(idoff + j, v0, v1, v2);
      OMC_ASSERT(v_orient[v0] == Sign::ZERO && v_orient[v1] == Sign::ZERO &&
                   v_orient[v2] == Sign::ZERO,
                 "Wrong base face.");
#endif
      base.push_back(idoff + j);
    }
  }
}

/**
 * @brief Check if a tetrahedron `tet_idoff` intersects a PLC face `face`.
 */
template <typename Traits>
bool FaceRecover<Traits>::tetIntersectsFace(index_t        tet_idoff,
                                            const PLCFace &face)
{
  auto isVtxBounding = [this](index_t vid) { return v_count[vid] > 0; };

  if (!tet_mesh.isFiniteTet(tet_idoff))
    return false;

  index_t tid = face.triangles[0];
  index_t t0 = plc.triVtx(tid, 0), t1 = plc.triVtx(tid, 1),
          t2 = plc.triVtx(tid, 2);

  for (int i = 0; i < 4; i++)
    for (int j = i + 1; j < 4; j++)
    {
      index_t ei = tet_mesh.tetNode(tet_idoff + i);
      index_t ej = tet_mesh.tetNode(tet_idoff + j);
      if (isVtxBounding(ei) || isVtxBounding(ej))
        continue; // the edge is coplanar

      Sign oei = orient3dCached(t0, t1, t2, ei);
      Sign oej = orient3dCached(t0, t1, t2, ej);
      if ((oei >= Sign::ZERO && oej >= Sign::ZERO) ||
          (oei <= Sign::ZERO && oej <= Sign::ZERO))
        continue; // the inner edge does not crosses the face

      if (segCrossesFace(ei, ej, face))
        return true; // the inner edge crosses one of the triangles
    }

  return false;
}

/**
 * @brief Compute orient3d for four vertices and cache the result.
 * Use a cached value if available to avoid redundant calculations.
 * @param v_0_1_2 the first three vertices define a plane
 * @param v3 the fourth query vertex
 * @return The orientation of v3 with respect to the plane defined by v_0_1_2.
 */
template <typename Traits>
Sign FaceRecover<Traits>::orient3dCached(index_t v0, index_t v1, index_t v2,
                                         index_t v3)
{
  if (v_orient[v3] != Sign::UNCERTAIN)
    return v_orient[v3];
  v_orient[v3] = Orient3D()(point(v0), point(v1), point(v2), point(v3));
  v_cached_orient.push_back(v3);
  return v_orient[v3];
}

/**
 * @brief Checks if a segment crosses any triangle within a given PLC face.
 * @param s0 The index of the first point of the segment.
 * @param s1 The index of the second point of the segment.
 * @param face The face containing the triangles to check for intersection.
 * @return True if the segment crosses any triangle in the face.
 * @note Assume that segment crosses the support plane of the PLC face.
 */
template <typename Traits>
bool FaceRecover<Traits>::segCrossesFace(index_t s0, index_t s1,
                                         const PLCFace &face) const
{
  const GPoint3 &p0 = point(s0), &p1 = point(s1);

  for (index_t tid : face.triangles)
  {
    const GPoint3 &v0 = point(plc.triVtx(tid, 0)),
                  &v1 = point(plc.triVtx(tid, 1)),
                  &v2 = point(plc.triVtx(tid, 2));
    // We have known that the segment crosses the support plane,
    // so, we skip checking orientation of p0, p1 w.r.t the plane.

    // We check whether the segment crosses the triangle.
    // s crosses t (borders included), if the signs of the three tetrahedra
    // obtained combining s with the three edges of t are all equal.
    Sign o1 = Orient3D()(p0, p1, v0, v1);
    Sign o2 = Orient3D()(p0, p1, v1, v2);
    if ((o1 > Sign::ZERO && o2 < Sign::ZERO) ||
        (o1 < Sign::ZERO && o2 > Sign::ZERO))
      continue;
    Sign o3 = Orient3D()(p0, p1, v2, v0);
    if ((o1 > Sign::ZERO && o3 < Sign::ZERO) ||
        (o1 < Sign::ZERO && o3 > Sign::ZERO))
      continue;
    if ((o2 > Sign::ZERO && o3 < Sign::ZERO) ||
        (o2 < Sign::ZERO && o3 > Sign::ZERO))
      continue;
    // OPT orientation can be cached locally

    return true;
  }
  return false;
}

} // namespace OMC
