#pragma once

#include "DetectClassifyTTI.h"

namespace OMC {

template <typename Traits>
DetectClassifyTTI<Traits>::DetectClassifyTTI(TriSoup &_ts, AuxStruct &_uniq_g,
                                             AuxStruct &_g,
                                             PntArena  &_pnt_arena,
                                             IdxArena  &_idx_arena)
  : ts(_ts)
  , uniq_g(_uniq_g)
  , g(_g)
  , pnt_arena(_pnt_arena)
  , idx_arena(_idx_arena)
{
}

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI(index_t ta, index_t tb)
{
	TTIHelper ha, hb;

	// vertices from two triangles `ta` and `tb`.
	ha.t_id = ta;
	ha.v    = {ts.triVertPtr(ta, 0), ts.triVertPtr(ta, 1), ts.triVertPtr(ta, 2)};
	ha.v_id = {ts.triVertID(ta, 0), ts.triVertID(ta, 1), ts.triVertID(ta, 2)};
	ha.e_id = {ts.triEdgeID(ta, 0), ts.triEdgeID(ta, 1), ts.triEdgeID(ta, 2)};

	hb.t_id = tb;
	hb.v    = {ts.triVertPtr(tb, 0), ts.triVertPtr(tb, 1), ts.triVertPtr(tb, 2)};
	hb.v_id = {ts.triVertID(tb, 0), ts.triVertID(tb, 1), ts.triVertID(tb, 2)};
	hb.e_id = {ts.triEdgeID(tb, 0), ts.triEdgeID(tb, 1), ts.triEdgeID(tb, 2)};

	// TODO get e_id only if necessary.

	// check if any triangle is degenerate
	OMC_EXPENSIVE_ASSERT(
	  CollinearPoints3D().misaligned(ha.v[0], ha.v[1], ha.v[2]) &&
	    CollinearPoints3D().misaligned(hb.v[0], hb.v[1], hb.v[2]),
	  "Detect degenerate triangle in check TTI.");

	// find vert correspondences
	// (we have removed duplicate vertices in arrangements pipeline, so, we can
	// check correspondences by index instead of geometry.)
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 3; j++)
			if (ha.v_id[i] == hb.v_id[j])
			{
				ha.v_in_vtx[i] = j;
				hb.v_in_vtx[j] = i;
			}

	// count number of coincident vertices in `ta` and `tb`
	size_t ta_count = (size_t)is_valid_idx(ha.v_in_vtx[0]) +
	                  (size_t)is_valid_idx(ha.v_in_vtx[1]) +
	                  (size_t)is_valid_idx(ha.v_in_vtx[2]);

	// `ta` and `tb` are coincident and do not intersect.
	// (although such case is impossible in arrangements pipeline because we have
	// removed duplicate triangles, we still keep this check for completeness.)
	OMC_EXPENSIVE_ASSERT(ta_count != 3, "Duplicate triangles are detected.");
	if (ta_count == 3)
		return;

	if (ta_count == 2) // `ta` and `tb` share an edge.
		check_TTI_share_edge(ha, hb);
	else if (ta_count == 1) // `ta` and `tb` share a vertex.
		check_TTI_share_vertex(ha, hb);
	else // `ta` and `tb` are separate
		check_TTI_separate(ha, hb);
}

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI_share_edge(TTIHelper &ha,
                                                     TTIHelper &hb)
{
	// Let e be the shared edge and { oppa, oppb } be the two vertices opposite to
	// e in `ta` and `tb`, respectively. If oppa and oppb lie at the same side of
	// e, the two triangles overlap. Otherwise they are edge-adjacent and form a
	// valid simplicial complex.
	index_t ea   = 0; // local indices of the shared edge in `ta`
	index_t eb   = 0; // local indices of the shared edge in `tb`
	index_t oppa = 0; // local index of the vertex opposite to e in `ta`
	index_t oppb = 0; // local index of the vertex opposite to e in `tb`

	for (size_t i = 0; i < 3; ++i)
	{
		if (!is_valid_idx(ha.v_in_vtx[i]))
		{
			oppa = i;
			ea   = (i + 1) % 3;
		}
		if (!is_valid_idx(hb.v_in_vtx[i]))
		{
			oppb = i;
			eb   = (i + 1) % 3;
		}
	}
	// check if `ta` and `tb` intersect.
	if (Orient3D()(ha.v[0], ha.v[1], ha.v[2], hb.v[oppb]) != Sign::ZERO)
		return; // do not intersect

	// now `ta` and `tb` are coplanar
	ha.t_nmax = MaxCompInTriNormal()(ha.v[0], ha.v[1], ha.v[2]);
	hb.t_nmax = ha.t_nmax;

	// check orientation of vertices oppa and oppb with respect to edge ea.
	Sign oppa_wrt_ea =
	  OrientOn2D()(ha.v[ea], ha.v[(ea + 1) % 3], ha.v[oppa], ha.t_nmax);
	Sign oppb_wrt_ea =
	  OrientOn2D()(ha.v[ea], ha.v[(ea + 1) % 3], hb.v[oppb], ha.t_nmax);

	OMC_EXPENSIVE_ASSERT(oppa_wrt_ea != Sign::ZERO && oppb_wrt_ea != Sign::ZERO,
	                     "degenerate triangle detected.");
	if (oppa_wrt_ea != oppb_wrt_ea)
		return; // do not intersect

	// cache the orientation result
	hb.v_wrt_seg[oppb][ea] = oppb_wrt_ea;
	if (ha.v_id[ea] == hb.v_id[eb])
	{
		OMC_EXPENSIVE_ASSERT(oppa_wrt_ea == get_v_wrt_seg(ha, oppa, hb, eb),
		                     "correspond orientation is wrong.");
		ha.v_wrt_seg[oppa][eb] = oppa_wrt_ea;
	}
	else
	{
		OMC_EXPENSIVE_ASSERT(oppa_wrt_ea ==
		                       reverse_sign(get_v_wrt_seg(ha, oppa, hb, eb)),
		                     "correspond orientation is wrong.");
		ha.v_wrt_seg[oppa][eb] = reverse_sign(oppa_wrt_ea);
	}

	// if they intersect, each triangle has one sharing edge and other two
	// coplanar edges.
	g.addCoplanarTriangles(ha.t_id, hb.t_id);
	uniq_g.setTriangleHasIntersections(ha.t_id);
	uniq_g.setTriangleHasIntersections(hb.t_id);

	// intersection list
	phmap::flat_hash_set<index_t>             inter_list;
	// only coplanar intersection points (edge-edge)
	phmap::flat_hash_set<CoplanarEEI, Hasher> coplanar_edge_crosses;

	// four checks for intersections between a coplanar edge and a triangle.
	// edges from `ta` (except the shared edge) with respect to `tb`
	classify_coplanar_edge_intersections(ha, (ea + 1) % 3, hb, inter_list,
	                                     coplanar_edge_crosses);
	classify_coplanar_edge_intersections(ha, (ea + 2) % 3, hb, inter_list,
	                                     coplanar_edge_crosses);
	// edges from `tb` (except the shared edge) with respect to `ta`
	classify_coplanar_edge_intersections(hb, (eb + 1) % 3, ha, inter_list,
	                                     coplanar_edge_crosses);
	classify_coplanar_edge_intersections(hb, (eb + 2) % 3, ha, inter_list,
	                                     coplanar_edge_crosses);

	OMC_EXPENSIVE_ASSERT(
	  inter_list.size() <= 1,
	  "more than 1 intersection points in coplanar triangles.");

	return; // end check.
}

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI_share_vertex(TTIHelper &ha,
                                                       TTIHelper &hb)
{
	// check if `ta` and `tb` intersect.

	// intersection list
	phmap::flat_hash_set<index_t>             inter_list;
	// only coplanar intersection points (edge-edge)
	phmap::flat_hash_set<CoplanarEEI, Hasher> coplanar_edge_crosses;

	index_t va = 0; // index of the shared vertex in `ta`
	index_t vb = 0; // index of the shared vertex in `tb`
	for (size_t i = 0; i < 3; ++i)
	{
		if (is_valid_idx(ha.v_in_vtx[i]))
			va = i;
		if (is_valid_idx(hb.v_in_vtx[i]))
			vb = i;
	}
	index_t ea = (va + 1) % 3; // index of the opposite edge in `ta`
	index_t eb = (vb + 1) % 3; // index of the opposite edge in `tb`

	index_t eav0 = ea, eav1 = (ea + 1) % 3; // two endpoints of ea
	index_t ebv0 = eb, ebv1 = (eb + 1) % 3; // two endpoints of eb

	// OPT most expensive
	Sign eav0_wrt_tb = Orient3D()(hb.v[0], hb.v[1], hb.v[2], ha.v[eav0]);
	Sign eav1_wrt_tb = Orient3D()(hb.v[0], hb.v[1], hb.v[2], ha.v[eav1]);

	if (eav0_wrt_tb != Sign::ZERO && eav0_wrt_tb == eav1_wrt_tb)
		return; // above or below `tb`, do not intersect

	if (eav0_wrt_tb == Sign::ZERO && eav1_wrt_tb == Sign::ZERO)
	{
		// ea is coplanar to `tb`, so `ta` and `tb` are coplanar
		auto &il  = inter_list;            // a short name :)
		auto &cec = coplanar_edge_crosses; // a short name :)

		ha.t_nmax = MaxCompInTriNormal()(ha.v[0], ha.v[1], ha.v[2]);
		hb.t_nmax = ha.t_nmax;

		size_t eav0_wrt_sector = get_vtx_wrt_sector(ha, eav0, hb, eb);
		if (eav0_wrt_sector == 0) // previous edge of ea intersects tb
			classify_coplanar_edge_intersections(ha, (ea + 2) % 3, hb, il, cec);
		size_t eav1_wrt_sector = get_vtx_wrt_sector(ha, eav1, hb, eb);
		if (eav1_wrt_sector == 0) // next edge of ea intersects tb
			classify_coplanar_edge_intersections(ha, (ea + 1) % 3, hb, il, cec);

		if ((eav0_wrt_sector == 1 || eav0_wrt_sector == 3) &&
		    (eav1_wrt_sector == 1 || eav1_wrt_sector == 3))
			return; // ta and tb do not intersect
		if ((eav0_wrt_sector == 2 || eav0_wrt_sector == 3) &&
		    (eav1_wrt_sector == 2 || eav1_wrt_sector == 3))
			return; // ta and tb do not intersect

		size_t ebv0_wrt_sector = get_vtx_wrt_sector(hb, ebv0, ha, ea);
		if (ebv0_wrt_sector == 0) // previous edge of eb intersects ta
			classify_coplanar_edge_intersections(hb, (eb + 2) % 3, ha, il, cec);
		size_t ebv1_wrt_sector = get_vtx_wrt_sector(hb, ebv1, ha, ea);
		if (ebv1_wrt_sector == 0) // next edge of eb intersects ta
			classify_coplanar_edge_intersections(hb, (eb + 1) % 3, ha, il, cec);

		if ((ebv0_wrt_sector == 1 || ebv0_wrt_sector == 3) &&
		    (ebv1_wrt_sector == 1 || ebv1_wrt_sector == 3))
			return; // ta and tb do not intersect
		if ((ebv0_wrt_sector == 2 || ebv0_wrt_sector == 3) &&
		    (ebv1_wrt_sector == 2 || ebv1_wrt_sector == 3))
			return; // ta and tb do not intersect

		classify_coplanar_edge_intersections(ha, ea, hb, il, cec);
		classify_coplanar_edge_intersections(hb, eb, ha, il, cec);

		if (!inter_list.empty())
		{
			g.addCoplanarTriangles(ha.t_id, hb.t_id);
			uniq_g.setTriangleHasIntersections(ha.t_id);
			uniq_g.setTriangleHasIntersections(hb.t_id);
		}

		OMC_EXPENSIVE_ASSERT(
		  inter_list.size() <= 3,
		  "more than 3 intersection points between coplanar triangles.");
		return;
	}

	// OPT most expensive
	Sign ebv0_wrt_ta = Orient3D()(ha.v[0], ha.v[1], ha.v[2], hb.v[ebv0]);
	Sign ebv1_wrt_ta = Orient3D()(ha.v[0], ha.v[1], ha.v[2], hb.v[ebv1]);

	if (ebv0_wrt_ta != Sign::ZERO && ebv0_wrt_ta == ebv1_wrt_ta)
		return; // above or below `tb`, do not intersect

	// go on classifying intersections on `tb`
	ha.t_nmax = MaxCompInTriNormal()(ha.v[0], ha.v[1], ha.v[2]);
	hb.t_nmax = MaxCompInTriNormal()(hb.v[0], hb.v[1], hb.v[2]);

	if (eav0_wrt_tb == Sign::POSITIVE && eav1_wrt_tb == Sign::NEGATIVE ||
	    eav0_wrt_tb == Sign::NEGATIVE && eav1_wrt_tb == Sign::POSITIVE)
	{ // ea cross the support plane of `tb`

		// the shared vertex is seen as an intersection point.
		inter_list.insert(ha.v_id[va]);

		// check if a new intersection point exists between ea and `tb`.
		if (noncoplanar_seg_tri_do_intersect(ha, ea, hb))
			classify_noncoplanar_edge_intersections(ha, ea, hb, inter_list);

		OMC_EXPENSIVE_ASSERT(
		  inter_list.size() <= 2,
		  "more than 2 intersection points between non-coplanar triangles.");
		// If we have detected more than one intersections in non-coplanar cases,
		// we end classification now and propagate undetected intersections later.
		if (inter_list.size() == 2)
		{
			// before ending, we record possible coplanar edges.
			if (ebv0_wrt_ta == Sign::ZERO || ebv1_wrt_ta == Sign::ZERO)
			{ // `tb` has a single coplanar edge to `ta`
				index_t copl_eb =
				  ebv0_wrt_ta == Sign::ZERO ? (eb + 2) % 3 : (eb + 1) % 3;
				index_t copl_vb = ebv0_wrt_ta == Sign::ZERO ? eb : (eb + 1) % 3;

				if (get_vtx_wrt_sector(hb, copl_vb, ha, ea) == 0)
					g.addCoplanarEdge(ha.t_id, hb.e_id[copl_eb]);
			}

			index_t v0_id = *(inter_list.begin());
			index_t v1_id = *(++inter_list.begin());
			add_symbolic_segment(v0_id, v1_id, ha.t_id, hb.t_id);
			uniq_g.setTriangleHasIntersections(ha.t_id);
			uniq_g.setTriangleHasIntersections(hb.t_id);

			return; // all possible intersections are found, return.
		}
	}
	else if (eav0_wrt_tb == Sign::ZERO || eav1_wrt_tb == Sign::ZERO)
	{ // one edge of `ta` on the support plane of `tb`
		index_t copl_ea = eav0_wrt_tb == Sign::ZERO ? (ea + 2) % 3 : (ea + 1) % 3;
		index_t copl_va = eav0_wrt_tb == Sign::ZERO ? ea : (ea + 1) % 3;

		if (get_vtx_wrt_sector(ha, copl_va, hb, eb) == 0)
		{
			classify_coplanar_edge_intersections(ha, copl_ea, hb, inter_list,
			                                     coplanar_edge_crosses);
			if (!inter_list.empty())
			{
				g.addCoplanarEdge(hb.t_id, ha.e_id[copl_ea]);
				uniq_g.setTriangleHasIntersections(ha.t_id);
				uniq_g.setTriangleHasIntersections(hb.t_id);
			}
		}
		return; // whether intersections are found or not, return.
	}
	else
	{
		OMC_ASSERT(false, "Impossible case happens.");
	}

	// go on classifying intersections on `ta`
	if (ebv0_wrt_ta == Sign::POSITIVE && ebv1_wrt_ta == Sign::NEGATIVE ||
	    ebv0_wrt_ta == Sign::NEGATIVE && ebv1_wrt_ta == Sign::POSITIVE)
	{ // `eb` cross the support plane of `ta`

		// the shared vertex is seen as an intersection point.
		inter_list.insert(hb.v_id[vb]);

		// check if a new intersection point exists between ea and `tb`.
		if (noncoplanar_seg_tri_do_intersect(hb, eb, ha))
			classify_noncoplanar_edge_intersections(hb, eb, ha, inter_list);

		OMC_EXPENSIVE_ASSERT(
		  inter_list.size() <= 2,
		  "more than 2 intersection points between non-coplanar triangles.");
		// If we have detected more than one intersections in non-coplanar cases,
		// we end classification now and propagate undetected intersections later.
		if (inter_list.size() == 2)
		{
			index_t v0_id = *(inter_list.begin());
			index_t v1_id = *(++inter_list.begin());
			add_symbolic_segment(v0_id, v1_id, ha.t_id, hb.t_id);
			uniq_g.setTriangleHasIntersections(ha.t_id);
			uniq_g.setTriangleHasIntersections(hb.t_id);
		}
	}
	else if (ebv0_wrt_ta == Sign::ZERO || ebv1_wrt_ta == Sign::ZERO)
	{ // one edge of `tb` on the support plane of `ta`
		index_t copl_eb = ebv0_wrt_ta == Sign::ZERO ? (eb + 2) % 3 : (eb + 1) % 3;
		index_t copl_vb = ebv0_wrt_ta == Sign::ZERO ? eb : (eb + 1) % 3;

		if (get_vtx_wrt_sector(hb, copl_vb, ha, ea) == 0)
		{
			classify_coplanar_edge_intersections(hb, copl_eb, ha, inter_list,
			                                     coplanar_edge_crosses);
			if (!inter_list.empty())
			{
				g.addCoplanarEdge(ha.t_id, hb.e_id[copl_eb]);
				uniq_g.setTriangleHasIntersections(ha.t_id);
				uniq_g.setTriangleHasIntersections(hb.t_id);
			}
		}
	}
	else
	{
		OMC_ASSERT(false, "Impossible case happens.");
	}
}

// clang-format off
inline bool _sameOrientation(const Sign o1, const Sign o2) { return o1 == o2; }
// true if all edges are coplanar to the triangle, return false otherwise
inline bool _allCoplanarEdges(const Sign o[])
{
	return (o[0] == Sign::ZERO && o[1] == Sign::ZERO && o[2] == Sign::ZERO);
}
// if there is a coplanar edge return true and set edge_id, return false
// otherwise
inline bool _singleCoplanarEdge(const Sign o[], index_t &edge_id)
{
	if (o[0] == Sign::ZERO && o[1] == Sign::ZERO && o[2] != Sign::ZERO) { edge_id = 0; return true; }
	if (o[1] == Sign::ZERO && o[2] == Sign::ZERO && o[0] != Sign::ZERO) { edge_id = 1; return true; }
	if (o[2] == Sign::ZERO && o[0] == Sign::ZERO && o[1] != Sign::ZERO) { edge_id = 2; return true; }
	edge_id = InvalidIndex;
	return false; // false
}
// if there is a vertex in the plane and the opposite edge doesn't intersect the
// plane return true and set the vtx id, return false otherwise
inline bool _vtxInPlaneAndOppositeEdgeOnSameSide(const Sign o[], index_t &vtx_id)
{
	if (o[0] == Sign::ZERO && o[1] == o[2] && o[1] != Sign::ZERO) { vtx_id = 0; return true; }
	if (o[1] == Sign::ZERO && o[0] == o[2] && o[0] != Sign::ZERO) { vtx_id = 1; return true; }
	if (o[2] == Sign::ZERO && o[0] == o[1] && o[0] != Sign::ZERO) { vtx_id = 2; return true; }
	vtx_id = InvalidIndex;
	return false;
}
// if there is a vertex in the plane and the opposite edge intersect the plane
// return true and set the vtx id, return false otherwise
inline bool _vtxInPlaneAndOppositeEdgeCrossPlane(const Sign o[], index_t &vtx_id)
{
	if (o[0] == Sign::ZERO && o[1] != o[2] && o[1] != Sign::ZERO && o[2] != Sign::ZERO) { vtx_id = 0; return true; }
	if (o[1] == Sign::ZERO && o[0] != o[2] && o[0] != Sign::ZERO && o[2] != Sign::ZERO) { vtx_id = 1; return true; }
	if (o[2] == Sign::ZERO && o[0] != o[1] && o[0] != Sign::ZERO && o[1] != Sign::ZERO) { vtx_id = 2; return true; }
	vtx_id = InvalidIndex;
	return false;
}
// if there is a vertex on one side and the opposite edge on the other return
// the relative informations, -1 otherwise
inline bool _vtxOnASideAndOppositeEdgeOnTheOther(const Sign o[], index_t &vtx_id)
{
	// one vtx on the plane
	if (o[0] == Sign::ZERO || o[1] == Sign::ZERO || o[2] == Sign::ZERO) { vtx_id = InvalidIndex; return false; }
	// all vtx on the same side of the plane
	if (o[0] == o[1] && o[1] == o[2]) { vtx_id = InvalidIndex; return false; }
	if (o[0] == o[1]) { vtx_id = 2; return true; }
	if (o[0] == o[2]) { vtx_id = 1; return true; }
	vtx_id = 0;
	return true;
}
// clang-format on

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI_separate(TTIHelper &ha, TTIHelper &hb)
{
	// `ta` and `tb` do not share any edge or vertex.

	// intersection list
	phmap::flat_hash_set<index_t>             inter_list;
	// only coplanar intersection points (edge-edge)
	phmap::flat_hash_set<CoplanarEEI, Hasher> coplanar_edge_crosses;

	// =========================================================================
	// check of A respect to B
	// =========================================================================

	index_t edge_id, vtx_id;

	Sign orAB[3]; // orientation of edge of `ta` with respect to `tb`
	orAB[0] = get_v_wrt_tri(ha, 0, hb);
	orAB[1] = get_v_wrt_tri(ha, 1, hb);
	orAB[2] = get_v_wrt_tri(ha, 2, hb);

	if (_sameOrientation(orAB[0], orAB[1]) &&
	    _sameOrientation(orAB[1], orAB[2]) && (orAB[0] != Sign::ZERO))
	{
		// CASE: no intersection found
		return;
	}

	if (_allCoplanarEdges(orAB))
	{
		// CASE: all edge of ta are coplanar to all edges of tb   (orAB: 0 0 0)
		hb.t_nmax = MaxCompInTriNormal()(hb.v[0], hb.v[1], hb.v[2]);
		ha.t_nmax = hb.t_nmax;

		if (coplanar_seg_tri_do_intersect(ha, 0, hb))
			classify_coplanar_edge_intersections(ha, 0, hb, inter_list,
			                                     coplanar_edge_crosses);
		if (coplanar_seg_tri_do_intersect(ha, 1, hb))
			classify_coplanar_edge_intersections(ha, 1, hb, inter_list,
			                                     coplanar_edge_crosses);
		if (coplanar_seg_tri_do_intersect(ha, 2, hb))
			classify_coplanar_edge_intersections(ha, 2, hb, inter_list,
			                                     coplanar_edge_crosses);

		if (coplanar_seg_tri_do_intersect(hb, 0, ha))
			classify_coplanar_edge_intersections(hb, 0, ha, inter_list,
			                                     coplanar_edge_crosses);
		if (coplanar_seg_tri_do_intersect(hb, 1, ha))
			classify_coplanar_edge_intersections(hb, 1, ha, inter_list,
			                                     coplanar_edge_crosses);
		if (coplanar_seg_tri_do_intersect(hb, 2, ha))
			classify_coplanar_edge_intersections(hb, 2, ha, inter_list,
			                                     coplanar_edge_crosses);

		if (!inter_list.empty())
		{
			g.addCoplanarTriangles(ha.t_id, hb.t_id);
			uniq_g.setTriangleHasIntersections(ha.t_id);
			uniq_g.setTriangleHasIntersections(hb.t_id);
		}

		return;
	}

	hb.t_nmax = MaxCompInTriNormal()(hb.v[0], hb.v[1], hb.v[2]);
	ha.t_nmax = MaxCompInTriNormal()(ha.v[0], ha.v[1], ha.v[2]);

	if (_singleCoplanarEdge(orAB, edge_id))
	{
		// CASE: a single edge of ta is coplanar to tb    (e.g. orAB: 1 0 0)
		if (coplanar_seg_tri_do_intersect(ha, edge_id, hb))
		{
			classify_coplanar_edge_intersections(ha, edge_id, hb, inter_list,
			                                     coplanar_edge_crosses);
			if (!inter_list.empty())
			{
				g.addCoplanarEdge(hb.t_id, ha.e_id[edge_id]);
				uniq_g.setTriangleHasIntersections(ha.t_id);
				uniq_g.setTriangleHasIntersections(hb.t_id);
			}
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeOnSameSide(orAB, vtx_id))
	{
		// CASE: a vertex of ta is coplanar to tb, and the opposite edge is on the
		// same side respect to tb  (e.g. orAB: 1 0 1)
		if (classify_coplanr_vtx_intersections(ha, vtx_id, hb, inter_list))
		{
			uniq_g.setTriangleHasIntersections(ha.t_id);
			uniq_g.setTriangleHasIntersections(hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeCrossPlane(orAB, vtx_id))
	{
		// CASE: a vertex of ta is coplanar to tb, and the opposite edge could
		// intersect tb (e.g. orAB: -1 0 1)
		classify_coplanr_vtx_intersections(ha, vtx_id, hb, inter_list);

		index_t opp_edge_id = (vtx_id + 1) % 3;

		if (noncoplanar_seg_tri_do_intersect(ha, opp_edge_id, hb))
			classify_noncoplanar_edge_intersections(ha, opp_edge_id, hb, inter_list);
		// go on checking B->A
	}
	else if (_vtxOnASideAndOppositeEdgeOnTheOther(orAB, vtx_id))
	{
		// CASE: a vertex of ta is on one side of the plane defined to tb, and the
		// opposite edge (always in ta) is in the other (e.g. orAB: -1 1 1)
		if (noncoplanar_seg_tri_do_intersect(ha, vtx_id, hb))
			classify_noncoplanar_edge_intersections(ha, vtx_id, hb, inter_list);
		if (noncoplanar_seg_tri_do_intersect(ha, (vtx_id + 2) % 3, hb))
			classify_noncoplanar_edge_intersections(ha, (vtx_id + 2) % 3, hb,
			                                        inter_list);
		// go on checking B->A
	}

	// =========================================================================
	// check of B respect to A
	// =========================================================================

	Sign orBA[3];

	orBA[0] = get_v_wrt_tri(hb, 0, ha);
	orBA[1] = get_v_wrt_tri(hb, 1, ha);
	orBA[2] = get_v_wrt_tri(hb, 2, ha);

	OMC_EXPENSIVE_ASSERT(
	  inter_list.size() <= 2,
	  "more than 2 intersection points between non-coplanar triangles.");
	if (inter_list.size() == 2)
	{
		// when detect more than one intersections in noncoplanar cases,
		// we do not continue classifying but propagate undetected intersections in
		// later stages.

		// Check if an edge from tA is coplanar to tB and intersects tB.
		if (_singleCoplanarEdge(orBA, edge_id)) // TODO remove orBA
		{
			if (coplanar_seg_tri_do_intersect(hb, edge_id, ha))
				g.addCoplanarEdge(ha.t_id, hb.e_id[edge_id]);
		}

		index_t v0_id = *(inter_list.begin());
		index_t v1_id = *(++inter_list.begin());

		add_symbolic_segment(v0_id, v1_id, ha.t_id, hb.t_id);
		uniq_g.setTriangleHasIntersections(ha.t_id);
		uniq_g.setTriangleHasIntersections(hb.t_id);

		return; // all possible intersections are found, return.
	}

	if (_sameOrientation(orBA[0], orBA[1]) &&
	    _sameOrientation(orBA[1], orBA[2]) && (orBA[0] != Sign::ZERO))
	{
		// CASE: no intersection found
		return;
	}

	if (_singleCoplanarEdge(orBA, edge_id))
	{
		// CASE: a single edge of tB is coplanar to tA    (e.g. orBA: 1 0 0)
		if (coplanar_seg_tri_do_intersect(hb, edge_id, ha))
		{
			classify_coplanar_edge_intersections(hb, edge_id, ha, inter_list,
			                                     coplanar_edge_crosses);
			if (!inter_list.empty())
			{
				g.addCoplanarEdge(ha.t_id, hb.e_id[edge_id]);
				uniq_g.setTriangleHasIntersections(ha.t_id);
				uniq_g.setTriangleHasIntersections(hb.t_id);
			}
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeOnSameSide(orBA, vtx_id))
	{
		// CASE: a vertex of tB is coplanar to tA, and the opposite edge is on the
		// same side respect to tA  (e.g. orBA: 1 0 1)
		if (classify_coplanr_vtx_intersections(hb, vtx_id, ha, inter_list))
		{
			uniq_g.setTriangleHasIntersections(ha.t_id);
			uniq_g.setTriangleHasIntersections(hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeCrossPlane(orBA, vtx_id))
	{
		// CASE: a vertex of tB is coplanar to tA, and the opposite edge could
		// intersect tA (e.g. orBA: -1 0 1)
		classify_coplanr_vtx_intersections(hb, vtx_id, ha, inter_list);

		index_t opp_edge_id = (vtx_id + 1) % 3;

		if (noncoplanar_seg_tri_do_intersect(hb, opp_edge_id, ha))
			classify_noncoplanar_edge_intersections(hb, opp_edge_id, ha, inter_list);
	}
	else if (_vtxOnASideAndOppositeEdgeOnTheOther(orBA, vtx_id))
	{
		// CASE: a vertex of tB is on one side of the plane defined to tA, and the
		// opposite edge (always in tB) is in the other (e.g. orBA: -1 1 1)
		if (noncoplanar_seg_tri_do_intersect(hb, vtx_id, ha))
			classify_noncoplanar_edge_intersections(hb, vtx_id, ha, inter_list);
		if (noncoplanar_seg_tri_do_intersect(hb, (vtx_id + 2) % 3, ha))
			classify_noncoplanar_edge_intersections(hb, (vtx_id + 2) % 3, ha,
			                                        inter_list);
	}

	OMC_EXPENSIVE_ASSERT(
	  inter_list.size() <= 2,
	  "more than 2 intersection points between noncoplanar traingles");

	if (inter_list.size() == 2)
	{
		index_t v0_id = *(inter_list.begin());
		index_t v1_id = *(++inter_list.begin());

		add_symbolic_segment(v0_id, v1_id, ha.t_id, hb.t_id);
		uniq_g.setTriangleHasIntersections(ha.t_id);
		uniq_g.setTriangleHasIntersections(hb.t_id);
	}
}

/**
 * @param ha helper for `ta`
 * @param va A valid local vertex index in `ta`
 * @param hb helper for `tb`
 * @param eb A valid local edge index in `tb`
 * @return true if `va` is strictly inside `eb`
 * @note we assume that `va` is on the support line of `eb`
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::get_v_in_seg(TTIHelper &ha, index_t va,
                                             TTIHelper &hb, index_t eb)
{
	const NT *p  = ha.v[va];
	const NT *s0 = hb.v[eb], *s1 = hb.v[(eb + 1) % 3];
	return ((p[0] > std::min(s0[0], s1[0]) && p[0] < std::max(s0[0], s1[0])) ||
	        (p[1] > std::min(s0[1], s1[1]) && p[1] < std::max(s0[1], s1[1])) ||
	        (p[2] > std::min(s0[2], s1[2]) && p[2] < std::max(s0[2], s1[2])));
}

/**
 * @param ha helper for `ta`.
 * @param va A valid local vertex index in `ta`.
 * @param hb helper for `tb`.
 * @return index_t the index of segment/edge of `tb` where `va` is in.
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::get_v_in_seg(TTIHelper &ha, index_t va,
                                                TTIHelper &hb)
{
	// check if cached
	if (TTIHelper::is_cached(ha.v_in_seg[va]))
		return ha.v_in_seg[va];

	// check if `va` is in any vertex of `tb`
	if (is_valid_idx(ha.v_in_vtx[va]))
	{
		ha.v_in_seg[va] = InvalidIndex;
		return InvalidIndex;
	}

	// find the edge
	for (index_t eb = 0; eb < 3; eb++)
	{ // traverse edges in `tb`
		if (get_v_wrt_seg(ha, va, hb, eb) == Sign::ZERO)
		{ // `va` is on the support line of edge ei
			const NT *p  = ha.v[va];
			const NT *s0 = hb.v[eb], *s1 = hb.v[(eb + 1) % 3];
			if ((p[0] > std::min(s0[0], s1[0]) && p[0] < std::max(s0[0], s1[0])) ||
			    (p[1] > std::min(s0[1], s1[1]) && p[1] < std::max(s0[1], s1[1])) ||
			    (p[2] > std::min(s0[2], s1[2]) && p[2] < std::max(s0[2], s1[2])))
			{
				// find the edge where `va` is in.
				ha.v_in_seg[va] = eb;
				return eb;
			}
		}
	}

	// can't find an edge where `va` is in.
	ha.v_in_seg[va] = InvalidIndex;
	return InvalidIndex;
}

/**
 * @param ha helper for `ta`.
 * @param va A valid local vertex index in `ta`.
 * @param hb helper for `tb`.
 * @return bool true if `va` strictly inside `tb`.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::get_v_in_tri(TTIHelper &ha, index_t va,
                                             TTIHelper &hb)
{
	// check if cached
	if (TTIHelper::is_cached(ha.v_in_tri[va]))
		return ha.v_in_tri[va];

	// check if `va` is in any vertex or segment of `tb`
	if (is_valid_idx(ha.v_in_vtx[va]) || is_valid_idx(get_v_in_seg(ha, va, hb)))
	{
		ha.v_in_tri[va] = false;
		return false;
	}

	Sign o0 = get_v_wrt_seg(ha, va, hb, 0);
	Sign o1 = get_v_wrt_seg(ha, va, hb, 1);
	Sign o2 = get_v_wrt_seg(ha, va, hb, 2);

	ha.v_in_tri[va] = ((o0 > Sign::ZERO && o1 > Sign::ZERO && o2 > Sign::ZERO) ||
	                   (o0 < Sign::ZERO && o1 < Sign::ZERO && o2 < Sign::ZERO));
	return ha.v_in_tri[va];
}

/**
 * @param ha helper for `ta`.
 * @param va A valid local vertex index in `ta`.
 * @param hb helper for `tb`.
 * @return Sign the orientation.
 */
template <typename Traits>
Sign DetectClassifyTTI<Traits>::get_v_wrt_tri(TTIHelper &ha, index_t va,
                                              TTIHelper &hb)
{
	// check if cached
	if (is_sign_reliable(ha.v_wrt_tri[va]))
		return ha.v_wrt_tri[va];
	// otherwise calculate it

	const NT *p  = ha.v[va];
	const NT *t0 = hb.v[0], *t1 = hb.v[1], *t2 = hb.v[2];

	Sign sign = Orient3D()(t0, t1, t2, p);

	ha.v_wrt_tri[va] = sign;
	return sign;
}

/**
 * @param ha helper for `ta`.
 * @param va A valid local vertex index in `ta`.
 * @param hb helper for `tb`.
 * @param eb A valid local edge index in `tb`.
 * @return Sign the orientation.
 */
template <typename Traits>
Sign DetectClassifyTTI<Traits>::get_v_wrt_seg(TTIHelper &ha, index_t va,
                                              TTIHelper &hb, index_t eb)
{
	// check if cached
	if (is_sign_reliable(ha.v_wrt_seg[va][eb]))
		return ha.v_wrt_seg[va][eb];
	// otherwise calculate it

	OMC_EXPENSIVE_ASSERT(hb.t_nmax != -1, "orthogonal plane not initialized.");

	const NT *p  = ha.v[va];
	const NT *s0 = hb.v[eb], *s1 = hb.v[(eb + 1) % 3];

	Sign sign = OrientOn2D()(s0, s1, p, hb.t_nmax);

	ha.v_wrt_seg[va][eb] = sign;
	return sign;
}

/**
 * @param ha helper for `ta`.
 * @param ea A valid local edge index in `ta`.
 * @param hb helper for `tb`.
 * @param eb A valid local edge index in `tb`.
 * @return Sign the orientation.
 */
template <typename Traits>
Sign DetectClassifyTTI<Traits>::get_seg_wrt_seg(TTIHelper &ha, index_t ea,
                                                TTIHelper &hb, index_t eb)
{
	// check if cached
	if (is_sign_reliable(ha.seg_wrt_seg[ea][eb]))
		return ha.seg_wrt_seg[ea][eb];
	// otherwise calculate it

	const NT *s0 = ha.v[ea], *s1 = ha.v[(ea + 1) % 3];
	const NT *t0 = hb.v[eb], *t1 = hb.v[(eb + 1) % 3];

	Sign sign = Orient3D()(s0, s1, t0, t1);

	ha.seg_wrt_seg[ea][eb] = sign; // seg_wrt_seg
	hb.seg_wrt_seg[eb][ea] = sign; // are symetric
	return sign;
}

/**
 * @brief
 * In case of that vertex `va` of triangle `ta` is on the support plane
 * of triangle `tb`, this func test the position of `va` with respect to a
 * sector formed by triangle `tb`.
 *
 * The sector is formed by the two adjacent edges of vertex `vb`, which is the
 * opposite vertex to the edge `eb` in triangle `tb`.
 * @details
 *    ____eb____
 *   	\        /
 *     \  0   /
 *  ebn \    / ebp
 *       \  /
 *  1     * vb    2
 *       / \
 *      /   \
 *     /  3  \
 *    /       \
 * @param ha helper for `ta`.
 * @param va A valid local vertex index in `ta`
 * @param hb helper for `tb`.
 * @param eb A valid local edge index in `tb`.
 * @return 0, 1, 2 and 3 represent different regions in the support plane.
 * see details for meanings of returned numbers.
 */
template <typename Traits>
size_t DetectClassifyTTI<Traits>::get_vtx_wrt_sector(TTIHelper &ha, index_t va,
                                                     TTIHelper &hb, index_t eb)
{
	Sign va_wrt_ebn = get_v_wrt_seg(ha, va, hb, (eb + 1) % 3);
	Sign va_wrt_ebp = get_v_wrt_seg(ha, va, hb, (eb + 2) % 3);
	if (va_wrt_ebn == Sign::POSITIVE && va_wrt_ebp == Sign::NEGATIVE)
		return 1;
	if (va_wrt_ebn == Sign::NEGATIVE && va_wrt_ebp == Sign::POSITIVE)
		return 2;

	Sign va_wrt_eb = get_v_wrt_seg(ha, va, hb, eb);
	if ((va_wrt_eb >= Sign::ZERO && va_wrt_ebn >= Sign::ZERO &&
	     va_wrt_ebp >= Sign::ZERO) ||
	    (va_wrt_eb <= Sign::ZERO && va_wrt_ebn <= Sign::ZERO &&
	     va_wrt_ebp <= Sign::ZERO))
	{
		return 0;
	}

	Sign vb_wrt_eb =
	  OrientOn2D()(hb.v[eb], hb.v[(eb + 1) % 3], hb.v[(eb + 2) % 3], hb.t_nmax);

	return va_wrt_eb != vb_wrt_eb ? 0 : 3;
}

/**
 * @param ha helper for `ta`.
 * @param ea A valid local edge index in `ta`.
 * @param hb helper for `tb`.
 * @param eb A valid local edge index in `tb`.
 * @return bool true if two edges intersect or overlap.
 * @note it is specially designed for classfication, so it is not general.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::seg_seg_do_intersect(TTIHelper &ha, index_t ea,
                                                     TTIHelper &hb, index_t eb,
                                                     Sign eb0_wrt_ea,
                                                     Sign eb1_wrt_ea)
{
	index_t ea0 = ea, ea1 = (ea + 1) % 3;

	OMC_UNUSED index_t eb0 = eb, eb1 = (eb + 1) % 3;

	Sign ea0_wrt_eb = get_v_wrt_seg(ha, ea0, hb, eb);
	Sign ea1_wrt_eb = get_v_wrt_seg(ha, ea1, hb, eb);

	if (ea0_wrt_eb != ea1_wrt_eb && eb0_wrt_ea != eb1_wrt_ea)
	{
		// two edges sharing an endpoint are impossible in classfication pipeline.
#if 0
		if (ha.v_in_vtx[ea0] == eb0 || ha.v_in_vtx[ea0] == eb1 ||
		    ha.v_in_vtx[ea1] == eb0 || ha.v_in_vtx[ea1] == eb1)
		 	return false;
#else
		OMC_EXPENSIVE_ASSERT(!(ha.v_in_vtx[ea0] == eb0 || ha.v_in_vtx[ea0] == eb1 ||
		                       ha.v_in_vtx[ea1] == eb0 || ha.v_in_vtx[ea1] == eb1),
		                     "impossible case.");
#endif

		// intersect at a point
		return true;
	}

// degenerate case: collinear segments.
// such case is impossible in classification pipeline.
#if 0
	if (ea0_wrt_eb == Sign::ZERO && ea1_wrt_eb == Sign::ZERO &&
	    eb0_wrt_ea == Sign::ZERO && eb1_wrt_ea == Sign::ZERO)
	{
	  // coincident segments are impossible in classification pipeline,
	  // but we keep this check for completeness.
	  if ((ha.v_in_vtx[ea0] == eb0 && ha.v_in_vtx[ea1] == eb1) ||
	      (ha.v_in_vtx[ea0] == eb1 && ha.v_in_vtx[ea1] == eb0))
	    return false;

	  size_t i = 0, j = 1; // ha.t_nmax == 2
	  if (ha.t_nmax == 0)
	  {
	    i = 1, j = 2;
	  }
	  else if (ha.t_nmax == 1)
	  {
	    i = 0, j = 2;
	  }

	  NT min_ea_i = ha.v[ea0][i] < ha.v[ea1][i] ? ha.v[ea0][i] : ha.v[ea1][i];
	  NT max_ea_i = ha.v[ea0][i] < ha.v[ea1][i] ? ha.v[ea1][i] : ha.v[ea0][i];
	  NT min_ea_j = ha.v[ea0][j] < ha.v[ea1][j] ? ha.v[ea0][j] : ha.v[ea1][j];
	  NT max_ea_j = ha.v[ea0][j] < ha.v[ea1][j] ? ha.v[ea1][j] : ha.v[ea0][j];

	  NT min_eb_i = hb.v[eb0][i] < hb.v[eb1][i] ? hb.v[eb0][i] : hb.v[eb1][i];
	  NT max_eb_i = hb.v[eb0][i] < hb.v[eb1][i] ? hb.v[eb1][i] : hb.v[eb0][i];
	  NT min_eb_j = hb.v[eb0][j] < hb.v[eb1][j] ? hb.v[eb0][j] : hb.v[eb1][j];
	  NT max_eb_j = hb.v[eb0][j] < hb.v[eb1][j] ? hb.v[eb1][j] : hb.v[eb0][j];

	  if ( // test s0 endpoints against s1 range
	    (ha.v[ea0][i] > min_eb_i && ha.v[ea0][i] < max_eb_i) ||
	    (ha.v[ea0][j] > min_eb_j && ha.v[ea0][j] < max_eb_j) ||
	    (ha.v[ea1][i] > min_eb_i && ha.v[ea1][i] < max_eb_i) ||
	    (ha.v[ea1][j] > min_eb_j && ha.v[ea1][j] < max_eb_j) ||
	    // test s1 endpoints against s0 range
	    (hb.v[eb0][i] > min_ea_i && hb.v[eb0][i] < max_ea_i) ||
	    (hb.v[eb0][j] > min_ea_j && hb.v[eb0][j] < max_ea_j) ||
	    (hb.v[eb1][i] > min_ea_i && hb.v[eb1][i] < max_ea_i) ||
	    (hb.v[eb1][j] > min_ea_j && hb.v[eb1][j] < max_ea_j))
	  {
	    return true; // overlap
	  }
	}
#else
	OMC_EXPENSIVE_ASSERT(!(ea0_wrt_eb == Sign::ZERO && ea1_wrt_eb == Sign::ZERO &&
	                       eb0_wrt_ea == Sign::ZERO && eb1_wrt_ea == Sign::ZERO),
	                     "impossible case.");
#endif
	return false; // do not intersect
}

/**
 * @param ha helper for `ta`
 * @param ea a valid index of edge in `ta`
 * @param hb helper for `tb`
 * @return true if they do intersect
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::coplanar_seg_tri_do_intersect(TTIHelper &ha,
                                                              index_t    ea,
                                                              TTIHelper &hb)
{
	index_t eav0 = ea, eav1 = (ea + 1) % 3;

	Sign tb_ori = OrientOn2D()(hb.v[0], hb.v[1], hb.v[2], hb.t_nmax);

	Sign eav0_wrt_eb0 = get_v_wrt_seg(ha, eav0, hb, 0);
	Sign eav1_wrt_eb0 = get_v_wrt_seg(ha, eav1, hb, 0);

	if (eav0_wrt_eb0 != tb_ori && eav1_wrt_eb0 != tb_ori)
		return false; // ea and tb are separated by eb0

	Sign eav0_wrt_eb1 = get_v_wrt_seg(ha, eav0, hb, 1);
	Sign eav1_wrt_eb1 = get_v_wrt_seg(ha, eav1, hb, 1);

	if (eav0_wrt_eb1 != tb_ori && eav1_wrt_eb1 != tb_ori)
		return false; // ea and tb are separated by eb1

	Sign eav0_wrt_eb2 = get_v_wrt_seg(ha, eav0, hb, 2);
	Sign eav1_wrt_eb2 = get_v_wrt_seg(ha, eav1, hb, 2);

	if (eav0_wrt_eb2 != tb_ori && eav1_wrt_eb2 != tb_ori)
		return false; // ea and tb are separated by eb2

	Sign vb0_wrt_ea, vb1_wrt_ea, vb2_wrt_ea;
	if (ha.t_nmax == hb.t_nmax)
	{
		vb0_wrt_ea = get_v_wrt_seg(hb, 0, ha, ea);
		vb1_wrt_ea = get_v_wrt_seg(hb, 1, ha, ea);
		vb2_wrt_ea = get_v_wrt_seg(hb, 2, ha, ea);
	}
	else
	{
		vb0_wrt_ea = OrientOn2D()(ha.v[eav0], ha.v[eav1], hb.v[0], hb.t_nmax);
		vb1_wrt_ea = OrientOn2D()(ha.v[eav0], ha.v[eav1], hb.v[1], hb.t_nmax);
		vb2_wrt_ea = OrientOn2D()(ha.v[eav0], ha.v[eav1], hb.v[2], hb.t_nmax);
	}

	if ((vb0_wrt_ea > Sign::ZERO && vb1_wrt_ea > Sign::ZERO &&
	     vb2_wrt_ea > Sign::ZERO) ||
	    (vb0_wrt_ea < Sign::ZERO && vb1_wrt_ea < Sign::ZERO &&
	     vb2_wrt_ea < Sign::ZERO))
		return false; // tb is strictly on the one side of ea

	return true;
}

/**
 * @brief Given an edge `ea` from triangle `ta`, it is known that `ea` crosses
 * the support plane of `ta`, then this func tests if `ea` intersects `ta` at a
 * point except the vertices of `ta`.
 * @param ha helper for triangle `ta`
 * @param ea a valid local index of edge in `ta`
 * @param hb helper for triangle `tb`
 * @return true if `ea` and `tb` intersect.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::noncoplanar_seg_tri_do_intersect(TTIHelper &ha,
                                                                 index_t    ea,
                                                                 TTIHelper &hb)
{
	if (hb.v_in_seg[0] == ea || hb.v_in_seg[1] == ea || hb.v_in_seg[2] == ea)
		return true;

	Sign vol_ea_eb0 = get_seg_wrt_seg(ha, ea, hb, 0);
	Sign vol_ea_eb1 = get_seg_wrt_seg(ha, ea, hb, 1);
	if ((vol_ea_eb0 > Sign::ZERO && vol_ea_eb1 < Sign::ZERO) ||
	    (vol_ea_eb0 < Sign::ZERO && vol_ea_eb1 > Sign::ZERO))
		return false; // DO NOT INTERSECT
	Sign vol_ea_eb2 = get_seg_wrt_seg(ha, ea, hb, 2);
	if ((vol_ea_eb1 > Sign::ZERO && vol_ea_eb2 < Sign::ZERO) ||
	    (vol_ea_eb1 < Sign::ZERO && vol_ea_eb2 > Sign::ZERO))
		return false; // DO NOT INTERSECT
	if ((vol_ea_eb2 > Sign::ZERO && vol_ea_eb0 < Sign::ZERO) ||
	    (vol_ea_eb2 < Sign::ZERO && vol_ea_eb0 > Sign::ZERO))
		return false; // DO NOT INTERSECT

	return true;
}

/**
 * @param ha helper for `ta`
 * @param va A valid local vertex index in `ta`
 * @param hb helper for `tb`
 * @param inter_list list to store intersection points.
 * @param inter_listlist to store intersection points.
 * @return true if intersections are detected.
 * @note this func does not assume that `va` and `tb` intersect, it will check
 * and then classify.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::classify_coplanr_vtx_intersections(
  TTIHelper &ha, index_t va, TTIHelper &hb,
  phmap::flat_hash_set<index_t> &inter_list)
{
	OMC_EXPENSIVE_ASSERT(!is_valid_idx(ha.v_in_vtx[va]), "impossible case");

	if (TTIHelper::is_cached(ha.v_in_seg[va]) && is_valid_idx(ha.v_in_seg[va]))
	{ // strictly inside an edge of `tb`
		g.addVertexInEdge(hb.e_id[ha.v_in_seg[va]], ha.v_id[va]);
		inter_list.insert(ha.v_id[va]);
		return true;
	}

	Sign ori_tb = OrientOn2D()(hb.v[0], hb.v[1], hb.v[2], hb.t_nmax);

	for (size_t i = 0; i < 3; i++)
	{
		Sign va_wrt_eb0 = get_v_wrt_seg(ha, va, hb, i);
		if (va_wrt_eb0 != Sign::ZERO && va_wrt_eb0 != ori_tb)
			return false; // strictly outside `tb`
		if (va_wrt_eb0 == Sign::ZERO)
		{ // on the support line of ebi
			if (get_v_in_seg(ha, va, hb, i))
			{ // strictly inside ebi
				g.addVertexInEdge(hb.e_id[i], ha.v_id[va]);
				inter_list.insert(ha.v_id[va]);
				ha.v_in_seg[va] = i;
				return true;
			}
			else
				return false;
		}
	}

	// strictly inside `tb`
	g.addVertexInTriangle(hb.t_id, ha.v_id[va]);
	inter_list.insert(ha.v_id[va]);
	return true;
}

/**
 * @param ha helper for `ta`.
 * @param ea A valid local edge index in `ta`.
 * @param hb helper for `tb`.
 * @param inter_list list to store intersection points.
 * @param copl_edge_crosses list to store coplanar intersection points.
 * @return true if intersections are detected.
 * @note this func does not assume that `ea` and `tb` intersect, it will check
 * and then classify.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::classify_coplanar_edge_intersections(
  TTIHelper &ha, index_t ea, TTIHelper &hb,
  phmap::flat_hash_set<index_t>             &inter_list,
  phmap::flat_hash_set<CoplanarEEI, Hasher> &copl_edge_crosses)
{
	index_t ev0 = ea, ev1 = (ea + 1) % 3;
	index_t ev0_in_vtx = ha.v_in_vtx[ev0];
	index_t ev1_in_vtx = ha.v_in_vtx[ev1];

	OMC_EXPENSIVE_ASSERT(!(is_valid_idx(ev0_in_vtx) && is_valid_idx(ev1_in_vtx)),
	                     "two coincident vertices.");
	if (is_valid_idx(ev0_in_vtx) &&
	    is_valid_idx(ev1_in_vtx)) // Although this case should not happen,
		return false;               // we keep this check for completeness.

	index_t ev0_in_seg = get_v_in_seg(ha, ev0, hb);
	if (is_valid_idx(ev0_in_seg))
	{
		g.addVertexInEdge(hb.e_id[ev0_in_seg], ha.v_id[ev0]);
		inter_list.insert(ha.v_id[ev0]);
	}
	index_t ev1_in_seg = get_v_in_seg(ha, ev1, hb);
	if (is_valid_idx(ev1_in_seg))
	{
		g.addVertexInEdge(hb.e_id[ev1_in_seg], ha.v_id[ev1]);
		inter_list.insert(ha.v_id[ev1]);
	}

	if (is_valid_idx(ev0_in_seg) && is_valid_idx(ev1_in_seg))
	{
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha.t_id, hb.t_id);
		return true;
	}
	if (is_valid_idx(ev0_in_seg) && is_valid_idx(ev1_in_vtx))
	{
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha.t_id, hb.t_id);
		return true;
	}
	if (is_valid_idx(ev1_in_seg) && is_valid_idx(ev0_in_vtx))
	{
		add_symbolic_segment(ha.v_id[ev1], ha.v_id[ev0], ha.t_id, hb.t_id);
		return true;
	}

	bool ev1_in_tri = get_v_in_tri(ha, ev1, hb);
	if (ev1_in_tri)
	{
		g.addVertexInTriangle(hb.t_id, ha.v_id[ev1]);
		inter_list.insert(ha.v_id[ev1]);
	}

	// v0 in a segment or vtx and v1 inside triangle
	if ((is_valid_idx(ev0_in_seg) || is_valid_idx(ev0_in_vtx)) && ev1_in_tri)
	{
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha.t_id, hb.t_id);
		return true;
	}

	bool ev0_in_tri = get_v_in_tri(ha, ev0, hb);
	if (ev0_in_tri)
	{
		g.addVertexInTriangle(hb.t_id, ha.v_id[ev0]);
		inter_list.insert(ha.v_id[ev0]);
	}

	// v1 in a segment or vtx and v0 inside triangle
	if ((is_valid_idx(ev1_in_seg) || is_valid_idx(ev1_in_vtx)) && ev0_in_tri)
	{
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha.t_id, hb.t_id);
		return true;
	}

	// v0 and v1 both inside the triangle
	if (ev0_in_tri && ev1_in_tri)
	{
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha.t_id, hb.t_id);
		return true;
	}

	// Edges cross checking

	auto is_vb_in_ea = [&](index_t vb)
	{
		// check if `vb` is in any vertex of `ea`
		if (ha.t_nmax == hb.t_nmax)
		{
			if (hb.v_in_vtx[vb] == ev0 || hb.v_in_vtx[vb] == ev1)
				return std::pair<bool, Sign>(false, get_v_wrt_seg(hb, vb, ha, ea));
			if (TTIHelper::is_cached(hb.v_in_seg[vb]))
			{
				if (hb.v_in_seg[vb] == ea)
					return std::pair<bool, Sign>(true, get_v_wrt_seg(hb, vb, ha, ea));
				else
					return std::pair<bool, Sign>(false, get_v_wrt_seg(hb, vb, ha, ea));
			}

			if (get_v_wrt_seg(hb, vb, ha, ea) == Sign::ZERO)
			{ // `vb` is on the support line of edge `ea`
				if (get_v_in_seg(hb, vb, ha, ea))
				{ // `vb` strictly inside `ea`
					hb.v_in_seg[vb] = ea;
					return std::pair<bool, Sign>(true, get_v_wrt_seg(hb, vb, ha, ea));
				}
			}
			return std::pair<bool, Sign>(false, get_v_wrt_seg(hb, vb, ha, ea));
		}
		else
		{
			Sign vb_wrt_ea = OrientOn2D()(ha.v[ev0], ha.v[ev1], hb.v[vb], hb.t_nmax);
			if (hb.v_in_vtx[vb] == ev0 || hb.v_in_vtx[vb] == ev1)
				return std::pair<bool, Sign>(false, vb_wrt_ea);
			if (vb_wrt_ea == Sign::ZERO)
			{ // `vb` is on the support line of edge `ea`
				const NT *p  = hb.v[vb];
				const NT *s0 = ha.v[ea], *s1 = ha.v[(ea + 1) % 3];
				if ((p[0] > std::min(s0[0], s1[0]) && p[0] < std::max(s0[0], s1[0])) ||
				    (p[1] > std::min(s0[1], s1[1]) && p[1] < std::max(s0[1], s1[1])) ||
				    (p[2] > std::min(s0[2], s1[2]) && p[2] < std::max(s0[2], s1[2])))
				{ // `vb` strictly inside `ea`
					hb.v_wrt_seg[vb][ea] = Sign::ZERO;
					hb.v_in_seg[vb]      = ea;
					return std::pair<bool, Sign>(true, vb_wrt_ea);
				}
			}
			return std::pair<bool, Sign>(false, vb_wrt_ea);
		}
	};

	// we check only if seg A cross seg B and not B cross A (we found the
	// intersection once)
	auto [vb0_in_ea, vb0_wrt_ea] = is_vb_in_ea(/*vb*/ 0);
	auto [vb1_in_ea, vb1_wrt_ea] = is_vb_in_ea(/*vb*/ 1);
	auto [vb2_in_ea, vb2_wrt_ea] = is_vb_in_ea(/*vb*/ 2);

	index_t seg0_cross = InvalidIndex, seg1_cross = InvalidIndex,
	        seg2_cross = InvalidIndex;

	index_t glob_ea  = ha.e_id[ea];
	index_t glob_eb0 = hb.e_id[/*eb*/ 0];
	index_t glob_eb1 = hb.e_id[/*eb*/ 1];
	index_t glob_eb2 = hb.e_id[/*eb*/ 2];

	if (vb0_in_ea)
	{
		g.addVertexInEdge(glob_ea, hb.v_id[/*vb*/ 0]);
		inter_list.insert(hb.v_id[/*vb*/ 0]);
	}
	if (vb1_in_ea)
	{
		g.addVertexInEdge(glob_ea, hb.v_id[/*vb*/ 1]);
		inter_list.insert(hb.v_id[/*vb*/ 1]);
	}
	if (vb2_in_ea)
	{
		g.addVertexInEdge(glob_ea, hb.v_id[/*vb*/ 2]);
		inter_list.insert(hb.v_id[/*vb*/ 2]);
	}

	if (
	  /* endpoints of ea are not on eb0 or eb0's endpoints */
	  ev0_in_vtx != /*vb*/ 0 && ev0_in_vtx != /*vb*/ 1 &&
	  ev1_in_vtx != /*vb*/ 0 && ev1_in_vtx != /*vb*/ 1 &&
	  ev0_in_seg != /*eb*/ 0 && ev1_in_seg != /*eb*/ 0 &&
	  /* endpoints of eb0 are not on ea */
	  !vb0_in_ea && !vb1_in_ea &&
	  /* ea and eb0 cross*/
	  seg_seg_do_intersect(ha, ea, hb, /*eb*/ 0, vb0_wrt_ea, vb1_wrt_ea))
	{
		// edge `ea` cross seg 0 in `tb`
		seg0_cross = add_edge_cross_coplanar_edge(glob_eb0, glob_ea, hb.t_id,
		                                          copl_edge_crosses);
		inter_list.insert(seg0_cross);

		if (is_valid_idx(ev0_in_vtx) || is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{
			add_symbolic_segment(ha.v_id[ev0], seg0_cross, ha.t_id, hb.t_id);
			return true;
		}
		if (is_valid_idx(ev1_in_vtx) || is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{
			add_symbolic_segment(ha.v_id[ev1], seg0_cross, ha.t_id, hb.t_id);
			return true;
		}
		if (vb2_in_ea)
		{
			add_symbolic_segment(hb.v_id[/*vb*/ 2], seg0_cross, hb.t_id, ha.t_id);
			return true;
		}
	}

	if (/* endpoints of ea are not on eb1 or eb1's endpoints */
	    ev0_in_vtx != /*vb*/ 1 && ev0_in_vtx != /*vb*/ 2 &&
	    ev1_in_vtx != /*vb*/ 1 && ev1_in_vtx != /*vb*/ 2 &&
	    ev0_in_seg != /*eb*/ 1 && ev1_in_seg != /*eb*/ 1 &&
	    /* endpoints of eb1 are not on ea */
	    !vb1_in_ea && !vb2_in_ea &&
	    /* ea and eb1 cross*/
	    seg_seg_do_intersect(ha, ea, hb, /*eb*/ 1, vb1_wrt_ea, vb2_wrt_ea))
	{
		// edge `ea` cross seg 1 in `tb`
		seg1_cross = add_edge_cross_coplanar_edge(glob_eb1, glob_ea, hb.t_id,
		                                          copl_edge_crosses);
		inter_list.insert(seg1_cross);

		if (is_valid_idx(ev0_in_vtx) || is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{
			add_symbolic_segment(ha.v_id[ev0], seg1_cross, ha.t_id, hb.t_id);
			return true;
		}
		if (is_valid_idx(ev1_in_vtx) || is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{
			add_symbolic_segment(ha.v_id[ev1], seg1_cross, ha.t_id, hb.t_id);
			return true;
		}
		if (vb0_in_ea)
		{
			add_symbolic_segment(hb.v_id[/*vb*/ 0], seg1_cross, hb.t_id, ha.t_id);
			return true;
		}
	}

	if (/* endpoints of ea are not on eb2 or eb2's endpoints */
	    ev0_in_vtx != /*vb*/ 0 && ev0_in_vtx != /*vb*/ 2 &&
	    ev1_in_vtx != /*vb*/ 0 && ev1_in_vtx != /*vb*/ 2 &&
	    ev0_in_seg != /*eb*/ 2 && ev1_in_seg != /*eb*/ 2 &&
	    /* endpoints of eb2 are not on ea */
	    !vb0_in_ea && !vb2_in_ea &&
	    /* ea and eb2 cross*/
	    seg_seg_do_intersect(ha, ea, hb, /*eb*/ 2, vb2_wrt_ea, vb0_wrt_ea))
	{
		// edge `ea` cross seg 0 in `tb`
		seg2_cross = add_edge_cross_coplanar_edge(glob_eb2, glob_ea, hb.t_id,
		                                          copl_edge_crosses);
		inter_list.insert(seg2_cross);

		if (is_valid_idx(ev0_in_vtx) || is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{
			add_symbolic_segment(ha.v_id[ev0], seg2_cross, ha.t_id, hb.t_id);
			return true;
		}
		if (is_valid_idx(ev1_in_vtx) || is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{
			add_symbolic_segment(ha.v_id[ev1], seg2_cross, ha.t_id, hb.t_id);
			return true;
		}
		if (vb1_in_ea)
		{
			add_symbolic_segment(hb.v_id[/*vb*/ 1], seg2_cross, ha.t_id, hb.t_id);
			return true;
		}
	}

	// final probably symbolic edges
	if (is_valid_idx(seg0_cross) && is_valid_idx(seg1_cross))
	{
		add_symbolic_segment(seg0_cross, seg1_cross, ha.t_id, hb.t_id);
		return true;
	}
	if (is_valid_idx(seg0_cross) && is_valid_idx(seg2_cross))
	{
		add_symbolic_segment(seg0_cross, seg2_cross, ha.t_id, hb.t_id);
		return true;
	}
	if (is_valid_idx(seg1_cross) && is_valid_idx(seg2_cross))
	{
		add_symbolic_segment(seg1_cross, seg2_cross, ha.t_id, hb.t_id);
		return true;
	}

	if (vb0_in_ea && vb1_in_ea)
	{
		add_symbolic_segment(hb.v_id[0], hb.v_id[1], hb.t_id, ha.t_id);
		return true;
	}
	if (vb1_in_ea && vb2_in_ea)
	{
		add_symbolic_segment(hb.v_id[1], hb.v_id[2], hb.t_id, ha.t_id);
		return true;
	}
	if (vb2_in_ea && vb0_in_ea)
	{
		add_symbolic_segment(hb.v_id[2], hb.v_id[0], hb.t_id, ha.t_id);
		return true;
	}

	if (vb0_in_ea)
	{
		if (is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{
			add_symbolic_segment(hb.v_id[0], ha.v_id[ev0], hb.t_id, ha.t_id);
			return true;
		}
		if (is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{
			add_symbolic_segment(hb.v_id[0], ha.v_id[ev1], hb.t_id, ha.t_id);
			return true;
		}
	}

	if (vb1_in_ea)
	{
		if (is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{
			add_symbolic_segment(hb.v_id[1], ha.v_id[ev0], hb.t_id, ha.t_id);
			return true;
		}
		if (is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{
			add_symbolic_segment(hb.v_id[1], ha.v_id[ev1], hb.t_id, ha.t_id);
			return true;
		}
	}

	if (vb2_in_ea)
	{
		if (is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{
			add_symbolic_segment(hb.v_id[2], ha.v_id[ev0], hb.t_id, ha.t_id);
			return true;
		}
		if (is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{
			add_symbolic_segment(hb.v_id[2], ha.v_id[ev1], hb.t_id, ha.t_id);
			return true;
		}
	}

	return false;
}

/**
 * @param ha helper for `ta`.
 * @param ea A valid local edge index in `ta`.
 * @param hb helper for `tb`.
 * @param inter_list list to store intersection points.
 * @param inter_listlist to store non-coplanar intersection points.
 * @return true if intersections are detected.
 * @note this func assumes that `ea` and `tb` definitely intersect at a point
 * except the vertices of `tb`
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::classify_noncoplanar_edge_intersections(
  TTIHelper &ha, index_t ea, TTIHelper &hb,
  phmap::flat_hash_set<index_t> &inter_list)
{
	for (size_t i = 0; i < 3; i++)
		if (hb.v_in_seg[i] == ea)
		{
			g.addVertexInEdge(ha.e_id[ea], hb.v_id[i]);
			inter_list.insert(hb.v_id[i]);
			return true;
		}

	Sign vol_ea_eb0 = get_seg_wrt_seg(ha, ea, hb, 0);
	Sign vol_ea_eb2 = get_seg_wrt_seg(ha, ea, hb, 2);

	// `ea` cross vb0
	if (vol_ea_eb0 == Sign::ZERO && vol_ea_eb2 == Sign::ZERO)
	{
		g.addVertexInEdge(ha.e_id[ea], hb.v_id[0]);
		inter_list.insert(hb.v_id[0]);
		hb.v_in_seg[0] = ea;
		return true;
	}
	Sign vol_ea_eb1 = get_seg_wrt_seg(ha, ea, hb, 1);
	// `ea` cross vb1
	if (vol_ea_eb0 == Sign::ZERO && vol_ea_eb1 == Sign::ZERO)
	{
		g.addVertexInEdge(ha.e_id[ea], hb.v_id[1]);
		inter_list.insert(hb.v_id[1]);
		hb.v_in_seg[1] = ea;
		return true;
	}
	// `ea` cross vb2
	if (vol_ea_eb1 == Sign::ZERO && vol_ea_eb2 == Sign::ZERO)
	{
		g.addVertexInEdge(ha.e_id[ea], hb.v_id[2]);
		inter_list.insert(hb.v_id[2]);
		hb.v_in_seg[2] = ea;
		return true;
	}
	// the edge intersects the tri in seg 0
	if (vol_ea_eb0 == Sign::ZERO)
	{
		index_t int_point =
		  add_edge_cross_noncoplanar_edge(ha.e_id[ea], hb.e_id[0]);
		inter_list.insert(int_point);
		return true;
	}
	// the edge intersects the tri in seg 1
	if (vol_ea_eb1 == Sign::ZERO)
	{
		index_t int_point =
		  add_edge_cross_noncoplanar_edge(ha.e_id[ea], hb.e_id[1]);
		inter_list.insert(int_point);
		return true;
	}
	// the edge intersects the tri in seg 2
	if (vol_ea_eb2 == Sign::ZERO)
	{
		index_t int_point =
		  add_edge_cross_noncoplanar_edge(ha.e_id[ea], hb.e_id[2]);
		inter_list.insert(int_point);
		return true;
	}

	// the edge intersects the inner triangle
	{
		index_t int_point = add_edge_cross_tri(ha.e_id[ea], hb.t_id);
		inter_list.insert(int_point);
	}
	return true;
}

/**
 * @param v0 global index of a vertex
 * @param v1 global index of a vertex
 * @param ta global index of a triangle
 * @param tb global index of a triangle
 */
template <typename Traits>
void DetectClassifyTTI<Traits>::add_symbolic_segment(index_t v0, index_t v1,
                                                     index_t ta, index_t tb)
{
	OMC_EXPENSIVE_ASSERT(v0 != v1, "trying to add a 0-lenght symbolic edge");
	UIPair seg = uniquePair(v0, v1);

	if (!ts.triContainsEdge(ta, v0, v1))
		g.addSegmentInTriangle(ta, seg);

	if (!ts.triContainsEdge(tb, v0, v1))
		g.addSegmentInTriangle(tb, seg);

	g.addTrianglesInSegment(seg, ta, tb);
}

/**
 * @param ea global index of an edge
 * @param eb global index of an edge
 * @param t global index of a triangle
 * @param copl_edge_crosses coplanar intersection points
 * @return index_t global index of the new implicit point
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_edge_cross_coplanar_edge(
  index_t ea, index_t eb, index_t t,
  phmap::flat_hash_set<CoplanarEEI, Hasher> &copl_edge_crosses)
{
	// find intersection point in existed points
	for (const CoplanarEEI &c : copl_edge_crosses)
		if (c.is_same(ea, eb))
			return c.p;
	// otherwise find it in global point set or create a new point

	IPoint_SSI *new_v = pnt_arena.emplace(
	  CreateSSI()(ts.edgeVert(ea, 0), ts.edgeVert(ea, 1), ts.edgeVert(eb, 0),
	              ts.edgeVert(eb, 1), ts.triPlane(t)));

	std::atomic<index_t> *idx_ptr = idx_arena.emplace(InvalidIndex);
	// index creation is deferred until succeesfully inserting point.

	// <index of point, succeed to insert?>
	std::pair<index_t, bool> ins =
	  uniq_g.addVertexInSortedList(new_v, idx_ptr, CreateIndex(ts, uniq_g));
	if (!ins.second)
	{
		OMC_ASSERT(ins.first != InvalidIndex, "");
		IPoint_SSI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
		idx_arena.recycle(idx_ptr);
	}

	g.addVertexInEdge(ea, ins.first);
	g.addVertexInEdge(eb, ins.first);
	copl_edge_crosses.insert(CoplanarEEI(ea, eb, ins.first));

	return ins.first;
}

/**
 * @param ea global index of an edge
 * @param eb global index of an edge
 * @return global index of the intersection point
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_edge_cross_noncoplanar_edge(index_t ea,
                                                                   index_t eb)
{
	// find a plane for two interseted segments
	int plane = -1;
	{
		static std::array<index_t, 12> _tri = {0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};

		std::array<const NT *, 4> _pnts = {
		  ts.edgeVert(ea, 0).data(), ts.edgeVert(ea, 1).data(),
		  ts.edgeVert(eb, 0).data(), ts.edgeVert(eb, 1).data()};

		for (size_t i = 0; i < 12; i += 3)
		{
			if (CollinearPoints3D().misaligned(_pnts[_tri[i]], _pnts[_tri[i + 1]],
			                                   _pnts[_tri[i + 2]]))
			{
				plane = MaxCompInTriNormal()(_pnts[_tri[i]], _pnts[_tri[i + 1]],
				                             _pnts[_tri[i + 2]]);
				break;
			}
		}
		OMC_ASSERT(plane != -1, "segments can't be projected to 2D.");
	}

	IPoint_SSI *new_v = pnt_arena.emplace(
	  CreateSSI()(ts.edgeVert(ea, 0), ts.edgeVert(ea, 1), ts.edgeVert(eb, 0),
	              ts.edgeVert(eb, 1), plane));

	std::atomic<index_t> *idx_ptr = idx_arena.emplace(InvalidIndex);
	// index creation is deferred until succeesfully inserting point.

	// <index of point, succeed to insert?>
	std::pair<index_t, bool> ins =
	  uniq_g.addVertexInSortedList(new_v, idx_ptr, CreateIndex(ts, uniq_g));

	if (!ins.second)
	{
		OMC_ASSERT(ins.first != InvalidIndex, "");
		IPoint_SSI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
		idx_arena.recycle(idx_ptr);
	}

	g.addVertexInEdge(ea, ins.first);
	g.addVertexInEdge(eb, ins.first);

	return ins.first;
}

template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_edge_cross_tri(index_t ea, index_t tb)
{
	IPoint_LPI *new_v = pnt_arena.emplace(
	  CreateLPI()(ts.edgeVert(ea, 0), ts.edgeVert(ea, 1), ts.triVert(tb, 0),
	              ts.triVert(tb, 1), ts.triVert(tb, 2)));
	std::atomic<index_t> *idx_ptr = idx_arena.emplace(InvalidIndex);
	// index creation is deferred until succeesfully inserting point.

	// <index of LPI point, succeed to insert?>
	std::pair<index_t, bool> ins =
	  uniq_g.addVertexInSortedList(new_v, idx_ptr, CreateIndex(ts, uniq_g));

	if (!ins.second)
	{
		OMC_ASSERT(ins.first != InvalidIndex, "");
		IPoint_LPI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
		idx_arena.recycle(idx_ptr);
	}

	g.addVertexInTriangle(tb, ins.first);
	g.addVertexInEdge(ea, ins.first);

	return ins.first;
}

} // namespace OMC