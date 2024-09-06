#pragma once

#include "DetectClassifyTTI.h"

namespace OMC {

#ifdef OMC_ARR_PROF_TTI
	#define OMC_ARR_PROF_TTI_INCT OMC_ARR_PROFILE_INC_TOTAL(ArrFuncNames::DC_TTI)
	#define OMC_ARR_PROF_TTI_INCR(branch) \
		OMC_ARR_PROFILE_INC_REACH(ArrFuncNames::DC_TTI, branch)
#else
	#define OMC_ARR_PROF_TTI_INCT
	#define OMC_ARR_PROF_TTI_INCR(branch)
#endif

template <typename Traits>
struct DetectClassifyTTI<Traits>::TTIHelper
{
	static const index_t UncachedIndex   = (index_t)(-2);
	static const char    UncachedBoolean = (char)(-1);

	static bool is_cached(index_t idx) { return idx != UncachedIndex; }
	static bool is_cached(char boolean) { return boolean != UncachedBoolean; }

	index_t t_id;
	int     t_nmax;

	// vectices from triangle
	std::array<EPoint, 3>     p;
	// pointers to vertices.
	std::array<const NT *, 3> v;

	// indices of vertices from triangle.
	std::array<index_t, 3> v_id;
	// indices of edges from triangle.
	std::array<index_t, 3> e_id;

	// position of sub-simplex of this triangle with respect to sub-simplex of
	// another triangle
	std::array<index_t, 3> v_in_vtx;
	std::array<index_t, 3> v_in_seg;
	std::array<char, 3>    v_in_tri;

	// cached orientations
	// Sign::UNCERTAIN means "not cached/calculated"
	std::array<std::array<Sign, 3>, 3> v_wrt_seg;
	std::array<std::array<Sign, 3>, 3> seg_wrt_seg;

	// v_in_vtx will always be set.
	void init_v_in_seg()
	{
		v_in_seg = {UncachedIndex, UncachedIndex, UncachedIndex};
	}
	void init_v_in_tri()
	{
		v_in_tri = {UncachedBoolean, UncachedBoolean, UncachedBoolean};
	}

	void init_v_wrt_seg()
	{
		v_wrt_seg = {
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN})};
	}
	void init_seg_wrt_seg()
	{
		seg_wrt_seg = {
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN})};
	}

	bool contains_vtx(index_t vidx)
	{
		return v_id[0] == vidx || v_id[1] == vidx || v_id[2] == vidx;
	}
};

template <typename Traits>
struct DetectClassifyTTI<Traits>::CoplanarEEI // Edge-Edge-Intersection
{
	index_t ea; // global index of edge a
	index_t eb; // global index of edge b
	index_t p;  // global index of the intersection point

	CoplanarEEI()  = default;
	~CoplanarEEI() = default;

	CoplanarEEI(index_t _ea, index_t _eb, index_t _p)
	  : ea(_ea)
	  , eb(_eb)
	  , p(_p)
	{
	}

	// used to check if two intersection points are same.
	// one is the stored intersection point, another is the intersection point
	// between query edges.
	bool is_same(index_t query_ea, index_t query_eb) const
	{
		return (query_ea == ea && query_eb == eb) ||
		       (query_ea == eb && query_eb == ea);
	}
};

template <typename Traits>
DetectClassifyTTI<Traits>::DetectClassifyTTI(TriSoup &_ts, PntArena &_pnt_arena)
  : ts(_ts)
  , pnt_arena(_pnt_arena)
{
}

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI(index_t ta, index_t tb)
{
	TTIHelper ha, hb;

	// vertices from two triangles `ta` and `tb`.
	ha.t_id     = ta;
	ha.v_id     = {ts.triVertID(ta, 0), ts.triVertID(ta, 1), ts.triVertID(ta, 2)};
	ha.p        = {AsEP()(ts.vert(ha.v_id[0])), AsEP()(ts.vert(ha.v_id[1])),
	               AsEP()(ts.vert(ha.v_id[2]))};
	ha.v        = {ha.p[0].data(), ha.p[1].data(), ha.p[2].data()};
	ha.e_id     = {InvalidIndex, InvalidIndex, InvalidIndex};
	ha.v_in_vtx = {InvalidIndex, InvalidIndex, InvalidIndex};

	hb.t_id     = tb;
	hb.v_id     = {ts.triVertID(tb, 0), ts.triVertID(tb, 1), ts.triVertID(tb, 2)};
	hb.p        = {AsEP()(ts.vert(hb.v_id[0])), AsEP()(ts.vert(hb.v_id[1])),
	               AsEP()(ts.vert(hb.v_id[2]))};
	hb.v        = {hb.p[0].data(), hb.p[1].data(), hb.p[2].data()};
	hb.e_id     = {InvalidIndex, InvalidIndex, InvalidIndex};
	hb.v_in_vtx = {InvalidIndex, InvalidIndex, InvalidIndex};

	// check if any triangle is degenerate
	OMC_EXPENSIVE_ASSERT(
	  CollinearPoints3D().misaligned(ha.v[0], ha.v[1], ha.v[2]) &&
	    CollinearPoints3D().misaligned(hb.v[0], hb.v[1], hb.v[2]),
	  "Detect degenerate triangle in check TTI.");
	OMC_ARR_PROF_TTI_INCT;

	// classify two triangles to three cases:
	// sharing edge, sharing vertex or sharing nothing.
	size_t correspond_count = find_vtx_correspondence(ha, hb);

	if (correspond_count == 3) // `ta` and `tb` are coincident.
		return; // (Although such case is impossible in arrangements pipeline
		        // because we have removed duplicate triangles, we still keep this
		        // check for completeness.)

	if (correspond_count == 2) // `ta` and `tb` share an edge.
		check_TTI_share_edge(ha, hb);
	else if (correspond_count == 1) // `ta` and `tb` share a vertex.
		check_TTI_share_vertex(ha, hb);
	else // `ta` and `tb` are separate
		check_TTI_separate(ha, hb);
}

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI_share_edge(TTIHelper &ha,
                                                     TTIHelper &hb)
{
	OMC_ARR_PROF_TTI_INCR(0);
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
	if (get_v_wrt_tri(hb, oppb, ha) != Sign::ZERO)
		return; // do not intersect

	OMC_ARR_PROF_TTI_INCR(1);

	// now `ta` and `tb` are coplanar
	ha.t_nmax = ts.triPlane(ha.t_id);
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

	// otherwise they intersect.
	// each triangle has one sharing edge and other two coplanar edges.
	ts.addCoplanarTriangles(ha.t_id, hb.t_id);
	ts.setTriangleHasIntersections(ha.t_id, hb.t_id);

	OMC_ARR_PROF_TTI_INCR(2);

	// cache the orientation result
	ha.init_v_wrt_seg();
	ha.init_v_in_seg();
	ha.init_v_in_tri();
	hb.init_v_wrt_seg();
	hb.init_v_in_seg();
	hb.init_v_in_tri();

	hb.v_wrt_seg[oppb][ea] = oppb_wrt_ea;
	ha.v_wrt_seg[oppa][eb] =
	  ha.v_id[ea] == hb.v_id[eb] ? oppa_wrt_ea : reverse_sign(oppa_wrt_ea);

	// only coplanar intersection points (edge-edge)
	CoplanarEEIList cec; // short name of coplanar_edge_crosses :)

	// four checks for intersections between a coplanar edge and a triangle.
	// edges from `ta` (except the shared edge) with respect to `tb`
	classify_coplanar_edge_intersections(ha, (ea + 1) % 3, hb, &cec);
	classify_coplanar_edge_intersections(ha, (ea + 2) % 3, hb, &cec);
	// edges from `tb` (except the shared edge) with respect to `ta`
	classify_coplanar_edge_intersections(hb, (eb + 1) % 3, ha, &cec);
	classify_coplanar_edge_intersections(hb, (eb + 2) % 3, ha, &cec);

	OMC_ARR_PROF_TTI_INCR(3);

	return; // end check.
}

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI_share_vertex(TTIHelper &ha,
                                                       TTIHelper &hb)
{
	OMC_ARR_PROF_TTI_INCR(4);
	// check if `ta` and `tb` intersect.

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

#ifdef OMC_ARR_DC_FILTER_O3D
	Sign eav0_wrt_tb = Orient3D().filter(hb.v[0], hb.v[1], hb.v[2], ha.v[eav0]);
	Sign eav1_wrt_tb = Orient3D().filter(hb.v[0], hb.v[1], hb.v[2], ha.v[eav1]);

	bool t_nmax_init    = false;
	bool v_wrt_seg_init = false;

	if (eav0_wrt_tb == Sign::UNCERTAIN && eav1_wrt_tb == Sign::UNCERTAIN)
	{ // it is possible that two triangles are coplanar.
		// we avoid expensive exact orient3d and filter out the no intersection case
		// by orient2d.
		ha.t_nmax   = ts.triPlane(ha.t_id);
		hb.t_nmax   = ts.triPlane(hb.t_id);
		t_nmax_init = true;
		if (ha.t_nmax == hb.t_nmax)
		{
			v_wrt_seg_init = true;
			if (!fast_check_on2d_share_vertex(ha, ea, hb, eb))
				return; // there is definitely no intersection
		}
	}
	// fail to filter out, go back to original pipeline
	if (eav0_wrt_tb == Sign::UNCERTAIN)
		eav0_wrt_tb = get_v_wrt_tri(ha, eav0, hb);
	if (eav1_wrt_tb == Sign::UNCERTAIN)
		eav1_wrt_tb = get_v_wrt_tri(ha, eav1, hb);
#else
	Sign eav0_wrt_tb = get_v_wrt_tri(ha, eav0, hb);
	Sign eav1_wrt_tb = get_v_wrt_tri(ha, eav1, hb);
#endif

	if (eav0_wrt_tb != Sign::ZERO && eav0_wrt_tb == eav1_wrt_tb)
		return; // above or below `tb`, do not intersect

	if (eav0_wrt_tb == Sign::ZERO && eav1_wrt_tb == Sign::ZERO)
	{
		OMC_ARR_PROF_TTI_INCR(5);

		{ // initialize cache
#ifdef OMC_ARR_DC_FILTER_O3D
			if (!t_nmax_init)
			{
				ha.t_nmax = ts.triPlane(ha.t_id);
				hb.t_nmax = ha.t_nmax;
			}
			if (!v_wrt_seg_init)
			{
				ha.init_v_wrt_seg();
				hb.init_v_wrt_seg();
			}
#else
			ha.t_nmax = ts.triPlane(ha.t_id);
			hb.t_nmax = ha.t_nmax;
			ha.init_v_wrt_seg();
			hb.init_v_wrt_seg();
#endif
		}

		size_t eav0_wrt_sector = get_vtx_wrt_sector(ha, eav0, hb, eb);
		size_t eav1_wrt_sector = get_vtx_wrt_sector(ha, eav1, hb, eb);

		if ((eav0_wrt_sector == 1 || eav0_wrt_sector == 3) &&
		    (eav1_wrt_sector == 1 || eav1_wrt_sector == 3))
			return; // ta and tb do not intersect
		if ((eav0_wrt_sector == 2 || eav0_wrt_sector == 3) &&
		    (eav1_wrt_sector == 2 || eav1_wrt_sector == 3))
			return; // ta and tb do not intersect
		OMC_ARR_PROF_TTI_INCR(6);

		size_t ebv0_wrt_sector = get_vtx_wrt_sector(hb, ebv0, ha, ea);
		size_t ebv1_wrt_sector = get_vtx_wrt_sector(hb, ebv1, ha, ea);

		if ((ebv0_wrt_sector == 1 || ebv0_wrt_sector == 3) &&
		    (ebv1_wrt_sector == 1 || ebv1_wrt_sector == 3))
			return; // ta and tb do not intersect
		if ((ebv0_wrt_sector == 2 || ebv0_wrt_sector == 3) &&
		    (ebv1_wrt_sector == 2 || ebv1_wrt_sector == 3))
			return; // ta and tb do not intersect
		OMC_ARR_PROF_TTI_INCR(7);

		// ea is coplanar to `tb`, so `ta` and `tb` are coplanar
		size_t ic = 0; // intersection count, excluding shared vertex

		CoplanarEEIList cec; // short name of coplanar_edge_crosses;

		{ // initialize cache
			ha.init_v_in_seg();
			ha.init_v_in_tri();
			hb.init_v_in_seg();
			hb.init_v_in_tri();
		}

		if (eav0_wrt_sector == 0) // previous edge of ea intersects tb
			ic += classify_coplanar_edge_intersections(ha, (ea + 2) % 3, hb, &cec);
		if (eav1_wrt_sector == 0) // next edge of ea intersects tb
			ic += classify_coplanar_edge_intersections(ha, (ea + 1) % 3, hb, &cec);
		if (ebv0_wrt_sector == 0) // previous edge of eb intersects ta
			ic += classify_coplanar_edge_intersections(hb, (eb + 2) % 3, ha, &cec);
		if (ebv1_wrt_sector == 0) // next edge of eb intersects ta
			ic += classify_coplanar_edge_intersections(hb, (eb + 1) % 3, ha, &cec);
		ic += classify_coplanar_edge_intersections(ha, ea, hb, &cec);
		ic += classify_coplanar_edge_intersections(hb, eb, ha, &cec);

		if (ic != 0)
		{
			OMC_ARR_PROF_TTI_INCR(8);
			ts.addCoplanarTriangles(ha.t_id, hb.t_id);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}

		return;
	}
	OMC_ARR_PROF_TTI_INCR(9);

	Sign ebv0_wrt_ta = get_v_wrt_tri(hb, ebv0, ha);
	Sign ebv1_wrt_ta = get_v_wrt_tri(hb, ebv1, ha);

	if (ebv0_wrt_ta != Sign::ZERO && ebv0_wrt_ta == ebv1_wrt_ta)
		return; // above or below `tb`, do not intersect
	OMC_ARR_PROF_TTI_INCR(10);

	// otherwise they cross each other's plane, possibly intersect.

	// go on classifying intersections on `tb`
	{ // initialize cache
#ifdef OMC_ARR_DC_FILTER_O3D
		if (!t_nmax_init)
		{
			ha.t_nmax = ts.triPlane(ha.t_id);
			hb.t_nmax = ts.triPlane(hb.t_id);
		}
		if (!v_wrt_seg_init)
		{
			ha.init_v_wrt_seg();
			hb.init_v_wrt_seg();
		}
#else
		ha.t_nmax = ts.triPlane(ha.t_id);
		hb.t_nmax = ts.triPlane(hb.t_id);
		ha.init_v_wrt_seg();
		hb.init_v_wrt_seg();
#endif
		ha.init_v_in_seg();
		ha.init_v_in_tri();
		ha.init_seg_wrt_seg();
		hb.init_v_in_seg();
		hb.init_v_in_tri();
		hb.init_seg_wrt_seg();
	}

	if ((eav0_wrt_tb == Sign::POSITIVE && eav1_wrt_tb == Sign::NEGATIVE) ||
	    (eav0_wrt_tb == Sign::NEGATIVE && eav1_wrt_tb == Sign::POSITIVE))
	{ // ea cross the support plane of `tb`
		OMC_ARR_PROF_TTI_INCR(11);

		// intersection list
		IntersectionPoints intersection_points;
		IntersectionTypes  intersection_types; // won't be use
		// the shared vertex is seen as an intersection point.
		intersection_points.push_back(ha.v_id[va]);

		// check if a new intersection point exists between `ea` and `tb`.
		if (noncoplanar_seg_tri_do_intersect(ha, ea, hb))
		{
			OMC_ARR_PROF_TTI_INCR(12);
			classify_noncoplanar_edge_intersections(ha, ea, hb, intersection_points,
			                                        intersection_types);
		}

		OMC_EXPENSIVE_ASSERT(
		  intersection_points.size() <= 2,
		  "more than 2 intersection points between non-coplanar triangles.");
		// If we have detected more than one intersections in non-coplanar cases,
		// we end classification now and propagate undetected intersections later.
		if (intersection_points.size() == 2)
		{
			index_t v0_id = intersection_points[0];
			index_t v1_id = intersection_points[1];

			// before ending, we record possible coplanar edges.
			if (ebv0_wrt_ta == Sign::ZERO || ebv1_wrt_ta == Sign::ZERO)
			{ // `tb` has a single coplanar edge to `ta`
				index_t copl_eb =
				  ebv0_wrt_ta == Sign::ZERO ? (eb + 2) % 3 : (eb + 1) % 3;

				OMC_EXPENSIVE_ASSERT_AUX_CODE(
				  index_t copl_vb = ebv0_wrt_ta == Sign::ZERO ? eb : (eb + 1) % 3);
				OMC_EXPENSIVE_ASSERT(get_vtx_wrt_sector(hb, copl_vb, ha, ea) == 0,
				                     "wrong coplanar edge");

				ts.addCoplanarEdge(ha.t_id, get_e_id(hb, copl_eb), v0_id, v1_id);
			}

			add_symbolic_segment(v0_id, v1_id, ha, InvalidIndex, hb, InvalidIndex);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);

			OMC_ARR_PROF_TTI_INCR(13);
			return; // all possible intersections are found, return.
		}
	}
	else if (eav0_wrt_tb == Sign::ZERO || eav1_wrt_tb == Sign::ZERO)
	{ // one edge of `ta` on the support plane of `tb`
		OMC_ARR_PROF_TTI_INCR(14);
		index_t copl_ea = eav0_wrt_tb == Sign::ZERO ? (ea + 2) % 3 : (ea + 1) % 3;
		index_t copl_va = eav0_wrt_tb == Sign::ZERO ? ea : (ea + 1) % 3;

		if (get_vtx_wrt_sector(ha, copl_va, hb, eb) == 0)
		{
			OMC_ARR_PROF_TTI_INCR(15);
			classify_coplanar_edge_intersections(ha, copl_ea, hb, nullptr);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else
	{
		OMC_ASSERT(false, "Impossible case happens.");
	}

	// go on classifying intersections on `ta`
	if ((ebv0_wrt_ta == Sign::POSITIVE && ebv1_wrt_ta == Sign::NEGATIVE) ||
	    (ebv0_wrt_ta == Sign::NEGATIVE && ebv1_wrt_ta == Sign::POSITIVE))
	{ // `eb` cross the support plane of `ta`
		OMC_ARR_PROF_TTI_INCR(16);

		// intersection list
		IntersectionPoints intersection_points;
		IntersectionTypes  intersection_types; // won't be use
		// the shared vertex is seen as an intersection point.
		intersection_points.push_back(hb.v_id[vb]);

		// check if a new intersection point exists between ea and `tb`.
		if (noncoplanar_seg_tri_do_intersect(hb, eb, ha))
		{
			OMC_ARR_PROF_TTI_INCR(17);
			classify_noncoplanar_edge_intersections(hb, eb, ha, intersection_points,
			                                        intersection_types);
		}

		OMC_EXPENSIVE_ASSERT(
		  intersection_points.size() <= 2,
		  "more than 2 intersection points between non-coplanar triangles.");
		// If we have detected more than one intersections in non-coplanar cases,
		// we end classification now and propagate undetected intersections later.
		if (intersection_points.size() == 2)
		{
			index_t v0_id = intersection_points[0];
			index_t v1_id = intersection_points[1];

			add_symbolic_segment(v0_id, v1_id, ha, InvalidIndex, hb, InvalidIndex);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else if (ebv0_wrt_ta == Sign::ZERO || ebv1_wrt_ta == Sign::ZERO)
	{ // one edge of `tb` on the support plane of `ta`
		OMC_ARR_PROF_TTI_INCR(18);
		index_t copl_eb = ebv0_wrt_ta == Sign::ZERO ? (eb + 2) % 3 : (eb + 1) % 3;
		index_t copl_vb = ebv0_wrt_ta == Sign::ZERO ? eb : (eb + 1) % 3;

		if (get_vtx_wrt_sector(hb, copl_vb, ha, ea) == 0)
		{
			OMC_ARR_PROF_TTI_INCR(19);
			classify_coplanar_edge_intersections(hb, copl_eb, ha, nullptr);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}
		return; // whether intersections are found or not, return.
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
template<typename Traits>
bool DetectClassifyTTI<Traits>::intersection_on_one_edge(
  const IntersectionTypes &intersection_types, index_t &edge_id)
{
	PointInSimplexType t0 = intersection_types[0];
	PointInSimplexType t1 = intersection_types[1];
	if ((t0 == PointInSimplexType::ON_VERT0 || t0 == PointInSimplexType::ON_VERT1 || t0 == PointInSimplexType::ON_EDGE0) &&
	    (t1 == PointInSimplexType::ON_VERT0 || t1 == PointInSimplexType::ON_VERT1 || t1 == PointInSimplexType::ON_EDGE0))
	{
		edge_id = 0;
		return true;
	}
	if ((t0 == PointInSimplexType::ON_VERT1 || t0 == PointInSimplexType::ON_VERT2 || t0 == PointInSimplexType::ON_EDGE1) &&
	    (t1 == PointInSimplexType::ON_VERT1 || t1 == PointInSimplexType::ON_VERT2 || t1 == PointInSimplexType::ON_EDGE1))
	{
		edge_id = 1;
		return true;
	}
	if ((t0 == PointInSimplexType::ON_VERT0 || t0 == PointInSimplexType::ON_VERT2 || t0 == PointInSimplexType::ON_EDGE2) &&
	    (t1 == PointInSimplexType::ON_VERT0 || t1 == PointInSimplexType::ON_VERT2 || t1 == PointInSimplexType::ON_EDGE2))
	{
		edge_id = 2;
		return true;
	}
	edge_id = InvalidIndex;
	return false;
}
// clang-format on

template <typename Traits>
void DetectClassifyTTI<Traits>::check_TTI_separate(TTIHelper &ha, TTIHelper &hb)
{
	// `ta` and `tb` do not share any edge or vertex.
	OMC_ARR_PROF_TTI_INCR(20);

	// intersection list
	IntersectionPoints intersection_points;
	IntersectionTypes  intersection_types;

	// =========================================================================
	// check of A respect to B
	// =========================================================================

	index_t edge_id, vtx_id;

	Sign orAB[3]; // orientation of edge of `ta` with respect to `tb`
#ifdef OMC_ARR_DC_FILTER_O3D
	orAB[0] = Orient3D().filter(hb.v[0], hb.v[1], hb.v[2], ha.v[0]);
	orAB[1] = Orient3D().filter(hb.v[0], hb.v[1], hb.v[2], ha.v[1]);
	orAB[2] = Orient3D().filter(hb.v[0], hb.v[1], hb.v[2], ha.v[2]);

	bool t_nmax_init    = false;
	bool v_wrt_seg_init = false;

	if (orAB[0] == Sign::UNCERTAIN && orAB[1] == Sign::UNCERTAIN &&
	    orAB[2] == Sign::UNCERTAIN)
	{ // it is possible that two triangles are coplanar.
		// we avoid expensive exact orient3d and filter out the no intersection case
		// by orient2d.
		t_nmax_init = true;
		ha.t_nmax   = ts.triPlane(ha.t_id);
		hb.t_nmax   = ts.triPlane(hb.t_id);
		if (ha.t_nmax == hb.t_nmax)
		{
			v_wrt_seg_init = true;
			if (!fast_check_on2d_separate(ha, hb))
				return; // there is definitely no intersection
		}
	}
	// go back to original pipeline
	if (orAB[0] == Sign::UNCERTAIN)
		orAB[0] = get_v_wrt_tri(ha, 0, hb);
	if (orAB[1] == Sign::UNCERTAIN)
		orAB[1] = get_v_wrt_tri(ha, 1, hb);
	if (orAB[2] == Sign::UNCERTAIN)
		orAB[2] = get_v_wrt_tri(ha, 2, hb);
#else
	orAB[0] = get_v_wrt_tri(ha, 0, hb);
	orAB[1] = get_v_wrt_tri(ha, 1, hb);
	orAB[2] = get_v_wrt_tri(ha, 2, hb);
#endif

	if (_sameOrientation(orAB[0], orAB[1]) &&
	    _sameOrientation(orAB[1], orAB[2]) && (orAB[0] != Sign::ZERO))
	{
		// CASE: no intersection found
		OMC_ARR_PROF_TTI_INCR(21);
		return;
	}

	if (_allCoplanarEdges(orAB))
	{
		OMC_ARR_PROF_TTI_INCR(22);
		// CASE: all edge of ta are coplanar to all edges of tb   (orAB: 0 0 0)

		{ // initialize cache
#ifdef OMC_ARR_DC_FILTER_O3D
			if (!t_nmax_init)
			{
				ha.t_nmax = ts.triPlane(ha.t_id);
				hb.t_nmax = ha.t_nmax;
			}
			if (!v_wrt_seg_init)
			{
				ha.init_v_wrt_seg();
				hb.init_v_wrt_seg();
			}
#else
			ha.t_nmax = ts.triPlane(ha.t_id);
			hb.t_nmax = ha.t_nmax;
			ha.init_v_wrt_seg();
			hb.init_v_wrt_seg();
#endif
			ha.init_v_in_seg();
			ha.init_v_in_tri();
			hb.init_v_in_seg();
			hb.init_v_in_tri();
		}
		// only coplanar intersection points (edge-edge)
		CoplanarEEIList cec;    // short name of coplanar_edge_crosses;
		size_t          ic = 0; // intersection count

		if (coplanar_seg_tri_do_intersect(ha, 0, hb))
			ic += classify_coplanar_edge_intersections(ha, 0, hb, &cec);
		if (coplanar_seg_tri_do_intersect(ha, 1, hb))
			ic += classify_coplanar_edge_intersections(ha, 1, hb, &cec);
		if (coplanar_seg_tri_do_intersect(ha, 2, hb))
			ic += classify_coplanar_edge_intersections(ha, 2, hb, &cec);

		if (coplanar_seg_tri_do_intersect(hb, 0, ha))
			ic += classify_coplanar_edge_intersections(hb, 0, ha, &cec);
		if (coplanar_seg_tri_do_intersect(hb, 1, ha))
			ic += classify_coplanar_edge_intersections(hb, 1, ha, &cec);
		if (coplanar_seg_tri_do_intersect(hb, 2, ha))
			ic += classify_coplanar_edge_intersections(hb, 2, ha, &cec);

		if (ic > 0)
		{
			OMC_ARR_PROF_TTI_INCR(23);
			ts.addCoplanarTriangles(ha.t_id, hb.t_id);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}

		return;
	}

	{ // initialize cache
#ifdef OMC_ARR_DC_FILTER_O3D
		if (!t_nmax_init)
		{
			ha.t_nmax = ts.triPlane(ha.t_id);
			hb.t_nmax = ts.triPlane(hb.t_id);
		}
		if (!v_wrt_seg_init)
		{
			ha.init_v_wrt_seg();
			hb.init_v_wrt_seg();
		}
#else
		ha.t_nmax = ts.triPlane(ha.t_id);
		hb.t_nmax = ts.triPlane(hb.t_id);
		ha.init_v_wrt_seg();
		hb.init_v_wrt_seg();
#endif
		ha.init_v_in_seg();
		ha.init_v_in_tri();
		ha.init_seg_wrt_seg();
		hb.init_v_in_seg();
		hb.init_v_in_tri();
		hb.init_seg_wrt_seg();
	}

	if (_singleCoplanarEdge(orAB, edge_id))
	{
		OMC_ARR_PROF_TTI_INCR(24);
		// CASE: a single coplanar edge of ta is coplanar to tb (e.g. orAB: 0 0 1).
		if (coplanar_seg_tri_do_intersect(ha, edge_id, hb))
		{
			OMC_ARR_PROF_TTI_INCR(25);
			classify_coplanar_edge_intersections(ha, edge_id, hb, nullptr);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeOnSameSide(orAB, vtx_id))
	{
		OMC_ARR_PROF_TTI_INCR(26);
		// CASE: a vertex of ta is coplanar to tb, and the opposite edge is on the
		// same side respect to tb  (e.g. orAB: 1 0 1)
		if (classify_coplanr_vtx_intersections(ha, vtx_id, hb, intersection_points,
		                                       intersection_types))
		{
			OMC_ARR_PROF_TTI_INCR(27);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeCrossPlane(orAB, vtx_id))
	{
		OMC_ARR_PROF_TTI_INCR(28);
		// CASE: a vertex of ta is coplanar to tb, and the opposite edge could
		// intersect tb (e.g. orAB: -1 0 1)
		classify_coplanr_vtx_intersections(ha, vtx_id, hb, intersection_points,
		                                   intersection_types);

		index_t opp_edge_id = (vtx_id + 1) % 3;

		if (noncoplanar_seg_tri_do_intersect(ha, opp_edge_id, hb))
		{
			OMC_ARR_PROF_TTI_INCR(29);
			classify_noncoplanar_edge_intersections(
			  ha, opp_edge_id, hb, intersection_points, intersection_types);
		}
		// go on checking B->A
	}
	else if (_vtxOnASideAndOppositeEdgeOnTheOther(orAB, vtx_id))
	{
		OMC_ARR_PROF_TTI_INCR(30);
		// CASE: a vertex of ta is on one side of the plane defined to tb, and the
		// opposite edge (always in ta) is in the other (e.g. orAB: -1 1 1)
		if (noncoplanar_seg_tri_do_intersect(ha, vtx_id, hb))
		{
			OMC_ARR_PROF_TTI_INCR(31);
			classify_noncoplanar_edge_intersections(
			  ha, vtx_id, hb, intersection_points, intersection_types);
		}
		if (noncoplanar_seg_tri_do_intersect(ha, (vtx_id + 2) % 3, hb))
		{
			OMC_ARR_PROF_TTI_INCR(32);
			classify_noncoplanar_edge_intersections(
			  ha, (vtx_id + 2) % 3, hb, intersection_points, intersection_types);
		}
		// go on checking B->A
	}

	OMC_EXPENSIVE_ASSERT(
	  intersection_points.size() <= 2,
	  "more than 2 intersection points between non-coplanar triangles.");
	if (intersection_points.size() == 2)
	{
		// when detect more than one intersections in noncoplanar cases,
		// we do not continue classifying but propagate undetected intersections in
		// later stages.
		index_t v0_id = intersection_points[0];
		index_t v1_id = intersection_points[1];

		OMC_EXPENSIVE_ASSERT(intersection_types.size() == 2, "mismatch size.");
		// Check if an edge from tB is coplanar to tA and intersects tA.
		if (intersection_on_one_edge(intersection_types, edge_id))
		{
			OMC_EXPENSIVE_ASSERT(coplanar_seg_tri_do_intersect(hb, edge_id, ha),
			                     "coplanar edge does not intersect triangle.");
			ts.addCoplanarEdge(ha.t_id, get_e_id(hb, edge_id), v0_id, v1_id);
		}

		add_symbolic_segment(v0_id, v1_id, ha, InvalidIndex, hb, InvalidIndex);
		ts.setTriangleHasIntersections(ha.t_id, hb.t_id);

		OMC_ARR_PROF_TTI_INCR(33);
		return; // all possible intersections are found, return.
	}
	// intersection_types is not used after here.

	// =========================================================================
	// check of B respect to A
	// =========================================================================

	Sign orBA[3];

	orBA[0] = get_v_wrt_tri(hb, 0, ha);
	orBA[1] = get_v_wrt_tri(hb, 1, ha);
	orBA[2] = get_v_wrt_tri(hb, 2, ha);

	if (_sameOrientation(orBA[0], orBA[1]) &&
	    _sameOrientation(orBA[1], orBA[2]) && (orBA[0] != Sign::ZERO))
	{
		// CASE: no intersection found
		OMC_ARR_PROF_TTI_INCR(34);
		return;
	}

	if (_singleCoplanarEdge(orBA, edge_id))
	{
		OMC_ARR_PROF_TTI_INCR(35);
		// CASE: a single edge of tB is coplanar to tA    (e.g. orBA: 1 0 0)
		if (coplanar_seg_tri_do_intersect(hb, edge_id, ha))
		{
			OMC_ARR_PROF_TTI_INCR(36);
			classify_coplanar_edge_intersections(hb, edge_id, ha, nullptr);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeOnSameSide(orBA, vtx_id))
	{
		OMC_ARR_PROF_TTI_INCR(37);
		// CASE: a vertex of tB is coplanar to tA, and the opposite edge is on the
		// same side respect to tA  (e.g. orBA: 1 0 1)
		if (classify_coplanr_vtx_intersections(hb, vtx_id, ha, intersection_points,
		                                       intersection_types))
		{
			OMC_ARR_PROF_TTI_INCR(38);
			ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
		}
		return; // whether intersections are found or not, return.
	}
	else if (_vtxInPlaneAndOppositeEdgeCrossPlane(orBA, vtx_id))
	{
		OMC_ARR_PROF_TTI_INCR(39);
		// CASE: a vertex of tB is coplanar to tA, and the opposite edge could
		// intersect tA (e.g. orBA: -1 0 1)
		classify_coplanr_vtx_intersections(hb, vtx_id, ha, intersection_points,
		                                   intersection_types);

		index_t opp_edge_id = (vtx_id + 1) % 3;

		if (noncoplanar_seg_tri_do_intersect(hb, opp_edge_id, ha))
		{
			OMC_ARR_PROF_TTI_INCR(40);
			classify_noncoplanar_edge_intersections(
			  hb, opp_edge_id, ha, intersection_points, intersection_types);
		}
	}
	else if (_vtxOnASideAndOppositeEdgeOnTheOther(orBA, vtx_id))
	{
		OMC_ARR_PROF_TTI_INCR(41);
		// CASE: a vertex of tB is on one side of the plane defined to tA, and the
		// opposite edge (always in tB) is in the other (e.g. orBA: -1 1 1)
		if (noncoplanar_seg_tri_do_intersect(hb, vtx_id, ha))
		{
			OMC_ARR_PROF_TTI_INCR(42);
			classify_noncoplanar_edge_intersections(
			  hb, vtx_id, ha, intersection_points, intersection_types);
		}
		if (noncoplanar_seg_tri_do_intersect(hb, (vtx_id + 2) % 3, ha))
		{
			OMC_ARR_PROF_TTI_INCR(43);
			classify_noncoplanar_edge_intersections(
			  hb, (vtx_id + 2) % 3, ha, intersection_points, intersection_types);
		}
	}

	std::sort(intersection_points.begin(), intersection_points.end());
	intersection_points.resize(
	  std::unique(intersection_points.begin(), intersection_points.end()) -
	  intersection_points.begin());

	OMC_EXPENSIVE_ASSERT(
	  intersection_points.size() <= 2,
	  "more than 2 intersection points between non-coplanar triangles");

	if (intersection_points.size() == 2)
	{
		index_t v0_id = intersection_points[0];
		index_t v1_id = intersection_points[1];

		add_symbolic_segment(v0_id, v1_id, ha, InvalidIndex, hb, InvalidIndex);
		ts.setTriangleHasIntersections(ha.t_id, hb.t_id);
	}
}

/**
 * @brief get index of an edge of triangle `ta`
 * @param ha helper for triangle `ta`
 * @param ea an edge of triangle `ta`
 * @return index_t index of the edge `ea`
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::get_e_id(TTIHelper &ha, index_t ea)
{
	if (is_valid_idx(ha.e_id[ea]))
		return ha.e_id[ea];

	ha.e_id[ea] = ts.triEdgeID(ha.t_id, ea, /*addIfNotFound*/ std::true_type());
	return ha.e_id[ea];
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
 * @return bool true if `va` is strictly inside `tb`.
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
	const NT *p  = ha.v[va];
	const NT *t0 = hb.v[0], *t1 = hb.v[1], *t2 = hb.v[2];

	Sign sign = Orient3D()(t0, t1, t2, p);

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
 * The sector is formed by the two adjacent edges of vertex `vb`.
 * `vb` is the opposite vertex to the edge `eb` in triangle `tb`.
 * @details
 *    ----eb----
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
		return 1; // excludes line ebn and ebp
	if (va_wrt_ebn == Sign::NEGATIVE && va_wrt_ebp == Sign::POSITIVE)
		return 2; // excludes line ebp and ebn

	Sign va_wrt_eb = get_v_wrt_seg(ha, va, hb, eb);
	if ((va_wrt_eb >= Sign::ZERO && va_wrt_ebn >= Sign::ZERO &&
	     va_wrt_ebp >= Sign::ZERO) ||
	    (va_wrt_eb <= Sign::ZERO && va_wrt_ebn <= Sign::ZERO &&
	     va_wrt_ebp <= Sign::ZERO))
		return 0; // includes line eb, ebn and ebp

	Sign vb_wrt_eb =
	  OrientOn2D()(hb.v[eb], hb.v[(eb + 1) % 3], hb.v[(eb + 2) % 3], hb.t_nmax);

	return va_wrt_eb != vb_wrt_eb ? 0 : 3;
	// 3 excludes line ebn and ebp
}

/**
 * @brief find vertex correspondences
 * @return size_t how many vertices correspond
 * @note we have removed duplicate vertices in arrangements pipeline, so, we can
 * check correspondences by index instead of geometry.
 */
template <typename Traits>
size_t DetectClassifyTTI<Traits>::find_vtx_correspondence(TTIHelper &ha,
                                                          TTIHelper &hb)
{
	size_t count = 0;
	// clang-format off
	if (ha.v_id[0] == hb.v_id[0]) { ha.v_in_vtx[0] = 0; hb.v_in_vtx[0] = 0; count += 1; }
	if (ha.v_id[0] == hb.v_id[1]) { ha.v_in_vtx[0] = 1; hb.v_in_vtx[1] = 0; count += 1; }
	if (ha.v_id[0] == hb.v_id[2]) { ha.v_in_vtx[0] = 2; hb.v_in_vtx[2] = 0; count += 1; }
	if (ha.v_id[1] == hb.v_id[0]) { ha.v_in_vtx[1] = 0; hb.v_in_vtx[0] = 1; count += 1; }
	if (ha.v_id[1] == hb.v_id[1]) { ha.v_in_vtx[1] = 1; hb.v_in_vtx[1] = 1; count += 1; }
	if (ha.v_id[1] == hb.v_id[2]) { ha.v_in_vtx[1] = 2; hb.v_in_vtx[2] = 1; count += 1; }
	if (ha.v_id[2] == hb.v_id[0]) { ha.v_in_vtx[2] = 0; hb.v_in_vtx[0] = 2; count += 1; }
	if (ha.v_id[2] == hb.v_id[1]) { ha.v_in_vtx[2] = 1; hb.v_in_vtx[1] = 2; count += 1; }
	if (ha.v_id[2] == hb.v_id[2]) { ha.v_in_vtx[2] = 2; hb.v_in_vtx[2] = 2; count += 1; }
	// clang-format on
	return count;
}

/**
 * @brief fast check if there is intersection on a 2d plane when two triangles
 * share a common vertex.
 * @param ha helper for triangle `ta`
 * @param ea the opposite edge to the shared vertex in triangle `ta`
 * @param hb helper for triangle `tb`
 * @param eb the opposite edge to the shared vertex in triangle `tb`
 * @return false if there is definitely no intersection.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::fast_check_on2d_share_vertex(TTIHelper &ha,
                                                             index_t    ea,
                                                             TTIHelper &hb,
                                                             index_t    eb)
{
	index_t eav0 = ea, eav1 = (ea + 1) % 3; // two endpoints of ea
	index_t ebv0 = eb, ebv1 = (eb + 1) % 3; // two endpoints of eb
	ha.init_v_wrt_seg();
	hb.init_v_wrt_seg();
	Sign ori_ta = OrientOn2D()(ha.v[0], ha.v[1], ha.v[2], ha.t_nmax);
	Sign ori_tb = OrientOn2D()(hb.v[0], hb.v[1], hb.v[2], hb.t_nmax);
	// use four lines adajcent to shared vertex to separate two triangles
	Sign eav0_wrt_line, eav1_wrt_line;
	eav0_wrt_line = get_v_wrt_seg(ha, eav0, hb, (eb + 1) % 3);
	eav1_wrt_line = get_v_wrt_seg(ha, eav1, hb, (eb + 1) % 3);
	if (eav0_wrt_line != Sign::ZERO && eav0_wrt_line == eav1_wrt_line &&
	    eav0_wrt_line != ori_tb)
		return false;
	eav0_wrt_line = get_v_wrt_seg(ha, eav0, hb, (eb + 2) % 3);
	eav1_wrt_line = get_v_wrt_seg(ha, eav1, hb, (eb + 2) % 3);
	if (eav0_wrt_line != Sign::ZERO && eav0_wrt_line == eav1_wrt_line &&
	    eav0_wrt_line != ori_tb)
		return false;
	Sign ebv0_wrt_line, ebv1_wrt_line;
	ebv0_wrt_line = get_v_wrt_seg(hb, ebv0, ha, (ea + 1) % 3);
	ebv1_wrt_line = get_v_wrt_seg(hb, ebv1, ha, (ea + 1) % 3);
	if (ebv0_wrt_line != Sign::ZERO && ebv0_wrt_line == ebv1_wrt_line &&
	    ebv0_wrt_line != ori_ta)
		return false;
	ebv0_wrt_line = get_v_wrt_seg(hb, ebv0, ha, (ea + 2) % 3);
	ebv1_wrt_line = get_v_wrt_seg(hb, ebv1, ha, (ea + 2) % 3);
	if (ebv0_wrt_line != Sign::ZERO && ebv0_wrt_line == ebv1_wrt_line &&
	    ebv0_wrt_line != ori_ta)
		return false;

	return true;
}

/**
 * @brief fast check if there is intersection on a 2d plane when two triangles
 * are separated.
 * @param ha helper for triangle `ta`
 * @param hb helper for triangle `tb`
 * @return false if there is definitely no intersection.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::fast_check_on2d_separate(TTIHelper &ha,
                                                         TTIHelper &hb)
{
	ha.init_v_wrt_seg();
	hb.init_v_wrt_seg();
	Sign ori_ta = OrientOn2D()(ha.v[0], ha.v[1], ha.v[2], ha.t_nmax);
	Sign ori_tb = OrientOn2D()(hb.v[0], hb.v[1], hb.v[2], hb.t_nmax);
	// use four lines adajcent to shared vertex to separate two triangles
	Sign eav0_wrt_line, eav1_wrt_line, eav2_wrt_line;
	eav0_wrt_line = get_v_wrt_seg(ha, /*vb*/ 0, hb, /*eb*/ 0);
	eav1_wrt_line = get_v_wrt_seg(ha, /*vb*/ 1, hb, /*eb*/ 0);
	eav2_wrt_line = get_v_wrt_seg(ha, /*vb*/ 2, hb, /*eb*/ 0);
	if (eav0_wrt_line != Sign::ZERO && eav0_wrt_line == eav1_wrt_line &&
	    eav1_wrt_line == eav2_wrt_line && eav0_wrt_line != ori_tb)
		return false;
	eav0_wrt_line = get_v_wrt_seg(ha, /*vb*/ 0, hb, /*eb*/ 1);
	eav1_wrt_line = get_v_wrt_seg(ha, /*vb*/ 1, hb, /*eb*/ 1);
	eav2_wrt_line = get_v_wrt_seg(ha, /*vb*/ 2, hb, /*eb*/ 1);
	if (eav0_wrt_line != Sign::ZERO && eav0_wrt_line == eav1_wrt_line &&
	    eav1_wrt_line == eav2_wrt_line && eav0_wrt_line != ori_tb)
		return false;
	eav0_wrt_line = get_v_wrt_seg(ha, /*vb*/ 0, hb, /*eb*/ 2);
	eav1_wrt_line = get_v_wrt_seg(ha, /*vb*/ 1, hb, /*eb*/ 2);
	eav2_wrt_line = get_v_wrt_seg(ha, /*vb*/ 2, hb, /*eb*/ 2);
	if (eav0_wrt_line != Sign::ZERO && eav0_wrt_line == eav1_wrt_line &&
	    eav1_wrt_line == eav2_wrt_line && eav0_wrt_line != ori_tb)
		return false;
	Sign ebv0_wrt_line, ebv1_wrt_line, ebv2_wrt_line;
	ebv0_wrt_line = get_v_wrt_seg(hb, /*vb*/ 0, ha, /*eb*/ 0);
	ebv1_wrt_line = get_v_wrt_seg(hb, /*vb*/ 1, ha, /*eb*/ 0);
	ebv2_wrt_line = get_v_wrt_seg(hb, /*vb*/ 2, ha, /*eb*/ 0);
	if (ebv0_wrt_line != Sign::ZERO && ebv0_wrt_line == ebv1_wrt_line &&
	    ebv1_wrt_line == ebv2_wrt_line && ebv0_wrt_line != ori_ta)
		return false;
	ebv0_wrt_line = get_v_wrt_seg(hb, /*vb*/ 0, ha, /*eb*/ 1);
	ebv1_wrt_line = get_v_wrt_seg(hb, /*vb*/ 1, ha, /*eb*/ 1);
	ebv2_wrt_line = get_v_wrt_seg(hb, /*vb*/ 2, ha, /*eb*/ 1);
	if (ebv0_wrt_line != Sign::ZERO && ebv0_wrt_line == ebv1_wrt_line &&
	    ebv1_wrt_line == ebv2_wrt_line && ebv0_wrt_line != ori_ta)
		return false;
	ebv0_wrt_line = get_v_wrt_seg(hb, /*vb*/ 0, ha, /*eb*/ 2);
	ebv1_wrt_line = get_v_wrt_seg(hb, /*vb*/ 1, ha, /*eb*/ 2);
	ebv2_wrt_line = get_v_wrt_seg(hb, /*vb*/ 2, ha, /*eb*/ 2);
	if (ebv0_wrt_line != Sign::ZERO && ebv0_wrt_line == ebv1_wrt_line &&
	    ebv1_wrt_line == ebv2_wrt_line && ebv0_wrt_line != ori_ta)
		return false;

	return true;
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

	OMC_EXPENSIVE_ASSERT_AUX_CODE(index_t eb0 = eb; index_t eb1 = (eb + 1) % 3;)

	Sign ea0_wrt_eb = get_v_wrt_seg(ha, ea0, hb, eb);
	Sign ea1_wrt_eb = get_v_wrt_seg(ha, ea1, hb, eb);

	if (ea0_wrt_eb != ea1_wrt_eb && eb0_wrt_ea != eb1_wrt_ea)
	{
		// two edges sharing an endpoint are impossible in classfication pipeline.
		OMC_EXPENSIVE_ASSERT(!(ha.v_in_vtx[ea0] == eb0 || ha.v_in_vtx[ea0] == eb1 ||
		                       ha.v_in_vtx[ea1] == eb0 || ha.v_in_vtx[ea1] == eb1),
		                     "impossible case.");

		// intersect at a point
		return true;
	}

	// degenerate case: collinear segments.
	// such case is impossible in classification pipeline.
	OMC_EXPENSIVE_ASSERT(!(ea0_wrt_eb == Sign::ZERO && ea1_wrt_eb == Sign::ZERO &&
	                       eb0_wrt_ea == Sign::ZERO && eb1_wrt_ea == Sign::ZERO),
	                     "impossible case.");
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
	OMC_EXPENSIVE_ASSERT(tb_ori != Sign::ZERO, "degenerate triangle");

	Sign eav0_wrt_eb0 = get_v_wrt_seg(ha, eav0, hb, 0);
	Sign eav1_wrt_eb0 = get_v_wrt_seg(ha, eav1, hb, 0);

	if (eav0_wrt_eb0 != Sign::ZERO && eav1_wrt_eb0 != Sign::ZERO &&
	    eav0_wrt_eb0 != tb_ori && eav1_wrt_eb0 != tb_ori)
		return false; // ea and tb are separated by eb0

	Sign eav0_wrt_eb1 = get_v_wrt_seg(ha, eav0, hb, 1);
	Sign eav1_wrt_eb1 = get_v_wrt_seg(ha, eav1, hb, 1);

	if (eav0_wrt_eb1 != Sign::ZERO && eav1_wrt_eb1 != Sign::ZERO &&
	    eav0_wrt_eb1 != tb_ori && eav1_wrt_eb1 != tb_ori)
		return false; // ea and tb are separated by eb1

	Sign eav0_wrt_eb2 = get_v_wrt_seg(ha, eav0, hb, 2);
	Sign eav1_wrt_eb2 = get_v_wrt_seg(ha, eav1, hb, 2);

	if (eav0_wrt_eb2 != Sign::ZERO && eav1_wrt_eb2 != Sign::ZERO &&
	    eav0_wrt_eb2 != tb_ori && eav1_wrt_eb2 != tb_ori)
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
 * @param intersection_points list to store intersection points.
 * @return PointInSimplexType return where the intersection point is in `tb`.
 * @note this func does not assume that `va` and `tb` intersect, it will check
 * and then classify.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::classify_coplanr_vtx_intersections(
  TTIHelper &ha, index_t va, TTIHelper &hb,
  IntersectionPoints &intersection_points,
  IntersectionTypes  &intersection_types)
{
	OMC_EXPENSIVE_ASSERT(!is_valid_idx(ha.v_in_vtx[va]), "impossible case");

	if (TTIHelper::is_cached(ha.v_in_seg[va]) && is_valid_idx(ha.v_in_seg[va]))
	{ // strictly inside an edge of `tb`
		// if v_in_seg is set, the vertex must have been added to the edge.
		// add_vertex_in_edge(hb, ha.v_in_seg[va], ha, va);
		intersection_points.push_back(ha.v_id[va]);
		intersection_types.push_back(ha.v_in_seg[va] == 0
		                               ? PointInSimplexType::ON_EDGE0
		                               : (ha.v_in_seg[va] == 1
		                                    ? PointInSimplexType::ON_EDGE1
		                                    : PointInSimplexType::ON_EDGE2));
		return true;
	}

	Sign ori_tb = OrientOn2D()(hb.v[0], hb.v[1], hb.v[2], hb.t_nmax);

	for (size_t i = 0; i < 3; i++)
	{
		Sign va_wrt_eb0 = get_v_wrt_seg(ha, va, hb, i);
		if (va_wrt_eb0 != Sign::ZERO && va_wrt_eb0 != ori_tb)
			return false;
		if (va_wrt_eb0 == Sign::ZERO)
		{ // on the support line of ebi
			if (get_v_in_seg(ha, va, hb, i))
			{ // strictly inside ebi
				ha.v_in_seg[va] = i;
				intersection_points.push_back(add_vertex_in_edge(hb, i, ha, va));
				intersection_types.push_back(
				  i == 0 ? PointInSimplexType::ON_EDGE0
				         : (i == 1 ? PointInSimplexType::ON_EDGE1
				                   : PointInSimplexType::ON_EDGE2));
				return true;
			}
			else
				return false;
		}
	}

	// strictly inside `tb`
	add_vertex_in_tri(hb, ha, va);
	intersection_points.push_back(ha.v_id[va]);
	intersection_types.push_back(PointInSimplexType::STRICTLY_INSIDE);
	return true;
}

/**
 * @param ha helper for `ta`.
 * @param ea A valid local edge index in `ta`.
 * @param hb helper for `tb`.
 * @param copl_edge_crosses list to store coplanar intersection points.
 * @return whether intersection points are detected (share vertex is excluded).
 * @note this func does not assume that `ea` and `tb` intersect, it will
 * check and then classify.
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::classify_coplanar_edge_intersections(
  TTIHelper &ha, index_t ea, TTIHelper &hb, CoplanarEEIList *copl_edge_crosses)
{
	index_t            ev0 = ea, ev1 = (ea + 1) % 3;
	IntersectionPoints ip; // intersection points, excluding shared vertices

	// check if `ev0` and `ev1` on any vertex of `tb`
	index_t ev0_in_vtx = ha.v_in_vtx[ev0];
	index_t ev1_in_vtx = ha.v_in_vtx[ev1];

	if (is_valid_idx(ev0_in_vtx) && is_valid_idx(ev1_in_vtx))
		return 0; // Although this case should not happen in arrangements,
		          // we keep this check for completeness.

	// check if `ev0` and `ev1` on any edge of `tb`
	index_t ev0_in_seg = get_v_in_seg(ha, ev0, hb);
	index_t ev1_in_seg = get_v_in_seg(ha, ev1, hb);

	if (is_valid_idx(ev0_in_seg))
		ip.push_back(add_vertex_in_edge(hb, ev0_in_seg, ha, ev0));
	if (is_valid_idx(ev1_in_seg))
		ip.push_back(add_vertex_in_edge(hb, ev1_in_seg, ha, ev1));

	if (is_valid_idx(ev0_in_seg) && is_valid_idx(ev1_in_seg))
	{
		if (ev0_in_seg != ev1_in_seg) // ea crosses tb
		{
			add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], ha.v_id[ev1]);
		}
		else // ea is totally contained in an edge of tb
		{
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, ev0_in_seg),
			                   ha.v_id[ev0], ha.v_id[ev1]);
		}
		return true;
	}
	if (is_valid_idx(ev0_in_seg) && is_valid_idx(ev1_in_vtx))
	{
		if (ev0_in_seg == (ev1_in_vtx + 1) % 3) // ea crosses tb
		{
			add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], ha.v_id[ev1]);
		}
		else // ea is totally contained in an edge of tb
		{
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, ev0_in_seg),
			                   ha.v_id[ev0], ha.v_id[ev1]);
		}
		return true;
	}
	if (is_valid_idx(ev1_in_seg) && is_valid_idx(ev0_in_vtx))
	{
		if (ev1_in_seg == (ev0_in_vtx + 1) % 3) // ea crosses tb
		{
			add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], ha.v_id[ev1]);
		}
		else // ea is totally contained in an edge of tb
		{
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, ev1_in_seg),
			                   ha.v_id[ev0], ha.v_id[ev1]);
		}
		return true;
	}

	// check if `ev0` and `ev1` inside `tb`
	bool ev1_in_tri = get_v_in_tri(ha, ev1, hb);
	if (ev1_in_tri)
		ip.push_back(add_vertex_in_tri(hb, ha, ev1));

	if ((is_valid_idx(ev0_in_seg) || is_valid_idx(ev0_in_vtx)) && ev1_in_tri)
	{ // v0 in a segment or vtx and v1 inside triangle
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha, ea, hb, InvalidIndex);
		ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], ha.v_id[ev1]);
		return true;
	}

	bool ev0_in_tri = get_v_in_tri(ha, ev0, hb);
	if (ev0_in_tri)
		ip.push_back(add_vertex_in_tri(hb, ha, ev0));

	if ((is_valid_idx(ev1_in_seg) || is_valid_idx(ev1_in_vtx)) && ev0_in_tri)
	{ // v1 in a segment or vtx and v0 inside triangle
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha, ea, hb, InvalidIndex);
		ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], ha.v_id[ev1]);
		return true;
	}

	if (ev0_in_tri && ev1_in_tri)
	{ // v0 and v1 are both inside the triangle, ea is totally inside tb
		add_symbolic_segment(ha.v_id[ev0], ha.v_id[ev1], ha, ea, hb, InvalidIndex);
		ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], ha.v_id[ev1]);
		return true;
	}

	// Check any vertex of `tb` is inside `ea`

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
				{ // `vb` is strictly inside `ea`
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
				{ // `vb` is strictly inside `ea`
					hb.v_wrt_seg[vb][ea] = Sign::ZERO;
					hb.v_in_seg[vb]      = ea;
					return std::pair<bool, Sign>(true, vb_wrt_ea);
				}
			}
			return std::pair<bool, Sign>(false, vb_wrt_ea);
		}
	};

	auto [vb0_in_ea, vb0_wrt_ea] = is_vb_in_ea(/*vb*/ 0);
	if (vb0_in_ea)
	{
		ip.push_back(add_vertex_in_edge(ha, ea, hb, /*vb*/ 0));

		if (ev0_in_seg == /*eb*/ 1 || ev0_in_tri)
		{ // ea crosses tb (also crosses vb0)
			add_symbolic_segment(hb.v_id[0], ha.v_id[ev0], hb, InvalidIndex, ha, ea);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], hb.v_id[0]);
			return true;
		}
		else if (ev0_in_seg == /*eb*/ 0 || ev0_in_vtx == /*vb*/ 1)
		{ // ea and eb0 overlap between vb0 and ev0
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 0), hb.v_id[0],
			                   ha.v_id[ev0]);
			return true;
		}
		else if (ev0_in_seg == /*eb*/ 2 || ev0_in_vtx == /*vb*/ 2)
		{ // ea and eb2 overlap between vb0 and ev0
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 2), hb.v_id[0],
			                   ha.v_id[ev0]);
			return true;
		}

		if (ev1_in_seg == /*eb*/ 1 || ev1_in_tri)
		{ // ea crosses tb (also crosses vb0)
			add_symbolic_segment(hb.v_id[0], ha.v_id[ev1], hb, InvalidIndex, ha, ea);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev1], hb.v_id[0]);
			return true;
		}
		else if (ev1_in_seg == /*eb*/ 0 || ev1_in_vtx == /*vb*/ 1)
		{ // ea and eb0 overlap between vb0 and ev1
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 0), hb.v_id[0],
			                   ha.v_id[ev1]);
			return true;
		}
		else if (ev1_in_seg == /*eb*/ 2 || ev1_in_vtx == /*vb*/ 2)
		{ // ea and eb2 overlap between vb0 and ev1
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 2), hb.v_id[0],
			                   ha.v_id[ev1]);
			return true;
		}
	}

	auto [vb1_in_ea, vb1_wrt_ea] = is_vb_in_ea(/*vb*/ 1);
	if (vb1_in_ea)
	{
		ip.push_back(add_vertex_in_edge(ha, ea, hb, /*vb*/ 1));

		if (ev0_in_seg == /*eb*/ 2 || ev0_in_tri)
		{ // ea crosses tb (also crosses vb1)
			add_symbolic_segment(hb.v_id[1], ha.v_id[ev0], hb, InvalidIndex, ha, ea);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], hb.v_id[1]);
			return true;
		}
		else if (ev0_in_seg == /*eb*/ 0 || ev0_in_vtx == /*vb*/ 0)
		{ // ea and eb0 overlap between vb1 and ev0
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 0), hb.v_id[1],
			                   ha.v_id[ev0]);
			return true;
		}
		else if (ev0_in_seg == /*eb*/ 1 || ev0_in_vtx == /*vb*/ 2)
		{ // ea and eb1 overlap between vb1 and ev0
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 1), hb.v_id[1],
			                   ha.v_id[ev0]);
			return true;
		}

		if (ev1_in_seg == /*eb*/ 2 || ev1_in_tri)
		{ // ea crosses tb (also crosses vb1)
			add_symbolic_segment(hb.v_id[1], ha.v_id[ev1], hb, InvalidIndex, ha, ea);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev1], hb.v_id[1]);
			return true;
		}
		else if (ev1_in_seg == /*eb*/ 0 || ev1_in_vtx == /*vb*/ 0)
		{ // ea and eb0 overlap between vb1 and ev1
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 0), hb.v_id[1],
			                   ha.v_id[ev1]);
			return true;
		}
		else if (ev1_in_seg == /*eb*/ 1 || ev1_in_vtx == /*vb*/ 2)
		{ // ea and eb1 overlap between vb1 and ev1
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 1), hb.v_id[1],
			                   ha.v_id[ev1]);
			return true;
		}
	}

	auto [vb2_in_ea, vb2_wrt_ea] = is_vb_in_ea(/*vb*/ 2);
	if (vb2_in_ea)
	{
		ip.push_back(add_vertex_in_edge(ha, ea, hb, /*vb*/ 2));

		if (ev0_in_seg == /*eb*/ 0 || ev0_in_tri)
		{ // ea crosses tb (also crosses vb2)
			add_symbolic_segment(hb.v_id[2], ha.v_id[ev0], hb, InvalidIndex, ha, ea);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], hb.v_id[2]);
			return true;
		}
		else if (ev0_in_seg == /*eb*/ 1 || ev0_in_vtx == /*vb*/ 1)
		{ // ea and eb1 overlap between vb2 and ev0
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 1), hb.v_id[2],
			                   ha.v_id[ev0]);
			return true;
		}
		else if (ev0_in_seg == /*eb*/ 2 || ev0_in_vtx == /*vb*/ 0)
		{ // ea and eb2 overlap between vb0 and ev0
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 2), hb.v_id[2],
			                   ha.v_id[ev0]);
			return true;
		}

		if (ev1_in_seg == /*eb*/ 0 || ev1_in_tri)
		{ // ea crosses tb (also crosses vb2)
			add_symbolic_segment(hb.v_id[2], ha.v_id[ev1], hb, InvalidIndex, ha, ea);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev1], hb.v_id[2]);
			return true;
		}
		else if (ev1_in_seg == /*eb*/ 1 || ev1_in_vtx == /*vb*/ 1)
		{ // ea and eb1 overlap between vb2 and ev1
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 1), hb.v_id[2],
			                   ha.v_id[ev1]);
			return true;
		}
		else if (ev1_in_seg == /*eb*/ 2 || ev1_in_vtx == /*vb*/ 0)
		{ // ea and eb2 overlap between vb0 and ev1
			ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, /*eb*/ 2), hb.v_id[2],
			                   ha.v_id[ev1]);
			return true;
		}
	}

	if (vb0_in_ea && vb1_in_ea)
	{ // eb0 is totally contained in ea
		ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, 0), hb.v_id[0],
		                   hb.v_id[1]);
		return true;
	}
	if (vb1_in_ea && vb2_in_ea)
	{ // eb1 is totally contained in ea
		ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, 1), hb.v_id[1],
		                   hb.v_id[2]);
		return true;
	}
	if (vb2_in_ea && vb0_in_ea)
	{ // eb2 is totally contained in ea
		ts.addColinearEdge(get_e_id(ha, ea), get_e_id(hb, 2), hb.v_id[2],
		                   hb.v_id[0]);
		return true;
	}

	// Edges cross check

	index_t seg0_cross = InvalidIndex;
	index_t seg1_cross = InvalidIndex;
	index_t seg2_cross = InvalidIndex;

	if (/* endpoints of ea are not on eb0 or eb0's endpoints */
	    ev0_in_vtx != /*vb*/ 0 && ev0_in_vtx != /*vb*/ 1 &&
	    ev1_in_vtx != /*vb*/ 0 && ev1_in_vtx != /*vb*/ 1 &&
	    ev0_in_seg != /*eb*/ 0 && ev1_in_seg != /*eb*/ 0 &&
	    /* endpoints of eb0 are not on ea */
	    !vb0_in_ea && !vb1_in_ea &&
	    /* ea and eb0 cross*/
	    seg_seg_do_intersect(ha, ea, hb, /*eb*/ 0, vb0_wrt_ea, vb1_wrt_ea))
	{
		// edge `ea` cross seg 0 in `tb`
		seg0_cross =
		  add_edge_cross_coplanar_edge(ha, ea, hb, /*eb*/ 0, copl_edge_crosses);
		ip.push_back(seg0_cross);

		if (is_valid_idx(ev0_in_vtx) || is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{ // cross eb0 and tb
			add_symbolic_segment(ha.v_id[ev0], seg0_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], seg0_cross);
			return true;
		}
		if (is_valid_idx(ev1_in_vtx) || is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{ // cross eb0 and tb
			add_symbolic_segment(ha.v_id[ev1], seg0_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev1], seg0_cross);
			return true;
		}
		if (vb2_in_ea)
		{ // cross eb0 and (opposite vertex) vb2
			add_symbolic_segment(hb.v_id[/*vb*/ 2], seg0_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), hb.v_id[2], seg0_cross);
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
		seg1_cross =
		  add_edge_cross_coplanar_edge(ha, ea, hb, /*eb*/ 1, copl_edge_crosses);
		ip.push_back(seg1_cross);

		if (is_valid_idx(ev0_in_vtx) || is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{ // cross eb1 and tb
			add_symbolic_segment(ha.v_id[ev0], seg1_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], seg1_cross);
			return true;
		}
		if (is_valid_idx(ev1_in_vtx) || is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{ // cross eb1 and tb
			add_symbolic_segment(ha.v_id[ev1], seg1_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev1], seg1_cross);
			return true;
		}
		if (vb0_in_ea)
		{ // cross eb1 and (opposite vertex) vb0
			add_symbolic_segment(hb.v_id[/*vb*/ 0], seg1_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), hb.v_id[0], seg1_cross);
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
		seg2_cross =
		  add_edge_cross_coplanar_edge(ha, ea, hb, /*eb*/ 2, copl_edge_crosses);
		ip.push_back(seg2_cross);

		if (is_valid_idx(ev0_in_vtx) || is_valid_idx(ev0_in_seg) || ev0_in_tri)
		{ // cross eb2 and tb
			add_symbolic_segment(ha.v_id[ev0], seg2_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev0], seg2_cross);
			return true;
		}
		if (is_valid_idx(ev1_in_vtx) || is_valid_idx(ev1_in_seg) || ev1_in_tri)
		{ // cross eb2 and tb
			add_symbolic_segment(ha.v_id[ev1], seg2_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), ha.v_id[ev1], seg2_cross);
			return true;
		}
		if (vb1_in_ea)
		{ // cross eb2 and (opposite vertex) vb1
			add_symbolic_segment(hb.v_id[/*vb*/ 1], seg2_cross, ha, ea, hb, InvalidIndex);
			ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), hb.v_id[1], seg2_cross);
			return true;
		}
	}

	// final probably symbolic edges: `ea` cross two edges of `tb`
	if (is_valid_idx(seg0_cross) && is_valid_idx(seg1_cross))
	{
		add_symbolic_segment(seg0_cross, seg1_cross, ha, ea, hb, InvalidIndex);
		ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), seg0_cross, seg1_cross);
		return true;
	}
	if (is_valid_idx(seg0_cross) && is_valid_idx(seg2_cross))
	{
		add_symbolic_segment(seg0_cross, seg2_cross, ha, ea, hb, InvalidIndex);
		ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), seg0_cross, seg2_cross);
		return true;
	}
	if (is_valid_idx(seg1_cross) && is_valid_idx(seg2_cross))
	{
		add_symbolic_segment(seg1_cross, seg2_cross, ha, ea, hb, InvalidIndex);
		ts.addCoplanarEdge(hb.t_id, get_e_id(ha, ea), seg1_cross, seg2_cross);
		return true;
	}

	// there remain two cases that cause intersection point but are neither
	// crossing triangle (add_symblic_segment) nor being colinear with an edge
	// (addColinearEdge).
	// They are:
	// 1. only ev0 or ev1 is inside one edge of tb (ev0_in_seg, ev1_in_seg)
	// 2. one vertex of tb is inside of ea (vb0_in_ea, vb1_in_ea, vb2_in_ea)
	// Intersection points are found in the two cases and are put in list.
	OMC_EXPENSIVE_ASSERT(ip.size() <= 1, "impossible case");
	return !ip.empty();
}

/**
 * @param ha helper for `ta`.
 * @param ea A valid local edge index in `ta`.
 * @param hb helper for `tb`.
 * @param intersection_points list to store intersection points.
 * @return PointInSimplexType return where the intersection point is in `tb`.
 * @note this func assumes that `ea` and `tb` definitely intersect at a point
 * except the vertices of `tb`
 */
template <typename Traits>
bool DetectClassifyTTI<Traits>::classify_noncoplanar_edge_intersections(
  TTIHelper &ha, index_t ea, TTIHelper &hb,
  IntersectionPoints &intersection_points,
  IntersectionTypes  &intersection_types)
{
	for (size_t i = 0; i < 3; i++)
	{
		if (hb.v_in_seg[i] == ea)
		{
			// if v_in_seg is set, the vertex must have been added to the edge.
			// add_vertex_in_edge(ha, ea, hb, /*vb*/ i);
			intersection_points.push_back(hb.v_id[i]);
			intersection_types.push_back(i == 0
			                               ? PointInSimplexType::ON_VERT0
			                               : (i == 1 ? PointInSimplexType::ON_VERT1
			                                         : PointInSimplexType::ON_VERT2));
			return true;
		}
	}

	Sign vol_ea_eb0 = get_seg_wrt_seg(ha, ea, hb, 0);
	Sign vol_ea_eb2 = get_seg_wrt_seg(ha, ea, hb, 2);

	// `ea` cross vb0
	if (vol_ea_eb0 == Sign::ZERO && vol_ea_eb2 == Sign::ZERO)
	{
		hb.v_in_seg[/*vb*/ 0] = ea;
		intersection_points.push_back(add_vertex_in_edge(ha, ea, hb, /*vb*/ 0));
		intersection_types.push_back(PointInSimplexType::ON_VERT0);
		return true;
	}
	Sign vol_ea_eb1 = get_seg_wrt_seg(ha, ea, hb, 1);
	// `ea` cross vb1
	if (vol_ea_eb0 == Sign::ZERO && vol_ea_eb1 == Sign::ZERO)
	{
		hb.v_in_seg[/*vb*/ 1] = ea;
		intersection_points.push_back(add_vertex_in_edge(ha, ea, hb, /*vb*/ 1));
		intersection_types.push_back(PointInSimplexType::ON_VERT1);
		return true;
	}
	// `ea` cross vb2
	if (vol_ea_eb1 == Sign::ZERO && vol_ea_eb2 == Sign::ZERO)
	{
		hb.v_in_seg[/*vb*/ 2] = ea;
		intersection_points.push_back(add_vertex_in_edge(ha, ea, hb, /*vb*/ 2));
		intersection_types.push_back(PointInSimplexType::ON_VERT2);
		return true;
	}
	// the edge intersects the tri in seg 0
	if (vol_ea_eb0 == Sign::ZERO)
	{
		intersection_points.push_back(
		  add_edge_cross_noncoplanar_edge(ha, ea, hb, /*eb*/ 0));
		intersection_types.push_back(PointInSimplexType::ON_EDGE0);
		return true;
	}
	// the edge intersects the tri in seg 1
	if (vol_ea_eb1 == Sign::ZERO)
	{
		intersection_points.push_back(
		  add_edge_cross_noncoplanar_edge(ha, ea, hb, /*eb*/ 1));
		intersection_types.push_back(PointInSimplexType::ON_EDGE1);
		return true;
	}
	// the edge intersects the tri in seg 2
	if (vol_ea_eb2 == Sign::ZERO)
	{
		intersection_points.push_back(
		  add_edge_cross_noncoplanar_edge(ha, ea, hb, /*eb*/ 2));
		intersection_types.push_back(PointInSimplexType::ON_EDGE2);
		return true;
	}

	// the edge intersects the inner triangle
	intersection_points.push_back(add_edge_cross_tri(ha, ea, hb));
	intersection_types.push_back(PointInSimplexType::STRICTLY_INSIDE);
	return true;
}

/**
 * @param v0 **GLOBAL** index of a vertex
 * @param v1 **GLOBAL** index of a vertex
 * @param ha helper for a triangle
 * @param ea ea is a valid index if segment comes from ea
 * @param hb helper for another triangle
 * @param eb eb is a valid index if segment comes from eb
 * @details symbolic segment is symmetric.
 */
template <typename Traits>
void DetectClassifyTTI<Traits>::add_symbolic_segment(index_t v0, index_t v1,
                                                     TTIHelper &ha, index_t ea,
                                                     TTIHelper &hb, index_t eb)
{
	OMC_EXPENSIVE_ASSERT(v0 != v1, "trying to add a 0-length symbolic edge");
	UIPair  seg = unique_pair(v0, v1);
	index_t seg_id;
	if (is_valid_idx(ea))
		seg_id = ts.getOrAddSegment(seg, get_e_id(ha, ea));
	else if (is_valid_idx(eb))
		seg_id = ts.getOrAddSegment(seg, get_e_id(hb, eb));
	else
		seg_id = ts.getOrAddSegment(seg, InvalidIndex);

	// check if (v0,v1) is an edge of ta
	// if (!(ha.contains_vtx(v0) && ha.contains_vtx(v1)))
	// check if (v0,v1) is an edge of ta
	if (!is_valid_idx(ea))
		// if not, add segment to ta
		ts.addSegmentInTriangle(ha.t_id, seg_id);

	// check if (v0,v1) is an edge of tb
	// if (!(hb.contains_vtx(v0) && hb.contains_vtx(v1)))
	// if not, add segment to tb
	if (!is_valid_idx(eb))
		ts.addSegmentInTriangle(hb.t_id, seg_id);

	// seg2tris will be constructed after classification.
}

/**
 * @param ha helper for triangle `ta`
 * @param hb helper for triangle `tb`
 * @param vb **LOCAL** index of a vertex of triangle `tb`
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_vertex_in_tri(TTIHelper &ha,
                                                     TTIHelper &hb, index_t vb)
{
	ts.addVertexInTriangle(ha.t_id, hb.v_id[vb]);
	return hb.v_id[vb];
}

/**
 * @param ha helper for triangle `ta`
 * @param ea **LOCAL** index of an edge of triangle `ta`
 * @param hb helper for triangle `tb`
 * @param vb **LOCAL** index of a vertex of triangle `tb`
 * @return index_t
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_vertex_in_edge(TTIHelper &ha, index_t ea,
                                                      TTIHelper &hb, index_t vb)
{
	index_t ea_id = get_e_id(ha, ea);
	index_t vb_id = hb.v_id[vb];

#ifndef OMC_ARR_GLOBAL_POINT_SET
	// get the mutex
	tbb::spin_mutex &ea_mutex = ts.getE2PMutex(ea_id);

	// lock until add or fix end.
	std::lock_guard<tbb::spin_mutex> lock(ea_mutex);

	index_t found_vid = ts.findVertexInEdge(ea_id, ts.vert(vb_id));
	if (is_valid_idx(found_vid))
	{
		if (ts.vert(found_vid).point_type() != GPoint::PointType::Explicit
		    /* || vb_id < found_vid */)
		{ // note: two explicit point won't be same in arrangements pipeline.
			// the only case is that a complicate implicit point is replaced by a
			// simple explicit point.
			ts.fixVertexInEdge(ea_id, /*old*/ found_vid, /*new*/ vb_id);
			// return vb_id;
		}
		// else
		//	return found_vid;
	}
	else
	{
		ts.addVertexInEdge(ea_id, vb_id);
		// return vb_id;
	}
#else
	ts.addVertexInEdge(ea_id, vb_id);
#endif
	return vb_id; // vb_id will always be added.
}

/**
 * @param ha helper for triangle `ta`
 * @param ea **LOCAL** index of an edge of triangle `ta`
 * @param hb helper for triangle `tb`
 * @param eb **LOCAL** index of an edge of triangle `tb`
 * @param copl_edge_crosses coplanar intersection points
 * @return index_t global index of the new implicit point
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_edge_cross_coplanar_edge(
  TTIHelper &ha, index_t ea, TTIHelper &hb, index_t eb,
  CoplanarEEIList *copl_edge_crosses)
{
	index_t ea_id = get_e_id(ha, ea);
	index_t eb_id = get_e_id(hb, eb);

	// find intersection point in existed points
	if (copl_edge_crosses)
	{
		for (const CoplanarEEI &c : *copl_edge_crosses)
			if (c.is_same(ea_id, eb_id))
				return c.p;
	}
	// otherwise find it or create a new point

#ifdef OMC_ARR_AUX_LPI
	// find a non-coplanar jolly point
	auto find_jolly = [this, &hb]()
	{
		const NT *v0 = hb.v[0], *v1 = hb.v[1], *v2 = hb.v[2];
		for (index_t i = 0; i < 5; i++)
		{
			const NT *jp = ts.jolly_points[i]->EXP().data();
			if (Orient3D()(v0, v1, v2, jp) != Sign::ZERO)
				return i;
		}
		OMC_ASSERT(false, "impossible.");
		return index_t(4);
	};
	index_t jp_idx = find_jolly();

	IPoint_LPI *new_v = pnt_arena.emplace(CreateLPI()(
	  ts.vert(ha.v_id[ea]), ts.vert(ha.v_id[(ea + 1) % 3]), ts.vert(hb.v_id[eb]),
	  ts.vert(hb.v_id[(eb + 1) % 3]), *ts.jolly_points[jp_idx]));

	// try to add the point
	auto [added_vid, new_vertex_created] = add_SSI(ea_id, eb_id, new_v);

	if (!new_vertex_created)
	{
		IPoint_LPI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
	}
#else
	// fix the vertex sequential of created implicit point
	// index_t min_eid, max_eid;
	index_t min_min_vid, min_max_vid, max_min_vid, max_max_vid;
	if (ea_id < eb_id)
	{
		// min_eid     = ea_id;
		min_min_vid = std::min(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
		min_max_vid = std::max(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
		// max_eid     = eb_id;
		max_min_vid = std::min(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
		max_max_vid = std::max(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
	}
	else
	{
		// min_eid     = eb_id;
		min_min_vid = std::min(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
		min_max_vid = std::max(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
		// max_eid     = ea_id;
		max_min_vid = std::min(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
		max_max_vid = std::max(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
	}

	IPoint_SSI *new_v = pnt_arena.emplace(CreateSSI()(
	  ts.vert(min_min_vid), ts.vert(min_max_vid), ts.vert(max_min_vid),
	  ts.vert(max_max_vid), int_to_OrPlane(hb.t_nmax)));

	// try to add the point
	auto [added_vid, new_vertex_created] = add_SSI(ea_id, eb_id, new_v);

	if (!new_vertex_created)
	{
		IPoint_SSI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
	}
#endif

	OMC_EXPENSIVE_ASSERT(is_valid_idx(added_vid), "invalid index");
	if (copl_edge_crosses)
		copl_edge_crosses->push_back(CoplanarEEI(ea_id, eb_id, added_vid));
	return added_vid;
}

/**
 * @brief edge `ea` of triangle `ta` crosses edge `eb` of triangle `tb`
 * @param ha helper for triangle `ta`
 * @param ea **LOCAL** index of an edge
 * @param hb helper for triangle `tb`
 * @param eb **LOCAL** index of an edge
 * @return global index of the intersection point
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_edge_cross_noncoplanar_edge(
  TTIHelper &ha, index_t ea, TTIHelper &hb, index_t eb)
{
#ifdef OMC_ARR_AUX_LPI
	// fix the vertex sequential of created implicit point
	index_t e_min_vid = std::min(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
	index_t e_max_vid = std::max(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);

	IPoint_LPI *new_v = pnt_arena.emplace(
	  CreateLPI()(ts.vert(e_min_vid), ts.vert(e_max_vid), ts.vert(hb.v_id[0]),
	              ts.vert(hb.v_id[1]), ts.vert(hb.v_id[2])));

	// try to add the point
	auto [added_vid, new_vertex_created] =
	  add_SSI(get_e_id(ha, ea), get_e_id(hb, eb), new_v);

	if (!new_vertex_created)
	{
		IPoint_LPI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
	}
#else
	// find a plane for two interseted segments
	int plane = -1;
	{
		static std::array<index_t, 12> _tri = {0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};

		std::array<const NT *, 4> _pnts = {ha.v[ea], ha.v[(ea + 1) % 3], hb.v[eb],
		                                   hb.v[(eb + 1) % 3]};

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

	// fix the vertex sequential of created implicit point
	index_t ea_id = get_e_id(ha, ea);
	index_t eb_id = get_e_id(hb, eb);
	// index_t min_eid, max_eid;
	index_t min_min_vid, min_max_vid, max_min_vid, max_max_vid;
	if (ea_id < eb_id)
	{
		// min_eid     = ea_id;
		min_min_vid = std::min(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
		min_max_vid = std::max(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
		// max_eid     = eb_id;
		max_min_vid = std::min(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
		max_max_vid = std::max(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
	}
	else
	{
		// min_eid     = eb_id;
		min_min_vid = std::min(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
		min_max_vid = std::max(hb.v_id[eb], hb.v_id[(eb + 1) % 3]);
		// max_eid     = ea_id;
		max_min_vid = std::min(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
		max_max_vid = std::max(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
	}

	IPoint_SSI *new_v = pnt_arena.emplace(
	  CreateSSI()(ts.vert(min_min_vid), ts.vert(min_max_vid),
	              ts.vert(max_min_vid), ts.vert(max_max_vid), plane));

	// try to add the point
	auto [added_vid, new_vertex_created] = add_SSI(ea_id, eb_id, new_v);

	if (!new_vertex_created)
	{
		IPoint_SSI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
	}
#endif

	OMC_EXPENSIVE_ASSERT(is_valid_idx(added_vid), "invalid index");
	return added_vid;
}

/**
 * @brief an edge `ea` of triangle `ta` crosses triangle `tb`
 * @param ha helper for triangle `ta`
 * @param ea **LOCAL** index of an edge of triangle `ta`
 * @param hb helper for triangle `tb`
 * @return index_t global index of the intersection point
 */
template <typename Traits>
index_t DetectClassifyTTI<Traits>::add_edge_cross_tri(TTIHelper &ha, index_t ea,
                                                      TTIHelper &hb)
{
	// fix the vertex sequential of created implicit point
	index_t e_min_vid = std::min(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);
	index_t e_max_vid = std::max(ha.v_id[ea], ha.v_id[(ea + 1) % 3]);

	IPoint_LPI *new_v = pnt_arena.emplace(
	  CreateLPI()(ts.vert(e_min_vid), ts.vert(e_max_vid), ts.vert(hb.v_id[0]),
	              ts.vert(hb.v_id[1]), ts.vert(hb.v_id[2])));

	// try to add the point
	auto [added_vid, new_vertex_created] =
	  add_LPI(get_e_id(ha, ea), hb.t_id, new_v);

	if (!new_vertex_created)
	{
		IPoint_LPI::gcv().remove(new_v);
		pnt_arena.recycle(new_v);
	}

	OMC_EXPENSIVE_ASSERT(is_valid_idx(added_vid), "invalid index");

	return added_vid;
}

/**
 * @brief Add a new SSI point.
 * @param ea_id **GLOBAL** index of an edge
 * @param eb_id **GLOBAL** index of an edge
 * @param new_v the new point
 * @return std::pair<index_t, bool> The first index is the index of the result
 * point (existed point or newly added point). The second boolean is true if
 * it succeeds to add the point, otherwise it fails because there is already an
 * existed point.
 */
template <typename Traits>
std::pair<index_t, bool>
DetectClassifyTTI<Traits>::add_SSI(index_t ea_id, index_t eb_id, GPoint *new_v)
{
#ifndef OMC_ARR_GLOBAL_POINT_SET
	// get mutexes of edges.
	index_t          min_eid   = std::min(ea_id, eb_id);
	index_t          max_eid   = std::max(ea_id, eb_id);
	tbb::spin_mutex &min_mutex = ts.getE2PMutex(min_eid);
	tbb::spin_mutex &max_mutex = ts.getE2PMutex(max_eid);

	// lock until add or fix end.
	std::lock_guard<tbb::spin_mutex> min_lock(min_mutex);
	std::lock_guard<tbb::spin_mutex> max_lock(max_mutex);

	index_t found_vid_in_ea = ts.findVertexInEdge(ea_id, *new_v);
	index_t found_vid_in_eb = ts.findVertexInEdge(eb_id, *new_v);

	bool    new_vertex_created = false;
	index_t suitable_vid = InvalidIndex, fix_vid = InvalidIndex;

	if (is_valid_idx(found_vid_in_ea) && is_valid_idx(found_vid_in_eb))
	{
		// get the suitable point with simpler type and smaller index
		if (found_vid_in_ea == found_vid_in_eb)
		{
			suitable_vid = found_vid_in_ea;
		}
		else
		{
			if (ts.vert(found_vid_in_ea).point_type() <
			    ts.vert(found_vid_in_eb).point_type())
			{
				suitable_vid = found_vid_in_ea;
				fix_vid      = found_vid_in_eb;
			}
			else if (ts.vert(found_vid_in_ea).point_type() >
			         ts.vert(found_vid_in_eb).point_type())
			{
				suitable_vid = found_vid_in_eb;
				fix_vid      = found_vid_in_ea;
			}
			else // same type
			{
				suitable_vid = std::min(found_vid_in_ea, found_vid_in_eb);
				fix_vid      = std::max(found_vid_in_ea, found_vid_in_eb);
			}
			index_t fix_edge = fix_vid == found_vid_in_ea ? ea_id : eb_id;
			ts.fixVertexInEdge(/*edge*/ fix_edge, /*old*/ fix_vid,
			                   /*new*/ suitable_vid);
		}
	}
	else if (is_valid_idx(found_vid_in_ea))
	{
		// get the suitable point with simpler type and smaller index
		if (ts.vert(found_vid_in_ea).point_type() <= new_v->point_type())
		{
			suitable_vid = found_vid_in_ea;
			ts.addVertexInEdge(eb_id, found_vid_in_ea);
		}
		else // ts.vert(found_vid_in_ea).point_type() > new_v->point_type()
		{
			suitable_vid = ts.addImplVert(new_v);
			fix_vid      = found_vid_in_ea;
			ts.addVertexInEdge(eb_id, suitable_vid);
			ts.fixVertexInEdge(/*edge*/ ea_id, /*old*/ fix_vid,
			                   /*new*/ suitable_vid);
			new_vertex_created = true;
		}
	}
	else if (is_valid_idx(found_vid_in_eb))
	{
		// get the suitable point with simpler type and smaller index
		if (ts.vert(found_vid_in_eb).point_type() <= new_v->point_type())
		{
			suitable_vid = found_vid_in_eb;
			ts.addVertexInEdge(ea_id, found_vid_in_eb);
		}
		else // ts.vert(found_vid_in_eb).point_type() > new_v->point_type()
		{
			suitable_vid = ts.addImplVert(new_v);
			fix_vid      = found_vid_in_eb;
			ts.addVertexInEdge(ea_id, suitable_vid);
			ts.fixVertexInEdge(/*edge*/ eb_id, /*old*/ fix_vid,
			                   /*new*/ suitable_vid);
			new_vertex_created = true;
		}
	}
	else
	{
		suitable_vid = ts.addImplVert(new_v);
		ts.addVertexInEdge(ea_id, suitable_vid);
		ts.addVertexInEdge(eb_id, suitable_vid);
		new_vertex_created = true;
	}

	return std::pair<index_t, bool>(suitable_vid, new_vertex_created);
#else
	auto [added_vid, new_vertex_created] = ts.addUniquePoint(*new_v);
	ts.addVertexInEdge(ea_id, added_vid);
	ts.addVertexInEdge(eb_id, added_vid);

	return std::pair<index_t, bool>(added_vid, new_vertex_created);
#endif
}

/**
 * @brief add an LPI point
 * @param e_id **GLOBAL** index of an edge
 * @param t_id **GLOBAL** index of a triangle
 * @param new_v the new point
 * @return std::pair<index_t, bool> The first index is the index of the result
 * point (existed point or newly added point). The second boolean is true if
 * it succeeds to add the point, otherwise it fails because there is already an
 * existed point.
 */
template <typename Traits>
std::pair<index_t, bool>
DetectClassifyTTI<Traits>::add_LPI(index_t e_id, index_t t_id, GPoint *new_v)
{
#ifndef OMC_ARR_GLOBAL_POINT_SET
	// get mutexes of edges.
	tbb::spin_mutex &mutex = ts.getE2PMutex(e_id);

	// lock until add or fix end.
	std::lock_guard<tbb::spin_mutex> lock(mutex);

	index_t found_vid = ts.findVertexInEdge(e_id, *new_v);

	bool    new_vertex_created = false;
	index_t suitable_vid = InvalidIndex, fix_vid = InvalidIndex;

	if (is_valid_idx(found_vid))
	{
		// get the suitable point with simpler type and smaller index
		if (ts.vert(found_vid).point_type() <= new_v->point_type())
		{
			suitable_vid = found_vid;
		}
		else // ts.vert(found_vid).point_type() > new_v->point_type()
		{
			suitable_vid = ts.addImplVert(new_v);
			fix_vid      = found_vid;
			ts.fixVertexInEdge(/*edge*/ e_id, /*old*/ fix_vid,
			                   /*new*/ suitable_vid);
			new_vertex_created = true;
		}
	}
	else // not found
	{
		suitable_vid = ts.addImplVert(new_v);
		ts.addVertexInEdge(e_id, suitable_vid);
		new_vertex_created = true;
	}

	ts.addVertexInTriangle(t_id, suitable_vid);

	return std::pair<index_t, bool>(suitable_vid, new_vertex_created);
#else
	auto [added_vid, new_vertex_created] = ts.addUniquePoint(*new_v);
	ts.addVertexInEdge(e_id, added_vid);
	ts.addVertexInTriangle(t_id, added_vid);

	return std::pair<index_t, bool>(added_vid, new_vertex_created);
#endif
}

} // namespace OMC