#pragma once

#include "ClassifyIntersections.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"
#include "OpenMeshCraft/Utils/Exception.h"

// #define OMC_ARR_CI_NO_CACHE

namespace OMC {

template <typename Traits>
ClassifyIntersections<Traits>::ClassifyIntersections(
  TriSoup &_ts, AuxStruct &_g, const std::vector<ShewchukCachePtr> &_cache,
  bool _parallel, MeshArrangements_Stats *_stats, bool _verbose)
  : ts(_ts)
  , pnt_arenas(*ts.pnt_arenas)
  , idx_arenas(*ts.idx_arenas)
  , uniq_g(_g)
  , cache(_cache)
  , parallel(_parallel)
  , verbose(_verbose)
  , stats(_stats)
{
	if (parallel)
	{
		concurrent_g =
		  std::vector<AuxStruct>(tbb::this_task_arena::max_concurrency());
		for (AuxStruct &g : concurrent_g)
			g.initialize(ts);

		thread_contexts.reserve(tbb::this_task_arena::max_concurrency());
		for (size_t i = 0; i < (size_t)tbb::this_task_arena::max_concurrency(); i++)
		{
			thread_contexts.emplace_back(pnt_arenas[i], idx_arenas[i],
			                             concurrent_g[i]);
		}
	}
	else
	{
		thread_contexts.emplace_back(pnt_arenas[0], idx_arenas[0], uniq_g);
	}
}

template <typename Traits>
void ClassifyIntersections<Traits>::checkTriangleTriangleIntersections()
{
	auto check_tt = [this](const UIPair &pair)
	{
		GPoint::clear_global_cached_values();

		index_t tA_id = pair.first, tB_id = pair.second;

		uniq_g.setTriangleHasIntersections(tA_id);
		uniq_g.setTriangleHasIntersections(tB_id);

		if (parallel)
		{ // we store results in coccurent aux struct
			int thread_id = tbb::this_task_arena::current_thread_index();
			checkTriTriInter(tA_id, tB_id, thread_contexts[thread_id]);
		}
		else
		{ // we store results in unique aux struct
			checkTriTriInter(tA_id, tB_id, thread_contexts[0]);
		}
	};

	GPoint::enable_global_cached_values(tbb::this_task_arena::max_concurrency());

	if (parallel)
		tbb::parallel_for_each(uniq_g.intersection_list.begin(),
		                       uniq_g.intersection_list.end(), check_tt);
	else
		std::for_each(uniq_g.intersection_list.begin(),
		              uniq_g.intersection_list.end(), check_tt);

	if (parallel) // merge results stored in concurrent aux struct
		mergeConcurrentAuxStructures();

	sortEdgePointsList();

	GPoint::disable_global_cached_values();
}

inline bool sameOrientation(const Sign o1, const Sign o2) { return o1 == o2; }

// true if all edges are coplanar to the triangle, return false otherwise
inline bool allCoplanarEdges(const Sign o[])
{
	return (o[0] == Sign::ZERO && o[1] == Sign::ZERO && o[2] == Sign::ZERO);
}

// if there is a coplanar edge return true and set edge_id, return false
// otherwise
inline bool singleCoplanarEdge(const Sign o[], index_t &edge_id)
{
	if (o[0] == Sign::ZERO && o[1] == Sign::ZERO && o[2] != Sign::ZERO)
	{
		edge_id = 0;
		return true;
	}
	if (o[1] == Sign::ZERO && o[2] == Sign::ZERO && o[0] != Sign::ZERO)
	{
		edge_id = 1;
		return true;
	}
	if (o[2] == Sign::ZERO && o[0] == Sign::ZERO && o[1] != Sign::ZERO)
	{
		edge_id = 2;
		return true;
	}
	edge_id = InvalidIndex;
	return false; // false
}

// if there is a vertex in the plane and the opposite edge doesn't intersect the
// plane return true and set the vtx id, return false otherwise
inline bool vtxInPlaneAndOppositeEdgeOnSameSide(const Sign o[], index_t &vtx_id)
{
	if (o[0] == Sign::ZERO && o[1] == o[2] && o[1] != Sign::ZERO)
	{
		vtx_id = 0;
		return true;
	}
	if (o[1] == Sign::ZERO && o[0] == o[2] && o[0] != Sign::ZERO)
	{
		vtx_id = 1;
		return true;
	}
	if (o[2] == Sign::ZERO && o[0] == o[1] && o[0] != Sign::ZERO)
	{
		vtx_id = 2;
		return true;
	}
	vtx_id = InvalidIndex;
	return false;
}

// if there is a vertex in the plane and the opposite edge intersect the plane
// return true and set the vtx id, return false otherwise
inline bool vtxInPlaneAndOppositeEdgeCrossPlane(const Sign o[], index_t &vtx_id)
{
	if (o[0] == Sign::ZERO && o[1] != o[2] && o[1] != Sign::ZERO &&
	    o[2] != Sign::ZERO)
	{
		vtx_id = 0;
		return true;
	}
	if (o[1] == Sign::ZERO && o[0] != o[2] && o[0] != Sign::ZERO &&
	    o[2] != Sign::ZERO)
	{
		vtx_id = 1;
		return true;
	}
	if (o[2] == Sign::ZERO && o[0] != o[1] && o[0] != Sign::ZERO &&
	    o[1] != Sign::ZERO)
	{
		vtx_id = 2;
		return true;
	}
	vtx_id = InvalidIndex;
	return false;
}

// if there is a vertex on one side and the opposite edge on the other return
// the relative informations, -1 otherwise
inline bool vtxOnASideAndOppositeEdgeOnTheOther(const Sign o[], index_t &vtx_id,
                                                index_t &opp_v0,
                                                index_t &opp_v1)
{
	if (o[0] == Sign::ZERO || o[1] == Sign::ZERO || o[2] == Sign::ZERO)
	{
		vtx_id = InvalidIndex;
		return false; // one vtx on the plane
	}

	if (o[0] == o[1] && o[1] == o[2])
	{
		vtx_id = InvalidIndex;
		return false; // all vtx on the same side of the plane
	}

	if (o[0] == o[1])
	{
		opp_v0 = 0;
		opp_v1 = 1;
		vtx_id = 2;
		return true;
	}

	if (o[0] == o[2])
	{
		opp_v0 = 0;
		opp_v1 = 2;
		vtx_id = 1;
		return true;
	}

	opp_v0 = 1;
	opp_v1 = 2;
	vtx_id = 0;
	return true;
}

template <typename Traits>
void ClassifyIntersections<Traits>::checkTriTriInter(index_t        tA_id,
                                                     index_t        tB_id,
                                                     ThreadContext &tc)
{
	// temporary vtx list for final symbolic edge creation
	phmap::flat_hash_set<size_t> vtx_list;
	// intersection list
	phmap::flat_hash_set<size_t> inter_list;

	bool coplanar_tris = false;

	/* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 *      check of tB respect to tA
	 * ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	Sign orBA[3];
#ifdef OMC_ARR_CI_NO_CACHE
	ShewchukCache tA_cache;
	Orient3D::get_minors(ts.triVertPtr(tA_id, 0), ts.triVertPtr(tA_id, 1),
	                     ts.triVertPtr(tA_id, 2), tA_cache.minor, tA_cache.perm);
#else
	const ShewchukCache &tA_cache = *cache[tA_id];
#endif
	orBA[0] = Orient3D::with_cached_minors(
	  ts.triVertPtr(tA_id, 0), ts.triVertPtr(tA_id, 1), ts.triVertPtr(tA_id, 2),
	  ts.triVertPtr(tB_id, 0), tA_cache.minor, tA_cache.perm);
	orBA[1] = Orient3D::with_cached_minors(
	  ts.triVertPtr(tA_id, 0), ts.triVertPtr(tA_id, 1), ts.triVertPtr(tA_id, 2),
	  ts.triVertPtr(tB_id, 1), tA_cache.minor, tA_cache.perm);
	orBA[2] = Orient3D::with_cached_minors(
	  ts.triVertPtr(tA_id, 0), ts.triVertPtr(tA_id, 1), ts.triVertPtr(tA_id, 2),
	  ts.triVertPtr(tB_id, 2), tA_cache.minor, tA_cache.perm);

	index_t edge_id, vtx_id, opp_v0, opp_v1;

	if (sameOrientation(orBA[0], orBA[1]) && sameOrientation(orBA[1], orBA[2]) &&
	    (orBA[0] != Sign::ZERO))
	{
		// CASE: no intersection found
		return;
	}
	else if (allCoplanarEdges(orBA))
	{
		// CASE: all edge of tB are coplanar to all edges of tA   (orBA: 0 0 0)
		tc.g.addCoplanarTriangles(tA_id, tB_id);
		coplanar_tris = true;

		checkSingleCoplanarEdgeIntersections(ts.triVertID(tB_id, 0),
		                                     ts.triVertID(tB_id, 1), tB_id, tA_id,
		                                     inter_list, tc);
		checkSingleCoplanarEdgeIntersections(ts.triVertID(tB_id, 1),
		                                     ts.triVertID(tB_id, 2), tB_id, tA_id,
		                                     inter_list, tc);
		checkSingleCoplanarEdgeIntersections(ts.triVertID(tB_id, 2),
		                                     ts.triVertID(tB_id, 0), tB_id, tA_id,
		                                     inter_list, tc);
	}
	else if (singleCoplanarEdge(orBA, edge_id))
	{
		// CASE: a single edge of tB is coplanar to tA    (e.g. orBA: 1 0 0)
		size_t e_v0_id = edge_id;
		size_t e_v1_id = (edge_id + 1) % 3;
		if (checkSingleCoplanarEdgeIntersections(ts.triVertID(tB_id, e_v0_id),
		                                         ts.triVertID(tB_id, e_v1_id),
		                                         tB_id, tA_id, inter_list, tc))
			tc.g.addCoplanarEdge(tA_id, ts.triEdgeID(tB_id, edge_id));
	}
	else if (vtxInPlaneAndOppositeEdgeOnSameSide(orBA, vtx_id))
	{
		// CASE: a vertex of tB is coplanar to tA, and the opposite edge is on the
		// same side respect to tA  (e.g. orBA: 1 0 1)
		checkVtxInTriangleIntersection(ts.triVertID(tB_id, vtx_id), tA_id, vtx_list,
		                               inter_list, tc);
	}
	else if (vtxInPlaneAndOppositeEdgeCrossPlane(orBA, vtx_id))
	{
		// CASE: a vertex of tB is coplanar to tA, and the opposite edge could
		// intersect tA (e.g. orBA: -1 0 1)
		index_t real_v_id = ts.triVertID(tB_id, vtx_id);
		checkVtxInTriangleIntersection(real_v_id, tA_id, vtx_list, inter_list, tc);

		index_t opp_edge_id =
		  ts.edgeOppositeToVert(tB_id, ts.triVertID(tB_id, vtx_id));
		checkSingleNoCoplanarEdgeIntersection(opp_edge_id, tA_id, vtx_list,
		                                      inter_list, tc);
	}
	else if (vtxOnASideAndOppositeEdgeOnTheOther(orBA, vtx_id, opp_v0, opp_v1))
	{
		// CASE: a vertex of tB is on one side of the plane defined to tA, and the
		// opposite edge (always in tB) is in the other (e.g. orBA: -1 1 1)
		index_t id_v      = ts.triVertID(tB_id, vtx_id);
		index_t id_opp_v0 = ts.triVertID(tB_id, opp_v0);
		index_t id_opp_v1 = ts.triVertID(tB_id, opp_v1);

		index_t edge_id0 = ts.edgeID(id_v, id_opp_v0);
		index_t edge_id1 = ts.edgeID(id_v, id_opp_v1);
		OMC_EXPENSIVE_ASSERT(is_valid_idx(edge_id0) && is_valid_idx(edge_id1),
		                     "edge invalid");

		checkSingleNoCoplanarEdgeIntersection(edge_id0, tA_id, vtx_list, inter_list,
		                                      tc);
		checkSingleNoCoplanarEdgeIntersection(edge_id1, tA_id, vtx_list, inter_list,
		                                      tc);
	}

	/* ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 *      check of A respect to B
	 * ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 */

	Sign orAB[3];

	// all edge of tA are coplanar to all edges of tB   (orAB: 0 0 0)
	if (coplanar_tris)
	{
		orAB[0] = Sign::ZERO;
		orAB[1] = Sign::ZERO;
		orAB[2] = Sign::ZERO;
	}
	else
	{
#ifdef OMC_ARR_CI_NO_CACHE
		ShewchukCache tB_cache;
		Orient3D::get_minors(ts.triVertPtr(tB_id, 0), ts.triVertPtr(tB_id, 1),
		                     ts.triVertPtr(tB_id, 2), tB_cache.minor,
		                     tB_cache.perm);
#else
		const ShewchukCache &tB_cache = *cache[tB_id];
#endif

		orAB[0] = Orient3D::with_cached_minors(
		  ts.triVertPtr(tB_id, 0), ts.triVertPtr(tB_id, 1), ts.triVertPtr(tB_id, 2),
		  ts.triVertPtr(tA_id, 0), tB_cache.minor, tB_cache.perm);
		orAB[1] = Orient3D::with_cached_minors(
		  ts.triVertPtr(tB_id, 0), ts.triVertPtr(tB_id, 1), ts.triVertPtr(tB_id, 2),
		  ts.triVertPtr(tA_id, 1), tB_cache.minor, tB_cache.perm);
		orAB[2] = Orient3D::with_cached_minors(
		  ts.triVertPtr(tB_id, 0), ts.triVertPtr(tB_id, 1), ts.triVertPtr(tB_id, 2),
		  ts.triVertPtr(tA_id, 2), tB_cache.minor, tB_cache.perm);
	}

	if (!coplanar_tris && inter_list.size() > 1)
	{
		// before goto final_check, check if an edge from tA is coplanar to tB and
		// intersects tB.
		if (singleCoplanarEdge(orAB, edge_id))
		{
			index_t e_v0_id  = edge_id;
			index_t e_v1_id  = (edge_id + 1) % 3;
			int     tB_plane = static_cast<int>(ts.triPlane(tB_id));
			if (Triangle3_Segment3_DoIntersect().intersection_type(
			      ts.triVertPtr(tB_id, 0), ts.triVertPtr(tB_id, 1),
			      ts.triVertPtr(tB_id, 2), ts.triVertPtr(tA_id, e_v0_id),
			      ts.triVertPtr(tA_id, e_v1_id), tB_plane, cache[tB_id]->minor,
			      cache[tB_id]->perm) >= SimplexIntersectionType::INTERSECT)
			{
				tc.g.addCoplanarEdge(tB_id, ts.triEdgeID(tA_id, edge_id));
			}
		}
		goto final_check; // sorry about that :(		// may god forgive you XD
	}

	if (sameOrientation(orAB[0], orAB[1]) && sameOrientation(orAB[1], orAB[2]) &&
	    (orAB[0] != Sign::ZERO))
	{
		// CASE: no intersection
		return;
	}
	else if (coplanar_tris)
	{
		// CASE: all edge of tA are coplanar to all edges of tB   (orAB: 0 0 0)
		checkSingleCoplanarEdgeIntersections(ts.triVertID(tA_id, 0),
		                                     ts.triVertID(tA_id, 1), tA_id, tB_id,
		                                     inter_list, tc);
		checkSingleCoplanarEdgeIntersections(ts.triVertID(tA_id, 1),
		                                     ts.triVertID(tA_id, 2), tA_id, tB_id,
		                                     inter_list, tc);
		checkSingleCoplanarEdgeIntersections(ts.triVertID(tA_id, 2),
		                                     ts.triVertID(tA_id, 0), tA_id, tB_id,
		                                     inter_list, tc);
	}
	else if (singleCoplanarEdge(orAB, edge_id))
	{
		// CASE: a single edge of tA is coplanar to tB    (e.g. orAB: 1 0 0)
		index_t e_v0_id = edge_id;
		index_t e_v1_id = (edge_id + 1) % 3;
		if (checkSingleCoplanarEdgeIntersections(ts.triVertID(tA_id, e_v0_id),
		                                         ts.triVertID(tA_id, e_v1_id),
		                                         tA_id, tB_id, inter_list, tc))
			tc.g.addCoplanarEdge(tB_id, ts.triEdgeID(tA_id, edge_id));
	}
	else if (vtxInPlaneAndOppositeEdgeOnSameSide(orAB, vtx_id))
	{
		// CASE: a vertex of tA is coplanar to tB, and the opposite edge is on the
		// same side respect to tB  (e.g. orAB: 1 0 1)
		checkVtxInTriangleIntersection(ts.triVertID(tA_id, vtx_id), tB_id, vtx_list,
		                               inter_list, tc);
	}
	else if (vtxInPlaneAndOppositeEdgeCrossPlane(orAB, vtx_id))
	{
		// CASE: a vertex of tA is coplanar to tB, and the opposite edge could
		// intersect tB (e.g. orAB: -1 0 1)
		index_t real_v_id = ts.triVertID(tA_id, vtx_id);
		checkVtxInTriangleIntersection(real_v_id, tB_id, vtx_list, inter_list, tc);

		index_t opp_edge_id =
		  ts.edgeOppositeToVert(tA_id, ts.triVertID(tA_id, vtx_id));
		checkSingleNoCoplanarEdgeIntersection(opp_edge_id, tB_id, vtx_list,
		                                      inter_list, tc);
	}
	else if (vtxOnASideAndOppositeEdgeOnTheOther(orAB, vtx_id, opp_v0, opp_v1))
	{
		// CASE: a vertex of tA is on one side of the plane defined to tB, and the
		// opposite edge (always in tA) is in the other (e.g. orBA: -1 1 1)
		index_t id_v      = ts.triVertID(tA_id, vtx_id);
		index_t id_opp_v0 = ts.triVertID(tA_id, opp_v0);
		index_t id_opp_v1 = ts.triVertID(tA_id, opp_v1);

		index_t edge_id0 = ts.edgeID(id_v, id_opp_v0);
		index_t edge_id1 = ts.edgeID(id_v, id_opp_v1);
		OMC_EXPENSIVE_ASSERT(is_valid_idx(edge_id0) && is_valid_idx(edge_id1),
		                     "edge invalid");

		checkSingleNoCoplanarEdgeIntersection(edge_id0, tB_id, vtx_list, inter_list,
		                                      tc);
		checkSingleNoCoplanarEdgeIntersection(edge_id1, tB_id, vtx_list, inter_list,
		                                      tc);
	}

final_check:

	if (coplanar_tris)
	{
		OMC_EXPENSIVE_ASSERT(
		  vtx_list.size() <= 3,
		  "more than 3 intersection points in coplanar triangles");
	}
	else
	{
		OMC_EXPENSIVE_ASSERT(
		  (!coplanar_tris && vtx_list.size() <= 2),
		  "more than 2 intersection points in 2 no-coplanar traingles");
	}

	if (vtx_list.size() == 2)
	{
		index_t v0_id = *(vtx_list.begin());
		index_t v1_id = *(++vtx_list.begin());

		addSymbolicSegment(v0_id, v1_id, tA_id, tB_id, tc);
	}
}

/**
 * @brief One edge ( \p ev0, \p e_v1 ) in triangle \p e_t_id is coplanar to
 * one triangle \p o_t_id, check intersection between them.
 * @param li intersection list
 * @return true if coplanar edge intersect triangle.
 */
template <typename Traits>
bool ClassifyIntersections<Traits>::checkSingleCoplanarEdgeIntersections(
  index_t e_v0, index_t e_v1, index_t e_t_id, index_t o_t_id,
  phmap::flat_hash_set<size_t> &li, ThreadContext &tc)
{
	bool    v0_in_vtx = false, v1_in_vtx = false;
	index_t v0_in_seg = InvalidIndex, v1_in_seg = InvalidIndex;
	bool    v0_in_tri = false, v1_in_tri = false;

	// e_v0 position
	PointInSimplexType v0_inters =
	  Triangle3_Point3_DoIntersect().intersection_type(
	    ts.triVertPtr(o_t_id, 0), ts.triVertPtr(o_t_id, 1),
	    ts.triVertPtr(o_t_id, 2), ts.vertPtr(e_v0), ts.triPlane(o_t_id));

	if (v0_inters == PointInSimplexType::ON_VERT0 ||
	    v0_inters == PointInSimplexType::ON_VERT1 ||
	    v0_inters == PointInSimplexType::ON_VERT2)
	{
		v0_in_vtx = true; // v0 in a vertex
		li.insert(e_v0);
	}
	else if (v0_inters == PointInSimplexType::ON_EDGE0)
		v0_in_seg = ts.triEdgeID(o_t_id, 0); // v0 in seg0
	else if (v0_inters == PointInSimplexType::ON_EDGE1)
		v0_in_seg = ts.triEdgeID(o_t_id, 1); // v0 in seg1
	else if (v0_inters == PointInSimplexType::ON_EDGE2)
		v0_in_seg = ts.triEdgeID(o_t_id, 2); // v0 in seg2
	else if (v0_inters == PointInSimplexType::STRICTLY_INSIDE)
		v0_in_tri = true; // v0 inside tri

	// e_v1 position
	PointInSimplexType v1_inters =
	  Triangle3_Point3_DoIntersect().intersection_type(
	    ts.triVertPtr(o_t_id, 0), ts.triVertPtr(o_t_id, 1),
	    ts.triVertPtr(o_t_id, 2), ts.vertPtr(e_v1), ts.triPlane(o_t_id));

	if (v1_inters == PointInSimplexType::ON_VERT0 ||
	    v1_inters == PointInSimplexType::ON_VERT1 ||
	    v1_inters == PointInSimplexType::ON_VERT2)
	{
		v1_in_vtx = true; // v1 in a vertex
		li.insert(e_v1);
	}
	else if (v1_inters == PointInSimplexType::ON_EDGE0)
		v1_in_seg = ts.triEdgeID(o_t_id, 0); // v1 in seg0
	else if (v1_inters == PointInSimplexType::ON_EDGE1)
		v1_in_seg = ts.triEdgeID(o_t_id, 1); // v1 in seg1
	else if (v1_inters == PointInSimplexType::ON_EDGE2)
		v1_in_seg = ts.triEdgeID(o_t_id, 2); // v1 in seg2
	else if (v1_inters == PointInSimplexType::STRICTLY_INSIDE)
		v1_in_tri = true; // v1 inside tri

	if (v0_in_vtx && v1_in_vtx)
		return false;

	if (is_valid_idx(v0_in_seg) && is_valid_idx(v1_in_seg))
	// edge in triangle composed by the link of two vtx in edge
	{
		tc.g.addVertexInEdge(v0_in_seg, e_v0);
		tc.g.addVertexInEdge(v1_in_seg, e_v1);
		li.insert(e_v0);
		li.insert(e_v1);

		addSymbolicSegment(e_v0, e_v1, e_t_id, o_t_id, tc);
		return true;
	}
	else if (is_valid_idx(v0_in_seg)) // only v0 is in a segment of T
	{
		tc.g.addVertexInEdge(v0_in_seg, e_v0);
		li.insert(e_v0);

		if (v1_in_vtx)
		{
			addSymbolicSegment(e_v0, e_v1, e_t_id, o_t_id, tc);
			return true;
		}
	}
	else if (is_valid_idx(v1_in_seg)) // only v1 is in a segment of T
	{
		tc.g.addVertexInEdge(v1_in_seg, e_v1);
		li.insert(e_v1);

		if (v0_in_vtx)
		{
			addSymbolicSegment(e_v1, e_v0, e_t_id, o_t_id, tc);
			return true;
		}
	}

	// v0 in a segment or vtx and v1 inside triangle
	if ((is_valid_idx(v0_in_seg) || v0_in_vtx) && v1_in_tri)
	{
		tc.g.addVertexInTriangle(o_t_id, e_v1);
		li.insert(e_v1);

		addSymbolicSegment(e_v0, e_v1, e_t_id, o_t_id, tc);
		return true;
	}

	// v1 in a segment or vtx and v0 inside triangle
	if ((is_valid_idx(v1_in_seg) || v1_in_vtx) && v0_in_tri)
	{
		tc.g.addVertexInTriangle(o_t_id, e_v0);
		li.insert(e_v0);

		addSymbolicSegment(e_v0, e_v1, e_t_id, o_t_id, tc);
		return true;
	}

	// v0 and v1 both inside the triangle
	if (v0_in_tri && v1_in_tri)
	{
		tc.g.addVertexInTriangle(o_t_id, e_v0);
		tc.g.addVertexInTriangle(o_t_id, e_v1);
		li.insert(e_v0);
		li.insert(e_v1);

		addSymbolicSegment(e_v0, e_v1, e_t_id, o_t_id, tc);
		return true;
	}

	if (v0_in_tri) // only v0 inside the triangle
	{
		tc.g.addVertexInTriangle(o_t_id, e_v0);
		li.insert(e_v0);
	}
	else if (v1_in_tri) // only v1 inside the triangle
	{
		tc.g.addVertexInTriangle(o_t_id, e_v1);
		li.insert(e_v1);
	}

	// Edges cross checking

	// we check only if seg A cross seg B and not B cross A (we found the
	// intersection once)
	index_t o_t_e0 = ts.triEdgeID(o_t_id, 0);
	index_t o_t_e1 = ts.triEdgeID(o_t_id, 1);
	index_t o_t_e2 = ts.triEdgeID(o_t_id, 2);

	bool tv0_in_edge =
	  Segment3_Point3_DoIntersect().intersection_type(
	    ts.vertPtr(e_v0), ts.vertPtr(e_v1), ts.triVertPtr(o_t_id, 0),
	    ts.triPlane(o_t_id)) != PointInSimplexType::STRICTLY_OUTSIDE;
	bool tv1_in_edge =
	  Segment3_Point3_DoIntersect().intersection_type(
	    ts.vertPtr(e_v0), ts.vertPtr(e_v1), ts.triVertPtr(o_t_id, 1),
	    ts.triPlane(o_t_id)) != PointInSimplexType::STRICTLY_OUTSIDE;
	bool tv2_in_edge =
	  Segment3_Point3_DoIntersect().intersection_type(
	    ts.vertPtr(e_v0), ts.vertPtr(e_v1), ts.triVertPtr(o_t_id, 2),
	    ts.triPlane(o_t_id)) != PointInSimplexType::STRICTLY_OUTSIDE;

	index_t seg0_cross = InvalidIndex, seg1_cross = InvalidIndex,
	        seg2_cross = InvalidIndex;
	index_t curr_e_id  = ts.edgeID(e_v0, e_v1);

	if (v0_in_seg != o_t_e0 && v1_in_seg != o_t_e0 && !tv0_in_edge &&
	    !tv1_in_edge &&
	    Segment3_Segment3_DoIntersect().intersection_type(
	      ts.vertPtr(e_v0), ts.vertPtr(e_v1), ts.triVertPtr(o_t_id, 0),
	      ts.triVertPtr(o_t_id, 1), ts.triPlane(o_t_id)) >=
	      SimplexIntersectionType::INTERSECT) // edge e cross seg 0
	{
		seg0_cross = addEdgeCrossCoplanarEdgeInters(o_t_e0, curr_e_id, o_t_id, tc);
		li.insert(seg0_cross);

		if (v0_in_vtx || is_valid_idx(v0_in_seg) || v0_in_tri)
		{
			addSymbolicSegment(e_v0, seg0_cross, e_t_id, o_t_id, tc);
			return true;
		}
		else if (v1_in_vtx || is_valid_idx(v1_in_seg) || v1_in_tri)
		{
			addSymbolicSegment(e_v1, seg0_cross, e_t_id, o_t_id, tc);
			return true;
		}
		else if (tv2_in_edge)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 2), seg0_cross, o_t_id, e_t_id,
			                   tc);
			tc.g.addVertexInEdge(curr_e_id, ts.triVertID(o_t_id, 2));
			li.insert(ts.triVertID(o_t_id, 2));
			return true;
		}
	}

	if (v0_in_seg != o_t_e1 && v1_in_seg != o_t_e1 && !tv1_in_edge &&
	    !tv2_in_edge &&
	    Segment3_Segment3_DoIntersect().intersection_type(
	      ts.vertPtr(e_v0), ts.vertPtr(e_v1), ts.triVertPtr(o_t_id, 1),
	      ts.triVertPtr(o_t_id, 2), ts.triPlane(o_t_id)) >=
	      SimplexIntersectionType::INTERSECT) // edge e cross seg 1
	{
		seg1_cross = addEdgeCrossCoplanarEdgeInters(o_t_e1, curr_e_id, o_t_id, tc);
		li.insert(seg1_cross);

		if (v0_in_vtx || is_valid_idx(v0_in_seg) || v0_in_tri)
		{
			addSymbolicSegment(e_v0, seg1_cross, e_t_id, o_t_id, tc);
			return true;
		}
		else if (v1_in_vtx || is_valid_idx(v1_in_seg) || v1_in_tri)
		{
			addSymbolicSegment(e_v1, seg1_cross, e_t_id, o_t_id, tc);
			return true;
		}
		else if (tv0_in_edge)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 0), seg1_cross, o_t_id, e_t_id,
			                   tc);
			tc.g.addVertexInEdge(curr_e_id, ts.triVertID(o_t_id, 0));
			li.insert(ts.triVertID(o_t_id, 0));
			return true;
		}
	}

	if (v0_in_seg != o_t_e2 && v1_in_seg != o_t_e2 && !tv2_in_edge &&
	    !tv0_in_edge &&
	    Segment3_Segment3_DoIntersect().intersection_type(
	      ts.vertPtr(e_v0), ts.vertPtr(e_v1), ts.triVertPtr(o_t_id, 2),
	      ts.triVertPtr(o_t_id, 0), ts.triPlane(o_t_id)) >=
	      SimplexIntersectionType::INTERSECT) // edge e cross seg2
	{
		seg2_cross = addEdgeCrossCoplanarEdgeInters(o_t_e2, curr_e_id, o_t_id, tc);
		li.insert(seg2_cross);

		if (v0_in_vtx || is_valid_idx(v0_in_seg) || v0_in_tri)
		{
			addSymbolicSegment(e_v0, seg2_cross, e_t_id, o_t_id, tc);
			return true;
		}
		else if (v1_in_vtx || is_valid_idx(v1_in_seg) || v1_in_tri)
		{
			addSymbolicSegment(e_v1, seg2_cross, e_t_id, o_t_id, tc);
			return true;
		}
		else if (tv1_in_edge)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 1), seg2_cross, o_t_id, e_t_id,
			                   tc);
			tc.g.addVertexInEdge(curr_e_id, ts.triVertID(o_t_id, 1));
			li.insert(ts.triVertID(o_t_id, 1));
			return true;
		}
	}

	bool add_sym_seg = false;
	// final probably symbolic edges
	if (is_valid_idx(seg0_cross) && is_valid_idx(seg1_cross))
	{
		addSymbolicSegment(seg0_cross, seg1_cross, e_t_id, o_t_id, tc);
		add_sym_seg = true;
	}
	else if (is_valid_idx(seg0_cross) && is_valid_idx(seg2_cross))
	{
		addSymbolicSegment(seg0_cross, seg2_cross, e_t_id, o_t_id, tc);
		add_sym_seg = true;
	}
	else if (is_valid_idx(seg1_cross) && is_valid_idx(seg2_cross))
	{
		addSymbolicSegment(seg1_cross, seg2_cross, e_t_id, o_t_id, tc);
		add_sym_seg = true;
	}

	if (tv0_in_edge)
	{
		if (is_valid_idx(v0_in_seg) || v0_in_tri)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 0), e_v0, o_t_id, e_t_id, tc);
			add_sym_seg = true;
		}
		else if (is_valid_idx(v1_in_seg) || v1_in_tri)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 0), e_v1, o_t_id, e_t_id, tc);
			add_sym_seg = true;
		}
	}

	if (tv1_in_edge)
	{
		if (is_valid_idx(v0_in_seg) || v0_in_tri)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 1), e_v0, o_t_id, e_t_id, tc);
			add_sym_seg = true;
		}
		else if (is_valid_idx(v1_in_seg) || v1_in_tri)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 1), e_v1, o_t_id, e_t_id, tc);
			add_sym_seg = true;
		}
	}

	if (tv2_in_edge)
	{
		if (is_valid_idx(v0_in_seg) || v0_in_tri)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 2), e_v0, o_t_id, e_t_id, tc);
			add_sym_seg = true;
		}
		else if (is_valid_idx(v1_in_seg) || v1_in_tri)
		{
			addSymbolicSegment(ts.triVertID(o_t_id, 2), e_v1, o_t_id, e_t_id, tc);
			add_sym_seg = true;
		}
	}

	return add_sym_seg;
}

/**
 * @brief One edge \p e_id is not coplanar to triangle \p t_id, check
 * intersection between them.
 * @param lv vertex list
 * @param li intersection list
 */
template <typename Traits>
void ClassifyIntersections<Traits>::checkSingleNoCoplanarEdgeIntersection(
  index_t e_id, index_t t_id, phmap::flat_hash_set<size_t> &lv,
  phmap::flat_hash_set<size_t> &li, ThreadContext &tc)
{
	const NT *s0 = ts.edgeVertPtr(e_id, 0);
	const NT *s1 = ts.edgeVertPtr(e_id, 1);
	const NT *t0 = ts.triVertPtr(t_id, 0);
	const NT *t1 = ts.triVertPtr(t_id, 1);
	const NT *t2 = ts.triVertPtr(t_id, 2);

	Sign vol_s_t01, vol_s_t12, vol_s_t20;
	vol_s_t01 = Orient3D()(s0, s1, t0, t1);
	vol_s_t12 = Orient3D()(s0, s1, t1, t2);
	vol_s_t20 = Orient3D()(s0, s1, t2, t0);

	if ((vol_s_t01 > Sign::ZERO && vol_s_t12 < Sign::ZERO) ||
	    (vol_s_t01 < Sign::ZERO && vol_s_t12 > Sign::ZERO))
		return; // SimplexIntersectionType::DO_NOT_INTERSECT
	if ((vol_s_t12 > Sign::ZERO && vol_s_t20 < Sign::ZERO) ||
	    (vol_s_t12 < Sign::ZERO && vol_s_t20 > Sign::ZERO))
		return; // SimplexIntersectionType::DO_NOT_INTERSECT
	if ((vol_s_t20 > Sign::ZERO && vol_s_t01 < Sign::ZERO) ||
	    (vol_s_t20 < Sign::ZERO && vol_s_t01 > Sign::ZERO))
		return; // SimplexIntersectionType::DO_NOT_INTERSECT;

	// From here, edge must cross triangle at a point.
	// next, we are going to find the intersection point.

	if ((vol_s_t01 == Sign::ZERO && vol_s_t20 == Sign::ZERO) || // edge cross t0
	    (vol_s_t01 == Sign::ZERO && vol_s_t12 == Sign::ZERO) || // edge cross t1
	    (vol_s_t12 == Sign::ZERO && vol_s_t20 == Sign::ZERO))   // edge cross t2
		return;                                                   // no intersection

	// the edge intersects the tri in seg 0
	if (vol_s_t01 == Sign::ZERO)
	{
		index_t e_id_in_tri = ts.triEdgeID(t_id, 0);
		index_t int_point =
		  addEdgeCrossNoCoplanarEdgeInters(e_id, e_id_in_tri, t_id, tc);
		li.insert(int_point);
		lv.insert(int_point);
		return;
	}

	// the edge intersects the tri in seg 1
	if (vol_s_t12 == Sign::ZERO)
	{
		index_t e_id_in_tri = ts.triEdgeID(t_id, 1);
		index_t int_point =
		  addEdgeCrossNoCoplanarEdgeInters(e_id, e_id_in_tri, t_id, tc);
		li.insert(int_point);
		lv.insert(int_point);
		return;
	}

	// the edge intersects the tri in seg 2
	if (vol_s_t20 == Sign::ZERO)
	{
		index_t e_id_in_tri = ts.triEdgeID(t_id, 2);
		index_t int_point =
		  addEdgeCrossNoCoplanarEdgeInters(e_id, e_id_in_tri, t_id, tc);
		li.insert(int_point);
		lv.insert(int_point);
		return;
	}

	// the edge intersects the inner triangle
	index_t int_point = addEdgeCrossTriInters(e_id, t_id, tc);
	li.insert(int_point);
	lv.insert(int_point);
}

/**
 * @brief One vertex \p v_id is coplanar to triangle \p t_id, check intersection
 * between them.
 * @param lv vertex list
 * @param li intersection list
 */
template <typename Traits>
void ClassifyIntersections<Traits>::checkVtxInTriangleIntersection(
  index_t v_id, index_t t_id, phmap::flat_hash_set<size_t> &lv,
  phmap::flat_hash_set<size_t> &li, ThreadContext &tc)
{
	PointInSimplexType inters = Triangle3_Point3_DoIntersect().intersection_type(
	  ts.triVertPtr(t_id, 0), ts.triVertPtr(t_id, 1), ts.triVertPtr(t_id, 2),
	  ts.vertPtr(v_id), ts.triPlane(t_id));

	switch (inters)
	{
	case PointInSimplexType::STRICTLY_OUTSIDE:
		break;

	case PointInSimplexType::ON_EDGE0:
	{
		index_t e_id = ts.triEdgeID(t_id, 0);
		tc.g.addVertexInEdge(e_id, v_id);
		li.insert(v_id);
		lv.insert(v_id);
	}
	break;

	case PointInSimplexType::ON_EDGE1:
	{
		index_t e_id = ts.triEdgeID(t_id, 1);
		tc.g.addVertexInEdge(e_id, v_id);
		li.insert(v_id);
		lv.insert(v_id);
	}
	break;

	case PointInSimplexType::ON_EDGE2:
	{
		index_t e_id = ts.triEdgeID(t_id, 2);
		tc.g.addVertexInEdge(e_id, v_id);
		li.insert(v_id);
		lv.insert(v_id);
	}
	break;

	case PointInSimplexType::STRICTLY_INSIDE:
	{
		tc.g.addVertexInTriangle(t_id, v_id);
		li.insert(v_id);
		lv.insert(v_id);
	}
	break;

	case PointInSimplexType::ON_VERT0:
	case PointInSimplexType::ON_VERT1:
	case PointInSimplexType::ON_VERT2:
	{
		lv.insert(v_id);
		li.insert(v_id);
	}
	break;

	default:
		break;
	}
}

/**
 * @brief Add constraint segment into triangles.
 * @param v01_id indices of two end points of the segment.
 * @param orig_v01_id edge where the constraint segment comes from.
 * @param tAB_id indices of triangles.
 */
template <typename Traits>
void ClassifyIntersections<Traits>::addSymbolicSegment(
  index_t v0_id, index_t v1_id, index_t tA_id, index_t tB_id, ThreadContext &tc)
{
	OMC_EXPENSIVE_ASSERT(v0_id != v1_id,
	                     "trying to add a 0-lenght symbolic edge");
	UIPair seg = uniquePair(v0_id, v1_id);

	if (!ts.triContainsEdge(tA_id, v0_id, v1_id))
		tc.g.addSegmentInTriangle(tA_id, seg);

	if (!ts.triContainsEdge(tB_id, v0_id, v1_id))
		tc.g.addSegmentInTriangle(tB_id, seg);

	tc.g.addTrianglesInSegment(seg, tA_id, tB_id);
}

template <typename Traits>
index_t ClassifyIntersections<Traits>::addEdgeCrossCoplanarEdgeInters(
  index_t e0_id, index_t e1_id, index_t t_id, ThreadContext &tc)
{
	IPoint_SSI *new_v = tc.pnt_arena.emplace(CreateSSI()(
	  ts.edgeVert(e0_id, 0), ts.edgeVert(e0_id, 1), ts.edgeVert(e1_id, 0),
	  ts.edgeVert(e1_id, 1), ts.triPlane(t_id)));

	std::atomic<index_t> *idx_ptr = tc.idx_arena.emplace(InvalidIndex);
	// index creation is deferred until succeesfully inserting point.

	// <index of point, succeed to insert?>
	std::pair<index_t, bool> ins =
	  uniq_g.addVertexInSortedList(new_v, idx_ptr, CreateIndex(ts, uniq_g));
	if (!ins.second)
	{
		OMC_ASSERT(ins.first != InvalidIndex, "");
		IPoint_SSI::gcv().remove(new_v);
		tc.pnt_arena.recycle(new_v);
		tc.idx_arena.recycle(idx_ptr);
	}

	tc.g.addVertexInEdge(e0_id, ins.first);
	tc.g.addVertexInEdge(e1_id, ins.first);

	return ins.first;
}

/**
 * @brief One edge \p e0_id cross another edge \p e1_id in triangle \p t_id,
 * create a new LPI for this case.
 * @return index_t index of new LPI point.
 */
template <typename Traits>
index_t ClassifyIntersections<Traits>::addEdgeCrossNoCoplanarEdgeInters(
  index_t e0_id, index_t e1_id, OMC_UNUSED index_t t_id, ThreadContext &tc)
{
	// find a plane for two interseted segments
	int plane = -1;
	{
		static std::array<index_t, 12> _tri = {0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};

		std::array<const NT *, 4> _pnts = {
		  ts.edgeVert(e0_id, 0).data(), ts.edgeVert(e0_id, 1).data(),
		  ts.edgeVert(e1_id, 0).data(), ts.edgeVert(e1_id, 1).data()};

		OMC_EXPENSIVE_ASSERT(Segment3_Segment3_DoIntersect().intersection_type(
		                       _pnts[0], _pnts[1], _pnts[2], _pnts[3]) >=
		                       SimplexIntersectionType::INTERSECT,
		                     "two edges do not intersect.");

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

	IPoint_SSI *new_v = tc.pnt_arena.emplace(
	  CreateSSI()(ts.edgeVert(e0_id, 0), ts.edgeVert(e0_id, 1),
	              ts.edgeVert(e1_id, 0), ts.edgeVert(e1_id, 1), plane));

	std::atomic<index_t> *idx_ptr = tc.idx_arena.emplace(InvalidIndex);
	// index creation is deferred until succeesfully inserting point.

	// <index of point, succeed to insert?>
	std::pair<index_t, bool> ins =
	  uniq_g.addVertexInSortedList(new_v, idx_ptr, CreateIndex(ts, uniq_g));

	if (!ins.second)
	{
		OMC_ASSERT(ins.first != InvalidIndex, "");
		IPoint_SSI::gcv().remove(new_v);
		tc.pnt_arena.recycle(new_v);
		tc.idx_arena.recycle(idx_ptr);
	}

	tc.g.addVertexInEdge(e0_id, ins.first);
	tc.g.addVertexInEdge(e1_id, ins.first);

	return ins.first;
}

/**
 * @brief One edge \p e_id cross one triangle \p t_id, create a new LPI for this
 * case.
 * @return index_t the index of the new LPI point.
 */
template <typename Traits>
index_t ClassifyIntersections<Traits>::addEdgeCrossTriInters(index_t e_id,
                                                             index_t t_id,
                                                             ThreadContext &tc)
{
	IPoint_LPI *new_v = tc.pnt_arena.emplace(
	  CreateLPI()(ts.edgeVert(e_id, 0), ts.edgeVert(e_id, 1), ts.triVert(t_id, 0),
	              ts.triVert(t_id, 1), ts.triVert(t_id, 2)));
	std::atomic<index_t> *idx_ptr = tc.idx_arena.emplace(InvalidIndex);
	// index creation is deferred until succeesfully inserting point.

	// <index of LPI point, succeed to insert?>
	std::pair<index_t, bool> ins =
	  uniq_g.addVertexInSortedList(new_v, idx_ptr, CreateIndex(ts, uniq_g));

	if (!ins.second)
	{
		OMC_ASSERT(ins.first != InvalidIndex, "");
		IPoint_LPI::gcv().remove(new_v);
		tc.pnt_arena.recycle(new_v);
		tc.idx_arena.recycle(idx_ptr);
	}

	tc.g.addVertexInTriangle(t_id, ins.first);
	tc.g.addVertexInEdge(e_id, ins.first);

	return ins.first;
}

template <typename Traits>
void ClassifyIntersections<Traits>::mergeConcurrentAuxStructures()
{
	auto mergeOneMember = [this](size_t idx)
	{
		switch (idx)
		{
		case 0: // merge coplanar_tris and coplanar_edges
		{
			for (const AuxStruct &cg : concurrent_g)
				for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
				{
					if (!cg.coplanar_tris[t_id].empty())
						uniq_g.coplanar_tris[t_id].insert(uniq_g.coplanar_tris[t_id].end(),
						                                  cg.coplanar_tris[t_id].begin(),
						                                  cg.coplanar_tris[t_id].end());
					if (!cg.coplanar_edges[t_id].empty())
						uniq_g.coplanar_edges[t_id].insert(
						  uniq_g.coplanar_edges[t_id].end(),
						  cg.coplanar_edges[t_id].begin(), cg.coplanar_edges[t_id].end());
				}
		}
		break;
		case 1: // merge tri2pts
		{
			for (const AuxStruct &cg : concurrent_g)
				for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
					for (index_t v_id : cg.trianglePointsList(t_id))
						uniq_g.addVertexInTriangle(t_id, v_id);
		}
		break;
		case 2: // merge edge2pts
		{
			for (const AuxStruct &cg : concurrent_g)
				for (index_t e_id = 0; e_id < ts.numEdges(); e_id++)
					for (index_t v_id : cg.edgePointsList(e_id))
						uniq_g.addVertexInEdge(e_id, v_id);
		}
		break;
		case 3: // merge tri2segs and seg2tris
		{
			for (const AuxStruct &cg : concurrent_g)
				for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
					for (const UIPair &seg : cg.triangleSegmentsList(t_id))
						uniq_g.addSegmentInTriangle(t_id, seg);
		}
		break;
		case 4: // merge seg2tris
		{
			for (const AuxStruct &cg : concurrent_g)
				for (const auto &[seg, tris] : cg.seg2tris)
					for (index_t t_id : tris)
						uniq_g.addTrianglesInSegment(seg, t_id, t_id);
		}
		break;
		}
	};

	tbb::parallel_for(size_t(0), size_t(5), mergeOneMember);

	concurrent_g.clear();
}

template <typename Traits>
void ClassifyIntersections<Traits>::propagateCoplanarTrianglesIntersections()
{
	std::vector<std::pair<index_t, index_t>> coplanar_tris;
	coplanar_tris.reserve(1024);
	for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
		for (index_t copl_t_id : uniq_g.coplanarTriangles(t_id))
			coplanar_tris.emplace_back(t_id, copl_t_id);

	std::vector<std::pair<index_t, index_t>> coplanar_tri_edge;
	coplanar_tri_edge.reserve(1024);
	for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
		for (index_t copl_e_id : uniq_g.coplanarEdges(t_id))
			coplanar_tri_edge.emplace_back(t_id, copl_e_id);

	// tri 2 pnts <copl_t_id, p_id>
	std::vector<std::vector<std::pair<index_t, index_t>>> new_tri2pnts;
	// tri 2 segs
	new_tri2pnts.resize(tbb::this_task_arena::max_concurrency());

	auto propagate_copl_tri =
	  [this, &new_tri2pnts](std::pair<index_t, index_t> copl_pair)
	{
		index_t t_id = copl_pair.first, copl_t_id = copl_pair.second;
		int     thread_id = tbb::this_task_arena::current_thread_index();

		GPoint::clear_global_cached_values();
		index_t e0_id = ts.triEdgeID(t_id, 0);
		index_t e1_id = ts.triEdgeID(t_id, 1);
		index_t e2_id = ts.triEdgeID(t_id, 2);

		bool inside = false;

		// intersection points inside triangle
		for (int i = 1; i < (int)uniq_g.edgePointsList(e0_id).size() - 1; i++)
		{ // skip two endpoints of edge
			index_t p_id = uniq_g.edgePointsList(e0_id)[i];
			if (!ts.triContainsVert(copl_t_id, p_id) &&
			    pointInsideTriangle(p_id, copl_t_id))
			{
				inside = true;
				new_tri2pnts[thread_id].emplace_back(copl_t_id, p_id);
			}
			else if (inside) // we step from inside to outside
				break;
		}

		inside = false;
		for (int i = 1; i < (int)uniq_g.edgePointsList(e1_id).size() - 1; i++)
		{
			index_t p_id = uniq_g.edgePointsList(e1_id)[i];
			if (!ts.triContainsVert(copl_t_id, p_id) &&
			    pointInsideTriangle(p_id, copl_t_id))
			{
				inside = true;
				new_tri2pnts[thread_id].emplace_back(copl_t_id, p_id);
			}
			else if (inside) // we step from inside to outside
				break;
		}

		inside = false;
		for (int i = 1; i < (int)uniq_g.edgePointsList(e2_id).size() - 1; i++)
		{
			index_t p_id = uniq_g.edgePointsList(e2_id)[i];
			if (!ts.triContainsVert(copl_t_id, p_id) &&
			    pointInsideTriangle(p_id, copl_t_id))
			{
				inside = true;
				new_tri2pnts[thread_id].emplace_back(copl_t_id, p_id);
			}
			else if (inside) // we step from inside to outside
				break;
		}
	};

	auto propagate_copl_edge =
	  [this, &new_tri2pnts](std::pair<index_t, index_t> copl_pair)
	{
		index_t t_id = copl_pair.first, copl_e_id = copl_pair.second;
		int     thread_id = tbb::this_task_arena::current_thread_index();

		GPoint::clear_global_cached_values();
		// intersection points inside triangle
		bool inside = false;
		for (int i = 1; i < (int)uniq_g.edgePointsList(copl_e_id).size() - 1; i++)
		{
			index_t p_id = uniq_g.edgePointsList(copl_e_id)[i];
			if (!ts.triContainsVert(t_id, p_id) && pointInsideTriangle(p_id, t_id))
			{
				inside = true;
				new_tri2pnts[thread_id].emplace_back(t_id, p_id);
			}
			else if (inside) // we step from inside to outside
				break;
		}
	};

	std::vector<tbb::spin_mutex> tri_mutexes =
	  std::vector<tbb::spin_mutex>(ts.numTris());

	auto collect = [this, &tri_mutexes, &new_tri2pnts](int thread_id)
	{
		for (const std::pair<index_t, index_t> &t_p : new_tri2pnts[thread_id])
		{
			std::lock_guard<tbb::spin_mutex> lock(tri_mutexes[t_p.first]);
			uniq_g.addVertexInTriangle(t_p.first, t_p.second);
		}
	};

	GPoint::enable_global_cached_values(tbb::this_task_arena::max_concurrency());

	tbb::parallel_for_each(coplanar_tris.begin(), coplanar_tris.end(),
	                       propagate_copl_tri);
	tbb::parallel_for_each(coplanar_tri_edge.begin(), coplanar_tri_edge.end(),
	                       propagate_copl_edge);
	tbb::parallel_for(0, tbb::this_task_arena::max_concurrency(), collect);

	GPoint::disable_global_cached_values();
}

template <typename Traits>
void ClassifyIntersections<Traits>::sortEdgePointsList()
{
	// sort points on edge
	tbb::parallel_for(
	  size_t(0), ts.numEdges(),
	  [this](index_t e_id)
	  {
		  if (uniq_g.edgePointsList(e_id).empty())
			  return;
		  // add two end points into list
		  index_t ev0 = ts.edge(e_id).first, ev1 = ts.edge(e_id).second;

		  std::array<double, 3> dim_diff = {
		    std::fabs(ts.vert(ev0).x() - ts.vert(ev1).x()),
		    std::fabs(ts.vert(ev0).y() - ts.vert(ev1).y()),
		    std::fabs(ts.vert(ev0).z() - ts.vert(ev1).z())};
		  size_t dim =
		    std::max_element(dim_diff.begin(), dim_diff.end()) - dim_diff.begin();

		  uniq_g.addVertexInEdge(e_id, ev0);
		  uniq_g.addVertexInEdge(e_id, ev1);
		  if (dim == 0)
			  uniq_g.sortEdgeList(
			    e_id,
			    [this](index_t a, index_t b) {
				    return LessThan3D().on_x(ts.vert(a), ts.vert(b)) == Sign::NEGATIVE;
			    });
		  else if (dim == 1)
			  uniq_g.sortEdgeList(
			    e_id,
			    [this](index_t a, index_t b) {
				    return LessThan3D().on_y(ts.vert(a), ts.vert(b)) == Sign::NEGATIVE;
			    });
		  else
			  uniq_g.sortEdgeList(
			    e_id,
			    [this](index_t a, index_t b) {
				    return LessThan3D().on_z(ts.vert(a), ts.vert(b)) == Sign::NEGATIVE;
			    });
	  });
}

template <typename Traits>
bool ClassifyIntersections<Traits>::pointInsideTriangle(index_t p_id,
                                                        index_t t_id)
{
	const GPoint &p   = ts.vert(p_id);
	const GPoint &tv0 = ts.triVert(t_id, 0);
	const GPoint &tv1 = ts.triVert(t_id, 1);
	const GPoint &tv2 = ts.triVert(t_id, 2);

	Sign ori = ts.triOrientation(t_id);

	// noexcept
	switch (ts.triPlane(t_id))
	{
	case XY:
	{
		return ((OrientOn2D().on_xy(tv0, tv1, p) == ori &&
		         OrientOn2D().on_xy(tv1, tv2, p) == ori &&
		         OrientOn2D().on_xy(tv2, tv0, p) == ori));
	}

	case YZ:
	{
		return ((OrientOn2D().on_yz(tv0, tv1, p) == ori &&
		         OrientOn2D().on_yz(tv1, tv2, p) == ori &&
		         OrientOn2D().on_yz(tv2, tv0, p) == ori));
	}

	default:
	{
		return ((OrientOn2D().on_zx(tv0, tv1, p) == ori &&
		         OrientOn2D().on_zx(tv1, tv2, p) == ori &&
		         OrientOn2D().on_zx(tv2, tv0, p) == ori));
	}
	}
}

} // namespace OMC