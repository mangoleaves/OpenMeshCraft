#pragma once

#include "DetectClassifyTTI.h"
#include "DetectClassifyTTIs.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"
#include "OpenMeshCraft/Utils/Exception.h"

// #define OMC_ARR_CI_NO_CACHE

namespace OMC {

template <typename Traits>
DetectClassifyTTIs<Traits>::DetectClassifyTTIs(TriSoup &_ts, AuxStruct &_g,
                                               bool _parallel,
                                               MeshArrangements_Stats *_stats,
                                               bool                    _verbose)
  : ts(_ts)
  , pnt_arenas(*ts.pnt_arenas)
  , idx_arenas(*ts.idx_arenas)
  , uniq_g(_g)
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
	}
}

template <typename Traits>
void DetectClassifyTTIs<Traits>::checkTriangleTriangleIntersections()
{
	auto check_tt = [this](const UIPair &pair)
	{
		GPoint::clear_global_cached_values();

		index_t tA_id = pair.first, tB_id = pair.second;

		if (parallel)
		{ // we store results in coccurent aux struct
			int thread_id = tbb::this_task_arena::current_thread_index();
			DetectClassifyTTI<Traits> dc(ts, uniq_g, concurrent_g[thread_id],
			                             pnt_arenas[thread_id],
			                             idx_arenas[thread_id]);
			dc.check_TTI(tA_id, tB_id);
		}
		else
		{ // we store results in unique aux struct
			DetectClassifyTTI<Traits> dc(ts, uniq_g, uniq_g, pnt_arenas[0],
			                             idx_arenas[0]);
			dc.check_TTI(tA_id, tB_id);
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

template <typename Traits>
void DetectClassifyTTIs<Traits>::mergeConcurrentAuxStructures()
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
void DetectClassifyTTIs<Traits>::propagateCoplanarTrianglesIntersections()
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
void DetectClassifyTTIs<Traits>::sortEdgePointsList()
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
bool DetectClassifyTTIs<Traits>::pointInsideTriangle(index_t p_id, index_t t_id)
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