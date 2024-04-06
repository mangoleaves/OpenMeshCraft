#pragma once

#include "DetectClassifyTTI.h"
#include "DetectClassifyTTIs.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"
#include "OpenMeshCraft/Utils/Exception.h"

// #define OMC_ARR_CI_NO_CACHE

namespace OMC {

template <typename Traits>
DetectClassifyTTIs<Traits>::DetectClassifyTTIs(
  TriSoup &_ts, std::vector<UIPair> &_intersection_list,
  MeshArrangements_Stats *_stats, bool _verbose)
  : ts(_ts)
  , pnt_arenas(*ts.pnt_arenas)
  , idx_arenas(*ts.idx_arenas)
  , intersection_list(_intersection_list)
  , verbose(_verbose)
  , stats(_stats)
{
	checkTriangleTriangleIntersections();
	propagateCoplanarTrianglesIntersections();
}

template <typename Traits>
void DetectClassifyTTIs<Traits>::checkTriangleTriangleIntersections()
{
	auto check_tt = [this](const UIPair &pair)
	{
		GPoint::clear_global_cached_values();

		index_t thread_id = tbb::this_task_arena::current_thread_index();

		index_t tA_id = pair.first, tB_id = pair.second;

		DetectClassifyTTI<Traits> dc(ts, pnt_arenas[thread_id],
		                             idx_arenas[thread_id]);
		dc.check_TTI(tA_id, tB_id);
	};

	GPoint::enable_global_cached_values(tbb::this_task_arena::max_concurrency());

#if 1
	tbb::parallel_for_each(intersection_list.begin(), intersection_list.end(),
	                       check_tt);
#else
	std::for_each(intersection_list.begin(), intersection_list.end(), check_tt);
#endif

	ts.removeAllDuplicates();
	ts.calcPlaneAndOrient();

	sortEdgePointsList();

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
		  if (ts.edgePointsList(e_id).empty())
			  return;
		  // add two end points into list
		  index_t ev0 = ts.edge(e_id).first, ev1 = ts.edge(e_id).second;

		  std::array<double, 3> dim_diff = {
		    std::fabs(ts.vert(ev0).x() - ts.vert(ev1).x()),
		    std::fabs(ts.vert(ev0).y() - ts.vert(ev1).y()),
		    std::fabs(ts.vert(ev0).z() - ts.vert(ev1).z())};
		  size_t dim =
		    std::max_element(dim_diff.begin(), dim_diff.end()) - dim_diff.begin();

		  ts.addVertexInEdge(e_id, ev0);
		  ts.addVertexInEdge(e_id, ev1);
		  if (dim == 0)
			  ts.sortEdgeList(e_id,
			                  [this](index_t a, index_t b) {
				                  return LessThan3D().on_x(ts.vert(a), ts.vert(b)) ==
				                         Sign::NEGATIVE;
			                  });
		  else if (dim == 1)
			  ts.sortEdgeList(e_id,
			                  [this](index_t a, index_t b) {
				                  return LessThan3D().on_y(ts.vert(a), ts.vert(b)) ==
				                         Sign::NEGATIVE;
			                  });
		  else
			  ts.sortEdgeList(e_id,
			                  [this](index_t a, index_t b) {
				                  return LessThan3D().on_z(ts.vert(a), ts.vert(b)) ==
				                         Sign::NEGATIVE;
			                  });
	  });
}

template <typename Traits>
void DetectClassifyTTIs<Traits>::propagateCoplanarTrianglesIntersections()
{
	std::vector<std::pair<index_t, index_t>> coplanar_tris;
	coplanar_tris.reserve(1024);
	for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
		for (index_t copl_t_id : ts.coplanarTriangles(t_id))
			coplanar_tris.emplace_back(t_id, copl_t_id);

	std::vector<std::pair<index_t, index_t>> coplanar_tri_edge;
	coplanar_tri_edge.reserve(1024);
	for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
		for (index_t copl_e_id : ts.coplanarEdges(t_id))
			coplanar_tri_edge.emplace_back(t_id, copl_e_id);

	// tri 2 pnts <t_id, p_id>
	std::vector<std::vector<std::pair<index_t, index_t>>> new_tri2pnts;
	new_tri2pnts.resize(tbb::this_task_arena::max_concurrency());
	// mutexes on tri
	std::vector<tbb::spin_mutex> tri_mutexes(ts.numTris());

	auto propagate_copl_tri =
	  [this, &new_tri2pnts](std::pair<index_t, index_t> copl_pair)
	{
		index_t t_id = copl_pair.first, copl_t_id = copl_pair.second;
		int     thread_id = tbb::this_task_arena::current_thread_index();

		index_t e0_id = ts.triEdgeID(t_id, 0);
		index_t e1_id = ts.triEdgeID(t_id, 1);
		index_t e2_id = ts.triEdgeID(t_id, 2);

		index_t cv0_id = ts.triVertID(copl_t_id, 0);
		index_t cv1_id = ts.triVertID(copl_t_id, 1);
		index_t cv2_id = ts.triVertID(copl_t_id, 2);

		// TODO // OPT
		// How often does an edge really intersect triangle in coplanar triangle
		// pair? Do we need to store exact coplanar tri-edge pair?
		for (index_t e_id : {e0_id, e1_id, e2_id})
			if (is_valid_idx(e_id))
			{
				bool inside = false;

				const tbb::concurrent_vector<index_t> &e2p = ts.edgePointsList(e_id);

				for (int i = 1; i < (int)e2p.size() - 1; i++)
				{ // skip two endpoints of edge
					index_t p_id = e2p[i];
					/*p_id is not on vertex and p_id insides triangle */
					if (p_id != cv0_id && p_id != cv1_id && p_id != cv2_id &&
					    pointInsideTriangle(p_id, copl_t_id))
					{
						inside = true;
						new_tri2pnts[thread_id].emplace_back(copl_t_id, p_id);
					}
					else if (inside) // we step from inside to outside
						break;
				}
			}
	};

	auto propagate_copl_edge =
	  [this, &new_tri2pnts](std::pair<index_t, index_t> copl_pair)
	{
		index_t t_id = copl_pair.first, copl_e_id = copl_pair.second;
		int     thread_id = tbb::this_task_arena::current_thread_index();

		index_t v0_id = ts.triVertID(t_id, 0);
		index_t v1_id = ts.triVertID(t_id, 1);
		index_t v2_id = ts.triVertID(t_id, 2);

		// intersection points inside triangle
		bool inside = false;

		const tbb::concurrent_vector<index_t> &e2p = ts.edgePointsList(copl_e_id);

		for (int i = 1; i < (int)e2p.size() - 1; i++)
		{
			index_t p_id = e2p[i];
			if (p_id != v0_id && p_id != v1_id && p_id != v2_id &&
			    pointInsideTriangle(p_id, t_id))
			{
				inside = true;
				new_tri2pnts[thread_id].emplace_back(t_id, p_id);
			}
			else if (inside) // we step from inside to outside
				break;
		}
	};

	auto collect = [this, &tri_mutexes, &new_tri2pnts](int thread_id)
	{
		for (const std::pair<index_t, index_t> &t_p : new_tri2pnts[thread_id])
		{
			std::lock_guard<tbb::spin_mutex> lock(tri_mutexes[t_p.first]);

			const tbb::concurrent_vector<index_t> &t2p =
			  ts.trianglePointsList(t_p.first);
			// check if point is in triangle
			if (std::find(t2p.begin(), t2p.end(), t_p.second) == t2p.end())
				// if not, add it
				ts.addVertexInTriangle(t_p.first, t_p.second);
		}
	};

#if 1
	tbb::parallel_for_each(coplanar_tris.begin(), coplanar_tris.end(),
	                       propagate_copl_tri);
	tbb::parallel_for_each(coplanar_tri_edge.begin(), coplanar_tri_edge.end(),
	                       propagate_copl_edge);
	tbb::parallel_for(0, tbb::this_task_arena::max_concurrency(), collect);
#else
	std::for_each(coplanar_tris.begin(), coplanar_tris.end(), propagate_copl_tri);
	std::for_each(coplanar_tri_edge.begin(), coplanar_tri_edge.end(),
	              propagate_copl_edge);
	std::ranges::for_each(
	  std::views::iota(0, tbb::this_task_arena::max_concurrency()), collect);
#endif
}

template <typename Traits>
bool DetectClassifyTTIs<Traits>::pointInsideTriangle(index_t p_id, index_t t_id)
{
	// TODO // OPT
	// merge this func to where it is called. avoid duplicate calls to ts.vert,
	// ts.triVert, ts.triOrientation and ts.triPlane.
	const GPoint &p   = ts.vert(p_id);
	const GPoint &tv0 = ts.triVert(t_id, 0);
	const GPoint &tv1 = ts.triVert(t_id, 1);
	const GPoint &tv2 = ts.triVert(t_id, 2);

	Sign ori = ts.triOrient(t_id);

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