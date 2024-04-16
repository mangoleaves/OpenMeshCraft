#pragma once

#include "Triangulation.h"

#include <execution>
#include <stack>

namespace OMC {

template <typename Traits>
struct Triangulation<Traits>::CreateIndex
{
	TriSoup              &ts;
	const GPoint         *pp;
	std::atomic<index_t> *ip;

	CreateIndex() = default;
	CreateIndex(TriSoup &_ts, const GPoint *_pp, std::atomic<index_t> *_ip)
	  : ts(_ts)
	  , pp(_pp)
	  , ip(_ip)
	{
	}

	index_t operator()()
	{
		index_t idx = InvalidIndex;
		{ // lock for new index
			std::lock_guard<tbb::spin_mutex> lock(ts.new_vertex_mutex);
			idx = ts.addImplVert(const_cast<GPoint *>(pp), ip);
		}
		ip->store(idx, std::memory_order_relaxed); // assign a valid index
		return idx;
	};
};

template <typename Traits>
Triangulation<Traits>::Triangulation(TriSoup              &_ts,
                                     std::vector<index_t> &new_tris,
                                     std::vector<Label>   &new_labels)
  : ts(_ts)
  , pnt_arenas(*ts.pnt_arenas)
  , idx_arenas(*ts.idx_arenas)
{
	new_labels.clear();
	new_tris.clear();
	new_tris.reserve(2 * 3 * ts.numTris());
	new_labels.reserve(2 * ts.numTris());
	ts.vertices.reserve(ts.numVerts() * 3);
	ts.indices.reserve(ts.numVerts() * 3);

	std::vector<index_t> tris_to_split;
	tris_to_split.reserve(ts.numTris());

	for (index_t t_id = 0; t_id < ts.numTris(); t_id++)
	{
		if (ts.triangleHasIntersections(t_id))
			tris_to_split.push_back(t_id);
		else
		{
			// triangle without intersections directly goes to the output list
			new_tris.push_back(ts.triVertID(t_id, 0));
			new_tris.push_back(ts.triVertID(t_id, 1));
			new_tris.push_back(ts.triVertID(t_id, 2));
			new_labels.push_back(ts.triLabel(t_id));
		}
	}

	// processing the triangles to split
	GPoint::enable_global_cached_values(tbb::this_task_arena::max_concurrency());
#if 1
	std::shuffle(tris_to_split.begin(), tris_to_split.end(),
	             std::mt19937(std::random_device()()));
	std::vector<FastTriMesh> subms(tbb::this_task_arena::max_concurrency());
	tbb::parallel_for_each(
	  tris_to_split.begin(), tris_to_split.end(),
	  [&](size_t t_id)
	  {
		  int thread_id = tbb::this_task_arena::current_thread_index();
		  subms[thread_id].initialize(&ts.triVert(t_id, 0), &ts.triVert(t_id, 1),
		                              &ts.triVert(t_id, 2), ts.tri(t_id),
		                              ts.triPlane(t_id), ts.triOrient(t_id));
		  triangulateSingleTriangle(t_id, subms[thread_id], new_tris, new_labels);
	  });
#else
	FastTriMesh subm;
	std::for_each(tris_to_split.begin(), tris_to_split.end(),
	              [&](size_t t_id)
	              {
		              subm.initialize(&ts.triVert(t_id, 0), &ts.triVert(t_id, 1),
		                              &ts.triVert(t_id, 2), ts.tri(t_id),
		                              ts.triPlane(t_id), ts.triOrient(t_id));
		              triangulateSingleTriangle(t_id, subm, new_tris, new_labels);
	              });
#endif
	GPoint::disable_global_cached_values();

	postFixIndices(new_tris, new_labels);
}

template <typename Traits>
void Triangulation<Traits>::triangulateSingleTriangle(
  index_t t_id, FastTriMesh &subm, std::vector<index_t> &new_tris,
  std::vector<Label> &new_labels)
{
	GPoint::clear_global_cached_values();
	subm.setMeshInfo(t_id);
	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 *                                  POINTS AND SEGMENTS RECOVERY
	 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	std::vector<index_t> t_points(ts.trianglePointsList(t_id).begin(),
	                              ts.trianglePointsList(t_id).end());

	std::vector<index_t> e0_points, e1_points, e2_points;

	if (index_t e0_id = ts.triEdgeID(t_id, 0); is_valid_idx(e0_id))
		sortedVertexListAlongSegment(ts.edgePointsList(e0_id), subm.vertInfo(0),
		                             subm.vertInfo(1), e0_points);
	if (index_t e1_id = ts.triEdgeID(t_id, 1); is_valid_idx(e1_id))
		sortedVertexListAlongSegment(ts.edgePointsList(e1_id), subm.vertInfo(1),
		                             subm.vertInfo(2), e1_points);
	if (index_t e2_id = ts.triEdgeID(t_id, 2); is_valid_idx(e2_id))
		sortedVertexListAlongSegment(ts.edgePointsList(e2_id), subm.vertInfo(2),
		                             subm.vertInfo(0), e2_points);

	SegmentsList t_segments;
	t_segments.reserve(ts.triangleSegmentsList(t_id).size());
	for (index_t seg_id : ts.triangleSegmentsList(t_id))
		t_segments.push_back(Segment(seg_id, ts.segment(seg_id)));

	size_t estimated_vert_num =
	  t_points.size() + e0_points.size() + e1_points.size() + e2_points.size();
	subm.preAllocateSpace(estimated_vert_num);

	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 *                                  TRIANGLE SPLIT
	 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	if (t_points.size() <= 10)
		splitSingleTriangle(subm, t_points);
	else
		splitSingleTriangleWithTree(subm, t_points);

	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 *                                  EDGE SPLIT
	 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	splitSingleEdge(subm, 0, 1, e0_points);
	splitSingleEdge(subm, 1, 2, e1_points);
	splitSingleEdge(subm, 2, 0, e2_points);

	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 *                           CONSTRAINT SEGMENT INSERTION
	 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	addConstraintSegmentsInSingleTriangle(subm, t_segments);

	/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	 *                      POCKETS IN COPLANAR TRIANGLES SOLVING
	 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

	if (!ts.coplanarTriangles(t_id).empty())
	{
		solvePocketsInCoplanarTriangle(subm, new_tris, new_labels,
		                               ts.triLabel(t_id));
	}
	else
	{
		/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		 *                     NEW TRIANGLE CREATION (for final mesh)
		 * :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

		{ // start critical section...
			std::lock_guard<tbb::spin_mutex> lock(ts.new_tris_mutex);
			for (index_t ti = 0; ti < subm.numTriangles(); ++ti)
			{
				const index_t *tri = subm.tri(ti);
				new_tris.push_back(subm.vertInfo(tri[0]));
				new_tris.push_back(subm.vertInfo(tri[1]));
				new_tris.push_back(subm.vertInfo(tri[2]));
				new_labels.push_back(ts.triLabel(t_id));
			}
		} // end critical section
	}
}

template <typename Traits>
void Triangulation<Traits>::sortedVertexListAlongSegment(
  const typename TriSoup::Edge2PntsSet &point_list, index_t v0_id,
  index_t v1_id, std::vector<index_t> &out_point_list)
{
	if (point_list.empty())
		return;
	out_point_list = std::vector<index_t>(point_list.begin(), point_list.end());
	if (out_point_list.front() == v1_id) // first element must be v0
		std::reverse(out_point_list.begin(), out_point_list.end());
	OMC_EXPENSIVE_ASSERT(out_point_list.front() == v0_id, "missing endpoint.");
}

template <typename Traits>
void Triangulation<Traits>::splitSingleTriangle(
  FastTriMesh &subm, const std::vector<size_t> &points)
{
	if (points.empty())
		return;

	std::vector<uint8_t> tri_visited(points.size() * 4, false);

	// traverse triangles in subm to find which on contains the point
	auto locatePointWalking =
	  [this, &subm, &tri_visited](const GPoint &p, index_t curr_tid)
	{
		std::fill(tri_visited.begin(), tri_visited.end(), false);
		Sign subm_ori = subm.Orientation();
		while (true)
		{
			index_t v0 = subm.triVertID(curr_tid, 0);
			index_t v1 = subm.triVertID(curr_tid, 1);
			index_t v2 = subm.triVertID(curr_tid, 2);
			index_t e0 = subm.triEdgeID(curr_tid, 0);
			index_t e1 = subm.triEdgeID(curr_tid, 1);
			index_t e2 = subm.triEdgeID(curr_tid, 2);
			// clang-format off
			Sign ori0 = OrientOn2D()(subm.vert(v0), subm.vert(v1), p, planeToInt(subm.refPlane()));
			Sign ori1 = OrientOn2D()(subm.vert(v1), subm.vert(v2), p, planeToInt(subm.refPlane()));
			Sign ori2 = OrientOn2D()(subm.vert(v2), subm.vert(v0), p, planeToInt(subm.refPlane()));
			// clang-format on

			if ((ori0 == Sign::ZERO || ori0 == subm_ori) &&
			    (ori1 == Sign::ZERO || ori1 == subm_ori) &&
			    (ori2 == Sign::ZERO || ori2 == subm_ori))
			{ // inside or on boundary
				if (ori0 == subm_ori && ori1 == subm_ori && ori2 == subm_ori)
				{ // inside
					return std::pair<index_t, bool>(curr_tid, true);
				}
				else
				{ // on edge
					if (ori0 == Sign::ZERO)
						return std::pair<index_t, bool>(e0, false);
					if (ori1 == Sign::ZERO)
						return std::pair<index_t, bool>(e1, false);
					if (ori2 == Sign::ZERO)
						return std::pair<index_t, bool>(e2, false);
				}
			}
			else
			{ // outside, walk to another triangle
				tri_visited[curr_tid] = true;
				if (ori0 == reverse_sign(subm_ori) && !subm.edgeIsBoundary(e0) &&
				    !tri_visited[subm.triOppToEdge(e0, curr_tid)])
					curr_tid = subm.triOppToEdge(e0, curr_tid);
				else if (ori1 == reverse_sign(subm_ori) && !subm.edgeIsBoundary(e1) &&
				         !tri_visited[subm.triOppToEdge(e1, curr_tid)])
					curr_tid = subm.triOppToEdge(e1, curr_tid);
				else if (ori2 == reverse_sign(subm_ori) && !subm.edgeIsBoundary(e2) &&
				         !tri_visited[subm.triOppToEdge(e2, curr_tid)])
					curr_tid = subm.triOppToEdge(e2, curr_tid);
				else
					break;
			}
		}

		// this should not happen
		OMC_ASSERT(false, "No containing triangle found!");
		return std::pair<index_t, bool>(InvalidIndex, true);
	};

	// add the first point
	auto    curr  = points.begin();
	index_t v_pos = subm.addVert(&ts.vert(*curr), *curr);
	subm.splitTri(0, v_pos);

	index_t last_tid = 0;

	// progressively add the other points, looking for the triangle
	// that contains them only among the newly generated triangles
	while (++curr != points.end())
	{
		v_pos = subm.addVert(&ts.vert(*curr), *curr);

		// find triangle from last found one.
		auto [cont_id, is_tri] = locatePointWalking(subm.vert(v_pos), last_tid);

		if (is_tri)
		{
			subm.splitTri(cont_id, v_pos);
			last_tid = cont_id;
		}
		else
			subm.splitEdge(cont_id, v_pos);
	}
}

template <typename Traits>
void Triangulation<Traits>::splitSingleTriangleWithTree(
  FastTriMesh &subm, const std::vector<size_t> &points)
{
	if (points.empty())
		return;

	// tree and the root node
	SplitTree tree(points.size() * 4);
	size_t    n_id = tree.addNode(0, 1, 2);
	subm.setTriInfo(0, n_id);

	// add the first point to split triangle
	auto    curr  = points.begin();
	index_t v_pos = subm.addVert(&ts.vert(*curr), *curr);
	subm.splitTri(0, v_pos, tree);

	// progressively add the other points, looking for the triangle
	// that contains them only among the newly generated triangles
	while (++curr != points.end())
	{
		v_pos = subm.addVert(&ts.vert(*curr), *curr);

		auto [cont_id, is_tri] = locatePointInTree(subm, v_pos, tree);

		OMC_EXPENSIVE_ASSERT(is_valid_idx(cont_id),
		                     "No containing triangle found!");

		if (is_tri)
			subm.splitTri(cont_id, v_pos, tree);
		else
			subm.splitEdge(cont_id, v_pos, tree);
	}
}

template <typename Traits>
void Triangulation<Traits>::splitSingleEdge(FastTriMesh &subm, index_t v0_id,
                                            index_t               v1_id,
                                            std::vector<index_t> &points)
{
	if (points.size() <= 2) // only two end points
		return;

	index_t e_id = subm.edgeID(v0_id, v1_id);
	OMC_EXPENSIVE_ASSERT(is_valid_idx(e_id), "invalid edge id");

	// new_vertices in mesh
	for (size_t p_pos = 1; p_pos < points.size() - 1; p_pos++)
	{
		size_t p_id   = points[p_pos];
		size_t v_pos  = subm.addVert(&ts.vert(p_id), p_id);
		points[p_pos] = v_pos;
	}

	points.front() = v0_id;
	points.back()  = v1_id;

	// avoid error caused by reallocating
	AuxVector4<index_t> e2t = subm.adjE2T(e_id);
	// make new_triangles
	for (auto i = points.begin(), j = i + 1; j < points.end(); ++i, ++j)
	{
		for (index_t t_id : e2t)
		{
			index_t opp = subm.triVertOppositeTo(t_id, v0_id, v1_id);

			if (subm.nextVertInTri(t_id, v0_id) != v1_id)
				subm.addTri(*j, *i, opp);
			else
				subm.addTri(*i, *j, opp);
		}
	}

	// remove the original edge and the tris attached to it
	subm.removeEdge(e_id);
}

template <typename Traits>
std::pair<index_t, bool>
Triangulation<Traits>::locatePointInTree(const FastTriMesh &subm, index_t p_id,
                                         const SplitTree &tree)
{
	return locatePointInTreeRecur(subm, subm.vert(p_id), tree, 0,
	                              UIPair(InvalidIndex, InvalidIndex));
}

template <typename Traits>
std::pair<index_t, bool> Triangulation<Traits>::locatePointInTreeRecur(
  const FastTriMesh &subm, const GPoint &p, const SplitTree &tree,
  index_t node_id, UIPair ev)
{
	const SplitTree::Node &nd = tree.getNode(node_id);

	Sign orientation = subm.Orientation();

	if (!is_valid_idx(nd.c[0])) // leaf node
	{
		if (is_valid_idx(ev.first))
			return std::pair<index_t, bool>(subm.edgeID(ev.first, ev.second),
			                                /*is triangle?*/ false);
		else
			return std::pair<index_t, bool>(subm.triID(nd.v[0], nd.v[1], nd.v[2]),
			                                /*is triangle?*/ true);
	}

	// check its children
	if (is_valid_idx(nd.c[2]))
	{ // three child nodes
		index_t ps = nd.ps;
		if (is_valid_idx(ev.first))
		{ // check which edge matches...
			for (int j = 0; j < 3; j++)
				if ((nd.v[j] == ev.first && nd.v[(j + 1) % 3] == ev.second) ||
				    (nd.v[j] == ev.second && nd.v[(j + 1) % 3] == ev.first))
					return locatePointInTreeRecur(subm, p, tree, nd.c[j], ev);
		}
		Sign ori[3] = {OrientOn2D()(subm.vert(ps), p, subm.vert(nd.v[0]),
		                            planeToInt(subm.refPlane())),
		               OrientOn2D()(subm.vert(ps), p, subm.vert(nd.v[1]),
		                            planeToInt(subm.refPlane())),
		               OrientOn2D()(subm.vert(ps), p, subm.vert(nd.v[2]),
		                            planeToInt(subm.refPlane()))};
		for (int j = 0; j < 3; j++)
		{
			if (ori[j] == Sign::ZERO && ori[(j + 1) % 3] == orientation &&
			    ori[(j + 2) % 3] == reverse_sign(orientation))
			{
				return locatePointInTreeRecur(subm, p, tree, nd.c[j],
				                              UIPair(ps, nd.v[j]));
			}
			else if (ori[j] == reverse_sign(orientation) &&
			         ori[(j + 1) % 3] == orientation)
			{
				return locatePointInTreeRecur(subm, p, tree, nd.c[j],
				                              UIPair(InvalidIndex, InvalidIndex));
			}
		}
	}
	else
	{ // two child nodes
		index_t ps = nd.ps;
		index_t vp = nd.v[0], vn = nd.v[1], vo = nd.v[2];
		index_t np = nd.c[0], nn = nd.c[1];
		if ((vo == ev.first && vp == ev.second) ||
		    (vo == ev.second && vp == ev.first))
			return locatePointInTreeRecur(subm, p, tree, np, ev);
		if ((vo == ev.first && vn == ev.second) ||
		    (vo == ev.second && vn == ev.first))
			return locatePointInTreeRecur(subm, p, tree, nn, ev);
		Sign ori = OrientOn2D()(subm.vert(ps), subm.vert(vo), p,
		                        planeToInt(subm.refPlane()));
		if (ori == Sign::ZERO)
			return locatePointInTreeRecur(subm, p, tree, np, UIPair(ps, vo));

		index_t nd_rec = ori == orientation ? np : nn;

		if ((vp == ev.first && vn == ev.second) ||
		    (vp == ev.second && vn == ev.first))
			ev = ori == orientation ? UIPair(ps, vp) : UIPair(ps, vn);
		else
			ev = UIPair(InvalidIndex, InvalidIndex);

		return locatePointInTreeRecur(subm, p, tree, nd_rec, ev);
	}

	OMC_ASSERT(false, "no containing triangle found");
	return std::pair<index_t, bool>(InvalidIndex, true); // warning killer
}

template <typename Traits>
void Triangulation<Traits>::addConstraintSegmentsInSingleTriangle(
  FastTriMesh &subm, SegmentsList &segment_list)
{
	// change segment's global index to local index
	for (Segment &seg : segment_list)
		seg.second = uniquePair(subm.vertLocalID(seg.second.first),
		                        subm.vertLocalID(seg.second.second));
	// add segments to map, map local segment to global segment id
	SubSegMap sub_segs_map;
	sub_segs_map.reserve(segment_list.size());
	for (const Segment &seg : segment_list)
		sub_segs_map[seg.second].insert(seg.first);
	// store segments adjacent to a TPI point in tpi2segs
	TPI2Segs tpi2segs;

	// add segments to triangle mesh
	while (segment_list.size() > 0)
	{
		Segment seg = segment_list.back();
		segment_list.pop_back();

		addConstraintSegment(subm, seg, segment_list, sub_segs_map, tpi2segs);
	}
}

template <typename Traits>
void Triangulation<Traits>::addConstraintSegment(FastTriMesh   &subm,
                                                 const Segment &seg,
                                                 SegmentsList  &segment_list,
                                                 SubSegMap     &sub_segs_map,
                                                 TPI2Segs      &tpi2segs)
{
	// get two local endpoints
	index_t v0_id = seg.second.first;
	index_t v1_id = seg.second.second;

	// check if edge already present in the mesh, just flag it as constraint
	if (index_t e_id = subm.edgeID(v0_id, v1_id); is_valid_idx(e_id))
	{
		subm.setEdgeConstr(e_id);
		return;
	}

	// for efficiency:
	// * it's better to start from the vert with lower valence.
	// * it's also better to start from a TPI vertex.
	index_t v_start = v0_id, v_stop = v1_id;
	if (!subm.vertFlag(v0_id) && !subm.vertFlag(v1_id) ||
	    subm.vertFlag(v0_id) && subm.vertFlag(v1_id))
	{
		if (subm.vertValence(v1_id) < subm.vertValence(v0_id))
			std::swap(v_start, v_stop);
	}
	else if (subm.vertFlag(v1_id))
	{
		std::swap(v_start, v_stop);
	}

	// find intersected edges and triangles of current segment
	std::vector<index_t> intersected_edges, intersected_tris;
	intersected_edges.reserve(64), intersected_tris.reserve(64);

	findIntersectingElements(subm, v_start, v_stop, intersected_edges,
	                         intersected_tris, segment_list, sub_segs_map,
	                         tpi2segs);

	if (intersected_edges.size() == 0)
		return;

	// walk along the border
	std::vector<index_t> h0, h1;
	h0.reserve(64), h1.reserve(64);
	boundaryWalker(subm, v_start, v_stop, intersected_tris.begin(),
	               intersected_edges.begin(), h0);
	boundaryWalker(subm, v_stop, v_start, intersected_tris.rbegin(),
	               intersected_edges.rbegin(), h1);

	OMC_EXPENSIVE_ASSERT(h0.size() >= 3, "insufficient edges of border");
	OMC_EXPENSIVE_ASSERT(h1.size() >= 3, "insufficient edges of border");

	// cut ears
	std::vector<index_t> new_tris;
	new_tris.reserve(64);
	earcutLinear(subm, h0, new_tris);
	earcutLinear(subm, h1, new_tris);

	// add new triangles and remove old triangles
	for (size_t i = 0; i < new_tris.size(); i += 3)
		subm.addTri(new_tris[i], new_tris[i + 1], new_tris[i + 2]);
	subm.removeTris(intersected_tris);

	// set segment as contrained.
	index_t e_id = subm.edgeID(v_start, v_stop);
	OMC_EXPENSIVE_ASSERT(is_valid_idx(e_id), "invalid edge");
	subm.setEdgeConstr(e_id); // edge marked as constr
}

template <typename Traits>
void Triangulation<Traits>::findIntersectingElements(
  FastTriMesh &subm, index_t &v_start, index_t &v_stop,
  std::vector<index_t> &intersected_edges,
  std::vector<index_t> &intersected_tris, SegmentsList &segment_list,
  SubSegMap &sub_segs_map, TPI2Segs &tpi2segs)
{
	using IdxPair  = std::pair<index_t, index_t>;
	using SignPair = std::pair<Sign, Sign>;

	/*******************************************************************/
	/*          Search the first intersecting element                  */
	/*******************************************************************/

	size_t link_size = subm.adjV2T(v_start).size();

	auto firstIntersectVertex = [&](index_t v_inter) -> void
	{
		// current segment (v_start, v_stop) is split in (v_start, v_inter) and
		// (v_inter, v_stop).

		// mark (v_start, v_inter) as constr
		index_t edge_id = subm.edgeID(v_start, v_inter);
		OMC_EXPENSIVE_ASSERT(is_valid_idx(edge_id), "invalid edge");
		subm.setEdgeConstr(edge_id);
		// split the original edge.
		splitSegmentInSubSegments(v_start, v_stop, v_inter, sub_segs_map);
		// push (v_inter, v_stop) to check list.
		// it is a sub-segment so we set its index to invalid index.
		segment_list.emplace_back(InvalidIndex, uniquePair(v_inter, v_stop));
		// check if current segment encounters a TPI point
		if (/*isTPI*/ subm.vertFlag(v_inter))
		{ // if v_inter is a TPI point, fix indices stored in related segments.
			const RefSegs &seg_ids    = sub_segs_map.at(uniquePair(v_start, v_stop));
			std::vector<index_t> &t2s = tpi2segs.at(v_inter);
			for (index_t seg_id : seg_ids)
			{
				if (std::find(t2s.begin(), t2s.end(), seg_id) != t2s.end())
					continue;
				IPoint_TPI *tpi_inter = (IPoint_TPI *)&subm.vert(v_inter);
				for (index_t other_seg_id : t2s)
				{
					index_t suitable_vid = fixTPI(seg_id, other_seg_id, tpi_inter);
					if (suitable_vid < subm.vertInfo(v_inter))
						subm.setVertInfo(v_inter, suitable_vid);
				}
				t2s.push_back(seg_id);
			}
		}
		// clear and stop this traversal.
		intersected_edges.clear();
		intersected_tris.clear();
	};

	/// @brief find the edge in link(seed) that intersect {A,B}
	/// @return true if found.
	auto sequencialSearchInterEdge = [&]() -> bool
	{
		const Sign orientation = subm.Orientation();
		const bool CCW         = orientation == Sign::POSITIVE;

		IdxPair  last_ev;
		SignPair last_ori;

		index_t  curr_tri = subm.adjV2T(v_start)[0];
		IdxPair  curr_ev{subm.nextVertInTri(curr_tri, v_start),
                    subm.prevVertInTri(curr_tri, v_start)};
		SignPair curr_ori{
		  OrientOn2D()(subm.vert(v_start), subm.vert(v_stop),
		               subm.vert(curr_ev.first), planeToInt(subm.refPlane())),
		  OrientOn2D()(subm.vert(v_start), subm.vert(v_stop),
		               subm.vert(curr_ev.second), planeToInt(subm.refPlane()))};

		// sequencially walk around v_start to find intersected vertex/edge.
		for (size_t i = 0; i < subm.adjV2T(v_start).size(); i++)
		{
			if (curr_ori.first == Sign::ZERO &&
			    pointInsideSegmentCollinear(subm, v_start, v_stop, curr_ev.first))
			{ // meet the collinear point, check if it is in the constraint segment
				firstIntersectVertex(curr_ev.first);
				return false;
			}
			if (curr_ori.second == Sign::ZERO &&
			    pointInsideSegmentCollinear(subm, v_start, v_stop, curr_ev.second))
			{ // meet the collinear point, check if it is in the constraint segment
				firstIntersectVertex(curr_ev.second);
				return false;
			}

			if (curr_ori.first != curr_ori.second)
			{ // meet edge that possibly crosses the contraint segment
				if (curr_ori.first == reverse_sign(orientation) &&
				    curr_ori.second == orientation)
				{
					intersected_edges.push_back(
					  subm.edgeID(curr_ev.first, curr_ev.second));
					intersected_tris.push_back(curr_tri);
					return true;
				}
			}

			bool step_dir;
			if (curr_ori.first <= Sign::ZERO && curr_ori.second <= Sign::ZERO)
				step_dir = CCW; // curr_tri is at negative side
			else if (curr_ori.first >= Sign::ZERO && curr_ori.second >= Sign::ZERO)
				step_dir = !CCW; // curr_tri is at positive side
			else
				step_dir = CCW; // curr_tri overlaps both sides, either direction is ok

			last_ev  = curr_ev;
			last_ori = curr_ori;

			// gonna step to next triangle
			curr_tri =
			  subm.rotateAroundVertex(v_start, curr_tri, /*step*/ 1, step_dir);
			curr_ev = IdxPair{subm.nextVertInTri(curr_tri, v_start),
			                  subm.prevVertInTri(curr_tri, v_start)};
			// clang-format off
			if      (curr_ev.first == last_ev.first)   curr_ori.first = last_ori.first;
			else if (curr_ev.first == last_ev.second)  curr_ori.first = last_ori.second;
			else    curr_ori.first = OrientOn2D()(subm.vert(v_start), subm.vert(v_stop), subm.vert(curr_ev.first), planeToInt(subm.refPlane()));

			if      (curr_ev.second == last_ev.first)  curr_ori.second = last_ori.first;
			else if (curr_ev.second == last_ev.second) curr_ori.second = last_ori.second;
			else    curr_ori.second = OrientOn2D()(subm.vert(v_start), subm.vert(v_stop), subm.vert(curr_ev.second), planeToInt(subm.refPlane()));
			// clang-format on
		}

		OMC_ASSERT(false, "can't find intersect vertex or edge.");
		return false;
	};

	/// @brief find the edge in link(seed) that intersect {A,B}
	/// @return true if found.
	auto binarySearchInterEdge = [&]() -> bool
	{
		const Sign orientation = subm.Orientation();
		const bool CCW         = orientation == Sign::POSITIVE;

		std::vector<index_t>  visited_tri(link_size, InvalidIndex);
		std::vector<IdxPair>  visited_ev(link_size,
		                                 IdxPair{InvalidIndex, InvalidIndex});
		std::vector<SignPair> calculated_ori(
		  link_size, SignPair{Sign::UNCERTAIN, Sign::UNCERTAIN});

		index_t  curr_tri = subm.adjV2T(v_start)[0];
		IdxPair  curr_ev{subm.nextVertInTri(curr_tri, v_start),
                    subm.prevVertInTri(curr_tri, v_start)};
		SignPair curr_ori{
		  OrientOn2D()(subm.vert(v_start), subm.vert(v_stop),
		               subm.vert(curr_ev.first), planeToInt(subm.refPlane())),
		  OrientOn2D()(subm.vert(v_start), subm.vert(v_stop),
		               subm.vert(curr_ev.second), planeToInt(subm.refPlane()))};

		// walk around v_start to find intersected vertex/edge.
		size_t local_idx = 0;
		size_t step      = link_size / 2;
		while (true)
		{
			if (curr_ori.first == Sign::ZERO &&
			    pointInsideSegmentCollinear(subm, v_start, v_stop, curr_ev.first))
			{ // meet the collinear point, check if it is in the constraint segment
				firstIntersectVertex(curr_ev.first);
				return false;
			}
			if (curr_ori.second == Sign::ZERO &&
			    pointInsideSegmentCollinear(subm, v_start, v_stop, curr_ev.second))
			{ // meet the collinear point, check if it is in the constraint segment
				firstIntersectVertex(curr_ev.second);
				return false;
			}

			if (curr_ori.first != curr_ori.second)
			{ // meet edge that possibly crosses the contraint segment
				if (curr_ori.first == reverse_sign(orientation) &&
				    curr_ori.second == orientation)
				{
					intersected_edges.push_back(
					  subm.edgeID(curr_ev.first, curr_ev.second));
					intersected_tris.push_back(curr_tri);
					return true;
				}
			}

			bool step_dir;
			if (curr_ori.first <= Sign::ZERO && curr_ori.second <= Sign::ZERO)
				step_dir = CCW; // curr_tri is at negative side
			else if (curr_ori.first >= Sign::ZERO && curr_ori.second >= Sign::ZERO)
				step_dir = !CCW; // curr_tri is at positive side
			else
				step_dir = CCW; // curr_tri overlaps both sides, either direction is ok

			visited_tri[local_idx]    = curr_tri;
			visited_ev[local_idx]     = curr_ev;
			calculated_ori[local_idx] = curr_ori;

			// gonna step to next triangle
			curr_tri  = subm.rotateAroundVertex(v_start, curr_tri, step, step_dir);
			curr_ev   = IdxPair{subm.nextVertInTri(curr_tri, v_start),
                        subm.prevVertInTri(curr_tri, v_start)};
			curr_ori  = SignPair{Sign::UNCERTAIN, Sign::UNCERTAIN};
			local_idx = CCW ? (local_idx + step) % link_size
			                : (local_idx + link_size - step) % link_size;
			step      = std::max(size_t(1), step / 2);
			index_t next_lid = (local_idx + 1) % link_size;
			index_t last_lid = (local_idx + link_size - 1) % link_size;
			// clang-format off
			if (is_valid_idx(visited_tri[next_lid]))
			{
				if      (curr_ev.first == visited_ev[next_lid].first)  curr_ori.first = calculated_ori[next_lid].first;
				else if (curr_ev.first == visited_ev[next_lid].second) curr_ori.first = calculated_ori[next_lid].second;

				if      (curr_ev.second == visited_ev[next_lid].first)  curr_ori.second = calculated_ori[next_lid].first;
				else if (curr_ev.second == visited_ev[next_lid].second) curr_ori.second = calculated_ori[next_lid].second;
			}
			if (is_valid_idx(visited_tri[last_lid]))
			{
				if      (curr_ev.first == visited_ev[last_lid].first)  curr_ori.first = calculated_ori[last_lid].first;
				else if (curr_ev.first == visited_ev[last_lid].second) curr_ori.first = calculated_ori[last_lid].second;

				if      (curr_ev.second == visited_ev[last_lid].first)  curr_ori.second = calculated_ori[last_lid].first;
				else if (curr_ev.second == visited_ev[last_lid].second) curr_ori.second = calculated_ori[last_lid].second;
			}
			if (curr_ori.first == Sign::UNCERTAIN)
				curr_ori.first = OrientOn2D()(subm.vert(v_start), subm.vert(v_stop), subm.vert(curr_ev.first), planeToInt(subm.refPlane()));
			if (curr_ori.second == Sign::UNCERTAIN)
				curr_ori.second = OrientOn2D()(subm.vert(v_start), subm.vert(v_stop), subm.vert(curr_ev.second), planeToInt(subm.refPlane()));
			// clang-format on
		}
	};

	if (link_size < 16)
	{
		if (!sequencialSearchInterEdge())
			return;
	}
	else
	{
		if (!binarySearchInterEdge())
			return;
	}

	OMC_ASSERT(intersected_edges.size() > 0, "empty intersected edges");

	/*******************************************************************/
	/*          Search remaining intersecting elements                 */
	/*******************************************************************/

	auto intersectVertex = [&](index_t v_inter, index_t t_id)
	{
		// current segment (v_start, v_stop) is split in (v_start, v_inter)
		// and (v_inter, v_stop).

		// split the original edge.
		splitSegmentInSubSegments(v_start, v_stop, v_inter, sub_segs_map);
		// put (v_inter, v_stop) in the segment_list to check later
		// it is a sub-segment so we set its index to invalid index.
		segment_list.emplace_back(InvalidIndex, uniquePair(v_inter, v_stop));
		// check if current segment encounters a TPI point
		if (/*isTPI*/ subm.vertFlag(v_inter))
		{ // if v_inter is a TPI point, fix indices stored in related segments.
			const RefSegs &seg_ids    = sub_segs_map.at(uniquePair(v_start, v_stop));
			std::vector<index_t> &t2s = tpi2segs.at(v_inter);
			for (index_t seg_id : seg_ids)
			{
				if (std::find(t2s.begin(), t2s.end(), seg_id) != t2s.end())
					continue;
				IPoint_TPI *tpi_inter = (IPoint_TPI *)&subm.vert(v_inter);
				for (index_t other_seg_id : t2s)
				{
					index_t suitable_vid = fixTPI(seg_id, other_seg_id, tpi_inter);
					if (suitable_vid < subm.vertInfo(v_inter))
						subm.setVertInfo(v_inter, suitable_vid);
				}
				t2s.push_back(seg_id);
			}
		}
		// output found edges/tris intersected with (v_start, v_inter).
		// adjust v_start and v_stop
		v_stop = v_inter;
		intersected_tris.push_back(t_id);
	};

	auto intersectEdge = [&](index_t ev0, index_t ev1, index_t t_id)
	{
		index_t e_id = subm.edgeID(ev0, ev1);
		OMC_EXPENSIVE_ASSERT(is_valid_idx(e_id), "invalid edge");
		intersected_edges.push_back(e_id);
		intersected_tris.push_back(t_id);
	};

	// walk along the topology to find the sorted list of edges and tris that
	// intersect {v_start, v_stop}
	Sign    orientation = subm.Orientation();
	index_t e_id        = intersected_edges.back();
	index_t ev0_id      = subm.edgeVertID(e_id, 0);
	index_t ev1_id      = subm.edgeVertID(e_id, 1);
	while (!subm.edgeIsConstr(e_id))
	{
		index_t t_id = subm.triOppToEdge(e_id, intersected_tris.back());
		OMC_EXPENSIVE_ASSERT(t_id != InvalidIndex, "invalid triangle index");
		index_t v2 = subm.triVertOppositeTo(t_id, ev0_id, ev1_id);
		if (v2 == v_stop)
		{ // converge.
			// append the last triangle (v_stop is in the last triangle)
			intersected_tris.push_back(t_id);
			OMC_EXPENSIVE_ASSERT(subm.triContainsVert(t_id, v_start) ||
			                       subm.triContainsVert(t_id, v_stop),
			                     "");
			return;
		}
		// walk to next triangle.
		Sign v2_ori = OrientOn2D()(subm.vert(v_start), subm.vert(v_stop),
		                           subm.vert(v2), planeToInt(subm.refPlane()));
		if (v2_ori == Sign::ZERO)
		{ // meet a vertex
			intersectVertex(v2, t_id);
			return;
		}
		else
		{ // meet an edge
			if (v2_ori == reverse_sign(orientation))
				intersectEdge(v2, subm.nextVertInTri(t_id, v2), t_id);
			else
				intersectEdge(v2, subm.prevVertInTri(t_id, v2), t_id);
		}
		// walk to next edge.
		e_id   = intersected_edges.back();
		ev0_id = subm.edgeVertID(e_id, 0);
		ev1_id = subm.edgeVertID(e_id, 1);
	}

	/*******************************************************************/
	/*  Intersect constrained segment, create TPI and split segments   */
	/*******************************************************************/

	// now, e_id is a constraint edge already present in the triangulation.
	// two constraint edges form a TPI point.
	const RefSegs &seg0_ids = sub_segs_map.at(uniquePair(v_start, v_stop));
	const RefSegs &seg1_ids = sub_segs_map.at(uniquePair(ev0_id, ev1_id));
	index_t        seg0_id  = *seg0_ids.begin();
	index_t        seg1_id  = *seg1_ids.begin();

	// TPI creation
	index_t global_tpi_id = createTPI(subm, seg0_id, seg1_id);

	// adding the TPI in the new_mesh
	index_t local_tpi_id = subm.addVert(&ts.vert(global_tpi_id), global_tpi_id);
	subm.setVertFlag(local_tpi_id, /*isTPI*/ true);

	// split edge
	subm.splitEdge(e_id, local_tpi_id);
	// set splitted edge as constraint
	index_t edge0_id = subm.edgeID(ev0_id, local_tpi_id);
	index_t edge1_id = subm.edgeID(local_tpi_id, ev1_id);
	OMC_EXPENSIVE_ASSERT(is_valid_idx(edge0_id), "invalid edge");
	OMC_EXPENSIVE_ASSERT(is_valid_idx(edge1_id), "invalid edge");
	subm.setEdgeConstr(edge0_id);
	subm.setEdgeConstr(edge1_id);

	// add two segments to tpi2segs
	OMC_EXPENSIVE_ASSERT(tpi2segs.find(local_tpi_id) == tpi2segs.end(),
	                     "the new tpi point has existed.");
	std::vector<index_t> &t2s = tpi2segs[local_tpi_id];
	t2s.insert(t2s.end(), seg0_ids.begin(), seg0_ids.end());
	t2s.insert(t2s.end(), seg1_ids.begin(), seg1_ids.end());

	if (t2s.size() > 2)
	{
		IPoint_TPI *new_tpi = (IPoint_TPI *)&subm.vert(local_tpi_id);
		for (index_t i = 0; i < t2s.size(); i++)
		{
			index_t seg_id = t2s[i];
			for (index_t j = i + 1; j < t2s.size(); j++)
			{
				index_t other_seg_id = t2s[j];
				index_t suitable_vid = fixTPI(seg_id, other_seg_id, new_tpi);
				if (suitable_vid < subm.vertInfo(local_tpi_id))
					subm.setVertInfo(local_tpi_id, suitable_vid);
			}
		}
	}

	// split two segments.
	// (note: seg0_ids and seg1_ids invalid after sub_segs_map changed)
	splitSegmentInSubSegments(v_start, v_stop, local_tpi_id, sub_segs_map);
	splitSegmentInSubSegments(ev0_id, ev1_id, local_tpi_id, sub_segs_map);

	// the original edge (v_start, v_stop) is split in (v_start, new_tpi) and
	// (new_tpi, v_stop).

	// put (new_tpi, v_stop) in the segment_list to check later
	segment_list.emplace_back(InvalidIndex, uniquePair(local_tpi_id, v_stop));

	// output found edges/tris intersected with (v_start, new_ip).
	if (intersected_tris.size() == 1)
	{
		// set (v_start, new_tpi) as constraint.
		index_t edge_id = subm.edgeID(v_start, local_tpi_id);
		OMC_EXPENSIVE_ASSERT(is_valid_idx(edge_id), "invalid edge");
		subm.setEdgeConstr(edge_id);
		// if (v_start, v_stop) only intersect one triangle,
		// after splitting, (v_start, new_ip) is an edge in subm.
		intersected_edges.clear();
		intersected_tris.clear();
	}
	else
	{
		// the last triangle and edge are splitted by new_ip,
		// we pop them, then...
		intersected_tris.pop_back();
		intersected_edges.pop_back();
		// ... push sub-triangle.
		index_t last_eid = intersected_edges.back();
		index_t last_tid = intersected_tris.back();
		index_t t_id     = subm.triOppToEdge(last_eid, last_tid);
		OMC_EXPENSIVE_ASSERT(is_valid_idx(t_id), "tri opposite to edge not found");
		intersected_tris.push_back(t_id);
		// adjust v_start and v_stop
		v_stop = local_tpi_id;
	}
}

template <typename Traits>
void Triangulation<Traits>::splitSegmentInSubSegments(index_t    v_start,
                                                      index_t    v_stop,
                                                      index_t    mid_point,
                                                      SubSegMap &sub_segs_map)
{
	UIPair orig_seg = uniquePair(v_start, v_stop);
	UIPair sub_seg0 = uniquePair(v_start, mid_point);
	UIPair sub_seg1 = uniquePair(v_stop, mid_point);

	// create two map before spliting and inserting
	sub_segs_map[sub_seg0];
	sub_segs_map[sub_seg1];

	const RefSegs &ref_seg_ids = sub_segs_map.find(orig_seg)->second;
	sub_segs_map.at(sub_seg0).insert(ref_seg_ids.begin(), ref_seg_ids.end());
	sub_segs_map.at(sub_seg1).insert(ref_seg_ids.begin(), ref_seg_ids.end());
}

template <typename Traits>
index_t Triangulation<Traits>::createTPI(FastTriMesh &subm, index_t seg0_id,
                                         index_t seg1_id)
{
	index_t t_id = subm.meshInfo();

	// find non-coplanar three planes (triangles) to create a TPI point.
	std::array<index_t, 3> planes;
	planes[0] = t_id;

	const concurrent_vector<index_t> &copl_tris = ts.coplanarTriangles(t_id);
	const concurrent_vector<index_t> &e0_s2t = ts.segmentTrianglesList(seg0_id);
	const concurrent_vector<index_t> &e1_s2t = ts.segmentTrianglesList(seg1_id);

	for (index_t ct_id : copl_tris)
	{
		if (std::binary_search(e0_s2t.begin(), e0_s2t.end(), ct_id) &&
		    std::binary_search(e1_s2t.begin(), e1_s2t.end(), ct_id))
		{
			planes[0] = ct_id;
			break;
		}
	}
	if (planes[0] > t_id)
		planes[0] = t_id;
	for (index_t st_id : e0_s2t)
	{
		if (st_id == t_id)
			continue;
		if (!std::binary_search(copl_tris.begin(), copl_tris.end(), st_id))
		{
			planes[1] = st_id;
			break;
		}
	}
	for (index_t st_id : e1_s2t)
	{
		if (st_id == t_id)
			continue;
		if (!std::binary_search(copl_tris.begin(), copl_tris.end(), st_id))
		{
			planes[2] = st_id;
			break;
		}
	}

	std::sort(planes.begin(), planes.end());

	auto [t0, t1, t2] = planes;

	std::array<const GPoint *, 3> t0_pnts = {
	  &ts.triVert(t0, 0), &ts.triVert(t0, 1), &ts.triVert(t0, 2)};
	std::array<const GPoint *, 3> t1_pnts = {
	  &ts.triVert(t1, 0), &ts.triVert(t1, 1), &ts.triVert(t1, 2)};
	std::array<const GPoint *, 3> t2_pnts = {
	  &ts.triVert(t2, 0), &ts.triVert(t2, 1), &ts.triVert(t2, 2)};

	int cur_thread_idx = tbb::this_task_arena::current_thread_index();

	IPoint_TPI *new_v = pnt_arenas[cur_thread_idx].emplace(
	  CreateTPI()(*t0_pnts[0], *t0_pnts[1], *t0_pnts[2], *t1_pnts[0], *t1_pnts[1],
	              *t1_pnts[2], *t2_pnts[0], *t2_pnts[1], *t2_pnts[2]));

	std::atomic<index_t> *new_idx =
	  idx_arenas[cur_thread_idx].emplace(InvalidIndex);

	auto [new_vertex_created, new_vid] =
	  addAndFixTPI(seg0_id, seg1_id, new_v, CreateIndex(ts, new_v, new_idx));

	if (!new_vertex_created)
	{
		IPoint_TPI::gcv().remove(new_v);
		pnt_arenas[cur_thread_idx].recycle(new_v);
		idx_arenas[cur_thread_idx].recycle(new_idx);
	}

	return new_vid;
}

template <typename Traits>
std::pair<bool, index_t>
Triangulation<Traits>::addAndFixTPI(index_t seg0_id, index_t seg1_id,
                                    IPoint_TPI *vtx, CreateIndex create_index)
{
	bool    new_vertex_created = false;
	index_t suitable_vid = InvalidIndex, fix_vid = InvalidIndex;

	index_t          min_segid = std::min(seg0_id, seg1_id);
	index_t          max_segid = std::max(seg0_id, seg1_id);
	tbb::spin_mutex &min_mutex = ts.getS2PMutex(min_segid);
	tbb::spin_mutex &max_mutex = ts.getS2PMutex(max_segid);

	// lock until add of fix end.
	std::lock_guard<tbb::spin_mutex> min_lock(min_mutex);
	std::lock_guard<tbb::spin_mutex> max_lock(max_mutex);

	index_t found_vid_in_e0 = ts.findVertexInSeg(seg0_id, AsGP()(*vtx));
	index_t found_vid_in_e1 = ts.findVertexInSeg(seg1_id, AsGP()(*vtx));

	if (is_valid_idx(found_vid_in_e0) && is_valid_idx(found_vid_in_e1))
	{
		// get the suitable point with simpler type and smaller index
		if (found_vid_in_e0 == found_vid_in_e1)
		{
			suitable_vid = found_vid_in_e0;
		}
		else
		{
			suitable_vid    = std::min(found_vid_in_e0, found_vid_in_e1);
			fix_vid         = std::max(found_vid_in_e0, found_vid_in_e1);
			index_t fix_seg = fix_vid == found_vid_in_e0 ? seg0_id : seg1_id;
			ts.fixVertexInSeg(fix_seg, fix_vid, suitable_vid);
		}
	}
	else if (is_valid_idx(found_vid_in_e0))
	{
		suitable_vid = found_vid_in_e0;
		ts.addVertexInSeg(seg1_id, found_vid_in_e0);
	}
	else if (is_valid_idx(found_vid_in_e1))
	{
		suitable_vid = found_vid_in_e1;
		ts.addVertexInSeg(seg0_id, found_vid_in_e1);
	}
	else
	{
		suitable_vid = create_index();
		ts.addVertexInSeg(seg0_id, suitable_vid);
		ts.addVertexInSeg(seg1_id, suitable_vid);
		new_vertex_created = true;
	}

	return std::pair<bool, index_t>(new_vertex_created, suitable_vid);
}

template <typename Traits>
index_t Triangulation<Traits>::fixTPI(index_t seg0_id, index_t seg1_id,
                                      IPoint_TPI *vtx)
{
	index_t suitable_vid = InvalidIndex, fix_vid = InvalidIndex;

	index_t          min_segid = std::min(seg0_id, seg1_id);
	index_t          max_segid = std::max(seg0_id, seg1_id);
	tbb::spin_mutex &min_mutex = ts.getS2PMutex(min_segid);
	tbb::spin_mutex &max_mutex = ts.getS2PMutex(max_segid);

	// lock until add of fix end.
	std::lock_guard<tbb::spin_mutex> min_lock(min_mutex);
	std::lock_guard<tbb::spin_mutex> max_lock(max_mutex);

	index_t found_vid_in_e0 = ts.findVertexInSeg(seg0_id, AsGP()(*vtx));
	index_t found_vid_in_e1 = ts.findVertexInSeg(seg1_id, AsGP()(*vtx));

	if (is_valid_idx(found_vid_in_e0) && is_valid_idx(found_vid_in_e1))
	{
		// get the suitable point with simpler type and smaller index
		if (found_vid_in_e0 != found_vid_in_e1)
		{
			suitable_vid    = std::min(found_vid_in_e0, found_vid_in_e1);
			fix_vid         = std::max(found_vid_in_e0, found_vid_in_e1);
			index_t fix_seg = fix_vid == found_vid_in_e0 ? seg0_id : seg1_id;
			ts.fixVertexInSeg(fix_seg, fix_vid, suitable_vid);
		}
		else
			suitable_vid = found_vid_in_e0; // any one is ok
	}
	else if (is_valid_idx(found_vid_in_e0))
	{
		suitable_vid = found_vid_in_e0;
		ts.addVertexInSeg(seg1_id, found_vid_in_e0);
	}
	else if (is_valid_idx(found_vid_in_e1))
	{
		suitable_vid = found_vid_in_e1;
		ts.addVertexInSeg(seg0_id, found_vid_in_e1);
	}
	else
	{
		OMC_ASSERT(false, "should not happen");
	}
	return suitable_vid;
}

template <typename Traits>
template <typename tri_iterator, typename edge_iterator>
void Triangulation<Traits>::boundaryWalker(const FastTriMesh &subm,
                                           index_t v_start, index_t v_stop,
                                           tri_iterator          curr_p,
                                           edge_iterator         curr_e,
                                           std::vector<index_t> &h)
{
	h.clear();
	h.push_back(v_start);
	do
	{
		index_t curr_v = h.back();
		index_t next_v = subm.nextVertInTri(*curr_p, curr_v);

		while (subm.edgeID(curr_v, next_v) == *curr_e)
		{
			++curr_p;
			if (subm.triContainsVert(*curr_p, v_stop))
			{
				h.push_back(v_stop);
				return;
			}

			++curr_e;

			OMC_EXPENSIVE_ASSERT(*curr_p < subm.numTriangles(), "index out of range");
			OMC_EXPENSIVE_ASSERT(is_valid_idx(*curr_e), "invalid edge");

			next_v = subm.nextVertInTri(*curr_p, curr_v);

			OMC_EXPENSIVE_ASSERT(is_valid_idx(subm.edgeID(curr_v, next_v)),
			                     "invalid edge");
		}

		h.push_back(next_v);
		++curr_p;

		if (subm.triContainsVert(*curr_p, v_stop))
		{
			h.push_back(v_stop);
			return;
		}

		++curr_e;
	} while (h.back() != v_stop);
}

template <typename Traits>
void Triangulation<Traits>::earcutLinear(const FastTriMesh          &subm,
                                         const std::vector<index_t> &poly,
                                         std::vector<index_t>       &tris)
{
	OMC_EXPENSIVE_ASSERT(poly.size() >= 3, "no valid poly dimension");

	Sign orientation = subm.Orientation();

	// doubly linked list for fat polygon inspection
	size_t               size = poly.size();
	std::vector<index_t> prev(size);
	std::vector<index_t> next(size);
	std::iota(prev.begin(), prev.end(), -1);
	std::iota(next.begin(), next.end(), 1);
	prev.front() = size - 1;
	next.back()  = 0;

	// keep a list of the ears to be cut
	std::vector<index_t> ears;
	ears.reserve(size);

	// this always has size |poly|, and keeps track of ears
	// (corners that were not ears at the beginning may become so later on)
	std::vector<bool> is_ear(size, false);

	// detect all safe ears in O(n).
	// This amounts to finding all convex vertices but the endpoints of the
	// constrained edge
	for (index_t curr = 1; curr < size - 1; ++curr)
	{
		// NOTE: the polygon may contain danging edges, prev!=next
		// avoids to even do the ear test for them

		const GPoint &p0 = subm.vert(poly[prev[curr]]);
		const GPoint &p1 = subm.vert(poly[curr]);
		const GPoint &p2 = subm.vert(poly[next[curr]]);

		Sign check = OrientOn2D()(p0, p1, p2, planeToInt(subm.refPlane()));

		if ((prev != next) && (check == orientation && check != Sign::ZERO))
		{
			ears.emplace_back(curr);
			is_ear.at(curr) = true;
		}
	}

	// progressively delete all ears, also updating the data structure
	size_t length = size;
	while (true)
	{
		index_t curr = ears.back();
		ears.pop_back();

		// make new_tri
		tris.push_back(poly[prev[curr]]);
		tris.push_back(poly[curr]);
		tris.push_back(poly[next[curr]]);

		// exclude curr from the polygon, connecting prev and next
		next[prev[curr]] = next[curr];
		prev[next[curr]] = prev[curr];

		// last triangle?
		if (--length < 3)
			return;

		// check if prev and next have become new_ears
		if (!is_ear[prev[curr]] && prev[curr] != 0)
		{
			const GPoint &p0 = subm.vert(poly[prev[prev[curr]]]);
			const GPoint &p1 = subm.vert(poly[prev[curr]]);
			const GPoint &p2 = subm.vert(poly[next[curr]]);

			Sign check = OrientOn2D()(p0, p1, p2, planeToInt(subm.refPlane()));

			if ((prev[prev[curr]] != next[curr]) &&
			    (check != Sign::ZERO && check == orientation))
			{
				ears.emplace_back(prev[curr]);
				is_ear.at(prev[curr]) = true;
			}
		}

		if (!is_ear[next[curr]] && next[curr] < size - 1)
		{
			const GPoint &p0 = subm.vert(poly[prev[curr]]);
			const GPoint &p1 = subm.vert(poly[next[curr]]);
			const GPoint &p2 = subm.vert(poly[next[next[curr]]]);

			Sign check = OrientOn2D()(p0, p1, p2, planeToInt(subm.refPlane()));

			if ((next[next[curr]] != prev[curr]) &&
			    (check != Sign::ZERO && check == orientation))
			{
				ears.emplace_back(next[curr]);
				is_ear.at(next[curr]) = true;
			}
		}
	}
}

template <typename Traits>
void Triangulation<Traits>::solvePocketsInCoplanarTriangle(
  const FastTriMesh &subm, std::vector<index_t> &new_tris,
  std::vector<Label> &new_labels, const Label &label)
{
	std::vector<std::vector<index_t>> tri_pockets;
	std::vector<std::set<index_t>>    polygons;

	findPocketsInTriangle(subm, tri_pockets, polygons);
	OMC_EXPENSIVE_ASSERT(tri_pockets.size() == polygons.size(), "size mismatch");

	std::vector<index_t> curr_polygon;
	for (index_t p_id = 0; p_id < polygons.size(); p_id++)
	{
		curr_polygon.clear();
		bool polygon_contain_tpi = false;
		for (index_t p : polygons[p_id])
		{
			if (subm.vertFlag(p))
				polygon_contain_tpi = true;
			// conversion from new to original vertices ids
			curr_polygon.push_back(subm.vertInfo(p));
		}

		remove_duplicates(curr_polygon);

		{ // critical section
			std::lock_guard<tbb::spin_mutex> lock(ts.new_tris_mutex);
			index_t pos = ts.addVisitedPolygonPocket(curr_polygon, new_labels.size(),
			                                         polygon_contain_tpi);

			if (!is_valid_idx(pos)) // pocket not added yet
			{
				const std::vector<index_t> &tri_list = tri_pockets[p_id];
				for (index_t t : tri_list)
				{
					const index_t *tri = subm.tri(t);
					new_tris.push_back(subm.vertInfo(tri[0]));
					new_tris.push_back(subm.vertInfo(tri[1]));
					new_tris.push_back(subm.vertInfo(tri[2]));
					new_labels.push_back(label);
				}
			}
			else // pocket already present
			{
				size_t num_tris = curr_polygon.size() - 2;

				for (size_t i = 0; i < num_tris; i++)
					new_labels[pos + i] |= label;
			}
		}
	}
}

template <typename Traits>
void Triangulation<Traits>::findPocketsInTriangle(
  const FastTriMesh &subm, std::vector<std::vector<index_t>> &tri_pockets,
  std::vector<std::set<index_t>> &tri_polygons)
{
	std::vector<bool> visited(subm.numTriangles(), false);

	for (size_t t_seed = 0; t_seed < subm.numTriangles(); ++t_seed)
	{
		if (visited[t_seed])
			continue;

		std::vector<index_t> curr_tri_pocket;
		std::set<index_t>    curr_tri_poly;
		std::stack<index_t>  stack;
		stack.push(t_seed);

		while (!stack.empty())
		{
			index_t curr_t = stack.top();
			stack.pop();

			if (visited[curr_t])
				continue;

			visited[curr_t] = true;
			curr_tri_pocket.push_back(curr_t);

			const std::array<index_t, 3> &t2e = subm.adjT2E(curr_t);

			for (index_t e_id : t2e)
			{
				if (subm.edgeIsConstr(e_id) || subm.edgeIsBoundary(e_id))
				{
					curr_tri_poly.insert(subm.edgeVertID(e_id, 0));
					curr_tri_poly.insert(subm.edgeVertID(e_id, 1));
				}
				else
				{
					for (index_t t : subm.adjE2T(e_id))
					{
						if (t != curr_t)
							stack.push(t);
					}
				}
			}
		}
		tri_pockets.push_back(curr_tri_pocket);
		tri_polygons.push_back(curr_tri_poly);
	}
}

template <typename Traits>
void Triangulation<Traits>::postFixIndices(std::vector<index_t> &new_tris,
                                           std::vector<Label>   &new_labels)
{
	if (!ts.any_index_fixed)
		return;

	Logger::info("TPI indices fixed.");

	// Fix indices stored in new_tris.
	// In this stage, no duplicate triangle will be generated.
	tbb::parallel_for_each(new_tris.begin(), new_tris.end(),
	                       [this](index_t &vid)
	                       { vid = ts.indices[vid]->load(); });

	// Fix indices stored in pockets with TPI points.
	// In this stage, duplicate triangles may be found, mark them as deleted
	// and then remove them.

	// we traversal all stored pockets, fix indices in them, check if them are
	// duplicate, and remove duplicates.
	typename TriSoup::PocketsMap &old_pockets_map = ts.pocketsMapWithTPI();
	typename TriSoup::PocketsMap  new_pockets_map;
	std::vector<uint8_t>          new_tri_deleted(new_labels.size(), false);

	if (old_pockets_map.empty())
		return;
	Logger::info("TPI pockets fixed.");

	for (std::pair<std::vector<index_t>, index_t> pocket : old_pockets_map)
	{
		std::vector<index_t> &polygon  = pocket.first;
		index_t               tris_pos = pocket.second;

		// fix indices stored in old pockets
		bool any_index_fixed = false;
		for (index_t &vid : polygon)
		{
			index_t new_vid = ts.indices[vid]->load();
			if (vid != new_vid)
			{
				vid             = new_vid;
				any_index_fixed = true;
			}
		}
		// add new pocket to new pockets map
		auto poly_iter = new_pockets_map.insert(std::make_pair(polygon, tris_pos));

		if (poly_iter.second) // pocket not added yet
		{
			// do nothing
		}
		else // pocket already present
		{
			index_t existed_tris_pos = poly_iter.first->second;
			size_t  num_tris         = polygon.size() - 2;
			// merge labels
			for (size_t i = 0; i < num_tris; i++)
				new_labels[existed_tris_pos + i] |= new_labels[tris_pos + i];
			// mark tris deleted
			for (size_t i = 0; i < num_tris; i++)
				new_tri_deleted[tris_pos + i] = true;
		}
	}

	// garbage collection
	index_t front = 0, back = new_labels.size() - 1;
	while (true)
	{
		while (!new_tri_deleted[front] && front < back)
			front++;
		while (new_tri_deleted[back] && front < back)
			back--;
		if (front >= back)
			break;
		// swap front and back
		std::swap(new_tri_deleted[front], new_tri_deleted[back]);
		// back is not used after swapping, so just update front.
		new_tris[front * 3]     = new_tris[back * 3];
		new_tris[front * 3 + 1] = new_tris[back * 3 + 1];
		new_tris[front * 3 + 2] = new_tris[back * 3 + 2];
		new_labels[front]       = new_labels[back];
	}
	size_t new_size = new_tri_deleted[back] ? (back) : (back + 1);
	new_tris.resize(new_size * 3);
	new_labels.resize(new_size);
}

template <typename Traits>
bool Triangulation<Traits>::pointInsideSegmentCollinear(const FastTriMesh &subm,
                                                        index_t ev0_id,
                                                        index_t ev1_id,
                                                        index_t p_id)
{
	return Segment3_Point3_DoIntersect().in_segment_collinear(
	         subm.vert(ev0_id), subm.vert(ev1_id), subm.vert(p_id)) ==
	       PointInType::STRICTLY_INSIDE;
}

} // namespace OMC