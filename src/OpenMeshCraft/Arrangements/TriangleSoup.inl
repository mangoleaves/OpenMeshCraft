#pragma once

#include "TriangleSoup.h"

namespace OMC {

/**
 * @brief Points on a coplanar edge are inside a coplanar triangle if and only
 * if they are between two intersection points.
 */
template <typename Traits>
struct TriangleSoup<Traits>::CCrEdgeInfo
{
	index_t e_id;  // index of the coplanar/colinear edge
	index_t v0_id; // index of one intersected point
	index_t v1_id; // index of another intersected point

	CCrEdgeInfo()  = default;
	~CCrEdgeInfo() = default;

	CCrEdgeInfo(index_t _e_id, index_t _v0_id, index_t _v1_id)
	  : e_id(_e_id)
	  , v0_id(_v0_id)
	  , v1_id(_v1_id)
	{
	}

	bool operator<(const CCrEdgeInfo &rhs) const { return e_id < rhs.e_id; }
	bool operator==(const CCrEdgeInfo &rhs) const { return e_id == rhs.e_id; }
};

template <typename Traits>
struct TriangleSoup<Traits>::EdgeComparator
{
	using is_transparent = std::true_type;
	// retrieve points in triangle soup.
	const TriangleSoup<Traits> &ts;
	// compare points on specific axis.
	uint32_t                    axis;

	EdgeComparator()  = default;
	~EdgeComparator() = default;

	EdgeComparator(TriangleSoup<Traits> &_ts, uint32_t _axis)
	  : ts(_ts)
	  , axis(_axis)
	{
	}

	bool operator()(index_t lhs, index_t rhs) const
	{
		return operator()(ts.vert(lhs), ts.vert(rhs));
	}
	bool operator()(index_t lhs, const GPoint &rhs) const
	{
		return operator()(ts.vert(lhs), rhs);
	}
	bool operator()(const GPoint &lhs, index_t rhs) const
	{
		return operator()(lhs, ts.vert(rhs));
	}
	bool operator()(const GPoint &lhs, const GPoint &rhs) const
	{
		if (!lhs.is_Explicit() && lhs.point_type() == rhs.point_type())
		{ // lhs and rhs are both implicit points, compare topology first.
			if (lhs.is_SSI())
			{
				const IPoint_SSI &_lhs = lhs.SSI(), &_rhs = rhs.SSI();
				if (&_lhs.A() == &_rhs.A() && &_lhs.B() == &_rhs.B() &&
				    &_lhs.P() == &_rhs.P() && &_lhs.Q() == &_rhs.Q())
					return false;
			}
			else if (lhs.is_LPI())
			{
				const IPoint_LPI &_lhs = lhs.LPI(), &_rhs = rhs.LPI();
				if (&_lhs.P() == &_rhs.P() && &_lhs.Q() == &_rhs.Q() &&
				    &_lhs.R() == &_rhs.R() && &_lhs.S() == &_rhs.S() &&
				    &_lhs.T() == &_rhs.T())
					return false;
			}
		}
		return LessThan3D().on(lhs, rhs, axis) == Sign::NEGATIVE;
	}
};

template <typename Traits>
void TriangleSoup<Traits>::initialize()
{
	num_orig_vtxs = static_cast<size_t>(vertices.size());
	num_orig_tris = static_cast<size_t>(triangles.size() / 3);

	edges.reserve(numOrigVerts());
	edge_map.resize(numOrigVerts());
	edge_mutexes = std::vector<tbb::spin_mutex>(numOrigVerts());

	tri_edges.resize(numOrigTris() * 3, InvalidIndex);

	any_index_fixed = false;

	tri_has_intersections.resize(numOrigTris(), false);
	coplanar_tris.resize(numOrigTris());
	coplanar_edges.resize(numOrigTris());

	tri2pts.resize(numOrigTris());
	edge2pts.reserve(numOrigVerts());
	tri2segs.resize(numOrigTris());

	int concurrency    = tbb::this_task_arena::max_concurrency();
	int num_hash_table = concurrency == 1 ? 1 : concurrency * 4;
	seg2tris.resize(num_hash_table);
	seg_mutexes = std::vector<tbb::spin_mutex>(num_hash_table);
}

template <typename Traits>
auto TriangleSoup<Traits>::vert(index_t v_id) const -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(v_id < numVerts(), "vtx id out of range");
	return *vertices[v_id];
}

template <typename Traits>
auto TriangleSoup<Traits>::vertPtr(index_t v_id) const -> const NT *
{
	OMC_EXPENSIVE_ASSERT(v_id < num_orig_vtxs,
	                     "vtx id out of range of original points");
	return vertices[v_id]->data();
}

template <typename Traits>
index_t TriangleSoup<Traits>::addImplVert(GPoint *pp, std::atomic<index_t> *ip)
{
	vertices.push_back(pp);
	indices.push_back(ip);
	return vertices.size() - 1;
}

template <typename Traits>
const UIPair &TriangleSoup<Traits>::edge(index_t e_id) const
{
	OMC_EXPENSIVE_ASSERT(e_id < numEdges(), "edge index out of range.");
	return edges[e_id];
}

template <typename Traits>
index_t TriangleSoup<Traits>::getOrAddEdge(index_t v0_id, index_t v1_id)
{
	Edge edge = uniquePair(v0_id, v1_id);

	phmap::flat_hash_map<index_t, index_t> &em = edge_map[edge.first];

	// lock until finding and/or adding operation end.
	std::lock_guard<tbb::spin_mutex> lock_edge_map(edge_mutexes[edge.first]);

	auto find_iter = em.find(edge.second);
	if (find_iter != em.end()) // edge exists, just return it.
		return find_iter->second;
	else // edge does not exist, add it
	{
		// 1. add a new edge and related instances.
		typename decltype(edges)::iterator    edges_iter;
		typename decltype(edge2pts)::iterator e2p_iter;

		// -- calculate the longest axis of the edge
		const NT *v0 = vertPtr(v0_id);
		const NT *v1 = vertPtr(v1_id);

		NT diff_x = fabs(v0[0] - v1[0]);
		NT diff_y = fabs(v0[1] - v1[1]);
		NT diff_z = fabs(v0[2] - v1[2]);

		uint32_t longest_axis =
		  diff_x > diff_y ? (diff_x > diff_z ? 0 : 2) : (diff_y > diff_z ? 1 : 2);

		OMC_EXPENSIVE_ASSERT(v0[longest_axis] - v1[longest_axis] != 0.,
		                     "degenerate edge.");

		// -- lock and add
		{
			// lock until all adding operations end.
			std::lock_guard<tbb::spin_mutex> lock_add_edge(new_edge_mutex);

#ifdef OMC_ARR_TS_PARA
			edges_iter = edges.push_back(edge);
			e2p_iter   = edge2pts.emplace_back(EdgeComparator(*this, longest_axis));
#else
			edges.push_back(edge);
			edge2pts.emplace_back(EdgeComparator(*this, longest_axis));
			edges_iter = edges.end() - 1;
			e2p_iter   = edge2pts.end() - 1;
#endif

			OMC_EXPENSIVE_ASSERT(e2p_iter - edge2pts.begin() ==
			                       edges_iter - edges.begin(),
			                     "size mismatch.");

			edge2pts_mutex.emplace_back();
			colinear_edges.emplace_back();
		}
		// 2. get new edge's index
		index_t edge_id = edges_iter - edges.begin();

		// 3. build map between new edge and its index
		auto insert_iter = em.insert({edge.second, edge_id});
		OMC_EXPENSIVE_ASSERT(insert_iter.second, "fail to add edge.");
		return edge_id;
	}
}

template <typename Traits>
const index_t *TriangleSoup<Traits>::tri(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return &triangles[3 * t_id];
}

template <typename Traits>
index_t TriangleSoup<Traits>::triVertID(index_t t_id, index_t off) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return triangles[3 * t_id + off];
}

template <typename Traits>
auto TriangleSoup<Traits>::triVert(index_t t_id,
                                   index_t off) const -> const GPoint &
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vert(triangles[3 * t_id + off]);
}

template <typename Traits>
auto TriangleSoup<Traits>::triVertPtr(index_t t_id,
                                      index_t off) const -> const NT *
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return vertPtr(triangles[3 * t_id + off]);
}

template <typename Traits>
index_t TriangleSoup<Traits>::triEdgeID(index_t t_id, index_t off)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	if (is_valid_idx(tri_edges[3 * t_id + off]))
		// edge index is cached, return it.
		return tri_edges[3 * t_id + off];
	// otherwise check if the edge exists
	index_t v0_id = triangles[3 * t_id + off];
	index_t v1_id = triangles[3 * t_id + ((off + 1) % 3)];
	Edge    edge  = uniquePair(v0_id, v1_id);
	if (edge_map[edge.first].find(edge.second) == edge_map[edge.first].end())
		// the edge does not exist, return invalid index.
		return InvalidIndex;
	// otherwise get the index of edge and return it.
	tri_edges[3 * t_id + off] = edge_map[edge.first].at(edge.second);
	return tri_edges[3 * t_id + off];
}

template <typename Traits>
index_t TriangleSoup<Traits>::triEdgeID(index_t t_id, index_t off,
                                        std::true_type)
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	// edge index is cached, return it.
	if (is_valid_idx(tri_edges[3 * t_id + off]))
		return tri_edges[3 * t_id + off];
	// otherwise get or add the edge
	index_t e_id              = getOrAddEdge(triangles[3 * t_id + off],
	                                         triangles[3 * t_id + ((off + 1) % 3)]);
	// get the index and return it
	tri_edges[3 * t_id + off] = e_id;
	return e_id;
}

template <typename Traits>
Label TriangleSoup<Traits>::triLabel(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numTris(), "t_id out of range");
	return tri_labels[t_id];
}

template <typename Traits>
void TriangleSoup<Traits>::buildVMap(Tree *tree)
{
	tree->clear_points();
	v_map = std::make_unique<AuxPointMap_Tree<Traits>>(tree);
}

template <typename Traits>
void TriangleSoup<Traits>::addVertexInTriangle(index_t t_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2pts.size(), "out of range");
	concurrent_vector<index_t> &points = tri2pts[t_id];
	if (points.empty())
		points.reserve(8);
	points.push_back(v_id);
}

template <typename Traits>
tbb::spin_mutex &TriangleSoup<Traits>::getE2PMutex(index_t e_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts_mutex.size(), "out of range");
	return edge2pts_mutex[e_id];
}

/**
 * @brief Find if a point in edge is geometrically concident to given point.
 * @return index_t Index of the found vertex, InvalidIndex if not found.
 * @note Before calling this func, lock the edge by edge2pts_mutex.
 */
template <typename Traits>
index_t TriangleSoup<Traits>::findVertexInEdge(index_t       e_id,
                                               const GPoint &pnt) const
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	const Edge2PntsSet &points = edge2pts[e_id];
	auto                iter   = points.find(pnt);
	return iter != points.cend() ? *iter : InvalidIndex;
}

/**
 * @brief Add a point in edge.
 * @note Before calling this func, lock the edge by edge2pts_mutex.
 * @note Make sure that the point to add does **NOT** exist in the edge.
 */
template <typename Traits>
void TriangleSoup<Traits>::addVertexInEdge(index_t e_id, index_t v_id)
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	Edge2PntsSet &points = edge2pts[e_id];
	auto [iter, succeed] = points.insert(v_id);
	OMC_EXPENSIVE_ASSERT(succeed, "fail to insert vertex in edge.");
}

/**
 * @brief Fix a point in edge to another point.
 * @note Before calling this func, lock the edge by edge2pts_mutex.
 * @note Make sure that the point to add **EXISTS** in the edge.
 */
template <typename Traits>
void TriangleSoup<Traits>::fixVertexInEdge(index_t e_id, index_t old_vid,
                                           index_t new_vid)
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	// fix index in edge2pts
	Edge2PntsSet &points = edge2pts[e_id];
	auto          iter   = points.find(old_vid);
#if 0 // std::set
	auto          hint   = iter;
	++hint;
	points.erase(iter);
	points.insert(hint, new_vid);
#else // boost::container::flat_set
	*iter = new_vid;
#endif
	// fix global index
	indices[old_vid]->store(new_vid, std::memory_order_relaxed);
	// fixing indices in tri2pts, tri2segs, seg2tris is delayed.
	any_index_fixed = true;
}

template <typename Traits>
void TriangleSoup<Traits>::addSegmentInTriangle(index_t t_id, const UIPair &seg)
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2segs.size(), "out of range");
	OMC_EXPENSIVE_ASSERT(isUnique(seg), "segment is not unique");
	concurrent_vector<UIPair> &segments = tri2segs[t_id];
	if (segments.empty())
		segments.reserve(8);
	segments.push_back(seg);
}

template <typename Traits>
void TriangleSoup<Traits>::addTrianglesInSegment(const UIPair &seg,
                                                 index_t       t_id)
{
	OMC_EXPENSIVE_ASSERT(isUnique(seg), "segment is not unique");
	index_t outer_map = (seg.first + seg.second) % seg2tris.size();
	{
		std::lock_guard<tbb::spin_mutex> lock(seg_mutexes[outer_map]);
		std::vector<index_t>            &tris = seg2tris[outer_map][seg];
		if (tris.empty())
			tris.reserve(8);
		tris.push_back(t_id);
	}
}

template <typename Traits>
void TriangleSoup<Traits>::addTrianglesInSegment(const UIPair &seg,
                                                 index_t tA_id, index_t tB_id)
{
	OMC_EXPENSIVE_ASSERT(isUnique(seg), "segment is not unique");
	index_t outer_map = (seg.first + seg.second) % seg2tris.size();
	{
		std::lock_guard<tbb::spin_mutex> lock(seg_mutexes[outer_map]);
		std::vector<index_t>            &tris = seg2tris[outer_map][seg];
		if (tris.empty())
			tris.reserve(8);
		tris.push_back(tA_id);
		tris.push_back(tB_id);
	}
}

template <typename Traits>
void TriangleSoup<Traits>::addCoplanarTriangles(index_t ta, index_t tb)
{
	OMC_EXPENSIVE_ASSERT(ta != tb, "same triangles");
	OMC_EXPENSIVE_ASSERT(ta < coplanar_tris.size() && tb < coplanar_tris.size(),
	                     "out of range");
	if (coplanar_tris[ta].empty())
		coplanar_tris[ta].reserve(8);
	if (coplanar_tris[tb].empty())
		coplanar_tris[tb].reserve(8);
	coplanar_tris[ta].push_back(tb);
	coplanar_tris[tb].push_back(ta);
}

template <typename Traits>
void TriangleSoup<Traits>::addCoplanarEdge(index_t t_id, index_t e_id,
                                           index_t v0_id, index_t v1_id)
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_edges.size(), "out of range");

	if (coplanar_edges[t_id].empty())
		coplanar_edges[t_id].reserve(8);
	coplanar_edges[t_id].emplace_back(e_id, v0_id, v1_id);
}

template <typename Traits>
void TriangleSoup<Traits>::addColinearEdge(index_t e0_id, index_t e1_id,
                                           index_t v0_id, index_t v1_id)
{
	OMC_EXPENSIVE_ASSERT(e0_id < colinear_edges.size(), "out of range");
	OMC_EXPENSIVE_ASSERT(e1_id < colinear_edges.size(), "out of range");

	index_t min_eid = std::min(e0_id, e1_id);
	index_t max_eid = std::max(e0_id, e1_id);
	if (colinear_edges[min_eid].empty())
		colinear_edges[min_eid].reserve(4);
	colinear_edges[min_eid].emplace_back(max_eid, v0_id, v1_id);
}

template <typename Traits>
const concurrent_vector<index_t> &
TriangleSoup<Traits>::coplanarTriangles(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_tris.size(), "out of range");
	return coplanar_tris[t_id];
}

template <typename Traits>
auto TriangleSoup<Traits>::coplanarEdges(index_t t_id) const
  -> const concurrent_vector<CCrEdgeInfo> &
{
	OMC_EXPENSIVE_ASSERT(t_id < coplanar_edges.size(), "out of range");
	return coplanar_edges[t_id];
}

template <typename Traits>
auto TriangleSoup<Traits>::colinearEdges(index_t e_id) const
  -> const concurrent_vector<CCrEdgeInfo> &
{
	OMC_EXPENSIVE_ASSERT(e_id < colinear_edges.size(), "out of range");
	return colinear_edges[e_id];
}

template <typename Traits>
void TriangleSoup<Traits>::setTriangleHasIntersections(index_t tA_id,
                                                       index_t tB_id)
{
	OMC_EXPENSIVE_ASSERT(tA_id < tri_has_intersections.size(), "out of range");
	OMC_EXPENSIVE_ASSERT(tB_id < tri_has_intersections.size(), "out of range");
	tri_has_intersections[tA_id] = true;
	tri_has_intersections[tB_id] = true;
}

template <typename Traits>
bool TriangleSoup<Traits>::triangleHasIntersections(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri_has_intersections.size(), "out of range");
	return tri_has_intersections[t_id];
}

template <typename Traits>
const concurrent_vector<index_t> &
TriangleSoup<Traits>::trianglePointsList(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2pts.size(), "out of range");
	return tri2pts[t_id];
}

template <typename Traits>
auto TriangleSoup<Traits>::edgePointsList(index_t e_id) const
  -> const Edge2PntsSet &
{
	OMC_EXPENSIVE_ASSERT(e_id < edge2pts.size(), "out of range");
	return edge2pts[e_id];
}

template <typename Traits>
const concurrent_vector<UIPair> &
TriangleSoup<Traits>::triangleSegmentsList(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < tri2segs.size(), "out of range");
	return tri2segs[t_id];
}

template <typename Traits>
const std::vector<index_t> &
TriangleSoup<Traits>::segmentTrianglesList(const UIPair &seg) const
{
	index_t outer_map = (seg.first + seg.second) % seg2tris.size();
	auto    res       = seg2tris[outer_map].find(seg);
	OMC_EXPENSIVE_ASSERT(res != seg2tris[outer_map].end(), "out of range");
	return res->second;
}

template <typename Traits>
Plane TriangleSoup<Traits>::triPlane(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numOrigTris(), "out of range.");
	OMC_EXPENSIVE_ASSERT(planeToInt(tri_plane[t_id]) >= 0, "invalid plane");
	return tri_plane[t_id];
}

template <typename Traits>
Sign TriangleSoup<Traits>::triOrient(index_t t_id) const
{
	OMC_EXPENSIVE_ASSERT(t_id < numOrigTris(), "out of range.");
	OMC_EXPENSIVE_ASSERT(is_sign_reliable(tri_orient[t_id]), "invalid sign");
	return tri_orient[t_id];
}

template <typename Traits>
std::pair<index_t, bool>
TriangleSoup<Traits>::addVertexInSortedList(const GPoint         *pp,
                                            std::atomic<index_t> *ip)
{
	auto [insert_ip, insert_res] =
	  static_cast<AuxPointMap_Tree<Traits> *>(v_map.get())->insert(pp, ip);
	return std::pair<index_t, bool>(insert_ip->load(), insert_res);
}

template <typename Traits>
template <typename GetIndex>
std::pair<index_t, bool> TriangleSoup<Traits>::addVertexInSortedList(
  const GPoint *pp, std::atomic<index_t> *ip, GetIndex get_idx)
{
	return static_cast<AuxPointMap_Tree<Traits> *>(v_map.get())
	  ->insert_F(pp, ip, get_idx);
}

template <typename Traits>
index_t TriangleSoup<Traits>::addVisitedPolygonPocket(
  const std::vector<index_t> &polygon, index_t pos)
{
	auto poly_it = pockets_map.insert(std::make_pair(polygon, pos));

	if (poly_it.second)
		return InvalidIndex; // polygon not present yet

	return poly_it.first->second;
}

template <typename Traits>
void TriangleSoup<Traits>::removeDuplicatesBeforeFix()
{
	tbb::parallel_for(size_t(0), numOrigTris(), [this](size_t t_id)
	                  { remove_duplicates(coplanar_tris[t_id]); });

	tbb::parallel_for(size_t(0), numEdges(), [this](index_t e_id)
	                  { remove_duplicates(colinear_edges[e_id]); });
}

template <typename Traits>
void TriangleSoup<Traits>::fixColinearEdgesIntersections()
{
	auto fix_colinear_edge = [this](index_t e_id)
	{
		// mutex to lock current edge
		tbb::spin_mutex    &e_mutex = getE2PMutex(e_id);
		const Edge2PntsSet &e2p     = edgePointsList(e_id);

		const concurrent_vector<CCrEdgeInfo> &ccr_edge_infos = colinearEdges(e_id);
		for (const CCrEdgeInfo &edge_info : ccr_edge_infos)
		{
			index_t coli_e_id = edge_info.e_id;
			index_t v0_id     = edge_info.v0_id;
			index_t v1_id     = edge_info.v1_id;

			OMC_EXPENSIVE_ASSERT(e_id < coli_e_id, "edge indices error");

			tbb::spin_mutex &coli_e_mutex = getE2PMutex(coli_e_id);

			// lock until fix end.
			std::lock_guard<tbb::spin_mutex> e_lock(e_mutex);
			std::lock_guard<tbb::spin_mutex> coli_e_lock(coli_e_mutex);

			bool overlap = false;
			for (index_t p_id : e2p)
			{
				if (!overlap)
				{
					if (p_id == v0_id || p_id == v1_id) // step into overlap interval
						overlap = true;
					continue;
				}

				// already inside overlap interval
				if (p_id == v0_id || p_id == v1_id) // step out of overlap interval
					break; // end traversal current colianar edge
				// else, current p_id is inside the overlap interval
				// check if point is already in coli_e2p
				index_t found_vid = findVertexInEdge(coli_e_id, vert(p_id));
				if (!is_valid_idx(found_vid))
					continue; // can't find a coincident point, continue to next.
				if (p_id == found_vid)
					continue; // find a same point with same index, continue to next.
				// two coincident points have different indices, need fix.
				if (vert(p_id).point_type() < vert(found_vid).point_type())
				{
					fixVertexInEdge(coli_e_id, found_vid, p_id);
				}
				else if (vert(found_vid).point_type() < vert(p_id).point_type())
				{
					fixVertexInEdge(e_id, p_id, found_vid);
				}
				else if (p_id < found_vid) // same type, but p_id smaller
				{
					fixVertexInEdge(coli_e_id, found_vid, p_id);
				}
				else // same type, but found_vid smaller
				{
					fixVertexInEdge(e_id, p_id, found_vid);
				}
			} // end for e2p
		} // end for ccr_edge_infos
	}; // end fix_colinear_edge

#ifdef OMC_ARR_TS_PARA
	tbb::parallel_for(size_t(0), numEdges(), fix_colinear_edge);
#else
	for (size_t e_id = 0; e_id < numEdges(); e_id++)
		fix_colinear_edge(e_id);
#endif
}

template <typename Traits>
void TriangleSoup<Traits>::fixAllIndices()
{
	// 1. fix global indices
	auto fix_vertex_idx = [this](index_t orig_vidx)
	{
		index_t fix_vidx = orig_vidx;
		index_t new_vidx = indices[fix_vidx]->load(std::memory_order_relaxed);
		while (fix_vidx != new_vidx)
		{
			fix_vidx = new_vidx;
			new_vidx = indices[fix_vidx]->load(std::memory_order_relaxed);
		}
		indices[orig_vidx]->store(fix_vidx, std::memory_order_relaxed);
	};

	if (any_index_fixed)
		tbb::parallel_for(numOrigVerts(), numVerts(), fix_vertex_idx);

#ifdef OMC_ENABLE_EXPENSIVE_ASSERT
	tbb::parallel_for(
	  numOrigVerts(), numVerts(),
	  [this](index_t idx)
	  {
		  index_t fixed_idx = indices[idx]->load(std::memory_order_relaxed);

		  OMC_ASSERT(fixed_idx ==
		               indices[fixed_idx]->load(std::memory_order_relaxed),
		             "unfinished fix.");
		  OMC_ASSERT(LessThan3D()(vert(idx), vert(fixed_idx)) == Sign::ZERO,
		             "vertex and fixed vertex do not coincident.");
	  });
#endif

	// 2. fix indices stored in tri2pts
	auto fix_tri2pts = [this](concurrent_vector<index_t> &t2p)
	{
		if (any_index_fixed)
		{
			for (index_t &i : t2p)
				i = indices[i]->load(std::memory_order_relaxed);
		}
		remove_duplicates(t2p);
	};
	tbb::parallel_for_each(tri2pts.begin(), tri2pts.end(), fix_tri2pts);

	// 3. fix indices stored in tri2segs
	auto fix_tri2segs = [this](index_t t_id)
	{
		concurrent_vector<UIPair> &t2s = tri2segs[t_id];

		if (any_index_fixed)
		{
			for (UIPair &seg : t2s)
			{
				seg = uniquePair(indices[seg.first]->load(std::memory_order_relaxed),
				                 indices[seg.second]->load(std::memory_order_relaxed));
			}
		}

		remove_duplicates(t2s);

		for (const UIPair &seg : t2s)
			addTrianglesInSegment(seg, t_id);
	};
	tbb::parallel_for(size_t(0), numOrigTris(), fix_tri2segs);

	// 4. fix indices stored in coplanar_edges
	auto fix_coplanar_edges = [this](index_t t_id)
	{
		tbb::concurrent_vector<CCrEdgeInfo> &ce = coplanar_edges[t_id];

		remove_duplicates(ce);
		if (any_index_fixed)
		{
			for (CCrEdgeInfo &edge_info : ce)
			{
				edge_info.v0_id =
				  indices[edge_info.v0_id]->load(std::memory_order_relaxed);
				edge_info.v1_id =
				  indices[edge_info.v1_id]->load(std::memory_order_relaxed);
			}
		}
	};
	tbb::parallel_for(size_t(0), numOrigTris(), fix_coplanar_edges);

#if defined(OMC_ENABLE_EXPENSIVE_ASSERT)
	struct Comparator
	{
		TriangleSoup<Traits> &ts;

		Comparator(TriangleSoup<Traits> &_ts)
		  : ts(_ts)
		{
		}

		bool operator()(index_t v0, index_t v1) const
		{
			return LessThan3D()(ts.vert(v0), ts.vert(v1)) == Sign::NEGATIVE;
		}
	};
	auto checkDuplOnTris = [this](index_t t_id)
	{
		std::set<index_t, Comparator> point_set =
		  std::set<index_t, Comparator>(Comparator(*this));
		size_t cnt = 0;
		point_set.insert(triVertID(t_id, 0));
		point_set.insert(triVertID(t_id, 1));
		point_set.insert(triVertID(t_id, 2));
		cnt += 3;
		for (size_t j = 0; j < 3; j++)
		{
			if (is_valid_idx(triEdgeID(t_id, j)))
			{
				const auto &e2p = edgePointsList(triEdgeID(t_id, j));
				OMC_ASSERT(e2p.size() >= 2, "miss two endpoints.");
				point_set.insert(e2p.begin(), e2p.end());
				cnt += e2p.size() - 2;
			}
		}
		point_set.insert(trianglePointsList(t_id).begin(),
		                 trianglePointsList(t_id).end());
		cnt += trianglePointsList(t_id).size();
		OMC_ASSERT(point_set.size() == cnt, "duplicate vertex detected");
	};
	for (size_t t_id = 0; t_id < numOrigTris(); t_id++)
		checkDuplOnTris(t_id);
#endif
}

template <typename Traits>
void TriangleSoup<Traits>::addEndPointsToE2P()
{
	auto add_end_points_to_edge2pts = [this](size_t e_id)
	{
		const Edge &e    = edge(e_id);
		index_t     ev0  = e.first;
		index_t     ev1  = e.second;
		uint32_t    axis = edge2pts[e_id].key_comp().axis;
		if (vert(ev0)[axis] < vert(ev1)[axis])
		{
			edge2pts[e_id].insert(edge2pts[e_id].begin(), ev0);
			edge2pts[e_id].insert(edge2pts[e_id].end(), ev1);
		}
		else
		{
			edge2pts[e_id].insert(edge2pts[e_id].begin(), ev1);
			edge2pts[e_id].insert(edge2pts[e_id].end(), ev0);
		}
	};

	tbb::parallel_for(size_t(0), numEdges(), add_end_points_to_edge2pts);
}

template <typename Traits>
void TriangleSoup<Traits>::calcPlaneAndOrient()
{
	tri_plane.resize(numOrigTris(), intToPlane(-1));
	tri_orient.resize(numOrigTris(), Sign::UNCERTAIN);

	auto calc_plane_and_orient = [this](index_t t_id)
	{
		if (triangleHasIntersections(t_id))
		{
			index_t v0_id = triVertID(t_id, 0);
			index_t v1_id = triVertID(t_id, 1);
			index_t v2_id = triVertID(t_id, 2);

			tri_plane[t_id] = intToPlane(
			  MaxCompInTriNormal()(vertPtr(v0_id), vertPtr(v1_id), vertPtr(v2_id)));
			tri_orient[t_id] = OrientOn2D()(vertPtr(v0_id), vertPtr(v1_id),
			                                vertPtr(v2_id), tri_plane[t_id]);
		}
	};

	tbb::parallel_for(size_t(0), numOrigTris(), calc_plane_and_orient);
}

} // namespace OMC