#pragma once

#include "FastTriMesh.h"
#include "Tree.h"
#include "TriangleSoup.h"

namespace OMC {

template <typename Traits>
class Triangulation
{
private:
	using NT         = typename Traits::NT;
	using EPoint     = typename Traits::EPoint;
	using GPoint     = typename Traits::GPoint;
	using IPoint_SSI = typename Traits::IPoint_SSI;
	using IPoint_LPI = typename Traits::IPoint_LPI;
	using IPoint_TPI = typename Traits::IPoint_TPI;
	using AsGP       = typename Traits::AsGP;
	using AsEP       = typename Traits::AsEP;
	using ToEP       = typename Traits::ToEP;
	using CreateSSI  = typename Traits::CreateSSI;
	using CreateLPI  = typename Traits::CreateLPI;
	using CreateTPI  = typename Traits::CreateTPI;

	using Orient3D   = typename Traits::Orient3D;
	using OrientOn2D = typename Traits::OrientOn2D;
	using LessThan3D = typename Traits::LessThan3D;

	using Segment3_Point3_DoIntersect =
	  typename Traits::Segment3_Point3_DoIntersect;

	using FastTriMesh = FastTriMesh<Traits>;
	using PntArena    = PointArena<Traits>;
	using TriSoup     = TriangleSoup<Traits>;

	// Segment, containing seg_id and seg's endpoints
	using Segment      = std::pair<index_t, UIPair>;
	// Collect segments on a triangle
	using SegmentsList = std::vector<Segment>;
	// Segment will be split to sub-segments by TPI points,
	// map sub-segments to its original segment id.
	using RefSegs =
	  boost::container::flat_set<index_t, std::less<index_t>,
	                             boost::container::small_vector<index_t, 4>>;
	using SubSegMap = phmap::flat_hash_map<UIPair, RefSegs>;
	// Store segments adajcent to TPI points in a triangle.
	using TPI2Segs  = phmap::flat_hash_map<index_t, std::vector<index_t>>;

	struct CreateIndex;

public:
	Triangulation(TriSoup &_ts, std::vector<index_t> &new_tris,
	              std::vector<Label> &new_labels);

private:
	void triangulateSingleTriangle(index_t t_id, FastTriMesh &subm,
	                               std::vector<index_t> &new_tris,
	                               std::vector<Label>   &new_labels);

	/* Split triangle and edges by constraint points ****************************/

	void
	sortedVertexListAlongSegment(const typename TriSoup::Edge2PntsSet &point_list,
	                             index_t v0_id, index_t v1_id,
	                             std::vector<index_t> &out_point_list);

	void splitSingleTriangle(FastTriMesh                &subm,
	                         const std::vector<index_t> &points);

	void splitSingleTriangleWithTree(FastTriMesh                &subm,
	                                 const std::vector<index_t> &points);

	void splitSingleEdge(FastTriMesh &subm, index_t v0_id, index_t v1_id,
	                     std::vector<index_t> &points);

	std::pair<index_t, bool> locatePointInTree(const FastTriMesh &subm,
	                                           index_t            p_id,
	                                           const SplitTree   &tree);

	std::pair<index_t, bool> locatePointInTreeRecur(const FastTriMesh &subm,
	                                                const GPoint      &p,
	                                                const SplitTree   &tree,
	                                                index_t node_id, UIPair ev);

	/* Split triangle by contraint segments *************************************/

	void addConstraintSegmentsInSingleTriangle(FastTriMesh  &subm,
	                                           SegmentsList &segment_list);

	void addConstraintSegment(FastTriMesh &subm, const Segment &seg,
	                          SegmentsList &segment_list, SubSegMap &sub_segs_map,
	                          TPI2Segs &tpi2segs);

	void findIntersectingElements(FastTriMesh &subm, index_t &v_start,
	                              index_t              &v_stop,
	                              std::vector<index_t> &intersected_edges,
	                              std::vector<index_t> &intersected_tris,
	                              SegmentsList         &segment_list,
	                              SubSegMap &sub_segs_map, TPI2Segs &tpi2segs);

	void splitSegmentInSubSegments(index_t v_start, index_t v_stop,
	                               index_t mid_point, SubSegMap &sub_segs_map);

	index_t createTPI(FastTriMesh &subm, index_t seg0_id, index_t seg1_id);

	std::pair<bool, index_t> addAndFixTPI(index_t seg0_id, index_t seg1_id,
	                                      IPoint_TPI *vtx,
	                                      CreateIndex create_index);

	index_t fixTPI(index_t seg0_id, index_t seg1_id, IPoint_TPI *vtx);

	template <typename tri_iterator, typename edge_iterator>
	void boundaryWalker(const FastTriMesh &subm, index_t v_start, index_t v_stop,
	                    tri_iterator curr_p, edge_iterator curr_e,
	                    std::vector<index_t> &h);

	void earcutLinear(const FastTriMesh &subm, const std::vector<index_t> &poly,
	                  std::vector<index_t> &tris);

	/* Solve pockets ************************************************************/

	void solvePocketsInCoplanarTriangle(const FastTriMesh    &subm,
	                                    std::vector<index_t> &new_tris,
	                                    std::vector<Label>   &new_labels,
	                                    const Label          &label);

	void findPocketsInTriangle(const FastTriMesh                 &subm,
	                           std::vector<std::vector<index_t>> &tri_pockets,
	                           std::vector<std::set<index_t>>    &tri_polygons);

	/* Postfix indices **********************************************************/

	void postFixIndices(std::vector<index_t> &new_tris,
	                    std::vector<Label>   &new_labels);

private:
	static bool pointInsideSegmentCollinear(const FastTriMesh &subm,
	                                        index_t ev0_id, index_t ev1_id,
	                                        index_t p_id);

private:
	TriSoup &ts;

	std::vector<PntArena> &pnt_arenas;
	std::vector<IdxArena> &idx_arenas;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Triangulation.inl"
#endif