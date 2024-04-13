#pragma once

#include "AuxStructure.h"
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

	using Segment      = std::pair<index_t, UIPair>; // seg_id, seg
	using SegmentsList = std::vector<Segment>;
	using SubSegMap    = phmap::flat_hash_map<UIPair, index_t>; // seg->seg_id

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
	                          SegmentsList &segment_list,
	                          SubSegMap    &sub_segs_map);

	void findIntersectingElements(FastTriMesh &subm, index_t &v_start,
	                              index_t              &v_stop,
	                              std::vector<index_t> &intersected_edges,
	                              std::vector<index_t> &intersected_tris,
	                              SegmentsList         &segment_list,
	                              SubSegMap            &sub_segs_map);

	void splitSegmentInSubSegments(index_t v_start, index_t v_stop,
	                               index_t mid_point, SubSegMap &sub_segs_map);

	index_t createTPI(FastTriMesh &subm, const UIPair &e0, const UIPair &e1,
	                  const SubSegMap &sub_segs_map);

	std::array<const GPoint *, 3> computeTriangleOfSegment(FastTriMesh &subm,
	                                                       index_t      seg_id);

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