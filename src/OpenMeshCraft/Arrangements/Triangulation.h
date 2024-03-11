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

	// clang-format onff
	using Segment3_Point3_DoIntersect =
	  typename Traits::Segment3_Point3_DoIntersect;
	using Segment3_Segment3_DoIntersect =
	  typename Traits::Segment3_Segment3_DoIntersect;
	using Triangle3_Point3_DoIntersect =
	  typename Traits::Triangle3_Point3_DoIntersect;
	using Triangle3_Segment3_DoIntersect =
	  typename Traits::Triangle3_Segment3_DoIntersect;
	// clang-format on

	using FastTriMesh = FastTriMesh<Traits>;
	using PntArena    = PointArena<Traits>;
	using TriSoup     = TriangleSoup<Traits>;
	using AuxStruct   = AuxiliaryStructure<Traits>;

public:
	Triangulation(TriSoup &_ts, AuxStruct &_g, std::vector<index_t> &new_tris,
	              std::vector<Label> &new_labels);

private:
	void triangulateSingleTriangle(index_t t_id, FastTriMesh &subm,
	                               std::vector<index_t> &new_tris,
	                               std::vector<Label>   &new_labels);

	/* Split triangle and edges by constraint points ****************************/

	void sortedVertexListAlongSegment(const std::vector<index_t> &point_list,
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

	void addConstraintSegmentsInSingleTriangle(FastTriMesh         &subm,
	                                           std::vector<UIPair> &segment_list);

	void addConstraintSegment(FastTriMesh &subm, const UIPair &seg,
	                          std::vector<UIPair>                  &segment_list,
	                          phmap::flat_hash_map<UIPair, UIPair> &sub_segs_map);

	void
	findIntersectingElements(FastTriMesh &subm, index_t &v_start, index_t &v_stop,
	                         std::vector<index_t> &intersected_edges,
	                         std::vector<index_t> &intersected_tris,
	                         std::vector<UIPair>  &segment_list,
	                         phmap::flat_hash_map<UIPair, UIPair> &sub_segs_map);

	void
	splitSegmentInSubSegments(index_t v_start, index_t v_stop, index_t mid_point,
	                          phmap::flat_hash_map<UIPair, UIPair> &sub_segs_map);

	index_t createTPI(FastTriMesh &subm, const UIPair &e0, const UIPair &e1,
	                  const phmap::flat_hash_map<UIPair, UIPair> &sub_segs_map);

	std::array<const GPoint *, 3> computeTriangleOfSegment(FastTriMesh  &subm,
	                                                       const UIPair &seg);

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
	static bool pointOnLine(const FastTriMesh &subm, index_t e_id, index_t p_id);

	static bool segmentsIntersectInside(const FastTriMesh &subm, index_t e00_id,
	                                    index_t e01_id, index_t e10_id,
	                                    index_t e11_id);

	static bool pointInsideSegment(const FastTriMesh &subm, index_t ev0_id,
	                               index_t ev1_id, index_t p_id);

	static bool pointInsideSegmentCollinear(const FastTriMesh &subm,
	                                        index_t ev0_id, index_t ev1_id,
	                                        index_t p_id);

private:
	TriSoup   &ts;
	AuxStruct &g;

	std::vector<PntArena> &pnt_arenas;
	std::vector<IdxArena> &idx_arenas;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Triangulation.inl"
#endif