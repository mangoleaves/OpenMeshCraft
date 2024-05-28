#pragma once

#include "TriangleSoup.h"
#include "Utils.h"

namespace OMC {

/// @brief Detect and Classify Triangle-Triangle-Intersection in one pass.
template <typename Traits>
class DetectClassifyTTI
{
public:
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

	using Orient3D           = typename Traits::Orient3D;
	using CollinearPoints3D  = typename Traits::CollinearPoints3D;
	using OrientOn2D         = typename Traits::OrientOn2D;
	using LessThan3D         = typename Traits::LessThan3D;
	using MaxCompInTriNormal = typename Traits::MaxCompInTriNormal;

	using PntArena = PointArena<Traits>;
	using TriSoup  = TriangleSoup<Traits>;

	DetectClassifyTTI(TriSoup &_ts, PntArena &_pnt_arena);

public:
	void check_TTI(index_t ta, index_t tb);

protected:
	/* Input data */
	TriSoup  &ts;
	PntArena &pnt_arena;

	struct TTIHelper;
	struct CoplanarEEI;

	using IntersectionPoints = AuxVector4<index_t>;
	using IntersectionTypes  = AuxVector4<PointInSimplexType>;
	using CoplanarEEIList    = AuxVector4<CoplanarEEI>;

protected:
	index_t get_e_id(TTIHelper &ha, index_t ea);

	bool get_v_in_seg(TTIHelper &ha, index_t va, TTIHelper &hb, index_t eb);

	index_t get_v_in_seg(TTIHelper &ha, index_t va, TTIHelper &hb);

	bool get_v_in_tri(TTIHelper &ha, index_t va, TTIHelper &hb);

	Sign get_v_wrt_seg(TTIHelper &ha, index_t va, TTIHelper &hb, index_t eb);

	Sign get_v_wrt_tri(TTIHelper &ha, index_t va, TTIHelper &hb);

	Sign get_seg_wrt_seg(TTIHelper &ha, index_t ea, TTIHelper &hb, index_t eb);

	size_t get_vtx_wrt_sector(TTIHelper &ha, index_t va, TTIHelper &hb,
	                          index_t eb);

	size_t find_vtx_correspondence(TTIHelper &ha, TTIHelper &hb);

	bool fast_check_on2d_share_vertex(TTIHelper &ha, index_t ea, TTIHelper &hb,
	                                  index_t eb);

	bool fast_check_on2d_separate(TTIHelper &ha, TTIHelper &hb);

	bool seg_seg_do_intersect(TTIHelper &ha, index_t ea, TTIHelper &hb,
	                          index_t eb, Sign eb0_wrt_ea, Sign eb1_wrt_ea);

	bool coplanar_seg_tri_do_intersect(TTIHelper &ha, index_t ea, TTIHelper &hb);

	bool noncoplanar_seg_tri_do_intersect(TTIHelper &ha, index_t ea,
	                                      TTIHelper &hb);

	bool intersection_on_one_edge(const IntersectionTypes &intersection_types,
	                              index_t                 &edge_id);

protected:
	void check_TTI_share_edge(TTIHelper &ha, TTIHelper &hb);

	void check_TTI_share_vertex(TTIHelper &ha, TTIHelper &hb);

	void check_TTI_separate(TTIHelper &ha, TTIHelper &hb);

	bool
	classify_coplanr_vtx_intersections(TTIHelper &ha, index_t va, TTIHelper &hb,
	                                   IntersectionPoints &intersection_points,
	                                   IntersectionTypes  &intersection_types);

	bool classify_coplanar_edge_intersections(
	  TTIHelper &ha, index_t ea, TTIHelper &hb,
	  CoplanarEEIList *copl_edge_crosses = nullptr);

	bool classify_noncoplanar_edge_intersections(
	  TTIHelper &ha, index_t ea, TTIHelper &hb,
	  IntersectionPoints &intersection_points,
	  IntersectionTypes  &intersection_types);

	void add_symbolic_segment(index_t v0, index_t v1, TTIHelper &ha,
	                          bool is_part_of_ea, TTIHelper &hb,
	                          bool is_part_of_eb);

	index_t add_vertex_in_tri(TTIHelper &ha, TTIHelper &hb, index_t vb);

	index_t add_vertex_in_edge(TTIHelper &ha, index_t ea, TTIHelper &hb,
	                           index_t vb);

	index_t add_edge_cross_coplanar_edge(TTIHelper &ha, index_t ea, TTIHelper &hb,
	                                     index_t          eb,
	                                     CoplanarEEIList *copl_edge_crosses);

	index_t add_edge_cross_noncoplanar_edge(TTIHelper &ha, index_t ea,
	                                        TTIHelper &hb, index_t eb);

	index_t add_edge_cross_tri(TTIHelper &ha, index_t ea, TTIHelper &hb);

	OMC_NODISCARD std::pair<index_t, bool> add_SSI(index_t ea_id, index_t eb_id,
	                                               GPoint *new_v);

	OMC_NODISCARD std::pair<index_t, bool> add_LPI(index_t e_id, index_t t_id,
	                                               GPoint *new_v);
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "DetectClassifyTTI.inl"
#endif