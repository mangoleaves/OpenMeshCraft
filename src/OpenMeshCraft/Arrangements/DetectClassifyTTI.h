/**
 * @file DetectClassifyTTI.h
 * @brief Detect and Classify Triangle-Triangle-Intersection in one pass.
 */
#pragma once

#include "AuxStructure.h"
#include "TriangleSoup.h"
#include "Utils.h"

namespace OMC {

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

	using PntArena  = PointArena<Traits>;
	using TriSoup   = TriangleSoup<Traits>;
	using AuxStruct = AuxiliaryStructure<Traits>;

	DetectClassifyTTI(TriSoup &_ts, AuxStruct &_uniq_g, AuxStruct &_g,
	                  PntArena &_pnt_arena, IdxArena &_idx_arena);

public:
	void check_TTI(index_t ta, index_t tb);

protected:
	/* Input data */
	TriSoup   &ts;
	AuxStruct &uniq_g;
	AuxStruct &g;
	PntArena  &pnt_arena;
	IdxArena  &idx_arena;

	struct TTIHelper
	{
		static const index_t UncachedIndex   = (index_t)(-2);
		static const char    UncachedBoolean = (char)(-1);

		static bool is_cached(index_t idx) { return idx != UncachedIndex; }
		static bool is_cached(char boolean) { return boolean != UncachedBoolean; }

		index_t t_id;
		int     t_nmax = -1;

		// indices of vertices from triangle.
		std::array<index_t, 3>    v_id;
		// indices of edges from triangle.
		std::array<index_t, 3>    e_id;
		// pointers to vertices from triangle.
		// TODO store coordinates for better cache performance?
		std::array<const NT *, 3> v;

		// position of sub-simplex of this triangle with respect to sub-simplex of
		// another triangle
		std::array<index_t, 3> v_in_vtx = {InvalidIndex, InvalidIndex,
		                                   InvalidIndex};
		std::array<index_t, 3> v_in_seg = {UncachedIndex, UncachedIndex,
		                                   UncachedIndex};
		std::array<char, 3>    v_in_tri = {UncachedBoolean, UncachedBoolean,
		                                   UncachedBoolean};

		// cached orientations
		// Sign::UNCERTAIN means "not cached/calculated"
		std::array<Sign, 3> v_wrt_tri = {Sign::UNCERTAIN, Sign::UNCERTAIN,
		                                 Sign::UNCERTAIN};
		std::array<std::array<Sign, 3>, 3> v_wrt_seg = {
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN})};
		std::array<std::array<Sign, 3>, 3> seg_wrt_seg = {
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN}),
		  std::array<Sign, 3>({Sign::UNCERTAIN, Sign::UNCERTAIN, Sign::UNCERTAIN})};
	};

	struct CoplanarEEI // Edge-Edge-Intersection
	{
		index_t ea; // global index of edge a
		index_t eb; // global index of edge b
		index_t p;  // global index of the intersection point

		CoplanarEEI(index_t _ea, index_t _eb, index_t _p)
		  : ea(_ea)
		  , eb(_eb)
		  , p(_p)
		{
		}

		// used in hash set
		bool operator==(const CoplanarEEI &rhs) const { return p == rhs.p; }
		// used to check if two intersection points are same.
		bool is_same(index_t query_ea, index_t query_eb) const
		{
			return (query_ea == ea && query_eb == eb) ||
			       (query_ea == eb && query_eb == ea);
		}
	};

	struct Hasher
	{
		index_t operator()(const CoplanarEEI &c) const
		{
			return std::hash<index_t>()(c.p);
		}
	};

	class CreateIndex
	{
	public:
		CreateIndex(TriSoup &_ts, AuxStruct &_g)
		  : ts(_ts)
		  , g(_g)
		{
		}

		void operator()(const GPoint *pp, std::atomic<index_t> *ip)
		{
			index_t idx = InvalidIndex;
			{ // lock for new index
				std::lock_guard<tbb::spin_mutex> lock(g.new_vertex_mutex);
				idx = ts.addImplVert(const_cast<GPoint *>(pp), ip);
			}
			ip->store(idx, std::memory_order_relaxed); // assign a valid index
		}

	private:
		TriSoup   &ts;
		AuxStruct &g;
	};

protected:

	bool get_v_in_seg(TTIHelper &ha, index_t va, TTIHelper &hb, index_t eb);

	index_t get_v_in_seg(TTIHelper &ha, index_t va, TTIHelper &hb);

	bool get_v_in_tri(TTIHelper &ha, index_t va, TTIHelper &hb);

	Sign get_v_wrt_seg(TTIHelper &ha, index_t va, TTIHelper &hb, index_t eb);

	Sign get_v_wrt_tri(TTIHelper &ha, index_t va, TTIHelper &hb);

	Sign get_seg_wrt_seg(TTIHelper &ha, index_t ea, TTIHelper &hb, index_t eb);

	size_t get_vtx_wrt_sector(TTIHelper &ha, index_t va, TTIHelper &hb,
	                          index_t eb);

	bool seg_seg_do_intersect(TTIHelper &ha, index_t ea, TTIHelper &hb,
	                          index_t eb);

	bool coplanar_seg_tri_do_intersect(TTIHelper &ha, index_t ea, TTIHelper &hb);

	bool noncoplanar_seg_tri_do_intersect(TTIHelper &ha, index_t ea,
	                                      TTIHelper &hb);

protected:
	void check_TTI_share_edge(TTIHelper &ha, TTIHelper &hb);

	void check_TTI_share_vertex(TTIHelper &ha, TTIHelper &hb);

	void check_TTI_separate(TTIHelper &ha, TTIHelper &hb);

	bool classify_coplanr_vtx_intersections(
	  TTIHelper &ha, index_t va, TTIHelper &hb,
	  phmap::flat_hash_set<index_t> &inter_list,
	  phmap::flat_hash_set<index_t> &vtx_list);

	bool classify_coplanar_edge_intersections(
	  TTIHelper &ha, index_t ea, TTIHelper &hb,
	  phmap::flat_hash_set<index_t>             &inter_list,
	  phmap::flat_hash_set<CoplanarEEI, Hasher> &copl_inter_list);

	bool classify_noncoplanar_edge_intersections(
	  TTIHelper &ha, index_t ea, TTIHelper &hb,
	  phmap::flat_hash_set<index_t> &inter_list,
	  phmap::flat_hash_set<index_t> &vtx_list);

	void add_symbolic_segment(index_t v0, index_t v1, index_t ta, index_t tb);

	index_t add_edge_cross_coplanar_edge(
	  index_t ea, index_t eb, index_t t,
	  phmap::flat_hash_set<CoplanarEEI, Hasher> &copl_inter_list);

	index_t add_edge_cross_noncoplanar_edge(index_t ea, index_t eb);

	index_t add_edge_cross_tri(index_t ea, index_t tb);
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "DetectClassifyTTI.inl"
#endif