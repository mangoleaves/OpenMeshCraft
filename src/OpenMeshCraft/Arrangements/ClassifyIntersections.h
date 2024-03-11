#pragma once

#include "AuxStructure.h"

namespace OMC {

template <typename Traits>
class ClassifyIntersections
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

	// clang-format off
	using Segment3_Point3_DoIntersect    = typename Traits::Segment3_Point3_DoIntersect;
	using Segment3_Segment3_DoIntersect  = typename Traits::Segment3_Segment3_DoIntersect;
	using Triangle3_Segment3_DoIntersect = typename Traits::Triangle3_Segment3_DoIntersect;
	using Triangle3_Point3_DoIntersect   = typename Traits::Triangle3_Point3_DoIntersect;
	// clang-format on

	using PntArena  = PointArena<Traits>;
	using TriSoup   = TriangleSoup<Traits>;
	using AuxStruct = AuxiliaryStructure<Traits>;

	ClassifyIntersections(TriSoup &_ts, AuxStruct &_g,
	                      const std::vector<ShewchukCachePtr> &_cache,
	                      bool _parallel, MeshArrangements_Stats *_stats,
	                      bool _verbose);

	void checkTriangleTriangleIntersections();

	void propagateCoplanarTrianglesIntersections();

protected:
	struct ThreadContext
	{
		ThreadContext(PntArena &_pa, IdxArena &_ia, AuxStruct &_g)
		  : pnt_arena(_pa)
		  , idx_arena(_ia)
		  , g(_g)
		{
		}

		PntArena  &pnt_arena;
		IdxArena  &idx_arena;
		AuxStruct &g;
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
	void checkTriTriInter(index_t tA_id, index_t tB_id, ThreadContext &tc);

	bool checkSingleCoplanarEdgeIntersections(index_t e_v0, index_t e_v1,
	                                          index_t e_t_id, index_t o_t_id,
	                                          phmap::flat_hash_set<size_t> &li,
	                                          ThreadContext                &tc);

	void checkSingleNoCoplanarEdgeIntersection(
	  index_t e_id, index_t t_id, phmap::flat_hash_set<size_t> &v_tmp,
	  phmap::flat_hash_set<size_t> &li, ThreadContext &tc);

	void checkVtxInTriangleIntersection(index_t v_id, index_t t_id,
	                                    phmap::flat_hash_set<size_t> &lv,
	                                    phmap::flat_hash_set<size_t> &li,
	                                    ThreadContext                &tc);

	void addSymbolicSegment(index_t v0_id, index_t v1_id, index_t tA_id,
	                        index_t tB_id, ThreadContext &tc);

	index_t addEdgeCrossCoplanarEdgeInters(index_t e0_id, index_t e1_id,
	                                       index_t t_id, ThreadContext &tc);

	index_t addEdgeCrossNoCoplanarEdgeInters(index_t e0_id, index_t e1_id,
	                                         index_t t_id, ThreadContext &tc);

	index_t addEdgeCrossTriInters(index_t e_id, index_t t_id, ThreadContext &tc);

	void mergeConcurrentAuxStructures();

	void sortEdgePointsList();

	bool pointInsideTriangle(index_t p_id, index_t t_id);

protected:
	TriSoup                   &ts;
	std::vector<PntArena>     &pnt_arenas;
	std::vector<IdxArena>     &idx_arenas;
	/// @brief used to store all aux data.
	AuxStruct                 &uniq_g;
	/// @brief concurrent_g is used to store data in multi-threading context.
	/// Final it will be merged in uniq_g.
	std::vector<AuxStruct>     concurrent_g;
	std::vector<ThreadContext> thread_contexts;

	/// cached data from DetectIntersections
	const std::vector<ShewchukCachePtr> &cache;

	bool                    parallel;
	bool                    verbose;
	MeshArrangements_Stats *stats;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ClassifyIntersections.inl"
#endif