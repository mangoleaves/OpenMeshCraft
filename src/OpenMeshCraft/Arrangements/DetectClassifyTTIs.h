#pragma once

#include "AuxStructure.h"

namespace OMC {

template <typename Traits>
class DetectClassifyTTIs
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

	DetectClassifyTTIs(TriSoup &_ts, AuxStruct &_g,
	                            bool _parallel, MeshArrangements_Stats *_stats,
	                            bool _verbose);

	void checkTriangleTriangleIntersections();

	void propagateCoplanarTrianglesIntersections();

protected:
	void mergeConcurrentAuxStructures();

	void sortEdgePointsList();

	bool pointInsideTriangle(index_t p_id, index_t t_id);

protected:
	TriSoup                   &ts;
	std::vector<PntArena>     &pnt_arenas;
	std::vector<IdxArena>     &idx_arenas;
	AuxStruct                 &uniq_g;

	/// @brief concurrent_g is used to store data in multi-threading context.
	/// Final it will be merged in uniq_g.
	std::vector<AuxStruct>     concurrent_g;

	bool                    parallel;
	bool                    verbose;
	MeshArrangements_Stats *stats;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "DetectClassifyTTIs.inl"
#endif