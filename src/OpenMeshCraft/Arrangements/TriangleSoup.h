#pragma once

#include "Utils.h"

// clang-format off
#include "OpenMeshCraft/Utils/DisableWarnings.h"
#include "parallel_hashmap/phmap.h"
#include "OpenMeshCraft/Utils/EnableWarnings.h"
// clang-format on

#include <bitset>
#include <map>
#include <vector>

namespace OMC {

template <typename Traits>
class TriangleSoup
{
public: /* Types **************************************************************/
	using NT     = typename Traits::NT;
	using GPoint = typename Traits::GPoint;

	using Orient3D           = typename Traits::Orient3D;
	using OrientOn2D         = typename Traits::OrientOn2D;
	using MaxCompInTriNormal = typename Traits::MaxCompInTriNormal;

	using Edge    = UIPair;
	using EdgeMap = phmap::flat_hash_map<Edge, index_t>;

	using PntArena = PointArena<Traits>;

public: /* Constructors *******************************************************/
	TriangleSoup() = default;

	~TriangleSoup() {}

	void initialize();

public: /* Size queries *******************************************************/
	size_t numVerts() const { return vertices.size(); }
	size_t numTris() const { return triangles.size() / 3; }
	size_t numEdges() const { return edges.size(); }

	size_t numOrigVertices() const { return num_orig_vtxs; }
	size_t numOrigTriangles() const { return num_orig_tris; }

public: /* Vertices ***********************************************************/
	const GPoint &vert(index_t v_id) const;

	const NT *vertPtr(index_t v_id) const;

	index_t addImplVert(GPoint *pp, std::atomic<index_t> *ip);

public: /* Edges **************************************************************/
	const Edge &edge(index_t e_id) const;

	index_t edgeID(index_t v0_id, index_t v1_id) const;

	const GPoint &edgeVert(index_t e_id, size_t off) const;

	const NT *edgeVertPtr(index_t e_id, size_t off) const;

	index_t edgeOppositeToVert(index_t t_id, index_t v_id) const;

	void addEdge(index_t v0_id, index_t v1_id);

	bool edgeContainsVert(index_t e_id, index_t v_id) const;

public: /* Triangles **********************************************************/
	const std::vector<index_t> &trisVector() const { return triangles; }

	const index_t *tri(index_t t_id) const;

	index_t triVertID(index_t t_id, size_t off) const;

	const GPoint &triVert(index_t t_id, size_t off) const;

	const NT *triVertPtr(index_t t_id, size_t off) const;

	index_t triEdgeID(index_t t_id, size_t off) const;

	Plane triPlane(index_t t_id) const;

	Sign triOrientation(index_t t_id) const;

	bool triContainsVert(index_t t_id, index_t v_id) const;

	bool triContainsEdge(const index_t t_id, index_t ev0_id,
	                     index_t ev1_id) const;

	Label triLabel(index_t t_id) const;

public:
	/// implement vertices
	tbb::concurrent_vector<GPoint *>               vertices;
	/// implement indices of vertices (used in v_map)
	tbb::concurrent_vector<std::atomic<index_t> *> indices;

	/// all explicit and implicit points
	std::vector<PntArena> *pnt_arenas = nullptr;
	/// local indices
	std::vector<IdxArena> *idx_arenas = nullptr;

	/// triangles and related data
	std::vector<index_t> triangles;
	std::vector<Label>   tri_labels;

protected:
	NT vertX(index_t v_id) const;
	NT vertY(index_t v_id) const;
	NT vertZ(index_t v_id) const;

private:
	/// edges (a pair of indices of end vertices)
	std::vector<Edge> edges;
	/// map edge to a unique index
	EdgeMap           edge_map;

	/// triangle planes
	std::vector<Plane> tri_planes;
	std::vector<Sign>  tri_orientations;

	size_t num_orig_vtxs;
	size_t num_orig_tris;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "TriangleSoup.inl"
#endif