#pragma once

#include "Utils.h"

#include "OpenMeshCraft/Utils/IndexDef.h"

namespace OMC {

enum class Sign;

class SplitTree;

#if 0	// OpenMeshCraft::InlinedVector, friendly for debug
template <typename T>
using AuxVector16 = InlinedVector<T, 16>;
template <typename T>
using AuxVector4 = InlinedVector<T, 4>;
#else	// absl::InlinedVector, possibly faster?
template <typename T>
using AuxVector16 = absl::InlinedVector<T, 16>;
template <typename T>
using AuxVector4 = absl::InlinedVector<T, 4>;
#endif

template <typename Traits>
class FastTriMesh
{
public:
	using Self = FastTriMesh<Traits>;

	using GPoint = typename Traits::GPoint;

	using Orient3D   = typename Traits::Orient3D;
	using OrientOn2D = typename Traits::OrientOn2D;

public:
	struct iVertex
	{
		iVertex() = delete;
		explicit iVertex(const GPoint *p, index_t id = 0) noexcept
		  : point(p)
		  , info(id)
		{
		}

		const GPoint *point;
		index_t       info;
	};

	struct iEdge
	{
		iEdge() = delete;
		explicit iEdge(index_t _v0, index_t _v1, const bool _b = false) noexcept
		  : verts(uniquePair(_v0, _v1))
		  , constr(_b)
		{
		}

		UIPair verts;
		bool   constr;

		bool operator<(const iEdge &rhs) const { return this->verts < rhs.verts; }
		bool operator<=(const iEdge &rhs) const { return this->verts <= rhs.verts; }
		bool operator==(const iEdge &rhs) const { return this->verts == rhs.verts; }
	};

	struct iTri
	{
		iTri() = delete;
		iTri(index_t _v0, index_t _v1, index_t _v2, index_t n = InvalidIndex)
		  : verts({_v0, _v1, _v2})
		  , info(n)
		{
		}

		std::array<index_t, 3> verts;
		std::array<index_t, 3> edges;
		index_t                info;
		Label                  label;
	};

	index_t mesh_info;

public: /* Constructors, Copy and Move ***************************************/
	FastTriMesh();

	FastTriMesh(const std::vector<GPoint *> &in_verts,
	            const std::vector<index_t>  &in_tris);

	FastTriMesh(const Self &rhs) noexcept { operator=(rhs); }

	FastTriMesh(Self &&rhs) noexcept { operator=(std::move(rhs)); }

	Self &operator=(Self &rhs) noexcept;

	Self &operator=(Self &&rhs) noexcept;

	~FastTriMesh();

	void initialize(const GPoint *tv0, const GPoint *tv1, const GPoint *tv2,
	                const index_t *tv_id, const Plane &ref_p, const Sign &ori);

	void preAllocateSpace(size_t estimated_num_verts);

	void clear();

public: /* Size queries ******************************************************/
	size_t numVerts() const;
	size_t numEdges() const;
	size_t numTriangles() const;

public: /* Info **************************************************************/
	// get index
	void resetVerticesInfo();
	void resetTrianglesInfo();

	index_t vertInfo(const index_t v_id) const;
	void    setVertInfo(const index_t v_id, const index_t info);

	index_t triInfo(index_t t_id) const;
	void    setTriInfo(index_t t_id, index_t val);

	index_t meshInfo() const;
	void    setMeshInfo(index_t val);

	const Label &triLabel(index_t t_id) const;
	void         setTriLabel(index_t t_id, const Label &label);

	index_t vertNewID(index_t orig_v_id) const;

public: /* Ajdacencies *******************************************************/
	bool edgeContainsVert(index_t e_id, index_t v_id) const;

	bool triContainsVert(index_t t_id, index_t v_id) const;

	const AuxVector4<index_t> &adjE2T(index_t e_id) const;

	const AuxVector16<index_t> &adjV2T(index_t v_id) const;

	const AuxVector16<index_t> &adjV2E(index_t v_id) const;

	const std::array<index_t, 3> &adjT2E(index_t t_id) const;

	index_t triVertOppositeTo(index_t t_id, index_t v0_id, index_t v1_id) const;

	index_t edgeOppToVert(index_t t_id, index_t v_id) const;

	index_t triOppToEdge(index_t e_id, index_t t_id) const;

	index_t nextVertInTri(index_t t_id, index_t curr_v_id) const;

	index_t prevVertInTri(index_t t_id, index_t curr_v_id) const;

public: /* Add / Remove ******************************************************/
	index_t addVert(const GPoint *p, index_t orig_v_id);

	index_t addEdge(index_t ev0_id, index_t ev1_id);

	index_t addTri(index_t tv0_id, index_t tv1_id, index_t tv2_id);

	void removeTri(index_t t_id);

	template <typename Container>
	void removeTris(const Container &t_ids);

	void removeEdge(index_t e_id);

public: /* Split **************************************************************/
	void splitEdge(const index_t e_id, index_t v_id);

	void splitEdge(const index_t e_id, index_t v_id, SplitTree &tree);

	void splitTri(index_t t_id, index_t v_id);

	void splitTri(index_t t_id, index_t v_id, SplitTree &tree);

	void flipTri(index_t t_id);

public: /* Primitive queries *************************************************/
	const GPoint &vert(index_t v_id) const;

	const GPoint &triVert(index_t t_id, size_t off) const;

	index_t triVertOffset(index_t t_id, index_t v_id) const;

	const index_t *tri(index_t t_id) const;

	index_t edgeID(index_t ev0_id, index_t ev1_id) const;

	index_t triID(index_t tv0_id, index_t tv1_id, index_t tv2_id) const;

	index_t triVertID(index_t t_id, size_t off) const;

	index_t triEdgeID(index_t t_id, index_t off) const;

	index_t edgeVertID(index_t e_id, index_t off) const;

public: /* Property queries **************************************************/
	Sign triOrientation(index_t t_id) const;

	void setEdgeConstr(index_t e_id);

	bool edgeIsConstr(index_t e_id) const;

	bool edgeIsBoundary(index_t e_id) const;

	bool edgeIsManifold(index_t e_id) const;

	size_t vertValence(index_t v_id) const;

	Plane refPlane() const;

	Sign Orientation() const;

public: /* Traversal meshes *************************************************/
	/* NOTE Below traversal functions are undefined if mesh is not manifold ***/

	index_t rotateAroundVertex(index_t center_v_id, index_t start_t_id,
	                           size_t step, bool CCW);

private: /* Auxiliary remove ************************************************/
	template <typename Container>
	void removeFromVec(Container &vec, index_t elem);
	template <typename Container>
	void updateInVec(Container &vec, index_t old_id, index_t new_id);

	void removeEdgeUnref(index_t e_id);
	void removeTriUnref(index_t t_id);

private:
	std::vector<iVertex> vertices;
	std::vector<iEdge>   edges;
	std::vector<iTri>    triangles;

	std::vector<AuxVector16<index_t>> v2e;
	std::vector<AuxVector16<index_t>> v2t;
	std::vector<AuxVector4<index_t>>  e2t;

	phmap::flat_hash_map<index_t, index_t> rev_vtx_map;

	Plane triangle_plane;
	Sign  orientation;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "FastTriMesh.inl"
#endif