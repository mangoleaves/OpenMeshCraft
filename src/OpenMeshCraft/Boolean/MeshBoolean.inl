#pragma once

#include "MeshBoolean.h"
// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
// Arrangements
#include "OpenMeshCraft/Arrangements/FastTriMesh.h"
#include "OpenMeshCraft/Arrangements/MeshArrangements.h"
#include "OpenMeshCraft/Arrangements/TriangleSoup.h"

namespace OMC {

template <typename Kernel, typename Traits>
class MeshBoolean<Kernel, Traits>::BooleanTraits
{
public:
	using K = Kernel;

	using NT         = typename K::NT;
	using EPoint     = typename K::EPoint3;
	using GPoint     = typename K::GPoint3;
	using IPoint_SSI = typename K::IPoint3T_SSI;
	using IPoint_LPI = typename K::IPoint3T_LPI;
	using IPoint_TPI = typename K::IPoint3T_TPI;

	using AsGP      = typename K::AsGP;
	using AsEP      = typename K::AsEP;
	using ToEP      = typename K::ToEP;
	using CreateSSI = typename K::CreateSSI3;
	using CreateLPI = typename K::CreateLPI;
	using CreateTPI = typename K::CreateTPI;

	using Segment     = typename K::Segment3;
	using Triangle    = typename K::Triangle3;
	using BoundingBox = typename K::BoundingBox3;

	using Orient2D           = typename K::Orient2D;
	using Orient3D           = typename K::Orient3D;
	using OrientOn2D         = typename K::OrientOn2D;
	using LessThan3D         = typename K::LessThan3D;
	using LongestAxis        = typename K::LongestAxis;
	using CollinearPoints3D  = typename K::CollinearPoints3D;
	using MaxCompInTriNormal = typename K::MaxCompInTriNormal;

	using CalcBbox = typename K::CalcBoundingBox3;

	// clang-format off
	using DoIntersect = typename K::DoIntersect;
	using Segment3_Point3_DoIntersect     = typename K::Segment3_Point3_DoIntersect;
	using Segment3_Segment3_DoIntersect   = typename K::Segment3_Segment3_DoIntersect;
	using Triangle3_Point3_DoIntersect    = typename K::Triangle3_Point3_DoIntersect;
	using Triangle3_Segment3_DoIntersect  = typename K::Triangle3_Segment3_DoIntersect;
	using Triangle3_Triangle3_DoIntersect = typename K::Triangle3_Triangle3_DoIntersect;
	// clang-format on
};

template <typename Traits>
class MeshBoolean_Impl
{
public: /* Traits ***********************************************************/
	// primitives
	using NT         = typename Traits::NT;
	using EPoint     = typename Traits::EPoint;
	using GPoint     = typename Traits::GPoint;
	using IPoint_SSI = typename Traits::IPoint_SSI;
	using IPoint_LPI = typename Traits::IPoint_LPI;
	using IPoint_TPI = typename Traits::IPoint_TPI;

	using AsGP      = typename Traits::AsGP;
	using AsEP      = typename Traits::AsEP;
	using ToEP      = typename Traits::ToEP;
	using CreateSSI = typename Traits::CreateSSI;
	using CreateLPI = typename Traits::CreateLPI;
	using CreateTPI = typename Traits::CreateTPI;

	using Segment  = typename Traits::Segment;
	using Triangle = typename Traits::Triangle;

	// predicates
	using Orient2D           = typename Traits::Orient2D;
	using Orient3D           = typename Traits::Orient3D;
	using OrientOn2D         = typename Traits::OrientOn2D;
	using LessThan3D         = typename Traits::LessThan3D;
	using LongestAxis        = typename Traits::LongestAxis;
	using CollinearPoints3D  = typename Traits::CollinearPoints3D;
	using MaxCompInTriNormal = typename Traits::MaxCompInTriNormal;

public: /* Auxiliary data structures *****************************************/
	// tree
	using Tree     = Arr_Tree_Intersection<Traits>;
	// fast triangle mesh
	using FTriMesh = FastTriMesh<Traits>;

	enum class IntersInfo
	{
		DISCARD,
		NO_INT,
		INT_IN_V0,
		INT_IN_V1,
		INT_IN_V2,
		INT_IN_EDGE01,
		INT_IN_EDGE12,
		INT_IN_EDGE20,
		INT_IN_TRI
	};

	struct Ray
	{
		Segment segment;
		index_t tv[3] = {InvalidIndex, InvalidIndex, InvalidIndex};
		char    dir   = 'X';
	};

#define LESS_THAN(V)                                                      \
	struct LessThanGPOn##V                                                  \
	{                                                                       \
		bool operator()(const std::pair<GPoint *, index_t> &p0,               \
		                const std::pair<GPoint *, index_t> &p1) const         \
		{                                                                     \
			return LessThan3D().on_##V(*p0.first, *p1.first) != Sign::POSITIVE; \
		}                                                                     \
		bool operator()(const std::pair<GPoint *, index_t> &p0,               \
		                const GPoint                       &p1) const                               \
		{                                                                     \
			return LessThan3D().on_##V(*p0.first, p1) == Sign::NEGATIVE;        \
		}                                                                     \
		bool operator()(const GPoint &p0, const GPoint &p1) const             \
		{                                                                     \
			return LessThan3D().on_##V(p0, p1) == Sign::NEGATIVE;               \
		}                                                                     \
	}

	LESS_THAN(x);
	LESS_THAN(y);
	LESS_THAN(z);

#undef LESS_THAN

public: /* Constructors ******************************************************/
	MeshBoolean_Impl(bool _verbose)
	  : arr_impl(nullptr)
	  // Get reference to data in arrangements
	  , in_coords(arr_impl.in_coords)
	  , in_tris(arr_impl.in_tris)
	  , in_labels(arr_impl.in_labels)
	  , arr_in_tris(arr_impl.arr_in_tris)
	  , arr_in_labels(arr_impl.arr_in_labels)
	  , arr_out_verts(arr_impl.arr_out_verts)
	  , arr_out_tris(arr_impl.arr_out_tris)
	  , arr_out_labels(arr_impl.arr_out_labels)
	  , tree(arr_impl.tree)
	  , dupl_triangles(arr_impl.dupl_triangles)
	  // behavior flags
	  , verbose(_verbose)
	{
	}

public: /* Pipelines ********************************************************/
	/// Complete boolean pipeline contains four steps:
	/// 1. Run arrangements on all triangles.
	/// 2. Label triangles.
	/// 3. apply boolean operation based on label.
	/// 4. compute explicit result and output.

	/// @brief This pipeline finish the first two steps of the boolean pipeline:
	/// 1. Run arrangements to all triangles.
	/// 2. Label each triangles.
	void computeLabelsPipeline();

	/// @brief This pipeline finish the third step of the boolean pipeline:
	/// 3. apply boolean operation on meshes based on label.
	void applyBooleanOp(MeshBooleanOp op);

	/// @brief This pipeline finish the last one step of the boolean pipeline:
	/// 4. compute explicit result and output.
	template <typename iPoint, typename iPoints, typename iTri, typename iTris>
	void
	computeFinalExplicitResult(iPoints &final_points, iTris &final_tris,
	                           std::vector<Label> *final_labels = nullptr) const;

private:
	/* Patches ****************************************************************/
	void computeSinglePatch(index_t seed, phmap::flat_hash_set<index_t> &patch);

	void computeAllPatches();

	/* Duplicate triangles ****************************************************/

	void addDuplicateTrisInfoInStructures();

	/* Ray intersection check utils *******************************************/

	static IntersInfo fast2DCheckIntersectionOnRay(const Ray    &ray,
	                                               const EPoint &tv0,
	                                               const EPoint &tv1,
	                                               const EPoint &tv2);

	static bool checkIntersectionInsideTriangle3D(const Ray    &ray,
	                                              const EPoint &tv0,
	                                              const EPoint &tv1,
	                                              const EPoint &tv2);

	static bool checkIntersectionInsideTriangle3DImplPoints(const Ray    &ray,
	                                                        const GPoint &tv0,
	                                                        const GPoint &tv1,
	                                                        const GPoint &tv2);
	/* Topology utils **********************************************************/

	static void findVertRingTris(index_t v_id, const Label &ref_label,
	                             const std::vector<index_t> &inters_tris,
	                             const std::vector<index_t> &in_tris,
	                             const std::vector<Label>   &in_labels,
	                             std::vector<index_t>       &one_ring);

	static void findEdgeTris(index_t ev0_id, index_t ev1_id,
	                         const Label                &ref_label,
	                         const std::vector<index_t> &inters_tris,
	                         const std::vector<index_t> &in_tris,
	                         const std::vector<Label>   &in_labels,
	                         std::vector<index_t>       &edge_tris);

	/* Compute inside/outside by ray casting method ****************************/

	Ray findRay(const phmap::flat_hash_set<index_t> &patch);

	std::vector<index_t> intersects_box(const Segment &rayAABB);

	Ray perturbRay(const Ray &ray, size_t i, size_t j, size_t dir);

	template <typename Pred>
	void sortIntersectedTrisAlongAxis(const Ray            &ray,
	                                  std::vector<index_t> &inters_tris);

	int perturbRayAndFindIntersTri(const Ray                  &ray,
	                               const std::vector<index_t> &tris_to_test);

	std::vector<index_t>
	pruneIntersectionsAndSortAlongRay(const Ray                  &ray,
	                                  const std::vector<index_t> &possible_inters,
	                                  const Label &patch_surface_label);

	Label analyzeSortedIntersections(const Ray                  &ray,
	                                 const std::vector<index_t> &sorted_inters);

	void
	propagateInnerLabelsOnPatch(const phmap::flat_hash_set<index_t> &patch_tris,
	                            const Label &patch_inner_label);

	void computeInOut();

	/* Apply boolean operation based on labels ******************************/

	size_t boolIntersection();
	size_t boolUnion();
	size_t boolSubtraction();
	size_t boolXOR();

public:
	MeshArrangements_Impl<Traits> arr_impl;

	/* Input data (reference to the same variable in MeshArrangement_Impl) */
	/// coordinates of all points
	std::vector<NT>      &in_coords;
	/// triangles of all meshes / triangle soups
	std::vector<index_t> &in_tris;
	/// labels of all triangles
	std::vector<size_t>  &in_labels;

	/* Arrangements data ****************************************************/
	/* Middle data */
	/// remove duplicate triangle in in_tris to get arr_in_tris
	std::vector<index_t> &arr_in_tris;
	/// convert label id to bitset to get arr_in_labels
	std::vector<Label>   &arr_in_labels;

	/* Output data */
	/// output vertices (pointers to points in arena)
	std::vector<GPoint *> &arr_out_verts;
	/// output triangles
	std::vector<index_t>  &arr_out_tris;
	/// output labels for all unique triangles
	Labels                &arr_out_labels;

	/* Auxiliary data */
	/// tree build on arr_in_tris (NOTE: not on in_tris)
	Tree                     &tree;
	/// information of removed duplicate triangles (maybe used again)
	std::vector<DuplTriInfo> &dupl_triangles;

	/* Boolean data **********************************************************/

	// triangle mesh used in boolean, build by arr_out_verts, arr_out_tris
	FTriMesh tm;

	// manifold patches in triangle mesh (tm) after arrangements
	std::vector<phmap::flat_hash_set<index_t>> patches;

	// absolute maximal coordinates of all points
	std::array<NT, 3> max_coords;

	/* Behavior control flags */
	bool verbose;
};

template <typename Traits>
void MeshBoolean_Impl<Traits>::computeLabelsPipeline()
{
	// Set config for mesh arrangements
	MeshArrangements_Config arr_config;
	arr_config.ignore_same_mesh = true;
	arr_config.verbose          = verbose;
	arr_impl.setConfig(arr_config);

	// Apply arrangments to all input meshes.
	arr_impl.meshArrangementsPipeline();

	auto start_init = OMC::Logger::elapse_reset();

	// Initialize triangle mesh used in ray-casting
	tm = FTriMesh(arr_out_verts, arr_out_tris);

	if (verbose)
		Logger::info(std::format("[OpenMeshCraft Boolean] Init Mesh time: {} s.",
		                         OMC::Logger::elapsed(start_init).count()));

	auto start_patch = OMC::Logger::elapse_reset();

	computeAllPatches();

	if (verbose)
		Logger::info(
		  std::format("[OpenMeshCraft Boolean] Compute patches time: {} s.",
		              OMC::Logger::elapsed(start_patch).count()));

	auto start_dupl = OMC::Logger::elapse_reset();

	// the informations about duplicated triangles (removed in arrangements) are
	// restored in the original structures
	addDuplicateTrisInfoInStructures();

	if (verbose)
		Logger::info(
		  std::format("[OpenMeshCraft Boolean] Add duplicate triangles time: {} s.",
		              OMC::Logger::elapsed(start_dupl).count()));

	auto start_inout = OMC::Logger::elapse_reset();

	max_coords = {tree.box().max_coord(0) + 0.5, tree.box().max_coord(1) + 0.5,
	              tree.box().max_coord(2) + 0.5};

	computeInOut();

	if (verbose)
		Logger::info(
		  std::format("[OpenMeshCraft Boolean] Compute Inside Outside time: {} s.",
		              OMC::Logger::elapsed(start_inout).count()));
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::applyBooleanOp(MeshBooleanOp op)
{
	// booleand operations
	if (op == MeshBooleanOp::INTERSECTION)
		boolIntersection();
	else if (op == MeshBooleanOp::UNION)
		boolUnion();
	else if (op == MeshBooleanOp::SUBTRACTION)
		boolSubtraction();
	else if (op == MeshBooleanOp::XOR)
		boolXOR();
	else
		OMC_THROW_DOMAIN_ERROR("Unregonized boolean operator.");
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::computeSinglePatch(
  index_t seed, phmap::flat_hash_set<index_t> &patch)
{
	OMC_UNUSED Label ref_l = arr_out_labels.surface[seed];

	// start from a seed triangle
	std::stack<index_t> tris_stack;
	tris_stack.push(seed);

	while (!tris_stack.empty())
	{
		// visit triangle on top of stack
		index_t curr_t = tris_stack.top();
		tris_stack.pop();
		// set triangle as visited, put it in patch
		tm.setTriInfo(curr_t, 1);
		patch.insert(curr_t);
		// check adjacent edges of current triangle
		for (index_t e_id : tm.adjT2E(curr_t))
		{
			if (tm.edgeIsManifold(e_id))
			{
				// if encounter a manifold edge,
				// step into the adjacent triangle.
				for (index_t t_id : tm.adjE2T(e_id))
				{
					if (t_id != curr_t && tm.triInfo(t_id) != 1)
					{
						OMC_EXPENSIVE_ASSERT(arr_out_labels.surface[t_id] == ref_l,
						                     "different label in computing single patch.");
						tris_stack.push(t_id);
					}
				}
			}
			else // e_id is not manifold -> stop flooding
			{
				// we set the vertices in the patch border with 1 (useful for ray
				// computation funcion)
				tm.setVertInfo(tm.edgeVertID(e_id, 0), 1);
				tm.setVertInfo(tm.edgeVertID(e_id, 1), 1);
			}
		}
	}
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::computeAllPatches()
{
	tm.resetVerticesInfo();
	tm.resetTrianglesInfo();

	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if (tm.triInfo(t_id) != 1)
		{
			patches.emplace_back();
			computeSinglePatch(t_id, patches.back());
		}
	}
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::addDuplicateTrisInfoInStructures()
{
	for (DuplTriInfo &item : dupl_triangles)
	{
		index_t v0_id    = arr_in_tris[3 * item.t_id];
		index_t v1_id    = arr_in_tris[3 * item.t_id + 1];
		index_t v2_id    = arr_in_tris[3 * item.t_id + 2];
		index_t new_t_id = arr_in_tris.size() / 3;

		if (item.w)
		{
			arr_in_tris.push_back(v0_id);
			arr_in_tris.push_back(v1_id);
			arr_in_tris.push_back(v2_id);
			tree.insert_triangle(arr_out_verts[v0_id], arr_out_verts[v1_id],
			                     arr_out_verts[v2_id], new_t_id);
		}
		else
		{
			arr_in_tris.push_back(v0_id);
			arr_in_tris.push_back(v2_id);
			arr_in_tris.push_back(v1_id);
			tree.insert_triangle(arr_out_verts[v0_id], arr_out_verts[v2_id],
			                     arr_out_verts[v1_id], new_t_id);
		}

		// we add the new_label to the new_triangle
		// we remove the dupl label from the orig triangle
		Label new_label;
		new_label[item.l_id] = true;
		arr_in_labels.push_back(new_label);
		arr_in_labels[item.t_id][item.l_id] = false;
	}
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::computeInOut()
{
	tbb::spin_mutex mutex;
	tbb::parallel_for(
	  (size_t)0, patches.size(),
	  [&](size_t p_id)
	  {
		  const phmap::flat_hash_set<index_t> &patch_tris = patches[p_id];

		  // label of the first triangle of the patch
		  const Label &patch_surface_label =
		    arr_out_labels.surface[*patch_tris.begin()];

		  Ray ray = findRay(patch_tris);

		  // find all the triangles having a bbox intersected by the ray
		  std::vector<index_t> possible_inters = intersects_box(ray.segment);

		  std::vector<index_t> sorted_inters = pruneIntersectionsAndSortAlongRay(
		    ray, possible_inters, patch_surface_label);

		  Label patch_inner_label = analyzeSortedIntersections(ray, sorted_inters);

		  propagateInnerLabelsOnPatch(patch_tris, patch_inner_label);
	  });
}

template <typename Traits>
bool MeshBoolean_Impl<Traits>::checkIntersectionInsideTriangle3DImplPoints(
  const Ray &ray, const GPoint &tv0, const GPoint &tv1, const GPoint &tv2)
{
	// clang-format off
	// we check the orientation of ray.segment.end() with respct to the planes
	// v0-v1-ray.segment.start(), v1-v2-ray.segment.start(), v2-v0-ray.segment.start()
	Sign or01f = Orient3D()(tv0, tv1, AsGP()(ray.segment.start()), AsGP()(ray.segment.end()));
	Sign or12f = Orient3D()(tv1, tv2, AsGP()(ray.segment.start()), AsGP()(ray.segment.end()));
	Sign or20f = Orient3D()(tv2, tv0, AsGP()(ray.segment.start()), AsGP()(ray.segment.end()));
	// clang-format on
	return (or01f != Sign::ZERO && or01f == or12f && or01f == or20f);
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::findVertRingTris(
  index_t v_id, const Label &ref_label, const std::vector<index_t> &inters_tris,
  const std::vector<index_t> &in_tris, const std::vector<Label> &in_labels,
  std::vector<index_t> &one_ring)
{
	for (index_t t_id : inters_tris)
	{
		if (in_labels[t_id] == ref_label && triContainsVert(t_id, v_id, in_tris))
			one_ring.push_back(t_id);
	}
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::findEdgeTris(
  index_t ev0_id, index_t ev1_id, const Label &ref_label,
  const std::vector<index_t> &inters_tris, const std::vector<index_t> &in_tris,
  const std::vector<Label> &in_labels, std::vector<index_t> &edge_tris)
{
	for (index_t t_id : inters_tris)
	{
		if (in_labels[t_id] == ref_label &&
		    triContainsVert(t_id, ev0_id, in_tris) &&
		    triContainsVert(t_id, ev1_id, in_tris))
			edge_tris.push_back(t_id);
	}

	OMC_EXPENSIVE_ASSERT(edge_tris.size() == 2,
	                     "problem in finding edge triangles");
	// always true in closed and manifold meshes
}

template <typename Traits>
auto MeshBoolean_Impl<Traits>::findRay(
  const phmap::flat_hash_set<index_t> &patch) -> Ray
{
	Ray     ray;
	// try to find an internal explicit point as start point of ray
	// (all operations with explicits are faster)
	index_t v_id = InvalidIndex;
	for (index_t t_id : patch)
	{
		const index_t tv[3] = {tm.triVertID(t_id, 0), tm.triVertID(t_id, 1),
		                       tm.triVertID(t_id, 2)};

		if (AsEP().is_explicit(tm.vert(tv[0])) && tm.vertInfo(tv[0]) == 0)
			v_id = tv[0];
		else if (AsEP().is_explicit(tm.vert(tv[1])) && tm.vertInfo(tv[1]) == 0)
			v_id = tv[1];
		else if (AsEP().is_explicit(tm.vert(tv[2])) && tm.vertInfo(tv[2]) == 0)
			v_id = tv[2];

		if (is_valid_idx(v_id))
		{
			const GPoint &v     = tm.vert(v_id);
			ray.segment.start() = AsEP()(v);
			ray.segment.end()   = EPoint(max_coords[0], v.y(), v.z());
			return ray;
		}
	}

	// parse triangles with all implicit points
	// clang-format off
	for (index_t t_id : patch)
	{
		EPoint a, b, c;
		a = ToEP()(tm.triVert(t_id, 0));
		b = ToEP()(tm.triVert(t_id, 1));
		c = ToEP()(tm.triVert(t_id, 2));

		if (!CollinearPoints3D().misaligned(a.data(), b.data(), c.data()))
			continue;

		int dir = MaxCompInTriNormal()(a.x(), a.y(), a.z(), b.x(), b.y(), b.z(), c.x(), c.y(), c.z());
		if (dir == 0) // dir = X
		{
			ray.segment.start() = EPoint(((a.x() + b.x() + c.x()) / 3.0) - 0.1, (a.y() + b.y() + c.y()) / 3.0, (a.z() + b.z() + c.z()) / 3.0);
			ray.segment.end() = EPoint(max_coords[0], ray.segment.start().y(), ray.segment.start().z());
			ray.dir = 'X';
		}
		else if (dir == 1) // dir = Y
		{
			ray.segment.start() = EPoint((a.x() + b.x() + c.x()) / 3.0, ((a.y() + b.y() + c.y()) / 3.0) - 0.1, (a.z() + b.z() + c.z()) / 3.0);
			ray.segment.end() = EPoint(ray.segment.start().x(), max_coords[1], ray.segment.start().z());
			ray.dir = 'Y';
		}
		else // dir = Z
		{
			ray.segment.start() = EPoint((a.x() + b.x() + c.x()) / 3.0, (a.y() + b.y() + c.y()) / 3.0, ((a.z() + b.z() + c.z()) / 3.0) - 0.1);
			ray.segment.end() = EPoint(ray.segment.start().x(), ray.segment.start().y(), max_coords[2]);
			ray.dir = 'Z';
		}

		Sign orf = Orient3D()(tm.triVert(t_id, 0), tm.triVert(t_id, 1), tm.triVert(t_id, 2), AsGP()(ray.segment.start()));
		Sign ors = Orient3D()(tm.triVert(t_id, 0), tm.triVert(t_id, 1), tm.triVert(t_id, 2), AsGP()(ray.segment.end()));

		if ((orf == Sign::NEGATIVE && ors == Sign::POSITIVE) ||
		    (orf == Sign::POSITIVE && ors == Sign::NEGATIVE))
		// the ray passes through the triangle
		{
			if (checkIntersectionInsideTriangle3DImplPoints(
			      ray, tm.triVert(t_id, 0), tm.triVert(t_id, 1), tm.triVert(t_id, 2)))
			// the ray passes inside the triangle
			{
				ray.tv[0] = tm.triVertID(t_id, 0);
				ray.tv[1] = tm.triVertID(t_id, 1);
				ray.tv[2] = tm.triVertID(t_id, 2);
				return ray;
			}
		}
	}
	// clang-format on

	std::cerr << "WARNING: the arrangement contains a fully implicit patch that "
	             "requires exact rationals for evaluation. This version of the "
	             "code does not support rationals, therefore the output result "
	             "may contain open boundaries."
	          << std::endl;
	return ray;
}

template <typename Traits>
std::vector<index_t>
MeshBoolean_Impl<Traits>::intersects_box(const Segment &rayAABB)
{
	std::vector<index_t> ids;
	tree.template all_intersections<Segment>(rayAABB, ids);
	return ids;
}

template <typename Traits>
auto MeshBoolean_Impl<Traits>::fast2DCheckIntersectionOnRay(
  const Ray &ray, const EPoint &tv0, const EPoint &tv1,
  const EPoint &tv2) -> IntersInfo
{
	NT v0[2]{}, v1[2]{}, v2[2]{}, vq[2]{};

	// clang-format off
	switch (ray.dir)
	{
	case 'X': // only YZ coordinates
	{
		v0[0] = tv0.y(); v0[1] = tv0.z();
		v1[0] = tv1.y(); v1[1] = tv1.z();
		v2[0] = tv2.y(); v2[1] = tv2.z();
		vq[0] = ray.segment.end().y();
		vq[1] = ray.segment.end().z();
	}
	break;

	case 'Y': // only xz coordinates
	{
		v0[0] = tv0.x(); v0[1] = tv0.z();
		v1[0] = tv1.x(); v1[1] = tv1.z();
		v2[0] = tv2.x(); v2[1] = tv2.z();
		vq[0] = ray.segment.end().x();
		vq[1] = ray.segment.end().z();
	}
	break;

	case 'Z': // only xy coordinates
	{
		v0[0] = tv0.x(); v0[1] = tv0.y();
		v1[0] = tv1.x(); v1[1] = tv1.y();
		v2[0] = tv2.x(); v2[1] = tv2.y();
		vq[0] = ray.segment.end().x();
		vq[1] = ray.segment.end().y();
	}
	break;
	}
	// clang-format on
	Sign or01 = Orient2D()(v0, v1, vq);
	Sign or12 = Orient2D()(v1, v2, vq);
	Sign or20 = Orient2D()(v2, v0, vq);

	if ((or01 >= Sign::ZERO && or12 >= Sign::ZERO && or20 >= Sign::ZERO) ||
	    (or01 <= Sign::ZERO && or12 <= Sign::ZERO && or20 <= Sign::ZERO))
	{
		// check if the the ray passes through a vert
		if (v0[0] == vq[0] && v0[1] == vq[1])
			return IntersInfo::INT_IN_V0;
		if (v1[0] == vq[0] && v1[1] == vq[1])
			return IntersInfo::INT_IN_V1;
		if (v2[0] == vq[0] && v2[1] == vq[1])
			return IntersInfo::INT_IN_V2;

		// check if the triangle is coplanar with the ray
		if (or01 == Sign::ZERO && or12 == Sign::ZERO)
			return IntersInfo::DISCARD;
		if (or12 == Sign::ZERO && or20 == Sign::ZERO)
			return IntersInfo::DISCARD;
		if (or20 == Sign::ZERO && or01 == Sign::ZERO)
			return IntersInfo::DISCARD;

		// check if the ray passes through an edge
		if (or01 == Sign::ZERO)
			return IntersInfo::INT_IN_EDGE01;
		if (or12 == Sign::ZERO)
			return IntersInfo::INT_IN_EDGE12;
		if (or20 == Sign::ZERO)
			return IntersInfo::INT_IN_EDGE20;

		return IntersInfo::INT_IN_TRI; // so the triangle intersect insede the
		                               // triangle area
	}

	return IntersInfo::NO_INT;
}

// offset is used to perturb the ray in all
// the possible directions
template <typename Traits>
auto MeshBoolean_Impl<Traits>::perturbRay(const Ray &ray, size_t i, size_t j,
                                          size_t dir) -> Ray
{
	Ray    new_ray = ray;
	size_t diri = (dir + 1) % 8, dirj = dir;
	NT     di = (diri > 3 ? -1. : 1.) * NT((diri % 4) < 3);
	NT     dj = (dirj > 3 ? -1. : 1.) * NT((dirj % 4) > 0);

	new_ray.segment.end()[i] =
	  std::nextafter(ray.segment.end()[i], (ray.segment.end()[i] + di));
	new_ray.segment.end()[j] =
	  std::nextafter(ray.segment.end()[j], (ray.segment.end()[j] + dj));

	return new_ray;
}

template <typename Traits>
bool MeshBoolean_Impl<Traits>::checkIntersectionInsideTriangle3D(
  const Ray &ray, const EPoint &tv0, const EPoint &tv1, const EPoint &tv2)
{
	// we check the orientation of ray.segment.end() with
	// respct to the planes v0-v1-ray.segment.start(),
	// v1-v2-ray.segment.start(), v2-v0-ray.segment.start()
	Sign or01f = Orient3D()(tv0.data(), tv1.data(), ray.segment.start().data(),
	                        ray.segment.end().data());
	Sign or12f = Orient3D()(tv1.data(), tv2.data(), ray.segment.start().data(),
	                        ray.segment.end().data());
	Sign or20f = Orient3D()(tv2.data(), tv0.data(), ray.segment.start().data(),
	                        ray.segment.end().data());

	return or01f != Sign::ZERO && or01f == or12f && or01f == or20f;
}

template <typename Traits>
template <typename Pred>
void MeshBoolean_Impl<Traits>::sortIntersectedTrisAlongAxis(
  const Ray &ray, std::vector<index_t> &inters_tris)
{
	// new LPI points
	std::vector<IPoint_LPI> arena_temp;
	arena_temp.reserve(inters_tris.size());

	// set of <t_id, impl_point>
	phmap::btree_set<std::pair<GPoint *, index_t>, Pred> inters_set;

	for (index_t t_id : inters_tris)
	{
		index_t v0_id = arr_in_tris[3 * t_id];
		index_t v1_id = arr_in_tris[3 * t_id + 1];
		index_t v2_id = arr_in_tris[3 * t_id + 2];

		std::pair<GPoint *, index_t> pair;
		pair.first  = &AsGP()(arena_temp.emplace_back(CreateLPI()(
      ray.segment.start(), ray.segment.end(), AsEP()(*arr_out_verts[v0_id]),
      AsEP()(*arr_out_verts[v1_id]), AsEP()(*arr_out_verts[v2_id]))));
		pair.second = t_id;

		inters_set.insert(pair);
	}

	inters_tris.clear();
	auto curr_int = inters_set.begin();

	// we discard the intersection before ray.first along axis
	if (ray.tv[0] != InvalidIndex) // the ray is generated
	{
		const GPoint &tv0 = *arr_out_verts[ray.tv[0]];
		const GPoint &tv1 = *arr_out_verts[ray.tv[1]];
		const GPoint &tv2 = *arr_out_verts[ray.tv[2]];

		if (Orient3D()(tv0, tv1, tv2, AsGP()(ray.segment.end())) == Sign::POSITIVE)
		{
			while (curr_int != inters_set.end() &&
			       Orient3D()(tv0, tv1, tv2, *curr_int->first) == Sign::NEGATIVE)
				curr_int++;
		}
		else
		{
			while (curr_int != inters_set.end() &&
			       Orient3D()(tv0, tv1, tv2, *curr_int->first) == Sign::POSITIVE)
				curr_int++;
		}
	}
	else // the ray is composed of 2 real explicit points
	{
		while (curr_int != inters_set.end() &&
		       Pred()(*curr_int->first, AsGP()(ray.segment.start())))
			curr_int++;
	}

	// we save all the intersecting triangles from ray.first to ray.second
	while (curr_int != inters_set.end())
	{
		inters_tris.push_back(curr_int->second);
		curr_int++;
	}
}

template <typename Traits>
int MeshBoolean_Impl<Traits>::perturbRayAndFindIntersTri(
  const Ray &ray, const std::vector<index_t> &tris_to_test)
{
	std::vector<index_t> inters_tris;
	Ray                  p_ray;

	for (int i = 0; i <= 7; i++)
	{
		if (ray.dir == 'X')
			p_ray = perturbRay(ray, 1, 2, i);
		else if (ray.dir == 'Y')
			p_ray = perturbRay(ray, 0, 2, i);
		else if (ray.dir == 'Z')
			p_ray = perturbRay(ray, 0, 1, i);

		for (index_t t_id : tris_to_test)
		{
			const EPoint tv0 = ToEP()(*arr_out_verts[arr_in_tris[3 * t_id]]);
			const EPoint tv1 = ToEP()(*arr_out_verts[arr_in_tris[3 * t_id + 1]]);
			const EPoint tv2 = ToEP()(*arr_out_verts[arr_in_tris[3 * t_id + 2]]);

			if (checkIntersectionInsideTriangle3D(p_ray, tv0, tv1, tv2))
				inters_tris.push_back(t_id);

			if (!inters_tris.empty())
				break;
		}
	}

	if (inters_tris.empty())
		return -1;

	if (ray.dir == 'X')
		sortIntersectedTrisAlongAxis<LessThanGPOnx>(p_ray, inters_tris);
	else if (ray.dir == 'Y')
		sortIntersectedTrisAlongAxis<LessThanGPOny>(p_ray, inters_tris);
	else
		sortIntersectedTrisAlongAxis<LessThanGPOnz>(p_ray, inters_tris);

	if (inters_tris.empty())
		return -1;

	// return the first triangle intersected
	return static_cast<int>(inters_tris[0]);
}

template <typename Traits>
std::vector<index_t>
MeshBoolean_Impl<Traits>::pruneIntersectionsAndSortAlongRay(
  const Ray &ray, const std::vector<index_t> &possible_inters,
  const Label &patch_surface_label)
{
	std::vector<index_t>          inters_tris;
	phmap::flat_hash_set<index_t> visited_tri;
	visited_tri.reserve(possible_inters.size() / 6);
	std::pair<phmap::flat_hash_set<index_t>::iterator, bool> ins;

	for (index_t t_id : possible_inters)
	{
		ins = visited_tri.insert(t_id);
		if (!ins.second)
			continue; // triangle already analyzed or in the one ring of a vert or
			          // in the adj of an edge

		const Label tested_tri_label = arr_in_labels[t_id];
		size_t      uint_tri_label   = LabelToIdx(tested_tri_label);
		if (patch_surface_label[uint_tri_label])
			continue; // <-- triangle of the same label of the tested patch

		const EPoint tv0 = ToEP()(*arr_out_verts[arr_in_tris[3 * t_id]]);
		const EPoint tv1 = ToEP()(*arr_out_verts[arr_in_tris[3 * t_id + 1]]);
		const EPoint tv2 = ToEP()(*arr_out_verts[arr_in_tris[3 * t_id + 2]]);

		IntersInfo ii = fast2DCheckIntersectionOnRay(ray, tv0, tv1, tv2);

		if (ii == IntersInfo::DISCARD || ii == IntersInfo::NO_INT)
			continue;

		if (ii == IntersInfo::INT_IN_TRI)
		{
			inters_tris.push_back(t_id);
		}
		else if (ii == IntersInfo::INT_IN_V0 || ii == IntersInfo::INT_IN_V1 ||
		         ii == IntersInfo::INT_IN_V2)
		{
			index_t v_id;
			if (ii == IntersInfo::INT_IN_V0)
				v_id = arr_in_tris[3 * t_id];
			else if (ii == IntersInfo::INT_IN_V1)
				v_id = arr_in_tris[3 * t_id + 1];
			else
				v_id = arr_in_tris[3 * t_id + 2];

			std::vector<index_t> vert_one_ring;
			findVertRingTris(v_id, tested_tri_label, possible_inters, arr_in_tris,
			                 arr_in_labels, vert_one_ring);
			// mark all the one ring as visited

			for (index_t t : vert_one_ring)
				visited_tri.insert(t);

			// the first inters triangle after ray perturbation
			int winner_tri = -1;
			winner_tri     = perturbRayAndFindIntersTri(ray, vert_one_ring);

			if (winner_tri != -1)
				inters_tris.push_back(winner_tri);
		}
		else if (ii == IntersInfo::INT_IN_EDGE01 ||
		         ii == IntersInfo::INT_IN_EDGE12 || ii == IntersInfo::INT_IN_EDGE20)
		{
			index_t ev0_id, ev1_id;
			if (ii == IntersInfo::INT_IN_EDGE01)
			{
				ev0_id = arr_in_tris[3 * t_id];
				ev1_id = arr_in_tris[3 * t_id + 1];
			}
			else if (ii == IntersInfo::INT_IN_EDGE12)
			{
				ev0_id = arr_in_tris[3 * t_id + 1];
				ev1_id = arr_in_tris[3 * t_id + 2];
			}
			else
			{
				ev0_id = arr_in_tris[3 * t_id + 2];
				ev1_id = arr_in_tris[3 * t_id];
			}

			std::vector<index_t> edge_tris;
			findEdgeTris(ev0_id, ev1_id, tested_tri_label, possible_inters,
			             arr_in_tris, arr_in_labels, edge_tris);

			for (index_t t : edge_tris)
				visited_tri.insert(t); // mark all the one ring as visited

			int winner_tri = -1;
			winner_tri     = perturbRayAndFindIntersTri(ray, edge_tris);

			if (winner_tri != -1)
				inters_tris.push_back(winner_tri);
		}
	}

	if (ray.dir == 'X')
		sortIntersectedTrisAlongAxis<LessThanGPOnx>(ray, inters_tris);
	else if (ray.dir == 'Y')
		sortIntersectedTrisAlongAxis<LessThanGPOny>(ray, inters_tris);
	else
		sortIntersectedTrisAlongAxis<LessThanGPOnz>(ray, inters_tris);

	return inters_tris;
}

template <typename Traits>
auto MeshBoolean_Impl<Traits>::analyzeSortedIntersections(
  const Ray &ray, const std::vector<index_t> &sorted_inters) -> Label
{
	Label patch_inner_label;
	Label visited_labels;

	for (index_t t_id : sorted_inters)
	{
		size_t t_label = LabelToIdx(arr_in_labels[t_id]);
		if (visited_labels[t_label])
			continue; // already visited patch

		const GPoint *tv0 = arr_out_verts[arr_in_tris[3 * t_id]];
		const GPoint *tv1 = arr_out_verts[arr_in_tris[3 * t_id + 1]];
		const GPoint *tv2 = arr_out_verts[arr_in_tris[3 * t_id + 2]];

		// checkOrientation -> 1 if inside, 0 if outside
		if (Orient3D()(tv0->data(), tv1->data(), tv2->data(),
		               ray.segment.end().data()) == Sign::POSITIVE)
			patch_inner_label[t_label] = true;

		visited_labels[t_label] = true;
	}

	return patch_inner_label;
}

template <typename Traits>
void MeshBoolean_Impl<Traits>::propagateInnerLabelsOnPatch(
  const phmap::flat_hash_set<index_t> &patch_tris,
  const Label                         &patch_inner_label)
{
	for (index_t t_id : patch_tris)
		arr_out_labels.inside[t_id] = patch_inner_label;
}

template <typename Traits>
size_t MeshBoolean_Impl<Traits>::boolIntersection()
{
	size_t num_tris_in_final_solution = 0;
	tm.resetTrianglesInfo();

	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if ((arr_out_labels.surface[t_id] ^ arr_out_labels.inside[t_id]).count() ==
		    arr_out_labels.num) // triangle to keep
		{
			tm.setTriInfo(t_id, 1);
			tm.setTriLabel(t_id, arr_out_labels.surface[t_id]);
			++num_tris_in_final_solution;
		}
	}

	return num_tris_in_final_solution;
}

template <typename Traits>
size_t MeshBoolean_Impl<Traits>::boolUnion()
{
	size_t num_tris_in_final_solution = 0;
	tm.resetTrianglesInfo();

	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if (arr_out_labels.inside[t_id].count() == 0) // triangle to keep
		{
			tm.setTriInfo(t_id, 1);
			tm.setTriLabel(t_id, arr_out_labels.surface[t_id]);
			++num_tris_in_final_solution;
		}
	}

	return num_tris_in_final_solution;
}

template <typename Traits>
size_t MeshBoolean_Impl<Traits>::boolSubtraction()
{
	size_t num_tris_in_final_solution = 0;
	tm.resetTrianglesInfo();

	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if (
		  // tri belongs to the first surface,
		  arr_out_labels.surface[t_id][0] &&
		  // and only belongs to the first surface.
		  arr_out_labels.surface[t_id].count() == 1 &&
		  // tri is not inside any surface.
		  arr_out_labels.inside[t_id].count() == 0)
		{
			tm.setTriInfo(t_id, 1);
			tm.setTriLabel(t_id, arr_out_labels.surface[t_id]);
			++num_tris_in_final_solution;
		}
		else if (
		  // tri does not belong to the first surface.
		  !arr_out_labels.surface[t_id][0] &&
		  // tri is inside the first surface,
		  arr_out_labels.inside[t_id][0] &&
		  // and only inside the first surface.
		  arr_out_labels.inside[t_id].count() == 1)
		{
			tm.setTriInfo(t_id, 1);
			tm.setTriLabel(t_id, arr_out_labels.surface[t_id]);
			++num_tris_in_final_solution;
		}
	}

	// fix triangles orientation
	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if (tm.triInfo(t_id) == 1 && arr_out_labels.surface[t_id][0] != 1)
			tm.flipTri(t_id);
	}

	return num_tris_in_final_solution;
}

template <typename Traits>
size_t MeshBoolean_Impl<Traits>::boolXOR()
{
	size_t num_tris_in_final_solution = 0;
	tm.resetTrianglesInfo();

	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if ((arr_out_labels.inside[t_id].count() == 0) ||
		    ((arr_out_labels.surface[t_id] ^ arr_out_labels.inside[t_id]).count() ==
		     arr_out_labels.num)) // triangle to keep
		{
			tm.setTriInfo(t_id, 1);
			tm.setTriLabel(t_id, arr_out_labels.surface[t_id]);
			++num_tris_in_final_solution;
		}
	}

	// fix triangles orientation
	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if (tm.triInfo(t_id) == 1 && arr_out_labels.inside[t_id].count() > 0)
			tm.flipTri(t_id);
	}

	return num_tris_in_final_solution;
}

template <typename Traits>
template <typename iPoint, typename iPoints, typename iTri, typename iTris>
void MeshBoolean_Impl<Traits>::computeFinalExplicitResult(
  iPoints &final_points, iTris &final_tris,
  std::vector<Label> *final_labels) const
{
	final_points.clear();
	final_tris.clear();
	final_tris.reserve(tm.numTriangles());
	if (final_labels)
	{
		final_labels->clear();
		final_labels->reserve(tm.numTriangles());
	}

	// loop over triangles and fix vertex indices
	OMC_UNUSED size_t    num_faces    = 0;
	size_t               num_vertices = 0;
	std::vector<index_t> vertex_index(tm.numVerts(), InvalidIndex);

	for (index_t t_id = 0; t_id < tm.numTriangles(); ++t_id)
	{
		if (tm.triInfo(t_id) == 0)
			continue; // triangle not included in final version

		const index_t         *triangle = tm.tri(t_id);
		std::array<index_t, 3> out_tris;
		for (size_t i = 0; i < 3; i++)
		{
			index_t old_vid = triangle[i];
			if (!is_valid_idx(vertex_index[old_vid]))
			{
				vertex_index[old_vid] = num_vertices++;
			}
			out_tris[i] = vertex_index[old_vid];
		}

		final_tris.push_back(iTri(out_tris[0], out_tris[1], out_tris[2]));
		if (final_labels)
			final_labels->push_back(tm.triLabel(t_id));
		num_faces++;
	}

	final_points.resize(num_vertices);

	// loop over vertices, fix index
	tbb::parallel_for(index_t(0), arr_out_verts.size(),
	                  [this, &vertex_index, &final_points](index_t v_id)
	                  {
		                  if (!is_valid_idx(vertex_index[v_id]))
			                  return;
		                  const GPoint *gp = arr_out_verts[v_id];
		                  EPoint        ep = ToEP()(*gp);
		                  final_points[vertex_index[v_id]] =
		                    iPoint(ep.x(), ep.y(), ep.z());
	                  });
}

/*****************************************************************************/
/*************** Implementations of interface class **************************/
/*****************************************************************************/

template <typename Kernel, typename Traits>
MeshBoolean<Kernel, Traits>::MeshBoolean(bool _verbose)
  : verbose(_verbose)
{
}

template <typename Kernel, typename Traits>
MeshBoolean<Kernel, Traits>::~MeshBoolean()
{
}

template <typename Kernel, typename Traits>
index_t
MeshBoolean<Kernel, Traits>::addTriMeshAsInput(const iPoints    &points,
                                               const iTriangles &triangles)
{
	input_meshes.emplace_back();
	input_meshes.back().points    = &points;
	input_meshes.back().triangles = &triangles;
	if (input_meshes.size() == NBIT)
		OMC_THROW_OUT_OF_RANGE("Input meshes for arrangments are too much, limit "
		                       "the number to less than {}",
		                       NBIT);
	// once the input meshes change, reset internal data
	m_impl = nullptr;
	return input_meshes.size() - 1;
}

template <typename Kernel, typename Traits>
void MeshBoolean<Kernel, Traits>::setTriMeshAsOutput(iPoints    &points,
                                                     iTriangles &triangles)
{
	output_points    = &points;
	output_triangles = &triangles;
}

template <typename Kernel, typename Traits>
void MeshBoolean<Kernel, Traits>::setOutputLabels(std::vector<Label> &labels)
{
	output_labels = &labels;
}

template <typename Kernel, typename Traits>
void MeshBoolean<Kernel, Traits>::clear()
{
	input_meshes.clear();
	output_points    = nullptr;
	output_triangles = nullptr;
	output_labels    = nullptr;
	m_impl           = nullptr;
}

template <typename Kernel, typename Traits>
void MeshBoolean<Kernel, Traits>::computeLabels()
{
	OMC_THROW_DOMAIN_ERROR_IF(input_meshes.size() <= 1,
	                          "need two meshes at less!");
	OMC_THROW_DOMAIN_ERROR_IF(output_points == nullptr ||
	                            output_triangles == nullptr,
	                          "output mesh is not set.");
	m_impl = std::make_unique<MeshBoolean_Impl<BooleanTraits>>(verbose);
	if (!loadMultipleMeshes()(input_meshes, m_impl->in_coords, m_impl->in_tris,
	                          m_impl->in_labels))
	{
		OMC_THROW_DOMAIN_ERROR("empty input meshes.");
		return;
	}
	m_impl->computeLabelsPipeline();
}

#define BOOL_OPERATION(var, bool_op)           \
	template <typename Kernel, typename Traits>  \
	void MeshBoolean<Kernel, Traits>::var()      \
	{                                            \
		MeshBooleanOp op = MeshBooleanOp::bool_op; \
		booleanPipeline(op);                       \
	}

BOOL_OPERATION(Union, UNION);
BOOL_OPERATION(Intersection, INTERSECTION);
BOOL_OPERATION(Subtraction, SUBTRACTION);
BOOL_OPERATION(Xor, XOR);

#undef BOOL_OPERATION

template <typename Kernel, typename Traits>
void MeshBoolean<Kernel, Traits>::booleanPipeline(MeshBooleanOp op)
{
	OMC_THROW_DOMAIN_ERROR_IF(input_meshes.size() <= 1,
	                          "need two meshes at less!");
	OMC_THROW_DOMAIN_ERROR_IF(output_points == nullptr ||
	                            output_triangles == nullptr,
	                          "output mesh is not set.");
	OMC_THROW_DOMAIN_ERROR_IF(m_impl == nullptr,
	                          "arrangements and label results are not computed.");

	m_impl->applyBooleanOp(op);
	m_impl->template computeFinalExplicitResult<iPoint, iPoints, iTriangle,
	                                            iTriangles>(
	  *output_points, *output_triangles, output_labels);
}

template <typename Kernel, typename Traits>
class MeshBoolean<Kernel, Traits>::loadMultipleMeshes
{
public:
	bool operator()(const std::vector<iTriSoup>      &input_meshes,
	                std::vector<typename Kernel::NT> &coords,
	                std::vector<index_t>             &tris,
	                std::vector<size_t>              &arr_out_labels)
	{
		for (size_t mesh_id = 0; mesh_id < input_meshes.size(); mesh_id++)
		{
			load(*input_meshes[mesh_id].points, *input_meshes[mesh_id].triangles,
			     mesh_id, coords, tris, arr_out_labels);
		}
		return !coords.empty() && !tris.empty();
	}
};

} // namespace OMC