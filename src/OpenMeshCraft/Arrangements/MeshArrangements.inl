#pragma once

#include "MeshArrangements.h"

// Data structures of arrangements
#include "FastTriMesh.h"
#include "TriangleSoup.h"

// Sub-algorithms of arrangements
#include "DetectClassifyTTIs.h"
#include "Triangulation.h"

// Geometry kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"

namespace OMC {

template <typename Kernel, typename Traits>
class MeshArrangements<Kernel, Traits>::ArrangementsTraits
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
	using CollinearPoints3D  = typename K::CollinearPoints3D;
	using LongestAxis        = typename K::LongestAxis;
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

/// @brief Implement class of MeshArrangments
template <typename Traits>
class MeshArrangements_Impl
{
public: /* Traits ************************************************************/
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
	using CollinearPoints3D  = typename Traits::CollinearPoints3D;
	using MaxCompInTriNormal = typename Traits::MaxCompInTriNormal;

public: /* Auxiliary data structures *****************************************/
	// tree
	using Tree     = Arr_Tree_Intersection<Traits>;
	// point arena
	using PntArena = PointArena<Traits>;
	// triangle soup
	using TriSoup  = TriangleSoup<Traits>;
	// fast triangle mesh
	using FastMesh = FastTriMesh<Traits>;

public: /* Constructors ******************************************************/
	MeshArrangements_Impl(MeshArrangements_Stats *_stats = nullptr)
	  : stats(_stats)
	{
	}

	/// @brief An experimental interface to set varying parameters.
	void setConfig(MeshArrangements_Config _config) { config = _config; }

public: /* Pipeline **********************************************************/
	/// @brief Detects intersections, classifies them, and triangulates
	/// all triangles with intersection points and constrained segments.
	/// @note All useful intermediate data will be stored until it is destroyed.
	/// @note The explicit result will not be output.
	void meshArrangementsPipeline();

	/// @brief Output explicit result (points and triangles).
	template <typename iPoint, typename iPoints, typename iTri, typename iTris>
	void computeExplicitResult(iPoints &final_points, iTris &final_tris,
	                           std::vector<Label> *final_labels = nullptr);

public: /* Preprocessing steps ***********************************************/
	void convertLabels();

	void mergeDuplicatedVertices();

	void removeDegenerateAndDuplicatedTriangles();

public: /* Routines before and after solving intersecions ********************/
	void initBeforeDetectClassify();

	void exitAfterTriangulation();

public:
	/* Input data */
	/// coordinates of all points
	std::vector<NT>      in_coords;
	/// triangles of all meshes / triangle soups
	std::vector<index_t> in_tris;
	/// labels of all triangles
	std::vector<size_t>  in_labels;

	/* Middle data */
	/// remove duplicate triangle in in_tris to get arr_in_tris
	std::vector<index_t> arr_in_tris;
	/// convert label id to bitset to get arr_in_labels
	std::vector<Label>   arr_in_labels;

	/* Output data */
	/// output vertices (pointers to points in arena)
	std::vector<GPoint *> arr_out_verts;
	/// output triangles
	std::vector<index_t>  arr_out_tris;
	/// output labels for all unique triangles
	Labels                arr_out_labels;

	/* Auxiliary data */
	/// tree build on arr_in_tris (NOTE: not on in_tris)
	Tree                     tree;
	/// information of removed duplicate triangles (maybe used again)
	std::vector<DuplTriInfo> dupl_triangles;

	/* Configuration */
	MeshArrangements_Config config;

	/* Statistics */
	MeshArrangements_Stats *stats;

private: /* Private middle data *******************************************/
	/// point arena for explicit points
	PntArena              exp_pnt_arena;
	/// All generated points in algorithm are stored in pnt_arena
	std::vector<PntArena> pnt_arenas;
	/// Triangle soup
	TriSoup               tri_soup;
};

template <typename Traits>
void MeshArrangements_Impl<Traits>::meshArrangementsPipeline()
{
	OMC_ASSERT(!in_coords.empty() && !in_tris.empty(), "empty input.");
	OMC_ASSERT(in_tris.size() % 3 == 0, "triangle size error.");
	OMC_ASSERT(in_tris.size() / 3 == in_labels.size(),
	           "size of triangles and labels mismatche.");

	/***** Preprocessing *****/

	OMC_ARR_START_ELAPSE(start_pp);

	// convert in_labels to arr_in_labels.
	// (in_labels will NOT be used again later).
	convertLabels();

	OMC_ARR_START_ELAPSE(start_mdv);

	// 1. merge duplicate vertices in in_coords and store them as explicit points
	// in arena.
	// 2. store pointers to explicit points in arr_out_verts.
	// 3. fix indices in in_tris to get arr_in_tris.
	// (in_coords and in_tris will NOT be used again later).
	mergeDuplicatedVertices();

	OMC_ARR_SAVE_ELAPSED(start_mdv, mdv_elapsed, "Merge duplicate vertices");
	OMC_ARR_START_ELAPSE(start_rdt);

	// 1. remove degenerate and duplicate triangles in arr_in_tris.
	// 2. merge labels of duplicate triangles to the unique one in arr_in_labels.
	// (so, multiple bits in arr_in_labels(std::bitset) are possibly set to true)
	removeDegenerateAndDuplicatedTriangles();

	// After removing degenerate triangles, there are possibly isolated vertices.
	// Remove them before arrangements pipeline if it is neccessary one day.
	// But even if isolated vertices are not removed here, they will be removed
	// when outputing explicit results. Because the unreferenced vertices (i.e.,
	// isolated vertices) won't be output.

	OMC_ARR_SAVE_ELAPSED(start_rdt, rdt_elapsed,
	                     "Remove degenerate and duplicate triangles");
	OMC_ARR_START_ELAPSE(start_tree);

	// initialize tree from triangle soup (vertices and triangles)
	tree.init_from_triangle_soup(arr_out_verts, arr_in_tris,
	                             arr_in_tris.size() / 3 + dupl_triangles.size(),
	                             config);

	OMC_ARR_SAVE_ELAPSED(start_tree, tree_elapsed, "Build tree");
	OMC_ARR_SAVE_ELAPSED(start_pp, pp_elapsed, "Preprocessing");

	/***** Arrangements *****/

	OMC_ARR_START_ELAPSE(start_ci);

	// build triangle soup
	initBeforeDetectClassify();

	// Classify intersections.
	DetectClassifyTTIs<Traits> DCI(tri_soup, tree, config, *stats);

	OMC_ARR_SAVE_ELAPSED(start_ci, ci_elapsed, "Classify intersection");
	OMC_ARR_START_ELAPSE(start_tr);

	// Triangulation.
	Triangulation<Traits> TR(tri_soup, arr_out_tris, arr_out_labels.surface,
	                         config, *stats);

	// collect vertices, triangles and labels.
	exitAfterTriangulation();

	OMC_ARR_SAVE_ELAPSED(start_tr, tr_elapsed, "Triangulation");
}

template <typename Traits>
void MeshArrangements_Impl<Traits>::convertLabels()
{
	arr_in_labels.resize(in_labels.size());
	Label mask;

	for (size_t i = 0; i < in_labels.size(); i++)
	{
		arr_in_labels[i][in_labels[i]] = true;
		mask[in_labels[i]]             = true;
	}
	arr_out_labels.num = mask.count();
}

template <typename Traits>
void MeshArrangements_Impl<Traits>::mergeDuplicatedVertices()
{
	bool parallel = in_coords.size() / 3 > 10000;
	arr_out_verts.reserve(in_coords.size() / 3);
	exp_pnt_arena.init.reserve(in_coords.size() / 3);
	arr_in_tris.reserve(in_tris.size());

	using vec3    = std::array<NT, 3>;
	vec3 *in_vecs = (vec3 *)in_coords.data();

	std::vector<index_t> sorted(in_coords.size() / 3);
	std::iota(sorted.begin(), sorted.end(), 0);
	size_t origin_num = sorted.size();

	if (parallel)
		tbb::parallel_sort(sorted.begin(), sorted.end(), [in_vecs](auto a, auto b)
		                   { return in_vecs[a] < in_vecs[b]; });
	else
		std::sort(sorted.begin(), sorted.end(),
		          [in_vecs](auto a, auto b) { return in_vecs[a] < in_vecs[b]; });

	std::vector<index_t> lookup(origin_num);
	for (size_t idx = 0; idx < origin_num; idx++)
	{
		if (idx == 0 || in_vecs[sorted[idx]] != in_vecs[sorted[idx - 1]])
		{
			const vec3 &v = in_vecs[sorted[idx]];
			// create explicit point and save it as generic point.
			arr_out_verts.push_back(
			  &AsGP()(exp_pnt_arena.init.emplace_back(v[0], v[1], v[2])));
		}
		lookup[sorted[idx]] = arr_out_verts.size() - 1;
	}

	arr_in_tris.resize(in_tris.size());
	if (parallel)
		std::transform(std::execution::par_unseq, in_tris.begin(), in_tris.end(),
		               arr_in_tris.begin(),
		               [&lookup](index_t idx) { return lookup[idx]; });
	else
		std::transform(std::execution::seq, in_tris.begin(), in_tris.end(),
		               arr_in_tris.begin(),
		               [&lookup](index_t idx) { return lookup[idx]; });

	if (config.verbose)
	{
		Logger::info(
		  std::format("[OpenMeshCraft Arrangements] removed {} duplicate vertices.",
		              std::to_string(in_coords.size() / 3 - arr_out_verts.size())));
	}
}

template <typename Traits>
void MeshArrangements_Impl<Traits>::removeDegenerateAndDuplicatedTriangles()
{
	using vec3i           = std::array<index_t, 3>;
	size_t num_orig_verts = arr_out_verts.size();
	size_t num_orig_tris  = arr_in_tris.size() / 3;
	vec3i *data_orig_tris = (vec3i *)arr_in_tris.data();
	size_t tri_off        = 0;

	// compute collinear
	std::vector<uint8_t> collinear_res(num_orig_tris, false);
	tbb::parallel_for((size_t)0, num_orig_tris,
	                  [this, data_orig_tris, &collinear_res](index_t t_id)
	                  {
		                  vec3i &t = data_orig_tris[t_id];
		                  collinear_res[t_id] =
		                    CollinearPoints3D()(arr_out_verts[t[0]]->data(),
		                                        arr_out_verts[t[1]]->data(),
		                                        arr_out_verts[t[2]]->data());
	                  });

	bool parallel = num_orig_tris > 10000;
	if (parallel)
	{
		// map: edge (two larger tri vertices) -> triangle_id
		using TriMap = phmap::flat_hash_map<std::pair<index_t, index_t>, index_t>;
		// the smallest tri vertex -> map
		std::vector<tbb::spin_mutex> tris_map_mutex(num_orig_verts);
		std::vector<TriMap>          tris_map(num_orig_verts);

		// unique index or compact index
		std::vector<index_t> tris_idx(num_orig_tris);

		bool exist_removed_tri = false;

		auto check_tri_is_unique = [this, &collinear_res, &tris_map,
		                            &tris_map_mutex, &tris_idx,
		                            &exist_removed_tri](index_t t_id)
		{
			if (collinear_res[t_id])
			{
				tris_idx[t_id]    = InvalidIndex; // means removed
				exist_removed_tri = true;
				return;
			}

			index_t v0_id = arr_in_tris[(3 * t_id)];
			index_t v1_id = arr_in_tris[(3 * t_id) + 1];
			index_t v2_id = arr_in_tris[(3 * t_id) + 2];

			std::array<index_t, 3> tri = {v0_id, v1_id, v2_id};
			std::sort(tri.begin(), tri.end());

			std::pair<index_t, index_t> edge = {tri[1], tri[2]};

			{ // critical section
				std::lock_guard<tbb::spin_mutex> lock(tris_map_mutex[tri[0]]);

				auto ins = tris_map[tri[0]].insert({edge, t_id});
				// if we meet this triangle for the first time, it is unique and
				// ins.first->second is same as t_id. otherwise, it is duplicate and
				// ins.first->second is the index of existed triangle.
				if (t_id < ins.first->second)
				{ // always save unique triangle with lower index
					tris_idx[ins.first->second] = t_id;
					tris_idx[t_id]              = t_id;
					tris_map[tri[0]][edge]      = t_id;
				}
				else
					tris_idx[t_id] = ins.first->second;
				if (!ins.second)
					exist_removed_tri = true;
			}
		};

		// parallel build map from triangle to unique tri
		tbb::parallel_for((size_t)0, num_orig_tris, check_tri_is_unique);

		if (!exist_removed_tri)
		{
			/* #region log */
			if (config.verbose)
			{
				Logger::info(std::format(
				  "[OpenMeshCraft Arrangements] No degenerate and duplicate "
				  "triangle found."));
			}
			/* #endregion */
			return;
		}
		// clear to save mem
		tris_map       = std::vector<TriMap>();
		tris_map_mutex = std::vector<tbb::spin_mutex>();

		// loop as before by use simpler way
		std::vector<index_t> compact_tris_idx(num_orig_tris);

		// the first traversal to map unique triangle to compact index
		for (index_t t_id = 0; t_id < num_orig_tris; t_id++)
		{
			if (!is_valid_idx(tris_idx[t_id]))
				continue;

			index_t v0_id = arr_in_tris[3 * t_id];
			index_t v1_id = arr_in_tris[3 * t_id + 1];
			index_t v2_id = arr_in_tris[3 * t_id + 2];
			Label   label = arr_in_labels[t_id];

			if (tris_idx[t_id] == t_id)
			{
				// triangle is unique, save it and its compact id.
				arr_in_tris[tri_off * 3]     = v0_id;
				arr_in_tris[tri_off * 3 + 1] = v1_id;
				arr_in_tris[tri_off * 3 + 2] = v2_id;
				arr_in_labels[tri_off]       = label;

				// update unique index to compact index
				compact_tris_idx[t_id] = tri_off;
				tri_off += 1;
			}
			else
			{
				// triangle is duplicate, save info abount duplication
				index_t unique_idx = t_id;
				while (unique_idx != tris_idx[unique_idx])
					unique_idx = tris_idx[unique_idx];
				index_t compact_idx = compact_tris_idx[unique_idx];

				arr_in_labels[compact_idx] |= label;
				size_t mesh_label = LabelToIdx(label);

				index_t curr_tri_verts[] = {v0_id, v1_id, v2_id};
				index_t uniq_tri_verts[] = {arr_in_tris[compact_idx * 3],
				                            arr_in_tris[compact_idx * 3 + 1],
				                            arr_in_tris[compact_idx * 3 + 2]};
				bool    w = consistentWinding(curr_tri_verts, uniq_tri_verts);
				dupl_triangles.push_back({/*compact triangle id*/ compact_idx,
				                          /*label of the actual triangle*/ mesh_label,
				                          /*winding with respect to the triangle stored
				                             in mesh (true -> same, false -> opposite)*/
				                          w});
			}
		}
	}
	else
	{
		// loop as before by use simpler way
		// map: tri_vertices -> tri_off
		phmap::flat_hash_map<std::array<index_t, 3>, size_t> tris_map;

		for (size_t t_id = 0; t_id < num_orig_tris; ++t_id)
		{
			if (collinear_res[t_id])
				continue;
			index_t v0_id = arr_in_tris[(3 * t_id)];
			index_t v1_id = arr_in_tris[(3 * t_id) + 1];
			index_t v2_id = arr_in_tris[(3 * t_id) + 2];
			Label   label = arr_in_labels[t_id];

			std::array<index_t, 3> tri = {v0_id, v1_id, v2_id};
			std::sort(tri.begin(), tri.end());

			auto ins = tris_map.insert({tri, tri_off});

			if (ins.second) // first time for tri v0, v1, v2
			{
				arr_in_tris[tri_off * 3]     = v0_id;
				arr_in_tris[tri_off * 3 + 1] = v1_id;
				arr_in_tris[tri_off * 3 + 2] = v2_id;
				arr_in_labels[tri_off]       = label;
				tri_off += 1;
			}
			else // triangle already present -> save info about duplicates
			{
				size_t orig_tri_off = ins.first->second;
				arr_in_labels[orig_tri_off] |= label; // label for duplicates

				size_t mesh_label = LabelToIdx(label);
				OMC_EXPENSIVE_ASSERT(mesh_label >= 0, "invalid label id");

				index_t curr_tri_verts[] = {v0_id, v1_id, v2_id};
				index_t orig_tri_verts[] = {arr_in_tris[orig_tri_off * 3],
				                            arr_in_tris[orig_tri_off * 3 + 1],
				                            arr_in_tris[orig_tri_off * 3 + 2]};

				bool w = consistentWinding(curr_tri_verts, orig_tri_verts);

				dupl_triangles.push_back({/*original triangle id*/ orig_tri_off,
				                          /*label of the actual triangle*/ mesh_label,
				                          /*winding with respect to the triangle stored
				                             in mesh (true -> same, false -> opposite)*/
				                          w});
			}
		}
	}

	/* #region log */
	if (config.verbose)
	{
		Logger::info(
		  std::format("[OpenMeshCraft Arrangements] removed {}  degenerate and "
		              "duplicate triangles.",
		              std::to_string(arr_in_labels.size() - tri_off)));
	}
	/* #endregion */

	arr_in_tris.resize(tri_off * 3);
	arr_in_labels.resize(tri_off);
}

template <typename Traits>
void MeshArrangements_Impl<Traits>::initBeforeDetectClassify()
{
	pnt_arenas = std::vector<PntArena>(tbb::this_task_arena::max_concurrency());

	for (GPoint *v : arr_out_verts)
		tri_soup.addImplVert(v);
	tri_soup.triangles  = arr_in_tris;   // copy, do not move
	tri_soup.tri_labels = arr_in_labels; // copy, do not move
	tri_soup.pnt_arenas = &pnt_arenas;
	tri_soup.initialize();
}

template <typename Traits>
void MeshArrangements_Impl<Traits>::exitAfterTriangulation()
{
	// initialize merged triangle soup and auxiliary structure
	arr_out_verts.resize(tri_soup.numVerts());
	std::copy(std::execution::par_unseq, tri_soup.vertices.begin(),
	          tri_soup.vertices.end(), arr_out_verts.begin());

	arr_out_labels.inside.resize(arr_out_tris.size() / 3);

#if defined(OMC_ARR_PROFILE)
	std::vector<uint8_t> vertex_used(arr_out_verts.size(), false);
	for (index_t vi : arr_out_tris)
		vertex_used[vi] = true;

	for (index_t vi = 0; vi < arr_out_verts.size(); vi++)
	{
		if (vertex_used[vi] && arr_out_verts[vi]->is_Implicit())
		{
			// count the number of different types of implicit points
			OMC_ARR_PROFILE_INC_TOTAL(ArrFuncNames::IP_CNT);
			OMC_ARR_PROFILE_INC_REACH(
			  ArrFuncNames::IP_CNT,
			  static_cast<size_t>(arr_out_verts[vi]->point_type()));

	#ifdef OMC_ARR_PROF_MAXVAR
			// save the order differences of maxvar
			NT indirect_maxvar, offset_maxvar;
			indirect_maxvar = arr_out_verts[vi]->getIndirectMaxVar();
			offset_maxvar   = arr_out_verts[vi]->getOffsetMaxVar();

			NT  order         = indirect_maxvar / offset_maxvar;
			int rounded_order = (int)std::floor(std::log10(order) * 10.);
			int rounded_offset_order =
			  std::min(std::max(rounded_order + 100, 0), 200);

			OMC_ARR_PROFILE_INC_TOTAL(ArrFuncNames::IP_MAXVAR_ORDER);
			OMC_ARR_PROFILE_INC_REACH(ArrFuncNames::IP_MAXVAR_ORDER,
			                          rounded_offset_order);
	#endif

	#ifdef OMC_ARR_PROF_TPI_LENGTH
			// get the length of TPI's expansion
			NT *lx = nullptr, *ly = nullptr, *lz = nullptr, *ld = nullptr;
			int lx_len = 0, ly_len = 0, lz_len = 0, d_len = 0;
			NT  bx, by, bz;
			arr_out_verts[vi]->getExpansionLambda(&lx, lx_len, &ly, ly_len, &lz,
			                                      lz_len, &ld, d_len, bx, by, bz);

			size_t point_type_name = static_cast<size_t>(ArrFuncNames::SSI_EXP_LEN);
			if (arr_out_verts[vi]->is_LPI())
				point_type_name = static_cast<size_t>(ArrFuncNames::LPI_EXP_LEN);
			if (arr_out_verts[vi]->is_TPI())
				point_type_name = static_cast<size_t>(ArrFuncNames::TPI_EXP_LEN);

			// clang-format off
			OMC_ARR_PROFILE_INC_TOTAL_CNT(static_cast<ArrFuncNames>(point_type_name), 4);
			OMC_ARR_PROFILE_INC_REACH(static_cast<ArrFuncNames>(point_type_name), lx_len);
			OMC_ARR_PROFILE_INC_REACH(static_cast<ArrFuncNames>(point_type_name), ly_len);
			OMC_ARR_PROFILE_INC_REACH(static_cast<ArrFuncNames>(point_type_name), lz_len);
			OMC_ARR_PROFILE_INC_REACH(static_cast<ArrFuncNames>(point_type_name), d_len);
			// clang-format on

			if (lx)
				FreeDoubles(lx);
			if (ly)
				FreeDoubles(ly);
			if (lz)
				FreeDoubles(lz);
			if (ld)
				FreeDoubles(ld);
	#endif
		}
	}

	#ifdef OMC_ARR_PROF_ES_PNTS
	std::fstream fout;
	fout.open("./data/test_output/arrangements/pnts_on_edges.txt", std::ios::out);
	if (fout.is_open())
	{
		size_t sum = 0;
		for (index_t ei = 0; ei < tri_soup.numEdges(); ei++)
		{
			auto &e2p = tri_soup.edgePointsList(ei);
			if (!e2p.empty())
			{
				sum += e2p.size();
				fout << e2p.size() << std::endl;
			}
		}
		if (sum != 0)
			fout << sum << std::endl;
	}
	fout.close();

	fout.open("./data/test_output/arrangements/pnts_on_segs.txt", std::ios::out);
	if (fout.is_open())
	{
		size_t sum = 0;
		for (const auto &s2p : tri_soup.seg2pts)
		{
			if (!s2p.empty())
			{
				sum += s2p.size();
				fout << s2p.size() << std::endl;
			}
		}
		if (sum != 0)
			fout << sum << std::endl;
	}
	fout.close();
	#endif
#endif
}

template <typename Traits>
template <typename iPoint, typename iPoints, typename iTri, typename iTris>
void MeshArrangements_Impl<Traits>::computeExplicitResult(
  iPoints &final_points, iTris &final_tris, std::vector<Label> *final_labels)
{
	final_points.clear();
	final_tris.clear();
	final_tris.resize(arr_out_tris.size() / 3);
	if (final_labels)
	{
		final_labels->clear();
		final_labels->resize(arr_out_tris.size() / 3);
	}

	// loop over triangles and fix vertex indices
	size_t               num_vertices = 0;
	std::vector<index_t> vertex_index(arr_out_verts.size(), InvalidIndex);
	for (index_t t_id = 0; t_id < arr_out_tris.size(); t_id += 3)
	{
		const index_t         *triangle = &arr_out_tris[t_id];
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

		final_tris[t_id / 3] = iTri(out_tris[0], out_tris[1], out_tris[2]);
		if (final_labels)
			(*final_labels)[t_id / 3] = arr_out_labels.surface[t_id / 3];
	}

	// loop over vertices
	final_points.resize(num_vertices);
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
MeshArrangements<Kernel, Traits>::MeshArrangements()
{
}

template <typename Kernel, typename Traits>
MeshArrangements<Kernel, Traits>::~MeshArrangements()
{
}

template <typename Kernel, typename Traits>
index_t
MeshArrangements<Kernel, Traits>::addTriMeshAsInput(const iPoints    &points,
                                                    const iTriangles &triangles)
{
	input_meshes.emplace_back();
	input_meshes.back().points    = &points;
	input_meshes.back().triangles = &triangles;
	if (input_meshes.size() == NBIT)
		OMC_THROW_OUT_OF_RANGE("Input meshes for arrangments are too much, limit "
		                       "the number to less than {}.",
		                       NBIT);
	return input_meshes.size() - 1;
}

template <typename Kernel, typename Traits>
void MeshArrangements<Kernel, Traits>::setTriMeshAsOutput(iPoints    &points,
                                                          iTriangles &triangles)
{
	output_points    = &points;
	output_triangles = &triangles;
}

template <typename Kernel, typename Traits>
void MeshArrangements<Kernel, Traits>::setOutputLabels(
  std::vector<Label> &labels)
{
	output_labels = &labels;
}

template <typename Kernel, typename Traits>
void MeshArrangements<Kernel, Traits>::clear()
{
	input_meshes.clear();
	output_points    = nullptr;
	output_triangles = nullptr;
	output_labels    = nullptr;
	m_impl           = nullptr;
}

template <typename Kernel, typename Traits>
void MeshArrangements<Kernel, Traits>::meshArrangements()
{
	m_impl =
	  std::make_unique<MeshArrangements_Impl<ArrangementsTraits>>(&arr_stats);

	if (!loadMultipleMeshes()(input_meshes, m_impl->in_coords, m_impl->in_tris,
	                          m_impl->in_labels))
	{
		OMC_THROW_DOMAIN_ERROR("Empty input meshes.");
		return;
	}

	m_impl->setConfig(config);
	m_impl->meshArrangementsPipeline();

	if (config.output_explicit_result)
	{
		OMC_THROW_DOMAIN_ERROR_IF(output_points == nullptr ||
		                            output_triangles == nullptr,
		                          "Output mesh is not set.");
		m_impl
		  ->template computeExplicitResult<iPoint, iPoints, iTriangle, iTriangles>(
		    *output_points, *output_triangles, output_labels);
		m_impl = nullptr;
	}
}

template <typename Kernel, typename Traits>
void MeshArrangements<Kernel, Traits>::setConfig(
  MeshArrangements_Config _config)
{
	config = _config;
}

template <typename Kernel, typename Traits>
MeshArrangements_Stats &MeshArrangements<Kernel, Traits>::stats()
{
	return arr_stats;
}

template <typename Kernel, typename Traits>
class MeshArrangements<Kernel, Traits>::loadMultipleMeshes
{
public:
	bool operator()(const std::vector<iTriSoup>      &input_meshes,
	                std::vector<typename Kernel::NT> &coords,
	                std::vector<index_t> &tris, std::vector<size_t> &labels)
	{
		for (size_t mesh_id = 0; mesh_id < input_meshes.size(); mesh_id++)
		{
			load(*input_meshes[mesh_id].points, *input_meshes[mesh_id].triangles,
			     mesh_id, coords, tris, labels);
		}
		return !coords.empty() && !tris.empty();
	}
};

} // namespace OMC