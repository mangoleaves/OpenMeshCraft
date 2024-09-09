#pragma once

#include "MeshArrangements.h"

// Data structures of arrangements
#include "FastTriMesh.h"
#include "TriangleSoup.h"

// Sub-algorithms of arrangements
#include "CleanMesh.h"
#include "DetectClassifyTTIs.h"
#include "Triangulation.h"

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
	using PntArena = ArrPointArena<Traits>;
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

public: /* Routines before and after solving intersecions ********************/
	void collectCleanResults(ArrCleanMesh<Traits> &CM);

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
	ArrLabels             arr_out_labels;

	/* Auxiliary data */
	/// tree build on arr_in_tris (NOTE: not on in_tris)
	Tree                        tree;
	/// information of removed duplicate triangles (maybe used again)
	std::vector<ArrDuplTriInfo> dupl_triangles;

	/* Configuration */
	MeshArrangements_Config config;

	/* Statistics */
	MeshArrangements_Stats *stats;

private: /* Private middle data *******************************************/
	/// explicit points and jolly points
	std::vector<EPoint> exp_pnt; // explicit points
#ifdef OMC_ARR_AUX_LPI
	std::vector<EPoint> jolly_pnt; // jolly points
#endif
	/// All generated implicit points in algorithm are stored in pnt_arena
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

	// clean input mesh
	ArrCleanMesh<Traits> CM(in_coords, in_tris, in_labels);
	CM.convertLabels();
	CM.mergeDuplicatedVertices();
	CM.removeDegenerateAndDuplicatedTriangles();
	CM.removeIsolatedVertices();
	collectCleanResults(CM);

	OMC_ARR_SAVE_ELAPSED(start_pp, pp_elapsed, "Preprocessing");

	/***** Build tree *****/

	OMC_ARR_START_ELAPSE(start_tree);

	// initialize tree from triangle soup (vertices and triangles)
	tree.init_from_triangle_soup(arr_out_verts, arr_in_tris,
	                             arr_in_tris.size() / 3 + dupl_triangles.size(),
	                             config);

	OMC_ARR_SAVE_ELAPSED(start_tree, tree_elapsed, "Build tree");

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
void MeshArrangements_Impl<Traits>::collectCleanResults(
  ArrCleanMesh<Traits> &CM)
{
	// collect vertices
	arr_out_verts.reserve(CM.out_coords.size() / 3);
	exp_pnt.reserve(CM.out_coords.size() / 3);
	for (index_t vi = 0; vi < CM.out_coords.size(); vi += 3)
		arr_out_verts.push_back(&AsGP()(exp_pnt.emplace_back(&CM.out_coords[vi])));
	// collect triangles
	arr_in_tris        = std::move(CM.out_tris);
	// collect labels
	arr_in_labels      = std::move(CM.out_labels);
	arr_out_labels.num = CM.num_labels;
	// collect info about duplicated triangles
	dupl_triangles     = std::move(CM.dupl_triangles);
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

#ifdef OMC_ARR_AUX_LPI
	// clang-format off
	jolly_pnt.reserve(5);
	tri_soup.jolly_points.push_back(&jolly_pnt.emplace_back(0.94280904158, 0.0, -0.333333333));
	tri_soup.jolly_points.push_back(&jolly_pnt.emplace_back(-0.47140452079, 0.81649658092, -0.333333333));
	tri_soup.jolly_points.push_back(&jolly_pnt.emplace_back(-0.47140452079, -0.81649658092, -0.333333333));
	tri_soup.jolly_points.push_back(&jolly_pnt.emplace_back(0.0, 0.0, 1.0));
	tri_soup.jolly_points.push_back(&jolly_pnt.emplace_back(1.0, 0.0, 0.0));
	// clang-format on
#endif
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
	if (input_meshes.size() == LABEL_NBIT)
		OMC_THROW_OUT_OF_RANGE("Input meshes for arrangments are too much, limit "
		                       "the number to less than {}.",
		                       LABEL_NBIT);
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

private:
	template <typename Points, typename Triangles, typename NT>
	void load(const Points &points, const Triangles &triangles,
	          const size_t label, std::vector<NT> &coords,
	          std::vector<index_t> &flat_tris, std::vector<size_t> &labels)
	{
		size_t p_off = coords.size() / 3; // prev num verts
		coords.resize(coords.size() + points.size() * 3);
		tbb::parallel_for(size_t(0), points.size(),
		                  [&coords, &points, &p_off](size_t p_id)
		                  {
			                  coords[(p_off + p_id) * 3]     = points[p_id][0];
			                  coords[(p_off + p_id) * 3 + 1] = points[p_id][1];
			                  coords[(p_off + p_id) * 3 + 2] = points[p_id][2];
		                  });

		size_t t_off = flat_tris.size() / 3; // prev num tris
		flat_tris.resize(flat_tris.size() + triangles.size() * 3);
		tbb::parallel_for(
		  size_t(0), triangles.size(),
		  [&flat_tris, &triangles, &t_off, &p_off](size_t t_id)
		  {
			  flat_tris[(t_off + t_id) * 3]     = p_off + triangles[t_id][0];
			  flat_tris[(t_off + t_id) * 3 + 1] = p_off + triangles[t_id][1];
			  flat_tris[(t_off + t_id) * 3 + 2] = p_off + triangles[t_id][2];
		  });

		size_t l_off = labels.size();
		labels.resize(labels.size() + triangles.size());
		std::fill(std::execution::par_unseq, labels.begin() + l_off, labels.end(),
		          label);
	}
};

} // namespace OMC