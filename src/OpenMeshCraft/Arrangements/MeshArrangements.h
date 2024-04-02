#pragma once

#include "OpenMeshCraft/Arrangements/Utils.h"

// Triangle soup
#include "OpenMeshCraft/Mesh/TriSoup.h"
// Utils
#include "OpenMeshCraft/Utils/IndexDef.h"
#include "OpenMeshCraft/Utils/Macros.h"

#include <memory>

namespace OMC {

class ExactIndirectPredicatesApproxConstructions;
using EIAC = ExactIndirectPredicatesApproxConstructions;

template <typename Traits>
class MeshArrangements_Impl;

template <typename Kernel, typename Traits>
class MeshArrangements
{
public: /* Traits about input data ********************************************/
	/** @brief Point type.
	 * E.g., Vec3d from any library or a wrapper of it.
	 */
	using iPoint    = typename Traits::PointT;
	/** @brief Triangle type. */
	using iTriangle = typename Traits::Triangle;
	/** @brief Number type in PointT.
	 * We also expect that PointT provides x(), y() and z().
	 */
	using iNT  = std::remove_reference_t<decltype(std::declval<iPoint>().x())>;
	/** @brief index in triangle*/
	using ti_t = std::remove_reference_t<decltype(std::declval<iTriangle>()[0])>;
	/** @brief An array that stores points for vertices.
	 * E.g., std::vector<PointT> or Eigen::VectorX<PointT>.
	 */
	using iPoints    = typename Traits::Points;
	/** @brief An array that stores indices of three vertices for all faces.
	 * E.g., std::vector<std::array<index_t,3>> or
	 * Eigen::VectorX<Eigen::Vector3i>.
	 */
	using iTriangles = typename Traits::Triangles;

public: /* Constructors ******************************************************/
	MeshArrangements(bool _verbose);
	~MeshArrangements();

	MeshArrangements(const MeshArrangements &) = delete;
	MeshArrangements(MeshArrangements &&)      = delete;

public: /* Interfaces ********************************************************/
	/**
	 * @brief add a triangle mesh (triangle soup) as one input.
	 * @param points points of the mesh.
	 * @param triangles triangles of the mesh.
	 * @return index_t The internal index of the just added mesh.
	 * @note We only store pointers to input, won't store a copy.
	 */
	index_t addTriMeshAsInput(const iPoints &points, const iTriangles &triangles);

	/**
	 * @brief Set the triangle mesh (triangle soup) as output destination.
	 * @param points points of the mesh.
	 * @param triangles triangles of the mesh.
	 * @note We store pointers to output and write to them without check.
	 */
	void setTriMeshAsOutput(iPoints &points, iTriangles &triangles);

	/**
	 * @brief Set the output labels.
	 * @details
	 * * Each input triangle mesh will be attached with an index internally
	 *   (the index is the return value of addTriMeshAsInput).
	 * * Each output triangle's label indicates the triangle mesh it locates on.
	 * * Labels are used to transfer attributes from input triangle meshes to
	 * output triangle meshes.
	 * @param labels labels of the output mesh.
	 */
	void setOutputLabels(std::vector<Label> &labels);

	/// @brief Clear input meshes and output mesh.
	void clear();

	/// @brief Apply arrangement operation on input meshes and save explicit
	/// result to output mesh.
	/// @param ignore_intersection_in_same_mesh If set to true, algorithm will
	/// ignore intersection between triangles in the same mesh. This is a feature
	/// used by mesh boolean.
	/// @param output_explicit_result If set to true, the explicit result (points
	/// and triangles) will be saved in output mesh set by setTriMeshAsOutput.
	/// All middle data will be cleared.
	void meshArrangements(bool ignore_intersection_in_same_mesh,
	                      bool output_explicit_result);

	/// @brief An experimental interface to set parameters.
	/// @note Don't modify it except you know it.
	void setConfig(MeshArrangements_Config _config);

	MeshArrangements_Stats &stats();

private:
	/// Input data
	struct iTriSoup
	{
		const iPoints    *points    = nullptr;
		const iTriangles *triangles = nullptr;
	};
	std::vector<iTriSoup> input_meshes;
	/// Output data
	iPoints              *output_points    = nullptr;
	iTriangles           *output_triangles = nullptr;

	std::vector<Label> *output_labels = nullptr;

	/* Parameters */
	MeshArrangements_Config config;

	/// behavior control flags
	bool                   verbose;
	MeshArrangements_Stats arr_stats;

private:
	class ArrangementsTraits;

	/// Implement class
	std::unique_ptr<MeshArrangements_Impl<ArrangementsTraits>> m_impl;

	class loadMultipleMeshes;
};

extern template class MeshArrangements<EIAC, TriSoupTraits>;

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "MeshArrangements.inl"
#endif