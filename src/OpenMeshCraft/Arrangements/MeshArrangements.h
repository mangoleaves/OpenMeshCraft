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
	MeshArrangements();
	~MeshArrangements();

	MeshArrangements(const MeshArrangements &) = delete;
	MeshArrangements(MeshArrangements &&)      = delete;

public: /* Interfaces ********************************************************/
	/**
	 * @brief Adds a triangle mesh (triangle soup) as an input.
	 * @param points The points of the mesh.
	 * @param triangles The triangles of the mesh.
	 * @return index_t The internal index of the newly added mesh.
	 * @note We only store pointers to the input; we do not store a copy.
	 */
	index_t addTriMeshAsInput(const iPoints &points, const iTriangles &triangles);

	/**
	 * @brief Sets the triangle mesh (triangle soup) as the output destination.
	 * @param points The points of the mesh.
	 * @param triangles The triangles of the mesh.
	 * @note We store pointers to the output and write to them without checks.
	 */
	void setTriMeshAsOutput(iPoints &points, iTriangles &triangles);

	/**
	 * @brief Sets the output labels.
	 * @details
	 * * Each input triangle mesh will be assigned an index internally
	 *   (the index is the return value of addTriMeshAsInput).
	 * * Each output triangle's label indicates the triangle mesh it belongs to.
	 * * Labels are used to transfer attributes from input triangle meshes to
	 *   output triangle meshes.
	 * @param labels The labels of the output mesh.
	 */
	void setOutputLabels(std::vector<Label> &labels);

	/// @brief Clear input meshes and output mesh.
	void clear();

	/// @brief Apply arrangement operation on input meshes and
	/// (optionally) save explicit result to output mesh.
	void meshArrangements();

	/// @brief An interface to set configuration (flags and parameters).
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

	/* Configuration */
	MeshArrangements_Config config;

	/* Statistics */
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