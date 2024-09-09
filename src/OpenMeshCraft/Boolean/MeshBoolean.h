#pragma once

// Definitions from mesh arrangements
#include "OpenMeshCraft/Arrangements/Utils.h"

namespace OMC {

enum class MeshBooleanOp
{
	UNION,
	INTERSECTION,
	SUBTRACTION,
	XOR,
	NONE
};

template <typename Traits>
class MeshBoolean_Impl;

template <typename Kernel, typename Traits>
class MeshBoolean
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
	MeshBoolean(bool _verbose);
	~MeshBoolean();

	MeshBoolean(const MeshBoolean &) = delete;
	MeshBoolean(MeshBoolean &&)      = delete;

public: /* Interfaces ********************************************************/
	/**
	 * @brief add a triangle mesh (triangle soup) as one input.
	 * @param points points of the mesh.
	 * @param triangles triangles of the mesh.
	 * @return index_t The internal index of the added mesh. It will used in
	 * labels of output triangle mesh.
	 * @note We only store pointers to input, won't store a copy, make sure the
	 * pointers are always valid during boolean.
	 */
	index_t addTriMeshAsInput(const iPoints &points, const iTriangles &triangles);

	/**
	 * @brief Set the triangle mesh (triangle soup) as output destination.
	 * @param points points of the output mesh.
	 * @param triangles triangles of the output mesh.
	 * @note We store pointers to output and write to them without check.
	 */
	void setTriMeshAsOutput(iPoints &points, iTriangles &triangles);

	/**
	 * @brief Set the output labels.
	 * @details
	 * * Each input triangle mesh will be attached with an index internally
	 *   (the index is the return value of addTriMeshAsInput).
	 * * Each output triangle's label indicates the triangle mesh it locates on.
	 * * ArrLabels are used to transfer attributes from input triangle meshes to
	 * output triangle meshes.
	 * @param labels labels of the output mesh.
	 */
	void setOutputLabels(std::vector<Label> &labels);

	///////////////////////////////////////////////////////////////////////
	/// Complete boolean pipeline contains four steps:                  ///
	/// (1) Run arrangements on all triangles.                          ///
	/// (2) Label triangles.                                            ///
	/// (3) apply boolean operation based on label.                     ///
	/// (4) compute explicit result and output.                         ///
	///////////////////////////////////////////////////////////////////////

	/// @brief Run (1)arrangements and (2)labeling. Store arrangements results and
	/// labels internally. (Then, users can apply different boolean operations and
	/// get final results faster.)
	void computeLabels();

	/// @brief apply (3)UNION and (4)put result to output mesh.
	void Union();
	/// @brief apply (3)INTERSECTION and (4)put result to output mesh.
	void Intersection();
	/// @brief apply (3)SUBTRACTION and (4)put result to output mesh.
	void Subtraction();
	/// @brief apply (3)XOR and (4)put result to output mesh.
	void Xor();

	/// @brief Clear input meshes, output mesh, and all internal data.
	/// Clearing is equivalent to reseting all to initial state.
	void clear();

private:
	struct iTriSoup
	{
		const iPoints    *points;
		const iTriangles *triangles;
	};
	std::vector<iTriSoup> input_meshes;

	iPoints    *output_points    = nullptr;
	iTriangles *output_triangles = nullptr;

	std::vector<Label> *output_labels = nullptr;

private: /* Internal boolean traits *******************************************/
	class BooleanTraits;
	/// Implement class
	std::unique_ptr<MeshBoolean_Impl<BooleanTraits>> m_impl;

private: /* Auxiliary classes and functions ***********************************/
	// Behavior control flags
	bool verbose;

	class loadMultipleMeshes;

	void booleanPipeline(MeshBooleanOp op);
};

// forward declaration
class ExactIndirectPredicatesApproxConstructions;
using EIAC = ExactIndirectPredicatesApproxConstructions;
class TriSoupTraits;

extern template class MeshBoolean<EIAC, TriSoupTraits>;

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "MeshBoolean.inl"
#endif