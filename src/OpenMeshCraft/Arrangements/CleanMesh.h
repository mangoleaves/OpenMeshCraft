#pragma once

#include "Utils.h"

namespace OMC {

/// @brief Repairs the input mesh using common routines to generate a clean
/// mesh. This serves as a pre-processing step in arrangements.
template <typename Traits>
class ArrCleanMesh
{
public: /* Traits ************************************************************/
	using NT = typename Traits::NT;

	using CollinearPoints3D  = typename Traits::CollinearPoints3D;

public: /* Constructors ******************************************************/
	ArrCleanMesh() = delete;

	ArrCleanMesh(const std::vector<NT>      &_in_coords,
	             const std::vector<index_t> &_in_tris,
	             const std::vector<size_t>  &_in_labels);

public: /* Data ***************************************************************/
	/* Input data */
	/// Coordinates of all points.
	const std::vector<NT>      &in_coords;
	/// Triangles of all meshes (triangle soups).
	const std::vector<index_t> &in_tris;
	/// Labels of all triangles.
	const std::vector<size_t>  &in_labels;

	/* Output data */
	/// Output coordinates of all points.
	std::vector<NT>      out_coords;
	/// Output triangles.
	std::vector<index_t> out_tris;
	/// Output labels of all triangles
	std::vector<Label>   out_labels;
	/// Number of output labels (i.e., number of different meshes)
	size_t               num_labels;

	/// Information about removed duplicate triangles (may be used again).
	std::vector<ArrDuplTriInfo> dupl_triangles;

public: /* Pipeline **********************************************************/
	void convertLabels();

	void mergeDuplicatedVertices();

	void removeDegenerateAndDuplicatedTriangles();

	void removeIsolatedVertices();

	static bool consistentWinding(const index_t *t0, const index_t *t1);
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "CleanMesh.inl"
#endif