#pragma once

#include "IOOptions.h"
#include "OpenMeshCraft/Mesh/TriSoup.h"

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/Macros.h"
#include "OpenMeshCraft/Utils/StringUtils.h"

#include <array>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

namespace OMC {

/**
 * @brief Write triangle soup to an OBJ file.
 * @tparam Traits Traits Define types for points, normals triangles.
 * OBJWriter will definen variables of these types to store results.
 */
template <typename Traits>
class OBJWriter
{
public:
	/** @brief Point type. */
	using PointT    = typename Traits::PointT;
	/** @brief Vector type. */
	using NormalT   = typename Traits::NormalT;
	/** @brief texture 3D. */
	using Tex3D     = typename Traits::Tex3D;
	/** @brief color. */
	using Color     = typename Traits::Color;
	/** @brief triangle. */
	using Triangle  = typename Traits::Triangle;
	/** @brief An array that stores points for vertices.
	 * E.g., std::vector<PointT> or Eigen::VectorX<PointT>.
	 */
	using Points    = typename Traits::Points;
	/** @brief An array that stores normals.
	 * E.g., std::vector<NormalT> or Eigen::VectorX<NormalT>.
	 */
	using Normals   = typename Traits::Normals;
	/** @brief An array that stores texture 3d coordinates. */
	using Tex3Ds    = typename Traits::Tex3Ds;
	/** @brief An array that stores colors . */
	using Colors    = typename Traits::Colors;
	/** @brief An array that stores indices of three vertices for all faces.
	 * E.g., std::vector<std::array<index_t,3>> or
	 * Eigen::VectorX<Eigen::Vector3i>.
	 */
	using Triangles = typename Traits::Triangles;

	using pn_t = std::remove_reference_t<decltype(std::declval<PointT>()[0])>;
	using nn_t = std::remove_reference_t<decltype(std::declval<NormalT>()[0])>;
	using tn_t = std::remove_reference_t<decltype(std::declval<Tex3D>()[0])>;
	using cn_t = std::remove_reference_t<decltype(std::declval<Color>()[0])>;

public:
	/**
	 * @brief Write triangle soup to file with given options and precison.
	 * @param opt Options to control which parts of triangle soup are written to
	 * file.
	 * @param precision precision of floating point numbers.
	 * @return true if succeed to write the triangle soup.
	 * @note Triangle soup must be stored in member variables of writer before
	 * writing.
	 */
	bool write(const std::string &filename, IOOptions &opt,
	           std::streamsize precision);

	/** @brief Clear data stored in writer */
	void clear();

	Points    m_points;    ///< point position
	Normals   m_normals;   ///< point normal
	Tex3Ds    m_texcoords; ///< point tex coordinates
	Triangles m_triangles; ///< triangle faces

private:
	std::string get_extensions() const { return "obj"; }

	bool save_to_stream(std::ostream &out);

private:
	std::string path, obj_name;

	IOOptions user_options;
};

extern template class OBJWriter<TriSoupTraits>;

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "OBJWriter.inl"
#endif