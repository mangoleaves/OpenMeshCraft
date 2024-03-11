#pragma once

#include "IOOptions.h"
#include "OpenMeshCraft/Mesh/TriSoup.h"

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/Macros.h"
#include "OpenMeshCraft/Utils/StringUtils.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

namespace OMC {

/**
 * @brief Read triangle soup from an OBJ file.
 * @tparam Traits Define types for points, normals triangles.
 * OBJReader will definen variables of these types to store results.
 */
template <typename Traits>
class OBJReader
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
	using ti_t = std::remove_reference_t<decltype(std::declval<Triangle>()[0])>;

public:
	OBJReader(){};
	~OBJReader(){};

	/**
	 * @brief Read triangle soup from the file. results will be stored in member
	 * variables of reader, e.g., \ref points, \ref faces.
	 * @param[in] filename path/name of the file.
	 * @param[inout] opt Options given by user, used to alternately read parts of
	 * file. File may not provide all options required by user, so, the common
	 * parts are returned.
	 * @return true if succeed to read the file.
	 */
	bool read(const std::string &filename, IOOptions &opt);

	Points    m_points;    ///< point position
	Normals   m_normals;   ///< point normal
	Colors    m_colors;    ///< point color
	Tex3Ds    m_texcoords; ///< point tex coordinates
	Triangles m_triangles; ///< triangle faces

private:
	std::string get_extensions() const { return "obj"; }

	/**
	 * @brief Check if reader is able to open/read the file.
	 * @param filename path/name of the file.
	 * @return true if reader can read the file.
	 */
	bool can_read(const std::string &filename) const;

	/**
	 * @brief After reading file as stream, load triangle soup from the stream.
	 * @param in input stream
	 * @param opt options given by user, used to alternately read part of triangle
	 * soup from file.
	 * @return true if succeed to read the triangle soup.
	 */
	bool load_from_stream(std::istream &in, IOOptions &opt);

	/**
	 * @brief Read vertices from the file, including points, texture coordinates,
	 * normals.
	 * @param in input stream
	 * @return true if succeed to read vertices.
	 */
	bool read_vertices(std::istream &in);

	/// @brief When key_word is v, process data store in stream.
	bool read_v(std::stringstream &stream);
	/// @brief When key_word is vt, process data store in stream.
	bool read_vt(std::stringstream &stream);
	/// @brief When key_word is vc, process data store in stream.
	bool read_vc(std::stringstream &stream);
	/// @brief When key_word is vn, process data store in stream.
	bool read_vn(std::stringstream &stream);

	/**
	 * @brief Read faces from the file. A line of face data contains multiple
	 * blocks (v0 v1 v2), each blocks may contains multiple components
	 * (vertex/texcoord/normal).
	 * @param in input stream.
	 * @return
	 */
	bool read_faces(std::istream &in);

	/**
	 * @brief Read a face data from a line of obj file.
	 * @param stream a line of obj file, containing multiple blocks
	 * (vertex/texCoord/normal).
	 */
	void read_a_face(std::stringstream &stream);

	/**
	 * @brief Process a block of face data.
	 * @param block block string
	 * @param comp_value the index, may be negative.
	 * @param component indicates which component the index belongs to.
	 * @return true if succeed to get current component.
	 */
	bool get_component_from_block(std::string &block, size_t &comp_idx,
	                              int &comp_value);

	/// @brief when process a vertex block data of face, process vertex component.
	index_t set_face_vertex(int comp_value);

	/// @brief when process a vertex block data of face, process texture
	/// coordinate component.
	void set_face_tex(int comp_value, index_t vert_idx);

	/// @brief when process a vertex block data of face, process normal component.
	void set_face_normal(int comp_value, index_t vert_idx);

	void postfix_when_no_faces();

private:
	std::string path;

	IOOptions user_options;
	IOOptions file_options;

	/// @name Used for temporarily store content of file.
	/// @{
	std::vector<PointT>  t_points;
	std::vector<NormalT> t_normals;
	std::vector<Color>   t_colors;
	std::vector<Tex3D>   t_texcoords;
	/// @}

	/// @name Used for building faces
	/// @{
	size_t n_point, n_tex, n_normal;
	/// @}
};

extern template class OBJReader<TriSoupTraits>;

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "OBJReader.inl"
#endif