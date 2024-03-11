#pragma once

#include "IOOptions.h"
#include "OpenMeshCraft/Mesh/TriSoup.h"

#include "OpenMeshCraft/Utils/Exception.h"
#include "OpenMeshCraft/Utils/Macros.h"
#include "OpenMeshCraft/Utils/StringUtils.h"

#include <array>
#include <fstream>
#include <string>
#include <vector>

namespace OMC {

/**
 * @brief Read triangle soup from an STL file.
 * @tparam Traits Define types for points, normals triangles.
 * STLReader will definen variables of these types to store results.
 */
template <typename Traits>
class STLReader
{
public:
	/** @brief Point type. */
	using PointT   = typename Traits::PointT;
	/** @brief Vector type.
	 * It is derived from PointT, representing the difference between two points.
	 */
	using NormalT  = typename Traits::NormalT;
	/** @brief triangle type */
	using Triangle = typename Traits::Triangle;
	/** @brief Number type in PointT.
	 * We also expect that PointT provides x(), y() and z().
	 */
	using NT     = std::remove_reference_t<decltype(std::declval<PointT>().x())>;
	/** @brief index used in triangle */
	using ti_t   = std::remove_reference_t<decltype(std::declval<Triangle>()[0])>;
	/** @brief An array that stores points for vertices.
	 * E.g., std::vector<PointT> or Eigen::VectorX<PointT>.
	 */
	using Points = typename Traits::Points;
	/** @brief An array that stores normals.
	 * E.g., std::vector<NormalT> or Eigen::VectorX<NormalT>.
	 */
	using Normals   = typename Traits::Normals;
	/** @brief An array that stores indices of three vertices for all faces.
	 * E.g., std::vector<std::array<index_t,3>> or
	 * Eigen::VectorX<Eigen::Vector3i>.
	 */
	using Triangles = typename Traits::Triangles;

public:
	STLReader(){};
	~STLReader(){};

	Points    m_points;    ///< vertex position
	Triangles m_triangles; ///< triangles
	Normals   m_normals;   ///< triangle normals

	/**
	 * @brief Read triangle soup from the file.
	 * @param[inout] opt Options given by user, used to alternately read parts of
	 * triangle soup from file. File may not provide all options required by user,
	 * so, the common parts are returned.
	 * @return true if succeed to read the triangle soup.
	 */
	bool read(const std::string &filename, IOOptions &opt);

private:
	std::string get_extensions() const { return "stl"; }

	bool read_stla(const std::string &filename, IOOptions &opt);
	bool read_stla(std::istream &_in, IOOptions &opt);
	bool read_stlb(const std::string &filename, IOOptions &opt);
	bool read_stlb(std::istream &_in, IOOptions &opt);

private:
	class CmpVec;

	int   read_int(FILE *_in, bool _swap = false);
	int   read_int(std::istream &_in, bool _swap = false);
	float read_float(FILE *_in, bool _swap = false);
	float read_float(std::istream &_in, bool _swap = false);

	enum class STL_Type
	{
		STLA,
		STLB,
		NONE
	};
	STL_Type check_stl_type(const std::string &filename);

	std::string path;

	IOOptions useroptions;
	IOOptions fileoptions;

	/// @name Used for vertices
	/// @{
	std::vector<index_t> vertices;
	std::vector<NormalT> normals;
	/// @}

	double eps_;
};

extern template class STLReader<TriSoupTraits>;

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "STLReader.inl"
#endif