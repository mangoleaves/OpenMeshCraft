#pragma once

#include "STLReader.h"

#include "OpenMeshCraft/Utils/Macros.h"

#include <cfloat>
#include <map>

namespace OMC {

template <typename Traits>
class STLReader<Traits>::CmpVec
{
public:
	CmpVec(double _eps = FLT_MIN)
	  : eps_(_eps)
	{
	}

	bool operator()(const PointT &_v0, const PointT &_v1) const
	{
		if (fabs(_v0[0] - _v1[0]) <= eps_)
		{
			if (fabs(_v0[1] - _v1[1]) <= eps_)
			{
				return (_v0[2] < _v1[2] - eps_);
			}
			else
				return (_v0[1] < _v1[1] - eps_);
		}
		else
			return (_v0[0] < _v1[0] - eps_);
	}

private:
	double eps_;
};

template <typename Traits>
int STLReader<Traits>::read_int(FILE *_in, bool _swap)
{
	union u2
	{
		int           i;
		unsigned char c[4];
	} ic;
	OMC_UNUSED size_t rc = fread((char *)ic.c, 1, 4, _in);
	if (_swap)
	{
		std::swap(ic.c[0], ic.c[3]);
		std::swap(ic.c[1], ic.c[2]);
	}
	return ic.i;
}

template <typename Traits>
int STLReader<Traits>::read_int(std::istream &_in, bool _swap)
{
	union u2
	{
		int           i;
		unsigned char c[4];
	} ic;
	_in.read((char *)ic.c, 4);
	if (_swap)
	{
		std::swap(ic.c[0], ic.c[3]);
		std::swap(ic.c[1], ic.c[2]);
	}
	return ic.i;
}

template <typename Traits>
float STLReader<Traits>::read_float(FILE *_in, bool _swap)
{
	union u3
	{
		float         f;
		unsigned char c[4];
	} fc;
	OMC_UNUSED size_t rc = fread((char *)fc.c, 1, 4, _in);
	if (_swap)
	{
		std::swap(fc.c[0], fc.c[3]);
		std::swap(fc.c[1], fc.c[2]);
	}
	return fc.f;
}

template <typename Traits>
float STLReader<Traits>::read_float(std::istream &_in, bool _swap)
{
	union u3
	{
		float         f;
		unsigned char c[4];
	} fc;
	_in.read((char *)fc.c, 4);
	if (_swap)
	{
		std::swap(fc.c[0], fc.c[3]);
		std::swap(fc.c[1], fc.c[2]);
	}
	return fc.f;
}

template <typename Traits>
auto STLReader<Traits>::check_stl_type(const std::string &filename) -> STL_Type
{
	// assume it's binary stl, then file size is known from #triangles
	// if size matches, it's really binary

	// open file
	FILE *in = fopen(filename.c_str(), "rb");
	if (!in)
		return STL_Type::NONE;

	// determine endian mode
	union
	{
		unsigned int  i;
		unsigned char c[4];
	} endian_test;
	endian_test.i = 1;
	bool swapFlag = (endian_test.c[3] == 1);

	// read number of triangles
	char              dummy[100];
	OMC_UNUSED size_t rc = fread(dummy, 1, 80, in);
	size_t            nT = read_int(in, swapFlag);

	// compute file size from nT
	size_t binary_size = 84 + nT * 50;

	// get actual file size
	size_t file_size(0);
	rewind(in);
	while (!feof(in))
		file_size += fread(dummy, 1, 100, in);
	fclose(in);

	// if sizes match -> it's STLB
	return (binary_size == file_size ? STL_Type::STLB : STL_Type::STLA);
}

template <typename Traits>
bool STLReader<Traits>::read(const std::string &filename, IOOptions &opt)
{
	std::fstream ifile(filename.c_str(), std::ios_base::in);

	OMC_THROW_INVALID_ARGUMENT_IF(!ifile.is_open() || !ifile.good(),
	                              "[STLReader] : cannot not open file {}",
	                              filename);

	{
#if defined(WIN32)
		std::string::size_type dot_pos = filename.find_last_of("\\/");
#else
		std::string::size_type dot_pos = filename.rfind("/");
#endif
		path = (dot_pos == std::string::npos)
		         ? "./"
		         : std::string(filename.substr(0, dot_pos + 1));
	}

	bool result = false;

	STL_Type file_type = STL_Type::NONE;

	if (filename.rfind(".stla") != std::string::npos)
		file_type = STL_Type::STLA;
	else if (filename.rfind(".stlb") != std::string::npos)
		file_type = STL_Type::STLB;
	else if (filename.rfind(".stl") != std::string::npos)
		file_type = check_stl_type(filename);

	switch (file_type)
	{
	case STL_Type::STLA:
	{
		result         = read_stla(filename, opt);
		opt.stl_binary = false;
		break;
	}
	case STL_Type::STLB:
	{
		result         = read_stlb(filename, opt);
		opt.stl_binary = true;
		break;
	}
	default:
	{
		result = false;
		break;
	}
	}

	return result;
}

template <typename Traits>
inline bool STLReader<Traits>::read_stla(const std::string &filename,
                                         IOOptions         &opt)
{
	std::fstream in(filename.c_str(), std::ios_base::in);

	bool res = read_stla(in, opt);

	if (in)
		in.close();

	return res;
}

template <typename Traits>
inline bool STLReader<Traits>::read_stla(std::istream         &_in,
                                         OMC_UNUSED IOOptions &opt)
{
	unsigned int i;
	PointT       v;
	NormalT      n;
	Points       unique_points;
	Triangles    faces_data;

	using PointIDMap     = std::map<PointT, std::size_t, CmpVec>;
	using PointIDMapIter = typename PointIDMap::iterator;
	CmpVec     comp(FLT_MIN);
	PointIDMap point_to_id(comp);

	std::string line;

	bool facet_normal(false);

	while (_in && !_in.eof())
	{
		// Get one line
		std::getline(_in, line);
		OMC_THROW_INVALID_ARGUMENT_IF(_in.bad(),
		                              "Warning! Could not read stream properly!");

		// Trim Both leading and trailing spaces
		OMC::trim_string(line);

		// Normal found?
		if (line.find("facet normal") != std::string::npos)
		{
			std::stringstream strstream(line);
			strstream.precision(16);

			std::string garbage;

			// facet
			strstream >> garbage;

			// normal
			strstream >> garbage;

			strstream >> n[0];
			strstream >> n[1];
			strstream >> n[2];

			facet_normal = true;
		}

		// Detected a triangle
		if ((line.find("outer") != std::string::npos) ||
		    (line.find("OUTER") != std::string::npos))
		{
			Triangle indices;

			for (i = 0; i < 3; ++i)
			{
				// Get one vertex
				std::getline(_in, line);
				OMC::trim_string(line);

				std::stringstream strstream;
				strstream.precision(16);

				strstream << line;

				std::string garbage;
				strstream >> garbage;

				strstream >> v[0];

				strstream >> v[1];

				strstream >> v[2];

				// has vector been referenced before?
				std::pair<PointIDMapIter, bool> is_insert_successful =
				  point_to_id.insert(std::make_pair(v, unique_points.size()));
				std::size_t id = is_insert_successful.first->second;

				if (id == unique_points.size())
					unique_points.push_back(v);
				/*else
				{
				  size_t pos = point_to_id[v];
				  if (v != unique_points[pos])
				  {
				    std::cout.precision(12);
				    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
				    std::cout << unique_points[pos][0] << " " << unique_points[pos][1]
				              << " " << unique_points[pos][2] << std::endl
				              << std::endl;
				  }
				}*/
				indices[i] = static_cast<ti_t>(id);
			}

			if (facet_normal)
			{
				m_normals.push_back(n);
			}

			m_triangles.push_back(indices);
		}
	}
	m_points = unique_points;

	return true;
}

template <typename Traits>
inline bool STLReader<Traits>::read_stlb(const std::string &filename,
                                         IOOptions         &opt)
{
	std::fstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);

	OMC_THROW_INVALID_ARGUMENT_IF(!in, "[STLReader] : cannot not open file {}",
	                              filename);

	bool res = read_stlb(in, opt);

	if (in)
		in.close();

	return res;
}

template <typename Traits>
inline bool STLReader<Traits>::read_stlb(std::istream         &_in,
                                         OMC_UNUSED IOOptions &opt)
{
	char         dummy[100];
	bool         swapFlag;
	unsigned int i, nT;
	PointT       v;
	NormalT      n;

	Points unique_points;

	using PointIDMap = std::map<PointT, std::size_t, CmpVec>;
	using PointIDMapIter =
	  typename std::map<PointT, std::size_t, CmpVec>::iterator;
	CmpVec         comp(FLT_MIN);
	PointIDMap     point_to_id(comp);
	PointIDMapIter vMapIt;

	// check size of types
	static_assert(sizeof(float) == 4 && sizeof(int) == 4);

	// determine endian mode
	union
	{
		unsigned int  i;
		unsigned char c[4];
	} endian_test;
	endian_test.i = 1;
	swapFlag      = (endian_test.c[3] == 1);

	// read number of triangles
	_in.read(dummy, 80);
	nT = read_int(_in, swapFlag);

	// read triangles
	while (nT)
	{
		Triangle indices;
		// read triangle normal
		n[0] = read_float(_in, swapFlag);
		n[1] = read_float(_in, swapFlag);
		n[2] = read_float(_in, swapFlag);

		// triangle's vertices
		for (i = 0; i < 3; ++i)
		{
			v[0] = read_float(_in, swapFlag);
			v[1] = read_float(_in, swapFlag);
			v[2] = read_float(_in, swapFlag);

			// has vector been referenced before?
			std::pair<PointIDMapIter, bool> is_insert_successful =
			  point_to_id.insert(std::make_pair(v, unique_points.size()));
			std::size_t id = is_insert_successful.first->second;

			if (id == unique_points.size())
				unique_points.push_back(v);
			/*else
			{
			  size_t pos = point_to_id[v];
			  if (v != unique_points[pos])
			  {
			    std::cout.precision(12);
			    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
			    std::cout << unique_points[pos][0] << " " << unique_points[pos][1]
			              << " " << unique_points[pos][2] << std::endl
			              << std::endl;
			  }

			}*/
			indices[i] = static_cast<ti_t>(id);
		}

		// Add face only if it is not degenerated
		m_triangles.push_back(indices);
		m_normals.push_back(n);

		_in.read(dummy, 2);
		--nT;
	}
	m_points = unique_points;

	return true;
}

} // namespace OMC