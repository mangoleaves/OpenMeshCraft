#pragma once

#include "STLWriter.h"

#include "OpenMeshCraft/Utils/Macros.h"

#include <execution>

namespace OMC {

template <typename Traits>
void STLWriter<Traits>::write_int(int _i, FILE *_out, bool _swap)
{
	union u2
	{
		int           i;
		unsigned char c[4];
	} ic;
	ic.i = _i;
	if (_swap)
	{
		std::swap(ic.c[0], ic.c[3]);
		std::swap(ic.c[1], ic.c[2]);
	}
	fwrite((char *)ic.c, 1, 4, _out);
}

template <typename Traits>
void STLWriter<Traits>::write_float(float _f, FILE *_out, bool _swap)
{
	union u3
	{
		float         f;
		unsigned char c[4];
	} fc;
	fc.f = _f;
	if (_swap)
	{
		std::swap(fc.c[0], fc.c[3]);
		std::swap(fc.c[1], fc.c[2]);
	}
	fwrite((char *)fc.c, 1, 4, _out);
}

template <typename Traits>
void STLWriter<Traits>::write_int(int _i, std::ostream &_out, bool _swap)
{
	union u2
	{
		int           i;
		unsigned char c[4];
	} ic;
	ic.i = _i;
	if (_swap)
	{
		std::swap(ic.c[0], ic.c[3]);
		std::swap(ic.c[1], ic.c[2]);
	}
	_out.write((char *)ic.c, 4);
}

template <typename Traits>
void STLWriter<Traits>::write_float(float _f, std::ostream &_out, bool _swap)
{
	union u3
	{
		float         f;
		unsigned char c[4];
	} fc;
	fc.f = _f;
	if (_swap)
	{
		std::swap(fc.c[0], fc.c[3]);
		std::swap(fc.c[1], fc.c[2]);
	}
	_out.write((char *)fc.c, 4);
}

template <typename Traits>
void STLWriter<Traits>::write_short(short int _i, FILE *_out, bool _swap)
{
	union u1
	{
		short int     s;
		unsigned char c[2];
	} sc;
	sc.s = _i;
	if (_swap)
		std::swap(sc.c[0], sc.c[1]);
	fwrite((char *)sc.c, 1, 2, _out);
}

template <typename Traits>
void STLWriter<Traits>::write_short(short int _i, std::ostream &_out,
                                    bool _swap)
{
	union u1
	{
		short int     s;
		unsigned char c[2];
	} sc;
	sc.s = _i;
	if (_swap)
		std::swap(sc.c[0], sc.c[1]);
	_out.write((char *)sc.c, 2);
}

template <typename Traits>
auto STLWriter<Traits>::calc_normal(const Triangle &t) -> NormalT
{
	const PointT &p0     = m_points[t[0]];
	const PointT &p1     = m_points[t[1]];
	const PointT &p2     = m_points[t[2]];
	NormalT       result = (p1 - p0).cross(p2 - p0);
	result.normalize();
	return result;
}

template <typename Traits>
bool STLWriter<Traits>::write(const std::string &filename, IOOptions &opt,
                              std::streamsize precision)
{
	// binary or ascii ?
	if (filename.rfind(".stla") != std::string::npos)
	{
		opt.stl_binary = false;
	}
	else if (filename.rfind(".stlb") != std::string::npos)
	{
		opt.stl_binary = true;
	}

	// open file
	std::fstream out(filename.c_str(),
	                 (opt.stl_binary ? std::ios_base::binary | std::ios_base::out
	                                 : std::ios_base::out));

	OMC_THROW_INVALID_ARGUMENT_IF(!out, "[STLWriter] : cannot open file {}",
	                              filename);

	user_options = opt;

	// Set precision on output stream. The default is set via IOManager and passed
	// through to all writers.
	// out.precision(precision);

	// Set fixed output to avoid problems with programs not reading scientific
	// notation correctly
	out << std::fixed;

	{
#if defined(WIN32)
		std::string::size_type dotposition = filename.find_last_of("\\/");
#else
		std::string::size_type dotposition = filename.rfind("/");
#endif

		if (dotposition == std::string::npos)
		{
			path     = "./";
			stl_name = filename;
		}
		else
		{
			path     = filename.substr(0, dotposition + 1);
			stl_name = filename.substr(dotposition + 1);
		}

		// remove the file extension
		dotposition = stl_name.find_last_of(".");

		if (dotposition != std::string::npos)
			stl_name = stl_name.substr(0, dotposition);
	}

	{
		// check points
		// OMC_THROW_Null_Value_if(m_points.empty(), "[STLWriter] : Warning
		// no points");

		// fix triangle normal if neccessary
		if (m_normals.empty() || m_normals.size() != m_triangles.size())
		{
			m_normals.resize(m_triangles.size());
			std::transform(std::execution::par_unseq, m_triangles.begin(),
			               m_triangles.end(), m_normals.begin(),
			               [this](const Triangle &t) { return calc_normal(t); });
		}
	}

	bool result = save_mesh_to_stream(out, opt, precision);

	out.close();
	return result;
}

template <typename Traits>
inline void STLWriter<Traits>::clear()
{
	m_points.clear();
	m_normals.clear();
	m_triangles.clear();
}

template <typename Traits>
bool STLWriter<Traits>::save_mesh_to_stream(std::ostream &out, IOOptions &opt,
                                            std::streamsize precision)
{
	if (!opt.stl_binary)
		out.precision(precision);

	if (opt.stl_binary)
		return write_stlb(out, opt);
	else
		return write_stla(out, opt);
}

template <typename Traits>
inline bool STLWriter<Traits>::write_stla(std::ostream         &out,
                                          OMC_UNUSED IOOptions &opt,
                                          std::streamsize       _precision)
{
	size_t i, nF(m_triangles.size()), nV;

	out.precision(_precision);

	// header
	out << "solid\n";

	// write face set
	for (i = 0; i < nF; ++i)
	{
		const auto &face_vertices = m_triangles[i];
		nV                        = face_vertices.size();

		if (nV == 3)
		{
			const PointT  &a = m_points[face_vertices[0]];
			const PointT  &b = m_points[face_vertices[1]];
			const PointT  &c = m_points[face_vertices[2]];
			const NormalT &n = m_normals[i];

			out << "facet normal " << n[0] << " " << n[1] << " " << n[2]
			    << "\nouter loop\n";
			out.precision(10);
			out << "vertex " << a[0] << " " << a[1] << " " << a[2] << "\n";
			out << "vertex " << b[0] << " " << b[1] << " " << b[2] << "\n";
			out << "vertex " << c[0] << " " << c[1] << " " << c[2] << "\n";
		}
		// else OMC_THROW_DOMAIN_ERROR("[STLWriter] : Warning non-triangle
		// data!");

		out << "\nendloop\nendfacet\n";
	}

	return true;
}

template <typename Traits>
inline bool STLWriter<Traits>::write_stlb(std::ostream         &out,
                                          OMC_UNUSED IOOptions &opt,
                                          std::streamsize       _precision)
{
	size_t i, nF(m_triangles.size()), nV;
	out.precision(_precision);

	// header
	const char header[80] =
	  "binary stl file"
	  "                                                                ";
	out.write(header, 80);

	// number of faces
	write_int(int(m_triangles.size()), out);

	// write face set
	for (i = 0; i < nF; ++i)
	{
		const auto &face_vertices = m_triangles[i];
		nV                        = face_vertices.size();

		if (nV == 3)
		{
			const PointT  &a = m_points[face_vertices[0]];
			const PointT  &b = m_points[face_vertices[1]];
			const PointT  &c = m_points[face_vertices[2]];
			const NormalT &n = m_normals[i];

			// face normal
			write_float((float)n[0], out);
			write_float((float)n[1], out);
			write_float((float)n[2], out);

			// face vertices
			write_float((float)a[0], out);
			write_float((float)a[1], out);
			write_float((float)a[2], out);

			write_float((float)b[0], out);
			write_float((float)b[1], out);
			write_float((float)b[2], out);

			write_float((float)c[0], out);
			write_float((float)c[1], out);
			write_float((float)c[2], out);

			// space filler
			write_short(0, out);
		}
		// else OMC_THROW_DOMAIN_ERROR("[STLWriter] : Warning: Skipped
		// non-triangle data!");
	}

	return true;
}

} // namespace OMC