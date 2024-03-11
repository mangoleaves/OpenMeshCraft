#pragma once

#include "OBJWriter.h"

namespace OMC {

template <typename Traits>
bool OBJWriter<Traits>::write(const std::string &filename, IOOptions &opt,
                              std::streamsize precision)
{
	std::fstream out(filename.c_str(), std::ios_base::out);

	OMC_THROW_INVALID_ARGUMENT_IF(!out, "[OBJWriter] : cannot open file {}",
	                              filename);

	user_options = opt;

	// Set precision on output stream. The default is set via IOManager and passed
	// through to all writers.
	out.precision(precision);

	// Set fixed output to avoid problems with programs not reading scientific
	// notation correctly
	out << std::fixed;

	// check for unsupported writer features
	if (user_options.vertex_has_color)
	{
		OMC_THROW_NOT_IMPLEMENTED();
	}

	{
#if defined(WIN32)
		std::string::size_type dotposition = filename.find_last_of("\\/");
#else
		std::string::size_type dotposition = filename.rfind("/");
#endif

		if (dotposition == std::string::npos)
		{
			path     = "./";
			obj_name = filename;
		}
		else
		{
			path     = filename.substr(0, dotposition + 1);
			obj_name = filename.substr(dotposition + 1);
		}

		// remove the file extension
		dotposition = obj_name.find_last_of(".");

		if (dotposition != std::string::npos)
			obj_name = obj_name.substr(0, dotposition);
	}

	bool result = save_to_stream(out);

	out.close();
	return result;
}

template <typename Traits>
void OBJWriter<Traits>::clear()
{
	m_points.clear();
	m_normals.clear();
	m_texcoords.clear();
	m_triangles.clear();
}

template <typename Traits>
bool OBJWriter<Traits>::save_to_stream(std::ostream &out)
{
	bool write_vertex_point = user_options.vertex_has_point && !m_points.empty();
	bool write_vertex_normal =
	  user_options.vertex_has_normal && !m_normals.empty();
	bool write_vertex_tex2d =
	  user_options.vertex_has_tex2d && !m_texcoords.empty();
	// TODO: create material file if needed

	// header
	out << "# " << m_points.size() << " vertices, " << m_triangles.size()
	    << " faces" << '\n';

	// TODO: support texture coordinates from faces.

	// vertex data (point, normals, texcoords)
	for (size_t i = 0; i < m_points.size(); ++i)
	{
		if (write_vertex_point)
		{
			const PointT &v = m_points[i];
			out << "v " << v[0] << " " << v[1] << " " << v[2] << '\n';
		}
		else
		{
			out << "v 0 0 0" << '\n'; // avoid generate invalid obj files
		}

		if (write_vertex_normal)
		{
			const NormalT &n = m_normals[i];
			out << "vn " << n[0] << " " << n[1] << " " << n[2] << '\n';
		}

		if (write_vertex_tex2d)
		{
			const Tex3D &t = m_texcoords[i];
			out << "vt " << t[0] << " " << t[1] << '\n';
		}
	}

	// we do not want to write seperators if we only write vertex indices
	bool only_vertices = !write_vertex_normal && !write_vertex_tex2d;

	// m_triangles (indices starting at 1 not 0)
	for (size_t i = 0; i < m_triangles.size(); ++i)
	{
		out << "f";

		auto output_idx = [&out, &only_vertices, &write_vertex_tex2d,
		                   &write_vertex_normal](index_t idx)
		{
			// Write vertex index
			out << " " << idx;

			if (!only_vertices)
			{
				// write separator
				out << "/";

				// write vertex texture coordinate index
				if (write_vertex_tex2d)
				{
					out << idx;
				}

				// write vertex normal index
				if (write_vertex_normal)
				{
					// write separator
					out << "/";
					out << idx;
				}
			}
		};

		output_idx(m_triangles[i][0] + 1);
		output_idx(m_triangles[i][1] + 1);
		output_idx(m_triangles[i][2] + 1);

		out << '\n';
	}

	return true;
}

} // namespace OMC