#pragma once

#include "OBJReader.h"

namespace OMC {

template <typename Traits>
bool OBJReader<Traits>::can_read(const std::string &filename) const
{
	// get file extension
	std::string            extension;
	std::string::size_type pos(filename.rfind("."));

	if (pos != std::string::npos)
		extension = filename.substr(pos + 1, filename.length() - pos - 1);
	else
		extension = filename; // check, if the whole filename defines the extension

	auto tolower_ = [](char c)
	{ return static_cast<char>(std::tolower(static_cast<int>(c))); };

	std::transform(extension.begin(), extension.end(), extension.begin(),
	               tolower_);

	// locate extension in extension string
	return (get_extensions().find(extension) != std::string::npos);
}

template <typename Traits>
bool OBJReader<Traits>::read(const std::string &filename, IOOptions &opt)
{
	std::fstream ifile(filename.c_str(), std::ios_base::in);

	OMC_THROW_INVALID_ARGUMENT_IF(!ifile.is_open() || !ifile.good(),
	                              "[OBJReader] : cannot not open file {}",
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

	bool result = load_from_stream(ifile, opt);

	ifile.close();
	return result;
}

template <typename Traits>
bool OBJReader<Traits>::load_from_stream(std::istream &in, IOOptions &opt)
{
	OMC_THROW_INVALID_ARGUMENT_IF(!in.good(),
	                              "[OBJReader] : cannot not use stream.");

	// Options supplied by the user
	user_options = opt;
	file_options.clear();

	// pass 1: read vertices
	if (!read_vertices(in))
		return false;

	// pass 2: read faces
	if (!read_faces(in))
		return false;

	// If we do not have any faces,
	// assume this is a point cloud and read the normals and colors directly
	if (m_triangles.empty())
	{
		postfix_when_no_faces();
	}

	// Return, what we actually read
	opt = user_options.intersection(file_options);

	return true;
}

template <typename Traits>
bool OBJReader<Traits>::read_vertices(std::istream &in)
{
	t_points.clear();
	t_normals.clear();
	t_colors.clear();
	t_texcoords.clear();
	// reset stream to begin position.
	in.clear();
	in.seekg(0, std::ios::beg);

	std::string       line, key_word;
	std::stringstream stream;

	while (in && !in.eof())
	{
		std::getline(in, line);
		OMC_THROW_INVALID_ARGUMENT_IF(
		  in.bad(), "[OBJReader] : Could not read file properly!");

		trim_string(line);

		// comment
		if (line.size() == 0 || line[0] == '#' || std::isspace(line[0]))
			continue;

		stream.str(line);
		stream.clear();
		stream >> key_word;

		if (key_word == "v") // vertex point
		{
			if (!read_v(stream))
				return false;
		}
		else if (key_word == "vt") // vertex texture
		{
			if (!read_vt(stream))
				return false;
		}
		else if (key_word == "vc") // vertex color
		{
			if (!read_vc(stream))
				return false;
		}
		else if (key_word == "vn") // vertex normal
		{
			if (!read_vn(stream))
				return false;
		}
	}

	m_points.clear();
	m_normals.clear();
	m_colors.clear();
	m_texcoords.clear();
	if (user_options.vertex_has_point && file_options.vertex_has_point)
	{
		m_points.resize(t_points.size());
		for (index_t i = 0; i < t_points.size(); i++)
			m_points[i] = t_points[i];
	}
	if (user_options.vertex_has_normal && file_options.vertex_has_normal)
		m_normals.resize(t_points.size());
	if (user_options.vertex_has_color && file_options.vertex_has_color)
		m_colors.resize(t_points.size());
	if ((user_options.vertex_has_tex2d && file_options.vertex_has_tex2d) ||
	    (user_options.vertex_has_tex3d && file_options.vertex_has_tex3d))
		m_texcoords.resize(t_points.size());
	return true;
}

template <typename Traits>
bool OBJReader<Traits>::read_v(std::stringstream &stream)
{
	pn_t x, y, z;
	cn_t r, g, b;
	stream >> x >> y >> z;
	file_options.vertex_has_point = true;
	if (!stream.fail())
	{
		if (user_options.vertex_has_point)
			t_points.push_back(PointT{x, y, z});
		else
			t_points.emplace_back();

		stream >> r >> g >> b;
		if (!stream.fail())
		{
			file_options.vertex_has_color = true;
			t_colors.push_back(Color{r, g, b});
		}
	}
	return true;
}

template <typename Traits>
bool OBJReader<Traits>::read_vt(std::stringstream &stream)
{
	tn_t u, v, w;
	stream >> u >> v;
	if (!stream.fail())
	{
		t_texcoords.push_back(Tex3D{u, v, tn_t(0.)});
		file_options.vertex_has_tex2d = true;

		stream >> w;
		if (!stream.fail())
		{
			file_options.vertex_has_tex3d = true;
			t_texcoords.back().z()        = w;
		}
		return true;
	}
	else
	{
		OMC_THROW_INVALID_ARGUMENT("[OBJReader] : Only single 2D or 3D "
		                           "texture coordinate per vertex allowed.");
	}
}

template <typename Traits>
bool OBJReader<Traits>::read_vc(std::stringstream &stream)
{
	cn_t r, g, b;
	stream >> r >> g >> b;
	if (!stream.fail())
	{
		t_colors.push_back(Color{r, g, b});
		file_options.vertex_has_color = true;
	}
	return true;
}

template <typename Traits>
bool OBJReader<Traits>::read_vn(std::stringstream &stream)
{
	nn_t x, y, z;
	stream >> x >> y >> z;
	if (!stream.fail())
	{
		t_normals.push_back(NormalT{x, y, z});
		file_options.vertex_has_normal = true;
	}
	return true;
}

template <typename Traits>
bool OBJReader<Traits>::read_faces(std::istream &in)
{
	// reset stream
	in.clear();
	in.seekg(0, std::ios::beg);

	std::string line, key_word;

	n_point = 0, n_tex = 0, n_normal = 0;
	m_triangles.clear();

	while (in && !in.eof())
	{
		std::getline(in, line);
		OMC_THROW_INVALID_ARGUMENT_IF(
		  in.bad(), "[OBJReader] : Could not read file properly.");

		trim_string(line);

		// comment
		if (line.size() == 0 || line[0] == '#' || isspace(line[0]))
		{
			continue;
		}

		std::stringstream stream;
		stream.str(line);
		stream.clear();
		stream >> key_word;

		// material file        "mtllib"
		// usemtl               "usemtl"
		// faces                "f"
		if (key_word == "mtllib")
		{
			// TODO
		}
		else if (key_word == "usemtl")
		{
			// TODO
		}
		// track current number of parsed vertex attributes, to allow for OBJs
		// negative indices
		else if (key_word == "v")
		{
			++n_point;
		}
		else if (key_word == "vt")
		{
			++n_tex;
		}
		else if (key_word == "vn")
		{
			++n_normal;
		}
		else if (key_word == "f")
		{
			// read full line after detecting a face
			std::string       face_line;
			std::stringstream line_data;

			std::getline(stream, face_line);
			line_data.str(face_line);
			line_data.clear();

			read_a_face(line_data);
		}
	}

	t_points.clear();
	t_normals.clear();
	t_colors.clear();
	t_texcoords.clear();
	return true;
}

template <typename Traits>
void OBJReader<Traits>::read_a_face(std::stringstream &stream)
{
	Triangle face;
	index_t  local_idx = 0;

	// work on the line until nothing left to read
	while (!stream.eof() && local_idx < 3)
	{
		// read one block from the line ( vertex/texCoord/normal )
		std::string block;
		stream >> block;

		index_t comp_idx = 0;
		int     comp_value;
		index_t vert_idx = InvalidIndex;
		do
		{
			// store the component ( each component is referenced by the index
			// here! )
			if (!get_component_from_block(block, comp_idx, comp_value))
				continue;

			// Calculation of index : -1 is the last vertex in the list
			// As obj counts from 1 and not zero add +1
			switch (comp_idx)
			{
			case 0: // vertex
				vert_idx          = set_face_vertex(comp_value);
				face[local_idx++] = static_cast<ti_t>(vert_idx);
				break;
			case 1: // texture coord
				OMC_THROW_INVALID_ARGUMENT_IF(local_idx == 0, "empty triangle face");
				set_face_tex(comp_value, vert_idx);
				break;
			case 2: // normal
				OMC_THROW_INVALID_ARGUMENT_IF(local_idx == 0, "empty triangle face");
				set_face_normal(comp_value, vert_idx);
				break;
			}
			// Prepare for reading next component
			++comp_idx;
			// Read until line does not contain any other info
		} while (!block.empty());
	}

	m_triangles.push_back(face);
}

template <typename Traits>
bool OBJReader<Traits>::get_component_from_block(std::string &block,
                                                 size_t      &comp_idx,
                                                 int         &comp_value)
{
	std::stringstream tmp;

	// get the component (vertex/texCoord/normal)
	size_t found = block.find("/");

	// parts are seperated by '/' So if no '/' found its the last
	// component
	if (found != std::string::npos)
	{
		// read the index value
		tmp.str(block.substr(0, found));
		tmp.clear();
		// If we get an empty string this property is undefined in the file
		if (block.substr(0, found).empty())
		{
			// Switch to next field
			block = block.substr(found + 1);
			// Now we are at the next component
			++comp_idx;
			// Skip further processing of this component
			return false;
		}
		// Read current value
		tmp >> comp_value;
		// remove the read part from the string
		block = block.substr(found + 1);
		return true;
	}
	else
	{
		// last component of the block, read it.
		tmp.str(block);
		tmp.clear();
		tmp >> comp_value;
		// Clear block after finished reading the line
		block = "";
		// Nothing to read here ( garbage at end of line )
		return !tmp.fail();
	}
}

template <typename Traits>
index_t OBJReader<Traits>::set_face_vertex(int comp_value)
{
	// Calculation of index :
	// -1 is the last vertex in the list
	// As obj counts from 1 and not zero add +1
	// Obj counts from 1 and not zero .. array counts from zero therefore -1
	size_t vert_idx =
	  (size_t)(comp_value < 0 ? n_point + comp_value + 1 : comp_value) - 1;

	if (file_options.vertex_has_color && user_options.vertex_has_color)
	{
		if (vert_idx < t_points.size() && vert_idx < t_colors.size())
			m_colors[vert_idx] = t_colors[vert_idx];
		else
			OMC_THROW_OUT_OF_RANGE(
			  "[OBJReader] : Error setting vertex color, index out of range.");
	}
	return vert_idx;
}

template <typename Traits>
void OBJReader<Traits>::set_face_tex(int comp_value, index_t vert_idx)
{
	size_t tex_idx =
	  (size_t)(comp_value < 0 ? n_tex + comp_value + 1 : comp_value) - 1;

	if (file_options.vertex_has_tex2d && user_options.vertex_has_tex2d)
	{
		if (tex_idx < t_texcoords.size())
			m_texcoords[vert_idx] =
			  Tex3D(t_texcoords[tex_idx].x(), t_texcoords[tex_idx].y(), tn_t(0));
		else
			OMC_THROW_OUT_OF_RANGE(
			  "[OBJReader] : Error setting Texture coordinates, index out of range.");
	}

	if (file_options.vertex_has_tex3d && user_options.vertex_has_tex3d)
	{
		if (tex_idx < t_texcoords.size())
			m_texcoords[vert_idx] = t_texcoords[tex_idx];
		else
			OMC_THROW_OUT_OF_RANGE(
			  "[OBJReader] : Error setting Texture coordinates, index out of range.");
	}
}

template <typename Traits>
void OBJReader<Traits>::set_face_normal(int comp_value, index_t vert_idx)
{
	size_t normal_idx =
	  (size_t)(comp_value < 0 ? n_normal + comp_value + 1 : comp_value) - 1;

	if (file_options.vertex_has_normal && user_options.vertex_has_normal)
	{
		if (normal_idx < t_normals.size())
			m_normals[vert_idx] = t_normals[normal_idx];
		else
			OMC_THROW_OUT_OF_RANGE(
			  "[OBJReader] : Error setting vertex normal, index out of range.");
	}
}

template <typename Traits>
void OBJReader<Traits>::postfix_when_no_faces()
{
	size_t vertnum = m_points.size();
	// add normal per vertex
	if (t_normals.size() == t_points.size())
	{
		if (user_options.vertex_has_normal && file_options.vertex_has_normal)
		{
			m_normals.resize(t_normals.size());
			for (size_t i = 0; i < vertnum; i++)
				m_normals[i] = t_normals[i];
		}
	}

	// add color per vertex
	if (t_colors.size() == t_points.size())
	{
		if (user_options.vertex_has_color && file_options.vertex_has_color)
		{
			m_colors.resize(t_colors.size());
			for (size_t i = 0; i < vertnum; i++)
				m_colors[i] = t_colors[i];
		}
	}
}

} // namespace OMC