#pragma once

namespace OMC {

/**
 * @brief Options about pre-described properties in mesh IO.
 */
class IOOptions
{
public:
	bool vertex_has_point;
	bool vertex_has_normal;
	bool vertex_has_tex2d;
	bool vertex_has_tex3d;
	bool vertex_has_color;

	bool face_has_normal;
	bool face_has_color;
	bool face_has_material;

	bool color_has_alpha;

	bool stl_binary;

	IOOptions()
	  : vertex_has_point(false)
	  , vertex_has_normal(false)
	  , vertex_has_tex2d(false)
	  , vertex_has_tex3d(false)
	  , vertex_has_color(false)
	  , face_has_normal(false)
	  , face_has_color(false)
	  , face_has_material(false)
	  , color_has_alpha(false)
	  , stl_binary(false)
	{
	}

	void clear()
	{
		vertex_has_point  = false;
		vertex_has_normal = false;
		vertex_has_tex2d  = false;
		vertex_has_tex3d  = false;
		vertex_has_color  = false;
		face_has_normal   = false;
		face_has_color    = false;
		face_has_material = false;
		color_has_alpha   = false;
		stl_binary        = false;
	}

	IOOptions &intersection(const IOOptions &rhs)
	{
		vertex_has_point  = vertex_has_point && rhs.vertex_has_point;
		vertex_has_normal = vertex_has_normal && rhs.vertex_has_normal;
		vertex_has_tex2d  = vertex_has_tex2d && rhs.vertex_has_tex2d;
		vertex_has_tex3d  = vertex_has_tex3d && rhs.vertex_has_tex3d;
		vertex_has_color  = vertex_has_color && rhs.vertex_has_color;
		face_has_normal   = face_has_normal && rhs.face_has_normal;
		face_has_color    = face_has_color && rhs.face_has_color;
		face_has_material = face_has_material && rhs.face_has_material;
		color_has_alpha   = color_has_alpha && rhs.color_has_alpha;
		stl_binary        = stl_binary && rhs.stl_binary;

		return *this;
	}
};

} // namespace OMC