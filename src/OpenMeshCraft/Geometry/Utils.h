#pragma once

namespace OMC {

/// @brief An enumeration represents three different orthogonal planes in 3D
/// space: YZ, ZX, and XY. Each plane is assigned a unique integer value.
enum OrPlane
{
	YZ = 0,   // normal is X(0)
	ZX = 1,   // normal is Y(1)
	XY = 2    // normal is Z(2)
};

inline OrPlane int_to_OrPlane(const int &norm)
{
	return static_cast<OrPlane>(norm);
}

inline int OrPlane_to_int(const OrPlane &p) { return static_cast<int>(p); }

} // namespace OMC