#pragma once

#include "OpenMeshCraft/Geometry/Primitives/Point2T.h"
#include "OpenMeshCraft/Geometry/Primitives/Point3T.h"
#include "OpenMeshCraft/Geometry/Primitives/Vector2T.h"
#include "OpenMeshCraft/Geometry/Primitives/Vector3T.h"

#include "OpenMeshCraft/Utils/IndexDef.h"

#include <vector>

namespace OMC {

class TriSoupTraits
{
public:
	using VecT      = Vec3T<double>;
	using NormalT   = Vec3T<double>;
	using PointT    = Point3T<double>;
	using Tex3D     = Point3T<double>;
	using Color     = Vec3T<float>;
	using Triangle  = Vec3T<index_t>;
	using Points    = std::vector<PointT>;
	using Normals   = std::vector<NormalT>;
	using Tex3Ds    = std::vector<Tex3D>;
	using Colors    = std::vector<Color>;
	using Triangles = std::vector<Triangle>;
};

} // namespace OMC