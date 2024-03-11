#pragma once
#include "OpenMeshCraft/Geometry/Primitives/Vector2T.h"
#include "OpenMeshCraft/Geometry/Primitives/Vector3T.h"
#include "OpenMeshCraft/Geometry/Primitives/VectorXT.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "Eigen/Dense"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "test_utils.h"

template <typename T, typename S, size_t dim>
class AreVecsEqual
{
public:
	::testing::AssertionResult operator()(const char *t_name, const char *s_name,
	                                      const T &t, const S &s)
	{
		for (size_t i = 0; i < dim; i++)
		{
			if (t[i] != s[i])
				return ::testing::AssertionFailure() << std::format(
				         "{}-th elements are not equal: {} in {} != {} in {}.", i, t[i],
				         t_name, s[i], s_name);
		}
		return ::testing::AssertionSuccess();
	}
};

using Vec3d  = OMC::Vec3d;
using tVec3d = OMC::VecXT<3, double>;
using eVec3d = Eigen::Vector3d;

using Vec2d  = OMC::Vec2d;
using tVec2d = OMC::VecXT<2, double>;
using eVec2d = Eigen::Vector2d;

using Vec2i  = OMC::Vec2i;
using tVec2i = OMC::VecXT<2, int>;
using eVec2i = Eigen::Vector2i;

using Vec3i  = OMC::Vec3i;
using tVec3i = OMC::VecXT<3, int>;
using eVec3i = Eigen::Vector3i;