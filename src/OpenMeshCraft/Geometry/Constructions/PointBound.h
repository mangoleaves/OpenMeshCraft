#pragma once

#include "OpenMeshCraft/Geometry/Primitives/ExplicitPoint2T.h"
#include "OpenMeshCraft/Geometry/Primitives/ExplicitPoint3T.h"
#include "OpenMeshCraft/Geometry/Primitives/GenericPoint2T.h"
#include "OpenMeshCraft/Geometry/Primitives/GenericPoint3T.h"

#include <type_traits>

namespace OMC {

/**
 * @brief Calculate lower/upper bound for various types of point
 * and output the result to an explicit point type.
 * @tparam GPointT Generic point type (e.g., implicit point type)
 * @tparam EPointT Explicit point type (storing coordinates explicitly)
 */
template <typename GPointT, typename EPointT>
class PointBound
{
public:
	template <typename GPT,
	          typename = std::enable_if_t<std::is_same_v<GPT, GPointT> &&
	                                      std::is_same_v<GPT, EPointT>>>
	const EPointT &lower_bound(const GPointT &gp)
	{
		return gp;
	}

	template <typename GPT,
	          typename = std::enable_if_t<std::is_same_v<GPT, GPointT> &&
	                                      std::is_same_v<GPT, EPointT>>>
	const EPointT &upper_bound(const GPointT &gp)
	{
		return gp;
	}
};

template <typename IT, typename ET>
class PointBound<GenericPoint2T<IT, ET>, ExplicitPoint2T<IT, ET>>
{
public:
	using EP = ExplicitPoint2T<IT, ET>;
	using GP = GenericPoint2T<IT, ET>;

	EP lower_bound(const GP &gp)
	{
		if (gp.is_Explicit())
		{
			return gp.to_Explicit();
		}
		else
		{
			IT lx, ly, d;
			gp.getIntervalLambda(lx, ly, d);
			IT x = lx / d, y = ly / d;
			return EP(x.inf(), y.inf());
		}
	}

	EP upper_bound(const GP &gp)
	{
		if (gp.is_Explicit())
		{
			return gp.to_Explicit();
		}
		else
		{
			IT lx, ly, d;
			gp.getIntervalLambda(lx, ly, d);
			IT x = lx / d, y = ly / d;
			return EP(x.sup(), y.sup());
		}
	}
};

template <typename IT, typename ET>
class PointBound<GenericPoint3T<IT, ET>, ExplicitPoint3T<IT, ET>>
{
public:
	using EP = ExplicitPoint3T<IT, ET>;
	using GP = GenericPoint3T<IT, ET>;

	EP lower_bound(const GP &gp)
	{
		if (gp.is_Explicit())
		{
			return gp.to_Explicit();
		}
		else
		{
#if defined(INDIRECT_PREDICATES)
			IT lx, ly, lz, d;
			gp.getIntervalLambda(lx, ly, lz, d);
			IT x = lx / d, y = ly / d, z = lz / d;
#else
			IT lx, ly, lz, d, bx, by, bz;
			gp.getIntervalLambda(lx, ly, lz, d, bx, by, bz);
			IT x = lx / d + bx, y = ly / d + by, z = lz / d + bz;
#endif
			return EP(x.inf(), y.inf(), z.inf());
		}
	}

	EP upper_bound(const GP &gp)
	{
		if (gp.is_Explicit())
		{
			return gp.to_Explicit();
		}
		else
		{
#if defined(INDIRECT_PREDICATES)
			IT lx, ly, lz, d;
			gp.getIntervalLambda(lx, ly, lz, d);
			IT x = lx / d, y = ly / d, z = lz / d;
#else
			IT lx, ly, lz, d, bx, by, bz;
			gp.getIntervalLambda(lx, ly, lz, d, bx, by, bz);
			IT x = lx / d + bx, y = ly / d + by, z = lz / d + bz;
#endif
			return EP(x.sup(), y.sup(), z.sup());
		}
	}
};

} // namespace OMC