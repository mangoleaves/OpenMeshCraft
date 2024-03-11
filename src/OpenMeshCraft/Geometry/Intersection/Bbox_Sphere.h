#pragma once

#include <cstddef>
#include <type_traits>

namespace OMC {

/**
 * @brief Check if BoundingBox and Sphere intersect.
 * @tparam Kernel.
 * @fixme This predicate is inexact in all kernels with inexact construction.
 * Because current implementation containts inexact construction.
 */
template <typename Kernel>
class Bbox_Sphere_Do_Intersect
{
public:
	using K = Kernel;

	using NT      = typename K::NT;
	using GPoint2 = typename K::GPoint2;
	using GPoint3 = typename K::GPoint3;
	using EPoint2 = typename K::EPoint2;
	using EPoint3 = typename K::EPoint3;
	using ToEP    = typename K::ToEP;
	using Bbox2   = typename K::BoundingBox2;
	using Bbox3   = typename K::BoundingBox3;
	using Sphere2 = typename K::Sphere2;
	using Sphere3 = typename K::Sphere3;

public:
	bool operator()(const Bbox2 &box, const Sphere2 &sphere) const
	{
		const EPoint2 &center   = sphere.center();
		NT             d        = NT(0.0);
		NT             distance = NT(0.0);
		// for x
		if (center.x() < box.min_bound().x())
		{
			d = box.min_bound().x() - center.x();
			distance += d * d;
		}
		else if (center.x() > box.max_bound().x())
		{
			d = center.x() - box.max_bound().x();
			distance += d * d;
		}
		if (distance > sphere.squared_radius())
			return false;
		// for y
		if (center.y() < box.min_bound().y())
		{
			d = box.min_bound().y() - center.y();
			distance += d * d;
		}
		else if (center.y() > box.max_bound().y())
		{
			d = center.y() - box.max_bound().y();
			distance += d * d;
		}
		if (distance > sphere.squared_radius())
			return false;
		return true;
	}

	bool operator()(const Bbox3 &box, const Sphere3 &sphere) const
	{
		const EPoint3 &center   = sphere.center();
		NT             d        = NT(0.0);
		NT             distance = NT(0.0);
		// for x
		if (center.x() < box.min_bound().x())
		{
			d = box.min_bound().x() - center.x();
			distance += d * d;
		}
		else if (center.x() > box.max_bound().x())
		{
			d = center.x() - box.max_bound().x();
			distance += d * d;
		}
		if (distance > sphere.squared_radius())
			return false;
		// for y
		if (center.y() < box.min_bound().y())
		{
			d = box.min_bound().y() - center.y();
			distance += d * d;
		}
		else if (center.y() > box.max_bound().y())
		{
			d = center.y() - box.max_bound().y();
			distance += d * d;
		}
		if (distance > sphere.squared_radius())
			return false;
		// for z
		if (center.z() < box.min_bound().z())
		{
			d = box.min_bound().z() - center.z();
			distance += d * d;
		}
		else if (center.z() > box.max_bound().z())
		{
			d = center.z() - box.max_bound().z();
			distance += d * d;
		}
		if (distance > sphere.squared_radius())
			return false;
		return true;
	}
};

} // namespace OMC