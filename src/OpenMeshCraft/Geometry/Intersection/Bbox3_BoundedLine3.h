#pragma once

#include "Bbox_Point.h"
#include "IntersectionUtils.h"

namespace OMC {

/**
 * @brief Check if BoundingBox3 and a bounded line intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Bbox3_BoundedLine3_Do_Intersect
{
public:
	using K      = Kernel;
	using NT     = typename K::NT;
	using GPoint = typename K::GPoint3;
	using EPoint = typename K::EPoint3;
	using Bbox   = typename K::BoundingBox3;

	using OrientOn2D = typename K::OrientOn2D;
	using LessThan3D = typename K::LessThan3D;

	using Box_Pnt_Do_Inter = Bbox_Point_Do_Intersect<Kernel>;

public:
	/// @brief The bounded line is given by two end points \p p and \p q.
	/// If the bounded flag is set to false, the line is unbounded at the
	/// corresponding end point.
	/// When both sides are bounded, the line is a segment.
	/// When one side is bounded, the line is a ray.
	/// When neither side is bounded, the line is an infinite line.
	bool do_intersect(const Bbox &box, const GPoint &p, bool p_bounded,
	                  const GPoint &q, bool q_bounded) const
	{
		if (Box_Pnt_Do_Inter()(box, p) || Box_Pnt_Do_Inter()(box, q))
			return true;

		Sign ret;

		auto less_x = [](const GPoint &lhs, const GPoint &rhs)
		{ return LessThan3D().on_x(lhs, rhs) == Sign::NEGATIVE; };
		auto less_y = [](const GPoint &lhs, const GPoint &rhs)
		{ return LessThan3D().on_y(lhs, rhs) == Sign::NEGATIVE; };
		auto less_z = [](const GPoint &lhs, const GPoint &rhs)
		{ return LessThan3D().on_z(lhs, rhs) == Sign::NEGATIVE; };
		auto leq_x = [](const GPoint &lhs, const GPoint &rhs)
		{ return LessThan3D().on_x(lhs, rhs) <= Sign::ZERO; };
		auto leq_y = [](const GPoint &lhs, const GPoint &rhs)
		{ return LessThan3D().on_y(lhs, rhs) <= Sign::ZERO; };
		auto leq_z = [](const GPoint &lhs, const GPoint &rhs)
		{ return LessThan3D().on_z(lhs, rhs) <= Sign::ZERO; };

		const GPoint &bmin = box.min_bound();
		const GPoint &bmax = box.max_bound();

		// First, we calculate the interval on the line that is inside between box
		// on each dimension. The interval is represented by two end points on the
		// line. The two end points are further represented by parameters of the
		// line: tmin and tmax.
		// For exactness, tmin and tmax are represented in rational form: tmin =
		// rmin/d, tmax = rmax/d, where d is distance between two end points.
		bool x_degenerate = LessThan3D().on_x(p, q) == Sign::ZERO;
		bool y_degenerate = LessThan3D().on_y(p, q) == Sign::ZERO;
		bool z_degenerate = LessThan3D().on_z(p, q) == Sign::ZERO;
		bool x_reverse, y_reverse, z_reverse;
		NT   rminx, rmaxx, rminy, rmaxy, rminz, rmaxz;

		// -----------------------------------
		// treat x coord
		// -----------------------------------
		if (leq_x(p, q)) // p.x <= q.x
		{
			if (p_bounded && less_x(bmax, p))
				// box.max.x < p.x <= q.x
				return false;
			if (q_bounded && less_x(q, bmin))
				// p.x <= q.x < box.min.x
				return false;
			rmaxx     = NT(bmax.x());
			rminx     = NT(bmin.x());
			// tmax = (rmaxx - p.x) / (q.x - p.x)
			// tmin = (rminx - p.x) / (q.x - p.x)
			x_reverse = false;
		}
		else // q.x < p.x
		{
			if (q_bounded && less_x(bmax, q))
				// bmax.x < q.x < p.x
				return false;
			if (p_bounded && less_x(p, bmin))
				// q.x < p.x < bmin
				return false;
			rmaxx     = NT(bmin.x());
			rminx     = NT(bmax.x());
			// tmax = (rmaxx - p.x) / (q.x - p.x)
			// tmin = (rminx - p.x) / (q.x - p.x)
			// To make denominator > 0, adjust to:
			// tmax = (p.x - rmaxx) / (p.x - q.x)
			// tmin = (p.x - rminx) / (p.x - q.x)
			x_reverse = true;
		}

		if (!p_bounded || !q_bounded)
		{
			// p.x == q.x, the line is vertical to x axis, check the line is
			// in the bounding box.
			if (less_x(bmax, p) || less_x(p, bmin))
				return false;
		}

		// -----------------------------------
		// treat y coord (totally same as x)
		// -----------------------------------
		if (leq_y(p, q)) // p.y <= q.y
		{
			if (p_bounded && less_y(bmax, p))
				// box.max.y < p.y <= q.y
				return false;
			if (q_bounded && less_y(q, bmin))
				// p.y <= q.y < box.min.y
				return false;
			rmaxy     = NT(bmax.y());
			rminy     = NT(bmin.y());
			y_reverse = false;
		}
		else // q.y < p.y
		{
			if (q_bounded && less_y(bmax, q))
				// bmax.y < q.y < p.y
				return false;
			if (p_bounded && less_y(p, bmin))
				// q.y < p.y < bmin
				return false;
			rmaxy     = NT(bmin.y());
			rminy     = NT(bmax.y());
			y_reverse = true;
		}

		if (!p_bounded || !q_bounded)
		{
			// p.y == q.y, the line is vertical to y axis, check the line is
			// in the bounding box.
			if (less_y(bmax, p) || less_y(p, bmin))
				return false;
		}

		// check if intervals on x and y onverlap.
		if (!x_degenerate || !y_degenerate)
		{
			// between x,y
			// (1) tmaxx < tminy ?
			//  => (rmaxx - p.x) / (q.x - p.x) < (rminy - p.y) / (q.y - p.y) ?
			//  => (rmaxx - p.x) * (q.y - p.y) < (rminy - p.y) * (q.x - p.x) ?
			//  => (rmaxx - p.x) * (q.y - p.y) - (rminy - p.y) * (q.x - p.x) < 0 ?
			//     (Reverse sign if x_reverse and y_reverse are not same)
			//  =>        | rmaxx rminy  1 |
			//      det = | q.x    q.y   1 |  < 0 ?
			//            | p.x    p.y   1 |
			//  => OrientOn2D.on_xy (p, q, (rmaxx, rminy, 1)) > 0 ?
			ret = OrientOn2D().on_xy(p, q, EPoint(rmaxx, rminy, 1.));
			if ((x_reverse == y_reverse && ret == Sign::POSITIVE) ||
			    (x_reverse != y_reverse && ret == Sign::NEGATIVE))
				return false;
			// (2) tmaxy < tminx ?
			//  => (rmaxy - p.y) / (q.y - p.y) < (rminx - p.x) / (q.x - p.x) ?
			//  => (rmaxy - p.y) * (q.x - p.x) < (rminx - p.x) * (q.y - p.y) ?
			//  => (rminx - p.x) * (q.y - p.y) - (rmaxy - p.y) * (q.x - p.x) > 0 ?
			//     (Reverse sign if x_reverse and y_reverse are not same)
			//  =>        | rminx rmaxy  1 |
			//      det = | q.x    q.y   1 |  > 0 ?
			//            | p.x    p.y   1 |
			//  => OrientOn2D.on_xy (p, q, (rminx, rmaxy, 1)) < 0 ?
			ret = OrientOn2D().on_xy(p, q, EPoint(rminx, rmaxy, 1.));
			if ((x_reverse == y_reverse && ret == Sign::NEGATIVE) ||
			    (x_reverse != y_reverse && ret == Sign::POSITIVE))
				return false;
		}

		// -----------------------------------
		// treat z coord (totally same as x and y)
		// -----------------------------------
		if (leq_z(p, q)) // p.z <= q.z
		{
			if (p_bounded && less_z(bmax, p))
				// box.max.z <= p.z < q.z
				return false;
			if (q_bounded && less_z(q, bmin))
				// p.z <= q.z < box.min.z
				return false;
			rmaxz     = NT(bmax.z());
			rminz     = NT(bmin.z());
			z_reverse = false;
		}
		else // q.z < p.z
		{
			if (q_bounded && less_z(bmax, q))
				// bmax.z < q.z < p.z
				return false;
			if (p_bounded && less_z(p, bmin))
				// q.z < p.z < bmin
				return false;
			rmaxz     = NT(bmin.z());
			rminz     = NT(bmax.z());
			z_reverse = true;
		}

		if (!p_bounded || !q_bounded)
		{
			// p.z == q.z, the line is vertical to z axis, check the line is
			// in the bounding box.
			if (less_z(bmax, p) || less_z(p, bmin))
				return false;
		}

		// check if intervals on x, y and z overlap each other.

		if (!y_degenerate || !z_degenerate)
		{
			ret = OrientOn2D().on_yz(p, q, EPoint(1., rmaxy, rminz));
			if ((y_reverse == z_reverse && ret == Sign::POSITIVE) ||
			    (y_reverse != z_reverse && ret == Sign::NEGATIVE))
				return false;

			ret = OrientOn2D().on_yz(p, q, EPoint(1., rminy, rmaxz));
			if ((y_reverse == z_reverse && ret == Sign::NEGATIVE) ||
			    (y_reverse != z_reverse && ret == Sign::POSITIVE))
				return false;
		}

		if (!x_degenerate || !z_degenerate)
		{
			ret = OrientOn2D().on_zx(p, q, EPoint(rmaxx, 1., rminz));
			if ((x_reverse == z_reverse && ret == Sign::NEGATIVE) ||
			    (x_reverse != z_reverse && ret == Sign::POSITIVE))
				return false;

			ret = OrientOn2D().on_zx(p, q, EPoint(rminx, 1., rmaxz));
			if ((x_reverse == z_reverse && ret == Sign::POSITIVE) ||
			    (x_reverse != z_reverse && ret == Sign::NEGATIVE))
				return false;
		}

		return true;
	}
};

} // namespace OMC