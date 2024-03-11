#pragma once

namespace OMC {

/**
 * @brief Check if BoundingBox3 and Line3 intersect.
 */
template <typename Kernel>
class Bbox3_Line3_Do_Intersect
{
public:
	using K      = Kernel;
	using NT     = typename K::NT;
	using Vec    = typename K::Vec3;
	using GPoint = typename K::GPoint3;
	using EPoint = typename K::EPoint3;
	using Line   = typename K::Line3;
	using Bbox   = typename K::BoundingBox3;

	using DotProductSignOn2D = typename K::DotProductSignOn2D;

	using Box_Pnt_Do_Inter = Bbox_Point_Do_Intersect<Kernel>;

public:
	bool operator()(const Bbox &box, const Line &line) const
	{
		const EPoint &start = line.start();
		const Vec    &dir   = line.direction();

		if (Box_Pnt_Do_Inter()(box, start))
			return true;

		DotProductSignOn2D dot_sign;
		Sign               ret;

		const EPoint &bmin = box.min_bound();
		const EPoint &bmax = box.max_bound();

		// First, we calculate the interval on the line that is inside between box
		// on each dimension. The interval is represented by two end points on the
		// line. The two end points are further represented by parameters of the
		// line: tmin and tmax.
		// For exactness, tmin and tmax are represented in rational form: tmin =
		// rmin/d, tmax = rmax/d, where d is distance between two end points.
		bool x_degenerate = dir.x() == NT(0);
		bool y_degenerate = dir.y() == NT(0);
		bool z_degenerate = dir.z() == NT(0);
		bool x_reverse, y_reverse, z_reverse;
		NT   rminx, rmaxx, rminy, rmaxy, rminz, rmaxz;

		// -----------------------------------
		// treat x coord
		// -----------------------------------
		if (dir.x() >= NT(0.))
		{
			rmaxx     = NT(bmax.x());
			rminx     = NT(bmin.x());
			// tmax = (rmaxx - start.x) / dir.x
			// tmin = (rminx - start.x) / dir.x
			x_reverse = false;
		}
		else // dir.x < 0
		{
			rmaxx     = NT(bmin.x());
			rminx     = NT(bmax.x());
			// tmax = (rmaxx - start.x) / dir.x
			// tmin = (rminx - start.x) / dir.x
			// To make denominator > 0, adjust to:
			// tmax = (start.x - rmaxx) / (-dir.x)
			// tmin = (start.x - rminx) / (-dir.x)
			x_reverse = true;
		}

		// -----------------------------------
		// treat y coord (totally same as x)
		// -----------------------------------
		if (dir.y() >= NT(0.))
		{
			rmaxy     = NT(bmax.y());
			rminy     = NT(bmin.y());
			// tmax = (rmaxy - start.y) / dir.y
			// tmin = (rminy - start.y) / dir.y
			y_reverse = false;
		}
		else // dir.y < 0
		{
			rmaxy     = NT(bmin.y());
			rminy     = NT(bmax.y());
			// tmax = (rmaxy - start.y) / dir.y
			// tmin = (rminy - start.y) / dir.y
			// To make denominator > 0, adjust to:
			// tmax = (start.y - rmaxy) / (-dir.y)
			// tmin = (start.y - rminy) / (-dir.y)
			y_reverse = true;
		}

		// check if intervals on x, y and z overlap each other.

		if (!x_degenerate || !y_degenerate)
		{
			// between x,y
			// (1) tmaxx < tminy ?
			//  => (rmaxx - start.x) / dir.x < (rminy - start.y) / dir.y ?
			//  => (rmaxx - start.x) * dir.y < (rminy - start.y) * dir.x ?
			//  => (rmaxx - start.x) * dir.y - (rminy - start.y) * dir.x < 0 ?
			//     (Reverse sign if x_reverse and y_reverse are not same)
			//  => DotProductSignOn2D.on_xy
			//     ((rmaxx, rminy, 0.), (dir.y, dir.x, 0.), start, (0.,0.,0.)) < 0 ?
			ret = dot_sign.on_xy(EPoint(rmaxx, rminy, NT(0.)),
			                     EPoint(dir.y(), -dir.x(), NT(0.)), start,
			                     EPoint(NT(0.), NT(0.), NT(0.)));
			if ((x_reverse == y_reverse && ret == Sign::NEGATIVE) ||
			    (x_reverse != y_reverse && ret == Sign::POSITIVE))
				return false;
			// (2) tmaxy < tminx ?
			//  => (rmaxy - start.y) / dir.y < (rminx - start.x) / dir.x ?
			//  => (rmaxy - start.y) * dir.x < (rminx - start.x) * dir.y ?
			//  => (rminx - start.x) * dir.y - (rmaxy - start.y) * dir.x > 0 ?
			//     (Reverse sign if x_reverse and y_reverse are not same)
			//  => DotProductSignOn2D.on_xy
			//     ((rminx, rmaxy, 0.), (dir.y, dir.x, 0.), start, (0.,0.,0.)) > 0 ?
			ret = dot_sign.on_xy(EPoint(rminx, rmaxy, NT(0.)),
			                     EPoint(dir.y(), -dir.x(), NT(0.)), start,
			                     EPoint(NT(0.), NT(0.), NT(0.)));
			if ((x_reverse == y_reverse && ret == Sign::POSITIVE) ||
			    (x_reverse != y_reverse && ret == Sign::NEGATIVE))
				return false;
		}

		// -----------------------------------
		// treat z coord (totally same as x and y)
		// -----------------------------------
		if (dir.z() >= NT(0.))
		{
			rmaxz     = NT(bmax.z());
			rminz     = NT(bmin.z());
			// tmax = (rmaxz - start.z) / dir.z
			// tmin = (rminz - start.z) / dir.z
			z_reverse = false;
		}
		else // dir.z < 0
		{
			rmaxz     = NT(bmin.z());
			rminz     = NT(bmax.z());
			// tmax = (rmaxz - start.z) / dir.z
			// tmin = (rminz - start.z) / dir.z
			// To make denominator > 0, adjust to:
			// tmax = (start.z - rmaxz) / (-dir.z)
			// tmin = (start.z - rminz) / (-dir.z)
			z_reverse = true;
		}

		if (!y_degenerate || !z_degenerate)
		{
			ret = dot_sign.on_yz(EPoint(NT(0.), rmaxy, rminz),
			                     EPoint(NT(0.), dir.z(), -dir.y()), start,
			                     EPoint(NT(0.), NT(0.), NT(0.)));
			if ((y_reverse == z_reverse && ret == Sign::NEGATIVE) ||
			    (y_reverse != z_reverse && ret == Sign::POSITIVE))
				return false;

			ret = dot_sign.on_yz(EPoint(NT(0.), rminy, rmaxz),
			                     EPoint(NT(0.), dir.z(), -dir.y()), start,
			                     EPoint(NT(0.), NT(0.), NT(0.)));
			if ((y_reverse == z_reverse && ret == Sign::POSITIVE) ||
			    (y_reverse != z_reverse && ret == Sign::NEGATIVE))
				return false;
		}

		if (!z_degenerate || !x_degenerate)
		{
			ret = dot_sign.on_zx(EPoint(rmaxx, NT(0.), rminz),
			                     EPoint(dir.z(), NT(0.), -dir.x()), start,
			                     EPoint(NT(0.), NT(0.), NT(0.)));
			if ((z_reverse == x_reverse && ret == Sign::NEGATIVE) ||
			    (z_reverse != x_reverse && ret == Sign::POSITIVE))
				return false;

			ret = dot_sign.on_zx(EPoint(rminx, NT(0.), rmaxz),
			                     EPoint(dir.z(), NT(0.), -dir.x()), start,
			                     EPoint(NT(0.), NT(0.), NT(0.)));
			if ((z_reverse == x_reverse && ret == Sign::POSITIVE) ||
			    (z_reverse != x_reverse && ret == Sign::NEGATIVE))
				return false;
		}

		return true;
	}
};

} // namespace OMC