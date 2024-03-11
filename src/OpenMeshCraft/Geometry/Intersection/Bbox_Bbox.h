#pragma once

namespace OMC {

/**
 * @brief Check if two BoundingBox intersect.
 * @tparam Kernel.
 */
template <typename Kernel>
class Bbox_Bbox_Do_Intersect
{
public:
	using K     = Kernel;
	using Bbox2 = typename K::BoundingBox2;
	using Bbox3 = typename K::BoundingBox3;

	// Bounding box must contains explicit points, so we directly call points
	// comparator.

	// using LessThan2D = typename K::LessThan2D;
	// using LessThan3D = typename K::LessThan3D;

public:
	bool operator()(const Bbox2 &lbox, const Bbox2 &rbox) const
	{
		return lbox.min_bound() <= rbox.max_bound() &&
		       rbox.min_bound() <= lbox.max_bound();
	}

	bool operator()(const Bbox3 &lbox, const Bbox3 &rbox) const
	{
		return lbox.min_bound() <= rbox.max_bound() &&
		       rbox.min_bound() <= lbox.max_bound();
	}
};

} // namespace OMC