#pragma once

namespace OMC {

/**
 * @brief Check if BoundingBox and Point intersect.
 * @tparam Kernel.
 */
template <typename Kernel>
class Bbox_Point_Do_Intersect
{
public:
	using K = Kernel;

	using GPoint2 = typename K::GPoint2;
	using GPoint3 = typename K::GPoint3;
	using EPoint2 = typename K::EPoint2;
	using EPoint3 = typename K::EPoint3;
	using Bbox2   = typename K::BoundingBox2;
	using Bbox3   = typename K::BoundingBox3;

	using LessThan2D = typename K::LessThan2D;
	using LessThan3D = typename K::LessThan3D;

public:
	template <typename GPT,
	          typename = std::enable_if_t<std::is_same_v<GPT, GPoint2> &&
	                                      !std::is_same_v<GPT, EPoint2>>>
	bool operator()(const Bbox2 &box, const GPT &point) const
	{
		return LessThan2D().on_x(box.min_bound(), point) <= Sign::ZERO &&
		       LessThan2D().on_y(box.min_bound(), point) <= Sign::ZERO &&
		       LessThan2D().on_x(point, box.max_bound()) <= Sign::ZERO &&
		       LessThan2D().on_y(point, box.max_bound()) <= Sign::ZERO;
	}

	template <typename GPT,
	          typename = std::enable_if_t<std::is_same_v<GPT, GPoint3> &&
	                                      !std::is_same_v<GPT, EPoint3>>>
	bool operator()(const Bbox3 &box, const GPT &point) const
	{
		return LessThan3D().on_x(box.min_bound(), point) <= Sign::ZERO &&
		       LessThan3D().on_y(box.min_bound(), point) <= Sign::ZERO &&
		       LessThan3D().on_z(box.min_bound(), point) <= Sign::ZERO &&
		       LessThan3D().on_x(point, box.max_bound()) <= Sign::ZERO &&
		       LessThan3D().on_y(point, box.max_bound()) <= Sign::ZERO &&
		       LessThan3D().on_z(point, box.max_bound()) <= Sign::ZERO;
	}

	bool operator()(const Bbox2 &box, const EPoint2 &point) const
	{
		return (box.min_bound() <= point) && (point <= box.max_bound());
	}

	bool operator()(const Bbox3 &box, const EPoint3 &point) const
	{
		return (box.min_bound() <= point) && (point <= box.max_bound());
	}
};

} // namespace OMC