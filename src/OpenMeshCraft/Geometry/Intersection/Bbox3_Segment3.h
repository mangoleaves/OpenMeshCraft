#pragma once

#include "Bbox3_BoundedLine3.h"

namespace OMC {

/**
 * @brief Check if BoundingBox3 and Segment3 intersect.
 * @tparam Kernel
 */
template <typename Kernel>
class Bbox3_Segment3_Do_Intersect
{
public:
	using K       = Kernel;
	using Segment = typename K::Segment3;
	using Bbox    = typename K::BoundingBox3;

	using Bbox3_BdLine3_DoInter = Bbox3_BoundedLine3_Do_Intersect<Kernel>;

public:
	bool operator()(const Bbox &box, const Segment &segment) const
	{
		return Bbox3_BdLine3_DoInter().do_intersect(
		  box, segment.start(), /*bounded start*/ true, segment.end(),
		  /*bounded end*/ true);
	}
};

} // namespace OMC