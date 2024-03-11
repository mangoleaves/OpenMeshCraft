#pragma once

#include "AdapOcTree.h"

#include <algorithm>

namespace OMC {

template <typename Traits>
void AdapOcTree<Traits>::calc_box_for_children(
  NodeRef nd, OrPointCRef c, std::array<Bbox, Degree> &child_boxes)
{
	OrPointCRef minb = nd.box().min_bound();
	OrPointCRef maxb = nd.box().max_bound();

	// clang-format off
	// child 0(000, -x-y-z)
	child_boxes[0] = Bbox(minb, c);
	// child 1(001, +x-y-z)
	child_boxes[1] = Bbox(OrPoint(c.x(), minb.y(), minb.z()), OrPoint(maxb.x(), c.y(), c.z()));
	// child 2(010, -x+y-z)
	child_boxes[2] = Bbox(OrPoint(minb.x(), c.y(), minb.z()), OrPoint(c.x(), maxb.y(), c.z()));
	// child 3(011, +x+y-z)
	child_boxes[3] = Bbox(OrPoint(c.x(), c.y(), minb.z()), OrPoint(maxb.x(), maxb.y(), c.z()));
	// child 4(100, -x-y+z)
	child_boxes[4] = Bbox(OrPoint(minb.x(), minb.y(), c.z()), OrPoint(c.x(), c.y(), maxb.z()));
	// child 5(101, +x-y+z)
	child_boxes[5] = Bbox(OrPoint(c.x(), minb.y(), c.z()), OrPoint(maxb.x(), c.y(), maxb.z()));
	// child 6(110, -x+y+z)
	child_boxes[6] = Bbox(OrPoint(minb.x(), c.y(), c.z()), OrPoint(c.x(), maxb.y(), maxb.z()));
	// child 7(111, +x+y+z)
	child_boxes[7] = Bbox(c, maxb);
	// clang-format on
}

} // namespace OMC