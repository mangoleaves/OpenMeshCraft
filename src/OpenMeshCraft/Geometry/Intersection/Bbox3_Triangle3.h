#pragma once

namespace OMC {

/**
 * @brief Check if BoundingBox3 and Triangle3 intersect.
 * @tparam Kernel
 * @todo Implement an standard box-triangle intersection test.
 * @fixme This predicate is inexact in all kernels with inexact construction.
 * Because current implementation containts inexact construction.
 */
template <typename Kernel>
class Bbox3_Triangle3_Do_Intersect
{
public:
	using K  = Kernel;
	using NT = typename K::NT;

	using EPointT   = typename K::EPoint3;
	using VecT      = typename K::Vec3;
	using BboxT     = typename K::BoundingBox3;
	using TriangleT = typename K::Triangle3;

public:
	bool operator()(const BboxT &box, const TriangleT &triangle) const;

	bool TestAxisEdges(const EPointT *v[], const VecT &e,
	                   const NT *box_length) const;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "Bbox3_Triangle3.inl"
#endif