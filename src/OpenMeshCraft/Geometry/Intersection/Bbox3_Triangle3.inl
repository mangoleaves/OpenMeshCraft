#pragma once

#include "Bbox3_Triangle3.h"
#include "Bbox_Point.h"

namespace OMC {

template <typename Kernel>
bool Bbox3_Triangle3_Do_Intersect<Kernel>::operator()(
  const BboxT &box, const TriangleT &triangle) const
{
	const EPointT &p0 = triangle.v0(), &p1 = triangle.v1(), &p2 = triangle.v2();

	// test three points
	EPointT min_p = p0;
	min_p.minimize(p1);
	min_p.minimize(p2);
	EPointT max_p = p0;
	max_p.maximize(p1);
	max_p.maximize(p2);
	if (!(min_p <= box.max_bound()) || !(max_p >= box.min_bound()))
		return false;

	// test plane
	VecT    tri_n = (p1 - p0).cross(p2 - p0);
	NT      dist  = -1.0 * tri_n.dot(p0.as_vec());
	EPointT P_min, P_max;
	for (int i = 0; i < 3; ++i)
	{
		if (tri_n[i] > NT(0))
		{
			P_min[i] = box.min_coord(i);
			P_max[i] = box.max_coord(i);
		}
		else
		{
			P_min[i] = box.max_coord(i);
			P_max[i] = box.min_coord(i);
		}
	}
	if (tri_n.dot(P_min.as_vec()) > -dist || tri_n.dot(P_max.as_vec()) < -dist)
		return false;

	// test axis edges
	NT box_length[3] = {(box.max_coord(0) - box.min_coord(0)) * 0.5,
	                    (box.max_coord(1) - box.min_coord(1)) * 0.5,
	                    (box.max_coord(2) - box.min_coord(2)) * 0.5};

	const EPointT center =
	  box.min_bound() + (box.max_bound() - box.min_bound()) * 0.5;

	EPointT points[3] = {EPointT(p0 - center), EPointT(p1 - center),
	                     EPointT(p2 - center)};

	const EPointT *u[3];

	// Edge e0
	const VecT e0 = points[1] - points[0];
	u[0]          = &points[0];
	u[1]          = &points[1];
	u[2]          = &points[2];
	if (!TestAxisEdges(u, e0, box_length))
		return false;

	// Edge e1
	const VecT e1 = points[2] - points[1];
	u[0]          = &points[1];
	u[1]          = &points[2];
	u[2]          = &points[0];
	if (!TestAxisEdges(u, e1, box_length))
		return false;

	// Edge e2
	const VecT e2 = points[2] - points[0];
	u[0]          = &points[2];
	u[1]          = &points[0];
	u[2]          = &points[1];
	if (!TestAxisEdges(u, e2, box_length))
		return false;

	return true;
}

#define Absolute(a) ((a) >= NT(0) ? (a) : -(a))
template <typename Kernel>
bool Bbox3_Triangle3_Do_Intersect<Kernel>::TestAxisEdges(
  const EPointT *v[], const VecT &e, const NT *box_length) const
{
	NT fex = Absolute(e[0]);
	NT fey = Absolute(e[1]);
	NT fez = Absolute(e[2]);

	// e and Axis X
	// n = (1,0,0)[cross]e = (0, -e.z, e.y)
	NT p0 = -e[2] * (*v[0])[1] + e[1] * (*v[0])[2];
	NT p2 = -e[2] * (*v[2])[1] + e[1] * (*v[2])[2];
	NT r  = fez * box_length[1] + fey * box_length[2];

	NT min_s = p0 < p2 ? p0 : p2;
	NT max_s = p0 < p2 ? p2 : p0;
	if (max_s < -r || min_s > r)
		return false;

	// e and Axis Y
	// n = (0,1,0)[cross]e = (e.z,0,-e.x)
	p0 = e[2] * (*v[0])[0] - e[0] * (*v[0])[2];
	p2 = e[2] * (*v[2])[0] - e[0] * (*v[2])[2];
	r  = fez * box_length[0] + fex * box_length[2];

	min_s = p0 < p2 ? p0 : p2;
	max_s = p0 < p2 ? p2 : p0;
	if (max_s < -r || min_s > r)
		return false;

	// e and Axis Z
	// n = (0,0,1)[cross]e = (-e.y,e.x,0)
	p0    = -e[1] * (*v[0])[0] + e[0] * (*v[0])[1];
	p2    = -e[1] * (*v[2])[0] + e[0] * (*v[2])[1];
	r     = fey * box_length[0] + fex * box_length[1];
	min_s = p0 < p2 ? p0 : p2;
	max_s = p0 < p2 ? p2 : p0;
	if (max_s < -r || min_s > r)
		return false;

	return true;
}
#undef Absolute

} // namespace OMC