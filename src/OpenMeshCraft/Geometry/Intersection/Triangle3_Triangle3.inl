#pragma once

#include "Triangle3_Triangle3.h"

#include <bitset>

namespace OMC {

template <typename Kernel>
bool Triangle3_Triangle3_Do_Intersect<Kernel>::operator()(
  const TriangleT &tri0, const TriangleT &tri1) const
{
	return intersection_type(tri0, tri1) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
bool Triangle3_Triangle3_Do_Intersect<Kernel>::operator()(
  const GPointT &t00, const GPointT &t01, const GPointT &t02,
  const GPointT &t10, const GPointT &t11, const GPointT &t12) const
{
	return intersection_type(t00, t01, t02, t10, t11, t12) >=
	       SimplexIntersectionType::SIMPLICIAL_COMPLEX;
}

template <typename Kernel>
SimplexIntersectionType
Triangle3_Triangle3_Do_Intersect<Kernel>::intersection_type(
  const TriangleT &tri0, const TriangleT &tri1) const
{
	return intersection_type(tri0.v0(), tri0.v1(), tri0.v2(), tri1.v0(),
	                         tri1.v1(), tri1.v2());
}

template <typename Kernel>
SimplexIntersectionType
Triangle3_Triangle3_Do_Intersect<Kernel>::intersection_type(
  const GPointT &t00, const GPointT &t01, const GPointT &t02,
  const GPointT &t10, const GPointT &t11, const GPointT &t12) const
{
	// triangle_is_degenerate_3d
	OMC_EXPENSIVE_ASSERT(CollinearPoints3D().misaligned(t00, t01, t02) &&
	                       CollinearPoints3D().misaligned(t10, t11, t12),
	                     "degenerate triangle");

	// binary flags to mark coincident vertices in t0 and t1
	std::bitset<3> t0_shared = {0b000};
	std::bitset<3> t1_shared = {0b000};

	// find vert correspondences
	// clang-format off
	if (LessThan3D().coincident(t00, t10)) { t0_shared[0] = true; t1_shared[0] = true; }
	if (LessThan3D().coincident(t00, t11)) { t0_shared[0] = true; t1_shared[1] = true; }
	if (LessThan3D().coincident(t00, t12)) { t0_shared[0] = true; t1_shared[2] = true; }
	if (LessThan3D().coincident(t01, t10)) { t0_shared[1] = true; t1_shared[0] = true; }
	if (LessThan3D().coincident(t01, t11)) { t0_shared[1] = true; t1_shared[1] = true; }
	if (LessThan3D().coincident(t01, t12)) { t0_shared[1] = true; t1_shared[2] = true; }
	if (LessThan3D().coincident(t02, t10)) { t0_shared[2] = true; t1_shared[0] = true; }
	if (LessThan3D().coincident(t02, t11)) { t0_shared[2] = true; t1_shared[1] = true; }
	if (LessThan3D().coincident(t02, t12)) { t0_shared[2] = true; t1_shared[2] = true; }
	// clang-format on

	// count number of coincident vertices in t0 and t1
	size_t t0_count = t0_shared.count();

	// either t0 and t1 are coincident
	if (t0_count == 3)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	// t0 and t1 share an edge. Let e be the shared edge and { opp0, opp1 } be the
	// two vertices opposite to e in t0 and t1, respectively. If opp0 and opp1 lie
	// at the same side of e, the two triangles overlap. Otherwise they are
	// edge-adjacent and form a valid simplicial complex
	if (t0_count == 2)
	{
		size_t e[2];      // indices of the shared vertices (in t0)
		size_t count = 0; // index for e (to fill it)
		size_t opp0  = 0; // index of the vertex opposite to e in t0
		size_t opp1  = 0; // index of the vertex opposite to e in t1
		for (size_t i = 0; i < 3; ++i)
		{
			if (!t0_shared[i])
				opp0 = i;
			else
				e[count++] = i;
			if (!t1_shared[i])
				opp1 = i;
		}

		const GPointT *t0[3] = {&t00, &t01, &t02};
		const GPointT *t1[3] = {&t10, &t11, &t12};
		Sign           opp0_wrt_e, opp1_wrt_e;

		// if they are not coplanar, then they form a valid complex
		if (Orient3D()(t00, t01, t02, *t1[opp1]) != Sign::ZERO)
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		opp0_wrt_e = OrientOn2D().on_xy(*t0[e[0]], *t0[e[1]], *t0[opp0]);
		opp1_wrt_e = OrientOn2D().on_xy(*t0[e[0]], *t0[e[1]], *t1[opp1]);
		if ((opp0_wrt_e == Sign::POSITIVE && opp1_wrt_e == Sign::NEGATIVE) ||
		    (opp0_wrt_e == Sign::NEGATIVE && opp1_wrt_e == Sign::POSITIVE))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		opp0_wrt_e = OrientOn2D().on_yz(*t0[e[0]], *t0[e[1]], *t0[opp0]);
		opp1_wrt_e = OrientOn2D().on_yz(*t0[e[0]], *t0[e[1]], *t1[opp1]);
		if ((opp0_wrt_e == Sign::POSITIVE && opp1_wrt_e == Sign::NEGATIVE) ||
		    (opp0_wrt_e == Sign::NEGATIVE && opp1_wrt_e == Sign::POSITIVE))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		opp0_wrt_e = OrientOn2D().on_zx(*t0[e[0]], *t0[e[1]], *t0[opp0]);
		opp1_wrt_e = OrientOn2D().on_zx(*t0[e[0]], *t0[e[1]], *t1[opp1]);
		if ((opp0_wrt_e == Sign::POSITIVE && opp1_wrt_e == Sign::NEGATIVE) ||
		    (opp0_wrt_e == Sign::NEGATIVE && opp1_wrt_e == Sign::POSITIVE))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		return SimplexIntersectionType::INTERSECT;
	}

	// t0 and t1 share a vertex. Let v be the shared vertex and { opp0 , opp1 } be
	// the two edges opposite to v in t0 and t1, respectively. If opp0 intersects
	// t1, or opp1 interects t0, the two triangles overlap. Otherwise they are
	// vertex-adjacent and form a valid simplicial complex
	if (t0_count == 1)
	{
		size_t v0 = 0; // index of the shared vertex in t0
		size_t v1 = 0; // index of the shared vertex in t1
		for (size_t i = 0; i < 3; ++i)
		{
			if (t0_shared[i])
				v0 = i;
			if (t1_shared[i])
				v1 = i;
		}

		const GPointT *t0[3]   = {&t00, &t01, &t02};
		const GPointT *t1[3]   = {&t10, &t11, &t12};
		const GPointT *opp0[2] = {t0[(v0 + 1) % 3], t0[(v0 + 2) % 3]};
		const GPointT *opp1[2] = {t1[(v1 + 1) % 3], t1[(v1 + 2) % 3]};

		// clang-format off
		if (Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, *opp0[0], *opp0[1]) >= SimplexIntersectionType::INTERSECT ||
		    Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, *opp1[0], *opp1[1]) >= SimplexIntersectionType::INTERSECT)
		{
			return SimplexIntersectionType::INTERSECT;
		}
		// clang-format on
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;
	}

	// t0 and t1 do not share sub-simplices. They can be fully disjoint,
	// intersecting at a single point, or overlapping

	// clang-format off
	if (Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, t00, t01) >= SimplexIntersectionType::INTERSECT ||
	    Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, t01, t02) >= SimplexIntersectionType::INTERSECT ||
	    Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, t02, t00) >= SimplexIntersectionType::INTERSECT ||
	    Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, t10, t11) >= SimplexIntersectionType::INTERSECT ||
	    Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, t11, t12) >= SimplexIntersectionType::INTERSECT ||
	    Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, t12, t10) >= SimplexIntersectionType::INTERSECT)
	{
		return SimplexIntersectionType::INTERSECT;
	}
	// clang-format on

	return SimplexIntersectionType::DO_NOT_INTERSECT;
}

template <typename Kernel>
SimplexIntersectionType
Triangle3_Triangle3_Do_Intersect<Kernel>::intersection_type(
  const NT *t00, const NT *t01, const NT *t02, const NT *t10, const NT *t11,
  const NT *t12, const NT *t0_min, const NT *t0_perm, const NT *t1_min,
  const NT *t1_perm) const
{
	// triangle_is_degenerate_3d
	OMC_EXPENSIVE_ASSERT(CollinearPoints3D().misaligned(t00, t01, t02) &&
	                       CollinearPoints3D().misaligned(t10, t11, t12),
	                     "degenerate triangle");
	OMC_INTER_PROFILE_INC_TOTAL(IntersectionNames::T3T3);

	// binary flags to mark coincident vertices in t0 and t1
	std::bitset<3> t0_shared = {0b000};
	std::bitset<3> t1_shared = {0b000};

	// find vert correspondences
	// clang-format off
	if (vec_equals_3d(t00, t10)) { t0_shared[0] = true; t1_shared[0] = true; }
	if (vec_equals_3d(t00, t11)) { t0_shared[0] = true; t1_shared[1] = true; }
	if (vec_equals_3d(t00, t12)) { t0_shared[0] = true; t1_shared[2] = true; }
	if (vec_equals_3d(t01, t10)) { t0_shared[1] = true; t1_shared[0] = true; }
	if (vec_equals_3d(t01, t11)) { t0_shared[1] = true; t1_shared[1] = true; }
	if (vec_equals_3d(t01, t12)) { t0_shared[1] = true; t1_shared[2] = true; }
	if (vec_equals_3d(t02, t10)) { t0_shared[2] = true; t1_shared[0] = true; }
	if (vec_equals_3d(t02, t11)) { t0_shared[2] = true; t1_shared[1] = true; }
	if (vec_equals_3d(t02, t12)) { t0_shared[2] = true; t1_shared[2] = true; }
	// clang-format on

	// count number of coincident vertices in t0 and t1
	size_t t0_count = t0_shared.count();

	// either t0 and t1 are coincident
	if (t0_count == 3)
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

	// t0 and t1 share an edge. Let e be the shared edge and { opp0, opp1 } be the
	// two vertices opposite to e in t0 and t1, respectively. If opp0 and opp1 lie
	// at the same side of e, the two triangles overlap. Otherwise they are
	// edge-adjacent and form a valid simplicial complex
	if (t0_count == 2)
	{
		size_t e[2];      // indices of the shared vertices (in t0)
		size_t count = 0; // index for e (to fill it)
		size_t opp0  = 0; // index of the vertex opposite to e in t0
		size_t opp1  = 0; // index of the vertex opposite to e in t1
		for (size_t i = 0; i < 3; ++i)
		{
			if (!t0_shared[i])
				opp0 = i;
			else
				e[count++] = i;
			if (!t1_shared[i])
				opp1 = i;
		}

		const NT *t0[3] = {t00, t01, t02};
		const NT *t1[3] = {t10, t11, t12};

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 0);

		// if they are not coplanar, then they form a valid complex
		if (t0_min != nullptr &&
		    Orient3D().with_cached_minors(t00, t01, t02, t1[opp1], t0_min,
		                                  t0_perm) != Sign::ZERO)
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;
		else if (Orient3D()(t00, t01, t02, t1[opp1]) != Sign::ZERO)
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 1);

		Sign opp0_wrt_e = OrientOn2D().on_xy(t0[e[0]], t0[e[1]], t0[opp0]);
		Sign opp1_wrt_e = OrientOn2D().on_xy(t0[e[0]], t0[e[1]], t1[opp1]);
		if ((opp0_wrt_e == Sign::POSITIVE && opp1_wrt_e == Sign::NEGATIVE) ||
		    (opp0_wrt_e == Sign::NEGATIVE && opp1_wrt_e == Sign::POSITIVE))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 2);

		opp0_wrt_e = OrientOn2D().on_yz(t0[e[0]], t0[e[1]], t0[opp0]);
		opp1_wrt_e = OrientOn2D().on_yz(t0[e[0]], t0[e[1]], t1[opp1]);
		if ((opp0_wrt_e == Sign::POSITIVE && opp1_wrt_e == Sign::NEGATIVE) ||
		    (opp0_wrt_e == Sign::NEGATIVE && opp1_wrt_e == Sign::POSITIVE))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 3);

		opp0_wrt_e = OrientOn2D().on_zx(t0[e[0]], t0[e[1]], t0[opp0]);
		opp1_wrt_e = OrientOn2D().on_zx(t0[e[0]], t0[e[1]], t1[opp1]);
		if ((opp0_wrt_e == Sign::POSITIVE && opp1_wrt_e == Sign::NEGATIVE) ||
		    (opp0_wrt_e == Sign::NEGATIVE && opp1_wrt_e == Sign::POSITIVE))
			return SimplexIntersectionType::SIMPLICIAL_COMPLEX;

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 4);

		return SimplexIntersectionType::INTERSECT;
	}

	// t0 and t1 share a vertex. Let v be the shared vertex and { opp0 , opp1 } be
	// the two edges opposite to v in t0 and t1, respectively. If opp0 intersects
	// t1, or opp1 interects t0, the two triangles overlap. Otherwise they are
	// vertex-adjacent and form a valid simplicial complex
	if (t0_count == 1)
	{
		size_t v0 = 0; // index of the shared vertex in t0
		size_t v1 = 0; // index of the shared vertex in t1
		for (size_t i = 0; i < 3; ++i)
		{
			if (t0_shared[i])
				v0 = i;
			if (t1_shared[i])
				v1 = i;
		}

		const NT *t0[3] = {t00, t01, t02};
		const NT *t1[3] = {t10, t11, t12};

		const NT *opp0[2] = {t0[(v0 + 1) % 3], t0[(v0 + 2) % 3]};
		const NT *opp1[2] = {t1[(v1 + 1) % 3], t1[(v1 + 2) % 3]};

		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 5);

		int n0_max = -1;
		int n1_max = -1;
		// clang-format off
		if (Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, opp0[0], opp0[1], n1_max, t1_min, t1_perm) >= SimplexIntersectionType::INTERSECT ||
		    Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, opp1[0], opp1[1], n0_max, t0_min, t0_perm) >= SimplexIntersectionType::INTERSECT)
		{
			OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 6);
			return SimplexIntersectionType::INTERSECT;
		}
		// clang-format on
		OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 7);
		return SimplexIntersectionType::SIMPLICIAL_COMPLEX;
	}

	// t0 and t1 do not share sub-simplices. They can be fully disjoint,
	// intersecting at a single point, or overlapping

	int n0_max = -1;
	int n1_max = -1;
	// clang-format off
	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 8);
	if (Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, t00, t01, n1_max, t1_min, t1_perm) >= SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 9);
	if (Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, t01, t02, n1_max, t1_min, t1_perm) >= SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 10);
	if (Triangle3_Segment3_DoInter().intersection_type(t10, t11, t12, t02, t00, n1_max, t1_min, t1_perm) >= SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 11);
	if (Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, t10, t11, n0_max, t0_min, t0_perm) >= SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 12);
	if (Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, t11, t12, n0_max, t0_min, t0_perm) >= SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 13);
	if (Triangle3_Segment3_DoInter().intersection_type(t00, t01, t02, t12, t10, n0_max, t0_min, t0_perm) >= SimplexIntersectionType::INTERSECT)
		return SimplexIntersectionType::INTERSECT;
	// clang-format on

	OMC_INTER_PROFILE_INC_REACH(IntersectionNames::T3T3, 14);
	return SimplexIntersectionType::DO_NOT_INTERSECT;
}

} // namespace OMC