#pragma once

#include "Bbox_Bbox.h"
#include "Bbox_Point.h"
#include "Bbox_Sphere.h"

#include "Segment2_Segment2.h"
#include "Triangle2_Point2.h"

#include "Bbox3_Line3.h"
#include "Bbox3_Ray3.h"
#include "Bbox3_Segment3.h"
#include "Bbox3_Triangle3.h"

#include "Segment3_Point3.h"
#include "Segment3_Segment3.h"

#include "Triangle3_Point3.h"
#include "Triangle3_Segment3.h"
#include "Triangle3_Triangle3.h"

#include "Tetrahedron3_Point3.h"
#include "Tetrahedron3_Segment3.h"
#include "Tetrahedron3_Triangle3.h"

namespace OMC {

template <typename Kernel>
class DoIntersectK
{
public:
	using K = Kernel;

	using GPoint2      = typename K::GPoint2;
	using EPoint2      = typename K::EPoint2;
	using BoundingBox2 = typename K::BoundingBox2;
	using Segment2     = typename K::Segment2;
	using Sphere2      = typename K::Sphere2;
	using Triangle2    = typename K::Triangle2;

	using GPoint3      = typename K::GPoint3;
	using EPoint3      = typename K::EPoint3;
	using BoundedLine3 = typename K::BoundedLine3;
	using BoundingBox3 = typename K::BoundingBox3;
	using Line3        = typename K::Line3;
	using Ray3         = typename K::Ray3;
	using Segment3     = typename K::Segment3;
	using Sphere3      = typename K::Sphere3;
	using Tetrahedron3 = typename K::Tetrahedron3;
	using Triangle3    = typename K::Triangle3;

public:
	// clang-format off
	template <typename GPT, typename = std::enable_if_t<std::is_same_v<GPT, GPoint2> && !std::is_same_v<GPT, EPoint2>>>
	bool operator()(const BoundingBox2 &box,  const GPT          &point)  const { return Bbox_Point_Do_Intersect<K>()(box, point); }
	bool operator()(const BoundingBox2 &box,  const EPoint2      &point)  const { return Bbox_Point_Do_Intersect<K>()(box, point); }
	bool operator()(const BoundingBox2 &box1, const BoundingBox2 &box2)   const { return Bbox_Bbox_Do_Intersect<K>()(box1, box2); }
	bool operator()(const BoundingBox2 &box,  const Sphere2      &sphere) const { return Bbox_Sphere_Do_Intersect<K>()(box, sphere); }
	bool operator()(const Segment2     &seg1, const Segment2     &seg2)   const { return Segment2_Segment2_Do_Intersect<K>()(seg1, seg2); }
	bool operator()(const Triangle2    &tri,  const GPoint2      &point)  const { return Triangle2_Point2_Do_Intersect<K>()(tri, point); }

	template <typename GPT, typename = std::enable_if_t<std::is_same_v<GPT, GPoint3> && !std::is_same_v<GPT, EPoint3>>>
	bool operator()(const BoundingBox3 &box,  const GPT          &point)    const { return Bbox_Point_Do_Intersect<K>()(box, point); }
	bool operator()(const BoundingBox3 &box,  const EPoint3      &point)    const { return Bbox_Point_Do_Intersect<K>()(box, point); }
	bool operator()(const BoundingBox3 &box1, const BoundingBox3 &box2)     const { return Bbox_Bbox_Do_Intersect<K>()(box1, box2); }
	bool operator()(const BoundingBox3 &box,  const Sphere3      &sphere)   const { return Bbox_Sphere_Do_Intersect<K>()(box, sphere); }
	bool operator()(const BoundingBox3 &box,  const Segment3     &segment)  const { return Bbox3_Segment3_Do_Intersect<K>()(box, segment); }
	bool operator()(const BoundingBox3 &box,  const BoundedLine3 &line)     const { return Bbox3_BoundedLine3_Do_Intersect<K>()( box, line.start(), line.start_bounded(), line.end(), line.end_bounded()); }
	bool operator()(const BoundingBox3 &box,  const Line3        &line)     const { return Bbox3_Line3_Do_Intersect<K>()(box, line); }
	bool operator()(const BoundingBox3 &box,  const Ray3         &ray)      const { return Bbox3_Ray3_Do_Intersect<K>()(box, ray); }
	bool operator()(const BoundingBox3 &box,  const Triangle3    &triangle) const { return Bbox3_Triangle3_Do_Intersect<K>()(box, triangle); }

	template <typename GPT, typename = std::enable_if_t<std::is_same_v<GPT, GPoint3> && !std::is_same_v<GPT, EPoint3>>>
	bool operator()(const Segment3 &seg,  const GPT      &point) const { return Segment3_Point3_Do_Intersect<K>()(seg, point); }
	bool operator()(const Segment3 &seg,  const EPoint3  &point) const { return Segment3_Point3_Do_Intersect<K>()(seg, point); }
	bool operator()(const Segment3 &seg1, const Segment3 &seg2)  const { return Segment3_Segment3_Do_Intersect<K>()(seg1, seg2); }

	template <typename GPT, typename = std::enable_if_t<std::is_same_v<GPT, GPoint3> && !std::is_same_v<GPT, EPoint3>>>
	bool operator()(const Triangle3 &tri,  const GPT       &point) const { return Triangle3_Point3_Do_Intersect<K>()(tri, point); }
	bool operator()(const Triangle3 &tri,  const EPoint3   &point) const { return Triangle3_Point3_Do_Intersect<K>()(tri, point); }
	bool operator()(const Triangle3 &tri,  const Segment3  &seg)   const { return Triangle3_Segment3_Do_Intersect<K>()(tri, seg); }
	bool operator()(const Triangle3 &tri1, const Triangle3 &tri2)  const { return Triangle3_Triangle3_Do_Intersect<K>()(tri1, tri2); }

	template <typename GPT, typename = std::enable_if_t<std::is_same_v<GPT, GPoint3> && !std::is_same_v<GPT, EPoint3>>>
	bool operator()(const Tetrahedron3 &tet,  const GPT       &point) const { return Tetrahedron3_Point3_Do_Intersect<K>()(tet, point); }
	bool operator()(const Tetrahedron3 &tet,  const EPoint3   &point) const { return Tetrahedron3_Point3_Do_Intersect<K>()(tet, point); }
	bool operator()(const Tetrahedron3 &tet,  const Segment3  &seg)   const { return Tetrahedron3_Segment3_Do_Intersect<K>()(tet, seg); }
	bool operator()(const Tetrahedron3 &tet,  const Triangle3 &tri)   const { return Tetrahedron3_Triangle3_Do_Intersect<K>()(tet, tri); }
	// clang-format on
};

} // namespace OMC