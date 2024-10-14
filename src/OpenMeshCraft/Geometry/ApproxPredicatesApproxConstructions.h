#pragma once

// Primitives
#include "Primitives/Vector2T.h"
#include "Primitives/Vector3T.h"
#include "Primitives/VectorXT.h"

#include "Primitives/Points.h"

#include "Primitives/RangeT.h"

#include "Primitives/BoundingBox2T.h"
#include "Primitives/BoundingBox3T.h"

#include "Primitives/Segment2T.h"
#include "Primitives/Sphere2T.h"
#include "Primitives/Triangle2T.h"

#include "Primitives/BoundedLine3T.h"
#include "Primitives/ImplicitPlane3T.h"
#include "Primitives/Line3T.h"
#include "Primitives/NormalCone3T.h"
#include "Primitives/Plane3T.h"
#include "Primitives/Ray3T.h"
#include "Primitives/Segment3T.h"
#include "Primitives/Sphere3T.h"
#include "Primitives/Tetrahedron3T.h"
#include "Primitives/Triangle3T.h"

// Predicates
#include "Predicates/GeneralPredicates.h"

// Intersection
#include "Intersection/DoIntersectK.h"

// Degeneration
#include "Predicates/CheckDegenerate3K.h"

// Constructions
#include "Constructions/Normal3.h"

// Projection
#include "Constructions/FastProjectPoint3K.h"
#include "Constructions/ProjectPoint3K.h"

// Bounding Volume Calculation
#include "Constructions/CalcBoundingBox3K.h"

namespace OMC {

/**
 * @brief This defines types for approximate predicates and approximate
 * construction.
 */
class ApproxPredicatesApproxConstructions
{
public:
	using Kernel = ApproxPredicatesApproxConstructions;

	/*******************************************/
	/********** Kernel properties **************/
	/*******************************************/
	static constexpr bool ArePredicatesExact    = false;
	static constexpr bool AreConstructionsExact = false;
	static constexpr bool ImplicitAndIndirect   = false;

public:
	/*******************************************/
	/************** Numbers ********************/
	/*******************************************/

	/// @name Number types
	/// @{
	using NT = double;
	/// @}

	/*****************************************/
	/************ Vectors ********************/
	/*****************************************/

	/// @name Vector types
	/// @{
	using Vec2 = Vec2T<NT>;
	using Vec3 = Vec3T<NT>;
	/// @}

	/********************************************/
	/************ Primitives ********************/
	/********************************************/

	/// @name Points
	/// @{
	using Point2     = Point2T<NT>;
	using Point3     = Point3T<NT>;
	using GPoint2    = Point2T<NT>;
	using GPoint3    = Point3T<NT>;
	using EPoint2    = Point2T<NT>;
	using EPoint3    = Point3T<NT>;
	// To make kernels based on implicit points are compatibal with kernels based
	// on regular points, we let implicit kernels be aware of implicit points
	// construction but unaware of underlying representation of implicit points.
	using AsGP       = AsGenericPoint_Ex<NT>;
	using AsEP       = AsExplicitPoint_Ex<NT>;
	using ToEP       = ToExplicitPoint_Ex<NT>;
	using CreateSSI2 = CreateImplicitSSI2_Ex<NT>;
	using CreateSSI3 = CreateImplicitSSI3_Ex<NT>;
	using CreateLPI  = CreateImplicitLPI_Ex<NT>;
	using CreateTPI  = CreateImplicitTPI_Ex<NT>;
	/// @}

	static_assert(std::is_trivially_copyable<Point2>::value);
	static_assert(std::is_trivially_copyable<Point3>::value);

	/// @name Primitive types
	/// @{
	using Range = RangeT<NT>;

	using BoundingBox2 = BoundingBox2T<NT, Vec2, EPoint2>;
	using BoundingBox3 = BoundingBox3T<NT, Vec3, EPoint3>;

	using Sphere2   = Sphere2T<NT, EPoint2>;
	using Segment2  = Segment2T<NT, EPoint2>;
	using Triangle2 = Triangle2T<NT, EPoint2>;

	using BoundedLine3   = BoundedLine3T<NT, EPoint3>;
	using ImplicitPlane3 = ImplicitPlane3T<NT, Vec3, EPoint3>;
	using Line3          = Line3T<NT, Vec3, EPoint3>;
	using NormalCone3    = NormalCone3T<NT, Vec3>;
	using Plane3         = Plane3T<NT, Vec3, EPoint3>;
	using Ray3           = Ray3T<NT, Vec3, EPoint3>;
	using Segment3       = Segment3T<NT, EPoint3>;
	using Sphere3        = Sphere3T<NT, EPoint3>;
	using Tetrahedron3   = Tetrahedron3T<NT, EPoint3>;
	using Triangle3      = Triangle3T<NT, EPoint3>;
	/// @}

	/********************************************/
	/************ Predicates ********************/
	/********************************************/

	/// @name Base predicates
	/// @{
	using DotProductSign2D   = DotProductSign2D_GNR<NT>;
	using DotProductSign3D   = DotProductSign3D_GNR<NT>;
	using DotProductSignOn2D = DotProductSignOn2D_GNR<NT>;

	using Orient2D          = Orient2D_GNR<NT>;
	using CollinearPoints2D = CollinearPoints2D_GNR<NT>;
	using LessThan2D        = LessThan2D_GNR<NT>;

	using Orient3D          = Orient3D_GNR<NT>;
	using OrientOn2D        = OrientOn2D_GNR<NT>;
	using CollinearPoints3D = CollinearPoints3D_GNR<NT>;
	using CollinearSort3D   = CollinearSort3D_GNR<NT>;
	using LessThan3D        = LessThan3D_GNR<NT>;
	/// @}

	/********************************************/
	/********* Basic constructions **************/
	/********************************************/
	using ConstructNormal3 = ConstructNormal3K<Kernel>;

	/**********************************************/
	/************ Degeneration ********************/
	/**********************************************/

	/// @name Degeneration
	/// @{
	using CheckDegenerate3 = CheckDegenerate3K<Kernel>;
	/// @}

	/**********************************************/
	/************ Intersection ********************/
	/**********************************************/

	/// @name Do Intersect
	/// @{
	// clang-format off
	using Segment2_Segment2_DoIntersect = Segment2_Segment2_Do_Intersect<Kernel>;
	using Triangle2_Point2_DoIntersect  = Triangle2_Point2_Do_Intersect<Kernel>;

	using Bbox_Bbox_DoIntersect          = Bbox_Bbox_Do_Intersect<Kernel>;
	using Bbox_Point_DoIntersect         = Bbox_Point_Do_Intersect<Kernel>;
	using Bbox_Sphere_DoIntersect        = Bbox_Sphere_Do_Intersect<Kernel>;
	using Bbox3_Segment3_DoIntersect     = Bbox3_Segment3_Do_Intersect<Kernel>;
	using Bbox3_BoundedLine3_DoIntersect = Bbox3_BoundedLine3_Do_Intersect<Kernel>;
	using Bbox3_Line3_DoIntersect        = Bbox3_Line3_Do_Intersect<Kernel>;
	using Bbox3_Ray3_DoIntersect         = Bbox3_Ray3_Do_Intersect<Kernel>;
	using Bbox3_Triangle3_DoIntersect    = Bbox3_Triangle3_Do_Intersect<Kernel>;

	using Segment3_Point3_DoIntersect   = Segment3_Point3_Do_Intersect<Kernel>;
	using Segment3_Segment3_DoIntersect = Segment3_Segment3_Do_Intersect<Kernel>;

	using Triangle3_Point3_DoIntersect    = Triangle3_Point3_Do_Intersect<Kernel>;
	using Triangle3_Segment3_DoIntersect  = Triangle3_Segment3_Do_Intersect<Kernel>;
	using Triangle3_Triangle3_DoIntersect = Triangle3_Triangle3_Do_Intersect<Kernel>;

	using Tetrahedron3_Point3_DoIntersect    = Tetrahedron3_Point3_Do_Intersect<Kernel>;
	using Tetrahedron3_Segment3_DoIntersect  = Tetrahedron3_Segment3_Do_Intersect<Kernel>;
	// clang-format on

	using DoIntersect = DoIntersectK<Kernel>;
	/// @}

	/***********************************************/
	/************ Constructions ********************/
	/***********************************************/

	/************ Projection ********************/

	/// @name Projection
	/// @{
	using ProjectPoint3     = ProjectPoint3K<Kernel>;
	using FastProjectPoint3 = FastProjectPoint3K<Kernel>;
	/// @}

	/********** Bounding Volume *****************/

	/// @name Axis-Aligned Bounding Box
	/// @{
	using CalcBoundingBox3 = CalcBoundingBox3K<Kernel>;
	/// @}
};

using APAC = ApproxPredicatesApproxConstructions;

} // namespace OMC