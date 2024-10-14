#pragma once

// Primitives
#include "Primitives/Vector2T.h"
#include "Primitives/Vector3T.h"
#include "Primitives/VectorXT.h"

#include "Primitives/ImplicitPoints.h"
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
#include "Primitives/Plane3T.h"
#include "Primitives/Ray3T.h"
#include "Primitives/Segment3T.h"
#include "Primitives/Sphere3T.h"
#include "Primitives/Tetrahedron3T.h"
#include "Primitives/Triangle3T.h"

// Predicates
#include "Predicates/IndirectPredicates.h"

// Degeneration
#include "Predicates/CheckDegenerate3K.h"

// Intersection
#include "Intersection/DoIntersectK.h"

// Constructions
#include "Constructions/Normal3.h"

// Projection
#include "Constructions/FastProjectPoint3K.h"
#include "Constructions/ProjectPoint3K.h"

// Bounding Volume Calculation
#include "Constructions/CalcBoundingBox3K.h"

namespace OMC {

class ExactIndirectPredicatesApproxConstructions
{
public:
	using Kernel = ExactIndirectPredicatesApproxConstructions;

	/*******************************************/
	/********** Kernel properties **************/
	/*******************************************/
	static constexpr bool ArePredicatesExact    = true;
	static constexpr bool AreConstructionsExact = false;
	static constexpr bool ImplicitAndIndirect   = true;

public:
	/*******************************************/
	/************** Numbers ********************/
	/*******************************************/

	/// @name Number types
	/// @{
	using FT = double;
	using NT = FT;
	using IT = IntervalNumber</*Protected*/ std::true_type>;
	// using ET = LazyNumber</*HasApprox*/ std::false_type, /*Protected*/
	// std::true_type, /*ThreadSafe*/ std::true_type>;
	using ET = BoostRational;
	// using ET = bigfloat;
	/// @}

	/*****************************************/
	/************ Vectors ********************/
	/*****************************************/

	/// @name Vector types
	/// @{
	using Vec2 = Vec2d;
	using Vec3 = Vec3d;
	/// @}

	/********************************************/
	/************ Primitives ********************/
	/********************************************/

	/// @name Points
	/// @{
	using Point2 = ExplicitPoint2T<IT, ET>;
	using Point3 = ExplicitPoint3T<IT, ET>;

	using GPoint2      = GenericPoint2T<IT, ET>;
	// Point2 and SSI inherit from GP2
	using EPoint2      = ExplicitPoint2T<IT, ET>;
	using IPoint2T_SSI = ImplicitPoint2T_SSI<IT, ET>;

	using GPoint3      = GenericPoint3T<IT, ET>;
	// Point3 and SSI/LPI/TPI inherit from GP3
	using EPoint3      = ExplicitPoint3T<IT, ET>;
	using IPoint3T_SSI = ImplicitPoint3T_SSI<IT, ET>;
	using IPoint3T_LPI = ImplicitPoint3T_LPI<IT, ET>;
	using IPoint3T_TPI = ImplicitPoint3T_TPI<IT, ET>;

	static_assert(std::is_trivially_copyable<GPoint2>::value);
	static_assert(std::is_trivially_copyable<GPoint3>::value);
	static_assert(std::is_trivially_copyable<EPoint2>::value);
	static_assert(std::is_trivially_copyable<EPoint3>::value);

	// To make kernels based on implicit points are compatibal with kernels based
	// on regular points, we let implicit kernels be aware of implicit points
	// construction but unaware of underlying representation of implicit points.
	using AsGP       = AsGenericPoint_Im<IT, ET>;
	using AsEP       = AsExplicitPoint_Im<IT, ET>;
	using ToEP       = ToExplicitPoint_Im<IT, ET>;
	using CreateSSI2 = CreateImplicitSSI2_Im<IT, ET>;
	using CreateSSI3 = CreateImplicitSSI3_Im<IT, ET>;
	using CreateLPI  = CreateImplicitLPI_Im<IT, ET>;
	using CreateTPI  = CreateImplicitTPI_Im<IT, ET>;
	/// @}

	/// @name Primitive types
	/// @{
	using Range = RangeT<NT>;

	using BoundingBox2 = BoundingBox2T<FT, Vec2, EPoint2>;
	using BoundingBox3 = BoundingBox3T<FT, Vec3, EPoint3>;

	using Sphere2   = Sphere2T<FT, EPoint2>;
	using Segment2  = Segment2T<FT, EPoint2>;
	using Triangle2 = Triangle2T<FT, EPoint2>;

	using BoundedLine3   = BoundedLine3T<FT, EPoint3>;
	using ImplicitPlane3 = ImplicitPlane3T<FT, Vec3, EPoint3>;
	using Line3          = Line3T<FT, Vec3, EPoint3>;
	using Plane3         = Plane3T<FT, Vec3, EPoint3>;
	using Ray3           = Ray3T<FT, Vec3, EPoint3>;
	using Segment3       = Segment3T<FT, EPoint3>;
	using Sphere3        = Sphere3T<FT, EPoint3>;
	using Tetrahedron3   = Tetrahedron3T<FT, EPoint3>;
	using Triangle3      = Triangle3T<FT, EPoint3>;
	/// @}

	/********************************************/
	/************ Predicates ********************/
	/********************************************/

	/// @name Base predicates
	/// @{
	using DotProductSign2D   = DotProductSign2D_Indirect<FT, IT, ET>;
	using DotProductSign3D   = DotProductSign3D_Indirect<FT, IT, ET>;
	using DotProductSignOn2D = DotProductSignOn2D_Indirect<FT, IT, ET>;

	using Orient2D          = Orient2D_Indirect<FT, IT, ET>;
	using LessThan2D        = LessThan2D_Indirect<FT, IT, ET>;
	using InCircle          = InCircle_Indirect<FT, IT, ET>;
	using CollinearPoints2D = CollinearPoints2D_Indirect<FT, IT, ET>;

	using LongestAxis        = LongestAxis_Indirect<FT, IT, ET>;
	using MaxCompInTriNormal = MaxComponentInTriangleNormal<FT, IT, ET>;

	using Orient3D          = Orient3D_Indirect<FT, IT, ET>;
	using OrientOn2D        = OrientOn2D_Indirect<FT, IT, ET>;
	using LessThan3D        = LessThan3D_Indirect<FT, IT, ET>;
	using CollinearPoints3D = CollinearPoints3D_Indirect<FT, IT, ET>;
	using CollinearSort3D   = CollinearSort3D_Indirect<FT, IT, ET>;
	using InSphere          = InSphere_Indirect<FT, IT, ET>;
	/// @}

	/********************************************/
	/********* Basic constructions **************/
	/********************************************/
	using ConstructNormal3 = ConstructNormal3K<Kernel>;

	/************ Degeneration ********************/

	/// @name Degeneration
	/// @{
	using CheckDegenerate3 = CheckDegenerate3K<Kernel>;
	/// @}

	/************ Intersection ********************/

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

using EIAC = ExactIndirectPredicatesApproxConstructions;

} // namespace OMC