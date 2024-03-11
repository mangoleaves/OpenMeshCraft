#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"

#include "AreVecsEqual.h"

#include "test_utils.h"

#include <type_traits>

using namespace OMC;

TEST(CheckDegenerate3, Simple)
{
	APAC::CheckDegenerate3 CheckDegenerate3;
	APAC::CheckDegenerate3::DgnType degeneration;
	APAC::Point3                 p0{0., 0., 0.}, p1(1., 0., 0.), p2(0.5, 0.5, 0.);

	APAC::Segment3 segment{p0, p0};
	degeneration = CheckDegenerate3(segment);
	ASSERT_EQ(std::holds_alternative<APAC::Point3>(degeneration), true);
	segment.end() = p1;
	degeneration = CheckDegenerate3(segment);
	ASSERT_EQ(std::holds_alternative<APAC::CheckDegenerate3::NoDgn>(degeneration), true);

	APAC::Triangle3 triangle{p0, p0, p0};
	degeneration = CheckDegenerate3(triangle);
	ASSERT_EQ(std::holds_alternative<APAC::Point3>(degeneration), true);
	triangle.v0() = p1;
	degeneration = CheckDegenerate3(triangle);
	ASSERT_EQ(std::holds_alternative<APAC::Segment3>(degeneration), true);
	triangle.v1() = p2;
	degeneration = CheckDegenerate3(triangle);
	ASSERT_EQ(std::holds_alternative<APAC::CheckDegenerate3::NoDgn>(degeneration), true);
}