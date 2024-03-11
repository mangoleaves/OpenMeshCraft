#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/intersection_3.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "test_utils.h"

#include <type_traits>

class test_Intersection_Ray3 : public testing::Test
{
protected:
	using Our_K  = OMC::APAC;
	using CGAL_K = CGAL::Epick;

	using Our_Pnt3  = Our_K::Point3;
	using Our_Bbox3 = Our_K::BoundingBox3;
	using Our_Seg3  = Our_K::Segment3;
	using Our_Line3 = Our_K::Line3;
	using Our_Ray3  = Our_K::Ray3;
	using Do_Inter  = Our_K::DoIntersect;
};

TEST_F(test_Intersection_Ray3, Do_Intersect_Triangle3)
{
#if 0
	CGAL_K::Point_3    cgal_ray_start(0.11, 0.09, 0.1);
	CGAL_K::Triangle_3 cgal_tri(CGAL_K::Point_3(1, 0, 0),
	                            CGAL_K::Point_3(0, 1, 0),
	                            CGAL_K::Point_3(0, 0, 1));

	size_t total     = 100000;
	size_t same_cnt  = 0;
	size_t inter_cnt = 0;
	for (size_t i = 0; i < total; i++)
	{
		double x = (double)rand() / RAND_MAX * 2.;
		double y = (double)rand() / RAND_MAX * 2.;
		double z = (double)rand() / RAND_MAX * 2.;

		CGAL_K::Ray_3 cgal_ray(cgal_ray_start, CGAL_K::Point_3(x, y, z));

		bool cgal_res = CGAL::do_intersect<CGAL_K>(cgal_tri, cgal_ray);
		inter_cnt += cgal_res;
	}

	EXPECT_EQ(same_cnt, total);
	std::cout << std::format(
	  "Test {} cases. Same cases: {}. Intersect cases: {}\n", total, same_cnt,
	  inter_cnt);
#endif
}
