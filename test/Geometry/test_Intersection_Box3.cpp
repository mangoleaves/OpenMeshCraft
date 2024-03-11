#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/intersection_3.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "test_utils.h"

#include <type_traits>

class test_Intersection_Box3 : public testing::Test
{
protected:
	using Our_K  = OMC::APAC;
	using CGAL_K = CGAL::Simple_cartesian<double>;
	// using CGAL_K = CGAL::Epick;

	using Our_Pnt3  = Our_K::Point3;
	using Our_Bbox3 = Our_K::BoundingBox3;
	using Our_Seg3  = Our_K::Segment3;
	using Our_Line3 = Our_K::Line3;
	using Our_Ray3  = Our_K::Ray3;
	using Do_Inter  = Our_K::DoIntersect;
};

TEST_F(test_Intersection_Box3, Do_Intersect_Point3)
{
	Do_Inter  do_inter;
	Our_Pnt3  our_p(1., 1., 1.);
	Our_Bbox3 our_box(Our_K::Point3(0., 0., 0.), Our_K::Point3(2., 2., 2.));

	ASSERT_EQ(do_inter(our_box, our_p), true);

	our_p.x() = 2.0;
	ASSERT_EQ(do_inter(our_box, our_p), true);
	our_p.y() = 2.0;
	ASSERT_EQ(do_inter(our_box, our_p), true);
	our_p.z() = 2.0;
	ASSERT_EQ(do_inter(our_box, our_p), true);
	our_p.x() = 0.0;
	ASSERT_EQ(do_inter(our_box, our_p), true);
	our_p.y() = 0.0;
	ASSERT_EQ(do_inter(our_box, our_p), true);
	our_p.z() = 0.0;
	ASSERT_EQ(do_inter(our_box, our_p), true);
	our_p.x() = std::nextafter(2.0, 3.0);
	ASSERT_EQ(do_inter(our_box, our_p), false);
	our_p.x() = std::nextafter(0.0, -1.0);
	ASSERT_EQ(do_inter(our_box, our_p), false);
}

TEST_F(test_Intersection_Box3, Do_Intersect_Segment3)
{
	Our_Seg3  our_seg;
	Our_Bbox3 our_box(Our_K::Point3(0.5, 0.5, 0.5), Our_K::Point3(1.5, 1.5, 1.5));
	our_seg.start() = Our_K::Point3(0.11, 0.09, 0.1);

	CGAL::Bbox_3    cgal_box = CGAL_K::Construct_bbox_3()(CGAL_K::Segment_3(
    CGAL_K::Point_3(0.5, 0.5, 0.5), CGAL_K::Point_3(1.5, 1.5, 1.5)));
	CGAL_K::Point_3 cgal_seg_start(0.11, 0.09, 0.1);

	size_t total     = 1000000;
	size_t same_cnt  = 0;
	size_t inter_cnt = 0;
	for (size_t i = 0; i < total; i++)
	{
		double x = (double)rand() / RAND_MAX * 2.;
		double y = (double)rand() / RAND_MAX * 2.;
		double z = (double)rand() / RAND_MAX * 2.;

		CGAL_K::Segment_3 cgal_seg(cgal_seg_start, CGAL_K::Point_3(x, y, z));
		our_seg.end() = Our_K::Point3(x, y, z);

		bool our_res  = Do_Inter()(our_box, our_seg);
		bool cgal_res = CGAL::do_intersect<CGAL_K>(cgal_box, cgal_seg);
		inter_cnt += cgal_res;
		same_cnt += our_res == cgal_res;
	}

	EXPECT_EQ(same_cnt, total);
	std::cout << std::format(
	  "Test {} cases. Same cases: {}. Intersect cases: {}\n", total, same_cnt,
	  inter_cnt);
}

TEST_F(test_Intersection_Box3, Do_Intersect_Line3)
{
	Our_Line3 our_line;
	Our_Bbox3 our_box(Our_K::Point3(0.5, 0.5, 0.5), Our_K::Point3(1.5, 1.5, 1.5));
	our_line.start() = Our_K::Point3(0.11, 0.09, 0.1);

	CGAL::Bbox_3    cgal_box = CGAL_K::Construct_bbox_3()(CGAL_K::Segment_3(
    CGAL_K::Point_3(0.5, 0.5, 0.5), CGAL_K::Point_3(1.5, 1.5, 1.5)));
	CGAL_K::Point_3 cgal_line_start(0.11, 0.09, 0.1);

	size_t total     = 1000000;
	size_t same_cnt  = 0;
	size_t inter_cnt = 0;
	for (size_t j = 0; j < total; j++)
	{
		double x = (double)rand() / RAND_MAX * 2.;
		double y = (double)rand() / RAND_MAX * 2.;
		double z = (double)rand() / RAND_MAX * 2.;

		CGAL_K::Line_3 cgal_line(cgal_line_start,
		                         CGAL_K::Point_3(x, y, z) - cgal_line_start);
		our_line.direction() = Our_K::Point3(x, y, z) - our_line.start();

		bool our_res  = Do_Inter()(our_box, our_line);
		bool cgal_res = CGAL::do_intersect<CGAL_K>(cgal_box, cgal_line);
		inter_cnt += cgal_res;
		same_cnt += our_res == cgal_res;
	}

	EXPECT_EQ(same_cnt, total);
	std::cout << std::format(
	  "Test {} cases. Same cases: {}. Intersect cases: {}\n", total, same_cnt,
	  inter_cnt);
}

TEST_F(test_Intersection_Box3, Do_Intersect_Ray3)
{
	Our_Ray3  our_ray;
	Our_Bbox3 our_box(Our_K::Point3(0.5, 0.5, 0.5), Our_K::Point3(1.5, 1.5, 1.5));
	our_ray.start() = Our_K::Point3(0.11, 0.09, 0.1);

	CGAL::Bbox_3    cgal_box = CGAL_K::Construct_bbox_3()(CGAL_K::Segment_3(
    CGAL_K::Point_3(0.5, 0.5, 0.5), CGAL_K::Point_3(1.5, 1.5, 1.5)));
	CGAL_K::Point_3 cgal_ray_start(0.11, 0.09, 0.1);

	size_t total     = 1000000;
	size_t same_cnt  = 0;
	size_t inter_cnt = 0;
	for (size_t k = 0; k < total; k++)
	{
		double x = (double)rand() / RAND_MAX * 2.;
		double y = (double)rand() / RAND_MAX * 2.;
		double z = (double)rand() / RAND_MAX * 2.;

		CGAL_K::Ray_3 cgal_ray(cgal_ray_start,
		                       CGAL_K::Point_3(x, y, z) - cgal_ray_start);
		our_ray.direction() = Our_K::Point3(x, y, z) - our_ray.start();

		bool our_res  = Do_Inter()(our_box, our_ray);
		bool cgal_res = CGAL::do_intersect<CGAL_K>(cgal_box, cgal_ray);
		inter_cnt += cgal_res;
		same_cnt += our_res == cgal_res;
	}

	EXPECT_EQ(same_cnt, total);
	std::cout << std::format(
	  "Test {} cases. Same cases: {}. Intersect cases: {}\n", total, same_cnt,
	  inter_cnt);
}