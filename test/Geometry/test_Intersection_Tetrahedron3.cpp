#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/Intersection/Tetrahedron3_Point3.h"
#include "OpenMeshCraft/Geometry/Intersection/Tetrahedron3_Segment3.h"
#include "OpenMeshCraft/Geometry/Predicates/IndirectPredicateDetails.h"
#include "OpenMeshCraft/NumberTypes/FPU.h"
#include "OpenMeshCraft/Utils/Macros.h"

#include "gtest/gtest.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/intersection_3.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include <execution>

// #define CHECK_EACH

class test_Intersection_Tetrahedron3 : public testing::Test
{
protected:
	using EIAC = OMC::EIAC;
	using Pnt2 = EIAC::Point2;
	using Pnt3 = EIAC::Point3;

	using CollinearPoints3D               = EIAC::CollinearPoints3D;
	using Tetrahedron3_Point3_DoIntersect = EIAC::Tetrahedron3_Point3_DoIntersect;
	using Tetrahedron3_Segment3_DoIntersect =
	  EIAC::Tetrahedron3_Segment3_DoIntersect;

	using FPnt2 = OMC::Point2T<double>;
	using FPnt3 = OMC::Point3T<double>;

	using EKIC = CGAL::Exact_predicates_inexact_constructions_kernel;
};

inline double randomUnitDouble() { return ((double)rand()) / RAND_MAX; }

TEST_F(test_Intersection_Tetrahedron3, pointTetrahedron)
{
	size_t       num_all_groups    = 100000;
	const double perc_degn         = 0.5;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all tetrahedron and point
	double *atp = new double[(num_random_groups + num_degn_groups) * 15];

	// Create vector of random tetrahedron and point
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 15; i++)
		random_groups[i] = randomUnitDouble();

	// Create vector of random tetrahedron and point on boundary
	double  group[15];
	double *degn_groups = atp + num_random_groups * 15;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		// clang-format off
		// the random tetrahedron
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = randomUnitDouble(); group[7] = randomUnitDouble(); group[8] = randomUnitDouble();
		group[9] = randomUnitDouble(); group[10] = randomUnitDouble(); group[11] = randomUnitDouble();
		// the coplanar point
		const double u = randomUnitDouble();
		const double v = 1 - u;
		group[12] = group[0] + (group[3] - group[0]) * u + (group[6] - group[0]) * v;
		group[13] = group[1] + (group[4] - group[1]) * u + (group[7] - group[1]) * v;
		group[14] = group[2] + (group[5] - group[2]) * u + (group[8] - group[2]) * v;
		// clang-format on
		for (size_t j = 0; j < 15; j++)
			degn_groups[i * 15 + j] = group[j];
	}

	std::vector<Pnt3>          our_points;
	std::vector<EKIC::Point_3> ref_points;
	our_points.reserve(num_all_groups * 5);
	ref_points.reserve(num_all_groups * 5);
	for (size_t i = 0; i < 5 * num_all_groups; ++i)
	{
		our_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
		ref_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
	}

	// Calculate predicates on them
	int  cgal_dummy = 0;
	auto start      = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += EKIC::Do_intersect_3()(
		  EKIC::Tetrahedron_3(ref_points[5 * i], ref_points[5 * i + 1],
		                      ref_points[5 * i + 2], ref_points[5 * i + 3]),
		  ref_points[5 * i + 4]);
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += Tetrahedron3_Point3_DoIntersect()(
		  our_points[5 * i], our_points[5 * i + 1], our_points[5 * i + 2],
		  our_points[5 * i + 3], our_points[5 * i + 4]);
	}

	std::cout << "OpenMeshCraft elapsed time: "
	          << OMC::Logger::elapsed(start).count() << "s\n";
	std::cout << "Dummy sum: " << ours_dummy << "\n";

	EXPECT_EQ(cgal_dummy, ours_dummy);

#ifdef CHECK_EACH
	std::vector indices(num_all_groups, 0);
	std::iota(indices.begin(), indices.end(), 0);
	std::for_each(
	  std::execution::seq, indices.begin(), indices.end(),
	  [&](size_t i)
	  {
		  bool ref_sign = EKIC::Do_intersect_3()(EKIC::Do_intersect_3()(
		    EKIC::Tetrahedron_3(ref_points[5 * i], ref_points[5 * i + 1],
		                        ref_points[5 * i + 2], ref_points[5 * i + 3]),
		    ref_points[5 * i + 4]));
		  bool our_sign = Tetrahedron3_Point3_DoIntersect()(
		    our_points[5 * i], our_points[5 * i + 1], our_points[5 * i + 2],
		    our_points[5 * i + 3], our_points[5 * i + 4]);

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}

TEST_F(test_Intersection_Tetrahedron3, segmentTetrahedron)
{
	size_t       num_all_groups    = 100000;
	const double perc_degn         = 0.5;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all tetrahedron and segment
	double *atp = new double[(num_random_groups + num_degn_groups) * 18];

	// Create vector of random tetrahedron and segment
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 18; i++)
		random_groups[i] = randomUnitDouble();

	double  group[18];
	double *degn_groups = atp + num_random_groups * 18;
	// Create vector of tetrahedron and coplanar segment
	for (size_t i = 0; i < size_t(num_degn_groups * 0.5); i++)
	{
		// clang-format off
		// the random tetrahedron
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = randomUnitDouble(); group[7] = randomUnitDouble(); group[8] = randomUnitDouble();
		group[9] = randomUnitDouble(); group[10] = randomUnitDouble(); group[11] = randomUnitDouble();
		// the coplanar segment
		double u = randomUnitDouble();
		double v = 1 - u;
		group[12] = group[0] + (group[3] - group[0]) * u + (group[6] - group[0]) * v;
		group[13] = group[1] + (group[4] - group[1]) * u + (group[7] - group[1]) * v;
		group[14] = group[2] + (group[5] - group[2]) * u + (group[8] - group[2]) * v;
		u = randomUnitDouble();
		v = 1 - u;
		group[15] = group[0] + (group[3] - group[0]) * u + (group[6] - group[0]) * v;
		group[16] = group[1] + (group[4] - group[1]) * u + (group[7] - group[1]) * v;
		group[17] = group[2] + (group[5] - group[2]) * u + (group[8] - group[2]) * v;
		// clang-format on
		for (size_t j = 0; j < 18; j++)
			degn_groups[i * 18 + j] = group[j];
	}
	// Create vector of tetrahedron and segment originating from one vertex
	for (size_t i = size_t(num_degn_groups * 0.5); i < num_degn_groups; i++)
	{
		// clang-format off
		// the random tetrahedron
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = randomUnitDouble(); group[7] = randomUnitDouble(); group[8] = randomUnitDouble();
		group[9] = randomUnitDouble(); group[10] = randomUnitDouble(); group[11] = randomUnitDouble();
		// the segment originating from the first vertex
		group[12] = group[0];
		group[13] = group[1];
		group[14] = group[2];
		group[15] = randomUnitDouble();
		group[16] = randomUnitDouble();
		group[17] = randomUnitDouble();
		// clang-format on
		for (size_t j = 0; j < 18; j++)
			degn_groups[i * 18 + j] = group[j];
	}

	std::vector<Pnt3>          our_points;
	std::vector<EKIC::Point_3> ref_points;
	our_points.reserve(num_all_groups * 6);
	ref_points.reserve(num_all_groups * 6);
	for (size_t i = 0; i < 6 * num_all_groups; ++i)
	{
		our_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
		ref_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
	}

	// Calculate predicates on them
	int  cgal_dummy = 0;
	auto start      = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += EKIC::Do_intersect_3()(
		  EKIC::Tetrahedron_3(ref_points[6 * i], ref_points[6 * i + 1],
		                      ref_points[6 * i + 2], ref_points[6 * i + 3]),
		  EKIC::Segment_3(ref_points[6 * i + 4], ref_points[6 * i + 5]));
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += Tetrahedron3_Segment3_DoIntersect()(
		  our_points[6 * i], our_points[6 * i + 1], our_points[6 * i + 2],
		  our_points[6 * i + 3], our_points[6 * i + 4], our_points[6 * i + 5]);
	}

	std::cout << "OpenMeshCraft elapsed time: "
	          << OMC::Logger::elapsed(start).count() << "s\n";
	std::cout << "Dummy sum: " << ours_dummy << "\n";

	EXPECT_EQ(cgal_dummy, ours_dummy);

#ifdef CHECK_EACH
	std::vector indices(num_all_groups, 0);
	std::iota(indices.begin(), indices.end(), 0);
	std::for_each(
	  std::execution::seq, indices.begin(), indices.end(),
	  [&](size_t i)
	  {
		  bool ref_sign = EKIC::Do_intersect_3()(
		    EKIC::Tetrahedron_3(ref_points[6 * i], ref_points[6 * i + 1],
		                        ref_points[6 * i + 2], ref_points[6 * i + 3]),
		    EKIC::Segment_3(ref_points[6 * i + 4], ref_points[6 * i + 5]));

		  bool our_sign = Tetrahedron3_Segment3_DoIntersect()(
		    our_points[6 * i], our_points[6 * i + 1], our_points[6 * i + 2],
		    our_points[6 * i + 3], our_points[6 * i + 4], our_points[6 * i + 5]);

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}

TEST_F(test_Intersection_Tetrahedron3, triangleTetrahedron)
{
	size_t       num_all_groups    = 100000;
	const double perc_degn         = 0.5;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all tetrahedron and segment
	double *atp = new double[(num_random_groups + num_degn_groups) * 21];

	// Create vector of random tetrahedron and segment
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 21; i++)
		random_groups[i] = randomUnitDouble();

	double  group[21];
	double *degn_groups = atp + num_random_groups * 21;
	// Create vector of tetrahedron and coplanar segment
	for (size_t i = 0; i < size_t(num_degn_groups * 0.33); i++)
	{
		// clang-format off
		// the random tetrahedron
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = randomUnitDouble(); group[7] = randomUnitDouble(); group[8] = randomUnitDouble();
		group[9] = randomUnitDouble(); group[10] = randomUnitDouble(); group[11] = randomUnitDouble();
		// the coplanar triangle
		double u = randomUnitDouble();
		double v = 1 - u;
		group[12] = group[0] + (group[3] - group[0]) * u + (group[6] - group[0]) * v;
		group[13] = group[1] + (group[4] - group[1]) * u + (group[7] - group[1]) * v;
		group[14] = group[2] + (group[5] - group[2]) * u + (group[8] - group[2]) * v;
		u = randomUnitDouble();
		v = 1 - u;
		group[15] = group[0] + (group[3] - group[0]) * u + (group[6] - group[0]) * v;
		group[16] = group[1] + (group[4] - group[1]) * u + (group[7] - group[1]) * v;
		group[17] = group[2] + (group[5] - group[2]) * u + (group[8] - group[2]) * v;
		u = randomUnitDouble();
		v = 1 - u;
		group[18] = group[0] + (group[3] - group[0]) * u + (group[6] - group[0]) * v;
		group[19] = group[1] + (group[4] - group[1]) * u + (group[7] - group[1]) * v;
		group[20] = group[2] + (group[5] - group[2]) * u + (group[8] - group[2]) * v;
		// clang-format on
		for (size_t j = 0; j < 21; j++)
			degn_groups[i * 21 + j] = group[j];
	}
	// Create vector of tetrahedron and segment originating from one vertex
	for (size_t i = size_t(num_degn_groups * 0.33);
	     i < size_t(num_degn_groups * 0.66); i++)
	{
		// clang-format off
		// the random tetrahedron
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = randomUnitDouble(); group[7] = randomUnitDouble(); group[8] = randomUnitDouble();
		group[9] = randomUnitDouble(); group[10] = randomUnitDouble(); group[11] = randomUnitDouble();
		// the triangle originating from the first vertex
		group[12] = group[0];
		group[13] = group[1];
		group[14] = group[2];
		group[15] = randomUnitDouble();
		group[16] = randomUnitDouble();
		group[17] = randomUnitDouble();
		group[18] = randomUnitDouble();
		group[19] = randomUnitDouble();
		group[20] = randomUnitDouble();
		// clang-format on
		for (size_t j = 0; j < 21; j++)
			degn_groups[i * 21 + j] = group[j];
	}
	// Create vector of tetrahedron and segment originating from one edge
	for (size_t i = size_t(num_degn_groups * 0.66); i < num_degn_groups; i++)
	{
		// clang-format off
		// the random tetrahedron
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = randomUnitDouble(); group[7] = randomUnitDouble(); group[8] = randomUnitDouble();
		group[9] = randomUnitDouble(); group[10] = randomUnitDouble(); group[11] = randomUnitDouble();
		// the triangle originating from the first vertex
		group[12] = group[0];
		group[13] = group[1];
		group[14] = group[2];
		group[15] = group[3];
		group[16] = group[4];
		group[17] = group[5];
		group[18] = randomUnitDouble();
		group[19] = randomUnitDouble();
		group[20] = randomUnitDouble();
		// clang-format on
		for (size_t j = 0; j < 21; j++)
			degn_groups[i * 21 + j] = group[j];
	}

	std::vector<Pnt3>          our_points;
	std::vector<EKIC::Point_3> ref_points;
	our_points.reserve(num_all_groups * 6);
	ref_points.reserve(num_all_groups * 6);
	for (size_t i = 0; i < 6 * num_all_groups; ++i)
	{
		our_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
		ref_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
	}

	// Calculate predicates on them
	int  cgal_dummy = 0;
	auto start      = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += EKIC::Do_intersect_3()(
		  EKIC::Tetrahedron_3(ref_points[6 * i], ref_points[6 * i + 1],
		                      ref_points[6 * i + 2], ref_points[6 * i + 3]),
		  EKIC::Segment_3(ref_points[6 * i + 4], ref_points[6 * i + 5]));
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += Tetrahedron3_Segment3_DoIntersect()(
		  our_points[6 * i], our_points[6 * i + 1], our_points[6 * i + 2],
		  our_points[6 * i + 3], our_points[6 * i + 4], our_points[6 * i + 5]);
	}

	std::cout << "OpenMeshCraft elapsed time: "
	          << OMC::Logger::elapsed(start).count() << "s\n";
	std::cout << "Dummy sum: " << ours_dummy << "\n";

	EXPECT_EQ(cgal_dummy, ours_dummy);

#ifdef CHECK_EACH
	std::vector indices(num_all_groups, 0);
	std::iota(indices.begin(), indices.end(), 0);
	std::for_each(
	  std::execution::seq, indices.begin(), indices.end(),
	  [&](size_t i)
	  {
		  bool ref_sign = EKIC::Do_intersect_3()(
		    EKIC::Tetrahedron_3(ref_points[6 * i], ref_points[6 * i + 1],
		                        ref_points[6 * i + 2], ref_points[6 * i + 3]),
		    EKIC::Segment_3(ref_points[6 * i + 4], ref_points[6 * i + 5]));

		  bool our_sign = Tetrahedron3_Segment3_DoIntersect()(
		    our_points[6 * i], our_points[6 * i + 1], our_points[6 * i + 2],
		    our_points[6 * i + 3], our_points[6 * i + 4], our_points[6 * i + 5]);

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}