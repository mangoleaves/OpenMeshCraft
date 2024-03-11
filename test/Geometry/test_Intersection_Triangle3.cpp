#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/Intersection/Triangle3_Point3.h"
#include "OpenMeshCraft/Geometry/Intersection/Triangle3_Segment3.h"
#include "OpenMeshCraft/Geometry/Intersection/Triangle3_Triangle3.h"
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

class test_Intersection_Triangle3 : public testing::Test
{
protected:
	using EIAC = OMC::EIAC;
	using Pnt2 = EIAC::Point2;
	using Pnt3 = EIAC::Point3;

	using CollinearPoints3D               = EIAC::CollinearPoints3D;
	using Triangle3_Point3_DoIntersect    = EIAC::Triangle3_Point3_DoIntersect;
	using Triangle3_Segment3_DoIntersect  = EIAC::Triangle3_Segment3_DoIntersect;
	using Triangle3_Triangle3_DoIntersect = EIAC::Triangle3_Triangle3_DoIntersect;

	using FPnt2 = OMC::Point2T<double>;
	using FPnt3 = OMC::Point3T<double>;

	using EKIC = CGAL::Exact_predicates_inexact_constructions_kernel;
};

inline double randomUnitDouble() { return ((double)rand()) / RAND_MAX; }

TEST_F(test_Intersection_Triangle3, pointTriangle)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.1;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 12];

	// Create vector of random 3D point tets
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 12; i++)
		random_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned tets
	double  group[12];
	double *degn_groups = atp + num_random_groups * 12;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		// clang-format off
		group[0] = randomUnitDouble();
		group[1] = randomUnitDouble();
		group[2] = randomUnitDouble();
		group[3] = randomUnitDouble();
		group[4] = randomUnitDouble();
		group[5] = randomUnitDouble();
		group[6] = randomUnitDouble();
		group[7] = randomUnitDouble();
		group[8] = randomUnitDouble();
		const double u = randomUnitDouble();
		const double v = 1 - u;
		group[9]  = group[0] + (group[3] - group[0]) * u + (group[6] - group[0]) *v;
		group[10] = group[1] + (group[4] - group[1]) * u + (group[7] - group[1]) *v;
		group[11] = group[2] + (group[5] - group[2]) * u + (group[8] - group[2]) *v;
		// clang-format on
		for (size_t j = 0; j < 12; j++)
			degn_groups[i * 12 + j] = group[j];
	}

	std::vector<Pnt3>          our_points;
	std::vector<EKIC::Point_3> ref_points;
	our_points.reserve(num_all_groups * 4);
	ref_points.reserve(num_all_groups * 4);
	for (size_t i = 0; i < 4 * num_all_groups; ++i)
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
		  EKIC::Triangle_3(ref_points[4 * i], ref_points[4 * i + 1],
		                   ref_points[4 * i + 2]),
		  ref_points[4 * i + 3]);
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += Triangle3_Point3_DoIntersect()(
		  our_points[4 * i], our_points[4 * i + 1], our_points[4 * i + 2],
		  our_points[4 * i + 3]);
	}

	std::cout << "OpenMeshCraft elapsed time: "
	          << OMC::Logger::elapsed(start).count() << "s\n";
	std::cout << "Dummy sum: " << ours_dummy << "\n";

	EXPECT_EQ(cgal_dummy, ours_dummy);

#ifdef CHECK_EACH
	std::vector indices(num_all_groups, 0);
	std::iota(indices.begin(), indices.end(), 0);
	std::for_each(std::execution::seq, indices.begin(), indices.end(),
	              [&](size_t i)
	              {
		              bool ref_sign = EKIC::Do_intersect_3()(
		                EKIC::Triangle_3(ref_points[4 * i], ref_points[4 * i + 1],
		                                 ref_points[4 * i + 2]),
		                ref_points[4 * i + 3]);
		              bool our_sign = Triangle3_Point3_DoIntersect()(
		                our_points[4 * i], our_points[4 * i + 1],
		                our_points[4 * i + 2], our_points[4 * i + 3]);

		              EXPECT_EQ(ref_sign, our_sign);
	              });
#endif
	delete[] atp;
}

TEST_F(test_Intersection_Triangle3, segmentTriangle)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.1;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 15];

	// Create vector of random 3D point tets
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 15; i++)
		random_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned tets
	double  group[15];
	double *degn_groups = atp + num_random_groups * 15;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		// clang-format off
		group[0] = randomUnitDouble();
		group[1] = randomUnitDouble();
		group[2] = randomUnitDouble();
		group[3] = randomUnitDouble();
		group[4] = randomUnitDouble();
		group[5] = randomUnitDouble();
		group[6] = randomUnitDouble();
		group[7] = randomUnitDouble();
		group[8] = randomUnitDouble();

		double u = randomUnitDouble();
		double v = 1 - u;
		group[9]  = group[0] + (group[3] - group[0]) * u + (group[3] - group[0]) * v;
		group[10] = group[1] + (group[4] - group[1]) * u + (group[4] - group[1]) * v;
		group[11] = group[2] + (group[5] - group[2]) * u + (group[5] - group[2]) * v;

		u = randomUnitDouble();
		v = 1 - u;
		group[12]  = group[0] + (group[3] - group[0]) * u + (group[3] - group[0]) * v;
		group[13] = group[1] + (group[4] - group[1]) * u + (group[4] - group[1]) * v;
		group[14] = group[2] + (group[5] - group[2]) * u + (group[5] - group[2]) * v;
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
		  EKIC::Triangle_3(ref_points[5 * i], ref_points[5 * i + 1],
		                   ref_points[5 * i + 2]),
		  EKIC::Segment_3(ref_points[5 * i + 3], ref_points[5 * i + 4]));
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += Triangle3_Segment3_DoIntersect()(
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
		  bool ref_sign = EKIC::Do_intersect_3()(
		    EKIC::Triangle_3(ref_points[5 * i], ref_points[5 * i + 1],
		                     ref_points[5 * i + 2]),
		    EKIC::Segment_3(ref_points[5 * i + 3], ref_points[5 * i + 4]));

		  bool our_sign = Triangle3_Segment3_DoIntersect()(
		    our_points[5 * i], our_points[5 * i + 1], our_points[5 * i + 2],
		    our_points[5 * i + 3], our_points[5 * i + 4]);

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}

TEST_F(test_Intersection_Triangle3, triangleTriangle)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.1;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 18];

	// Create vector of random 3D point tets
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 18; i++)
		random_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned tets
	double  group[18];
	double *degn_groups = atp + num_random_groups * 18;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		// clang-format off
		group[0] = randomUnitDouble();
		group[1] = randomUnitDouble();
		group[2] = randomUnitDouble();

		group[3] = randomUnitDouble();
		group[4] = randomUnitDouble();
		group[5] = randomUnitDouble();

		group[6] = randomUnitDouble();
		group[7] = randomUnitDouble();
		group[8] = randomUnitDouble();

		double u = randomUnitDouble();
		double v = 1 - u;
		group[9]  = group[0] + (group[3] - group[0]) * u + (group[3] - group[0]) * v;
		group[10] = group[1] + (group[4] - group[1]) * u + (group[4] - group[1]) * v;
		group[11] = group[2] + (group[5] - group[2]) * u + (group[5] - group[2]) * v;

		u = randomUnitDouble();
		v = 1 - u;
		group[12] = group[0] + (group[3] - group[0]) * u + (group[3] - group[0]) * v;
		group[13] = group[1] + (group[4] - group[1]) * u + (group[4] - group[1]) * v;
		group[14] = group[2] + (group[5] - group[2]) * u + (group[5] - group[2]) * v;

		u = randomUnitDouble();
		v = 1 - u;
		group[15] = group[0] + (group[3] - group[0]) * u + (group[3] - group[0]) * v;
		group[16] = group[1] + (group[4] - group[1]) * u + (group[4] - group[1]) * v;
		group[17] = group[2] + (group[5] - group[2]) * u + (group[5] - group[2]) * v;
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
		if (EKIC::Collinear_3()(ref_points[6 * i], ref_points[6 * i + 1],
		                        ref_points[6 * i + 2]) ||
		    EKIC::Collinear_3()(ref_points[6 * i + 3], ref_points[6 * i + 4],
		                        ref_points[6 * i + 5]))
			continue;
		cgal_dummy += EKIC::Do_intersect_3()(
		  EKIC::Triangle_3(ref_points[6 * i], ref_points[6 * i + 1],
		                   ref_points[6 * i + 2]),
		  EKIC::Triangle_3(ref_points[6 * i + 3], ref_points[6 * i + 4],
		                   ref_points[6 * i + 5]));
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		if (CollinearPoints3D()(our_points[6 * i], our_points[6 * i + 1],
		                        our_points[6 * i + 2]) ||
		    CollinearPoints3D()(our_points[6 * i + 3], our_points[6 * i + 4],
		                        our_points[6 * i + 5]))
			continue;
		ours_dummy += Triangle3_Triangle3_DoIntersect()(
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
		  if (!EKIC::Collinear_3()(ref_points[6 * i], ref_points[6 * i + 1],
		                           ref_points[6 * i + 2]) &&
		      !EKIC::Collinear_3()(ref_points[6 * i + 3], ref_points[6 * i + 4],
		                           ref_points[6 * i + 5]))
		  {
			  bool ref_sign = EKIC::Do_intersect_3()(
			    EKIC::Triangle_3(ref_points[6 * i], ref_points[6 * i + 1],
			                     ref_points[6 * i + 2]),
			    EKIC::Triangle_3(ref_points[6 * i + 3], ref_points[6 * i + 4],
			                     ref_points[6 * i + 5]));

			  bool our_sign = Triangle3_Triangle3_DoIntersect()(
			    our_points[6 * i], our_points[6 * i + 1], our_points[6 * i + 2],
			    our_points[6 * i + 3], our_points[6 * i + 4], our_points[6 * i + 5]);

			  EXPECT_EQ(ref_sign, our_sign);
		  }
	  });
#endif
	delete[] atp;
}