#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/Intersection/Segment3_Point3.h"
#include "OpenMeshCraft/Geometry/Intersection/Segment3_Segment3.h"
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

class test_Intersection_Segment3 : public testing::Test
{
protected:
	using EIAC = OMC::EIAC;
	using Pnt2 = EIAC::Point2;
	using Pnt3 = EIAC::Point3;

	using Segment3_Point3_DoIntersect = OMC::Segment3_Point3_Do_Intersect<EIAC>;
	using Segment3_Segment3_DoIntersect =
	  OMC::Segment3_Segment3_Do_Intersect<EIAC>;

	using FPnt2 = OMC::Point2T<double>;
	using FPnt3 = OMC::Point3T<double>;

	using EKIC = CGAL::Exact_predicates_inexact_constructions_kernel;
};

inline double randomUnitDouble() { return ((double)rand()) / RAND_MAX; }

TEST_F(test_Intersection_Segment3, pointSegment)
{
	size_t       num_all_groups    = 100000;
	const double perc_degn         = 0.5;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	// reset random seed for this test everytime
	srand(0);

	// Create vector of all segment and point
	double *atp = new double[(num_random_groups + num_degn_groups) * 9];

	// Create vector of random segment and point
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 9; i++)
		random_groups[i] = randomUnitDouble();

	// Create vector of random segment and (not strictly) colinear point
	double  group[9];
	double *degn_groups = atp + num_random_groups * 9;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		// clang-format off
		// random segment
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		// random colinear point
		double delta = randomUnitDouble();
		group[6]     = group[0] + (group[3] - group[0]) * delta * 2.0;
		group[7]     = group[1] + (group[4] - group[1]) * delta * 2.0;
		group[8]     = group[2] + (group[5] - group[2]) * delta * 2.0;
		// clang-format on
		for (size_t j = 0; j < 9; j++)
			degn_groups[i * 9 + j] = group[j];
	}

	std::vector<Pnt3>          our_points;
	std::vector<EKIC::Point_3> ref_points;
	our_points.reserve(num_all_groups * 3);
	ref_points.reserve(num_all_groups * 3);
	for (size_t i = 0; i < 3 * num_all_groups; ++i)
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
		  EKIC::Segment_3(ref_points[3 * i], ref_points[3 * i + 1]),
		  ref_points[3 * i + 2]);
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += Segment3_Point3_DoIntersect()(
		  our_points[3 * i + 1], our_points[3 * i + 2], our_points[3 * i]);
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
		                EKIC::Segment_3(ref_points[3 * i], ref_points[3 * i + 1]),
		                ref_points[3 * i + 2]);
		              bool our_sign = Segment3_Point3_DoIntersect()(
		                our_points[3 * i + 1], our_points[3 * i + 2],
		                our_points[3 * i]);

		              EXPECT_EQ(ref_sign, our_sign);
	              });
#endif
	delete[] atp;
}

TEST_F(test_Intersection_Segment3, segmentSegment)
{
	size_t       num_all_groups    = 100000;
	const double perc_degn         = 0.5;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	// reset random seed for this test everytime
	srand(0);

	// Create vector of all segments
	double *atp = new double[(num_random_groups + num_degn_groups) * 12];

	// Create two random segments
	double *random_groups = atp;
	for (size_t i = 0; i < num_random_groups * 12; i++)
		random_groups[i] = randomUnitDouble();

	// Create two random and coplanar segments
	double  group[12];
	double *degn_groups = atp + num_random_groups * 12;
	for (size_t i = 0; i < size_t(num_degn_groups * 0.5); i++)
	{
		// clang-format off
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = randomUnitDouble(); group[7] = randomUnitDouble(); group[8] = randomUnitDouble();
		const double delta = randomUnitDouble();
		group[9]  = group[6] + (group[3] - group[6]) * delta + (group[0] - group[6]) * delta;
		group[10] = group[7] + (group[4] - group[7]) * delta + (group[1] - group[7]) * delta;
		group[11] = group[8] + (group[5] - group[8]) * delta + (group[2] - group[8]) * delta;
		// clang-format on

		for (size_t j = 0; j < 12; j++)
			degn_groups[i * 12 + j] = group[j];
	}
	// Create two random segments originating from the same point
	for (size_t i = size_t(num_degn_groups * 0.5);
	     i < size_t(num_degn_groups * 0.75); i++)
	{
		// clang-format off
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		group[6] = group[0];           group[7] = group[1];           group[8] = group[2];
		const double delta = randomUnitDouble() + 1e-6;
		group[9]  = group[6] + (group[3] - group[6]) * delta + (group[0] - group[6]) * delta;
		group[10] = group[7] + (group[4] - group[7]) * delta + (group[1] - group[7]) * delta;
		group[11] = group[8] + (group[5] - group[8]) * delta + (group[2] - group[8]) * delta;
		// clang-format on

		for (size_t j = 0; j < 12; j++)
			degn_groups[i * 12 + j] = group[j];
	}
	// Create two random and (not strictly) colinear segments
	for (size_t i = size_t(num_degn_groups * 0.75); i < num_degn_groups; i++)
	{
		// clang-format off
		group[0] = randomUnitDouble(); group[1] = randomUnitDouble(); group[2] = randomUnitDouble();
		group[3] = randomUnitDouble(); group[4] = randomUnitDouble(); group[5] = randomUnitDouble();
		double delta = randomUnitDouble() + 1e-6;
		group[6] = group[0] + (group[3] - group[0]) * delta;
		group[7] = group[1] + (group[4] - group[1]) * delta;
		group[8] = group[2] + (group[5] - group[2]) * delta;
		delta = randomUnitDouble() + 1e-6;
		group[9]  = group[6] + (group[3] - group[0]) * delta;
		group[10] = group[7] + (group[4] - group[1]) * delta;
		group[11] = group[8] + (group[5] - group[2]) * delta;
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
		  EKIC::Segment_3(ref_points[4 * i], ref_points[4 * i + 1]),
		  EKIC::Segment_3(ref_points[4 * i + 2], ref_points[4 * i + 3]));
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start          = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += Segment3_Segment3_DoIntersect()(
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
	std::for_each(
	  std::execution::seq, indices.begin(), indices.end(),
	  [&](size_t i)
	  {
		  bool ref_sign = EKIC::Do_intersect_3()(
		    EKIC::Segment_3(ref_points[4 * i], ref_points[4 * i + 1]),
		    EKIC::Segment_3(ref_points[4 * i + 2], ref_points[4 * i + 3]));
		  bool our_sign = Segment3_Segment3_DoIntersect()(
		    our_points[4 * i], our_points[4 * i + 1], our_points[4 * i + 2],
		    our_points[4 * i + 3]);

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}