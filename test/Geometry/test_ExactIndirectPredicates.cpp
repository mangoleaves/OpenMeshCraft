#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/Predicates/IndirectPredicateDetails.h"
#include "OpenMeshCraft/NumberTypes/FPU.h"
#include "OpenMeshCraft/Utils/Macros.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Simple_cartesian.h"
#include "CGAL/intersection_3.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "gtest/gtest.h"

#include <execution>

// #define CHECK_EACH

class test_ExactIndirectPredicates : public testing::Test
{
protected:
	using EIAC               = OMC::EIAC;
	using Pnt2               = EIAC::Point2;
	using Pnt3               = EIAC::Point3;
	using Orient2D           = EIAC::Orient2D;
	using Orient3D           = EIAC::Orient3D;
	using OrientOn2D         = EIAC::OrientOn2D;
	using CollinearPoints3D  = EIAC::CollinearPoints3D;
	using LessThan3D         = EIAC::LessThan3D;
	using MaxCompInTriNormal = EIAC::MaxCompInTriNormal;
	using InCircle           = EIAC::InCircle;
	using InSphere           = EIAC::InSphere;

	using FPnt2 = OMC::Point2T<double>;
	using FPnt3 = OMC::Point3T<double>;

	using EKIC = CGAL::Exact_predicates_inexact_constructions_kernel;
};

inline double randomUnitDouble() { return ((double)rand()) / RAND_MAX; }

TEST_F(test_ExactIndirectPredicates, PlainPredicates2D)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.05;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point triplets
	double *atp = new double[(num_random_groups + num_degn_groups) * 6];

	// Create vector of random 2D point triplets
	double *num_groups = atp;
	for (size_t i = 0; i < num_random_groups * 6; i++)
		num_groups[i] = randomUnitDouble();

	// Create vector of 2D point aligned triplets
	double *degn_groups = atp + num_random_groups * 6;

	for (size_t i = 0; i < num_degn_groups; i++)
	{
		double p0x   = randomUnitDouble();
		double p0y   = randomUnitDouble();
		double p1x   = randomUnitDouble();
		double p1y   = randomUnitDouble();
		double delta = randomUnitDouble();
		double p2x   = p1x + (p1x - p0x) * delta;
		double p2y   = p1y + (p1y - p0y) * delta;

		degn_groups[i * 6]     = p0x;
		degn_groups[i * 6 + 1] = p0y;
		degn_groups[i * 6 + 2] = p1x;
		degn_groups[i * 6 + 3] = p1y;
		degn_groups[i * 6 + 4] = p2x;
		degn_groups[i * 6 + 5] = p2y;
	}

	std::vector<Pnt2>          our_points;
	std::vector<EKIC::Point_2> ref_points;
	our_points.reserve(3 * num_all_groups);
	ref_points.reserve(3 * num_all_groups);
	for (size_t i = 0; i < 3 * num_all_groups; ++i)
	{
		our_points.emplace_back(atp[2 * i], atp[2 * i + 1]);
		ref_points.emplace_back(atp[2 * i], atp[2 * i + 1]);
	}

	// Calculate predicates on them
	int  cgal_dummy = 0;
	auto start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += EKIC::Orientation_2()(ref_points[3 * i], ref_points[3 * i + 1],
		                               ref_points[3 * i + 2]);
	}
	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		OMC::Sign sign = Orient2D()(our_points[3 * i], our_points[3 * i + 1],
		                            our_points[3 * i + 2]);
		ours_dummy += static_cast<int>(sign);
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
		  int ref_sign = EKIC::Orientation_2()(
		    ref_points[3 * i], ref_points[3 * i + 1], ref_points[3 * i + 2]);
		  int our_sign = static_cast<int>(Orient2D()(
		    our_points[3 * i], our_points[3 * i + 1], our_points[3 * i + 2]));

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}

TEST_F(test_ExactIndirectPredicates, PlainPredicates3D)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.05;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 12];

	// Create vector of random 3D point tets
	double *num_groups = atp;
	for (size_t i = 0; i < num_random_groups * 12; i++)
		num_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned tets
	double  group[12];
	double *degn_groups = atp + num_random_groups * 12;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		group[0]       = randomUnitDouble();
		group[1]       = randomUnitDouble();
		group[2]       = randomUnitDouble();
		group[3]       = randomUnitDouble();
		group[4]       = randomUnitDouble();
		group[5]       = randomUnitDouble();
		group[6]       = randomUnitDouble();
		group[7]       = randomUnitDouble();
		group[8]       = randomUnitDouble();
		double delta = randomUnitDouble();
		group[9]       = group[3] + (group[3] - group[0]) * delta;
		group[10]      = group[4] + (group[4] - group[1]) * delta;
		group[11]      = group[5] + (group[5] - group[2]) * delta;
		for (int j = 0; j < 12; j++)
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
	auto start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy +=
		  EKIC::Orientation_3()(ref_points[4 * i], ref_points[4 * i + 1],
		                        ref_points[4 * i + 2], ref_points[4 * i + 3]);
	}
	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		OMC::Sign sign = Orient3D()(our_points[4 * i], our_points[4 * i + 1],
		                            our_points[4 * i + 2], our_points[4 * i + 3]);
		ours_dummy += static_cast<int>(sign);
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
		              int ref_sign = EKIC::Orientation_3()(
		                ref_points[4 * i], ref_points[4 * i + 1],
		                ref_points[4 * i + 2], ref_points[4 * i + 3]);
		              int our_sign = static_cast<int>(
		                Orient3D()(our_points[4 * i], our_points[4 * i + 1],
		                           our_points[4 * i + 2], our_points[4 * i + 3]));

		              EXPECT_EQ(ref_sign, our_sign);
	              });
#endif
	delete[] atp;
}

TEST_F(test_ExactIndirectPredicates, Coincident)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.05;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 6];

	// Create vector of random 3D point tets
	double *num_groups = atp;
	for (size_t i = 0; i < num_random_groups * 6; i++)
		num_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned pairs
	double  group[6];
	double *degn_groups = atp + num_random_groups * 6;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		group[0] = randomUnitDouble();
		group[1] = randomUnitDouble();
		group[2] = randomUnitDouble();

		group[3] = group[0];
		group[4] = group[1];
		group[5] = group[2];

		for (size_t j = 0; j < 6; j++)
			degn_groups[i * 6 + j] = group[j];
	}

	std::vector<Pnt3>          our_points;
	std::vector<EKIC::Point_3> ref_points;
	our_points.reserve(2 * num_all_groups);
	ref_points.reserve(2 * num_all_groups);
	for (size_t i = 0; i < 2 * num_all_groups; ++i)
	{
		our_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
		ref_points.emplace_back(atp[3 * i], atp[3 * i + 1], atp[3 * i + 2]);
	}

	// Calculate predicates on them
	int  cgal_dummy = 0;
	auto start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += EKIC::Equal_3()(ref_points[2 * i], ref_points[2 * i + 1]);
	}
	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	// Calculate predicates on them
	int ours_dummy = 0;
	start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += LessThan3D().coincident(our_points[2 * i], our_points[2 * i + 1]);
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
		              bool ref_sign =
		                EKIC::Equal_3()(ref_points[2 * i], ref_points[2 * i + 1]);
		              bool our_sign = LessThan3D().coincident(
		                our_points[2 * i], our_points[2 * i + 1]);

		              EXPECT_EQ(ref_sign, our_sign);
	              });
#endif
	delete[] atp;
}

TEST_F(test_ExactIndirectPredicates, Misaligned)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.05;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 9];

	// Create vector of random 3D point tets
	double *num_groups = atp;
	for (size_t i = 0; i < num_random_groups * 9; i++)
		num_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned tets
	double  group[9];
	double *degn_groups = atp + num_random_groups * 9;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		group[0] = randomUnitDouble();
		group[1] = randomUnitDouble();
		group[2] = randomUnitDouble();
		group[3] = randomUnitDouble();
		group[4] = randomUnitDouble();
		group[5] = randomUnitDouble();

		const double delta = randomUnitDouble();

		group[6] = group[3] + (group[3] - group[0]) * delta;
		group[7] = group[4] + (group[4] - group[1]) * delta;
		group[8] = group[5] + (group[5] - group[2]) * delta;

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
	int cgal_dummy      = 0;
	auto start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += !EKIC::Collinear_3()(ref_points[3 * i], ref_points[3 * i + 1],
		                              ref_points[3 * i + 2]);
	}
	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += CollinearPoints3D().misaligned(
		  our_points[3 * i], our_points[3 * i + 1], our_points[3 * i + 2]);
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
		  bool ref_sign = !EKIC::Collinear_3()(
		    ref_points[3 * i], ref_points[3 * i + 1], ref_points[3 * i + 2]);
		  bool our_sign = CollinearPoints3D().misaligned(
		    our_points[3 * i], our_points[3 * i + 1], our_points[3 * i + 2]);

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}

TEST_F(test_ExactIndirectPredicates, InCircle)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.05;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 15];

	// Create vector of random 3D point tets
	double *num_groups = atp;
	for (size_t i = 0; i < num_random_groups * 8; i++)
		num_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned tets
	double  group[15];
	double *degn_groups = atp + num_random_groups * 8;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		group[0] = randomUnitDouble();
		group[1] = randomUnitDouble();

		group[2] = randomUnitDouble();
		group[3] = randomUnitDouble();

		group[4] = randomUnitDouble();
		group[5] = randomUnitDouble();

		// not real degenerate
		const double delta = randomUnitDouble();
		group[6]             = group[4] + (group[4] - group[1]) * delta;
		group[7]             = group[5] + (group[5] - group[2]) * delta;

		for (size_t j = 0; j < 8; j++)
			degn_groups[i * 8 + j] = group[j];
	}

	std::vector<Pnt2>          our_points;
	std::vector<EKIC::Point_2> ref_points;
	our_points.reserve(num_all_groups * 4);
	ref_points.reserve(num_all_groups * 4);
	for (size_t i = 0; i < 4 * num_all_groups; ++i)
	{
		our_points.emplace_back(atp[2 * i], atp[2 * i + 1]);
		ref_points.emplace_back(atp[2 * i], atp[2 * i + 1]);
	}

	// Calculate predicates on them
	int  cgal_dummy = 0;
	auto start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += EKIC::Side_of_oriented_circle_2()(
		  ref_points[4 * i], ref_points[4 * i + 1], ref_points[4 * i + 2],
		  ref_points[4 * i + 3]);
	}

	std::cout << "CGAL elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += static_cast<int>(
		  InCircle()(our_points[4 * i], our_points[4 * i + 1],
		             our_points[4 * i + 2], our_points[4 * i + 3]));
	}

	std::cout << "OpenMeshCraft IPreds elapsed time: "
	          << OMC::Logger::elapsed(start).count() << "s\n";
	std::cout << "Dummy sum: " << ours_dummy << "\n";

	EXPECT_EQ(cgal_dummy, ours_dummy);

#ifdef CHECK_EACH
	std::vector indices(num_all_groups, 0);
	std::iota(indices.begin(), indices.end(), 0);
	std::for_each(std::execution::seq, indices.begin(), indices.end(),
	              [&](size_t i)
	              {
		              int ref_sign = EKIC::Side_of_oriented_circle_2()(
		                ref_points[4 * i], ref_points[4 * i + 1],
		                ref_points[4 * i + 2], ref_points[4 * i + 3]);

		              int our_sign = static_cast<int>(
		                InCircle()(our_points[4 * i], our_points[4 * i + 1],
		                           our_points[4 * i + 2], our_points[4 * i + 3]));

		              EXPECT_EQ(ref_sign, our_sign);
	              });
#endif
	delete[] atp;
}

TEST_F(test_ExactIndirectPredicates, InSphere)
{
	size_t       num_all_groups    = 1000000;
	const double perc_degn         = 0.05;
	const size_t num_random_groups = (size_t)(num_all_groups * (1.0 - perc_degn));
	const size_t num_degn_groups   = (size_t)(num_all_groups * perc_degn);
	num_all_groups                 = num_random_groups + num_degn_groups;
	srand(0);

	// Create vector of all the point tets
	double *atp = new double[(num_random_groups + num_degn_groups) * 15];

	// Create vector of random 3D point tets
	double *num_groups = atp;
	for (size_t i = 0; i < num_random_groups * 15; i++)
		num_groups[i] = randomUnitDouble();

	// Create vector of 3D point aligned tets
	double  group[15];
	double *degn_groups = atp + num_random_groups * 15;
	for (size_t i = 0; i < num_degn_groups; i++)
	{
		group[0] = randomUnitDouble();
		group[1] = randomUnitDouble();
		group[2] = randomUnitDouble();

		group[3] = randomUnitDouble();
		group[4] = randomUnitDouble();
		group[5] = randomUnitDouble();

		group[6] = randomUnitDouble();
		group[7] = randomUnitDouble();
		group[8] = randomUnitDouble();

		const double delta = randomUnitDouble();
		group[9]             = group[3] + (group[3] - group[0]) * delta;
		group[10]            = group[4] + (group[4] - group[1]) * delta;
		group[11]            = group[5] + (group[5] - group[2]) * delta;

		group[12] = randomUnitDouble();
		group[13] = randomUnitDouble();
		group[14] = randomUnitDouble();
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
	int cgal_dummy      = 0;
	auto start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		cgal_dummy += EKIC::Side_of_oriented_sphere_3()(
		  ref_points[5 * i], ref_points[5 * i + 1], ref_points[5 * i + 2],
		  ref_points[5 * i + 3], ref_points[5 * i + 4]);
	}

	std::cout << "IPreds elapsed time: " << OMC::Logger::elapsed(start).count()
	          << "s\n";
	std::cout << "Dummy sum: " << cgal_dummy << "\n";

	int ours_dummy = 0;
	start = OMC::Logger::elapse_reset();
	for (size_t i = 0; i < num_all_groups; ++i)
	{
		ours_dummy += static_cast<int>(InSphere()(
		  our_points[5 * i], our_points[5 * i + 1], our_points[5 * i + 2],
		  our_points[5 * i + 3], our_points[5 * i + 4]));
	}

	std::cout << "OpenMeshCraft IPreds elapsed time: "
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
		  int ref_sign = EKIC::Side_of_oriented_sphere_3()(
		    ref_points[5 * i], ref_points[5 * i + 1], ref_points[5 * i + 2],
		    ref_points[5 * i + 3], ref_points[5 * i + 4]);

		  int our_sign = static_cast<int>(InSphere()(
		    our_points[5 * i], our_points[5 * i + 1], our_points[5 * i + 2],
		    our_points[5 * i + 3], our_points[5 * i + 4]));

		  EXPECT_EQ(ref_sign, our_sign);
	  });
#endif
	delete[] atp;
}