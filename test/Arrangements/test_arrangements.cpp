#include <chrono>
#include <iostream>

// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
// Arrangements
#include "OpenMeshCraft/Arrangements/MeshArrangements.h"

#include "test_utils.h"

/*
 * All tests about Mesh Arrangements are put here,
 * until it's better to separate tests :D.
 */

class test_Arrangements : public testing::Test
{
protected:
	void SetUp() override {}

	void TearDown() override {}

protected:
	using index_t       = OMC::index_t;
	using EIAC          = OMC::EIAC;
	using TriSoupTraits = OMC::TriSoupTraits;

	using Points    = typename TriSoupTraits::Points;
	using Triangles = typename TriSoupTraits::Triangles;

	using Arrangements = OMC::MeshArrangements<EIAC, TriSoupTraits>;
};

/**
 * @brief check if it will crash.
 */
TEST_F(test_Arrangements, TestIfCrash)
{
	TEST_OUTPUT_DIRECTORY(Arrangements, TestIfComplexArrangementCrash);
	TEST_GET_CONFIG(Arrangements, TestIfComplexArrangementCrash);

	// Define IO
	IOOptions io_options;
	io_options.vertex_has_point = true;

	// Define mesh
	Points    input_points, result_points;
	Triangles input_triangles, result_triangles;

	// read mesh
	std::string dir      = config.get<std::string>("dir");
	std::string filename = config.get<std::string>("filename");

	read_mesh(dir + filename, input_points, input_triangles, io_options);

	auto start = OMC::Logger::elapse_reset();

	Arrangements arrangements(/*verbose*/ config.get<bool>("verbose"));

	arrangements.addTriMeshAsInput(input_points, input_triangles);
	arrangements.setTriMeshAsOutput(result_points, result_triangles);
	arrangements.meshArrangements(false, true);

	std::cout << "complex mesh arrangement uses "
	          << OMC::Logger::elapsed(start).count() << " s\n";

	write_mesh(outdir + filename, result_points, result_triangles, io_options);
}
