#include <chrono>
#include <iostream>

// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
// Boolean
#include "OpenMeshCraft/Boolean/MeshBoolean.h"

#include "test_utils.h"

/*
 * All tests about Mesh Boolean are put here,
 * until it's better to separate tests :D.
 */

class test_Boolean : public testing::Test
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

	using Boolean = OMC::MeshBoolean<EIAC, TriSoupTraits>;
};

/**
 * @brief Test all boolean operations to check if it will crash.
 */
TEST_F(test_Boolean, TestIfCrash)
{
	TEST_OUTPUT_DIRECTORY(Boolean, TestIfCrash);
	TEST_GET_CONFIG(Boolean, TestIfCrash);

	// Define IO
	IOOptions io_options;
	io_options.vertex_has_point = true;

	Points    points1, points2, result_points;
	Triangles triangles1, triangles2, result_triangles;

	// read mesh
	read_mesh(config.get<std::string>("mesh_path_1"), points1,
	                         triangles1, io_options);
	read_mesh(config.get<std::string>("mesh_path_2"), points2,
	                         triangles2, io_options);

	Boolean boolean(/*verbose*/ true);
	boolean.addTriMeshAsInput(points1, triangles1);
	boolean.addTriMeshAsInput(points2, triangles2);
	boolean.setTriMeshAsOutput(result_points, result_triangles);
	boolean.computeLabels();

	{
		boolean.Union();
		write_mesh(outdir + "Union.obj", result_points,
		                          result_triangles, io_options);
	}

	{
		boolean.Intersection();
		write_mesh(outdir + "Intersection.obj", result_points,
		                          result_triangles, io_options);
	}

	{
		boolean.Xor();
		write_mesh(outdir + "Xor.obj", result_points,
		                          result_triangles, io_options);
	}

	{
		boolean.Subtraction();
		write_mesh(outdir + "Subtraction.obj", result_points,
		                          result_triangles, io_options);
	}
}
