#include "OpenMeshCraft/Boolean/MeshBoolean.h"
#include "OpenMeshCraft/Mesh/IO/OBJReader.h"
#include "OpenMeshCraft/Mesh/IO/OBJWriter.h"
#include "OpenMeshCraft/Mesh/IO/STLReader.h"
#include "OpenMeshCraft/Mesh/IO/STLWriter.h"

#include "example_utils.h"

#include <Eigen/Dense>

// geometry kernel, hidden for users, just use it. :D
using EIAC          = OMC::EIAC;
// Triangle soup traits
using TriSoupTraits = OMC::TriSoupTraits;

// points and triangles from traits
using Points    = typename TriSoupTraits::Points;
using Triangles = typename TriSoupTraits::Triangles;

// Define reader and writer
using OBJReader = OMC::OBJReader<TriSoupTraits>;
using OBJWriter = OMC::OBJWriter<TriSoupTraits>;
using STLReader = OMC::STLReader<TriSoupTraits>;
using STLWriter = OMC::STLWriter<TriSoupTraits>;
using IOOptions = OMC::IOOptions;

// Define boolean
using MeshBoolean = OMC::MeshBoolean<EIAC, TriSoupTraits>;

void read_mesh(const std::string &filename, Points &points,
               Triangles &triangles, IOOptions &io_options)
{
	if (OMC::ends_with(filename, ".obj"))
	{
		OBJReader reader;
		reader.read(filename, io_options);

		points    = std::move(reader.m_points);
		triangles = std::move(reader.m_triangles);
	}
	else if (OMC::ends_with(filename, ".stl"))
	{
		STLReader reader;
		reader.read(filename, io_options);

		points    = std::move(reader.m_points);
		triangles = std::move(reader.m_triangles);
	}
}

void write_mesh(const std::string &filename, const Points &points,
                const Triangles &triangles, IOOptions &io_options)
{
	OBJWriter writer;
	writer.m_points    = std::move(points);
	writer.m_triangles = std::move(triangles);
	writer.write(filename, io_options, 10);
};

void example_MeshBoolean()
{
	OUTPUT_DIRECTORY(Boolean, Examples);

	// Define IO
	IOOptions io_options;
	io_options.vertex_has_point = true;

	Points    points1, points2, result_points;
	Triangles triangles1, triangles2, result_triangles;

	// read mesh
	read_mesh("./data/boolean_data/bunny25k.obj", points1, triangles1,
	          io_options);
	read_mesh("./data/boolean_data/cow.obj", points2, triangles2, io_options);

	MeshBoolean boolean(/*verbose*/ true);
	// set input meshes
	boolean.addTriMeshAsInput(points1, triangles1);
	boolean.addTriMeshAsInput(points2, triangles2);
	// set output mesh
	boolean.setTriMeshAsOutput(result_points, result_triangles);

	// run arrangements and labeling.
	// the results will be temporarily cached, until be cleared or destroyed.
	boolean.computeLabels();

	// after arrangements and labeling done, we can
	// apply differenct boolean operations faster.
	// (NOTE: must call computeLabels before boolean operation)
	{
		boolean.Union();
		make_file_writable(outdir + "OMC_Union.obj");
		write_mesh(outdir + "OMC_Union.obj", result_points, result_triangles,
		           io_options);
	}

	{
		boolean.Intersection();
		make_file_writable(outdir + "OMC_Intersection.obj");
		write_mesh(outdir + "OMC_Intersection.obj", result_points, result_triangles,
		           io_options);
	}

	{
		boolean.Xor();
		make_file_writable(outdir + "OMC_Xor.obj");
		write_mesh(outdir + "OMC_Xor.obj", result_points, result_triangles,
		           io_options);
	}

	{
		boolean.Subtraction();
		make_file_writable(outdir + "OMC_Subtraction.obj");
		write_mesh(outdir + "OMC_Subtraction.obj", result_points, result_triangles,
		           io_options);
	}
}
