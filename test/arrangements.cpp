#include <chrono>
#include <iostream>

// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
// Arrangements
#include "OpenMeshCraft/Arrangements/MeshArrangements.h"

#include "OpenMeshCraft/Utils/StringUtils.h"

#include "test_utils.h"

using index_t       = OMC::index_t;
using EIAC          = OMC::EIAC;
using TriSoupTraits = OMC::TriSoupTraits;

using Points    = typename TriSoupTraits::Points;
using Triangles = typename TriSoupTraits::Triangles;

using Arrangements = OMC::MeshArrangements<EIAC, TriSoupTraits>;

int main(int argc, char *argv[])
{
	std::string filename;
	bool        output_stats  = false;
	bool        output_result = false;
	bool        verbose       = false;

	auto print_help = []()
	{
		std::cout << "*.exe input_file -v -s -r\n"
		             "-v verbose\n"
		             "-s output stats to default file\n"
		             "-r output result model to default file\n";
		exit(1);
	};

	if (argc > 1)
		filename = argv[1];
	else
		print_help();
	for (int i = 2; i < argc; i++)
	{
		if (std::string(argv[i]) == "-v")
			verbose = true;
		else if (std::string(argv[i]) == "-s")
			output_stats = true;
		else if (std::string(argv[i]) == "-r")
			output_result = true;
		else
			print_help();
	}

	// Define IO
	IOOptions io_options;
	io_options.vertex_has_point = true;

	// Define mesh
	Points    input_points, result_points;
	Triangles input_triangles, result_triangles;

	read_mesh(filename, input_points, input_triangles, io_options);
#if 0
	filename = filename.substr(filename.find_last_of("/\\") + 1);
	std::fstream fout;
	if (output_stats)
	{
		fout.open("./ours_time.txt", std::ios::out | std::ios::app);
		fout << filename << ",";
	}

	auto start = OMC::Logger::elapse_reset();

	Arrangements arrangements(verbose);

	OMC::MeshArrangements_Stats &stats = arrangements.stats();

	arrangements.addTriMeshAsInput(input_points, input_triangles);
	arrangements.setTriMeshAsOutput(result_points, result_triangles);
	arrangements.meshArrangements(false, true);

	double time = OMC::Logger::elapsed(start).count();

	if (output_stats)
	{
		fout << stats.pp_elapsed << "," << stats.di_elapsed << ","
		     << stats.cn_elapsed << "," << stats.ci_elapsed << ","
		     << stats.tr_elapsed << "," << time << ","
		     << OMC::getPeakMegabytesUsed() << "," << result_points.size() << ","
		     << result_triangles.size() << "\n";
		std::cout << filename << ": " << time << "s, "
		          << OMC::getPeakMegabytesUsed() << " MB\n";
		fout.close();
	}

	if (output_result)
		write_mesh(filename, result_points, result_triangles, io_options);
#endif
	return 0;
}