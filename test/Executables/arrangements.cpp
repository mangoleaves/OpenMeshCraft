/**
 * @file arrangements.cpp
 * @brief An independent executable file to run mesh arrangements.
 */
#include <chrono>
#include <iostream>

// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
// Traits
#include "OpenMeshCraft/Mesh/TriSoup.h"
// Arrangements
#include "OpenMeshCraft/Arrangements/MeshArrangements.h"

#include "OpenMeshCraft/Utils/StringUtils.h"

#include "test_utils.h"

using EIAC          = OMC::EIAC;
using TriSoupTraits = OMC::TriSoupTraits;

using Points    = typename TriSoupTraits::Points;
using Triangles = typename TriSoupTraits::Triangles;

using Arrangements       = OMC::MeshArrangements<EIAC, TriSoupTraits>;
using ArrangementsConfig = OMC::MeshArrangements_Config;

boost::property_tree::ptree omc_test_config;

int main(int argc, char *argv[])
{
	std::string filename;
	bool        output_stats  = false;

	ArrangementsConfig arr_config;

	auto print_help = []()
	{
		std::cout
		  << "exe input_file -v -s -r -p=n --tree_leaf=n --tree_adaptive=f\n"
		     "-v verbose\n"
		     "-s output stats to default file\n"
		     "-r output result model to default file\n"
		     "-p=n parallel threads number, n is an integer\n"
		     "--tree_leaf=n tree's leaf node size, n is an integer\n"
		     "--tree_adaptive=f tree's adaptive threshold, f is an float number\n";
		exit(1);
	};

	if (argc > 1)
		filename = argv[1];
	else
		print_help();
	for (int i = 2; i < argc; i++)
	{
		std::string param(argv[i]);
		if (param == "-v")
			arr_config.verbose = true;
		else if (param == "-s")
			output_stats = true;
		else if (param == "-r")
			arr_config.output_explicit_result = true;
		else if (OMC::starts_with(param, "-p="))
		{
			int                 threads_num = std::stoi(param.substr(3));
			tbb::global_control tbb_gc(tbb::global_control::max_allowed_parallelism,
			                           threads_num);
		}
		else if (OMC::starts_with(param, "--tree_leaf="))
			arr_config.tree_split_size_thres = std::stoi(param.substr(12));
		else if (OMC::starts_with(param, "--tree_adaptive="))
			arr_config.tree_adaptive_thres = std::stof(param.substr(16));
		else
			print_help();
	}

	// Define IO
	IOOptions io_options;
	io_options.vertex_has_point = true;
	io_options.stl_binary       = true;

	// Define mesh
	Points    input_points, result_points;
	Triangles input_triangles, result_triangles;

	// read mesh
	read_mesh(filename, input_points, input_triangles, io_options);
	// prepare to output stats
	filename = filename.substr(filename.find_last_of("/\\") + 1);
	std::fstream fout;
	if (output_stats)
	{
		fout.open("./ours_time.txt", std::ios::out | std::ios::app);
		fout << filename << ",";
	}

	// set config, input and output
	Arrangements arrangements;
	arrangements.setConfig(arr_config);
	arrangements.addTriMeshAsInput(input_points, input_triangles);
	arrangements.setTriMeshAsOutput(result_points, result_triangles);
	// get stats
	OMC::MeshArrangements_Stats &stats = arrangements.stats();
	// run the algorithm
	auto start = OMC::Logger::elapse_reset();
	arrangements.meshArrangements();
	double time = OMC::Logger::elapsed(start).count();
	// output stats
	if (output_stats)
	{
		fout << std::fixed;
		fout << stats.pp_elapsed << "," << stats.tree_elapsed << ","
		     << stats.ci_elapsed << "," << stats.tr_elapsed << "," << time << ","
		     << OMC::getPeakMegabytesUsed() << "," << result_points.size() << ","
		     << result_triangles.size() << "\n";
		std::cout << filename << ": " << time << "s, "
		          << OMC::getPeakMegabytesUsed() << " MB\n";
		fout.close();
	}
	// output result
	if (arr_config.output_explicit_result)
		write_mesh(filename, result_points, result_triangles, io_options);
	return 0;
}