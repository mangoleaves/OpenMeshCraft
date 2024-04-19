#include <chrono>
#include <iostream>

// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
// Arrangements
#include "OpenMeshCraft/Arrangements/MeshArrangements.h"

#include "OpenMeshCraft/Utils/StringUtils.h"

#include "test_utils.h"

using EIAC          = OMC::EIAC;
using TriSoupTraits = OMC::TriSoupTraits;

using Points    = typename TriSoupTraits::Points;
using Triangles = typename TriSoupTraits::Triangles;

using Arrangements = OMC::MeshArrangements<EIAC, TriSoupTraits>;

boost::property_tree::ptree omc_test_config;

int main(int argc, char *argv[])
{
	std::string filename;
	bool        output_stats  = false;
	bool        output_result = false;
	bool        verbose       = false;
	bool        has_config    = false;

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
		std::string param(argv[i]);
		if (param == "-v")
			verbose = true;
		else if (param == "-s")
			output_stats = true;
		else if (param == "-r")
			output_result = true;
		else if (OMC::starts_with(param, "--config="))
		{
			std::string config_path = param.substr(param.find_first_of('=') + 1);
			if (boost::filesystem::is_regular_file(
			      boost::filesystem::path(config_path)))
			{
				boost::property_tree::read_json(config_path, omc_test_config);
				has_config = true;
			}
			else
			{
				OMC_ASSERT(false, "test config file error.");
			}
		}
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

	read_mesh(filename, input_points, input_triangles, io_options);
	filename = filename.substr(filename.find_last_of("/\\") + 1);
	std::fstream fout;
	if (output_stats)
	{
		fout.open("./ours_time.txt", std::ios::out | std::ios::app);
		fout << filename << ",";
	}

	Arrangements arrangements(verbose);
	arrangements.addTriMeshAsInput(input_points, input_triangles);
	arrangements.setTriMeshAsOutput(result_points, result_triangles);

	if (has_config)
	{
		TEST_GET_CONFIG(Arrangements, TestDataSet);
		bool set_parameter = config.get<bool>("set_parameter", false);
		OMC::MeshArrangements_Config arr_config;
		if (set_parameter)
		{
			boost::property_tree::ptree &parameters = config.get_child("parameters");
			arr_config.tree_enlarge_ratio =
			  parameters.get<double>("tree_enlarge_ratio");
			arr_config.tree_adaptive_thres =
			  parameters.get<double>("tree_adaptive_thres");
			arr_config.tree_split_size_thres =
			  parameters.get<size_t>("tree_split_size_thres");

			arrangements.setConfig(arr_config);
		}
	}

	OMC::MeshArrangements_Stats &stats = arrangements.stats();

	auto start = OMC::Logger::elapse_reset();

	arrangements.meshArrangements(false, true);

	double time = OMC::Logger::elapsed(start).count();

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

	if (output_result)
		write_mesh(filename, result_points, result_triangles, io_options);
	return 0;
}