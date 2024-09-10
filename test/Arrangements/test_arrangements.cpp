#include <chrono>
#include <iostream>

// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"
// triangle soup traits
#include "OpenMeshCraft/Mesh/TriSoup.h"
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
	using Config       = OMC::MeshArrangements_Config;
};

/**
 * @brief check if it will crash.
 */
TEST_F(test_Arrangements, TestIfCrash)
{
	TEST_OUTPUT_DIRECTORY(Arrangements, TestIfCrash);
	TEST_GET_CONFIG(Arrangements, TestIfCrash);

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

	// set config for multi-threads
	tbb::global_control tbb_gc(
	  tbb::global_control::max_allowed_parallelism,
	  config.get<size_t>("thread_num", tbb::this_task_arena::max_concurrency()));

	// set config for arrangements
	Config arr_config;
	arr_config.verbose                = config.get<bool>("verbose");
	arr_config.output_explicit_result = config.get<bool>("num_vf");
	arr_config.ignore_same_mesh       = false;

	// Run the algorithm
	Arrangements arrangements;
	arrangements.setConfig(arr_config);
	arrangements.addTriMeshAsInput(input_points, input_triangles);
	arrangements.setTriMeshAsOutput(result_points, result_triangles);

	auto start = OMC::Logger::elapse_reset();
	arrangements.meshArrangements();
	double total_time = OMC::Logger::elapsed(start).count();

	// report & log
	std::cout << std::format("mesh arrangement uses {} s\n", total_time);
	std::cout << std::format("result vertices {} result faces {}\n",
	                         result_points.size(), result_triangles.size());
	// write result
	if (config.get<bool>("write"))
		write_mesh(outdir + filename, result_points, result_triangles, io_options);
}

/**
 * @brief Test on datasets to check if it will crash.
 */
TEST_F(test_Arrangements, TestDataSet)
{
	TEST_OUTPUT_DIRECTORY(Arrangements, TestDataSet);
	TEST_GET_CONFIG(Arrangements, TestDataSet);

	// Define IO
	OMC::IOOptions io_options;
	io_options.vertex_has_point = true;

	// Define mesh
	Points    input_points, result_points;
	Triangles input_triangles, result_triangles;

	// open file for logging
	std::string  log_path = outdir + "stats.txt";
	std::fstream log_file;
	log_file.open(log_path, std::ios::out | std::ios::app);

	log_file << "filename,pp,tree,ci,tr,time\n";

	// get directory of models
	std::string             models_dir = config.get<std::string>("models_dir");
	boost::filesystem::path model_dir_path(models_dir);
	boost::filesystem::directory_iterator endIter;

	// configure the arrangements
	bool verbose = config.get<bool>("verbose");
	bool num_vf  = config.get<bool>("num_vf");

	bool set_parameter = config.get<bool>("set_parameter", false);

	OMC::MeshArrangements_Config arr_config;
	arr_config.verbose                = verbose;
	arr_config.output_explicit_result = num_vf;
	arr_config.ignore_same_mesh       = false;
	if (set_parameter)
	{
		boost::property_tree::ptree &parameters = config.get_child("parameters");
		arr_config.tree_enlarge_ratio =
		  parameters.get<double>("tree_enlarge_ratio");
		arr_config.tree_adaptive_thres =
		  parameters.get<double>("tree_adaptive_thres");
		arr_config.tree_split_size_thres =
		  parameters.get<size_t>("tree_split_size_thres");
	}

	// configure the multi-threads
	tbb::global_control tbb_gc(
	  tbb::global_control::max_allowed_parallelism,
	  config.get<size_t>("thread_num", tbb::this_task_arena::max_concurrency()));

	size_t skip_step   = config.get<size_t>("skip_step", 0);
	size_t process_cnt = 0;
	for (boost::filesystem::directory_iterator iter(model_dir_path);
	     iter != endIter; iter++)
	{
		if (boost::filesystem::is_directory(*iter))
		{
			// do nothing continue
		}
		else
		{
			process_cnt += 1;
			if (skip_step != 0 && process_cnt % skip_step != 0)
				continue;

			std::cout << "processing " << iter->path().filename().string()
			          << std::endl;
			// read mesh
			read_mesh(iter->path().string(), input_points, input_triangles,
			          io_options);
			// run the algorithm
			Arrangements arrangements;
			arrangements.setConfig(arr_config);
			arrangements.addTriMeshAsInput(input_points, input_triangles);
			arrangements.setTriMeshAsOutput(result_points, result_triangles);
			OMC::MeshArrangements_Stats &stats = arrangements.stats();

			auto start = OMC::Logger::elapse_reset();
			arrangements.meshArrangements();
			double total_time = OMC::Logger::elapsed(start).count();
			// report & log
			std::cout << total_time << " s\n";

			log_file << std::fixed;
			log_file << iter->path().filename().string();
			log_file << "," << stats.pp_elapsed << "," << stats.tree_elapsed << ","
			         << stats.ci_elapsed << "," << stats.tr_elapsed << ","
			         << total_time << "," << result_points.size() << ","
			         << result_triangles.size() << std::endl;
		}
	}
	log_file.close();
}
