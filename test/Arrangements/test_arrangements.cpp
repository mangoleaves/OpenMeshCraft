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
 * @brief Check if the mesh arrangement algorithm will crash.
 */
TEST_F(test_Arrangements, TestIfCrash)
{
	TEST_OUTPUT_DIRECTORY(Arrangements, TestIfCrash);
	TEST_GET_CONFIG(Arrangements, TestIfCrash);

	// Define IO options
	IOOptions io_options;
	io_options.vertex_has_point = true;

	// Define mesh containers
	Points    input_points, result_points;
	Triangles input_triangles, result_triangles;

	// Read mesh from file
	std::string dir      = config.get<std::string>("dir");
	std::string filename = config.get<std::string>("filename");

	read_mesh(dir + filename, input_points, input_triangles, io_options);

	// Set configuration for multi-threading
	tbb::global_control tbb_gc(
	  tbb::global_control::max_allowed_parallelism,
	  config.get<size_t>("thread_num", tbb::this_task_arena::max_concurrency()));

	// Set configuration for arrangements
	Config arr_config;
	arr_config.verbose                = config.get<bool>("verbose");
	arr_config.output_explicit_result = config.get<bool>("explicit");
	arr_config.ignore_same_mesh       = false;

	// Run the mesh arrangement algorithm
	Arrangements arrangements;
	arrangements.setConfig(arr_config);
	arrangements.addTriMeshAsInput(input_points, input_triangles);
	arrangements.setTriMeshAsOutput(result_points, result_triangles);

	auto start = OMC::Logger::elapse_reset();
	arrangements.meshArrangements();
	double total_time = OMC::Logger::elapsed(start).count();

	// Report and log the results
	std::cout << std::format("Mesh arrangement took {} seconds\n", total_time);
	std::cout << std::format("Result vertices: {} Result faces: {}\n",
	                         result_points.size(), result_triangles.size());
	// Write the result to a file if configured to do so
	if (config.get<bool>("write"))
		write_mesh(outdir + filename, result_points, result_triangles, io_options);
}

/**
 * @brief Test the mesh arrangement algorithm on a dataset to check for crashes.
 */
TEST_F(test_Arrangements, TestDataSet)
{
	TEST_OUTPUT_DIRECTORY(Arrangements, TestDataSet);
	TEST_GET_CONFIG(Arrangements, TestDataSet);

	// Define IO options
	OMC::IOOptions io_options;
	io_options.vertex_has_point = true;

	// Define mesh containers
	Points    input_points, result_points;
	Triangles input_triangles, result_triangles;

	// Open file for logging
	std::string  log_path = outdir + "stats.txt";
	std::fstream log_file;
	log_file.open(log_path, std::ios::out | std::ios::app);

	log_file << "filename,pp,tree,ci,tr,time\n";

	// Get directory of models
	std::string           models_dir = config.get<std::string>("models_dir");
	std::filesystem::path model_dir_path(models_dir);
	std::filesystem::directory_iterator endIter;

	// Configure the arrangements
	OMC::MeshArrangements_Config arr_config;
	arr_config.verbose                = config.get<bool>("verbose");
	arr_config.output_explicit_result = config.get<bool>("explicit");
	arr_config.ignore_same_mesh       = false;
	if (config.get<bool>("set_parameter", false))
	{
		boost::property_tree::ptree &parameters = config.get_child("parameters");
		arr_config.tree_enlarge_ratio =
		  parameters.get<double>("tree_enlarge_ratio");
		arr_config.tree_adaptive_thres =
		  parameters.get<double>("tree_adaptive_thres");
		arr_config.tree_split_size_thres =
		  parameters.get<size_t>("tree_split_size_thres");
	}

	// Configure multi-threading
	tbb::global_control tbb_gc(
	  tbb::global_control::max_allowed_parallelism,
	  config.get<size_t>("thread_num", tbb::this_task_arena::max_concurrency()));

	size_t skip_step   = config.get<size_t>("skip_step", 0);
	size_t process_cnt = 0;
	for (std::filesystem::directory_iterator iter(model_dir_path);
	     iter != endIter; iter++)
	{
		if (std::filesystem::is_directory(*iter))
		{
			// Skip directories
			continue;
		}
		else
		{
			process_cnt += 1;
			if (skip_step != 0 && process_cnt % skip_step != 0)
				continue;

			std::cout << "Processing " << iter->path().filename().string()
			          << std::endl;
			// Read mesh from file
			read_mesh(iter->path().string(), input_points, input_triangles,
			          io_options);
			// Run the mesh arrangement algorithm
			Arrangements arrangements;
			arrangements.setConfig(arr_config);
			arrangements.addTriMeshAsInput(input_points, input_triangles);
			arrangements.setTriMeshAsOutput(result_points, result_triangles);
			OMC::MeshArrangements_Stats &stats = arrangements.stats();

			auto start = OMC::Logger::elapse_reset();
			arrangements.meshArrangements();
			double total_time = OMC::Logger::elapsed(start).count();
			// Report and log the results
			std::cout << total_time << " seconds\n";

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
