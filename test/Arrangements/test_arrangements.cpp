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

	auto start = OMC::Logger::elapse_reset();

	Arrangements arrangements(/*verbose*/ config.get<bool>("verbose"));

	arrangements.addTriMeshAsInput(input_points, input_triangles);
	arrangements.setTriMeshAsOutput(result_points, result_triangles);
	arrangements.meshArrangements(false, true);

	std::cout << "complex mesh arrangement uses "
	          << OMC::Logger::elapsed(start).count() << " s\n";

	write_mesh(outdir + filename, result_points, result_triangles, io_options);
}

/**
 * @brief Test all boolean operations to check if it will crash.
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

	std::string  log_path = outdir + "stats.txt";
	std::fstream log_file;
	log_file.open(log_path, std::ios::out | std::ios::app);

	log_file << "filename,pp,tree,di,cn,ci,tr,time\n";

	std::string             models_dir = config.get<std::string>("models_dir");
	boost::filesystem::path model_dir_path(models_dir);
	boost::filesystem::directory_iterator endIter;

	bool verbose = config.get<bool>("verbose");

	bool set_parameter = config.get<bool>("set_parameter", false);

	OMC::MeshArrangements_Config arr_config;
	arr_config.tree_enlarge_ratio    = 1.1;
	arr_config.tree_adaptive_thres   = 0.1;
	arr_config.tree_parallel_scale   = 10000;
	arr_config.tree_split_size_thres = 1000;
	if (set_parameter)
	{
		boost::property_tree::ptree &parameters = config.get_child("parameters");
		arr_config.tree_enlarge_ratio =
		  parameters.get<double>("tree_enlarge_ratio");
		arr_config.tree_adaptive_thres =
		  parameters.get<double>("tree_adaptive_thres");
		arr_config.tree_parallel_scale =
		  parameters.get<size_t>("tree_parallel_scale");
		arr_config.tree_split_size_thres =
		  parameters.get<size_t>("tree_split_size_thres");
	}

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
			read_mesh(iter->path().string(), input_points, input_triangles,
			          io_options);

			Arrangements arrangements(/*verbose*/ verbose);
			arrangements.addTriMeshAsInput(input_points, input_triangles);
			arrangements.setTriMeshAsOutput(result_points, result_triangles);

			if (set_parameter)
			{
				arrangements.setConfig(arr_config);
			}

			OMC::MeshArrangements_Stats &stats = arrangements.stats();

			auto start = OMC::Logger::elapse_reset();

			arrangements.meshArrangements(false, false);

			double total_time = OMC::Logger::elapsed(start).count();
			std::cout << total_time << " s\n";

			log_file << std::fixed;
			log_file << iter->path().filename().string();
			log_file << "," << stats.pp_elapsed << "," << stats.tree_elapsed << ","
			         << stats.di_elapsed << "," << stats.cn_elapsed << ","
			         << stats.ci_elapsed << "," << stats.tr_elapsed << ","
			         << total_time << std::endl;
		}
	}
	log_file.close();
}
