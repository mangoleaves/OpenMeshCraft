#include "OpenMeshCraft/Arrangements/Utils.h"
#include "OpenMeshCraft/Geometry/Intersection/IntersectionUtils.h"
#include "OpenMeshCraft/Geometry/Predicates/IndirectDefinitions.h"
#include "OpenMeshCraft/Utils/Logger.h"

#include "test_utils.h"

boost::property_tree::ptree omc_test_config;

int main(int argc, char **argv)
{
	OMC_PRED_PROFILE_INIT;
	OMC_INTER_PROFILE_INIT;
	OMC_ARR_PROFILE_INIT;

	::testing::InitGoogleTest(&argc, argv);

	for (int i = 1; i < argc; i++)
	{
		std::string param = std::string(argv[i]);
		if (OMC::starts_with(param, "--config="))
		{
			std::string config_path = param.substr(param.find_first_of('=') + 1);
			if (std::filesystem::is_regular_file(std::filesystem::path(config_path)))
			{
				boost::property_tree::read_json(config_path, omc_test_config);
			}
			else
			{
				OMC_ASSERT(false, "test config file error.");
			}
			break;
		}
	}

	int ret = RUN_ALL_TESTS();

	OMC_PRED_PROFILE_PRINT;
	OMC_INTER_PROFILE_PRINT;
	OMC_ARR_PROFILE_PRINT;

	return ret;
}