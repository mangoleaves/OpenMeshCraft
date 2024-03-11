#include "OpenMeshCraft/Utils/Logger.h"
#include "test_utils.h"

boost::property_tree::ptree omc_test_config;

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);

	for (int i = 1; i < argc; i++)
	{
		std::string param = std::string(argv[i]);
		if (OMC::starts_with(param, "--config="))
		{
			std::string config_path = param.substr(param.find_first_of('=') + 1);
			if (boost::filesystem::is_regular_file(
			      boost::filesystem::path(config_path)))
			{
				boost::property_tree::read_json(config_path,
				                                omc_test_config);
			}
			else
			{
				OMC_ASSERT(false, "test config file error.");
			}
			break;
		}
	}

	return RUN_ALL_TESTS();
}