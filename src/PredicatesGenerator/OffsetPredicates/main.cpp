#include "converter.h"

#include <iostream>
#include <format>

int main(int argc, char *argv[])
{
	// insufficient parameters
	if (argc < 2)
		help_info();

	// try open the file
	std::ifstream file(argv[1]);
	if (!file.is_open())
		error("Error: Can't open file", 0);

	// parse program parameters
	bool append                    = false;
	bool output_filtered_function  = false;
	bool output_interval_function  = false;
	bool output_exact_function     = false;
	bool output_expansion_function = false;

	for (int i = 2; i < argc; i++)
	{
		if (std::string(argv[i]) == "-a")
			append = true;
		else if (std::string(argv[i]) == "-ft")
			output_filtered_function = true;
		else if (std::string(argv[i]) == "-it")
			output_interval_function = true;
		else if (std::string(argv[i]) == "-et")
			output_exact_function = true;
		else if (std::string(argv[i]) == "-ex")
			output_expansion_function = true;
		else
		{
			std::cerr << std::format("unrecognized parameter {}\n", argv[i]);
			help_info();
		}
	}

	_controlfp(_RC_UP, _MCW_RC);

	Predicate predicate(append, output_filtered_function,
	                    output_interval_function, output_exact_function,
	                    output_expansion_function);

	std::string fullname(argv[1]);
	size_t      sls_idx = fullname.find_last_of("/\\");
	if (sls_idx == std::string::npos)
		sls_idx = size_t(-1);
	size_t      ext_idx   = fullname.find_last_of(".");
	std::string func_name = fullname.substr(sls_idx + 1, ext_idx - sls_idx - 1);

	if (func_name.substr(0, 6) == "lambda")
		predicate.is_lambda = true;

	int ln = 1;

	while (predicate.parseLine(file, ln))
		ln++;

	// Print all
	predicate.produceAllCode(func_name);

	predicate.printErrorBounds();

	_controlfp(_RC_NEAR, _MCW_RC);

	return 0;
}
