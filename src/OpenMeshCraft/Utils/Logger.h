#pragma once

#include "Macros.h"

#include <chrono>
#include <cstdio>
#include <format>
#include <iostream>
#include <string>

namespace OMC {

/*****************************************************************************/
/* Logger ********************************************************************/
/*****************************************************************************/

class Logger
{
public:
	enum class Level
	{
		TRACE = 1,
		DEBUG = 2,
		INFO  = 3,
		WARN  = 4,
		FATAL = 5
	};

public:
	static void trace(const std::string& msg);
	static void debug(const std::string& msg);
	static void info(const std::string& msg);
	static void warn(const std::string& msg);
	static void fatal(const std::string& msg);

	/// @brief return the current time point.
	static OMC_NODISCARD std::chrono::steady_clock::time_point elapse_reset();

	/// @brief return the eplased time between current and last time point.
	static OMC_NODISCARD std::chrono::duration<double>
	elapsed(const std::chrono::steady_clock::time_point &last);
};

/*****************************************************************************/
/* Memory ********************************************************************/
/*****************************************************************************/

double getPeakMegabytesUsed();
double getMegabytesUsed();

} // namespace OMC