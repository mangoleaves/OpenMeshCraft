#pragma once

#include "Macros.h"

#include <string>
#include <vector>

namespace OMC {

std::vector<std::string> split_string(const std::string &str, char delim);

void trim_string(std::string &str);

bool starts_with(const std::string &big_str, const std::string &small_str);

bool ends_with(const std::string &big_str, const std::string &small_str);

std::string replace_first(const std::string &str, const std::string &orig,
                          const std::string &rep);

std::string replace_last(const std::string &str, const std::string &orig,
                         const std::string &rep);

std::string replace_all(const std::string &str, const std::string &orig,
                        const std::string &rep);

} // namespace OMC