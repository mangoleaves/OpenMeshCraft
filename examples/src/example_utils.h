#pragma once

#include "boost/filesystem.hpp"

#include <fstream>
#include <stack>
#include <string>
#include <vector>

inline bool make_file_writable(const std::string &filename)
{
	namespace bfs = boost::filesystem;

	bfs::path filepath(filename);

	if (bfs::is_regular_file(filepath))
		return true;

	bfs::path parent_path = filepath.parent_path();
	if (bfs::is_directory(parent_path))
		return true;

	if (!bfs::create_directories(parent_path))
		throw std::logic_error("can't write file " + filename);

	return false;
}

inline bool make_dir_writable(const std::string &dirname)
{
	namespace bfs = boost::filesystem;

	bfs::path dirpath(dirname);

	if (bfs::is_directory(dirpath))
		return true;

	if (!bfs::create_directories(dirpath))
		throw std::logic_error("can't write dir " + dirname);

	return false;
}

#define OUTPUT_DIRECTORY(TEST_SUIT_NAME, TEST_NAME)           \
	std::string outdir =                                        \
	  "./data/test_output/" #TEST_SUIT_NAME "/" #TEST_NAME "/"; \
	make_dir_writable(outdir)
