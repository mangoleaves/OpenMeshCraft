#pragma once

/* OpenMeshCraft */

#include "OpenMeshCraft/Mesh/IO/OBJReader.h"
#include "OpenMeshCraft/Mesh/IO/OBJWriter.h"
#include "OpenMeshCraft/Mesh/IO/STLReader.h"
#include "OpenMeshCraft/Mesh/IO/STLWriter.h"

#include "OpenMeshCraft/Mesh/TriSoup.h"

#include "OpenMeshCraft/Utils/StringUtils.h"

/* External libs */

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"

#include "gtest/gtest.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include <filesystem>
#include <fstream>
#include <stack>
#include <string>
#include <vector>

extern boost::property_tree::ptree omc_test_config;

/**
 * @brief Ensures the file's directory is writable.
 *
 * @param filename The name of the file.
 * @return true if the file's directory is writable.
 * @return false if the file's directory is not writable.
 */
inline bool make_file_writable(const std::string &filename)
{
	namespace fs = std::filesystem;

	fs::path filepath(filename);

	if (fs::is_regular_file(filepath))
		return true;

	fs::path parent_path = filepath.parent_path();
	if (fs::is_directory(parent_path))
		return true;

	if (!fs::create_directories(parent_path))
		throw std::logic_error("can't write file " + filename);

	return false;
}

/**
 * @brief Ensures the directory is writable.
 *
 * @param dirname The name of the directory.
 * @return true if the directory is writable.
 * @return false if the directory is not writable.
 */
inline bool make_dir_writable(const std::string &dirname)
{
	namespace fs = std::filesystem;

	fs::path dirpath(dirname);

	if (fs::is_directory(dirpath))
		return true;

	if (!fs::create_directories(dirpath))
		throw std::logic_error("can't write dir " + dirname);

	return false;
}

using index_t       = OMC::index_t;
using TriSoupTraits = OMC::TriSoupTraits;

using Points    = typename TriSoupTraits::Points;
using Triangles = typename TriSoupTraits::Triangles;

// Define reader and writer
using OBJReader = OMC::OBJReader<TriSoupTraits>;
using OBJWriter = OMC::OBJWriter<TriSoupTraits>;
using STLReader = OMC::STLReader<TriSoupTraits>;
using STLWriter = OMC::STLWriter<TriSoupTraits>;

using IOOptions = OMC::IOOptions;

/**
 * @brief Reads a mesh from a file.
 *
 * @param filename The name of the file.
 * @param points The points of the mesh.
 * @param triangles The triangles of the mesh.
 * @param io_options The IO options.
 */
inline void read_mesh(const std::string &filename, Points &points,
                      Triangles &triangles, IOOptions &io_options)
{
	if (OMC::ends_with(filename, ".obj") || OMC::ends_with(filename, ".OBJ"))
	{
		OBJReader reader;
		reader.read(filename, io_options);

		points    = std::move(reader.m_points);
		triangles = std::move(reader.m_triangles);
	}
	else if (OMC::ends_with(filename, ".stl") || OMC::ends_with(filename, ".STL"))
	{
		STLReader reader;
		reader.read(filename, io_options);

		points    = std::move(reader.m_points);
		triangles = std::move(reader.m_triangles);
	}
	else
	{
		throw std::runtime_error("unsupported file type.");
	}
}

/**
 * @brief Writes a mesh to a file.
 *
 * @param filename The name of the file.
 * @param points The points of the mesh.
 * @param triangles The triangles of the mesh.
 * @param io_options The IO options.
 */
inline void write_mesh(const std::string &filename, const Points &points,
                       const Triangles &triangles, IOOptions &io_options)
{
	if (OMC::ends_with(filename, ".obj") || OMC::ends_with(filename, ".OBJ"))
	{
		OBJWriter writer;
		writer.m_points    = std::move(points);
		writer.m_triangles = std::move(triangles);
		writer.write(filename, io_options, DBL_DIG);
	}
	else if (OMC::ends_with(filename, ".stl") || OMC::ends_with(filename, ".STL"))
	{
		/* Precision is not enough, change to obj */
		std::string new_fn = OMC::replace_last(filename, ".stl", ".obj");
		OBJWriter   writer;
		writer.m_points       = std::move(points);
		writer.m_triangles    = std::move(triangles);
		io_options.stl_binary = true;
		writer.write(new_fn, io_options, DBL_DIG);
	}
	else
	{
		throw std::runtime_error("unsupported file type.");
	}
};

#define TEST_OUTPUT_DIRECTORY(TEST_SUIT_NAME, TEST_NAME)                     \
	std::string outdir = "./data/test_output/" #TEST_SUIT_NAME "/" #TEST_NAME; \
	make_dir_writable(outdir)

#define TEST_GET_CONFIG(TEST_SUIT_NAME, TEST_NAME) \
	boost::property_tree::ptree &config =            \
	  omc_test_config.get_child(#TEST_SUIT_NAME).get_child(#TEST_NAME)