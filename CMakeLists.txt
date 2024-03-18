# ##############################################################################
# cmake settings
# ##############################################################################

cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0074 NEW)

project(OpenMeshCraft VERSION 1.0.0)

# ##############################################################################
# build and deploy settings
# ##############################################################################

if(NOT DEFINED OMC_MASTER_PROJECT)
  set(OMC_MASTER_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(OMC_MASTER_PROJECT ON)
  endif()
endif()

# cmake-format: off
option(OMC_CMAKE_ENABLE_EXCEPTION "Enable exception in OpenMeshCraft" ON)
option(OMC_CMAKE_ENABLE_ASSERT "Enable assertion in OpenMeshCraft" ON)
option(OMC_CMAKE_ENABLE_EXPENSIVE_ASSERT "Enable expensive assertion in OpenMeshCraft" OFF)
option(OMC_CMAKE_ENABLE_SSE2 "Enable SSE2 instruction set in OpenMeshCraft" ON)
option(OMC_CMAKE_ENABLE_AVX "Enable AVX instruction set in OpenMeshCraft" ON)
option(OMC_CMAKE_ENABLE_DEBUG_INFO "Generate debug information for OpenMeshCraft" ON)

option(OMC_CMAKE_LOG_TRACE "Log runtime trace messages." OFF)
option(OMC_CMAKE_LOG_DEBUG "Log runtime debug messages." OFF)
option(OMC_CMAKE_LOG_INFO "Log runtime info messages." ON)
option(OMC_CMAKE_LOG_WARN "Log runtime warn messages." ON)
option(OMC_CMAKE_LOG_FATAL "Log runtime fatal messages." ON)

option(OMC_CMAKE_PRED_PROFILE "Enable profiling predicates in OpenMeshCraft." OFF)

option(OMC_BUILD_TEST "Build test project of OpenMeshCraft" ${OMC_MASTER_PROJECT})
option(OMC_BUILD_EXAMPLE "Build example project of OpenMeshCraft" ${OMC_MASTER_PROJECT})
# cmake-format: on

# ##############################################################################
# build and deploy settings
# ##############################################################################

macro(set_dir_with_configs _DIR_NAME _DIR_PATH)

  set(${_DIR_NAME} ${_DIR_PATH})

  set(${_DIR_NAME}_RELEASE ${_DIR_PATH}/Release)
  set(${_DIR_NAME}_DEBUG ${_DIR_PATH}/Debug)

endmacro()

# global directory of static libraries
set_dir_with_configs(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# global directory of shared libararies
set_dir_with_configs(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# global directory of executable files
set_dir_with_configs(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# global install prefix
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)

# global module path.
#
# if user sets CMAKE_MODULE_PATH in command line, CMAKE_MODULE_PATH won't change
# in this project.
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
                      "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# global prefix path.
#
# if user sets CMAKE_PREFIX_PATH in command line, CMAKE_PREFIX_PATH won't change
# in this project.
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}"
                      "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ##############################################################################
# compile options and definitions for third party libs
# ##############################################################################

# block useless warnings
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # The file contains a character that cannot be represented in the current code
  # page (number). Save the file in Unicode format to prevent data loss.
  add_compile_options("/wd 4819")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # pass
endif()

# ##############################################################################
# add external libraries
# ##############################################################################

add_subdirectory(external)

# ##############################################################################
# global compile options and definitions
# ##############################################################################

# global compile definitions
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_compile_definitions(_USE_MATH_DEFINES)
  add_compile_definitions(NOMINMAX)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # pass
endif()

# ##############################################################################
# add subdirectory/sub-projects
# ##############################################################################

add_subdirectory(src)

if(OMC_BUILD_TEST)
  add_subdirectory(test)
endif()

if(OMC_BUILD_EXAMPLE)
  add_subdirectory(examples)
endif()
