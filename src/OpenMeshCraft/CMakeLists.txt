cmake_minimum_required(VERSION 3.11)
# set CMP0074 to not ignore XXX_ROOT variables.
cmake_policy(SET CMP0074 NEW)
set(CMAKE_CXX_STANDARD 20)

# set the project name
project(
  OpenMeshCraft
  VERSION 1.0.0
  LANGUAGES C CXX)

# ##############################################################################
# Configure OpenMeshCraft
# ##############################################################################

# get source files
# ##############################################################################

file(GLOB_RECURSE SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB_RECURSE SRC_CPP "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB_RECURSE SRC_INL "${CMAKE_CURRENT_SOURCE_DIR}/*.inl")

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SRC_H} ${SRC_INL}
                                                    ${SRC_CPP})

# add the library
# ##############################################################################

add_library(${PROJECT_NAME} STATIC ${SRC_H} ${SRC_INL} ${SRC_CPP})

# include directory
# ##############################################################################

# BUILD_INTERFACE is used for building, INSTALL_INTERFACE is used for
# installing.
target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                         $<INSTALL_INTERFACE:include>)

# link libraries
# ##############################################################################

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC # Eigen
         Eigen3::Eigen
         # GMP and MPFR
         GMP
         MPFR
         # parallel hashmap
         phmap
         # shewchuk predicates
         shewchuk_predicates
         # oneTBB
         tbb)

set(OMC_CONFIG_TARGET ${PROJECT_NAME})
set(OMC_CONFIGDEP_BOOST ON)
include(OpenMeshCraft_ConfigureDeps)

# compile definitions and options
# ##############################################################################
target_compile_definitions(
  ${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:OMC_HAS_IMPL>
  PUBLIC $<INSTALL_INTERFACE:OMC_HAS_IMPL>)

include(OpenMeshCraft_ConfigureMacros)

include(OpenMeshCraft_ConfigureCompiler)
