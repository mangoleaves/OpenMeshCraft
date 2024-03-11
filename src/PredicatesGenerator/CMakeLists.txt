cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0074 NEW)
set(CMAKE_CXX_STANDARD 20)

# set the project name
project(PredicatesGenerator)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES converter.cpp)

# add the library
add_executable(${PROJECT_NAME} converter.cpp converter.h main.cpp)

# include directory
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC)

set(OMC_CONFIG_TARGET ${PROJECT_NAME})
include(OpenMeshCraft_ConfigureCompiler)