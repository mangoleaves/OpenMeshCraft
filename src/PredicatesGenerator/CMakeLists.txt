cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0074 NEW)
set(CMAKE_CXX_STANDARD 20)

# set the project name
project(PredicatesGenerator)

# #############################################################################
# Indirect predicates generator
# #############################################################################

file(GLOB IND_SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/IndirectPredicates/*.h")
file(GLOB IND_SRC_CPP "${CMAKE_CURRENT_SOURCE_DIR}/IndirectPredicates/*.cpp")

# add the library
add_executable(Indirect${PROJECT_NAME} ${IND_SRC_H} ${IND_SRC_CPP})

# include directory
target_include_directories(Indirect${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link library
target_link_libraries(Indirect${PROJECT_NAME} PUBLIC)

# configure compiler
set(OMC_CONFIG_TARGET Indirect${PROJECT_NAME})
include(OpenMeshCraft_ConfigureCompiler)

# #############################################################################
# Offset predicates generator
# #############################################################################

file(GLOB OFF_SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/OffsetPredicates/*.h")
file(GLOB OFF_SRC_CPP "${CMAKE_CURRENT_SOURCE_DIR}/OffsetPredicates/*.cpp")

# add the library
add_executable(Offset${PROJECT_NAME} ${OFF_SRC_H} ${OFF_SRC_CPP})

# include directory
target_include_directories(Offset${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link library
target_link_libraries(Offset${PROJECT_NAME} PUBLIC)

# configure compiler
set(OMC_CONFIG_TARGET Offset${PROJECT_NAME})
include(OpenMeshCraft_ConfigureCompiler)
