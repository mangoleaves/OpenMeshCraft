# ##############################################################################
# cmake settings
# ##############################################################################

cmake_minimum_required(VERSION 3.11)
set(CMAKE_CXX_STANDARD 20)

# ##############################################################################
# add subdirectory/sub-projects
# ##############################################################################

add_subdirectory(OpenMeshCraft)

add_subdirectory(PredicatesGenerator EXCLUDE_FROM_ALL)
