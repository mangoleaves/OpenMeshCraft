# ##############################################################################
# cmake settings
# ##############################################################################

cmake_minimum_required(VERSION 3.11)
set(CMAKE_CXX_STANDARD 20)

# ##############################################################################
# add subdirectory/sub-projects
# ##############################################################################

add_subdirectory(OpenMeshCraft)

if(OMC_MASTER_PROJECT AND OMC_BUILD_PREDICATES_GENERATOR)
  add_subdirectory(PredicatesGenerator)
endif()
