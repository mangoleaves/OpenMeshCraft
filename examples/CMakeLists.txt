cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0074 NEW)
set(CMAKE_CXX_STANDARD 20)

# set the project name
project(OpenMeshCraft-Example)

# find Boost
# ##############################################################################
# We only use filesystem in example_utils.
set(Boost_NO_WARN_NEW_VERSIONS ON)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS filesystem)

# Configure examples
# ##############################################################################

# search all source files in current directory
aux_source_directory("${CMAKE_CURRENT_SOURCE_DIR}/src" Example_SRC)

# add example
add_executable(${PROJECT_NAME} ${Example_SRC})

target_include_directories(${PROJECT_NAME}
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCraft Boost::headers
                                             Boost::filesystem)
