cmake_minimum_required(VERSION 3.11)
cmake_policy(SET CMP0074 NEW)
set(CMAKE_CXX_STANDARD 20)

# set the project name
project(OpenMeshCraft-Test)

# ##############################################################################
# Configure test library
# ##############################################################################

# search all source files in current directory
# ##############################################################################
# cmake-format: off
file(GLOB_RECURSE TEST_Arrangements   "${CMAKE_CURRENT_SOURCE_DIR}/Arrangements/*")
file(GLOB_RECURSE TEST_Boolean        "${CMAKE_CURRENT_SOURCE_DIR}/Boolean/*")
file(GLOB_RECURSE TEST_Geometry       "${CMAKE_CURRENT_SOURCE_DIR}/Geometry/*")
file(GLOB_RECURSE TEST_Mesh           "${CMAKE_CURRENT_SOURCE_DIR}/Mesh/*")
file(GLOB_RECURSE TEST_NumberTypes    "${CMAKE_CURRENT_SOURCE_DIR}/NumberTypes/*")
file(GLOB_RECURSE TEST_Utils          "${CMAKE_CURRENT_SOURCE_DIR}/Utils/*")
# cmake-format: on

source_group(
  TREE ${CMAKE_CURRENT_SOURCE_DIR}
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/test_utils.h
        ${TEST_Arrangements}
        ${TEST_Boolean}
        ${TEST_Geometry}
        ${TEST_Mesh}
        ${TEST_NumberTypes}
        ${TEST_Utils})

# add test
# ##############################################################################
add_executable(
  ${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_utils.h
  ${TEST_Arrangements}
  ${TEST_Boolean}
  ${TEST_Geometry}
  ${TEST_Mesh}
  ${TEST_NumberTypes}
  ${TEST_Utils})

# include directories
# ##############################################################################
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link libraries
# ##############################################################################
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCraft CGAL::CGAL
                                             gtest_main)

# general configuration of current target
# ##############################################################################

set(OMC_CONFIG_TARGET ${PROJECT_NAME})

set(OMC_CONFIGDEP_BOOST ON)
set(OMC_CONFIGDEP_BOOST_COMPONENTS filesystem)
# finally configure...
include(OpenMeshCraft_ConfigureDeps)

# compile definitions and options
# ##############################################################################

# pass

# deploy thirdparty libraries after building
# ##############################################################################
# It is necessary only on Windows to deploy runtime libs.

if(CMAKE_HOST_WIN32)

  # ## abseil-cpp
  get_target_property(_IS_ABSL_IMPORTED abseil_dll IMPORTED)
  if(_IS_ABSL_IMPORTED)
    # cmake-format: off
    install(
      IMPORTED_RUNTIME_ARTIFACTS abseil_dll
      EXCLUDE_FROM_ALL
      COMPONENT ${PROJECT_NAME}
      RUNTIME DESTINATION $<TARGET_FILE_DIR:${PROJECT_NAME}>)
    # cmake-format: on
  else()
    install(
      TARGETS abseil_dll
      EXCLUDE_FROM_ALL
      COMPONENT ${PROJECT_NAME}
      RUNTIME DESTINATION $<TARGET_FILE_DIR:${PROJECT_NAME}>)
  endif()

  # ## GMP and MPFR
  # cmake-format: off
  install(
    IMPORTED_RUNTIME_ARTIFACTS GMP MPFR
    EXCLUDE_FROM_ALL
    COMPONENT ${PROJECT_NAME}
    RUNTIME DESTINATION $<TARGET_FILE_DIR:${PROJECT_NAME}>)
  # cmake-format: on

  # ## TBB
  get_target_property(_IS_TBB_IMPORTED tbb IMPORTED)
  if(_IS_TBB_IMPORTED)
    # cmake-format: off
    install(
      IMPORTED_RUNTIME_ARTIFACTS tbb
      EXCLUDE_FROM_ALL
      COMPONENT ${PROJECT_NAME}
      RUNTIME DESTINATION $<TARGET_FILE_DIR:${PROJECT_NAME}>)
    # cmake-format: on
  else()
    install(
      TARGETS tbb
      EXCLUDE_FROM_ALL
      COMPONENT ${PROJECT_NAME}
      RUNTIME DESTINATION $<TARGET_FILE_DIR:${PROJECT_NAME}>)
  endif()
  unset(_IS_TBB_IMPORTED)

  add_custom_command(
    TARGET ${PROJECT_NAME}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR} --component
            ${PROJECT_NAME} --config $<CONFIG>)
endif()

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})

# ##############################################################################
# ##############################################################################
# ##############################################################################

# ##############################################################################
# Configure an executable to run dataset for arrangements.
# ##############################################################################

add_executable(OpenMeshCraft-Arrangements paper_arrangements.cpp)

# include directories
# ##############################################################################
target_include_directories(OpenMeshCraft-Arrangements
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# link libraries
# ##############################################################################
target_link_libraries(OpenMeshCraft-Arrangements PUBLIC OpenMeshCraft gtest_main)

# general configuration of current target
# ##############################################################################

set(OMC_CONFIG_TARGET OpenMeshCraft-Arrangements)

set(OMC_CONFIGDEP_BOOST ON)
set(OMC_CONFIGDEP_BOOST_COMPONENTS filesystem)
# finally configure...
include(OpenMeshCraft_ConfigureDeps)

add_custom_command(
  TARGET OpenMeshCraft-Arrangements
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR} --component
          ${PROJECT_NAME} --config $<CONFIG>)
