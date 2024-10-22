# policy settings #############################################################
cmake_policy(SET CMP0074 NEW) # use XXX_ROOT variable in find_package
cmake_policy(SET CMP0167 OLD) # use BOOST_ROOT in find_package(boost)

# GMP ##########################################################################
# ##############################################################################
# The target GMP can be linked to other targets.

# set dir to variable, it will be used by find_package()
if(CMAKE_HOST_WIN32)
  set(GMP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gmp/windows")
  list(APPEND CMAKE_PREFIX_PATH
       "${CMAKE_CURRENT_SOURCE_DIR}/gmp/windows/lib/cmake")
elseif(CMAKE_HOST_UNIX)
  set(GMP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gmp/unix")
  list(APPEND CMAKE_PREFIX_PATH
       "${CMAKE_CURRENT_SOURCE_DIR}/gmp/unix/lib/cmake")
else()
  message(FATAL_ERROR "[OpenMeshCraft] Unsupported host for GMP:"
                      "${CMAKE_HOST_SYSTEM_NAME}.")
endif()

find_package(GMP REQUIRED NO_MODULE)

# MPFR #########################################################################
# ##############################################################################
# Target MPFR can be linked to other targets.

# set dir to variable, it will be used by find_package()
if(CMAKE_HOST_WIN32)
  set(MPFR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/mpfr/windows")
  list(APPEND CMAKE_PREFIX_PATH
       "${CMAKE_CURRENT_SOURCE_DIR}/mpfr/windows/lib/cmake")
elseif(CMAKE_HOST_UNIX)
  set(MPFR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/mpfr/unix")
  list(APPEND CMAKE_PREFIX_PATH
       "${CMAKE_CURRENT_SOURCE_DIR}/mpfr/unix/lib/cmake")
else()
  message(FATAL_ERROR "[OpenMeshCraft] Unsupported host for MPFR:"
                      "${CMAKE_HOST_SYSTEM_NAME}.")
endif()

find_package(MPFR REQUIRED NO_MODULE)

# CGAL #########################################################################
# ##############################################################################
# The target CGAL::CGAL can be linked to other targets.

# CGAL_DIR will be used by find_pacakge(CGAL ...)
set(CGAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cgal)

# find boost before finding CGAL
set(Boost_DEBUG OFF)
set(Boost_NO_WARN_NEW_VERSIONS ON)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.78.0 REQUIRED)

# find CGAL. This procedure should be executed after finding boost.
set(CGAL_Boost_USE_STATIC_LIBS ON)
find_package(CGAL)
if(CGAL_FOUND)
  message(STATUS "[OpenMeshCraft] Found CGAL.")
else()
  message(FATAL_ERROR "[OpenMeshCraft] Can't find CGAL.")
endif()

# Eigen ########################################################################
# ##############################################################################
# Header-only library, set EIGEN3_INCLUDE_DIR for further find_package. We also
# add its FindEigen3.cmake to module path. If an installed shared Eigen3 lib is
# provided, add the lib/cmake to CMAKE_PREFIX_PATH.

set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)

find_package(Eigen3 3.4.0 REQUIRED)

# Google test ##################################################################
# ##############################################################################
# Target gtest_main can be linked to other targets.
#
# After linking it, run below commands:
#
# * include(GoogleTest)
# * gtest_discover_tests(YOUR_TARGET)

if(CMAKE_HOST_WIN32)
  set(gtest_force_shared_crt
      ON
      CACHE BOOL "")
endif()
set(BUILD_GMOCK OFF)
set(INSTALL_GTEST OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/googletest EXCLUDE_FROM_ALL)

# parallel-hashmap #############################################################
# ##############################################################################
# Target phmap can be linked to other targets.

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/parallel-hashmap EXCLUDE_FROM_ALL)

# shewchuk-predicates ##########################################################
# ##############################################################################
# A very simple C lib. Target shewchuk_predicates can be linked to other
# targets.

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/shewchuk-predicates
                 EXCLUDE_FROM_ALL)

# TBB ##########################################################################
# ##############################################################################
# Target tbb can be linked to other targets.

# cmake-format: off
set(TBB_TEST OFF CACHE BOOL "")
set(TBBMALLOC_BUILD OFF CACHE BOOL "")
# cmake-format: on

if(WIN32 AND "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # When compile tbb by gcc (mingw64) on Windows, compiler complains a lot of warnings.
  # It's ok to compile tbb by MSVC on Windows or by gcc on Unix* OS.
  set(TBB_WARNING_SUPPRESS -Wno-stringop-overflow -Wno-unused-value
                           -Wno-array-bounds)
endif()

# always build TBB to shared lib

if(NOT BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/oneTBB EXCLUDE_FROM_ALL)
  set(BUILD_SHARED_LIBS OFF)
else()
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/oneTBB EXCLUDE_FROM_ALL)
endif()

# prevent generating versioned .so files.
set_property(TARGET tbb PROPERTY VERSION)
set_property(TARGET tbb PROPERTY SOVERSION)

# Pass variabls to parent scope
set(CMAKE_PREFIX_PATH
    ${CMAKE_PREFIX_PATH}
    PARENT_SCOPE)
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    PARENT_SCOPE)
