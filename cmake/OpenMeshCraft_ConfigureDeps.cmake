# ##############################################################################
#
# What does this script do?
#
# It configures external libs for a target of OpenMeshCraft.
#
# How to use this script?
#
# 1. Set OMC_CONFIG_TARGET before including this script. The OMC_CONFIG_TARGET
#    will be used in target_include_directories and other similar commands.
# 2. Set specific variables for each dependency before including the script. See
#    details in the description of each dependency.
#
# Dependency list:
#
# * Boost
#
# ##############################################################################

# ##############################################################################
# Boost
#
# How to use:
#
# * Set OMC_CONFIGDEP_BOOST to ON to configure boost as a dependency.
# * Set OMC_CONFIGDEP_BOOST_COMPONENTS to find specified components, e.g.,
#   filesystem.
# * It calls find_package to find boost, users may provide hints to
#   find_package.
# * If boost is found, boost::headers and required components will be linked to
#   OMC_CONFIG_TARGET.
#
# ##############################################################################

if(DEFINED OMC_CONFIGDEP_BOOST AND OMC_CONFIGDEP_BOOST)

  # set options for boost
  set(Boost_DEBUG OFF)
  set(Boost_NO_WARN_NEW_VERSIONS ON)
  set(Boost_USE_STATIC_LIBS ON)

  if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(Boost_ARCHITECTURE "-x64")
  endif()

  # At least 1.78
  find_package(Boost 1.78.0 REQUIRED
               COMPONENTS ${OMC_CONFIGDEP_BOOST_COMPONENTS})

  # check if boost is found and report
  if(Boost_FOUND)
    if(OMC_CMAKE_VERBOSE OR OMC_CMAKE_DEBUG)
      message(
        STATUS "[OpenMeshCraft] Boost include dirs: ${Boost_INCLUDE_DIRS}")
      message(STATUS "[OpenMeshCraft] Boost libraries: ${Boost_LIBRARIES}")
    endif()
  else()
    message(FATAL_ERROR "[OpenMeshCraft] Can't find Boost.")
  endif()

  # link boost targets
  target_link_libraries(${OMC_CONFIG_TARGET} PUBLIC Boost::headers)
  foreach(_COMPONENT ${OMC_CONFIGDEP_BOOST_COMPONENTS})
    target_link_libraries(${OMC_CONFIG_TARGET} PUBLIC Boost::${_COMPONENT})
  endforeach()

endif()
