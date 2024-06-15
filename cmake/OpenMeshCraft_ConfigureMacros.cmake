# enable runtime Exception
if(OMC_CMAKE_ENABLE_EXCEPTION)
  message(STATUS "[OpenMeshCraft] Enable runtime exception")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_ENABLE_EXCEPTION)
else()
  message(STATUS "[OpenMeshCraft] Disable runtime exception")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_DISABLE_EXCEPTION)
endif()

# enable expensive assertion
if(OMC_CMAKE_ENABLE_ASSERT)
  message(STATUS "[OpenMeshCraft] Enable assert")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_ENABLE_ASSERT)
else()
  message(STATUS "[OpenMeshCraft] Disable assert")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_DISABLE_ASSERT)
endif()

# enable expensive assertion
if(OMC_CMAKE_ENABLE_EXPENSIVE_ASSERT)
  message(STATUS "[OpenMeshCraft] Enable expensive assert")
  target_compile_definitions(${OMC_CONFIG_TARGET}
                             PUBLIC OMC_ENABLE_EXPENSIVE_ASSERT)
else()
  message(STATUS "[OpenMeshCraft] Disable expensive assert")
  target_compile_definitions(${OMC_CONFIG_TARGET}
                             PUBLIC OMC_DISABLE_EXPENSIVE_ASSERT)
endif()

# enable SSE2 vectorization set
if(OMC_CMAKE_ENABLE_SSE2)
  message(STATUS "[OpenMeshCraft] Enable SSE2")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_ENABLE_SSE2)
else()
  message(STATUS "[OpenMeshCraft] Disable SSE2")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_DISABLE_SSE2)
endif()

# enable AVX vectorization set
if(OMC_CMAKE_ENABLE_AVX)
  message(STATUS "[OpenMeshCraft] Enable AVX")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_ENABLE_AVX)
else()
  message(STATUS "[OpenMeshCraft] Disable AVX")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_DISABLE_AVX)
endif()

# enable AVX2 vectorization set
if(OMC_CMAKE_ENABLE_AVX2)
  message(STATUS "[OpenMeshCraft] Enable AVX2")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_ENABLE_AVX2)
else()
  message(STATUS "[OpenMeshCraft] Disable AVX2")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_DISABLE_AVX2)
endif()

# enable AVX2 vectorization set
if(OMC_CMAKE_ENABLE_FMA)
  message(STATUS "[OpenMeshCraft] Enable FMA")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_ENABLE_FMA)
else()
  message(STATUS "[OpenMeshCraft] Disable FMA")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_DISABLE_FMA)
endif()

# enable log message at some level
if(OMC_CMAKE_LOG_TRACE)
  message(STATUS "[OpenMeshCraft] Log runtime trace message.")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_LOG_TRACE)
endif()
if(OMC_CMAKE_LOG_DEBUG)
  message(STATUS "[OpenMeshCraft] Log runtime debug message.")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_LOG_DEBUG)
endif()
if(OMC_CMAKE_LOG_INFO)
  message(STATUS "[OpenMeshCraft] Log runtime information message.")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_LOG_INFO)
endif()
if(OMC_CMAKE_LOG_WARN)
  message(STATUS "[OpenMeshCraft] Log runtime warning message.")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_LOG_WARN)
endif()
if(OMC_CMAKE_LOG_FATAL)
  message(STATUS "[OpenMeshCraft] Log runtime fatal error message.")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_LOG_FATAL)
endif()

# enable profiling predicates
if(OMC_CMAKE_PRED_PROFILE)
  message(STATUS "[OpenMeshCraft] Enable profiling predicates")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_PRED_PROFILE)
else()
  # message(STATUS "[OpenMeshCraft] Disable profiling predicates")
endif()

# enable profiling intersections
if(OMC_CMAKE_INTER_PROFILE)
  message(STATUS "[OpenMeshCraft] Enable profiling intersections")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_INTER_PROFILE)
else()
  # message(STATUS "[OpenMeshCraft] Disable profiling intersections")
endif()

# enable profiling arrangements
if(OMC_CMAKE_ARR_PROFILE)
  message(STATUS "[OpenMeshCraft] Enable profiling arrangements")
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_ARR_PROFILE)
else()
  # message(STATUS "[OpenMeshCraft] Disable profiling arrangements")
endif()

# generate debug information
if(OMC_CMAKE_ENABLE_DEBUG_INFO)

  message(STATUS "[OpenMeshCraft] Enable DebugInfo")

  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -g)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # generate debug information
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /Zi)
    # turn on Debug information.
    target_link_options(${OMC_CONFIG_TARGET} PUBLIC /Debug)
  endif()

endif()
