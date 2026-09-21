#pragma once

#include <algorithm>
#include <numeric>

#if defined(__has_include)
  #if __has_include(<execution>)
    #include <execution>
  #endif
#elif defined(_MSC_VER)
  #include <execution>
#endif

#if defined(__cpp_lib_parallel_algorithm) && \
  __cpp_lib_parallel_algorithm >= 201603L
  #define OMC_IF_EXEC(...) __VA_ARGS__,
#else
  #define OMC_IF_EXEC(...)
#endif
