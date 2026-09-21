#pragma once

#ifdef OMC_HAS_STD_FORMAT
  #include <format>

namespace OMC {
using std::format;
}
#else
  #ifndef FMT_CONSTEVAL
    #define FMT_CONSTEVAL
  #endif
  #ifndef FMT_HEADER_ONLY
    #define FMT_HEADER_ONLY
  #endif
  #include <fmt/format.h>

namespace OMC {
using fmt::format;
}
#endif
