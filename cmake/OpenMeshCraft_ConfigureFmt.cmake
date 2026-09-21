if(NOT DEFINED OMC_USE_STD_FORMAT)
  if(MSVC)
    set(OMC_USE_STD_FORMAT ON)
  else()
    set(OMC_USE_STD_FORMAT OFF)
  endif()
endif()

set(_OMC_HAS_STD_FORMAT OFF)
if(OMC_USE_STD_FORMAT)
  include(CheckCXXSourceCompiles)
  set(_OMC_SAVED_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  if(MSVC)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} /std:c++20")
  else()
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++20")
  endif()
  check_cxx_source_compiles(
    "
    #include <format>
    int main() {
      (void)std::format(\"{}\", 1);
      return 0;
    }
    "
    OMC_STD_FORMAT_WORKS)
  set(CMAKE_REQUIRED_FLAGS "${_OMC_SAVED_REQUIRED_FLAGS}")
  if(OMC_STD_FORMAT_WORKS)
    set(_OMC_HAS_STD_FORMAT ON)
  endif()
endif()

if(_OMC_HAS_STD_FORMAT)
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC OMC_HAS_STD_FORMAT)
  return()
endif()

if(DEFINED OMC_FMT_SOURCE_DIR AND
   EXISTS "${OMC_FMT_SOURCE_DIR}/include/fmt/format.h")
  set(_OMC_FMT_INCLUDE "${OMC_FMT_SOURCE_DIR}/include")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../external/fmt/include/fmt/format.h")
  set(_OMC_FMT_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../../external/fmt/include")
else()
  include(FetchContent)
  FetchContent_Declare(
    omc_fmt
    URL https://github.com/fmtlib/fmt/archive/refs/tags/10.2.1.zip
    URL_HASH
      SHA256=D368F9C39A33A3AEF800F5BE372EC1DF1C12AD57ADA1F60ADC62F24C0E348469)
  FetchContent_GetProperties(omc_fmt)
  if(NOT omc_fmt_POPULATED)
    FetchContent_Populate(omc_fmt)
  endif()
  set(_OMC_FMT_INCLUDE "${omc_fmt_SOURCE_DIR}/include")
endif()

target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC FMT_HEADER_ONLY
                                                       FMT_CONSTEVAL=)
target_include_directories(
  ${OMC_CONFIG_TARGET}
  PUBLIC $<BUILD_INTERFACE:${_OMC_FMT_INCLUDE}>)
