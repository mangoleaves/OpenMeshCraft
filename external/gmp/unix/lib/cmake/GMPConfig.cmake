if(NOT TARGET GMP)

  set(GMP_INCLUDE_DIR
      "${CMAKE_CURRENT_LIST_DIR}/../../include"
      CACHE INTERNAL "Path to GMP include directory.")
  set(GMP_LIBRARIES_DIR "${CMAKE_CURRENT_LIST_DIR}/..")

  set(GMP_LIBRARIES ${GMP_LIBRARIES_DIR}/libgmp.so.10.4.1)

  set(GMP_SONAME libgmp.so.10)

  message(STATUS "GMP include dir: ${GMP_INCLUDE_DIR}")
  message(STATUS "GMP libraries: ${GMP_LIBRARIES}")

  message(STATUS "Add GMP as imported shared library.")
  add_library(GMP SHARED IMPORTED GLOBAL)
  set_target_properties(
    GMP
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
               IMPORTED_LOCATION "${GMP_LIBRARIES}"
               IMPORTED_SONAME "${GMP_SONAME}")
endif()
