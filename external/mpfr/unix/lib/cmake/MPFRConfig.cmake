if(NOT TARGET MPFR)

  set(MPFR_INCLUDE_DIR
      "${CMAKE_CURRENT_LIST_DIR}/../../include"
      CACHE INTERNAL "Path to MPFR include directory.")
  set(MPFR_LIBRARIES_DIR "${CMAKE_CURRENT_LIST_DIR}/..")

  set(MPFR_LIBRARIES ${MPFR_LIBRARIES_DIR}/libmpfr.so.6.1.0)

  set(MPFR_SONAME libmpfr.so.6)

  message(STATUS "MPFR include dir: ${MPFR_INCLUDE_DIR}")
  message(STATUS "MPFR libraries: ${MPFR_LIBRARIES}")

  message(STATUS "Add MPFR as imported shared library.")
  add_library(MPFR SHARED IMPORTED GLOBAL)
  set_target_properties(
    MPFR
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR}"
               IMPORTED_LOCATION "${MPFR_LIBRARIES}"
               IMPORTED_SONAME "${MPFR_SONAME}")
endif()
