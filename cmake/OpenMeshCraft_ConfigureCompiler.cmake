# Compiler-specific options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # -fno-unsafe-math-optimizations
  # ############################################################################
  # Allow optimizations for floating-point arithmetic that (a) assume that
  # arguments and results are valid and (b) may violate IEEE or ANSI standards.
  # When used at link time, it may include libraries or startup files that
  # change the default FPU control word or other similar optimizations.
  # This option is not turned on by any -O option besides -Ofast since it can
  # result in incorrect output for programs that depend on an exact implementation
  # of IEEE or ISO rules/specifications for math functions. It may, however,
  # yield faster code for programs that do not require the guarantees of these
  # specifications. Enables -fno-signed-zeros, -fno-trapping-math, -fassociative-math
  # and -freciprocal-math.
  # The default is -fno-unsafe-math-optimizations.
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -fno-unsafe-math-optimizations)

  # -frounding-math
  # ############################################################################
  # "Disable transformations and optimizations that assume default floating
  # point rounding behavior. This is round-to-zero for all floating point to
  # integer conversions, and round-to-nearest for all other arithmetic
  # truncations. This option should be specified for programs that change the FP
  # rounding mode dynamically, or that may be executed with a non-default
  # rounding mode. This option disables constant folding of floating point
  # expressions at compile-time (which may be affected by rounding mode) and
  # arithmetic transformations that are unsafe in the presence of sign-dependent
  # rounding modes."
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -frounding-math)

  # -fsignaling-nans
  # ############################################################################
  # Compile code assuming that IEEE signaling NaNs may generate user-visible
  # traps during floating-point operations. Setting this option disables
  # optimizations that may change the number of exceptions visible with
  # signaling NaNs. This option implies -ftrapping-math.  This option causes the
  # preprocessor macro __SUPPORT_SNAN__ to be defined.
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -fsignaling-nans)

  # handle warnings
  target_compile_options(
    ${OMC_CONFIG_TARGET}
    PUBLIC # raise up warning level
           -Wall
           -Wextra
           # block comment error
           -Wno-comment
           # pass pragmas only on Windows
           -Wno-unknown-pragmas
           # treat all warnings as errors
           -Werror)

  # handle vectorization
  if (OMC_CMAKE_ENABLE_AVX2)
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -mavx2)
  elseif(OMC_CMAKE_ENABLE_AVX)
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -mavx)
  elseif(OMC_CMAKE_ENABLE_SSE2)
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -msse2)
  endif()

  if (OMC_CMAKE_ENABLE_FMA)
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -mfma)
  endif()

  # reserve enough stack size
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC -Wl,-z,stacksize=8421376)

elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # /fp:strict
  # ############################################################################
  # the compiler preserves the source ordering and rounding properties of
  # floating-point code when it generates and optimizes object code for the
  # target machine, and observes the standard when handling special values. The
  # program may also safely access or modify the floating-point environment at
  # runtime.
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /fp:strict)

  # optimize speed and use intrinsic functions
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /Oi)

  # turn off annoying warnings
  target_compile_options(${OMC_CONFIG_TARGET}
                         PUBLIC "/D _CRT_SECURE_NO_WARNINGS")

  # turn on multiprocessor compile.
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /MP)

  # raise up warning level and treat all warnings as errors.
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /W4 /WX)

  # handle vectorization
  if (OMC_CMAKE_ENABLE_AVX2)
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /arch:AVX2)
  elseif(OMC_CMAKE_ENABLE_AVX)
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /arch:AVX)
  elseif(OMC_CMAKE_ENABLE_SSE2)
    target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /arch:SSE2)
  endif()

  # enable big obj file
  target_compile_options(${OMC_CONFIG_TARGET} PUBLIC /bigobj)

  # reserve enough stack size
  target_link_options(${OMC_CONFIG_TARGET} PUBLIC "/STACK:8421376")

  # _USE_MATH_DEFINES for PI, NOMINMAX for std::min and std::max
  target_compile_definitions(${OMC_CONFIG_TARGET} PUBLIC _USE_MATH_DEFINES
                                                         NOMINMAX)
endif()
