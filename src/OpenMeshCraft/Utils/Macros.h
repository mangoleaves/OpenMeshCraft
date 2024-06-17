#pragma once

#include <cstring>

// NOTE When you add new macros to this file,
//      make sure the new macros are put into right region.

/*********************************************************/
/* Macros that operate on macros *************************/
/*********************************************************/
#pragma region Macros_that_operate_on_macros

#ifndef OMC_CONCATE
	#define _OMC_CONCATE_(A, B) A##B
	#define OMC_CONCATE(A, B) _OMC_CONCATE_(A, B)
#endif

#pragma endregion Macros_that_operate_on_macros

/*********************************************************/
/* Macros that abbreviate others *************************/
/*********************************************************/
#pragma region Macros_that_abbreviate_others

#ifndef OMC_UNUSED
	#define OMC_UNUSED [[maybe_unused]]
#endif

#ifndef OMC_NODISCARD
	#define OMC_NODISCARD [[nodiscard]]
#endif

#ifndef OMC_DEPRECATED
	#define OMC_DEPRECATED [[deprecated]]
#endif

#ifndef OMC_FALLTHROUGH
	#define OMC_FALLTHROUGH [[fallthrough]]
#endif

#pragma endregion Macros_that_abbreviate_others

/*********************************************************/
/* Macros about vectorization    *************************/
/*********************************************************/
#pragma region Macros_about_vectorization

// In MSVC:
//  ----------------------------------------------------------------
// 	  _M_IX86_FP  |  compile options
//       0        |    /arch:IA32
//       1        |    /arch:SSE
//       2        |    /arch:SSE2,/arch:AVX,/arch:AVX2,/arch:AVX512
//  -----------------------------------------------------------------

#if !defined(OMC_SSE2) && defined(OMC_ENABLE_SSE2)
// __SSE2__ is a general and explicit macro to define SSE2.
	#if defined(__SSE2__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2) || \
	  defined(_M_X64)
		#define OMC_SSE2
	#endif
#endif

#if !defined(OMC_AVX) && defined(OMC_ENABLE_AVX)
	#if defined(__AVX__)
		#define OMC_AVX
	#endif
#endif

#if !defined(OMC_AVX2) && defined(OMC_ENABLE_AVX2)
	#if defined(__AVX2__)
		#define OMC_AVX2
	#endif
#endif

#if !defined(OMC_FMA) && defined(OMC_ENABLE_FMA)
	#if defined(_MSC_VER) && defined(__AVX2__)
		#define OMC_FMA
	#elif defined(__FMA__)	// GCC or Clang
		#define OMC_FMA
	#endif
#endif

#pragma endregion Macros_about_vectorization

/*********************************************************/
/* Macros that control Debug and Exception ***************/
/*********************************************************/
#pragma region Macros_that_control_Debug_and_Exception

// Enable exception by default
#if !defined(OMC_ENABLE_EXCEPTION) && !defined(OMC_DISABLE_EXCEPTION)
	#define OMC_ENABLE_EXCEPTION
#endif

// Enable assert by default
#if !defined(OMC_ENABLE_ASSERT) && !defined(OMC_DISABLE_ASSERT)
	#define OMC_ENABLE_ASSERT
#endif

// Disable expensive assert by default
#if !defined(OMC_ENABLE_EXPENSIVE_ASSERT) && \
  !defined(OMC_DISABLE_EXPENSIVE_ASSERT)
	#define OMC_DISABLE_EXPENSIVE_ASSERT
#endif

#pragma endregion Macros_that_control_Debug_and_Exception

/*********************************************************/
/* Macros that convert one to another ********************/
/*********************************************************/

// conver a value to a string
#define OMC_VALUE(x) #x
#define OMC_STRING(x) OMC_VALUE(x)

// convert absolute file path to file name
#if defined(_WIN32)
	#define __FILENAME__ \
		(strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
#elif defined(__unix)
	#define __FILENAME__ \
		(strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#else
	#define __FILENAME__ __FILE__
#endif

/*********************************************************/
/* Macros that control behavior of predicates ************/
/*********************************************************/

// #define OMC_INDIRECT_PRED
#define OMC_OFFSET_PRED

#define OMC_CACHE_SSF