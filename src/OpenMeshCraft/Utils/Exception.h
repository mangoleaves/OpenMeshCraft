#pragma once

#include "Logger.h"

#include <cassert>
#include <iostream>
#include <stdexcept>

namespace OMC {

/*************************************************************************/
/*                     Define Error Types                                */
/*************************************************************************/

#define DEFINE_ERROR_TYPE(err, base)         \
	class err : public base                    \
	{                                          \
	public:                                    \
		explicit err(const std::string &message) \
		  : base(message)                        \
		{                                        \
		}                                        \
	};

// clang-format off
DEFINE_ERROR_TYPE(DomainError      ,std::domain_error)
DEFINE_ERROR_TYPE(InvalidArgument  ,std::invalid_argument)
DEFINE_ERROR_TYPE(LogicError       ,std::logic_error)
DEFINE_ERROR_TYPE(NotImplemented   ,std::runtime_error)
DEFINE_ERROR_TYPE(OutOfRange       ,std::out_of_range)
DEFINE_ERROR_TYPE(OverflowError    ,std::overflow_error)
DEFINE_ERROR_TYPE(RuntimeError     ,std::runtime_error)
DEFINE_ERROR_TYPE(UnderflowError   ,std::underflow_error)
// clang-format on

#undef DEFINE_ERROR_TYPE

/*************************************************************************/
/*                      Throw Error Types                                */
/*************************************************************************/

#if defined(OMC_ENABLE_EXCEPTION)

	#define OMC_THROW(exception, message, ...)                   \
		{                                                          \
			std::string fmt_msg = std::format(message, __VA_ARGS__); \
			OMC::Logger::fatal(fmt_msg);                             \
			throw OMC::exception(fmt_msg);                           \
		}

	#define OMC_THROW_IF(condition, exception, message, ...)     \
		if (condition)                                             \
		{                                                          \
			std::string fmt_msg = std::format(message, __VA_ARGS__); \
			OMC::Logger::fatal(fmt_msg);                             \
			throw OMC::exception(fmt_msg);                           \
		}

#else

	#define OMC_THROW(exception, message, ...)
	#define OMC_THROW_IF(condition, exception, message, ...)

#endif

#if !defined(NOEXCEPTION)

  /* Unconditonal Throw */
	#define OMC_THROW_DOMAIN_ERROR(message, ...) \
		OMC_THROW(DomainError, message, ##__VA_ARGS__)
	#define OMC_THROW_INVALID_ARGUMENT(message, ...) \
		OMC_THROW(InvalidArgument, message, ##__VA_ARGS__)
	#define OMC_THROW_LOGIC_ERROR(message, ...) \
		OMC_THROW(LogicError, message, ##__VA_ARGS__)
	#define OMC_THROW_OUT_OF_RANGE(message, ...) \
		OMC_THROW(OutOfRange, message, ##__VA_ARGS__)
	#define OMC_THROW_OVERFLOW_ERROR(message, ...) \
		OMC_THROW(OverflowError, message, ##__VA_ARGS__)
	#define OMC_THROW_RUNTIME_ERROR(message, ...) \
		OMC_THROW(RuntimeError, message, ##__VA_ARGS__)
	#define OMC_THROW_UNDERFLOW_ERROR(message, ...) \
		OMC_THROW(UnderflowError, message, ##__VA_ARGS__)

  /* Conditonal Throw */
	#define OMC_THROW_DOMAIN_ERROR_IF(cond, message, ...) \
		OMC_THROW_IF(cond, DomainError, message, ##__VA_ARGS__)
	#define OMC_THROW_INVALID_ARGUMENT_IF(cond, message, ...) \
		OMC_THROW_IF(cond, InvalidArgument, message, ##__VA_ARGS__)
	#define OMC_THROW_LOGIC_ERROR_IF(cond, message, ...) \
		OMC_THROW_IF(cond, LogicError, message, ##__VA_ARGS__)
	#define OMC_THROW_OUT_OF_RANGE_IF(cond, message, ...) \
		OMC_THROW_IF(cond, OutOfRange, message, ##__VA_ARGS__)
	#define OMC_THROW_OVERFLOW_ERROR_IF(cond, message, ...) \
		OMC_THROW_IF(cond, OverflowError, message, ##__VA_ARGS__)
	#define OMC_THROW_RUNTIME_ERROR_IF(cond, message, ...) \
		OMC_THROW_IF(cond, RuntimeError, message, ##__VA_ARGS__)
	#define OMC_THROW_UNDERFLOW_ERROR_IF(cond, message, ...) \
		OMC_THROW_IF(cond, UnderflowError, message, ##__VA_ARGS__)

#else

  /* Unconditonal Throw */
	#define OMC_THROW_DOMAIN_ERROR(message, ...)
	#define OMC_THROW_INVALID_ARGUMENT(message, ...)
	#define OMC_THROW_LOGIC_ERROR(message, ...)
	#define OMC_THROW_OUT_OF_RANGE(message, ...)
	#define OMC_THROW_OVERFLOW_ERROR(message, ...)
	#define OMC_THROW_RUNTIME_ERROR(message, ...)
	#define OMC_THROW_UNDERFLOW_ERROR(message, ...)

  /* Conditonal Throw */
	#define OMC_THROW_DOMAIN_ERROR_IF(cond, message, ...)
	#define OMC_THROW_INVALID_ARGUMENT_IF(cond, message, ...)
	#define OMC_THROW_LOGIC_ERROR_IF(cond, message, ...)
	#define OMC_THROW_OUT_OF_RANGE_IF(cond, message, ...)
	#define OMC_THROW_OVERFLOW_ERROR_IF(cond, message, ...)
	#define OMC_THROW_RUNTIME_ERROR_IF(cond, message, ...)
	#define OMC_THROW_UNDERFLOW_ERROR_IF(cond, message, ...)

#endif

#define OMC_THROW_NOT_IMPLEMENTED() \
	throw OMC::NotImplemented(std::string(__func__));

} // namespace OMC

namespace OMC {
/*************************************************************************/
/*                     Define Debug Assertion                            */
/*************************************************************************/

#if defined(OMC_ENABLE_ASSERT)
	#define OMC_ASSERT_AUX_CODE(code) code

  /// Exit program
	#define OMC_EXIT(message, ...)                                         \
		{                                                                    \
			OMC::Logger::fatal(std::format("[{}] [{}] " message, __FILENAME__, \
			                               __LINE__, ##__VA_ARGS__));          \
			std::exit(1);                                                      \
		}

  /// Assert
	#define OMC_ASSERT(condition, message, ...) \
		OMC_THROW_RUNTIME_ERROR_IF(!(condition), message, ##__VA_ARGS__)

#else
	#define OMC_ASSERT_AUX_CODE(code)
	#define OMC_EXIT(message, ...)
	#define OMC_ASSERT(condition, message, ...)
#endif

#ifdef OMC_ENABLE_EXPENSIVE_ASSERT

	#define OMC_EXPENSIVE_ASSERT_AUX_CODE(code) code

	#define OMC_EXPENSIVE_ASSERT(condition, message, ...) \
		OMC_THROW_RUNTIME_ERROR_IF(!(condition), message, ##__VA_ARGS__)

#else
	#define OMC_EXPENSIVE_ASSERT_AUX_CODE(code)
	#define OMC_EXPENSIVE_ASSERT(condition, message, ...)
#endif

} // namespace OMC