#pragma once

#include <type_traits>

namespace OMC {
/* SFINAE */
/* Substitution failure is not an error. */

// Some useful websites:
// 1. How to enable/disable member functions
// https://www.foonathan.net/2016/12/conditionally-removing-functions/

/// @brief  Check if a class has specified function.
/// @tparam Args Parameter types passed to the function.
/// @note You need to pass function parameters to this template by yourself.
#define CLASS_HAS_FUNC(ClassName, FuncName)                               \
	template <typename... Args>                                             \
	class SFINAE_##ClassName##_has_##FuncName                               \
	{                                                                       \
		template <typename C, typename = decltype(std::declval<C>().FuncName( \
		                        std::declval<Args>()...))>                    \
		static std::true_type test(int);                                      \
		template <typename C>                                                 \
		static std::false_type test(...);                                     \
                                                                          \
	public:                                                                 \
		static constexpr bool value = decltype(test<ClassName>(0))::value;    \
	};

/// @brief  Check if a class has specified const function.
/// @tparam Args Parameter types passed to the function.
/// @note You need to pass function parameters to this template by yourself.
#define CLASS_HAS_CONST_FUNC(ClassName, FuncName)                        \
	template <typename... Args>                                            \
	class SFINAE_##ClassName##_has_##FuncName                              \
	{                                                                      \
		template <typename C,                                                \
		          typename = decltype(std::declval<std::add_const_t<C>>()    \
		                                .FuncName(std::declval<Args>()...))> \
		static std::true_type test(int);                                     \
		template <typename C>                                                \
		static std::false_type test(...);                                    \
                                                                         \
	public:                                                                \
		static constexpr bool value = decltype(test<ClassName>(0))::value;   \
	};

/// @brief  check if a class has specified function.
/// @details
/// A non-static member function receives an instance of the class as first
/// parameter. So we need to explicitly declare the function type of "func".
/// "typename F" is the type of function, "F" is the instance of function.
/// @note You need to pass return type and paramters type to this template by
/// yourself.
#define CLASS_HAS_STATIC_FUNC(ClassName, FuncName)                           \
	template <typename RetT, typename... Args>                                 \
	class SFINAE_##ClassName##_has_##FuncName                                  \
	{                                                                          \
		template <typename F, F>                                                 \
		struct helper                                                            \
		{                                                                        \
		};                                                                       \
		template <typename C>                                                    \
		static constexpr std::true_type                                          \
		test(helper<RetT (*)(Args...), C::FuncName> *);                          \
		template <typename C>                                                    \
		static constexpr std::false_type test(...);                              \
                                                                             \
	public:                                                                    \
		static constexpr bool value = decltype(test<ClassName>(nullptr))::value; \
	};

/// @brief  check if a class has specified type.
#define CLASS_HAS_TYPE(ClassName, TypeName)                                  \
	class SFINAE_##ClassName##_has_##TypeName                                  \
	{                                                                          \
		template <typename C>                                                    \
		static constexpr std::true_type test(typename C::TypeName *);            \
		template <typename C>                                                    \
		static constexpr std::false_type test(...);                              \
                                                                             \
	public:                                                                    \
		static constexpr bool value = decltype(test<ClassName>(nullptr))::value; \
	};

/// @brief  Check if a class has specified type then get the type,
/// otherwise get the default type.
/// @param ClassName The class where to find type.
/// @param TypeName The name of the type.
/// @param DefaultType If doesn't find the type, the return \p DefaultType.
/// @param Result Found result will be stored in \p Result.
#define GET_TYPE_OTHERWISE_DEFAULT(ClassName, TypeName, DefaultType, Result) \
	class SFINAE_get_##TypeName##_from_##ClassName                             \
	{                                                                          \
		template <typename C>                                                    \
		static constexpr auto test(typename C::TypeName *) ->                    \
		  typename C::TypeName;                                                  \
		template <typename C>                                                    \
		static constexpr auto test(...) -> DefaultType;                          \
                                                                             \
	public:                                                                    \
		typedef decltype(test<ClassName>(nullptr)) type;                         \
	};                                                                         \
	using Result = typename SFINAE_get_##TypeName##_from_##ClassName::type;

/// @brief  Check if a class has specified type then get the type,
/// otherwise get the `void` type.
/// @param ClassName The class where to find type.
/// @param TypeName The name of the type.
/// @param Result Found result will be stored in \p Result.
#define GET_TYPE_OTHERWISE_VOID(ClassName, TypeName, Result) \
	GET_TYPE_OTHERWISE_DEFAULT(ClassName, TypeName, void, Result)

/// @brief  check if a class has specified value.
#define CLASS_HAS_VALUE(ClassName, ValueName, Result)                        \
	class SFINAE_##ClassName##_has_##ValueName                                 \
	{                                                                          \
		template <typename C>                                                    \
		static constexpr std::true_type test(typename decltype(C::ValueName) *); \
		template <typename C>                                                    \
		static constexpr std::false_type test(...);                              \
                                                                             \
	public:                                                                    \
		static constexpr bool value = decltype(test<ClassName>(nullptr))::value; \
	};                                                                         \
	static constexpr bool Result = SFINAE_##ClassName##_has_##ValueName::value;

#define GET_VALUE_OTHERWISE_DEFAULT(ClassName, ValueType, ValueName, \
                                    DefaultValue, Result)            \
	class SFINAE_get_##ValueName##_from_##ClassName                    \
	{                                                                  \
		template <typename C>                                            \
		static constexpr ValueType test(decltype(C::ValueName) *)        \
		{                                                                \
			return C::ValueName;                                           \
		}                                                                \
		template <typename C>                                            \
		static constexpr ValueType test(...)                             \
		{                                                                \
			return DefaultValue;                                           \
		}                                                                \
                                                                     \
	public:                                                            \
		static constexpr ValueType value = test<ClassName>(nullptr);     \
	};                                                                 \
	static constexpr ValueType Result =                                \
	  SFINAE_get_##ValueName##_from_##ClassName::value;

} // namespace OMC