#ifndef CHEAPR_R_CONCEPTS_H
#define CHEAPR_R_CONCEPTS_H

#include <cheapr/internal/r_types.h>

namespace cheapr {
// Compile-time type check `is<>`
template<typename T, typename U>
    inline constexpr bool is = std::same_as<std::remove_cvref_t<T>, std::remove_cvref_t<U>>;

// Concepts to enable R type templates

template<typename T>
concept RIntegerType = std::same_as<std::remove_cvref_t<T>, r_bool_t> ||
std::same_as<std::remove_cvref_t<T>, r_int_t> ||
std::same_as<std::remove_cvref_t<T>, r_int64_t>;

template<typename T>
concept CppIntegerType = std::is_integral_v<std::remove_cvref_t<T>>;

template<typename T>
concept IntegerType = RIntegerType<T> || CppIntegerType<T>;

template<typename T>
concept RMathType = RIntegerType<T> || std::same_as<std::remove_cvref_t<T>, r_double_t>;

template<typename T>
concept CppMathType = std::is_arithmetic_v<std::remove_cvref_t<T>>;

template<typename T>
concept MathType = RMathType<T> || CppMathType<T>;

template<typename T, typename U>
concept AtLeastOneRMathType =
(RMathType<T> || RMathType<U>) && (MathType<T> && MathType<U>);

template<typename T>
concept RType = RMathType<T> ||
std::same_as<std::remove_cvref_t<T>, r_complex_t> ||
std::same_as<std::remove_cvref_t<T>, r_string_t> ||
std::same_as<std::remove_cvref_t<T>, r_byte_t> ||
std::same_as<std::remove_cvref_t<T>, r_symbol_t> ||
std::same_as<std::remove_cvref_t<T>, sexp_t>;

template<typename T>
concept CppType = !RType<T>;

template<typename T, typename U>
concept AtLeastOneRType = (RType<T> || RType<U>);

template <typename T>
concept RPtrWritableType = RMathType<T> ||
is<T, r_complex_t> ||
is<T, r_byte_t>;

template <class... T>
inline constexpr bool always_false = false;

template <class... T>
inline constexpr bool always_true = true;

}

#endif
