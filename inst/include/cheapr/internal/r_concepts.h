#ifndef CHEAPR_R_CONCEPTS_H
#define CHEAPR_R_CONCEPTS_H

#include <cheapr/internal/r_types.h>

namespace cheapr {

template <class... T>
inline constexpr bool always_false = false;

template <class... T>
inline constexpr bool always_true = true;

// Compile-time type check `is<>`
template<typename T, typename U>
    inline constexpr bool is = std::same_as<std::remove_cvref_t<T>, std::remove_cvref_t<U>>;


template<typename T, typename... Args>
inline constexpr bool any = (is<T, Args> || ...);

// Concepts to enable R type templates

template<typename T>
concept RIntegerType = any<T, r_lgl, r_int, r_int64>;

template<typename T>
concept CppIntegerType = std::is_integral_v<std::remove_cvref_t<T>>;

template<typename T>
concept IntegerType = RIntegerType<T> || CppIntegerType<T>;

template<typename T>
concept RMathType = RIntegerType<T> || std::same_as<std::remove_cvref_t<T>, r_dbl>;

template<typename T>
concept CppMathType = std::is_arithmetic_v<std::remove_cvref_t<T>>;

template<typename T>
concept MathType = RMathType<T> || CppMathType<T>;

template<typename T, typename U>
concept AtLeastOneRMathType =
(RMathType<T> || RMathType<U>) && (MathType<T> && MathType<U>);

template<typename T>
concept RType = RMathType<T> || any<T, r_cplx, r_str, r_raw, r_sym, r_sexp>;

template<typename T>
concept CppType = !RType<T>;

template<typename T, typename U>
concept AtLeastOneRType = (RType<T> || RType<U>);

template <typename T>
concept RPtrWritableType = RMathType<T> || any<T, r_cplx, r_raw>;

// template <RType T>
// SEXPTYPE sexp_type(T x){
// static_assert(always_false<T>, "Unsupported type for `sexp_type`");
// }

}

#endif
