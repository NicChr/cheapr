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
concept RScalar = RMathType<T> || any<T, r_cplx, r_str, r_raw, r_sym, r_sexp>;

// template<typename T, typename U>
// concept AtLeastOneRScalar = (RScalar<T> || RScalar<U>);

template <typename T>
concept RPtrWritableType = RMathType<T> || any<T, r_cplx, r_raw>;

// Forward declare structs to define concepts now
template<RScalar T>
struct r_vec;

struct r_df;
struct r_factors;
struct r_dates;
struct r_posixcts;

namespace internal {

template<typename T>
struct is_r_vector : std::false_type {};

template<typename T>
struct is_r_vector<r_vec<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_r_vector_v = is_r_vector<std::remove_cvref_t<T>>::value;

}

template<typename T>
concept RVector = internal::is_r_vector_v<T> || is<T, r_dates> || is<T, r_posixcts>;

template <typename T> 
concept RObject = any<T, r_sexp, r_factors, r_df> || RVector<T>;

template <typename T>
concept CppScalar = std::is_scalar_v<T> && !RObject<T> && !RScalar<T>;

template<typename T, typename U>
concept AtLeastOneRScalar = 
(RScalar<T> && RScalar<U>) ||
(RScalar<T> && CppScalar<U>) ||
(CppScalar<T> && RScalar<U>);

}

#endif
