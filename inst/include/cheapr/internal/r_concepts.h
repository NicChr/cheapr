#ifndef CHEAPR_R_CONCEPTS_H
#define CHEAPR_R_CONCEPTS_H

#include <cheapr/internal/r_setup.h>

namespace cheapr {

// Forward declare structs to define concepts now
struct r_lgl;
struct r_int;
struct r_int64;
struct r_dbl;
struct r_str;
struct r_cplx;
struct r_raw;
struct r_sym;
struct r_sexp;

template <class... T>
inline constexpr bool always_false = false;

template <class... T>
inline constexpr bool always_true = true;

// Compile-time type check `is<>`
template <typename T, typename U>
inline constexpr bool is = std::same_as<std::remove_cvref_t<T>, std::remove_cvref_t<U>>;

template <typename T, typename... Args>
inline constexpr bool any = (is<T, Args> || ...);

// Concepts to enable R type templates

template <typename T>
concept RIntegerType = any<T, r_lgl, r_int, r_int64>;

template <typename T>
concept CppIntegerType = std::is_integral_v<std::remove_cvref_t<T>>;

template <typename T>
concept IntegerType = RIntegerType<T> || CppIntegerType<T>;

template <typename T>
concept RMathType = RIntegerType<T> || is<T, r_dbl>;

template<typename T>
concept CppMathType = std::is_arithmetic_v<std::remove_cvref_t<T>>;

template <typename T>
concept MathType = RMathType<T> || CppMathType<T>;

template <typename T, typename U>
concept AtLeastOneRMathType =
(RMathType<T> || RMathType<U>) && (MathType<T> && MathType<U>);

template <typename T>
concept RScalar = RMathType<T> || any<T, r_cplx, r_str, r_raw, r_sym>;

// RVal is anything that can be stored in `r_vec<>`
template <typename T>
concept RVal = RScalar<T> || is<T, r_sexp>;

// template <typename T, typename U>
// concept AtLeastOneRVal = (RVal<T> || RVal<U>);

template <typename T>
concept RPtrWritableType = RMathType<T> || any<T, r_cplx, r_raw>;

// Forward declare structs to define concepts now
template<RVal T>
struct r_vec;

struct r_dates; // Inherits from r_vec<r_int>
struct r_posixcts; // Inherits from r_vec<r_dbl>
struct r_factors;
struct r_df;

namespace internal {

template <typename T>
struct is_r_vector : std::false_type {};

template <typename T>
struct is_r_vector<r_vec<T>> : std::true_type {};

template <typename T>
inline constexpr bool is_r_vector_v = is_r_vector<std::remove_cvref_t<T>>::value;

}

template <typename T>
concept RVector = internal::is_r_vector_v<T> || is<T, r_dates> || is<T, r_posixcts>;

// RObject is any object that can be represented in R - it excludes internal R types like CHARSXP
// Also, these are all implicitly convertible to `SEXP`
template <typename T> 
concept RObject = any<T, r_sexp, r_factors, r_df, r_sym> || RVector<T>;

template <typename T>
concept CppScalar = std::is_scalar_v<T>;

template <typename T, typename U>
concept AtLeastOneRVal = 
(RVal<T> && RVal<U>) ||
(RVal<T> && CppScalar<U>) ||
(CppScalar<T> && RVal<U>);

template <typename T>
concept Scalar = CppScalar<T> || RVal<T>;

template <typename T>
inline constexpr bool is_sexp = any<T, SEXP, r_sexp>;

template <typename From, typename... ToTypes>
concept constructible_to_any = (std::is_constructible_v<ToTypes, From> || ...);

template <typename T>
concept constructible_to_rval = 
    constructible_to_any<T, r_lgl, r_int, r_int64, r_dbl, r_cplx, r_str, r_raw, r_sym, r_sexp>;


// Helper to convert C/C++ typenamesto RVal typenames
// template <typename T>
// struct to_r_val {
//     static_assert(always_false<T>, "No defined C/C++ to RVal typename mapping defined");
// };

// template<>
// struct to_r_val<int> {
//     using type = r_int;
// };

// template<>
// struct to_r_val<double> {
//     using type = r_dbl;
// };

}

#endif
