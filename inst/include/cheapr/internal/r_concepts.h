#ifndef CHEAPR_R_CONCEPTS_H
#define CHEAPR_R_CONCEPTS_H

#include <cheapr/internal/r_setup.h>

namespace cheapr {

// Forward declare structs to enable defining concepts now
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

// A `SEXP` which we can write data to directly via a pointer
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

template <typename T>
concept RFactor = is<T, r_factors>;

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
concept Scalar = CppScalar<T> || RScalar<T>;

template <typename T>
inline constexpr bool is_sexp = any<T, SEXP, r_sexp>;

namespace internal {

template<typename T>
inline consteval bool can_definitely_be_int(){

    constexpr int max_int = std::numeric_limits<int>::max();

    using xt = std::remove_cvref_t<T>;

    if constexpr (!CppIntegerType<xt>) return false;

    if constexpr (sizeof(xt) <= sizeof(int)) {
        if constexpr (std::is_unsigned_v<xt>) {
             return std::cmp_less_equal(std::numeric_limits<xt>::max(), max_int);
        }
        return true;
    }
    return false;
}

template<typename T>
inline consteval bool can_definitely_be_int64(){

    constexpr int64_t max_int64 = std::numeric_limits<int64_t>::max();

    using xt = std::remove_cvref_t<T>;

    if constexpr (!CppIntegerType<xt>) return false;

    if constexpr (sizeof(xt) <= sizeof(int64_t)) {
        if constexpr (std::is_unsigned_v<xt>) {
             return std::cmp_less_equal(std::numeric_limits<xt>::max(), max_int64);
        }
        return true; 
    }
    return false;
}

// C/C++ -> RVal typenames
// While these aren't the only ways of constructing RVals, they are many-to-one and non-ambiguous

// This is essentially a map of defined non-RVal to RVal conversion operators
// Allowing any of these to be cast to an RVal via static_cast<>
template <typename T>
struct r_val_mapping_impl {};

template <RVal T>
struct r_val_mapping_impl<T> { using type = T; };

template<> struct r_val_mapping_impl<bool>          { using type = r_lgl; };
template<> struct r_val_mapping_impl<int>           { using type = r_int; };
template<> struct r_val_mapping_impl<int64_t>       { using type = r_int64; };
template<> struct r_val_mapping_impl<double>        { using type = r_dbl; };
template<> struct r_val_mapping_impl<const char*>   { using type = r_str; };
template<> struct r_val_mapping_impl<Rcomplex>      { using type = r_cplx; };
template<> struct r_val_mapping_impl<Rbyte>         { using type = r_raw; };
template<> struct r_val_mapping_impl<SEXP>          { using type = r_sexp; };

// R vectors & other containers
template <RVector T>
struct r_val_mapping_impl<T> { using type = r_sexp; };
template <RFactor T>
struct r_val_mapping_impl<T> { using type = r_sexp; };

template<CppMathType T>
struct r_val_mapping_impl<T> {
    using type = std::conditional_t<
    can_definitely_be_int<T>(), 
    r_int,
    std::conditional_t<
        can_definitely_be_int64<T>(), 
        r_int64,
        r_dbl
    >
>;

};

}

template <typename T>
using to_r_val_t = typename internal::r_val_mapping_impl<std::decay_t<T>>::type;

template <typename T>
concept CastableToRVal = requires {
    typename to_r_val_t<T>;
};

// Rules for determining math type promotion in binary operators

namespace internal {

template <RMathType T>
consteval uint8_t r_math_rank() {
    if constexpr (is<T, r_lgl>)   return 0;
    if constexpr (is<T, r_int>)   return 1;
    if constexpr (is<T, r_int64>) return 2;
    if constexpr (is<T, r_dbl>)   return 3;
    return std::numeric_limits<uint8_t>::max();
}

template <RMathType T, RMathType U>
struct common_r_math_impl {
    static constexpr uint8_t rank_t = r_math_rank<T>();
    static constexpr uint8_t rank_u = r_math_rank<U>();
    
    using type = std::conditional_t<(rank_t >= rank_u), T, U>;
};

// If we decide to always upgrade r_lgl to r_int
// template <RMathType T, RMathType U>
// struct common_r_math_impl {
//     static constexpr uint8_t rank_t = r_math_rank<T>();
//     static constexpr uint8_t rank_u = r_math_rank<U>();
    
//     static constexpr bool both_r_lgl = is<T, r_lgl> && is<U, r_lgl>;
//     using type = std::conditional_t<
//         both_r_lgl, r_int, std::conditional_t<
//             (rank_t >= rank_u), T, U
//         >
//     >;
// };

}

template <RMathType T, RMathType U>
using common_r_math_t = typename internal::common_r_math_impl<T, U>::type;

}

#endif
