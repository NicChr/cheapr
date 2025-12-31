#ifndef CHEAPR_R_NAS
#define CHEAPR_R_NAS

#include <limits>
#include <r_types.h>

namespace cheapr {
    
// NAs

namespace na {
inline constexpr r_bool_t logical = r_na;
inline constexpr r_int_t integer = r_int_t{std::numeric_limits<int>::min()};
inline constexpr r_int64_t integer64 = r_int64_t(std::numeric_limits<int64_t>::min());
inline const r_double_t real = r_double_t{NA_REAL};
inline const r_complex_t complex = r_complex_t{real, real};
inline constexpr r_byte_t raw = r_byte_t{0};
inline const r_string_t string = r_string_t{NA_STRING};
inline const sexp_t nil = r_null;
}

namespace internal {

template<typename T>
inline constexpr T na_value_impl() {
  static_assert(
    always_false<T>,
    "Unimplemented `na_value` specialisation"
  );
  return T{};
}

template<>
inline constexpr r_bool_t na_value_impl<r_bool_t>(){
  return na::logical;
}

template<>
inline constexpr r_int_t na_value_impl<r_int_t>(){
  return na::integer;
}

template<>
inline r_double_t na_value_impl<r_double_t>(){
  return na::real;
}
template<>
inline constexpr int na_value_impl<int>(){
  return na::integer.value;
}

template<>
inline double na_value_impl<double>(){
  return na::real.value;
}

template<>
inline constexpr r_int64_t na_value_impl<r_int64_t>(){
  return na::integer64;
}

template<>
inline r_complex_t na_value_impl<r_complex_t>(){
  return na::complex;
}

template<>
inline constexpr r_byte_t na_value_impl<r_byte_t>(){
  return r_byte_t{0};
}

template<>
inline r_string_t na_value_impl<r_string_t>(){
  return na::string;
}

template<>
inline SEXP na_value_impl<SEXP>(){
  return r_null;
}

template<>
inline sexp_t na_value_impl<sexp_t>(){
  return r_null;
}

}

template<typename T>
inline constexpr auto na_value() {
  return internal::na_value_impl<std::remove_cvref_t<T>>();
}

namespace internal {
template<typename T>
inline constexpr bool is_r_na_impl(T x) {
    if constexpr (RType<T>){
        return x == na_value<T>();
    } else {
        return false;
    }
}


template<>
inline constexpr bool is_r_na_impl<r_double_t>(r_double_t x){
  return x.value != x.value;
}
template<>
inline constexpr bool is_r_na_impl<double>(double x){
  return x != x;
}

template<>
inline constexpr bool is_r_na_impl<r_complex_t>(r_complex_t x){
  return is_r_na<r_double_t>(x.re()) || is_r_na<r_double_t>(x.im());
}

template<>
inline constexpr bool is_r_na_impl<r_byte_t>(r_byte_t x){
  return false;
}

// NULL is treated as NA of general R objects
template<>
inline bool is_r_na_impl<SEXP>(SEXP x){
  return is_null(x);
}

}

template<typename T>
inline constexpr bool is_r_na(const T x) {
    return is_r_na_impl<std::remove_cvref_t<T>>(x);
}

}

#endif
