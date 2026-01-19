#ifndef CHEAPR_R_NAS_H
#define CHEAPR_R_NAS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <limits>

namespace cheapr {
    
// NAs

namespace na {
inline constexpr r_lgl logical = r_na;
inline constexpr r_int integer = r_int{std::numeric_limits<int>::min()};
inline constexpr r_int64 integer64 = r_int64(std::numeric_limits<int64_t>::min());
inline const r_dbl real = r_dbl{NA_REAL};
inline const r_cplx complex = r_cplx{real, real};
inline constexpr r_raw raw = r_raw{0};
inline const r_str string = r_str{NA_STRING};
inline const r_sexp nil = r_null;
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
inline constexpr r_lgl na_value_impl<r_lgl>(){
  return na::logical;
}

template<>
inline constexpr r_int na_value_impl<r_int>(){
  return na::integer;
}

template<>
inline r_dbl na_value_impl<r_dbl>(){
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
inline constexpr r_int64 na_value_impl<r_int64>(){
  return na::integer64;
}

template<>
inline r_cplx na_value_impl<r_cplx>(){
  return na::complex;
}

template<>
inline constexpr r_raw na_value_impl<r_raw>(){
  return r_raw{0};
}

template<>
inline r_str na_value_impl<r_str>(){
  return na::string;
}

template<>
inline SEXP na_value_impl<SEXP>(){
  return r_null;
}

template<>
inline r_sexp na_value_impl<r_sexp>(){
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
  if constexpr (RScalar<T>){
      return x.value == na_value<T>().value;
  } else {
      return false;
  }
}


template<>
inline constexpr bool is_r_na_impl<r_dbl>(r_dbl x){
  return x.value != x.value;
}

template<>
inline constexpr bool is_r_na_impl<r_cplx>(r_cplx x){
  return is_r_na_impl<r_dbl>(x.re()) || is_r_na_impl<r_dbl>(x.im());
}

template<>
inline constexpr bool is_r_na_impl<r_raw>(r_raw x){
  return false;
}
 
// NULL is treated as NA of general R objects
template<>
inline bool is_r_na_impl<SEXP>(SEXP x){
  return x == R_NilValue;
}

}

template<typename T>
inline constexpr bool is_r_na(const T x) {
    return internal::is_r_na_impl<std::remove_cvref_t<T>>(x);
}

}

#endif
