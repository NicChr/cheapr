#ifndef CHEAPR_VECTOR_UTILS_H
#define CHEAPR_VECTOR_UTILS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_types.h>


namespace cheapr {

namespace internal {

template <typename T>
using vec_ptr_t = std::conditional_t<std::is_const_v<T>, const unwrapped_t<T>*, unwrapped_t<T>*>;

inline r_sexp new_vec(SEXPTYPE type, r_size_t n){
  return r_sexp(cpp11::safe[Rf_allocVector](type, n));
}
// Not used, to be removed later
inline r_sexp coerce_vec(SEXP x, SEXPTYPE type){
  return r_sexp(cpp11::safe[Rf_coerceVector](x, type));
}

template<RPtrWritableType T>
inline vec_ptr_t<T> vector_ptr(SEXP x) {
  static_assert(
    always_false<T>,
    "Unsupported type for vector_ptr"
  );
  return nullptr;
}

template<>
inline int* vector_ptr<r_lgl>(SEXP x) {
  return LOGICAL(x);
}

template<>
inline const int* vector_ptr<const r_lgl>(SEXP x) {
  return LOGICAL_RO(x);
}
template<>
inline int* vector_ptr<r_int>(SEXP x) {
  return INTEGER(x);
}
template<>
inline const int* vector_ptr<const r_int>(SEXP x) {
  return INTEGER_RO(x);
}

template<>
inline double* vector_ptr<r_dbl>(SEXP x) {
  return REAL(x);
}
template<>
inline const double* vector_ptr<const r_dbl>(SEXP x) {
  return REAL_RO(x);
}

template<>
inline int64_t* vector_ptr<r_int64>(SEXP x) {
  return reinterpret_cast<int64_t*>(REAL(x));
}
template<>
inline const int64_t* vector_ptr<const r_int64>(SEXP x) {
  return reinterpret_cast<const int64_t*>(REAL_RO(x));
}


template<>
inline Rcomplex* vector_ptr<r_cplx>(SEXP x) {
  return COMPLEX(x);
}
template<>
inline const Rcomplex* vector_ptr<const r_cplx>(SEXP x) {
  return COMPLEX_RO(x);
}

template<>
inline Rbyte* vector_ptr<r_raw>(SEXP x) {
  return RAW(x);
}
template<>
inline const Rbyte* vector_ptr<const r_raw>(SEXP x) {
  return RAW_RO(x);
}

// Internal vec constructor
template <RVal T>
inline r_sexp new_vec_impl(r_size_t n) {
  static_assert(
    always_false<T>,
    "Unimplemented `new_vec_impl` specialisation"
  );
  return r_null;
}

template <>
inline r_sexp new_vec_impl<r_lgl>(r_size_t n) {
  return internal::new_vec(LGLSXP, n);
}
template <>
inline r_sexp new_vec_impl<r_int>(r_size_t n){
  return internal::new_vec(INTSXP, n);
}
template <>
inline r_sexp new_vec_impl<r_dbl>(r_size_t n){
  return internal::new_vec(REALSXP, n);
}
template <>
inline r_sexp new_vec_impl<r_int64>(r_size_t n){
  r_sexp out = r_sexp(internal::new_vec(REALSXP, n));
  r_sexp cls = r_sexp(Rf_ScalarString(r_str("integer64")));
  Rf_setAttrib(out, R_ClassSymbol, cls);
  return out;
}
template <>
inline r_sexp new_vec_impl<r_str>(r_size_t n){
  return internal::new_vec(STRSXP, n);
}
template <>
inline r_sexp new_vec_impl<r_cplx>(r_size_t n){
  return internal::new_vec(CPLXSXP, n);
}
template <>
inline r_sexp new_vec_impl<r_raw>(r_size_t n){
  return internal::new_vec(RAWSXP, n);
}
template <>
inline r_sexp new_vec_impl<r_sexp>(r_size_t n){
  return internal::new_vec(VECSXP, n);
}


template <RVal T>
inline r_sexp new_scalar_vec(T default_value) {
  static_assert(
    always_false<T>,
    "Unimplemented `new_scalar_vec` specialisation"
  );
  return r_null;
}

template <>
inline r_sexp new_scalar_vec<r_lgl>(r_lgl default_value) {
  return r_sexp(Rf_ScalarLogical(unwrap(default_value)));
}
template <>
inline r_sexp new_scalar_vec<r_int>(r_int default_value){
  return r_sexp(Rf_ScalarInteger(unwrap(default_value)));
}
template <>
inline r_sexp new_scalar_vec<r_dbl>(r_dbl default_value){
  return r_sexp(Rf_ScalarReal(unwrap(default_value)));
}
template <>
inline r_sexp new_scalar_vec<r_int64>(r_int64 default_value){
  r_sexp out = new_vec_impl<r_int64>(1);
  REAL(out.value)[0] = default_value;
  return out;
}
template <>
inline r_sexp new_scalar_vec<r_str>(r_str default_value){
  return r_sexp(Rf_ScalarString(unwrap(default_value)));
}
template <>
inline r_sexp new_scalar_vec<r_cplx>(r_cplx default_value){
  return r_sexp(Rf_ScalarComplex(unwrap(default_value)));
}
template <>
inline r_sexp new_scalar_vec<r_raw>(r_raw default_value){
  return r_sexp(Rf_ScalarRaw(unwrap(default_value)));
}
template <>
inline r_sexp new_scalar_vec<r_sexp>(r_sexp default_value){
  r_sexp out = new_vec_impl<r_sexp>(1);
  SET_VECTOR_ELT(out.value, 0, unwrap(default_value));
  return out;
}

}

}

#endif
