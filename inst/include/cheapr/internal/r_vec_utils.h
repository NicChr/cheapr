#ifndef CHEAPR_VECTOR_UTILS_H
#define CHEAPR_VECTOR_UTILS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_types.h>


namespace cheapr {

namespace internal {

inline r_sexp new_vec(SEXPTYPE type, r_size_t n){
  return r_sexp(cpp11::safe[Rf_allocVector](type, n));
}
// Not used, to be removed later
inline r_sexp coerce_vec(SEXP x, SEXPTYPE type){
  return r_sexp(cpp11::safe[Rf_coerceVector](x, type));
}

template<RPtrWritableType T>
inline T* vector_ptr(SEXP x) {
  static_assert(
    always_false<T>,
    "Unsupported type for vector_ptr"
  );
  return nullptr;
}

template<>
inline r_lgl* vector_ptr<r_lgl>(SEXP x) {
  return reinterpret_cast<r_lgl*>(LOGICAL(x));
}

template<>
inline const r_lgl* vector_ptr<const r_lgl>(SEXP x) {
  return reinterpret_cast<const r_lgl*>(LOGICAL_RO(x));
}
template<>
inline r_int* vector_ptr<r_int>(SEXP x) {
  return reinterpret_cast<r_int*>(INTEGER(x));
}
template<>
inline const r_int* vector_ptr<const r_int>(SEXP x) {
  return reinterpret_cast<const r_int*>(INTEGER_RO(x));
}

template<>
inline r_dbl* vector_ptr<r_dbl>(SEXP x) {
  return reinterpret_cast<r_dbl*>(REAL(x));
}
template<>
inline const r_dbl* vector_ptr<const r_dbl>(SEXP x) {
  return reinterpret_cast<const r_dbl*>(REAL_RO(x));
}

template<>
inline r_int64* vector_ptr<r_int64>(SEXP x) {
  return reinterpret_cast<r_int64*>(REAL(x));
}
template<>
inline const r_int64* vector_ptr<const r_int64>(SEXP x) {
  return reinterpret_cast<const r_int64*>(REAL_RO(x));
}


template<>
inline r_cplx* vector_ptr<r_cplx>(SEXP x) {
  return reinterpret_cast<r_cplx*>(COMPLEX(x));
}
template<>
inline const r_cplx* vector_ptr<const r_cplx>(SEXP x) {
  return reinterpret_cast<const r_cplx*>(COMPLEX_RO(x));
}

template<>
inline r_raw* vector_ptr<r_raw>(SEXP x) {
  return reinterpret_cast<r_raw*>(RAW(x));
}
template<>
inline const r_raw* vector_ptr<const r_raw>(SEXP x) {
  return reinterpret_cast<const r_raw*>(RAW_RO(x));
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
  return r_sexp(Rf_ScalarLogical(default_value.value));
}
template <>
inline r_sexp new_scalar_vec<r_int>(r_int default_value){
  return r_sexp(Rf_ScalarInteger(default_value.value));
}
template <>
inline r_sexp new_scalar_vec<r_dbl>(r_dbl default_value){
  return r_sexp(Rf_ScalarReal(default_value.value));
}
template <>
inline r_sexp new_scalar_vec<r_int64>(r_int64 default_value){
  r_sexp out = new_vec_impl<r_int64>(1);
  REAL(out.value)[0] = default_value;
  return out;
}
template <>
inline r_sexp new_scalar_vec<r_str>(r_str default_value){
  return r_sexp(Rf_ScalarString(default_value.value));
}
template <>
inline r_sexp new_scalar_vec<r_cplx>(r_cplx default_value){
  return r_sexp(Rf_ScalarComplex(default_value.value));
}
template <>
inline r_sexp new_scalar_vec<r_raw>(r_raw default_value){
  return r_sexp(Rf_ScalarRaw(default_value.value));
}
template <>
inline r_sexp new_scalar_vec<r_sexp>(r_sexp default_value){
  r_sexp out = new_vec_impl<r_sexp>(1);
  SET_VECTOR_ELT(out.value, 0, default_value);
  return out;
}

}

}

#endif
