#ifndef CHEAPR_VECTOR_UTILS_H
#define CHEAPR_VECTOR_UTILS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>


namespace cheapr {

namespace internal {

inline r_sexp new_vec(SEXPTYPE type, r_size_t n){
  return r_sexp(cpp11::safe[Rf_allocVector](type, n));
}

inline r_sexp coerce_vec(SEXP x, SEXPTYPE type){
  return r_sexp(cpp11::safe[Rf_coerceVector](x, type));
}


template<RType T>
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

template<>
inline const r_str* vector_ptr<const r_str>(SEXP x) {
  return reinterpret_cast<const r_str*>(STRING_PTR_RO(x));
}

template<>
inline const r_sexp* vector_ptr<const r_sexp>(SEXP x) {
  return reinterpret_cast<const r_sexp*>(VECTOR_PTR_RO(x));
}

// R vector getters + setters

template<RType T>
inline T get_value(SEXP x, const r_size_t i){
  return vector_ptr<T>(x)[i];
}

template<RType T>
inline T get_value(T *p_x, const r_size_t i){
  return p_x[i];
}

template<RType T>
inline const T get_value(const T *p_x, const r_size_t i){
  return p_x[i];
}


template<typename T>
requires (RType<T> || is<T, const char *>)
inline void set_value(T *p_x, const r_size_t i, T val){
  p_x[i] = val;
}

template<typename T>
requires (RType<T> || is<T, const char *>)
inline void set_value(SEXP x, const r_size_t i, T val){
  vector_ptr<T>(x)[i] = val;
}
template<>
inline void set_value<r_cplx>(r_cplx* p_x, const r_size_t i, r_cplx val){
  p_x[i].re() = val.re();
  p_x[i].im() = val.im();
}
template<>
inline void set_value<r_cplx>(SEXP x, const r_size_t i, r_cplx val){
  auto *p_x = vector_ptr<r_cplx>(x);
  set_value(p_x, i, val);
}
template<>
inline void set_value<r_str>(SEXP x, const r_size_t i, r_str val){
  SET_STRING_ELT(x, i, static_cast<SEXP>(val));
}
template<>
inline void set_value<const char *>(SEXP x, const r_size_t i, const char* val){
  set_value<r_str>(x, i, r_str(val));
}

// Never use the pointer here to assign
template<>
inline void set_value<r_sexp>(SEXP x, const r_size_t i, r_sexp val){
  SET_VECTOR_ELT(x, i, val);
}

// One-parameter template version
template <RType T>
inline r_sexp new_vector_impl(r_size_t n) {
  static_assert(
    always_false<T>,
    "Unimplemented `new_vector_impl` specialisation"
  );
  return r_null;
}

template <>
inline r_sexp new_vector_impl<r_lgl>(r_size_t n) {
  return internal::new_vec(LGLSXP, n);
}
template <>
inline r_sexp new_vector_impl<r_int>(r_size_t n){
  return internal::new_vec(INTSXP, n);
}
template <>
inline r_sexp new_vector_impl<r_dbl>(r_size_t n){
  return internal::new_vec(REALSXP, n);
}
template <>
inline r_sexp new_vector_impl<r_int64>(r_size_t n){
  r_sexp out = r_sexp(internal::new_vec(REALSXP, n));
  r_sexp cls = r_sexp(Rf_ScalarString(r_str("integer64")));
  Rf_setAttrib(out, R_ClassSymbol, cls);
  return out;
}
template <>
inline r_sexp new_vector_impl<r_str>(r_size_t n){
  return internal::new_vec(STRSXP, n);
}
template <>
inline r_sexp new_vector_impl<r_cplx>(r_size_t n){
  return internal::new_vec(CPLXSXP, n);
}
template <>
inline r_sexp new_vector_impl<r_raw>(r_size_t n){
  return internal::new_vec(RAWSXP, n);
}
template <>
inline r_sexp new_vector_impl<r_sexp>(r_size_t n){
  return internal::new_vec(VECSXP, n);
}

template <RType T>
inline r_sexp new_scalar_vector(T default_value) {
    r_sexp out = new_vector_impl<T>(1);
    set_value<T>(out, 0, default_value);
    return out;
}

}

}

#endif
