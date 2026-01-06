#ifndef CHEAPR_VECTOR_UTILS_H
#define CHEAPR_VECTOR_UTILS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_utf8.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>


namespace cheapr {

namespace internal {

inline r_int* integer_ptr(SEXP x){
  return reinterpret_cast<r_int*>(INTEGER(x));
}
inline const r_int* integer_ptr_ro(SEXP x){
  return reinterpret_cast<const r_int*>(INTEGER_RO(x));
}
inline r_lgl* logical_ptr(SEXP x){
  return reinterpret_cast<r_lgl*>(integer_ptr(x));
}
inline const r_lgl* logical_ptr_ro(SEXP x){
  return reinterpret_cast<const r_lgl*>(integer_ptr_ro(x));
}
inline r_dbl* real_ptr(SEXP x){
  return reinterpret_cast<r_dbl*>(REAL(x));
}
inline const r_dbl* real_ptr_ro(SEXP x){
  return reinterpret_cast<const r_dbl*>(REAL_RO(x));
}
inline r_int64* integer64_ptr(SEXP x){
  return reinterpret_cast<r_int64*>(real_ptr(x));
}
inline const r_int64* integer64_ptr_ro(SEXP x){
  return reinterpret_cast<const r_int64*>(real_ptr_ro(x));
}
inline r_cplx* complex_ptr(SEXP x){
  return reinterpret_cast<r_cplx*>(COMPLEX(x));
}
inline const r_cplx* complex_ptr_ro(SEXP x){
  return reinterpret_cast<const r_cplx*>(COMPLEX_RO(x));
}
inline r_raw* raw_ptr(SEXP x){
  return reinterpret_cast<r_raw*>(RAW(x));
}
inline const r_raw* raw_ptr_ro(SEXP x){
  return reinterpret_cast<const r_raw*>(RAW_RO(x));
}
inline const r_str* string_ptr_ro(SEXP x){
  return reinterpret_cast<const r_str*>(STRING_PTR_RO(x));
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
inline T get_value(SEXP x, const R_xlen_t i){
  return vector_ptr<T>(x)[i];
}

template<RType T>
inline T get_value(T *p_x, const R_xlen_t i){
  return p_x[i];
}

template<RType T>
inline const T get_value(const T *p_x, const R_xlen_t i){
  return p_x[i];
}


template<typename T>
requires (RType<T> || is<T, const char *>)
inline void set_value(T *p_x, const R_xlen_t i, T val){
  p_x[i] = val;
}

template<typename T>
requires (RType<T> || is<T, const char *>)
inline void set_value(SEXP x, const R_xlen_t i, T val){
  vector_ptr<T>(x)[i] = val;
}
template<>
inline void set_value<r_cplx>(r_cplx* p_x, const R_xlen_t i, r_cplx val){
  p_x[i].re() = val.re();
  p_x[i].im() = val.im();
}
template<>
inline void set_value<r_cplx>(SEXP x, const R_xlen_t i, r_cplx val){
  auto *p_x = vector_ptr<r_cplx>(x);
  set_value(p_x, i, val);
}
template<>
inline void set_value<r_str>(SEXP x, const R_xlen_t i, r_str val){
  SET_STRING_ELT(x, i, static_cast<SEXP>(val));
}
template<>
inline void set_value<const char *>(SEXP x, const R_xlen_t i, const char* val){
  set_value<r_str>(x, i, static_cast<r_str>(internal::make_utf8_charsxp(val)));
}

// Never use the pointer here to assign
template<>
inline void set_value<r_sexp>(SEXP x, const R_xlen_t i, r_sexp val){
  SET_VECTOR_ELT(x, i, val);
}

// One-parameter template version
template <RType T>
inline SEXP new_vector_impl(R_xlen_t n) {
  static_assert(
    always_false<T>,
    "Unimplemented `new_vector_impl` specialisation"
  );
  return r_null;
}

template <>
inline SEXP new_vector_impl<r_lgl>(R_xlen_t n) {
  return internal::new_vec(LGLSXP, n);
}
template <>
inline SEXP new_vector_impl<r_int>(R_xlen_t n){
  return internal::new_vec(INTSXP, n);
}
template <>
inline SEXP new_vector_impl<r_dbl>(R_xlen_t n){
  return internal::new_vec(REALSXP, n);
}
template <>
inline SEXP new_vector_impl<r_int64>(R_xlen_t n){
  SEXP out = Rf_protect(internal::new_vec(REALSXP, n));
  SEXP cls = Rf_protect(internal::make_utf8_strsxp("integer64"));
  Rf_setAttrib(out, R_ClassSymbol, cls);
  Rf_unprotect(2);
  return out;
}
template <>
inline SEXP new_vector_impl<r_str>(R_xlen_t n){
  return internal::new_vec(STRSXP, n);
}
template <>
inline SEXP new_vector_impl<r_cplx>(R_xlen_t n){
  return internal::new_vec(CPLXSXP, n);
}
template <>
inline SEXP new_vector_impl<r_raw>(R_xlen_t n){
  return internal::new_vec(RAWSXP, n);
}
template <>
inline SEXP new_vector_impl<r_sexp>(R_xlen_t n){
  return internal::new_vec(VECSXP, n);
}

template <RType T>
inline SEXP new_scalar_vector(T default_value) {
    SEXP out = SHIELD(new_vector_impl<T>(1));
    set_value<T>(out, 0, default_value);
    YIELD(1);
    return out;
}

}

}

#endif
