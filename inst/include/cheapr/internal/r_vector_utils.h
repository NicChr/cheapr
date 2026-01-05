#ifndef CHEAPR_VECTOR_UTILS_H
#define CHEAPR_VECTOR_UTILS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_utf8.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_attrs.h>


namespace cheapr {

namespace internal {

inline r_int_t* integer_ptr(SEXP x){
  return reinterpret_cast<r_int_t*>(INTEGER(x));
}
inline const r_int_t* integer_ptr_ro(SEXP x){
  return reinterpret_cast<const r_int_t*>(INTEGER_RO(x));
}
inline r_bool_t* logical_ptr(SEXP x){
  return reinterpret_cast<r_bool_t*>(integer_ptr(x));
}
inline const r_bool_t* logical_ptr_ro(SEXP x){
  return reinterpret_cast<const r_bool_t*>(integer_ptr_ro(x));
}
inline r_double_t* real_ptr(SEXP x){
  return reinterpret_cast<r_double_t*>(REAL(x));
}
inline const r_double_t* real_ptr_ro(SEXP x){
  return reinterpret_cast<const r_double_t*>(REAL_RO(x));
}
inline r_int64_t* integer64_ptr(SEXP x){
  return reinterpret_cast<r_int64_t*>(real_ptr(x));
}
inline const r_int64_t* integer64_ptr_ro(SEXP x){
  return reinterpret_cast<const r_int64_t*>(real_ptr_ro(x));
}
inline r_complex_t* complex_ptr(SEXP x){
  return reinterpret_cast<r_complex_t*>(COMPLEX(x));
}
inline const r_complex_t* complex_ptr_ro(SEXP x){
  return reinterpret_cast<const r_complex_t*>(COMPLEX_RO(x));
}
inline r_byte_t* raw_ptr(SEXP x){
  return reinterpret_cast<r_byte_t*>(RAW(x));
}
inline const r_byte_t* raw_ptr_ro(SEXP x){
  return reinterpret_cast<const r_byte_t*>(RAW_RO(x));
}
inline const r_string_t* string_ptr_ro(SEXP x){
  return reinterpret_cast<const r_string_t*>(STRING_PTR_RO(x));
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
inline r_bool_t* vector_ptr<r_bool_t>(SEXP x) {
  return reinterpret_cast<r_bool_t*>(LOGICAL(x));
}

template<>
inline const r_bool_t* vector_ptr<const r_bool_t>(SEXP x) {
  return reinterpret_cast<const r_bool_t*>(LOGICAL_RO(x));
}
template<>
inline r_int_t* vector_ptr<r_int_t>(SEXP x) {
  return reinterpret_cast<r_int_t*>(INTEGER(x));
}
template<>
inline const r_int_t* vector_ptr<const r_int_t>(SEXP x) {
  return reinterpret_cast<const r_int_t*>(INTEGER_RO(x));
}

template<>
inline r_double_t* vector_ptr<r_double_t>(SEXP x) {
  return reinterpret_cast<r_double_t*>(REAL(x));
}
template<>
inline const r_double_t* vector_ptr<const r_double_t>(SEXP x) {
  return reinterpret_cast<const r_double_t*>(REAL_RO(x));
}

template<>
inline r_int64_t* vector_ptr<r_int64_t>(SEXP x) {
  return reinterpret_cast<r_int64_t*>(REAL(x));
}
template<>
inline const r_int64_t* vector_ptr<const r_int64_t>(SEXP x) {
  return reinterpret_cast<const r_int64_t*>(REAL_RO(x));
}


template<>
inline r_complex_t* vector_ptr<r_complex_t>(SEXP x) {
  return reinterpret_cast<r_complex_t*>(COMPLEX(x));
}
template<>
inline const r_complex_t* vector_ptr<const r_complex_t>(SEXP x) {
  return reinterpret_cast<const r_complex_t*>(COMPLEX_RO(x));
}

template<>
inline r_byte_t* vector_ptr<r_byte_t>(SEXP x) {
  return reinterpret_cast<r_byte_t*>(RAW(x));
}
template<>
inline const r_byte_t* vector_ptr<const r_byte_t>(SEXP x) {
  return reinterpret_cast<const r_byte_t*>(RAW_RO(x));
}

template<>
inline const r_string_t* vector_ptr<const r_string_t>(SEXP x) {
  return reinterpret_cast<const r_string_t*>(STRING_PTR_RO(x));
}

template<>
inline const sexp_t* vector_ptr<const sexp_t>(SEXP x) {
  return reinterpret_cast<const sexp_t*>(VECTOR_PTR_RO(x));
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
inline void set_value<r_complex_t>(r_complex_t* p_x, const R_xlen_t i, r_complex_t val){
  p_x[i].re() = val.re();
  p_x[i].im() = val.im();
}
template<>
inline void set_value<r_complex_t>(SEXP x, const R_xlen_t i, r_complex_t val){
  auto *p_x = vector_ptr<r_complex_t>(x);
  set_value(p_x, i, val);
}
template<>
inline void set_value<r_string_t>(SEXP x, const R_xlen_t i, r_string_t val){
  SET_STRING_ELT(x, i, static_cast<SEXP>(val));
}
template<>
inline void set_value<const char *>(SEXP x, const R_xlen_t i, const char* val){
  set_value<r_string_t>(x, i, static_cast<r_string_t>(internal::make_utf8_charsxp(val)));
}

// Never use the pointer here to assign
template<>
inline void set_value<sexp_t>(SEXP x, const R_xlen_t i, sexp_t val){
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
inline SEXP new_vector_impl<r_bool_t>(R_xlen_t n) {
  return internal::new_vec(LGLSXP, n);
}
template <>
inline SEXP new_vector_impl<r_int_t>(R_xlen_t n){
  return internal::new_vec(INTSXP, n);
}
template <>
inline SEXP new_vector_impl<r_double_t>(R_xlen_t n){
  return internal::new_vec(REALSXP, n);
}
template <>
inline SEXP new_vector_impl<r_int64_t>(R_xlen_t n){
  SEXP out = Rf_protect(internal::new_vec(REALSXP, n));
  SEXP cls = Rf_protect(internal::make_utf8_strsxp("integer64"));
  Rf_setAttrib(out, R_ClassSymbol, cls);
  Rf_unprotect(2);
  return out;
}
template <>
inline SEXP new_vector_impl<r_string_t>(R_xlen_t n){
  return internal::new_vec(STRSXP, n);
}
template <>
inline SEXP new_vector_impl<r_complex_t>(R_xlen_t n){
  return internal::new_vec(CPLXSXP, n);
}
template <>
inline SEXP new_vector_impl<r_byte_t>(R_xlen_t n){
  return internal::new_vec(RAWSXP, n);
}
template <>
inline SEXP new_vector_impl<sexp_t>(R_xlen_t n){
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
