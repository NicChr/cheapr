#ifndef CHEAPR_INTERNAL_HELPERS
#define CHEAPR_INTERNAL_HELPERS

#include <core.h>

namespace cheapr {

// Math helpers

template <typename T>
inline int sign(T x) {
  return (T(0) < x) - (x < T(0));
}

template<typename T>
inline T negate(T x){
  return -x;
}

// trunc but eliminates negative zeroes
template<typename T>
inline double trunc(T x){
  return std::trunc(x) + 0.0;
}

inline SEXP CHEAPR_CORES = R_NilValue;

inline int num_cores(){
  if (is_null(CHEAPR_CORES)){
    CHEAPR_CORES = install_utf8("cheapr.cores");
  }
  int n_cores = Rf_asInteger(Rf_GetOption1(CHEAPR_CORES));
  return n_cores >= 1 ? n_cores : 1;
}

inline R_xlen_t r_length(SEXP x){
  SEXP length_fn = SHIELD(find_pkg_fun("length", "base", false));
  SEXP expr = SHIELD(Rf_lang2(length_fn, x));
  SEXP r_length = SHIELD(Rf_eval(expr, R_BaseEnv));
  R_xlen_t out = TYPEOF(r_length) == INTSXP ? INTEGER(r_length)[0] : REAL(r_length)[0];
  YIELD(3);
  return out;
}

inline bool address_equal(SEXP x, SEXP y){
  return address(x) == address(y);
}

template<typename T>
inline constexpr bool is_r_integerable(T x){
  return between<T>(x, INTEGER_MIN, INTEGER_MAX);
}

inline bool is_int64(SEXP x){
  return Rf_isReal(x) && Rf_inherits(x, "integer64");
}

inline bool is_df(SEXP x){
  return Rf_inherits(x, "data.frame");
}

// Definition of simple vector is one in which
// - It is a vector or it is a list with no class
// - Attributes are data-independent
//
// Care must be taken when combining different simple vectors
// as attributes may only be applicable within a single vector
// e.g. for factors (different levels) and POSIXct (different timezones)
//

inline bool is_bare_atomic(SEXP x){
  return !Rf_isObject(x) && Rf_isVectorAtomic(x);
}

inline bool is_simple_atomic_vec2(SEXP x){
  return is_simple_atomic_vec(x) || is_int64(x);
}

inline bool is_simple_vec2(SEXP x){
  return is_simple_vec(x) || is_int64(x);
}

inline bool is_bare_df(SEXP x){
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  return Rf_length(cls) == 1 &&
    std::strcmp(CHAR(STRING_ELT(cls, 0)), "data.frame") == 0;
}

inline bool is_bare_tbl(SEXP x){
  SEXP xclass = Rf_getAttrib(x, R_ClassSymbol);

  return Rf_length(xclass) == 3 &&
    std::strcmp(CHAR(STRING_ELT(xclass, 0)), "tbl_df") == 0 &&
    std::strcmp(CHAR(STRING_ELT(xclass, 1)), "tbl") == 0 &&
    std::strcmp(CHAR(STRING_ELT(xclass, 2)), "data.frame") == 0;
}

// string paste helper
inline void str_paste(std::string &x, const std::string &sep, const std::string &y){
  x += sep;
  x += y;
}

} // end of cheapr namespace

#endif
