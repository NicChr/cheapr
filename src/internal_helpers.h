#ifndef CHEAPR_INTERNAL_HELPERS_H
#define CHEAPR_INTERNAL_HELPERS_H

#include <core.h>

namespace cheapr {

namespace internal {

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

inline SEXP CHEAPR_CORES = r_null;

inline int num_cores(){
  if (vec::is_null(CHEAPR_CORES)){
    CHEAPR_CORES = install_utf8("cheapr.cores");
  }
  int n_cores = Rf_asInteger(Rf_GetOption1(CHEAPR_CORES));
  return n_cores >= 1 ? n_cores : 1;
}

template<typename... Args>
inline SEXP eval_pkg_fun(const char* fn, const char* pkg, SEXP envir, Args... args){
  SEXP r_fn = SHIELD(find_pkg_fun(fn, pkg, true));
  SEXP out = SHIELD(eval_fun(r_fn, envir, args...));
  YIELD(2);
  return out;
}

inline R_xlen_t r_length(SEXP x){
  SEXP r_length = SHIELD(eval_pkg_fun("length", "base", base_env, x));
  R_xlen_t out = TYPEOF(r_length) == INTSXP ? INTEGER(r_length)[0] : REAL(r_length)[0];
  YIELD(1);
  return out;
}

inline bool address_equal(SEXP x, SEXP y){
  return vec::address(x) == vec::address(y);
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
  return !vec::is_object(x) && Rf_isVectorAtomic(x);
}

inline bool cheapr_is_simple_atomic_vec(SEXP x){
  return (
      Rf_isVectorAtomic(x) && (
          !vec::is_object(x) || (
              Rf_inherits(x, "Date") || Rf_inherits(x, "factor") ||
              Rf_inherits(x, "POSIXct")
          )
      )
  );
}

inline bool is_bare_list(SEXP x){
  return (!vec::is_object(x) && TYPEOF(x) == VECSXP);
}

inline bool cheapr_is_simple_vec(SEXP x){
  return (cheapr_is_simple_atomic_vec(x) || is_bare_list(x));
}

inline bool cheapr_is_simple_atomic_vec2(SEXP x){
  return cheapr_is_simple_atomic_vec(x) || is_int64(x);
}

inline bool cheapr_is_simple_vec2(SEXP x){
  return cheapr_is_simple_vec(x) || is_int64(x);
}

inline bool is_bare_df(SEXP x){
  SEXP cls = cheapr::vec::get_attrib(x, R_ClassSymbol);
  return Rf_length(cls) == 1 &&
    std::strcmp(CHAR(STRING_ELT(cls, 0)), "data.frame") == 0;
}

inline bool is_bare_tbl(SEXP x){
  SEXP xclass = cheapr::vec::get_attrib(x, R_ClassSymbol);

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

} // end of internal namespace

} // end of cheapr namespace

#endif
