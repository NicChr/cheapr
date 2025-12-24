#ifndef CHEAPR_INTERNAL_HELPERS_H
#define CHEAPR_INTERNAL_HELPERS_H

#include <cheapr/internal/c_core.h>

namespace cheapr {

namespace internal {

template<typename... Args>
inline SEXP eval_pkg_fun(const char* fn, const char* pkg, SEXP envir, Args... args){
  SEXP r_fn = SHIELD(fn::find_pkg_fun(fn, pkg, true));
  SEXP out = SHIELD(fn::eval_fn(r_fn, envir, args...));
  YIELD(2);
  return out;
}

inline R_xlen_t r_length(SEXP x){
  SEXP r_length = SHIELD(eval_pkg_fun("length", "base", env::base_env, x));
  R_xlen_t out = TYPEOF(r_length) == INTSXP ? internal::integer_ptr(r_length)[0] : internal::real_ptr(r_length)[0];
  YIELD(1);
  return out;
}

inline bool address_equal(SEXP x, SEXP y){
  return address(x) == address(y);
}

inline bool is_int64(SEXP x){
  return Rf_isReal(x) && inherits1(x, "integer64");
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
  return !vec::is_object(x) && vec::is_atomic(x);
}

inline bool cheapr_is_simple_atomic_vec(SEXP x){
  return (
      vec::is_atomic(x) && (
          !vec::is_object(x) || (
              inherits1(x, "Date") || inherits1(x, "factor") ||
              inherits1(x, "POSIXct")
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
  SEXP cls = cheapr::attr::get_old_class(x);
  return Rf_length(cls) == 1 &&
    std::strcmp(CHAR(STRING_ELT(cls, 0)), "data.frame") == 0;
}

inline bool is_bare_tbl(SEXP x){
  SEXP xclass = cheapr::attr::get_old_class(x);

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
