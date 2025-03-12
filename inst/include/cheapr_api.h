#ifndef CHEAPR_API_H
#define CHEAPR_API_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// -----------------------------------------------------------------------------

namespace cheapr {

static inline R_xlen_t
vec_length(SEXP x) {
  typedef R_xlen_t fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_vec_length");
  return fn(x);
}

static inline SEXP
r_address(SEXP x) {
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_r_address");
  return fn(x);
}

static inline bool
is_compact_seq(SEXP x) {
  typedef bool fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_is_compact_seq");
  return fn(x);
}

static inline SEXP
compact_seq_data(SEXP x) {
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_compact_seq_data");
  return fn(x);
}

static inline SEXP
create_df_row_names(int x) {
  typedef SEXP fn_t(int);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_create_df_row_names");
  return fn(x);
}

static inline SEXP
shallow_copy(SEXP x) {
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_shallow_copy");
  return fn(x);
}

static inline SEXP
set_add_attrs(SEXP x, SEXP attributes, bool add) {
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_set_add_attrs");
  return fn(x, attributes, add);
}

static inline SEXP
set_rm_attrs(SEXP x) {
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_set_rm_attrs");
  return fn(x);
}

static inline SEXP
exclude_locs(SEXP exclude, R_xlen_t xn) {
  typedef SEXP fn_t(SEXP, R_xlen_t);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_exclude_locs");
  return fn(exclude, xn);
}

}

// -----------------------------------------------------------------------------

#endif
