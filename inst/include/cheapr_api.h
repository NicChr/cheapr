#ifndef CHEAPR_API_H
#define CHEAPR_API_H

#include <stdbool.h>
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

static inline R_xlen_t
unnested_length(SEXP x){
  typedef R_xlen_t fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_unnested_length");
  return fn(x);
}

static inline SEXP
drop_null(SEXP l, bool always_shallow_copy) {
  typedef SEXP fn_t(SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_drop_null");
  return fn(l, always_shallow_copy);
}

static inline SEXP
list_as_df(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_list_as_df");
  return fn(x);
}

static inline SEXP
lengths(SEXP x, bool names){
  typedef SEXP fn_t(SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_lengths");
  return fn(x, names);
}

static inline SEXP
df_slice(SEXP x, SEXP indices, bool check){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_slice");
  return fn(x, indices, check);
}

static inline SEXP
df_select(SEXP x, SEXP locs){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_select");
  return fn(x, locs);
}

static inline SEXP
sset(SEXP x, SEXP indices){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_sset");
  return fn(x, indices);
}

static inline SEXP
val_find(SEXP x, SEXP value, bool invert){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_val_find");
  return fn(x, value, invert);
}

static inline SEXP
sequence(SEXP size, SEXP from, SEXP by){
  typedef SEXP fn_t(SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_sequence");
  return fn(size, from, by);
}

static inline SEXP
seq_len(R_xlen_t n){
  typedef SEXP fn_t(R_xlen_t);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_seq_len");
  return fn(n);
}

static inline bool
is_simple_atomic_vec(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_is_simple_atomic_vec");
  return fn(x);
}

static inline SEXP
rep_len(SEXP x, int length){
  typedef SEXP fn_t(SEXP, int);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_rep_len");
  return fn(x, length);
}

static inline SEXP
recycle(SEXP x, SEXP length){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_recycle");
  return fn(x, length);
}

static inline SEXP
c(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_c");
  return fn(x);
}

static inline SEXP
name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep){
  typedef SEXP fn_t(SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_name_repair");
  return fn(names, dup_sep, empty_sep);
}

static inline SEXP
unique(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_unique");
  return fn(x);
}

static inline SEXP
setdiff(SEXP x, SEXP y){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_setdiff");
  return fn(x, y);
}

static inline SEXP
get_ptype(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_get_ptype");
  return fn(x);
}

}

// -----------------------------------------------------------------------------

#endif
