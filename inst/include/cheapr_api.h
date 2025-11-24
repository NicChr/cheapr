#ifndef CHEAPR_API_H
#define CHEAPR_API_H

// Core constants and small functions
#include <core.h>
// R tools for exporting C fns
#include <R_ext/Rdynload.h>

// -----------------------------------------------------------------------------

namespace cheapr {

inline bool
is_compact_seq(SEXP x) {
  typedef bool fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_is_compact_seq");
  return fn(x);
}

inline SEXP
set_add_attrs(SEXP x, SEXP attributes, bool add) {
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_set_add_attrs");
  return fn(x, attributes, add);
}

inline SEXP
set_rm_attrs(SEXP x) {
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_set_rm_attrs");
  return fn(x);
}

inline SEXP
exclude_locs(SEXP exclude, R_xlen_t xn) {
  typedef SEXP fn_t(SEXP, R_xlen_t);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_exclude_locs");
  return fn(exclude, xn);
}

inline R_xlen_t
unnested_length(SEXP x){
  typedef R_xlen_t fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_unnested_length");
  return fn(x);
}

inline SEXP
drop_null(SEXP list_of_vecs, bool always_shallow_copy) {
  typedef SEXP fn_t(SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_drop_null");
  return fn(list_of_vecs, always_shallow_copy);
}

inline SEXP
lengths(SEXP list_of_vecs, bool names){
  typedef SEXP fn_t(SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_lengths");
  return fn(list_of_vecs, names);
}

inline SEXP
new_list(R_xlen_t length, SEXP default_value){
  typedef SEXP fn_t(R_xlen_t, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_new_list");
  return fn(length, default_value);
}

inline SEXP
list_assign(SEXP x, SEXP values){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_list_assign");
  return fn(x, values);
}

inline SEXP
sset(SEXP x, SEXP indices, bool check){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_sset");
  return fn(x, indices, check);
}

inline SEXP
sset_vec(SEXP x, SEXP indices, bool check){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_sset_vec");
  return fn(x, indices, check);
}

inline SEXP
val_find(SEXP x, SEXP value, bool invert){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_val_find");
  return fn(x, value, invert);
}

inline SEXP
val_remove(SEXP x, SEXP value){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_val_remove");
  return fn(x, value);
}

inline SEXP
sequence(SEXP size, SEXP from, SEXP by, bool as_list, bool add_id){
  typedef SEXP fn_t(SEXP, SEXP, SEXP, bool, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_sequence");
  return fn(size, from, by, as_list, add_id);
}

inline SEXP
seq_len(R_xlen_t n){
  typedef SEXP fn_t(R_xlen_t);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_seq_len");
  return fn(n);
}

inline SEXP
rep_len(SEXP x, int length){
  typedef SEXP fn_t(SEXP, int);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_rep_len");
  return fn(x, length);
}

inline SEXP
rep(SEXP x, SEXP times){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_rep");
  return fn(x, times);
}

inline SEXP
rep_each(SEXP x, SEXP each){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_rep_each");
  return fn(x, each);
}

inline SEXP
recycle(SEXP list_of_vecs, SEXP length){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_recycle");
  return fn(list_of_vecs, length);
}

inline SEXP
c(SEXP list_of_vecs){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_c");
  return fn(list_of_vecs);
}

inline SEXP
list_c(SEXP list_of_vecs){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_list_c");
  return fn(list_of_vecs);
}

inline SEXP
name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep){
  typedef SEXP fn_t(SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_name_repair");
  return fn(names, dup_sep, empty_sep);
}

inline SEXP
unique(SEXP x, bool names){
  typedef SEXP fn_t(SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_unique");
  return fn(x, names);
}

inline SEXP
setdiff(SEXP x, SEXP y, bool unique){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_setdiff");
  return fn(x, y, unique);
}

inline SEXP
intersect(SEXP x, SEXP y, bool unique){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_intersect");
  return fn(x, y, unique);
}

inline SEXP
get_ptype(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_get_ptype");
  return fn(x);
}

// Data frame functions

inline SEXP
create_df_row_names(int x) {
  typedef SEXP fn_t(int);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_create_df_row_names");
  return fn(x);
}

inline SEXP
list_as_df(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_list_as_df");
  return fn(x);
}

inline SEXP
new_df(SEXP list_of_vecs, SEXP nrows, bool recycle, bool name_repair){
  typedef SEXP fn_t(SEXP, SEXP, bool, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_new_df");
  return fn(list_of_vecs, nrows, recycle, name_repair);
}

inline SEXP
df_slice(SEXP x, SEXP indices, bool check){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_slice");
  return fn(x, indices, check);
}

inline SEXP
df_select(SEXP x, SEXP locs){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_select");
  return fn(x, locs);
}

inline SEXP
df_subset(SEXP x, SEXP i, SEXP j, bool check){
  typedef SEXP fn_t(SEXP, SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_subset");
  return fn(x, i, j, check);
}

inline SEXP
df_assign_cols(SEXP x, SEXP cols){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_assign_cols");
  return fn(x, cols);
}

inline SEXP
df_col_c(SEXP list_of_vecs, bool recycle, bool name_repair){
  typedef SEXP fn_t(SEXP, bool, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_col_c");
  return fn(list_of_vecs, recycle, name_repair);
}

inline SEXP
str_coalesce(SEXP list_of_vecs){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_str_coalesce");
  return fn(list_of_vecs);
}

inline SEXP
rebuild(SEXP x, SEXP source, bool shallow_copy){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_rebuild");
  return fn(x, source, shallow_copy);
}

inline SEXP
semi_copy(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_semi_copy");
  return fn(x);
}

inline SEXP
paste(SEXP x, SEXP sep, SEXP collapse){
  typedef SEXP fn_t(SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_paste");
  return fn(x, sep, collapse);
}

inline void
replace(SEXP x, SEXP where, SEXP with, bool quiet){
  typedef void fn_t(SEXP, SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_replace");
  return fn(x, where, with, quiet);
}

inline SEXP
if_else(SEXP condition, SEXP yes, SEXP no, SEXP na){
  typedef SEXP fn_t(SEXP, SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_if_else");
  return fn(condition, yes, no, na);
}

inline SEXP
gcd(SEXP x, double tol, bool na_rm){
  typedef SEXP fn_t(SEXP, double, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_gcd");
  return fn(x, tol, na_rm);
}

inline SEXP
clean_indices(SEXP locs, SEXP x){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_clean_indices");
  return fn(locs, x);
}

// Deprecated, use cheapr::vector_length
inline R_xlen_t
vec_length(SEXP x){
  return cheapr::vector_length(x);
}

// Deprecated, use cheapr::address
inline SEXP
r_address(SEXP x){
  return cheapr::address(x);
}

// Deprecated, use cheapr::replace
inline SEXP
loc_set_replace(SEXP x, SEXP where, SEXP what){
  typedef SEXP fn_t(SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_loc_set_replace");
  return fn(x, where, what);
}

// R vector constructor from C++ types and SEXP
// This is also defined in variadic.h but difficult to export with one definition
template<typename... Args>
SEXP new_r_vec(Args... args) {
  SEXP out = cheapr::SHIELD(cheapr::new_r_list(args...));
  cheapr::SHIELD(out = cheapr::c(out));
  cheapr::YIELD(2);
  return out;
}

template<typename... Args>
SEXP new_r_df(Args... args) {
  SEXP out = cheapr::SHIELD(cheapr::new_r_list(args...));
  cheapr::SHIELD(out = new_df(out, R_NilValue, true, true));
  cheapr::YIELD(2);
  return out;
}

template<typename... Args>
inline SEXP r_paste(SEXP sep, SEXP collapse, Args... args){
  SEXP objs = cheapr::SHIELD(cheapr::new_r_list(args...));
  SEXP out = cheapr::SHIELD(cheapr::paste(objs, sep, collapse));
  cheapr::YIELD(2);
  return out;
}

}

// -----------------------------------------------------------------------------

#endif
