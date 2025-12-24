#ifndef CHEAPR_C_API_H
#define CHEAPR_C_API_H

#if defined(CHEAPR_API_SELECTED)
#if CHEAPR_API_SELECTED != 2
#error "Only one cheapr API header may be included"
#endif
#else
#define CHEAPR_API_SELECTED 2 // C API
#endif

// Core constants and small functions
#include <cheapr/internal/c_core.h>
// R tools for exporting C fns
#include <R_ext/Rdynload.h>

// -----------------------------------------------------------------------------

namespace cheapr {


namespace altrep {
inline bool
is_compact_seq(SEXP x) {
  typedef bool fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_is_compact_seq");
  return fn(x);
}
}


namespace lst {
inline R_xlen_t
unnested_length(SEXP x){
  typedef R_xlen_t fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_unnested_length");
  return fn(x);
}

inline SEXP
drop_null(SEXP list_of_vecs) {
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_drop_null");
  return fn(list_of_vecs);
}


inline SEXP
lengths(SEXP list_of_vecs, bool names){
  typedef SEXP fn_t(SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_lengths");
  return fn(list_of_vecs, names);
}

inline SEXP
modify(SEXP x, SEXP values){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_list_assign");
  return fn(x, values);
}
}

namespace internal {
inline SEXP
sset_vec(SEXP x, SEXP indices, bool check){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_sset_vec");
  return fn(x, indices, check);
}
}

namespace vec {

inline SEXP
sset(SEXP x, SEXP indices){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_sset");
  return fn(x, indices, true);
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
rep_len(SEXP x, R_xlen_t length){
  typedef SEXP fn_t(SEXP, R_xlen_t);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_rep_len2");
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
semi_copy(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_semi_copy");
  return fn(x);
}

}

inline SEXP
name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep){
  typedef SEXP fn_t(SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_name_repair");
  return fn(names, dup_sep, empty_sep);
}

namespace df {

inline SEXP
new_df(SEXP list_of_vecs, SEXP nrows, bool recycle, bool name_repair){
  typedef SEXP fn_t(SEXP, SEXP, bool, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_new_df");
  return fn(list_of_vecs, nrows, recycle, name_repair);
}

// Data frame functions

template<typename... Args>
SEXP make_df(Args... args) {
  SEXP out = SHIELD(cheapr::vec::make_list(args...));
  SHIELD(out = new_df(out, r_null, true, true));
  YIELD(2);
  return out;
}

inline SEXP
list_as_df(SEXP x){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_list_as_df");
  return fn(x);
}


inline SEXP
slice(SEXP x, SEXP indices, bool check){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_slice");
  return fn(x, indices, check);
}

inline SEXP
select(SEXP x, SEXP locs){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_select");
  return fn(x, locs);
}

inline SEXP
modify(SEXP x, SEXP cols){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_assign_cols");
  return fn(x, cols);
}
}

inline SEXP
rebuild(SEXP x, SEXP source, bool shallow_copy){
  typedef SEXP fn_t(SEXP, SEXP, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_rebuild");
  return fn(x, source, shallow_copy);
}

inline SEXP
clean_indices(SEXP locs, SEXP x){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_clean_indices");
  return fn(locs, x);
}


namespace collapse {

inline SEXP
list_combine(SEXP list_of_vecs){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_list_c");
  return fn(list_of_vecs);
}

inline SEXP
str_coalesce(SEXP list_of_vecs){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_str_coalesce");
  return fn(list_of_vecs);
}

inline SEXP
recycle(SEXP list_of_vecs, SEXP length){
  typedef SEXP fn_t(SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_recycle");
  return fn(list_of_vecs, length);
}

inline SEXP
combine(SEXP list_of_vecs){
  typedef SEXP fn_t(SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_c");
  return fn(list_of_vecs);
}
inline SEXP
paste(SEXP list_of_chars, SEXP sep, SEXP collapse){
  typedef SEXP fn_t(SEXP, SEXP, SEXP);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_paste");
  return fn(list_of_chars, sep, collapse);
}
inline SEXP
col_combine(SEXP list_of_vecs, bool recycle, bool name_repair){
  typedef SEXP fn_t(SEXP, bool, bool);
  static fn_t *fn = (fn_t*) R_GetCCallable("cheapr", "api_df_col_c");
  return fn(list_of_vecs, recycle, name_repair);
}
}

namespace vec {
// R vector constructor from C++ types and SEXP
// This is also defined in variadic.h but difficult to export with one definition
template<typename... Args>
SEXP combine(Args... args) {
  SEXP out = SHIELD(cheapr::vec::make_list(args...));
  SHIELD(out = cheapr::collapse::combine(out));
  YIELD(2);
  return out;
}

template<typename... Args>
SEXP make_vec(Args... args) {
  return combine(args...);
}

template<typename... Args>
inline SEXP r_paste(SEXP sep, SEXP collapse, Args... args){
  SEXP objs = SHIELD(cheapr::vec::make_list(args...));
  SEXP out = SHIELD(cheapr::collapse::paste(objs, sep, collapse));
  YIELD(2);
  return out;
}

template<typename... Args>
SEXP recycle(Args... args){
  SEXP out = SHIELD(cheapr::vec::make_list(args...));
  SHIELD(out = collapse::recycle(out, r_null));
  YIELD(2);
  return out;
}

}

}

// -----------------------------------------------------------------------------

#endif
