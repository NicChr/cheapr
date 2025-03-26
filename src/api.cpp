#include "cheapr.h"
#include <cpp11/R.hpp>
#include <cpp11/protect.hpp>
#include <R_ext/Rdynload.h> // For DllInfo on R 3.3
#include <stdbool.h>

// -----------------------------------------------------------------------------

R_xlen_t
api_vec_length(SEXP x) {
  try {
    return vec_length(x);
  } catch (...) {
    return 0;
  }
}

SEXP
api_r_address(SEXP x) {
  try {
    return r_address(x);
  } catch (...) {
    return R_NilValue;
  }
}

bool
api_is_compact_seq(SEXP x) {
  try {
    return is_compact_seq(x);
  } catch (...) {
    return false;
  }
}

SEXP
api_compact_seq_data(SEXP x) {
  try {
    return compact_seq_data(x);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_create_df_row_names(int n) {
  try {
    return create_df_row_names(n);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_shallow_copy(SEXP x) {
  try {
    return shallow_copy(x);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_set_add_attrs(SEXP x, SEXP attributes, bool add){
  try {
    return cpp_set_add_attributes(x, attributes, add);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_set_rm_attrs(SEXP x) {
  try {
    return cpp_set_rm_attributes(x);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_exclude_locs(SEXP exclude, R_xlen_t xn) {
  try {
    return exclude_locs(exclude, xn);
  } catch (...) {
    return R_NilValue;
  }
}

R_xlen_t
api_unnested_length(SEXP x){
  try {
    return unnested_length(x);
  } catch (...) {
    return (R_xlen_t) 0;
  }
}

SEXP
api_drop_null(SEXP l, bool always_shallow_copy) {
  try {
    return cpp_drop_null(l, always_shallow_copy);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_list_as_df(SEXP x){
  try {
    return cpp_list_as_df(x);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_new_df(SEXP x, SEXP nrows, bool recycle, bool name_repair){
  try {
    return cpp_new_df(x, nrows, recycle, name_repair);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_lengths(SEXP x, bool names){
  try {
    return cpp_lengths(x, names);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_df_slice(SEXP x, SEXP indices, bool check){
  try {
    return cpp_df_slice(x, indices, check);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_df_select(SEXP x, SEXP locs){
  try {
    return cpp_df_select(x, locs);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_sset(SEXP x, SEXP indices){
  try {
    return cpp_sset(x, indices);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_val_find(SEXP x, SEXP value, bool invert){
  try {
    return cpp_which_val(x, value, invert);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_sequence(SEXP size, SEXP from, SEXP by){
  try {
    return cpp_sequence(size, from, by);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_seq_len(R_xlen_t n){
  try {
    return cpp_seq_len(n);
  } catch (...) {
    return R_NilValue;
  }
}

bool
api_is_simple_atomic_vec(SEXP x){
  try {
    return is_simple_atomic_vec(x);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_rep_len(SEXP x, int length){
  try {
    return cpp_rep_len(x, length);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_recycle(SEXP x, SEXP length){
  try {
    return cpp_recycle(x, length);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_c(SEXP x){
  try {
    return cpp_c(x);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep){
  try {
    return cpp_name_repair(names, dup_sep, empty_sep);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_unique(SEXP x){
  try {
    return cpp_unique(x);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_setdiff(SEXP x, SEXP y){
  try {
    return cpp_setdiff(x, y);
  } catch (...) {
    return R_NilValue;
  }
}

SEXP
api_get_ptype(SEXP x){
  try {
    return get_ptype(x);
  } catch (...) {
    return R_NilValue;
  }
}

// -----------------------------------------------------------------------------

[[cpp11::init]]
void api_init(DllInfo* dll){
  R_RegisterCCallable("cheapr", "api_vec_length",    (DL_FUNC)api_vec_length);
  R_RegisterCCallable("cheapr", "api_r_address",    (DL_FUNC)api_r_address);
  R_RegisterCCallable("cheapr", "api_is_compact_seq",    (DL_FUNC)api_is_compact_seq);
  R_RegisterCCallable("cheapr", "api_compact_seq_data",    (DL_FUNC)api_compact_seq_data);
  R_RegisterCCallable("cheapr", "api_create_df_row_names",    (DL_FUNC)api_create_df_row_names);
  R_RegisterCCallable("cheapr", "api_shallow_copy",    (DL_FUNC)api_shallow_copy);
  R_RegisterCCallable("cheapr", "api_set_add_attrs",    (DL_FUNC)api_set_add_attrs);
  R_RegisterCCallable("cheapr", "api_set_rm_attrs",    (DL_FUNC)api_set_rm_attrs);
  R_RegisterCCallable("cheapr", "api_exclude_locs",    (DL_FUNC)api_exclude_locs);
  R_RegisterCCallable("cheapr", "api_unnested_length",    (DL_FUNC)api_unnested_length);
  R_RegisterCCallable("cheapr", "api_drop_null",    (DL_FUNC)api_drop_null);
  R_RegisterCCallable("cheapr", "api_list_as_df",    (DL_FUNC)api_list_as_df);
  R_RegisterCCallable("cheapr", "api_new_df",    (DL_FUNC)api_new_df);
  R_RegisterCCallable("cheapr", "api_lengths",    (DL_FUNC)api_lengths);
  R_RegisterCCallable("cheapr", "api_df_slice",    (DL_FUNC)api_df_slice);
  R_RegisterCCallable("cheapr", "api_df_select",    (DL_FUNC)api_df_select);
  R_RegisterCCallable("cheapr", "api_sset",    (DL_FUNC)api_sset);
  R_RegisterCCallable("cheapr", "api_val_find",    (DL_FUNC)api_val_find);
  R_RegisterCCallable("cheapr", "api_sequence",    (DL_FUNC)api_sequence);
  R_RegisterCCallable("cheapr", "api_seq_len",    (DL_FUNC)api_seq_len);
  R_RegisterCCallable("cheapr", "api_is_simple_atomic_vec",    (DL_FUNC)api_is_simple_atomic_vec);
  R_RegisterCCallable("cheapr", "api_recycle",    (DL_FUNC)api_recycle);
  R_RegisterCCallable("cheapr", "api_rep_len",    (DL_FUNC)api_rep_len);
  R_RegisterCCallable("cheapr", "api_c",    (DL_FUNC)api_c);
  R_RegisterCCallable("cheapr", "api_name_repair",    (DL_FUNC)api_name_repair);
  R_RegisterCCallable("cheapr", "api_unique",    (DL_FUNC)api_unique);
  R_RegisterCCallable("cheapr", "api_setdiff",    (DL_FUNC)api_setdiff);
  R_RegisterCCallable("cheapr", "api_get_ptype",    (DL_FUNC)api_get_ptype);
}
