#ifndef CHEAPR_DECLARATIONS_H
#define CHEAPR_DECLARATIONS_H

#include "types.h"
#include <cpp11/function.hpp>

// R fns

inline cpp11::function cheapr_sset = cpp11::package("cheapr")["cheapr_sset"];
inline cpp11::function cheapr_is_na = cpp11::package("cheapr")["is_na"];
inline cpp11::function base_rep = cpp11::package("base")["rep"];
inline cpp11::function cheapr_fast_match = cpp11::package("cheapr")["fast_match"];
inline cpp11::function cheapr_fast_unique = cpp11::package("cheapr")["fast_unique"];
inline cpp11::function cheapr_rebuild = cpp11::package("cheapr")["rebuild"];
inline cpp11::function base_assign = cpp11::package("cheapr")["base_assign_at"];
inline cpp11::function base_length = cpp11::package("base")["length"];

// Make function definitions visible to all C++ files
SEXP cpp_which_(SEXP x, bool invert);
SEXP cpp_missing_row(SEXP x, double threshold, bool threshold_is_prop);
SEXP xlen_to_r(R_xlen_t x);
R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive);
SEXP cpp_list_as_df(SEXP x);
SEXP cpp_new_df(SEXP x, SEXP nrows, bool recycle, bool name_repair);
SEXP cpp_is_na(SEXP x);
SEXP cpp_which_na(SEXP x);
SEXP cpp_which_not_na(SEXP x);
SEXP check_transform_altrep(SEXP x);
SEXP altrep_materialise(SEXP x);
SEXP compact_seq_data(SEXP x);
bool is_compact_seq(SEXP x);
R_xlen_t na_count(SEXP x, bool recursive);
bool cpp_any_na(SEXP x, bool recursive);
SEXP cpp_int64_to_double(SEXP x);
SEXP cpp_numeric_to_int64(SEXP x);
SEXP cpp_int64_to_numeric(SEXP x);
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add);
SEXP cpp_set_rm_attributes(SEXP x);
bool implicit_na_coercion(SEXP x, SEXP target);
SEXP cpp_val_find(SEXP x, SEXP value, bool invert, SEXP n_values);
SEXP cpp_set_divide(SEXP x, SEXP y);
SEXP cpp_val_remove(SEXP x, SEXP value);
SEXP cpp_seq_len(R_xlen_t n);
SEXP create_df_row_names(int n);
SEXP cpp_shallow_copy(SEXP x);
SEXP exclude_locs(SEXP exclude, R_xlen_t xn);
R_xlen_t unnested_length(SEXP x);
SEXP cpp_drop_null(SEXP l, bool always_shallow_copy);
SEXP cpp_lengths(SEXP x, bool names);
SEXP sset_vec(SEXP x, SEXP indices, bool check);
SEXP cpp_sset(SEXP x, SEXP indices, bool check);
SEXP cpp_df_slice(SEXP x, SEXP indices, bool check);
SEXP cpp_df_select(SEXP x, SEXP locs);
SEXP cpp_df_subset(SEXP x, SEXP i, SEXP j, bool check);
SEXP cpp_which_val(SEXP x, SEXP value, bool invert);
SEXP cpp_sequence(SEXP size, SEXP from, SEXP by, bool as_list, bool add_id);
SEXP cpp_rep_len(SEXP x, int length);
SEXP cpp_rep(SEXP x, SEXP times);
SEXP cpp_rep_each(SEXP x, SEXP each);
SEXP cpp_recycle(SEXP x, SEXP length);
SEXP cpp_c(SEXP x);
SEXP cpp_list_c(SEXP x);
SEXP cpp_loc_set_replace(SEXP x, SEXP where, SEXP what);
SEXP cpp_name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep);
SEXP cpp_unique(SEXP x, bool names);
SEXP cpp_setdiff(SEXP x, SEXP y, bool unique);
SEXP cpp_intersect(SEXP x, SEXP y, bool unique);
SEXP get_ptype(SEXP x);
SEXP get_list_element(SEXP list, SEXP str);
SEXP rebuild(SEXP x, SEXP source, bool shallow_copy);
SEXP cpp_df_assign_cols(SEXP x, SEXP cols);
SEXP cpp_df_col_c(SEXP x, bool recycle, bool name_repair);
SEXP cpp_list_assign(SEXP x, SEXP values);
SEXP slice_loc(SEXP x, R_xlen_t i);
double cpp_sum(SEXP x);
double cpp_min(SEXP x);
SEXP cpp_str_coalesce(SEXP x);
SEXP cpp_na_init(SEXP x, int n);
SEXP new_list(R_xlen_t length, SEXP default_value);
void set_list_as_df(SEXP x);
SEXP cpp_semi_copy(SEXP x);
void clear_attributes(SEXP x);
uint_fast64_t null_count(SEXP x);
SEXP compact_seq_len(R_xlen_t n);
SEXP clean_indices(SEXP indices, SEXP x, bool count);
SEXP cpp_lgl_count(SEXP x);
SEXP cpp_lgl_locs(SEXP x, R_xlen_t n_true, R_xlen_t n_false,
                  bool include_true, bool include_false, bool include_na);
SEXP cpp_cast_common(SEXP x);
SEXP factor_as_character(SEXP x);
SEXP cpp_val_replace(SEXP x, SEXP value, SEXP replace, bool recursive);
SEXP character_as_factor(SEXP x, SEXP levels);
SEXP match(SEXP y, SEXP x, int no_match);
SEXP cpp_replace(SEXP x, SEXP where, SEXP with, bool in_place, bool quiet);
SEXP clean_locs(SEXP locs, SEXP x);
SEXP cpp_common_template(SEXP x);
SEXP combine_internal(SEXP x, const R_xlen_t out_size, SEXP vec_template);
SEXP cpp_as_df(SEXP x);
void recycle_in_place(SEXP x, R_xlen_t n);
R_xlen_t length_common(SEXP x);
SEXP cpp_paste(SEXP x, SEXP sep, SEXP collapse);
cheapr::r_type r_common_type(SEXP x);
SEXP cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na);
SEXP cpp_gcd(SEXP x, double tol, bool na_rm, bool break_early, bool round);
bool cpp_all_na(SEXP x, bool return_true_on_empty, bool recursive);
SEXP get_vec_names(SEXP x);
SEXP set_vec_names(SEXP x, SEXP names);
bool vec_has_names(SEXP x);

#endif
