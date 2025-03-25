// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// altrep.cpp
bool is_compact_seq(SEXP x);
extern "C" SEXP _cheapr_is_compact_seq(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(is_compact_seq(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// altrep.cpp
SEXP compact_seq_data(SEXP x);
extern "C" SEXP _cheapr_compact_seq_data(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(compact_seq_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// attrs.cpp
SEXP cpp_set_rm_attributes(SEXP x);
extern "C" SEXP _cheapr_cpp_set_rm_attributes(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_rm_attributes(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// attrs.cpp
SEXP cpp_set_add_attr(SEXP x, SEXP which, SEXP value);
extern "C" SEXP _cheapr_cpp_set_add_attr(SEXP x, SEXP which, SEXP value) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_add_attr(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(which), cpp11::as_cpp<cpp11::decay_t<SEXP>>(value)));
  END_CPP11
}
// attrs.cpp
SEXP cpp_set_rm_attr(SEXP x, SEXP which);
extern "C" SEXP _cheapr_cpp_set_rm_attr(SEXP x, SEXP which) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_rm_attr(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(which)));
  END_CPP11
}
// attrs.cpp
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add);
extern "C" SEXP _cheapr_cpp_set_add_attributes(SEXP x, SEXP attributes, SEXP add) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_add_attributes(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(attributes), cpp11::as_cpp<cpp11::decay_t<bool>>(add)));
  END_CPP11
}
// attrs.cpp
void cpp_shallow_duplicate_attrs(SEXP source, SEXP target);
extern "C" SEXP _cheapr_cpp_shallow_duplicate_attrs(SEXP source, SEXP target) {
  BEGIN_CPP11
    cpp_shallow_duplicate_attrs(cpp11::as_cpp<cpp11::decay_t<SEXP>>(source), cpp11::as_cpp<cpp11::decay_t<SEXP>>(target));
    return R_NilValue;
  END_CPP11
}
// attrs.cpp
void cpp_copy_most_attrs(SEXP source, SEXP target);
extern "C" SEXP _cheapr_cpp_copy_most_attrs(SEXP source, SEXP target) {
  BEGIN_CPP11
    cpp_copy_most_attrs(cpp11::as_cpp<cpp11::decay_t<SEXP>>(source), cpp11::as_cpp<cpp11::decay_t<SEXP>>(target));
    return R_NilValue;
  END_CPP11
}
// combine.cpp
SEXP cpp_rep_len(SEXP x, int length);
extern "C" SEXP _cheapr_cpp_rep_len(SEXP x, SEXP length) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_rep_len(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<int>>(length)));
  END_CPP11
}
// combine.cpp
SEXP cpp_recycle(SEXP x, SEXP length);
extern "C" SEXP _cheapr_cpp_recycle(SEXP x, SEXP length) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_recycle(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(length)));
  END_CPP11
}
// combine.cpp
SEXP cpp_setdiff(SEXP x, SEXP y);
extern "C" SEXP _cheapr_cpp_setdiff(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_setdiff(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y)));
  END_CPP11
}
// combine.cpp
SEXP get_ptypes(SEXP x);
extern "C" SEXP _cheapr_get_ptypes(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(get_ptypes(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// combine.cpp
SEXP cpp_combine_levels(SEXP x);
extern "C" SEXP _cheapr_cpp_combine_levels(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_combine_levels(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// combine.cpp
SEXP cpp_combine_factors(SEXP x);
extern "C" SEXP _cheapr_cpp_combine_factors(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_combine_factors(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// combine.cpp
SEXP cpp_c(SEXP x);
extern "C" SEXP _cheapr_cpp_c(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_c(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// gcd.cpp
double cpp_gcd2(double x, double y, double tol, bool na_rm);
extern "C" SEXP _cheapr_cpp_gcd2(SEXP x, SEXP y, SEXP tol, SEXP na_rm) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_gcd2(cpp11::as_cpp<cpp11::decay_t<double>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(y), cpp11::as_cpp<cpp11::decay_t<double>>(tol), cpp11::as_cpp<cpp11::decay_t<bool>>(na_rm)));
  END_CPP11
}
// gcd.cpp
double cpp_lcm2(double x, double y, double tol, bool na_rm);
extern "C" SEXP _cheapr_cpp_lcm2(SEXP x, SEXP y, SEXP tol, SEXP na_rm) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lcm2(cpp11::as_cpp<cpp11::decay_t<double>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(y), cpp11::as_cpp<cpp11::decay_t<double>>(tol), cpp11::as_cpp<cpp11::decay_t<bool>>(na_rm)));
  END_CPP11
}
// gcd.cpp
SEXP cpp_gcd(SEXP x, double tol, bool na_rm, bool break_early, bool round);
extern "C" SEXP _cheapr_cpp_gcd(SEXP x, SEXP tol, SEXP na_rm, SEXP break_early, SEXP round) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_gcd(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(tol), cpp11::as_cpp<cpp11::decay_t<bool>>(na_rm), cpp11::as_cpp<cpp11::decay_t<bool>>(break_early), cpp11::as_cpp<cpp11::decay_t<bool>>(round)));
  END_CPP11
}
// gcd.cpp
SEXP cpp_lcm(SEXP x, double tol, bool na_rm);
extern "C" SEXP _cheapr_cpp_lcm(SEXP x, SEXP tol, SEXP na_rm) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lcm(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(tol), cpp11::as_cpp<cpp11::decay_t<bool>>(na_rm)));
  END_CPP11
}
// gcd.cpp
SEXP cpp_gcd2_vectorised(SEXP x, SEXP y, double tol, bool na_rm);
extern "C" SEXP _cheapr_cpp_gcd2_vectorised(SEXP x, SEXP y, SEXP tol, SEXP na_rm) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_gcd2_vectorised(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y), cpp11::as_cpp<cpp11::decay_t<double>>(tol), cpp11::as_cpp<cpp11::decay_t<bool>>(na_rm)));
  END_CPP11
}
// gcd.cpp
SEXP cpp_lcm2_vectorised(SEXP x, SEXP y, double tol, bool na_rm);
extern "C" SEXP _cheapr_cpp_lcm2_vectorised(SEXP x, SEXP y, SEXP tol, SEXP na_rm) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lcm2_vectorised(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y), cpp11::as_cpp<cpp11::decay_t<double>>(tol), cpp11::as_cpp<cpp11::decay_t<bool>>(na_rm)));
  END_CPP11
}
// int64.cpp
SEXP cpp_int64_to_int(SEXP x);
extern "C" SEXP _cheapr_cpp_int64_to_int(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_int64_to_int(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// int64.cpp
SEXP cpp_int64_to_double(SEXP x);
extern "C" SEXP _cheapr_cpp_int64_to_double(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_int64_to_double(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// int64.cpp
SEXP cpp_int64_to_numeric(SEXP x);
extern "C" SEXP _cheapr_cpp_int64_to_numeric(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_int64_to_numeric(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// int64.cpp
SEXP cpp_numeric_to_int64(SEXP x);
extern "C" SEXP _cheapr_cpp_numeric_to_int64(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_numeric_to_int64(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// int64.cpp
SEXP cpp_format_numeric_as_int64(SEXP x);
extern "C" SEXP _cheapr_cpp_format_numeric_as_int64(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_format_numeric_as_int64(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// lag.cpp
SEXP cpp_lag(SEXP x, R_xlen_t k, SEXP fill, bool set, bool recursive);
extern "C" SEXP _cheapr_cpp_lag(SEXP x, SEXP k, SEXP fill, SEXP set, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lag(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<R_xlen_t>>(k), cpp11::as_cpp<cpp11::decay_t<SEXP>>(fill), cpp11::as_cpp<cpp11::decay_t<bool>>(set), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// lag.cpp
SEXP cpp_lag2(SEXP x, SEXP lag, SEXP order, SEXP run_lengths, SEXP fill, bool recursive);
extern "C" SEXP _cheapr_cpp_lag2(SEXP x, SEXP lag, SEXP order, SEXP run_lengths, SEXP fill, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lag2(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(lag), cpp11::as_cpp<cpp11::decay_t<SEXP>>(order), cpp11::as_cpp<cpp11::decay_t<SEXP>>(run_lengths), cpp11::as_cpp<cpp11::decay_t<SEXP>>(fill), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// lists.cpp
SEXP cpp_unnested_length(SEXP x);
extern "C" SEXP _cheapr_cpp_unnested_length(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_unnested_length(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// lists.cpp
SEXP cpp_lengths(SEXP x, bool names);
extern "C" SEXP _cheapr_cpp_lengths(SEXP x, SEXP names) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lengths(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(names)));
  END_CPP11
}
// lists.cpp
SEXP cpp_new_list(SEXP size, SEXP default_value);
extern "C" SEXP _cheapr_cpp_new_list(SEXP size, SEXP default_value) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_new_list(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size), cpp11::as_cpp<cpp11::decay_t<SEXP>>(default_value)));
  END_CPP11
}
// lists.cpp
SEXP shallow_copy(SEXP x);
extern "C" SEXP _cheapr_shallow_copy(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(shallow_copy(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// lists.cpp
SEXP cpp_drop_null(SEXP l, bool always_shallow_copy);
extern "C" SEXP _cheapr_cpp_drop_null(SEXP l, SEXP always_shallow_copy) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_drop_null(cpp11::as_cpp<cpp11::decay_t<SEXP>>(l), cpp11::as_cpp<cpp11::decay_t<bool>>(always_shallow_copy)));
  END_CPP11
}
// lists.cpp
SEXP cpp_list_assign(SEXP x, SEXP values);
extern "C" SEXP _cheapr_cpp_list_assign(SEXP x, SEXP values) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_list_assign(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(values)));
  END_CPP11
}
// lists.cpp
SEXP cpp_list_as_df(SEXP x);
extern "C" SEXP _cheapr_cpp_list_as_df(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_list_as_df(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// lists.cpp
SEXP cpp_df_assign_cols(SEXP x, SEXP cols);
extern "C" SEXP _cheapr_cpp_df_assign_cols(SEXP x, SEXP cols) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_df_assign_cols(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(cols)));
  END_CPP11
}
// lists.cpp
SEXP cpp_df_reconstruct(SEXP data, SEXP from, bool keep_attrs);
extern "C" SEXP _cheapr_cpp_df_reconstruct(SEXP data, SEXP from, SEXP keep_attrs) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_df_reconstruct(cpp11::as_cpp<cpp11::decay_t<SEXP>>(data), cpp11::as_cpp<cpp11::decay_t<SEXP>>(from), cpp11::as_cpp<cpp11::decay_t<bool>>(keep_attrs)));
  END_CPP11
}
// nas.cpp
SEXP cpp_num_na(SEXP x, bool recursive);
extern "C" SEXP _cheapr_cpp_num_na(SEXP x, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_num_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// nas.cpp
bool cpp_any_na(SEXP x, bool recursive);
extern "C" SEXP _cheapr_cpp_any_na(SEXP x, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_any_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// nas.cpp
bool cpp_all_na(SEXP x, bool return_true_on_empty, bool recursive);
extern "C" SEXP _cheapr_cpp_all_na(SEXP x, SEXP return_true_on_empty, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_all_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(return_true_on_empty), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// nas.cpp
SEXP cpp_is_na(SEXP x);
extern "C" SEXP _cheapr_cpp_is_na(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_is_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_df_row_na_counts(SEXP x);
extern "C" SEXP _cheapr_cpp_df_row_na_counts(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_df_row_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_df_col_na_counts(SEXP x);
extern "C" SEXP _cheapr_cpp_df_col_na_counts(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_df_col_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_col_any_na(SEXP x, bool names);
extern "C" SEXP _cheapr_cpp_col_any_na(SEXP x, SEXP names) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_col_any_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(names)));
  END_CPP11
}
// nas.cpp
SEXP cpp_col_all_na(SEXP x, bool names);
extern "C" SEXP _cheapr_cpp_col_all_na(SEXP x, SEXP names) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_col_all_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(names)));
  END_CPP11
}
// nas.cpp
SEXP cpp_matrix_row_na_counts(SEXP x);
extern "C" SEXP _cheapr_cpp_matrix_row_na_counts(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_matrix_row_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_matrix_col_na_counts(SEXP x);
extern "C" SEXP _cheapr_cpp_matrix_col_na_counts(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_matrix_col_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_row_na_counts(SEXP x, bool names);
extern "C" SEXP _cheapr_cpp_row_na_counts(SEXP x, SEXP names) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_row_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(names)));
  END_CPP11
}
// nas.cpp
SEXP cpp_col_na_counts(SEXP x, bool names);
extern "C" SEXP _cheapr_cpp_col_na_counts(SEXP x, SEXP names) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_col_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(names)));
  END_CPP11
}
// scalars.cpp
SEXP cpp_count_val(SEXP x, SEXP value, bool recursive);
extern "C" SEXP _cheapr_cpp_count_val(SEXP x, SEXP value, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_count_val(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(value), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// scalars.cpp
SEXP cpp_val_replace(SEXP x, SEXP value, SEXP replace, bool recursive);
extern "C" SEXP _cheapr_cpp_val_replace(SEXP x, SEXP value, SEXP replace, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_val_replace(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(value), cpp11::as_cpp<cpp11::decay_t<SEXP>>(replace), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// scalars.cpp
SEXP cpp_val_set_replace(SEXP x, SEXP value, SEXP replace, bool recursive);
extern "C" SEXP _cheapr_cpp_val_set_replace(SEXP x, SEXP value, SEXP replace, SEXP recursive) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_val_set_replace(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(value), cpp11::as_cpp<cpp11::decay_t<SEXP>>(replace), cpp11::as_cpp<cpp11::decay_t<bool>>(recursive)));
  END_CPP11
}
// scalars.cpp
SEXP cpp_loc_set_replace(SEXP x, SEXP where, SEXP what);
extern "C" SEXP _cheapr_cpp_loc_set_replace(SEXP x, SEXP where, SEXP what) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_loc_set_replace(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(where), cpp11::as_cpp<cpp11::decay_t<SEXP>>(what)));
  END_CPP11
}
// scalars.cpp
SEXP cpp_val_remove(SEXP x, SEXP value);
extern "C" SEXP _cheapr_cpp_val_remove(SEXP x, SEXP value) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_val_remove(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(value)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_int_sequence(SEXP size, SEXP from, SEXP by);
extern "C" SEXP _cheapr_cpp_int_sequence(SEXP size, SEXP from, SEXP by) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_int_sequence(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size), cpp11::as_cpp<cpp11::decay_t<SEXP>>(from), cpp11::as_cpp<cpp11::decay_t<SEXP>>(by)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_dbl_sequence(SEXP size, SEXP from, SEXP by);
extern "C" SEXP _cheapr_cpp_dbl_sequence(SEXP size, SEXP from, SEXP by) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_dbl_sequence(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size), cpp11::as_cpp<cpp11::decay_t<SEXP>>(from), cpp11::as_cpp<cpp11::decay_t<SEXP>>(by)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_sequence(SEXP size, SEXP from, SEXP by);
extern "C" SEXP _cheapr_cpp_sequence(SEXP size, SEXP from, SEXP by) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_sequence(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size), cpp11::as_cpp<cpp11::decay_t<SEXP>>(from), cpp11::as_cpp<cpp11::decay_t<SEXP>>(by)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_window_sequence(SEXP size, double k, bool partial, bool ascending);
extern "C" SEXP _cheapr_cpp_window_sequence(SEXP size, SEXP k, SEXP partial, SEXP ascending) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_window_sequence(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(k), cpp11::as_cpp<cpp11::decay_t<bool>>(partial), cpp11::as_cpp<cpp11::decay_t<bool>>(ascending)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_lag_sequence(SEXP size, double k, bool partial);
extern "C" SEXP _cheapr_cpp_lag_sequence(SEXP size, SEXP k, SEXP partial) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lag_sequence(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(k), cpp11::as_cpp<cpp11::decay_t<bool>>(partial)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_lead_sequence(SEXP size, double k, bool partial);
extern "C" SEXP _cheapr_cpp_lead_sequence(SEXP size, SEXP k, SEXP partial) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lead_sequence(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size), cpp11::as_cpp<cpp11::decay_t<double>>(k), cpp11::as_cpp<cpp11::decay_t<bool>>(partial)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_sequence_id(SEXP size);
extern "C" SEXP _cheapr_cpp_sequence_id(SEXP size) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_sequence_id(cpp11::as_cpp<cpp11::decay_t<SEXP>>(size)));
  END_CPP11
}
// sequences.cpp
SEXP cpp_fixed_width_breaks(double start, double end, double n, bool pretty, bool expand_min, bool expand_max);
extern "C" SEXP _cheapr_cpp_fixed_width_breaks(SEXP start, SEXP end, SEXP n, SEXP pretty, SEXP expand_min, SEXP expand_max) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_fixed_width_breaks(cpp11::as_cpp<cpp11::decay_t<double>>(start), cpp11::as_cpp<cpp11::decay_t<double>>(end), cpp11::as_cpp<cpp11::decay_t<double>>(n), cpp11::as_cpp<cpp11::decay_t<bool>>(pretty), cpp11::as_cpp<cpp11::decay_t<bool>>(expand_min), cpp11::as_cpp<cpp11::decay_t<bool>>(expand_max)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_abs(SEXP x);
extern "C" SEXP _cheapr_cpp_set_abs(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_abs(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_floor(SEXP x);
extern "C" SEXP _cheapr_cpp_set_floor(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_floor(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_ceiling(SEXP x);
extern "C" SEXP _cheapr_cpp_set_ceiling(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_ceiling(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_trunc(SEXP x);
extern "C" SEXP _cheapr_cpp_set_trunc(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_trunc(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_change_sign(SEXP x);
extern "C" SEXP _cheapr_cpp_set_change_sign(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_change_sign(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_exp(SEXP x);
extern "C" SEXP _cheapr_cpp_set_exp(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_exp(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_sqrt(SEXP x);
extern "C" SEXP _cheapr_cpp_set_sqrt(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_sqrt(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_add(SEXP x, SEXP y);
extern "C" SEXP _cheapr_cpp_set_add(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_add(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_subtract(SEXP x, SEXP y);
extern "C" SEXP _cheapr_cpp_set_subtract(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_subtract(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_multiply(SEXP x, SEXP y);
extern "C" SEXP _cheapr_cpp_set_multiply(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_multiply(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_divide(SEXP x, SEXP y);
extern "C" SEXP _cheapr_cpp_set_divide(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_divide(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_pow(SEXP x, SEXP y);
extern "C" SEXP _cheapr_cpp_set_pow(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_pow(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_log(SEXP x, SEXP base);
extern "C" SEXP _cheapr_cpp_set_log(SEXP x, SEXP base) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_log(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(base)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_set_round(SEXP x, SEXP digits);
extern "C" SEXP _cheapr_cpp_set_round(SEXP x, SEXP digits) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_round(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(digits)));
  END_CPP11
}
// set_math.cpp
SEXP cpp_int_sign(SEXP x);
extern "C" SEXP _cheapr_cpp_int_sign(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_int_sign(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// sset.cpp
SEXP clean_indices(SEXP indices, SEXP x);
extern "C" SEXP _cheapr_clean_indices(SEXP indices, SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(clean_indices(cpp11::as_cpp<cpp11::decay_t<SEXP>>(indices), cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// sset.cpp
SEXP cpp_sset(SEXP x, SEXP indices);
extern "C" SEXP _cheapr_cpp_sset(SEXP x, SEXP indices) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_sset(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(indices)));
  END_CPP11
}
// sset.cpp
SEXP cpp_rev(SEXP x, bool set);
extern "C" SEXP _cheapr_cpp_rev(SEXP x, SEXP set) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_rev(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(set)));
  END_CPP11
}
// sset.cpp
SEXP cpp_df_select(SEXP x, SEXP locs);
extern "C" SEXP _cheapr_cpp_df_select(SEXP x, SEXP locs) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_df_select(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(locs)));
  END_CPP11
}
// sset.cpp
SEXP cpp_df_slice(SEXP x, SEXP indices, bool check);
extern "C" SEXP _cheapr_cpp_df_slice(SEXP x, SEXP indices, SEXP check) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_df_slice(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(indices), cpp11::as_cpp<cpp11::decay_t<bool>>(check)));
  END_CPP11
}
// sset.cpp
SEXP cpp_df_subset(SEXP x, SEXP i, SEXP j, bool keep_attrs);
extern "C" SEXP _cheapr_cpp_df_subset(SEXP x, SEXP i, SEXP j, SEXP keep_attrs) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_df_subset(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(i), cpp11::as_cpp<cpp11::decay_t<SEXP>>(j), cpp11::as_cpp<cpp11::decay_t<bool>>(keep_attrs)));
  END_CPP11
}
// utils.cpp
SEXP cpp_is_simple_atomic_vec(SEXP x);
extern "C" SEXP _cheapr_cpp_is_simple_atomic_vec(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_is_simple_atomic_vec(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
SEXP cpp_is_simple_vec(SEXP x);
extern "C" SEXP _cheapr_cpp_is_simple_vec(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_is_simple_vec(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
SEXP cpp_vector_length(SEXP x);
extern "C" SEXP _cheapr_cpp_vector_length(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_vector_length(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
SEXP cpp_address(SEXP x);
extern "C" SEXP _cheapr_cpp_address(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_address(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
SEXP r_copy(SEXP x);
extern "C" SEXP _cheapr_r_copy(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(r_copy(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
double var_sum_squared_diff(SEXP x, double mu);
extern "C" SEXP _cheapr_var_sum_squared_diff(SEXP x, SEXP mu) {
  BEGIN_CPP11
    return cpp11::as_sexp(var_sum_squared_diff(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(mu)));
  END_CPP11
}
// utils.cpp
SEXP cpp_bin(SEXP x, SEXP breaks, bool codes, bool right, bool include_lowest, bool include_oob);
extern "C" SEXP _cheapr_cpp_bin(SEXP x, SEXP breaks, SEXP codes, SEXP right, SEXP include_lowest, SEXP include_oob) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_bin(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(breaks), cpp11::as_cpp<cpp11::decay_t<bool>>(codes), cpp11::as_cpp<cpp11::decay_t<bool>>(right), cpp11::as_cpp<cpp11::decay_t<bool>>(include_lowest), cpp11::as_cpp<cpp11::decay_t<bool>>(include_oob)));
  END_CPP11
}
// utils.cpp
SEXP cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na);
extern "C" SEXP _cheapr_cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_if_else(cpp11::as_cpp<cpp11::decay_t<SEXP>>(condition), cpp11::as_cpp<cpp11::decay_t<SEXP>>(yes), cpp11::as_cpp<cpp11::decay_t<SEXP>>(no), cpp11::as_cpp<cpp11::decay_t<SEXP>>(na)));
  END_CPP11
}
// utils.cpp
SEXP cpp_lgl_count(SEXP x);
extern "C" SEXP _cheapr_cpp_lgl_count(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lgl_count(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
void cpp_set_copy_elements(SEXP source, SEXP target);
extern "C" SEXP _cheapr_cpp_set_copy_elements(SEXP source, SEXP target) {
  BEGIN_CPP11
    cpp_set_copy_elements(cpp11::as_cpp<cpp11::decay_t<SEXP>>(source), cpp11::as_cpp<cpp11::decay_t<SEXP>>(target));
    return R_NilValue;
  END_CPP11
}
// utils.cpp
SEXP cpp_set_or(SEXP x, SEXP y);
extern "C" SEXP _cheapr_cpp_set_or(SEXP x, SEXP y) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_set_or(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(y)));
  END_CPP11
}
// utils.cpp
SEXP cpp_growth_rate(SEXP x);
extern "C" SEXP _cheapr_cpp_growth_rate(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_growth_rate(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
SEXP cpp_name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep);
extern "C" SEXP _cheapr_cpp_name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_name_repair(cpp11::as_cpp<cpp11::decay_t<SEXP>>(names), cpp11::as_cpp<cpp11::decay_t<SEXP>>(dup_sep), cpp11::as_cpp<cpp11::decay_t<SEXP>>(empty_sep)));
  END_CPP11
}
// which.cpp
SEXP cpp_which_(SEXP x, bool invert);
extern "C" SEXP _cheapr_cpp_which_(SEXP x, SEXP invert) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_which_(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(invert)));
  END_CPP11
}
// which.cpp
SEXP cpp_which_val(SEXP x, SEXP value, bool invert);
extern "C" SEXP _cheapr_cpp_which_val(SEXP x, SEXP value, SEXP invert) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_which_val(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<SEXP>>(value), cpp11::as_cpp<cpp11::decay_t<bool>>(invert)));
  END_CPP11
}
// which.cpp
SEXP cpp_which_na(SEXP x);
extern "C" SEXP _cheapr_cpp_which_na(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_which_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// which.cpp
SEXP cpp_which_not_na(SEXP x);
extern "C" SEXP _cheapr_cpp_which_not_na(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_which_not_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// which.cpp
SEXP cpp_lgl_locs(SEXP x, R_xlen_t n_true, R_xlen_t n_false, bool include_true, bool include_false, bool include_na);
extern "C" SEXP _cheapr_cpp_lgl_locs(SEXP x, SEXP n_true, SEXP n_false, SEXP include_true, SEXP include_false, SEXP include_na) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lgl_locs(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<R_xlen_t>>(n_true), cpp11::as_cpp<cpp11::decay_t<R_xlen_t>>(n_false), cpp11::as_cpp<cpp11::decay_t<bool>>(include_true), cpp11::as_cpp<cpp11::decay_t<bool>>(include_false), cpp11::as_cpp<cpp11::decay_t<bool>>(include_na)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_cheapr_clean_indices",               (DL_FUNC) &_cheapr_clean_indices,               2},
    {"_cheapr_compact_seq_data",            (DL_FUNC) &_cheapr_compact_seq_data,            1},
    {"_cheapr_cpp_address",                 (DL_FUNC) &_cheapr_cpp_address,                 1},
    {"_cheapr_cpp_all_na",                  (DL_FUNC) &_cheapr_cpp_all_na,                  3},
    {"_cheapr_cpp_any_na",                  (DL_FUNC) &_cheapr_cpp_any_na,                  2},
    {"_cheapr_cpp_bin",                     (DL_FUNC) &_cheapr_cpp_bin,                     6},
    {"_cheapr_cpp_c",                       (DL_FUNC) &_cheapr_cpp_c,                       1},
    {"_cheapr_cpp_col_all_na",              (DL_FUNC) &_cheapr_cpp_col_all_na,              2},
    {"_cheapr_cpp_col_any_na",              (DL_FUNC) &_cheapr_cpp_col_any_na,              2},
    {"_cheapr_cpp_col_na_counts",           (DL_FUNC) &_cheapr_cpp_col_na_counts,           2},
    {"_cheapr_cpp_combine_factors",         (DL_FUNC) &_cheapr_cpp_combine_factors,         1},
    {"_cheapr_cpp_combine_levels",          (DL_FUNC) &_cheapr_cpp_combine_levels,          1},
    {"_cheapr_cpp_copy_most_attrs",         (DL_FUNC) &_cheapr_cpp_copy_most_attrs,         2},
    {"_cheapr_cpp_count_val",               (DL_FUNC) &_cheapr_cpp_count_val,               3},
    {"_cheapr_cpp_dbl_sequence",            (DL_FUNC) &_cheapr_cpp_dbl_sequence,            3},
    {"_cheapr_cpp_df_assign_cols",          (DL_FUNC) &_cheapr_cpp_df_assign_cols,          2},
    {"_cheapr_cpp_df_col_na_counts",        (DL_FUNC) &_cheapr_cpp_df_col_na_counts,        1},
    {"_cheapr_cpp_df_reconstruct",          (DL_FUNC) &_cheapr_cpp_df_reconstruct,          3},
    {"_cheapr_cpp_df_row_na_counts",        (DL_FUNC) &_cheapr_cpp_df_row_na_counts,        1},
    {"_cheapr_cpp_df_select",               (DL_FUNC) &_cheapr_cpp_df_select,               2},
    {"_cheapr_cpp_df_slice",                (DL_FUNC) &_cheapr_cpp_df_slice,                3},
    {"_cheapr_cpp_df_subset",               (DL_FUNC) &_cheapr_cpp_df_subset,               4},
    {"_cheapr_cpp_drop_null",               (DL_FUNC) &_cheapr_cpp_drop_null,               2},
    {"_cheapr_cpp_fixed_width_breaks",      (DL_FUNC) &_cheapr_cpp_fixed_width_breaks,      6},
    {"_cheapr_cpp_format_numeric_as_int64", (DL_FUNC) &_cheapr_cpp_format_numeric_as_int64, 1},
    {"_cheapr_cpp_gcd",                     (DL_FUNC) &_cheapr_cpp_gcd,                     5},
    {"_cheapr_cpp_gcd2",                    (DL_FUNC) &_cheapr_cpp_gcd2,                    4},
    {"_cheapr_cpp_gcd2_vectorised",         (DL_FUNC) &_cheapr_cpp_gcd2_vectorised,         4},
    {"_cheapr_cpp_growth_rate",             (DL_FUNC) &_cheapr_cpp_growth_rate,             1},
    {"_cheapr_cpp_if_else",                 (DL_FUNC) &_cheapr_cpp_if_else,                 4},
    {"_cheapr_cpp_int64_to_double",         (DL_FUNC) &_cheapr_cpp_int64_to_double,         1},
    {"_cheapr_cpp_int64_to_int",            (DL_FUNC) &_cheapr_cpp_int64_to_int,            1},
    {"_cheapr_cpp_int64_to_numeric",        (DL_FUNC) &_cheapr_cpp_int64_to_numeric,        1},
    {"_cheapr_cpp_int_sequence",            (DL_FUNC) &_cheapr_cpp_int_sequence,            3},
    {"_cheapr_cpp_int_sign",                (DL_FUNC) &_cheapr_cpp_int_sign,                1},
    {"_cheapr_cpp_is_na",                   (DL_FUNC) &_cheapr_cpp_is_na,                   1},
    {"_cheapr_cpp_is_simple_atomic_vec",    (DL_FUNC) &_cheapr_cpp_is_simple_atomic_vec,    1},
    {"_cheapr_cpp_is_simple_vec",           (DL_FUNC) &_cheapr_cpp_is_simple_vec,           1},
    {"_cheapr_cpp_lag",                     (DL_FUNC) &_cheapr_cpp_lag,                     5},
    {"_cheapr_cpp_lag2",                    (DL_FUNC) &_cheapr_cpp_lag2,                    6},
    {"_cheapr_cpp_lag_sequence",            (DL_FUNC) &_cheapr_cpp_lag_sequence,            3},
    {"_cheapr_cpp_lcm",                     (DL_FUNC) &_cheapr_cpp_lcm,                     3},
    {"_cheapr_cpp_lcm2",                    (DL_FUNC) &_cheapr_cpp_lcm2,                    4},
    {"_cheapr_cpp_lcm2_vectorised",         (DL_FUNC) &_cheapr_cpp_lcm2_vectorised,         4},
    {"_cheapr_cpp_lead_sequence",           (DL_FUNC) &_cheapr_cpp_lead_sequence,           3},
    {"_cheapr_cpp_lengths",                 (DL_FUNC) &_cheapr_cpp_lengths,                 2},
    {"_cheapr_cpp_lgl_count",               (DL_FUNC) &_cheapr_cpp_lgl_count,               1},
    {"_cheapr_cpp_lgl_locs",                (DL_FUNC) &_cheapr_cpp_lgl_locs,                6},
    {"_cheapr_cpp_list_as_df",              (DL_FUNC) &_cheapr_cpp_list_as_df,              1},
    {"_cheapr_cpp_list_assign",             (DL_FUNC) &_cheapr_cpp_list_assign,             2},
    {"_cheapr_cpp_loc_set_replace",         (DL_FUNC) &_cheapr_cpp_loc_set_replace,         3},
    {"_cheapr_cpp_matrix_col_na_counts",    (DL_FUNC) &_cheapr_cpp_matrix_col_na_counts,    1},
    {"_cheapr_cpp_matrix_row_na_counts",    (DL_FUNC) &_cheapr_cpp_matrix_row_na_counts,    1},
    {"_cheapr_cpp_name_repair",             (DL_FUNC) &_cheapr_cpp_name_repair,             3},
    {"_cheapr_cpp_new_list",                (DL_FUNC) &_cheapr_cpp_new_list,                2},
    {"_cheapr_cpp_num_na",                  (DL_FUNC) &_cheapr_cpp_num_na,                  2},
    {"_cheapr_cpp_numeric_to_int64",        (DL_FUNC) &_cheapr_cpp_numeric_to_int64,        1},
    {"_cheapr_cpp_recycle",                 (DL_FUNC) &_cheapr_cpp_recycle,                 2},
    {"_cheapr_cpp_rep_len",                 (DL_FUNC) &_cheapr_cpp_rep_len,                 2},
    {"_cheapr_cpp_rev",                     (DL_FUNC) &_cheapr_cpp_rev,                     2},
    {"_cheapr_cpp_row_na_counts",           (DL_FUNC) &_cheapr_cpp_row_na_counts,           2},
    {"_cheapr_cpp_sequence",                (DL_FUNC) &_cheapr_cpp_sequence,                3},
    {"_cheapr_cpp_sequence_id",             (DL_FUNC) &_cheapr_cpp_sequence_id,             1},
    {"_cheapr_cpp_set_abs",                 (DL_FUNC) &_cheapr_cpp_set_abs,                 1},
    {"_cheapr_cpp_set_add",                 (DL_FUNC) &_cheapr_cpp_set_add,                 2},
    {"_cheapr_cpp_set_add_attr",            (DL_FUNC) &_cheapr_cpp_set_add_attr,            3},
    {"_cheapr_cpp_set_add_attributes",      (DL_FUNC) &_cheapr_cpp_set_add_attributes,      3},
    {"_cheapr_cpp_set_ceiling",             (DL_FUNC) &_cheapr_cpp_set_ceiling,             1},
    {"_cheapr_cpp_set_change_sign",         (DL_FUNC) &_cheapr_cpp_set_change_sign,         1},
    {"_cheapr_cpp_set_copy_elements",       (DL_FUNC) &_cheapr_cpp_set_copy_elements,       2},
    {"_cheapr_cpp_set_divide",              (DL_FUNC) &_cheapr_cpp_set_divide,              2},
    {"_cheapr_cpp_set_exp",                 (DL_FUNC) &_cheapr_cpp_set_exp,                 1},
    {"_cheapr_cpp_set_floor",               (DL_FUNC) &_cheapr_cpp_set_floor,               1},
    {"_cheapr_cpp_set_log",                 (DL_FUNC) &_cheapr_cpp_set_log,                 2},
    {"_cheapr_cpp_set_multiply",            (DL_FUNC) &_cheapr_cpp_set_multiply,            2},
    {"_cheapr_cpp_set_or",                  (DL_FUNC) &_cheapr_cpp_set_or,                  2},
    {"_cheapr_cpp_set_pow",                 (DL_FUNC) &_cheapr_cpp_set_pow,                 2},
    {"_cheapr_cpp_set_rm_attr",             (DL_FUNC) &_cheapr_cpp_set_rm_attr,             2},
    {"_cheapr_cpp_set_rm_attributes",       (DL_FUNC) &_cheapr_cpp_set_rm_attributes,       1},
    {"_cheapr_cpp_set_round",               (DL_FUNC) &_cheapr_cpp_set_round,               2},
    {"_cheapr_cpp_set_sqrt",                (DL_FUNC) &_cheapr_cpp_set_sqrt,                1},
    {"_cheapr_cpp_set_subtract",            (DL_FUNC) &_cheapr_cpp_set_subtract,            2},
    {"_cheapr_cpp_set_trunc",               (DL_FUNC) &_cheapr_cpp_set_trunc,               1},
    {"_cheapr_cpp_setdiff",                 (DL_FUNC) &_cheapr_cpp_setdiff,                 2},
    {"_cheapr_cpp_shallow_duplicate_attrs", (DL_FUNC) &_cheapr_cpp_shallow_duplicate_attrs, 2},
    {"_cheapr_cpp_sset",                    (DL_FUNC) &_cheapr_cpp_sset,                    2},
    {"_cheapr_cpp_unnested_length",         (DL_FUNC) &_cheapr_cpp_unnested_length,         1},
    {"_cheapr_cpp_val_remove",              (DL_FUNC) &_cheapr_cpp_val_remove,              2},
    {"_cheapr_cpp_val_replace",             (DL_FUNC) &_cheapr_cpp_val_replace,             4},
    {"_cheapr_cpp_val_set_replace",         (DL_FUNC) &_cheapr_cpp_val_set_replace,         4},
    {"_cheapr_cpp_vector_length",           (DL_FUNC) &_cheapr_cpp_vector_length,           1},
    {"_cheapr_cpp_which_",                  (DL_FUNC) &_cheapr_cpp_which_,                  2},
    {"_cheapr_cpp_which_na",                (DL_FUNC) &_cheapr_cpp_which_na,                1},
    {"_cheapr_cpp_which_not_na",            (DL_FUNC) &_cheapr_cpp_which_not_na,            1},
    {"_cheapr_cpp_which_val",               (DL_FUNC) &_cheapr_cpp_which_val,               3},
    {"_cheapr_cpp_window_sequence",         (DL_FUNC) &_cheapr_cpp_window_sequence,         4},
    {"_cheapr_get_ptypes",                  (DL_FUNC) &_cheapr_get_ptypes,                  1},
    {"_cheapr_is_compact_seq",              (DL_FUNC) &_cheapr_is_compact_seq,              1},
    {"_cheapr_r_copy",                      (DL_FUNC) &_cheapr_r_copy,                      1},
    {"_cheapr_shallow_copy",                (DL_FUNC) &_cheapr_shallow_copy,                1},
    {"_cheapr_var_sum_squared_diff",        (DL_FUNC) &_cheapr_var_sum_squared_diff,        2},
    {NULL, NULL, 0}
};
}

void api_init(DllInfo* dll);

extern "C" attribute_visible void R_init_cheapr(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  api_init(dll);
  R_forceSymbols(dll, TRUE);
}
