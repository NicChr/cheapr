// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

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
SEXP cpp_which_na(SEXP x);
extern "C" SEXP _cheapr_cpp_which_na(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_which_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_which_not_na(SEXP x);
extern "C" SEXP _cheapr_cpp_which_not_na(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_which_not_na(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
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
SEXP cpp_row_na_counts(SEXP x);
extern "C" SEXP _cheapr_cpp_row_na_counts(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_row_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_col_na_counts(SEXP x);
extern "C" SEXP _cheapr_cpp_col_na_counts(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_col_na_counts(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// nas.cpp
SEXP cpp_missing_row(SEXP x, double threshold, bool threshold_is_prop);
extern "C" SEXP _cheapr_cpp_missing_row(SEXP x, SEXP threshold, SEXP threshold_is_prop) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_missing_row(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(threshold_is_prop)));
  END_CPP11
}
// nas.cpp
SEXP cpp_missing_col(SEXP x, double threshold, bool threshold_is_prop);
extern "C" SEXP _cheapr_cpp_missing_col(SEXP x, SEXP threshold, SEXP threshold_is_prop) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_missing_col(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(threshold_is_prop)));
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
SEXP cpp_matrix_missing_row(SEXP x, double threshold, bool threshold_is_prop);
extern "C" SEXP _cheapr_cpp_matrix_missing_row(SEXP x, SEXP threshold, SEXP threshold_is_prop) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_matrix_missing_row(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(threshold_is_prop)));
  END_CPP11
}
// nas.cpp
SEXP cpp_matrix_missing_col(SEXP x, double threshold, bool threshold_is_prop);
extern "C" SEXP _cheapr_cpp_matrix_missing_col(SEXP x, SEXP threshold, SEXP threshold_is_prop) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_matrix_missing_col(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<double>>(threshold), cpp11::as_cpp<cpp11::decay_t<bool>>(threshold_is_prop)));
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
// utils.cpp
SEXP cpp_r_unnested_length(SEXP x);
extern "C" SEXP _cheapr_cpp_r_unnested_length(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_r_unnested_length(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
SEXP cpp_lengths(SEXP x);
extern "C" SEXP _cheapr_cpp_lengths(SEXP x) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_lengths(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x)));
  END_CPP11
}
// utils.cpp
SEXP cpp_new_list(R_xlen_t size, SEXP default_value);
extern "C" SEXP _cheapr_cpp_new_list(SEXP size, SEXP default_value) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_new_list(cpp11::as_cpp<cpp11::decay_t<R_xlen_t>>(size), cpp11::as_cpp<cpp11::decay_t<SEXP>>(default_value)));
  END_CPP11
}
// utils.cpp
SEXP cpp_list_rm_null(SEXP l);
extern "C" SEXP _cheapr_cpp_list_rm_null(SEXP l) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_list_rm_null(cpp11::as_cpp<cpp11::decay_t<SEXP>>(l)));
  END_CPP11
}
// which.cpp
SEXP cpp_which_(SEXP x, bool invert);
extern "C" SEXP _cheapr_cpp_which_(SEXP x, SEXP invert) {
  BEGIN_CPP11
    return cpp11::as_sexp(cpp_which_(cpp11::as_cpp<cpp11::decay_t<SEXP>>(x), cpp11::as_cpp<cpp11::decay_t<bool>>(invert)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_cheapr_cpp_all_na",               (DL_FUNC) &_cheapr_cpp_all_na,               3},
    {"_cheapr_cpp_any_na",               (DL_FUNC) &_cheapr_cpp_any_na,               2},
    {"_cheapr_cpp_col_na_counts",        (DL_FUNC) &_cheapr_cpp_col_na_counts,        1},
    {"_cheapr_cpp_dbl_sequence",         (DL_FUNC) &_cheapr_cpp_dbl_sequence,         3},
    {"_cheapr_cpp_gcd",                  (DL_FUNC) &_cheapr_cpp_gcd,                  5},
    {"_cheapr_cpp_gcd2",                 (DL_FUNC) &_cheapr_cpp_gcd2,                 4},
    {"_cheapr_cpp_gcd2_vectorised",      (DL_FUNC) &_cheapr_cpp_gcd2_vectorised,      4},
    {"_cheapr_cpp_int_sequence",         (DL_FUNC) &_cheapr_cpp_int_sequence,         3},
    {"_cheapr_cpp_is_na",                (DL_FUNC) &_cheapr_cpp_is_na,                1},
    {"_cheapr_cpp_lag_sequence",         (DL_FUNC) &_cheapr_cpp_lag_sequence,         3},
    {"_cheapr_cpp_lcm",                  (DL_FUNC) &_cheapr_cpp_lcm,                  3},
    {"_cheapr_cpp_lcm2",                 (DL_FUNC) &_cheapr_cpp_lcm2,                 4},
    {"_cheapr_cpp_lcm2_vectorised",      (DL_FUNC) &_cheapr_cpp_lcm2_vectorised,      4},
    {"_cheapr_cpp_lead_sequence",        (DL_FUNC) &_cheapr_cpp_lead_sequence,        3},
    {"_cheapr_cpp_lengths",              (DL_FUNC) &_cheapr_cpp_lengths,              1},
    {"_cheapr_cpp_list_rm_null",         (DL_FUNC) &_cheapr_cpp_list_rm_null,         1},
    {"_cheapr_cpp_matrix_col_na_counts", (DL_FUNC) &_cheapr_cpp_matrix_col_na_counts, 1},
    {"_cheapr_cpp_matrix_missing_col",   (DL_FUNC) &_cheapr_cpp_matrix_missing_col,   3},
    {"_cheapr_cpp_matrix_missing_row",   (DL_FUNC) &_cheapr_cpp_matrix_missing_row,   3},
    {"_cheapr_cpp_matrix_row_na_counts", (DL_FUNC) &_cheapr_cpp_matrix_row_na_counts, 1},
    {"_cheapr_cpp_missing_col",          (DL_FUNC) &_cheapr_cpp_missing_col,          3},
    {"_cheapr_cpp_missing_row",          (DL_FUNC) &_cheapr_cpp_missing_row,          3},
    {"_cheapr_cpp_new_list",             (DL_FUNC) &_cheapr_cpp_new_list,             2},
    {"_cheapr_cpp_num_na",               (DL_FUNC) &_cheapr_cpp_num_na,               2},
    {"_cheapr_cpp_r_unnested_length",    (DL_FUNC) &_cheapr_cpp_r_unnested_length,    1},
    {"_cheapr_cpp_row_na_counts",        (DL_FUNC) &_cheapr_cpp_row_na_counts,        1},
    {"_cheapr_cpp_which_",               (DL_FUNC) &_cheapr_cpp_which_,               2},
    {"_cheapr_cpp_which_na",             (DL_FUNC) &_cheapr_cpp_which_na,             1},
    {"_cheapr_cpp_which_not_na",         (DL_FUNC) &_cheapr_cpp_which_not_na,         1},
    {"_cheapr_cpp_window_sequence",      (DL_FUNC) &_cheapr_cpp_window_sequence,      4},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_cheapr(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
