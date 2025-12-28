#include "cheapr.h"

// NA handling functions
// Author: Nick Christofides


#define CHEAPR_ANY_NA                                            \
for (R_xlen_t i = 0; i < n; ++i){                                \
  if (is_r_na(p_x[i])){                                          \
    out = true;                                                  \
    break;                                                       \
  }                                                              \
}                                                                \

#define CHEAPR_ALL_NA                                            \
for (R_xlen_t i = 0; i < n; ++i){                                \
  if (!is_r_na(p_x[i])){                                         \
    out = false;                                                 \
    break;                                                       \
  }                                                              \
}

#define CHEAPR_IS_NA                                                          \
if (n_threads > 1){                                                           \
  OMP_PARALLEL_FOR_SIMD(n_threads)                                            \
  for (R_xlen_t i = 0; i < n; ++i){                                           \
    set_value(p_out, i, static_cast<r_bool_t>(is_r_na(p_x[i])));                \
  }                                                                           \
} else {                                                                      \
  OMP_SIMD                                                                    \
  for (R_xlen_t i = 0; i < n; ++i){                                           \
    set_value(p_out, i, static_cast<r_bool_t>(is_r_na(p_x[i])));                \
  }                                                                           \
}


bool vec_any(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  const r_bool_t *p_x = logical_ptr_ro(x);

  bool out = false;
  for (R_xlen_t i = 0; i < n; ++i){
    if (is_r_true(p_x[i])){
     out = true;
      break;
    }
  }
  return out;
}

bool vec_all(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  const r_bool_t *p_x = logical_ptr_ro(x);

  bool out = true;

  for (R_xlen_t i = 0; i < n; ++i){
    if (p_x[i] == r_false){
      out = false;
      break;
    }
  }
  return out;
}

R_xlen_t na_count(SEXP x, bool recursive){
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t count = 0;
  int32_t NP = 0;
  int n_threads = calc_threads(n);
  visit_vector(x, [&](auto p_x) {

    using data_t = std::remove_const_t<std::remove_pointer_t<std::decay_t<decltype(p_x)>>>;

    auto default_scalar_count = [&] {
      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, x)); ++NP;
      SEXP scalar_true = SHIELD(as_vector(r_true)); ++NP;
      count = scalar_count(is_missing, scalar_true, true);
    };

    if constexpr (std::is_same_v<data_t, std::nullptr_t>){
      if (is_null(x)){
        return;
      } else {
        default_scalar_count();
      }
    } else if constexpr (std::is_same_v<data_t, SEXP>) {
      if (recursive){
        R_CheckStack(); // Check C Stack size isn't close to the limit
        for (R_xlen_t i = 0; i < n; ++i){
          count += na_count(p_x[i], true);
        }
      } else {
        default_scalar_count();
      }
    } else {
      if (n_threads > 1){
        _Pragma("omp parallel for simd num_threads(n_threads) reduction(+:count)")
        for (R_xlen_t i = 0; i < n; ++i) count += is_r_na(p_x[i]);
      } else {
        _Pragma("omp simd reduction(+:count)")
        for (R_xlen_t i = 0; i < n; ++i) count += is_r_na(p_x[i]);
      }
    }
  });
  YIELD(NP);
  return count;
}

[[cpp11::register]]
SEXP cpp_num_na(SEXP x, bool recursive){
  return as_vector(na_count(x, recursive));
}


[[cpp11::register]]
bool cpp_any_na(SEXP x, bool recursive){
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);
  bool out = false;
  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    return out;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_x = integer_ptr(x);
    CHEAPR_ANY_NA;
    break;
  }
  case CHEAPR_INT64SXP: {
    int64_t *p_x = integer64_ptr(x);
    CHEAPR_ANY_NA;
    break;
  }
  case REALSXP: {
    double *p_x = real_ptr(x);
    CHEAPR_ANY_NA;
    break;
  }
  case STRSXP: {
    const r_string_t *p_x = string_ptr_ro(x);
    CHEAPR_ANY_NA;
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    r_complex_t *p_x = complex_ptr(x);
    CHEAPR_ANY_NA;
    break;
  }
  case VECSXP: {
    if (recursive){
    R_CheckStack(); // Check C Stack size isn't close to the limit
    for (int i = 0; i < n; ++i){
      out = cpp_any_na(VECTOR_ELT(x, i), true);
      if (out) break;
    }
    break;
  }
  }
  default: {
    SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, x)); ++NP;
    out = vec_any(is_missing);
    break;
  }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
bool cpp_all_na(SEXP x, bool return_true_on_empty, bool recursive){
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);
  bool out = true;
  if (n == 0){
    if (return_true_on_empty){
      return true;
    } else {
      return false;
    }
  }
  switch ( CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    int *p_x = integer_ptr(x);
    CHEAPR_ALL_NA;
    break;
  }
  case CHEAPR_INT64SXP: {
    int64_t *p_x = integer64_ptr(x);
    CHEAPR_ALL_NA;
    break;
  }
  case REALSXP: {
    double *p_x = real_ptr(x);
    CHEAPR_ALL_NA;
    break;
  }
  case STRSXP: {
    const r_string_t *p_x = string_ptr_ro(x);
    CHEAPR_ALL_NA;
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    r_complex_t *p_x = complex_ptr(x);
    CHEAPR_ALL_NA;
    break;
  }
  case VECSXP: {
    if (recursive){
    R_CheckStack(); // Check C Stack size isn't close to the limit
    for (int i = 0; i < n; ++i){
      out = cpp_all_na(VECTOR_ELT(x, i), return_true_on_empty, true);
      if (!out) break;
    }
    break;
  }
  }
  default: {
    SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, x)); ++NP;
    out = vec_all(is_missing);
    break;
  }
  }
  YIELD(NP);
  return out;
}

// A multi-threaded version of `is.na()`
// lists are handled differently in that each element
// must contain only NA in all nested elements to be regarded as NA

[[cpp11::register]]
SEXP cpp_is_na(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out;
  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    out = SHIELD(new_vector<r_bool_t>(0));
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = SHIELD(new_vector<r_bool_t>(n));
    r_bool_t* RESTRICT p_out = logical_ptr(out);
    const int *p_x = integer_ptr(x);
    CHEAPR_IS_NA
    break;
  }
  case CHEAPR_INT64SXP: {
    out = SHIELD(new_vector<r_bool_t>(n));
    r_bool_t* RESTRICT p_out = logical_ptr(out);
    const int64_t *p_x = integer64_ptr_ro(x);
    CHEAPR_IS_NA
    break;
  }
  case REALSXP: {
    out = SHIELD(new_vector<r_bool_t>(n));
    r_bool_t* RESTRICT p_out = logical_ptr(out);
    const double *p_x = real_ptr(x);
    CHEAPR_IS_NA
    break;
  }
  case STRSXP: {
    out = SHIELD(new_vector<r_bool_t>(n));
    r_bool_t* RESTRICT p_out = logical_ptr(out);
    const r_string_t *p_x = string_ptr_ro(x);
    CHEAPR_IS_NA
    break;
  }
  case RAWSXP: {
    out = SHIELD(new_vector<r_bool_t>(n, r_false));
    break;
  }
  case CPLXSXP: {
    out = SHIELD(new_vector<r_bool_t>(n));
    r_bool_t* RESTRICT p_out = logical_ptr(out);
    const r_complex_t *p_x = complex_ptr(x);
    CHEAPR_IS_NA
    break;
  }
  case VECSXP: {
    if (!vec::is_object(x)){
    out = SHIELD(new_vector<r_bool_t>(n, r_false));
    break;
  }
  }
  default: {
    return eval_pkg_fun("is.na", "base", env::base_env, x);
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_df_row_na_counts(SEXP x){
  if (!is_df(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = list_ptr_ro(x);
  int32_t NP = 0;
  int num_col = Rf_length(x);
  int num_row = df::nrow(x);
  SEXP out = SHIELD(vec::new_vector<int>(num_row, 0)); ++NP;
  int* RESTRICT p_out = integer_ptr(out);
  for (int j = 0; j < num_col; ++j){
    switch ( CHEAPR_TYPEOF(p_x[j]) ){
    case LGLSXP:
    case INTSXP: {
      const int *p_xj = integer_ptr_ro(p_x[j]);
      OMP_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case CHEAPR_INT64SXP: {
      const int64_t *p_xj = integer64_ptr_ro(p_x[j]);
      OMP_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
     break;
    }
    case REALSXP: {
      const double *p_xj = real_ptr_ro(p_x[j]);
      OMP_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case STRSXP: {
      const r_string_t *p_xj = string_ptr_ro(p_x[j]);
      OMP_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case RAWSXP: {
      break;
    }
    case CPLXSXP: {
      const r_complex_t *p_xj = complex_ptr_ro(p_x[j]);
      OMP_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case VECSXP: {
      if (vec::is_object(p_x[j])){

      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, p_x[j])); ++NP;
      if (Rf_length(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_length(is_missing); ++NP;
        SEXP names = SHIELD(get_old_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(get_value<r_string_t>(names, j)), element_length, int_nrows);
      }
      const r_bool_t* RESTRICT p_is_missing = logical_ptr_ro(is_missing);
      for (int k = 0; k < num_row; ++k){
        p_out[k] += p_is_missing[k];
      }
    } else {
      const SEXP *p_xj = list_ptr_ro(p_x[j]);
      for (int i = 0; i < num_row; ++i){
        p_out[i] += cpp_all_na(p_xj[i], false, true);
      }
    }
    break;
    }
    default: {
      YIELD(NP);
      Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(p_x[j])));
    }
    }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_df_col_na_counts(SEXP x){
  if (!is_df(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = list_ptr_ro(x);
  int num_col = Rf_length(x);
  int32_t NP = 0;
  int num_row = df::nrow(x);
  SEXP out = SHIELD(vec::new_vector<int>(num_col, 0)); ++NP;
  int *p_out = integer_ptr(out);
  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (vec::is_object(p_x[j])){
      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, p_x[j])); ++NP;
      if (Rf_length(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_length(is_missing); ++NP;
        SEXP names = SHIELD(get_old_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(get_value<r_string_t>(names, j)), element_length, int_nrows);
      }
      int *p_is_missing = integer_ptr(is_missing);
      for (int k = 0; k < num_row; ++k){
        p_out[j] += p_is_missing[k];
      }
    } else {
      for (int i = 0; i < num_row; ++i){
        p_out[j] += cpp_all_na(VECTOR_ELT(p_x[j], i), false, true);
      }
    }
    break;
    }
    default: {
      p_out[j] = na_count(p_x[j], false);
      break;
    }
    }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_col_any_na(SEXP x, bool names){
  if (!is_df(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = list_ptr_ro(x);

  int32_t NP = 0;
  int num_row = df::nrow(x);
  int num_col = Rf_length(x);

  SEXP out = SHIELD(new_vector<r_bool_t>(num_col)); ++NP;
  int *p_out = integer_ptr(out);

  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (vec::is_object(p_x[j])){
      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, p_x[j])); ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = SHIELD(get_old_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(get_value<r_string_t>(names, j)), element_length, int_nrows);
      }
      p_out[j] = static_cast<int>(vec_any(is_missing));
    } else {
      bool any_na = false;
      bool all_na;
      for (int i = 0; i < num_row; ++i){
        // x[[i]] is only 'NA' if all nested elements are NA
        all_na = cpp_all_na(VECTOR_ELT(p_x[j], i), false, true);
        if (all_na){
          any_na = true;
          break;
        }
      }
      p_out[j] = any_na;
    }
    break;
    }
    default: {
      p_out[j] = cpp_any_na(p_x[j], false);
      break;
    }
    }
  }
  SEXP x_names = SHIELD(get_old_names(x)); ++NP;
  if (names){
    set_old_names(out, x_names);
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_col_all_na(SEXP x, bool names){
  if (!is_df(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = list_ptr_ro(x);

  int32_t NP = 0;
  int num_row = df::nrow(x);
  int num_col = Rf_length(x);

  SEXP out = SHIELD(new_vector<r_bool_t>(num_col)); ++NP;
  int *p_out = integer_ptr(out);

  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (vec::is_object(p_x[j])){
      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, p_x[j])); ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = SHIELD(get_old_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(get_value<r_string_t>(names, j)), element_length, int_nrows);
      }
      p_out[j] = static_cast<int>(vec_all(is_missing));
    } else {
      bool all_na2 = true;
      bool all_na;
      for (int i = 0; i < num_row; ++i){
        // x[[i]] is only 'NA' if all nested elements are NA
        all_na = cpp_all_na(VECTOR_ELT(p_x[j], i), false, true);
        if (!all_na){
          all_na2 = false;
          break;
        }
      }
      p_out[j] = all_na2;
    }
    break;
    }
    default: {
      p_out[j] = cpp_all_na(p_x[j], true, false);
      break;
    }
    }
  }
  SEXP x_names = SHIELD(get_old_names(x)); ++NP;
  if (names){
    set_old_names(out, x_names);
  }
  YIELD(NP);
  return out;
}

R_xlen_t cpp_clean_threshold(double threshold, bool threshold_is_prop, R_xlen_t n){
  if (is_r_na(threshold)){
    Rf_error("threshold cannot be NA");
  }
  R_xlen_t out = threshold;
  if (threshold_is_prop){
    if (threshold < 0){
      out = 0;
    } else if (threshold == r_limits::r_pos_inf){
      out = n + 1;
    } else {
      out = std::floor( (threshold * n) + 0.0000000001);
    }
  } else {
    if (threshold < 0){
      out = 0;
    }
    if (threshold == r_limits::r_pos_inf){
      out = n + 1;
    }
  }
  return out;
}

// Matrix methods
// The methods for matrices are substantially different

[[cpp11::register]]
SEXP cpp_matrix_row_na_counts(SEXP x){
  if (!Rf_isMatrix(x)){
    Rf_error("x must be a matrix");
  }
  R_xlen_t num_row = Rf_nrows(x);
  R_xlen_t num_col = Rf_ncols(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = SHIELD(vec::new_vector<int>(num_row, 0));
  int *p_out = integer_ptr(out);
  if (num_row > 0 && num_col > 0){
    switch ( CHEAPR_TYPEOF(x) ){
    case LGLSXP:
    case INTSXP: {
      int *p_x = integer_ptr(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    case REALSXP: {
      double *p_x = real_ptr(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    case CHEAPR_INT64SXP: {
      int64_t *p_x = integer64_ptr(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    case STRSXP: {
      const r_string_t *p_x = string_ptr_ro(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    case RAWSXP: {
      break;
    }
    case CPLXSXP: {
      r_complex_t *p_x = complex_ptr(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    default: {
      YIELD(1);
      Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
    }
    }
  }
  YIELD(1);
  return out;
}
[[cpp11::register]]
SEXP cpp_matrix_col_na_counts(SEXP x){
  if (!Rf_isMatrix(x)){
    Rf_error("x must be a matrix");
  }
  R_xlen_t num_row = Rf_nrows(x);
  R_xlen_t num_col = Rf_ncols(x);
  R_xlen_t n = Rf_xlength(x);
  bool new_col;
  SEXP out = SHIELD(vec::new_vector<int>(num_col, 0));
  int *p_out = integer_ptr(out);
  if (num_row > 0 && num_col > 0){
    switch ( CHEAPR_TYPEOF(x) ){
    case LGLSXP:
    case INTSXP: {
      int *p_x = integer_ptr(x);
      for (R_xlen_t i = 0, coli = 0, rowi = 0; i < n; ++rowi, ++i){
        new_col = rowi == num_row;
        if (new_col){
          ++coli;
          rowi = 0;
        }
        p_out[coli] += is_r_na(p_x[i]);
      }
      break;
    }
    case REALSXP: {
      double *p_x = real_ptr(x);
      for (R_xlen_t i = 0, coli = 0, rowi = 0; i < n; ++rowi, ++i){
        new_col = rowi == num_row;
        if (new_col){
          ++coli;
          rowi = 0;
        }
        p_out[coli] += is_r_na(p_x[i]);
      }
      break;
    }
    case CHEAPR_INT64SXP: {
      int64_t *p_x = integer64_ptr(x);
      for (R_xlen_t i = 0, coli = 0, rowi = 0; i < n; ++rowi, ++i){
        new_col = rowi == num_row;
        if (new_col){
          ++coli;
          rowi = 0;
        }
        p_out[coli] += is_r_na(p_x[i]);
      }
      break;
    }
    case STRSXP: {
      const r_string_t *p_x = string_ptr_ro(x);
      for (R_xlen_t i = 0, coli = 0, rowi = 0; i < n; ++rowi, ++i){
        new_col = rowi == num_row;
        if (new_col){
          ++coli;
          rowi = 0;
        }
        p_out[coli] += is_r_na(p_x[i]);
      }
      break;
    }
    case RAWSXP: {
      break;
    }
    case CPLXSXP: {
      r_complex_t *p_x = complex_ptr(x);
      for (R_xlen_t i = 0, coli = 0, rowi = 0; i < n; ++rowi, ++i){
        new_col = rowi == num_row;
        if (new_col){
          ++coli;
          rowi = 0;
        }
        p_out[coli] += is_r_na(p_x[i]);
      }
      break;
    }
    default: {
      YIELD(1);
      Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
    }
    }
  }
  YIELD(1);
  return out;
}

// Helpers to get matrix row and col names

SEXP matrix_rownames(SEXP x) {
  SEXP dimnames = SHIELD(get_attr(x, symbol::dim_names_sym));

  if (is_null(dimnames) ||
      TYPEOF(dimnames) != VECSXP ||
      Rf_length(dimnames) != 2){
    YIELD(1);
    return r_null;
  }
  YIELD(1);
  return VECTOR_ELT(dimnames, 0);
}

SEXP matrix_colnames(SEXP x) {
  SEXP dimnames = SHIELD(get_attr(x, symbol::dim_names_sym));

  if (is_null(dimnames) ||
      TYPEOF(dimnames) != VECSXP ||
      Rf_length(dimnames) != 2){
    YIELD(1);
    return r_null;
  }
  YIELD(1);
  return VECTOR_ELT(dimnames, 1);
}

[[cpp11::register]]
SEXP cpp_row_na_counts(SEXP x, bool names){
  bool is_matrix = Rf_isMatrix(x);
  bool is_data_frame = is_df(x);

  if (!is_matrix && !is_data_frame){
    Rf_error("x must be a matrix or data frame");
  }

  int32_t NP = 0;

  SEXP out;
  if (is_matrix){
    out = SHIELD(cpp_matrix_row_na_counts(x)); ++NP;
    if (names){
      SEXP row_names = SHIELD(vec::deep_copy(matrix_rownames(x))); ++NP;
      set_old_names(out, row_names);
    }
  } else {
    out = SHIELD(cpp_df_row_na_counts(x)); ++NP;
    if (names){
      SEXP row_names = SHIELD(vec::deep_copy(get_attr(x, symbol::row_names_sym))); ++NP;
      set_old_names(out, row_names);
    }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_col_na_counts(SEXP x, bool names){
  bool is_matrix = Rf_isMatrix(x);
  bool is_data_frame = is_df(x);

  if (!is_matrix && !is_data_frame){
    Rf_error("x must be a matrix or data frame");
  }

  int32_t NP = 0;

  SEXP out;
  if (is_matrix){
    out = SHIELD(cpp_matrix_col_na_counts(x)); ++NP;
    if (names){
      SEXP col_names = SHIELD(vec::deep_copy(matrix_colnames(x))); ++NP;
      set_old_names(out, col_names);
    }
  } else {
    out = SHIELD(cpp_df_col_na_counts(x)); ++NP;
    SEXP x_names = SHIELD(get_old_names(x)); ++NP;
    if (names){
      set_old_names(out, x_names);
    }
  }
  YIELD(NP);
  return out;
}
