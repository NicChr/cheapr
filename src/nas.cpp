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

#define CHEAPR_IS_NA                                               \
if (n_cores > 1){                                                  \
  OMP_PARALLEL_FOR_SIMD                                            \
  for (R_xlen_t i = 0; i < n; ++i){                                \
    p_out[i] = is_r_na(p_x[i]);                                    \
  }                                                                \
} else {                                                           \
  OMP_FOR_SIMD                                                     \
  for (R_xlen_t i = 0; i < n; ++i){                                \
    p_out[i] = is_r_na(p_x[i]);                                    \
  }                                                                \
}

#define CHEAPR_NA_COUNT                                                   \
if (do_parallel){                                                         \
  _Pragma("omp parallel for simd num_threads(n_cores) reduction(+:count)")\
  for (R_xlen_t i = 0; i < n; ++i) count += is_r_na(p_x[i]);              \
} else {                                                                  \
  OMP_FOR_SIMD                                                            \
  for (R_xlen_t i = 0; i < n; ++i) count += is_r_na(p_x[i]);              \
}


bool vec_any(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  const r_bool_t *p_x = BOOLEAN_RO(x);

  bool out = false;
  for (R_xlen_t i = 0; i < n; ++i){
    if (p_x[i] == r_true){
     out = true;
      break;
    }
  }
  return out;
}

bool vec_all(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  const r_bool_t *p_x = BOOLEAN_RO(x);

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
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  bool do_parallel = n_cores > 1;
  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    const int *p_x = INTEGER_RO(x);
    CHEAPR_NA_COUNT
    break;
  }
  case CHEAPR_INT64SXP: {
    const int64_t *p_x = INTEGER64_PTR_RO(x);
    CHEAPR_NA_COUNT
    break;
  }
  case REALSXP: {
    const double *p_x = REAL_RO(x);
    CHEAPR_NA_COUNT
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_NA_COUNT
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    const Rcomplex *p_x = COMPLEX_RO(x);
    CHEAPR_NA_COUNT
    break;
  }
  case VECSXP: {
    // We use a recursive method if recursive is true
    // Otherwise we skip to the default section below
    if (recursive){
    const SEXP *p_x = LIST_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      count += na_count(p_x[i], true);
    }
    break;
  } else {
    const SEXP *p_x = LIST_PTR_RO(x);
    CHEAPR_NA_COUNT
    break;
  }
  }
  default: {
    SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), x)); ++NP;
    SEXP r_true = SHIELD(as_r_vec(true)); ++NP;
    count = scalar_count(is_missing, r_true, true);
    break;
  }
  }
  YIELD(NP);
  return count;
}

[[cpp11::register]]
SEXP cpp_num_na(SEXP x, bool recursive){
  return as_r_vec(na_count(x, recursive));
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
    int *p_x = INTEGER(x);
    CHEAPR_ANY_NA;
    break;
  }
  case CHEAPR_INT64SXP: {
    int64_t *p_x = INTEGER64_PTR(x);
    CHEAPR_ANY_NA;
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    CHEAPR_ANY_NA;
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_ANY_NA;
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    CHEAPR_ANY_NA;
    break;
  }
  case VECSXP: {
    if (recursive){
    for (int i = 0; i < n; ++i){
      out = cpp_any_na(VECTOR_ELT(x, i), true);
      if (out) break;
    }
    break;
  }
  }
  default: {
    SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), x)); ++NP;
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
    int *p_x = INTEGER(x);
    CHEAPR_ALL_NA;
    break;
  }
  case CHEAPR_INT64SXP: {
    int64_t *p_x = INTEGER64_PTR(x);
    CHEAPR_ALL_NA;
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    CHEAPR_ALL_NA;
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_ALL_NA;
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    CHEAPR_ALL_NA;
    break;
  }
  case VECSXP: {
    if (recursive){
    for (int i = 0; i < n; ++i){
      out = cpp_all_na(VECTOR_ELT(x, i), return_true_on_empty, true);
      if (!out) break;
    }
    break;
  }
  }
  default: {
    SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), x)); ++NP;
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
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  SEXP out;
  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    out = SHIELD(vec::new_vec(LGLSXP, 0));
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = SHIELD(vec::new_vec(LGLSXP, n));
    int* RESTRICT p_out = LOGICAL(out);
    const int *p_x = INTEGER(x);
    CHEAPR_IS_NA
    break;
  }
  case CHEAPR_INT64SXP: {
    out = SHIELD(vec::new_vec(LGLSXP, n));
    int* RESTRICT p_out = LOGICAL(out);
    const int64_t *p_x = INTEGER64_PTR_RO(x);
    CHEAPR_IS_NA
    break;
  }
  case REALSXP: {
    out = SHIELD(vec::new_vec(LGLSXP, n));
    int* RESTRICT p_out = LOGICAL(out);
    const double *p_x = REAL(x);
    CHEAPR_IS_NA
    break;
  }
  case STRSXP: {
    out = SHIELD(vec::new_vec(LGLSXP, n));
    int* RESTRICT p_out = LOGICAL(out);
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_IS_NA
    break;
  }
  case RAWSXP: {
    out = SHIELD(vec::new_vec(LGLSXP, n));
    int* RESTRICT p_out = LOGICAL(out);
    std::fill(p_out, p_out + n, 0);
    break;
  }
  case CPLXSXP: {
    out = SHIELD(vec::new_vec(LGLSXP, n));
    int* RESTRICT p_out = LOGICAL(out);
    const Rcomplex *p_x = COMPLEX(x);
    CHEAPR_IS_NA
    break;
  }
  case VECSXP: {
    if (!vec::is_object(x)){
    out = SHIELD(vec::new_vec(LGLSXP, n));
    int* RESTRICT p_out = LOGICAL(out);
    const SEXP *p_x = LIST_PTR_RO(x);
    CHEAPR_IS_NA
    break;
  }
  }
  default: {
    return eval_pkg_fun("is.na", "base", R_GetCurrentEnv(), x);
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
  const SEXP *p_x = LIST_PTR_RO(x);
  int32_t NP = 0;
  int num_col = Rf_length(x);
  int num_row = df::nrow(x);
  SEXP out = SHIELD(vec::new_vec(INTSXP, num_row)); ++NP;
  int* RESTRICT p_out = INTEGER(out);
  std::fill(p_out, p_out + num_row, 0);
  for (int j = 0; j < num_col; ++j){
    switch ( CHEAPR_TYPEOF(p_x[j]) ){
    case LGLSXP:
    case INTSXP: {
      const int *p_xj = INTEGER_RO(p_x[j]);
      OMP_FOR_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case CHEAPR_INT64SXP: {
      const int64_t *p_xj = INTEGER64_PTR_RO(p_x[j]);
      OMP_FOR_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
     break;
    }
    case REALSXP: {
      const double *p_xj = REAL_RO(p_x[j]);
      OMP_FOR_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case STRSXP: {
      const SEXP *p_xj = STRING_PTR_RO(p_x[j]);
      OMP_FOR_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case RAWSXP: {
      break;
    }
    case CPLXSXP: {
      const Rcomplex *p_xj = COMPLEX_RO(p_x[j]);
      OMP_FOR_SIMD
      for (int i = 0; i < num_row; ++i){
        p_out[i] += is_r_na(p_xj[i]);
      }
      break;
    }
    case VECSXP: {
      if (vec::is_object(p_x[j])){

      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), p_x[j])); ++NP;
      if (Rf_length(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_length(is_missing); ++NP;
        SEXP names = SHIELD(get_r_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(STRING_ELT(names, j)), element_length, int_nrows);
      }
      const int* RESTRICT p_is_missing = LOGICAL_RO(is_missing);
      for (int k = 0; k < num_row; ++k){
        p_out[k] += p_is_missing[k];
      }
    } else {
      const SEXP *p_xj = LIST_PTR_RO(p_x[j]);
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
  const SEXP *p_x = LIST_PTR_RO(x);
  int num_col = Rf_length(x);
  int32_t NP = 0;
  int num_row = df::nrow(x);
  SEXP out = SHIELD(vec::new_vec(INTSXP, num_col)); ++NP;
  int *p_out = INTEGER(out);
  std::fill(p_out, p_out + num_col, 0);
  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (vec::is_object(p_x[j])){
      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), p_x[j])); ++NP;
      if (Rf_length(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_length(is_missing); ++NP;
        SEXP names = SHIELD(get_r_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(STRING_ELT(names, j)), element_length, int_nrows);
      }
      int *p_is_missing = LOGICAL(is_missing);
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
  const SEXP *p_x = LIST_PTR_RO(x);

  int32_t NP = 0;
  int num_row = df::nrow(x);
  int num_col = Rf_length(x);

  SEXP out = SHIELD(vec::new_vec(LGLSXP, num_col)); ++NP;
  int *p_out = INTEGER(out);

  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (vec::is_object(p_x[j])){
      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), p_x[j])); ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = SHIELD(get_r_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(STRING_ELT(names, j)), element_length, int_nrows);
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
  SEXP x_names = SHIELD(get_r_names(x)); ++NP;
  if (names){
    set_r_names(out, x_names);
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_col_all_na(SEXP x, bool names){
  if (!is_df(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = LIST_PTR_RO(x);

  int32_t NP = 0;
  int num_row = df::nrow(x);
  int num_col = Rf_length(x);

  SEXP out = SHIELD(vec::new_vec(LGLSXP, num_col)); ++NP;
  int *p_out = INTEGER(out);

  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (vec::is_object(p_x[j])){
      SEXP is_missing = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), p_x[j])); ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = SHIELD(get_r_names(x));
        YIELD(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 utf8_char(STRING_ELT(names, j)), element_length, int_nrows);
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
  SEXP x_names = SHIELD(get_r_names(x)); ++NP;
  if (names){
    set_r_names(out, x_names);
  }
  YIELD(NP);
  return out;
}

R_xlen_t cpp_clean_threshold(double threshold, bool threshold_is_prop, R_xlen_t n){
  if (threshold != threshold){
    Rf_error("threshold cannot be NA");
  }
  R_xlen_t out = threshold;
  if (threshold_is_prop){
    if (threshold < 0){
      out = 0;
    } else if (threshold == R_PosInf){
      out = n + 1;
    } else {
      out = std::floor( (threshold * n) + 0.0000000001);
    }
  } else {
    if (threshold < 0){
      out = 0;
    }
    if (threshold == R_PosInf){
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
  SEXP out = SHIELD(vec::new_vec(INTSXP, num_row));
  int *p_out = INTEGER(out);
  std::fill(p_out, p_out + num_row, 0);
  if (num_row > 0 && num_col > 0){
    switch ( CHEAPR_TYPEOF(x) ){
    case LGLSXP:
    case INTSXP: {
      int *p_x = INTEGER(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    case REALSXP: {
      double *p_x = REAL(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    case CHEAPR_INT64SXP: {
      int64_t *p_x = INTEGER64_PTR(x);
      for (R_xlen_t i = 0, rowi = 0; i < n; ++rowi, ++i){
        if (rowi == num_row) rowi = 0;
        p_out[rowi] += is_r_na(p_x[i]);
      }
      break;
    }
    case STRSXP: {
      const SEXP *p_x = STRING_PTR_RO(x);
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
      Rcomplex *p_x = COMPLEX(x);
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
  SEXP out = SHIELD(vec::new_vec(INTSXP, num_col));
  int *p_out = INTEGER(out);
  std::fill(p_out, p_out + num_col, 0);
  if (num_row > 0 && num_col > 0){
    switch ( CHEAPR_TYPEOF(x) ){
    case LGLSXP:
    case INTSXP: {
      int *p_x = INTEGER(x);
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
      double *p_x = REAL(x);
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
      int64_t *p_x = INTEGER64_PTR(x);
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
      const SEXP *p_x = STRING_PTR_RO(x);
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
      Rcomplex *p_x = COMPLEX(x);
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
  SEXP dimnames = SHIELD(get_attrib(x, R_DimNamesSymbol));

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
  SEXP dimnames = SHIELD(get_attrib(x, R_DimNamesSymbol));

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
      set_r_names(out, row_names);
    }
  } else {
    out = SHIELD(cpp_df_row_na_counts(x)); ++NP;
    if (names){
      SEXP row_names = SHIELD(vec::deep_copy(get_attrib(x, R_RowNamesSymbol))); ++NP;
      set_r_names(out, row_names);
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
      set_r_names(out, col_names);
    }
  } else {
    out = SHIELD(cpp_df_col_na_counts(x)); ++NP;
    SEXP x_names = SHIELD(get_r_names(x)); ++NP;
    if (names){
      set_r_names(out, x_names);
    }
  }
  YIELD(NP);
  return out;
}
