#include "cheapr_cpp.h"

#define CHEAPR_COUNT_NA(_ISNA_)                                \
for (R_xlen_t i = 0; i < n; ++i){                              \
  count += _ISNA_(p_x[i]);                                     \
}                                                              \

#define CHEAPR_ANY_NA(_ISNA_)                                  \
for (R_xlen_t i = 0; i < n; ++i){                              \
  if (_ISNA_(p_x[i])){                                         \
    out = true;                                                \
    break;                                                     \
  }                                                            \
}                                                              \

#define CHEAPR_ALL_NA(_ISNA_)                                  \
for (R_xlen_t i = 0; i < n; ++i){                              \
  if (!_ISNA_(p_x[i])){                                        \
    out = false;                                               \
    break;                                                     \
  }                                                            \
}                                                              \

#define CHEAPR_VEC_IS_NA(_ISNA_)                               \
for (R_xlen_t i = 0; i < n; ++i){                              \
  p_out[i] = _ISNA_(p_x[i]);                                   \
}                                                              \


R_xlen_t na_count(SEXP x, bool recursive){
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t count = 0;
  int NP = 0;
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  bool do_parallel = n_cores > 1;
  switch ( TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      CHEAPR_COUNT_NA(cheapr_is_na_int);
    } else {
      OMP_FOR_SIMD
      CHEAPR_COUNT_NA(cheapr_is_na_int);
    }
    break;
  }
  case REALSXP: {
    if (is_int64(x)){
    long long *p_x = INTEGER64_PTR(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      CHEAPR_COUNT_NA(cheapr_is_na_int64);
    } else {
      OMP_FOR_SIMD
      CHEAPR_COUNT_NA(cheapr_is_na_int64);
    }
  } else {
    double *p_x = REAL(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      CHEAPR_COUNT_NA(cheapr_is_na_dbl);
    } else {
      OMP_FOR_SIMD
      CHEAPR_COUNT_NA(cheapr_is_na_dbl);
    }
  }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      CHEAPR_COUNT_NA(cheapr_is_na_str);
    } else {
      OMP_FOR_SIMD
      CHEAPR_COUNT_NA(cheapr_is_na_str);
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      CHEAPR_COUNT_NA(cheapr_is_na_cplx);
    } else {
      OMP_FOR_SIMD
      CHEAPR_COUNT_NA(cheapr_is_na_cplx);
    }
    break;
  }
  case VECSXP: {
    // We use a recursive method if recursive is true
    // Otherwise we skip to the default section below
    if (recursive){
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      count += na_count(p_x[i], true);
    }
    break;
  }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x)); ++NP;
    SEXP r_true = Rf_protect(Rf_ScalarLogical(true)); ++NP;
    count = scalar_count(is_missing, r_true, true);
    break;
  }
  }
  Rf_unprotect(NP);
  return count;
}

[[cpp11::register]]
SEXP cpp_num_na(SEXP x, bool recursive){
  return xlen_to_r(na_count(x, recursive));
}


[[cpp11::register]]
bool cpp_any_na(SEXP x, bool recursive){
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);
  bool out = false;
  switch ( TYPEOF(x) ){
  case NILSXP: {
    return out;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    CHEAPR_ANY_NA(cheapr_is_na_int);
    break;
  }
  case REALSXP: {
    if (is_int64(x)){
    long long *p_x = INTEGER64_PTR(x);
    CHEAPR_ANY_NA(cheapr_is_na_int64);
  } else {
    double *p_x = REAL(x);
    CHEAPR_ANY_NA(cheapr_is_na_dbl);
  }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_ANY_NA(cheapr_is_na_str);
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    CHEAPR_ANY_NA(cheapr_is_na_cplx);
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
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x)); ++NP;
    SEXP any_missing = Rf_protect(cpp11::package("base")["any"](is_missing)); ++NP;
    out = Rf_asLogical(any_missing);
    break;
  }
  }
  Rf_unprotect(NP);
  return out;
}

[[cpp11::register]]
bool cpp_all_na(SEXP x, bool return_true_on_empty, bool recursive){
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);
  bool out = true;
  if (n == 0){
    if (return_true_on_empty){
      return true;
    } else {
      return false;
    }
  }
  switch ( TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    CHEAPR_ALL_NA(cheapr_is_na_int);
    break;
  }
  case REALSXP: {
    if (is_int64(x)){
    long long *p_x = INTEGER64_PTR(x);
    CHEAPR_ALL_NA(cheapr_is_na_int64);
  } else {
    double *p_x = REAL(x);
    CHEAPR_ALL_NA(cheapr_is_na_dbl);
  }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_ALL_NA(cheapr_is_na_str);
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    CHEAPR_ALL_NA(cheapr_is_na_cplx);
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
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x)); ++NP;
    SEXP all_missing = Rf_protect(cpp11::package("base")["all"](is_missing)); ++NP;
    out = Rf_asLogical(all_missing);
    break;
  }
  }
  Rf_unprotect(NP);
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
  switch ( TYPEOF(x) ){
  case NILSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, 0));
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    int *p_x = INTEGER(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_VEC_IS_NA(cheapr_is_na_int);
    } else {
      OMP_FOR_SIMD
      CHEAPR_VEC_IS_NA(cheapr_is_na_int);
    }

    break;
  }
  case REALSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    if (is_int64(x)){
      long long *p_x = INTEGER64_PTR(x);
      if (n_cores > 1){
        OMP_PARALLEL_FOR_SIMD
        CHEAPR_VEC_IS_NA(cheapr_is_na_int64);
      } else {
        OMP_FOR_SIMD
        CHEAPR_VEC_IS_NA(cheapr_is_na_int64);
      }
    } else {
      double *p_x = REAL(x);
      if (n_cores > 1){
        OMP_PARALLEL_FOR_SIMD
        CHEAPR_VEC_IS_NA(cheapr_is_na_dbl);
      } else {
        OMP_FOR_SIMD
        CHEAPR_VEC_IS_NA(cheapr_is_na_dbl);
      }
    }

    break;
  }
  case STRSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    const SEXP *p_x = STRING_PTR_RO(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_VEC_IS_NA(cheapr_is_na_str);
    } else {
      OMP_FOR_SIMD
      CHEAPR_VEC_IS_NA(cheapr_is_na_str);
    }
    break;
  }
  case RAWSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    memset(p_out, 0, n * sizeof(int));
    break;
  }
  case CPLXSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    Rcomplex *p_x = COMPLEX(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_VEC_IS_NA(cheapr_is_na_cplx);
    } else {
      OMP_FOR_SIMD
      CHEAPR_VEC_IS_NA(cheapr_is_na_cplx);
    }
    break;
  }
  case VECSXP: {
    if (!Rf_isObject(x)){
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = cpp_all_na(p_x[i], false, true);
    }
    break;
  }
  }
  default: {
    out = Rf_protect(cpp11::package("base")["is.na"](x));
    break;
  }
  }
  Rf_unprotect(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_df_row_na_counts(SEXP x){
  if (!Rf_isFrame(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  int num_col = Rf_length(x);
  int NP = 0;
  R_xlen_t num_row = cpp_df_nrow(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, num_row)); ++NP;
  int *p_out = INTEGER(out);
  memset(p_out, 0, num_row * sizeof(int));
  int n_cores = num_row >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  bool do_parallel = n_cores > 1;
  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case LGLSXP:
    case INTSXP: {
      int *p_xj = INTEGER(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_out[i] += cheapr_is_na_int(p_xj[i]);
      }
      break;
    }
    case REALSXP: {
      if (is_int64(p_x[j])){
      long long *p_xj = (long long *) REAL(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_out[i] += cheapr_is_na_int64(p_xj[i]);
      }
    } else {
      double *p_xj = REAL(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_out[i] += cheapr_is_na_dbl(p_xj[i]);
      }
    }
      break;
    }
    case STRSXP: {
      const SEXP *p_xj = STRING_PTR_RO(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_out[i] += cheapr_is_na_str(p_xj[i]);
      }
      break;
    }
    case RAWSXP: {
      break;
    }
    case CPLXSXP: {
      Rcomplex *p_xj = COMPLEX(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_out[i] += cheapr_is_na_cplx(p_xj[i]);
      }
      break;
    }
    case VECSXP: {
      if (Rf_isObject(p_x[j])){
      SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](p_x[j])); ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
        Rf_unprotect(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 CHAR(STRING_ELT(names, j)), element_length, int_nrows);
      }
      int *p_is_missing = LOGICAL(is_missing);
      for (R_xlen_t k = 0; k < num_row; ++k){
        p_out[k] += p_is_missing[k];
      }
    } else {
      const SEXP *p_xj = VECTOR_PTR_RO(p_x[j]);
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_out[i] += cpp_all_na(p_xj[i], false, true);
      }
    }
    break;
    }
    default: {
      Rf_unprotect(NP);
      Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(p_x[j])));
    }
    }
  }
  Rf_unprotect(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_df_col_na_counts(SEXP x){
  if (!Rf_isFrame(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  int num_col = Rf_length(x);
  int NP = 0;
  R_xlen_t num_row = cpp_df_nrow(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, num_col)); ++NP;
  int *p_out = INTEGER(out);
  memset(p_out, 0, num_col * sizeof(int));
  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (Rf_isObject(p_x[j])){
      SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](p_x[j]));
      ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
        Rf_unprotect(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 CHAR(STRING_ELT(names, j)), element_length, int_nrows);
      }
      int *p_is_missing = LOGICAL(is_missing);
      for (R_xlen_t k = 0; k < num_row; ++k){
        p_out[j] += p_is_missing[k];
      }
    } else {
      for (R_xlen_t i = 0; i < num_row; ++i){
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
  Rf_unprotect(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_col_any_na(SEXP x, bool names){
  if (!Rf_isFrame(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);

  int NP = 0;
  int num_row = cpp_df_nrow(x);
  int num_col = Rf_length(x);

  SEXP out = Rf_protect(Rf_allocVector(LGLSXP, num_col)); ++NP;
  int *p_out = INTEGER(out);

  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (Rf_isObject(p_x[j])){
      SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](p_x[j]));
      cpp11::function r_any = cpp11::package("base")["any"];
      ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
        Rf_unprotect(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 CHAR(STRING_ELT(names, j)), element_length, int_nrows);
      }
      SEXP r_any_true = Rf_protect(r_any(is_missing)); ++NP;
      p_out[j] = Rf_asLogical(r_any_true);
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
  if (names) cpp_copy_names(x, out, true);
  Rf_unprotect(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_col_all_na(SEXP x, bool names){
  if (!Rf_isFrame(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);

  int NP = 0;
  int num_row = cpp_df_nrow(x);
  int num_col = Rf_length(x);

  SEXP out = Rf_protect(Rf_allocVector(LGLSXP, num_col)); ++NP;
  int *p_out = INTEGER(out);

  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (Rf_isObject(p_x[j])){
      SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](p_x[j]));
      cpp11::function r_all = cpp11::package("base")["all"];
      ++NP;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing); ++NP;
        SEXP names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
        Rf_unprotect(NP);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 CHAR(STRING_ELT(names, j)), element_length, int_nrows);
      }
      SEXP r_all_true = Rf_protect(r_all(is_missing)); ++NP;
      p_out[j] = Rf_asLogical(r_all_true);
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
  if (names) cpp_copy_names(x, out, true);
  Rf_unprotect(NP);
  return out;
}

// INCOMPLETE
//TO-DO - Rewrite this but for all_na()
// More likely to be faster for the all_na case
// SEXP cpp_row_any_na(SEXP x, bool names){
//   if (!Rf_isFrame(x)){
//     Rf_error("x must be a data frame");
//   }
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//
//   int NP = 0;
//   int num_row = cpp_df_nrow(x);
//   int num_col = Rf_length(x);
//
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, num_row)); ++NP;
//   int *p_out = INTEGER(out);
//
//   bool row_has_na = false;
//
//   for (int i = 0; i < num_row; ++i){
//     row_has_na = false;
//     for (int j = 0; j < num_col; ++j){
//       switch ( TYPEOF(p_x[j]) ){
//       case LGLSXP:
//       case INTSXP: {
//         if (cheapr_is_na_int(INTEGER(p_x[j])[i])){
//         row_has_na = true;
//         break;
//       }
//         break;
//       }
//       case REALSXP: {
//         if (is_int64(p_x[j])){
//         long long *p_xj = (long long *) REAL(p_x[j]);
//         if (cheapr_is_na_int64(p_xj[i])){
//           row_has_na = true;
//           break;
//         }
//       } else {
//         if (cheapr_is_na_dbl(REAL(p_x[j])[i])){
//           row_has_na = true;
//           break;
//         }
//       }
//         break;
//       }
//       case STRSXP: {
//         if (cheapr_is_na_str(STRING_PTR_RO(p_x[j])[i])){
//         row_has_na = true;
//         break;
//       }
//         break;
//       }
//       default: {
//         Rf_unprotect(NP);
//         Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(p_x[j])));
//       }
//       }
//     }
//     p_out[i] = row_has_na;
//   }
//   if (names) cpp_copy_names(x, out, true);
//   Rf_unprotect(NP);
//   return out;
// }

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
  int num_row = Rf_nrows(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, num_row));
  int *p_out = INTEGER(out);
  memset(p_out, 0, num_row * sizeof(int));
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  bool do_parallel = n_cores > 1;
#pragma omp parallel num_threads(n_cores) if(do_parallel)
  switch ( TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += cheapr_is_na_int(p_x[i]);
    }
    break;
  }
  case REALSXP: {
    if (is_int64(x)){
    long long *p_x = INTEGER64_PTR(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += cheapr_is_na_int64(p_x[i]);
    }
  } else {
    double *p_x = REAL(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += cheapr_is_na_dbl(p_x[i]);
    }
  }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += cheapr_is_na_str(p_x[i]);
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += cheapr_is_na_cplx(p_x[i]);
    }
    break;
  }
  default: {
    Rf_unprotect(1);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(1);
  return out;
}
[[cpp11::register]]
SEXP cpp_matrix_col_na_counts(SEXP x){
  if (!Rf_isMatrix(x)){
    Rf_error("x must be a matrix");
  }
  int num_row = Rf_nrows(x);
  int num_col = Rf_ncols(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, num_col));
  int *p_out = INTEGER(out);
  memset(p_out, 0, num_col * sizeof(int));
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  bool do_parallel = n_cores > 1;
#pragma omp parallel num_threads(n_cores) if(do_parallel)
  switch ( TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += cheapr_is_na_int(p_x[i]);
    }
    break;
  }
  case REALSXP: {
    if (is_int64(x)){
    long long *p_x = INTEGER64_PTR(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += cheapr_is_na_int64(p_x[i]);
    }
  } else {
    double *p_x = REAL(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += cheapr_is_na_dbl(p_x[i]);
    }
  }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += cheapr_is_na_str(p_x[i]);
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += cheapr_is_na_cplx(p_x[i]);
    }
    break;
  }
  default: {
    Rf_unprotect(1);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(1);
  return out;
}

// Helpers to get matrix row and col names

SEXP matrix_rownames(SEXP x) {
  SEXP dimnames = Rf_protect(Rf_getAttrib(x, R_DimNamesSymbol));

  if (Rf_isNull(dimnames) ||
      TYPEOF(dimnames) != VECSXP ||
      Rf_length(dimnames) != 2){
    Rf_unprotect(1);
    return R_NilValue;
  }
  Rf_unprotect(1);
  return VECTOR_ELT(dimnames, 0);
}
SEXP matrix_colnames(SEXP x) {
  SEXP dimnames = Rf_protect(Rf_getAttrib(x, R_DimNamesSymbol));

  if (Rf_isNull(dimnames) ||
      TYPEOF(dimnames) != VECSXP ||
      Rf_length(dimnames) != 2){
    Rf_unprotect(1);
    return R_NilValue;
  }
  Rf_unprotect(1);
  return VECTOR_ELT(dimnames, 1);
}

[[cpp11::register]]
SEXP cpp_row_na_counts(SEXP x, bool names){
  bool is_matrix = Rf_isMatrix(x);
  bool is_data_frame = Rf_isFrame(x);

  if (!is_matrix && !is_data_frame){
    Rf_error("x must be a matrix or data frame");
  }

  int NP = 0;

  SEXP out;
  if (is_matrix){
    out = Rf_protect(cpp_matrix_row_na_counts(x)); ++NP;
    if (names){
      SEXP row_names = Rf_protect(Rf_duplicate(matrix_rownames(x))); ++NP;
      Rf_setAttrib(out, R_NamesSymbol, row_names);
    }
  } else {
    out = Rf_protect(cpp_df_row_na_counts(x)); ++NP;
    if (names){
      SEXP row_names = Rf_protect(Rf_duplicate(Rf_getAttrib(x, R_RowNamesSymbol))); ++NP;
      Rf_setAttrib(out, R_NamesSymbol, row_names);
    }
  }
  Rf_unprotect(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_col_na_counts(SEXP x, bool names){
  bool is_matrix = Rf_isMatrix(x);
  bool is_data_frame = Rf_isFrame(x);

  if (!is_matrix && !is_data_frame){
    Rf_error("x must be a matrix or data frame");
  }

  int NP = 0;

  SEXP out;
  if (is_matrix){
    out = Rf_protect(cpp_matrix_col_na_counts(x)); ++NP;
    if (names){
      SEXP col_names = Rf_protect(Rf_duplicate(matrix_colnames(x))); ++NP;
      Rf_setAttrib(out, R_NamesSymbol, col_names);
    }
  } else {
    out = Rf_protect(cpp_df_col_na_counts(x)); ++NP;
    if (names) cpp_copy_names(x, out, true);
  }
  Rf_unprotect(NP);
  return out;
}
