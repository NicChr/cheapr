#include "cheapr_cpp.h"

int int_div(int x, int y){
  return x / y;
}

[[cpp11::register]]
R_xlen_t cpp_vec_length(SEXP x){
  if (Rf_isFrame(x)){
    return cpp_df_nrow(x);
    // Is x a list?
  } else if (Rf_isVectorList(x)){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return cpp_vec_length(VECTOR_ELT(x, 0));
    } else if (Rf_inherits(x, "POSIXlt")){
      return Rf_xlength(VECTOR_ELT(x, 0));
    } else if (Rf_isObject(x)){
      return Rf_asReal(cpp11::package("base")["length"](x));
    } else {
      return Rf_xlength(x);
    }
    // Catch-all
  } else {
    return Rf_xlength(x);
  }
}

int num_cores(){
  SEXP num_cores = Rf_protect(Rf_GetOption1(Rf_installChar(Rf_mkChar("cheapr.cores"))));
  int out = Rf_asInteger(num_cores);
  Rf_unprotect(1);
  return out >= 1 ? out : 1;
}

SEXP xlen_to_r(R_xlen_t x){
  return x > integer_max_ ? Rf_ScalarReal(x) : Rf_ScalarInteger(x);
}

R_xlen_t cpp_df_nrow(SEXP x){
  return Rf_xlength(Rf_getAttrib(x, R_RowNamesSymbol));
}

R_xlen_t cpp_unnested_length(SEXP x){
  if (!Rf_isVectorList(x)){
    return Rf_xlength(x);
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t out = 0;
  for (R_xlen_t i = 0; i < n; ++i){
    out += Rf_isVectorList(p_x[i]) ? cpp_unnested_length(p_x[i]) : Rf_xlength(p_x[i]);
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_r_unnested_length(SEXP x){
  return xlen_to_r(cpp_unnested_length(x));
}

[[cpp11::register]]
SEXP cpp_lengths(SEXP x) {
  R_xlen_t n = Rf_xlength(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
  int *p_out = INTEGER(out);
  if (!Rf_isVectorList(x)){
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = 1;
    }
    Rf_unprotect(1);
    return out;
  } else {
    const SEXP* p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = cpp_vec_length(p_x[i]);
    }
    Rf_unprotect(1);
    return out;
  }
}

[[cpp11::register]]
SEXP cpp_new_list(R_xlen_t size, SEXP default_value) {
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, size));
  if (!Rf_isNull(default_value)){
    for (R_xlen_t i = 0; i < size; ++i) {
      SET_VECTOR_ELT(out, i, default_value);
    }
  }
  Rf_unprotect(1);
  return out;
}

// Remove NULL elements from list

[[cpp11::register]]
SEXP cpp_drop_null(SEXP l, bool always_shallow_copy) {
  Rf_protect(l = Rf_coerceVector(l, VECSXP));
  const SEXP *p_l = VECTOR_PTR_RO(l);
  int n = Rf_length(l);
  int n_null = 0;
  for (int i = 0; i < n; ++i) {
    n_null += (p_l[i] == R_NilValue);
  }
  if (n_null == 0 && !always_shallow_copy){
    Rf_unprotect(1);
    return l;
  }
  int n_keep = n - n_null;
  int whichj = 0;
  int j = 0;

  // Which list elements should we keep?

  SEXP keep = Rf_protect(Rf_allocVector(INTSXP, n_keep));
  int *p_keep = INTEGER(keep);
  while (whichj < n_keep){
    p_keep[whichj] = j;
    whichj += (p_l[j++] != R_NilValue);
  }

  // Subset on both the list and names of the list

  SEXP out = Rf_protect(Rf_allocVector(VECSXP, n_keep));
  SEXP names = Rf_protect(Rf_getAttrib(l, R_NamesSymbol));
  bool has_names = !Rf_isNull(names);
  if (has_names){
    SEXP *p_names = STRING_PTR(names);
    SEXP out_names = Rf_protect(Rf_allocVector(STRSXP, n_keep));
    for (int k = 0; k < n_keep; ++k) {
      SET_STRING_ELT(out_names, k, p_names[p_keep[k]]);
      SET_VECTOR_ELT(out, k, p_l[p_keep[k]]);
    }
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    Rf_unprotect(5);
    return out;
  } else {
    for (int k = 0; k < n_keep; ++k) {
      SET_VECTOR_ELT(out, k, p_l[p_keep[k]]);
    }
    Rf_unprotect(4);
    return out;
  }
}

[[cpp11::register]]
SEXP cpp_list_as_df(SEXP x) {
  SEXP out = Rf_protect(cpp_drop_null(x, true));
  int N; // Number of rows
  if (Rf_inherits(x, "data.frame")){
    N = cpp_df_nrow(x);
  } else if (Rf_length(out) == 0){
    N = 0;
  } else {
    N = cpp_vec_length(VECTOR_ELT(out, 0));
  }
  SEXP df_str = Rf_protect(Rf_ScalarString(Rf_mkChar("data.frame")));
  if (N > 0){
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 2));
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -N;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    Rf_classgets(out, df_str);
    Rf_unprotect(3);
    return out;
  } else {
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    Rf_classgets(out, df_str);
    Rf_unprotect(3);
    return out;
  }
}

SEXP r_address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return Rf_mkChar(buf);
}

[[cpp11::register]]
SEXP r_copy(SEXP x){
  return Rf_duplicate(x);
}

// Potentially useful for rolling calculations
// Computes the rolling number of true values in a given
// series of consecutive true values

// SEXP cpp_run_id(SEXP x, bool left_to_right){
//   if (!Rf_isLogical(x)){
//     Rf_error("x must be a logical vector");
//   }
//   R_xlen_t count = 0;
//   R_xlen_t n = Rf_xlength(x);
//   SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
//   int *p_out = INTEGER(out);
//   int *p_x = LOGICAL(x);
//   if (left_to_right){
//     for (R_xlen_t i = 0; i < n; ++i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   } else {
//     for (R_xlen_t i = n - 1; i >= 0; --i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }

// Would use data.table as it is very efficient, but would require extra dependency
// Here x must be an integer vector

// SEXP cpp_between(SEXP x, SEXP lower, SEXP upper){
//   int *p_x = INTEGER(x);
//   int n = Rf_length(x);
//   if (Rf_length(lower) != 1){
//     Rf_error("lower must be of length 1");
//   }
//   if (Rf_length(upper) != 1){
//     Rf_error("upper must be of length 1");
//   }
//   Rf_protect(lower = Rf_coerceVector(lower, INTSXP));
//   Rf_protect(upper = Rf_coerceVector(upper, INTSXP));
//   int lo = Rf_asInteger(lower);
//   int hi = Rf_asInteger(upper);
//   if (hi < lo){
//     int hi2 = hi;
//     hi = lo;
//     lo = hi2;
//   }
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
//   int *p_out = LOGICAL(out);
//   bool do_parallel = n >= 100000;
//   int n_cores = do_parallel ? num_cores() : 1;
//   if (lo == NA_INTEGER || hi == NA_INTEGER){
// #pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
//     for (int i = 0; i < n; ++i){
//       p_out[i] = NA_LOGICAL;
//     }
//   } else {
//     unsigned int rng = hi - lo;
// #pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
//     for (int i = 0; i < n; ++i){
//       int xi = p_x[i];
//       p_out[i] = (xi != NA_INTEGER) ? unsigned(xi - lo) <= rng : NA_LOGICAL;
//     }
//   }
//   Rf_unprotect(3);
//   return out;
// }
