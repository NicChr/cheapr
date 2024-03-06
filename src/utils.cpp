#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>

int int_div(int x, int y){
  return x / y;
}

R_xlen_t cpp_vec_length(SEXP x){
  return Rf_inherits(x, "vctrs_rcrd") ? cpp_vec_length(VECTOR_ELT(x, 0)) : Rf_xlength(x);
}

int num_cores(){
  int out = Rf_asInteger(Rf_GetOption1(Rf_installChar(Rf_mkChar("cheapr.cores"))));
  if (out >= 1){
    return out;
  } else {
    return 1;
  }
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

// bool cpp_is_list_df_like(SEXP x){
//   if (!Rf_isVectorList(x)){
//     Rf_error("x must be a list.");
//   }
//   if (Rf_isFrame(x)) return true;
//   R_xlen_t n = Rf_xlength(x);
//   bool out = true;
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//   R_xlen_t init = cpp_vec_length(p_x[0]);
//   for (R_xlen_t i = 1; i < n; ++i) {
//     if (cpp_vec_length(p_x[i]) != init){
//      out = false;
//       break;
//     }
//   }
//   return out;
// }

//
// bool list_has_list(SEXP x){
//   Rf_protect(x = Rf_coerceVector(x, VECSXP));
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//   R_xlen_t n = Rf_xlength(x);
//   bool out = false;
//   for (R_xlen_t i = 0; i < n; ++i){
//     if (Rf_isVectorList(p_x[i])){
//       out = true;
//       break;
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }

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
