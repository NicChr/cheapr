#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>

int int_div(int x, int y){
  return x / y;
}

int num_cores(){
  int out = Rf_asInteger(Rf_GetOption1(Rf_installChar(Rf_mkChar("cheapr.cores"))));
  if (out >= 1){
    return out;
  } else {
    return 1;
  }
}

[[cpp11::register]]
R_xlen_t cpp_vector_size(SEXP x){
  if (Rf_isFrame(x)){
    return Rf_xlength(Rf_getAttrib(x, R_RowNamesSymbol));
  } else if (Rf_isVectorList(x)){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return cpp_vector_size(VECTOR_ELT(x, 0));
    } else {
      // return Rf_xlength(x);
      int n = Rf_length(x);
      if (n == 0){
        return 0;
      } else {
        R_xlen_t init = cpp_vector_size(VECTOR_ELT(x, 0));
        for (int i = 1; i < n; ++i) {
          if (cpp_vector_size(VECTOR_ELT(x, i)) != init){
            Rf_error("All list elements must be of equal length");
          }
        }
        return init;
      }
    }
  } else {
    return Rf_xlength(x);
  }
}

[[cpp11::register]]
int cpp_vector_width(SEXP x){
  if (Rf_isFrame(x)){
    return Rf_length(Rf_getAttrib(x, R_NamesSymbol));
  } else if (Rf_isVectorList(x)){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return Rf_length(x);
    } else {
      int n = Rf_length(x);
      if (n == 0){
        return 0;
      } else {
        const SEXP *p_x = VECTOR_PTR_RO(x);
        R_xlen_t init = cpp_vector_size(p_x[0]);
        for (int i = 1; i < n; ++i) {
          if (cpp_vector_size(p_x[i]) != init){
            Rf_error("All list elements must be of equal length");
          }
        }
        return n;
      }
    }
  } else {
    return 0;
  }
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
