#include "cheapr_cpp.h"

// A more memory-efficient which()
// Author: Nick Christofides

// Count the number of true values

R_xlen_t count_true(int *px, R_xlen_t n){
  R_xlen_t size = 0;
  if (n >= 100000){
#pragma omp parallel for simd num_threads(num_cores()) reduction(+:size)
    for (R_xlen_t j = 0; j < n; ++j) size += (px[j] == TRUE);
    return size;
  } else {
#pragma omp for simd
    for (R_xlen_t j = 0; j < n; ++j) size += (px[j] == TRUE);
    return size;
  }
}

[[cpp11::register]]
SEXP cpp_which_(SEXP x, bool invert){
  R_xlen_t n = Rf_xlength(x);
  int *p_x = LOGICAL(x);
  bool is_long = (n > integer_max_);
  if (invert){
    if (is_long){
      R_xlen_t size = count_true(p_x, n);
      R_xlen_t out_size = n - size;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] != TRUE);
      }
      Rf_unprotect(1);
      return out;
    } else {
      int size = count_true(p_x, n);
      int out_size = n - size;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] != TRUE);
      }
      Rf_unprotect(1);
      return out;
    }
  } else {
    if (is_long){
      R_xlen_t size = count_true(p_x, n);
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] == TRUE);
      }
      Rf_unprotect(1);
      return out;
    } else {
      int size = count_true(p_x, n);
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] == TRUE);
      }
      Rf_unprotect(1);
      return out;
    }
  }
}

[[cpp11::register]]
SEXP cpp_which_val(SEXP x, SEXP value, bool invert){
  int n_protections = 0;
  R_xlen_t n = Rf_xlength(x);
  bool is_long = (n > integer_max_);
  if (Rf_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  SEXP val_is_na = Rf_protect(cpp_is_na(value));
  ++n_protections;
  if (Rf_asLogical(val_is_na)){
    Rf_unprotect(n_protections);
    if (invert){
      return cpp_which_not_na(x);
    } else {
      return cpp_which_na(x);
    }
  }
#define WHICH_VAL(_val_)                                           \
  if (invert){                                                     \
    while (whichi < out_size){                                     \
      p_out[whichi] = i + 1;                                       \
      whichi += (p_x[i++] != _val_);                               \
    }                                                              \
  } else {                                                         \
    while (whichi < out_size){                                     \
      p_out[whichi] = i + 1;                                       \
      whichi += (p_x[i++] == _val_);                               \
    }                                                              \
  }
  R_xlen_t n_vals = scalar_count(x, value, false);
  R_xlen_t out_size = invert ? n - n_vals : n_vals;
  R_xlen_t whichi = 0;
  R_xlen_t i = 0;
  switch ( TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++n_protections;
    Rf_protect(value = Rf_coerceVector(value, INTSXP));
    ++n_protections;
    int val = Rf_asInteger(value);
    int *p_x = INTEGER(x);
    if (is_long){
      double *p_out = REAL(out);
      WHICH_VAL(val);
    } else {
      int *p_out = INTEGER(out);
      WHICH_VAL(val);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case REALSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++n_protections;
    Rf_protect(value = Rf_coerceVector(value, REALSXP));
    ++n_protections;
    double val = Rf_asReal(value);
    double *p_x = REAL(x);
    if (is_long){
      double *p_out = REAL(out);
      WHICH_VAL(val);
    } else {
      int *p_out = INTEGER(out);
      WHICH_VAL(val);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case STRSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++n_protections;
    Rf_protect(value = Rf_coerceVector(value, STRSXP));
    ++n_protections;
    SEXP val = Rf_protect(Rf_asChar(value));
    ++n_protections;
    SEXP *p_x = STRING_PTR(x);
    if (is_long){
      double *p_out = REAL(out);
      WHICH_VAL(val);
    } else {
      int *p_out = INTEGER(out);
      WHICH_VAL(val);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  default: {
    Rf_unprotect(n_protections);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }

}

// 2 more which() alternatives
// list cpp_which2(SEXP x){
//   int n = Rf_xlength(x);
//   int *p_x = LOGICAL(x);
//   // std::vector<int> out;
//   // out.reserve(n);
//   // for (int i = 0; i < n; ++i){
//   //   if (p_x[i] == TRUE){
//   //     out.push_back(i + 1);
//   //   }
//   // }
//   int k = 0;
//   std::vector<int> out(n);
//   for (int i = 0; i < n; ++i){
//     if (p_x[i] == TRUE){
//       out[k++] = i + 1;
//     } else {
//       out.pop_back();
//     }
//   }
//   return writable::list({
//     "out"_nm = out
//   });
// }
//
// SEXP cpp_which3(SEXP x){
//   int n = Rf_xlength(x);
//   int *p_x = LOGICAL(x);
//   int size = 0;
//   int j;
//   for (j = 0; j < n; ++j) size += (p_x[j] == TRUE);
//   SEXP out = Rf_protect(Rf_allocVector(INTSXP, size));
//   int *p_out = INTEGER(out);
//   int k = 0;
//   for (int i = 0; i < j; ++i){
//     if (p_x[i] == TRUE){
//       p_out[k++] = i + 1;
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }
