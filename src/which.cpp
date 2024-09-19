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

#define CHEAPR_WHICH_VAL(_val_)                                \
while (whichi < out_size){                                     \
  p_out[whichi] = i + 1;                                       \
  whichi += (p_x[i++] == _val_);                               \
}                                                              \

#define CHEAPR_WHICH_VAL_INVERTED(_val_)                       \
while (whichi < out_size){                                     \
  p_out[whichi] = i + 1;                                       \
  whichi += (p_x[i++] != _val_);                               \
}                                                              \

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
      CHEAPR_WHICH_VAL_INVERTED(TRUE);
      Rf_unprotect(1);
      return out;
    } else {
      int size = count_true(p_x, n);
      int out_size = n - size;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL_INVERTED(TRUE);
      Rf_unprotect(1);
      return out;
    }
  } else {
    if (is_long){
      R_xlen_t out_size = count_true(p_x, n);
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL(TRUE);
      Rf_unprotect(1);
      return out;
    } else {
      int out_size = count_true(p_x, n);
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL(TRUE);
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
  if (Rf_isVectorList(x)){
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
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
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int *p_out = INTEGER(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
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
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int *p_out = INTEGER(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
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
    const SEXP *p_x = STRING_PTR_RO(x);
    if (is_long){
      double *p_out = REAL(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int *p_out = INTEGER(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
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

// Memory-efficient which(is.na(x))

[[cpp11::register]]
SEXP cpp_which_na(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  bool is_short = (n <= integer_max_);
  switch ( TYPEOF(x) ){
  case NILSXP: {
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_unprotect(1);
    return out;
  }
  case LGLSXP:
  case INTSXP: {
    R_xlen_t count = na_count(x, true);
    int *p_x = INTEGER(x);
    if (is_short){
      int out_size = count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL(NA_INTEGER);
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL(NA_INTEGER);
      Rf_unprotect(1);
      return out;
    }
  }
  case REALSXP: {
    R_xlen_t count = na_count(x, true);
    if (Rf_inherits(x, "integer64")){
      long long *p_x = INTEGER64_PTR(x);
      if (is_short){
        int out_size = count;
        SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
        int *p_out = INTEGER(out);
        int whichi = 0;
        int i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i++] == NA_INTEGER64);
        }
        Rf_unprotect(1);
        return out;
      } else {
        R_xlen_t out_size = count;
        SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
        double *p_out = REAL(out);
        R_xlen_t whichi = 0;
        R_xlen_t i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i++] == NA_INTEGER64);
        }
        Rf_unprotect(1);
        return out;
      }
    } else {
      double *p_x = REAL(x);
      if (is_short){
        int out_size = count;
        SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
        int *p_out = INTEGER(out);
        int whichi = 0;
        int i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i] != p_x[i]);
          ++i;
        }
        Rf_unprotect(1);
        return out;
      } else {
        R_xlen_t out_size = count;
        SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
        double *p_out = REAL(out);
        R_xlen_t whichi = 0;
        R_xlen_t i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i] != p_x[i]);
          ++i;
        }
        Rf_unprotect(1);
        return out;
      }
    }
  }
  case STRSXP: {
    R_xlen_t count = na_count(x, true);
    const SEXP *p_x = STRING_PTR_RO(x);
    if (is_short){
      int out_size = count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL(NA_STRING);
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL(NA_STRING);
      Rf_unprotect(1);
      return out;
    }
  }
  case RAWSXP: {
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_unprotect(1);
    return out;
  }
  case CPLXSXP: {
    R_xlen_t count = na_count(x, true);
    Rcomplex *p_x = COMPLEX(x);
    if (is_short){
      int out_size = count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += cheapr_is_na_cplx(p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += cheapr_is_na_cplx(p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x));
    SEXP out = Rf_protect(cpp_which_(is_missing, false));
    Rf_unprotect(2);
    return out;
  }
  }
}

[[cpp11::register]]
SEXP cpp_which_not_na(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  bool is_short = (n <= integer_max_);
  switch ( TYPEOF(x) ){
  case NILSXP: {
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_unprotect(1);
    return out;
  }
  case LGLSXP:
  case INTSXP: {
    R_xlen_t count = na_count(x, true);
    int *p_x = INTEGER(x);
    if (is_short){
      int out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL_INVERTED(NA_INTEGER);
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL_INVERTED(NA_INTEGER);
      Rf_unprotect(1);
      return out;
    }
  }
  case REALSXP: {
    R_xlen_t count = na_count(x, true);
    if (Rf_inherits(x, "integer64")){
      long long *p_x = INTEGER64_PTR(x);
      if (is_short){
        int out_size = n - count;
        SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
        int *p_out = INTEGER(out);
        int whichi = 0;
        int i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i++] != NA_INTEGER64);
        }
        Rf_unprotect(1);
        return out;
      } else {
        R_xlen_t out_size = n - count;
        SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
        double *p_out = REAL(out);
        R_xlen_t whichi = 0;
        R_xlen_t i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i++] != NA_INTEGER64);
        }
        Rf_unprotect(1);
        return out;
      }
    } else {
      double *p_x = REAL(x);
      if (is_short){
        int out_size = n - count;
        SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
        int *p_out = INTEGER(out);
        int whichi = 0;
        int i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i] == p_x[i]);
          ++i;
        }
        Rf_unprotect(1);
        return out;
      } else {
        R_xlen_t out_size = n - count;
        SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
        double *p_out = REAL(out);
        R_xlen_t whichi = 0;
        R_xlen_t i = 0;
        while (whichi < out_size){
          p_out[whichi] = i + 1;
          whichi += (p_x[i] == p_x[i]);
          ++i;
        }
        Rf_unprotect(1);
        return out;
      }
    }
  }
  case STRSXP: {
    R_xlen_t count = na_count(x, true);
    const SEXP *p_x = STRING_PTR_RO(x);
    if (is_short){
      int out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL_INVERTED(NA_STRING);
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL_INVERTED(NA_STRING);
      Rf_unprotect(1);
      return out;
    }
  }
  case RAWSXP: {
    if (is_short){
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
    int *p_out = INTEGER(out);
    for (int i = 0; i < n; ++i){
      p_out[i] = i + 1;
    }
    Rf_unprotect(1);
    return out;
  } else {
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
    double *p_out = REAL(out);
    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = i + 1;
    }
    Rf_unprotect(1);
    return out;
  }
  }
  case CPLXSXP: {
    R_xlen_t count = na_count(x, true);
    Rcomplex *p_x = COMPLEX(x);
    if (is_short){
      int out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += !cheapr_is_na_cplx(p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += !cheapr_is_na_cplx(p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x));
    SEXP out = Rf_protect(cpp_which_(is_missing, true));
    Rf_unprotect(2);
    return out;
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
