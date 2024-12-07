#include "cheapr.h"

// A more memory-efficient which()
// Author: Nick Christofides

// Count the number of true values

R_xlen_t count_true(int *px, R_xlen_t n){
  R_xlen_t size = 0;
  if (n >= CHEAPR_OMP_THRESHOLD){
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
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);
  bool is_long = (n > integer_max_);
  if (Rf_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  if (Rf_isVectorList(x)){
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  SEXP val_is_na = Rf_protect(cpp_is_na(value)); ++NP;
  if (Rf_asLogical(val_is_na)){
    Rf_unprotect(NP);
    if (invert){
      return cpp_which_not_na(x);
    } else {
      return cpp_which_na(x);
    }
  }

  if (implicit_na_coercion(value, x)){
    Rf_unprotect(NP);
    Rf_error("Value has been implicitly converted to NA, please check");
  }

  R_xlen_t n_vals = scalar_count(x, value, false);
  R_xlen_t out_size = invert ? n - n_vals : n_vals;
  R_xlen_t whichi = 0;
  R_xlen_t i = 0;
  switch ( CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
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
    Rf_unprotect(NP);
    return out;
  }
  case REALSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, REALSXP)); ++NP;
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
    Rf_unprotect(NP);
    return out;
  }
  case CHEAPR_INT64SXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, CHEAPR_INT64SXP)); ++NP;
    long long int val = INTEGER64_PTR(value)[0];
    long long int *p_x = INTEGER64_PTR(x);
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
    Rf_unprotect(NP);
    return out;
  }
  case STRSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, STRSXP)); ++NP;
    SEXP val = Rf_protect(Rf_asChar(value)); ++NP;
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
    Rf_unprotect(NP);
    return out;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}

// Memory-efficient which(is.na(x))

[[cpp11::register]]
SEXP cpp_which_na(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  bool is_short = (n <= integer_max_);
  switch ( CHEAPR_TYPEOF(x) ){
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
  case CHEAPR_INT64SXP: {
    R_xlen_t count = na_count(x, true);
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
  }
  case REALSXP: {
    R_xlen_t count = na_count(x, true);
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
  switch ( CHEAPR_TYPEOF(x) ){
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
  case CHEAPR_INT64SXP: {
    R_xlen_t count = na_count(x, true);
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
  }
  case REALSXP: {
    R_xlen_t count = na_count(x, true);
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

// Return the locations of T, F, and NA in one pass
// Must provide the correct num of T and F as args

[[cpp11::register]]
SEXP cpp_lgl_locs(SEXP x, R_xlen_t n_true, R_xlen_t n_false,
                  bool include_true, bool include_false, bool include_na){
  R_xlen_t n = Rf_xlength(x);
  int *p_x = LOGICAL(x);

  if (n > integer_max_){
    SEXP true_locs = Rf_protect(Rf_allocVector(REALSXP, include_true ? n_true : 0));
    SEXP false_locs = Rf_protect(Rf_allocVector(REALSXP, include_false ? n_false : 0));
    SEXP na_locs = Rf_protect(Rf_allocVector(REALSXP, include_na ? (n - n_true - n_false) : 0));

    double *p_true = REAL(true_locs);
    double *p_false = REAL(false_locs);
    double *p_na = REAL(na_locs);

    R_xlen_t k1 = 0;
    R_xlen_t k2 = 0;
    R_xlen_t k3 = 0;

    for (R_xlen_t i = 0; i < n; ++i){
      if (include_true && p_x[i] == TRUE){
        p_true[k1++] = i + 1;
      } else if (include_false && p_x[i] == FALSE){
        p_false[k2++] = i + 1;
      } else if (include_na && p_x[i] == NA_LOGICAL){
        p_na[k3++] = i + 1;
      }
    }
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, true_locs);
    SET_VECTOR_ELT(out, 1, false_locs);
    SET_VECTOR_ELT(out, 2, na_locs);

    SEXP names = Rf_protect(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("true"));
    SET_STRING_ELT(names, 1, Rf_mkChar("false"));
    SET_STRING_ELT(names, 2, Rf_mkChar("na"));
    Rf_setAttrib(out, R_NamesSymbol, names);

    Rf_unprotect(5);
    return out;
  } else {
    SEXP true_locs = Rf_protect(Rf_allocVector(INTSXP, include_true ? n_true : 0));
    SEXP false_locs = Rf_protect(Rf_allocVector(INTSXP, include_false ? n_false : 0));
    SEXP na_locs = Rf_protect(Rf_allocVector(INTSXP, include_na ? (n - n_true - n_false) : 0));

    int *p_true = INTEGER(true_locs);
    int *p_false = INTEGER(false_locs);
    int *p_na = INTEGER(na_locs);

    int k1 = 0;
    int k2 = 0;
    int k3 = 0;

    for (int i = 0; i < n; ++i){
      if (include_true && p_x[i] == TRUE){
        p_true[k1++] = i + 1;
      } else if (include_false && p_x[i] == FALSE){
        p_false[k2++] = i + 1;
      } else if (include_na && p_x[i] == NA_LOGICAL){
        p_na[k3++] = i + 1;
      }
    }
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, true_locs);
    SET_VECTOR_ELT(out, 1, false_locs);
    SET_VECTOR_ELT(out, 2, na_locs);

    SEXP names = Rf_protect(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("true"));
    SET_STRING_ELT(names, 1, Rf_mkChar("false"));
    SET_STRING_ELT(names, 2, Rf_mkChar("na"));
    Rf_setAttrib(out, R_NamesSymbol, names);

    Rf_unprotect(5);
    return out;
  }
}

// like cpp_which_val but an optional n_values arg to skip counting
SEXP cpp_val_find(SEXP x, SEXP value, bool invert, SEXP n_values){
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);
  bool is_long = (n > integer_max_);
  if (Rf_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  if (Rf_isVectorList(x)){
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }

  bool val_is_na = cpp_any_na(value, true);

  if (implicit_na_coercion(value, x)){
    Rf_unprotect(NP);
    Rf_error("Value has been implicitly converted to NA, please check");
  }

  R_xlen_t n_vals = Rf_isNull(n_values) ? scalar_count(x, value, false) : Rf_asReal(n_values);
  R_xlen_t out_size = invert ? n - n_vals : n_vals;
  R_xlen_t whichi = 0;
  R_xlen_t i = 0;
  switch ( CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
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
    Rf_unprotect(NP);
    return out;
  }
  case REALSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, REALSXP)); ++NP;
    double val = Rf_asReal(value);
    double *p_x = REAL(x);
    if (val_is_na){
      if (is_long){
        double *p_out = REAL(out);
        if (invert){
          while (whichi < out_size){
            p_out[whichi] = i + 1;
            whichi += !cheapr_is_na_dbl(p_x[i]);
            ++i;
          }
        } else {
          while (whichi < out_size){
            p_out[whichi] = i + 1;
            whichi += cheapr_is_na_dbl(p_x[i]);
            ++i;
          }
        }
      } else {
        int *p_out = INTEGER(out);
        if (invert){
          CHEAPR_WHICH_VAL_INVERTED(val);
        } else {
          CHEAPR_WHICH_VAL(val);
        }
      }
    } else {
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
    }
    Rf_unprotect(NP);
    return out;
  }
  case CHEAPR_INT64SXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, CHEAPR_INT64SXP)); ++NP;
    long long int val = INTEGER64_PTR(value)[0];
    long long int *p_x = INTEGER64_PTR(x);
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
    Rf_unprotect(NP);
    return out;
  }
  case STRSXP: {
    SEXP out = Rf_protect(Rf_allocVector(is_long ? REALSXP : INTSXP, out_size));
    ++NP;
    Rf_protect(value = coerce_vector(value, STRSXP)); ++NP;
    SEXP val = Rf_protect(Rf_asChar(value)); ++NP;
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
    Rf_unprotect(NP);
    return out;
  }
  default: {
    Rf_unprotect(NP);
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
