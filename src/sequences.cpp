#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>

double r_sum(SEXP x, bool na_rm = false){
  cpp11::function base_sum = cpp11::package("base")["sum"];
  return Rf_asReal(base_sum(x, cpp11::named_arg("na.rm") = na_rm));
}

double r_min(SEXP x){
  cpp11::function base_min = cpp11::package("base")["min"];
  double out = R_PosInf;
  if (Rf_xlength(x) > 0){
    out = Rf_asReal(base_min(x));
  }
  return out;
}

// My version of base::sequence()

[[cpp11::register]]
SEXP cpp_int_sequence(SEXP size, SEXP from, SEXP by) {
  int size_n = Rf_length(size);
  int from_n = Rf_length(from);
  int by_n = Rf_length(by);
  if (size_n > 0 && (from_n <= 0 || by_n <= 0)){
    Rf_error("from and by must both have length > 0");
  }
  double out_size = r_sum(size, false);
  double min_size = r_min(size);
  if (!(out_size == out_size)){
    Rf_error("size must not contain NA values");
  }
  if (min_size < 0){
    Rf_error("size must be a vector of non-negative integers");
  }
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
  int *p_out = INTEGER(out);
  R_xlen_t index = 0;
  int fj;
  int bj;
  int start;
  int increment;
  int seq_size;
  double seq_end;
  if (size_n > 0){
    int *p_size = INTEGER(size);
    int *p_from = INTEGER(from);
    int *p_by = INTEGER(by);
    for (int j = 0; j < size_n; ++j){
      seq_size = p_size[j];
      fj = j % from_n;
      bj = j % by_n;
      start = p_from[fj];
      increment = p_by[bj];
      // Throw error if integer overflow
      seq_end = ( (std::fmax(seq_size - 1, 0.0)) * increment ) + start;
      if (std::fabs(seq_end) > integer_max_){
        Rf_unprotect(1);
        Rf_error("Integer overflow value of %g in sequence %d", seq_end, j + 1);
      }
      if (start == NA_INTEGER){
        Rf_unprotect(1);
        Rf_error("from contains NA values");
      }
      if (increment == NA_INTEGER){
        Rf_unprotect(1);
        Rf_error("by contains NA values");
      }
      for (int i = 0; i < seq_size; ++i){
        p_out[index] = start;
        start += increment;
        ++index;
      }
    }
  }
  Rf_unprotect(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_dbl_sequence(SEXP size, SEXP from, SEXP by) {
  int size_n = Rf_length(size);
  int from_n = Rf_length(from);
  int by_n = Rf_length(by);
  if (size_n > 0 && (from_n <= 0 || by_n <= 0)){
    Rf_error("from and by must both have length > 0");
  }
  // To recycle we would need to do sum * remainder of the sum over n
  double out_size = r_sum(size, false);
  double min_size = r_min(size);
  if (!(out_size == out_size)){
    Rf_error("size must not contain NA values");
  }
  if (min_size < 0){
    Rf_error("size must be a vector of non-negative integers");
  }
  SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
  double *p_out = REAL(out);
  R_xlen_t index = 0;
  int fj;
  int bj;
  int seq_size;
  double start;
  double increment;
  if (size_n > 0){
    int *p_size = INTEGER(size);
    double *p_from = REAL(from);
    double *p_by = REAL(by);
    for (int j = 0; j < size_n; ++j){
      seq_size = p_size[j];
      fj = j % from_n;
      bj = j % by_n;
      start = p_from[fj];
      increment = p_by[bj];
      if (!(start == start)){
        Rf_unprotect(1);
        Rf_error("from contains NA values");
      }
      if (!(increment == increment)){
        Rf_unprotect(1);
        Rf_error("by contains NA values");
      }
      for (int i = 0; i < seq_size; ++i){
        p_out[index] = ( start + (i * increment) );
        ++index;
      }
    }
  }
  Rf_unprotect(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_sequence(SEXP size, SEXP from, SEXP by) {
  int size_n = Rf_length(size);
  int from_n = Rf_length(from);
  int by_n = Rf_length(by);
  switch (TYPEOF(from)){
  case INTSXP: {
    switch (TYPEOF(by)){
  case INTSXP: {
    int n = std::max(std::max(size_n, from_n), by_n);
    double seq_end;
    bool out_is_integer = true;
    int *p_size = INTEGER(size);
    int *p_from = INTEGER(from);
    int *p_by = INTEGER(by);

    // Checking that the sequence values are integers
    // Only do the loop if vectors are not zero-length
    if (size_n > 0 && from_n > 0 && by_n > 0){
      for (int i = 0; i < n; ++i){
        seq_end = (std::fmax(p_size[i % size_n] - 1.0, 0.0) * p_by[i % by_n]) + p_from[i % from_n];
        if (seq_end > integer_max_){
          out_is_integer = false;
          break;
        }
      }
    }
    // If all sequence values are < 2^31 then we can safely use cpp_int_sequence
    if (out_is_integer){
      return cpp_int_sequence(size, from, by);
    } else {
      Rf_protect(from = Rf_coerceVector(from, REALSXP));
      Rf_protect(by = Rf_coerceVector(by, REALSXP));
      SEXP out = Rf_protect(cpp_dbl_sequence(size, from, by));
      Rf_unprotect(3);
      return out;
    }
  }
  case REALSXP: {
    Rf_protect(from = Rf_coerceVector(from, REALSXP));
    SEXP out = Rf_protect(cpp_dbl_sequence(size, from, by));
    Rf_unprotect(2);
    return out;
  }
  default: {
    Rf_error("by must have type integer or double in %s", __func__);
  }
  }
    break;
  }
  case REALSXP: {
    switch (TYPEOF(by)){
  case INTSXP: {
    Rf_protect(by = Rf_coerceVector(by, REALSXP));
    SEXP out = Rf_protect(cpp_dbl_sequence(size, from, by));
    Rf_unprotect(2);
    return out;
  }
  case REALSXP: {
    return cpp_dbl_sequence(size, from, by);
  }
  default: {
    Rf_error("by must have type integer or double in %s", __func__);
  }
  }
  }
  default: {
    Rf_error("from must have type integer or double in %s", __func__);
  }
  }
}

[[cpp11::register]]
SEXP cpp_window_sequence(SEXP size,
                         double k,
                         bool partial = true,
                         bool ascending = true) {
  int size_n = Rf_length(size);
  SEXP size_sexp = Rf_protect(Rf_coerceVector(size, INTSXP));
  if (r_min(size_sexp) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  k = std::fmax(k, 0);
  R_xlen_t N = r_sum(size_sexp, false);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, N));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size_sexp);
  R_xlen_t index = 0;
  if (ascending){
    // right aligned window sequences
    if (partial){
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          if (i < k){
            p_out[index] = i + 1;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    } else {
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          if (i < (k - 1)){
            p_out[index] = NA_INTEGER;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    }
  } else {
    // left aligned window sequences
    int idiff;
    if (partial){
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          idiff = p_size[j] - i - 1;
          if (idiff < k){
            p_out[index] = idiff + 1;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    } else {
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          idiff = p_size[j] - i - 1;
          if (idiff < (k - 1)){
            p_out[index] = NA_INTEGER;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    }
  }
  Rf_unprotect(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_lag_sequence(SEXP size, double k, bool partial = false) {
  Rf_protect(size = Rf_coerceVector(size, INTSXP));
  if (r_min(size) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  int size_n = Rf_length(size);
  k = std::fmax(k, 0);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, r_sum(size, false)));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size);
  R_xlen_t index = 0;
  if (partial){
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        if (i < k){
          p_out[index] = i;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  } else {
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        if (i < k){
          p_out[index] = NA_INTEGER;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  }
  Rf_unprotect(2);
  return out;
}
[[cpp11::register]]
SEXP cpp_lead_sequence(SEXP size, double k, bool partial = false) {
  Rf_protect(size = Rf_coerceVector(size, INTSXP));
  if (r_min(size) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  int size_n = Rf_length(size);
  k = std::fmax(k, 0);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, r_sum(size, false)));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size);
  R_xlen_t index = 0;
  int idiff;
  if (partial){
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        idiff = p_size[j] - i - 1;
        if (idiff < k){
          p_out[index] = idiff;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  } else {
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        idiff = p_size[j] - i - 1;
        if (idiff < k){
          p_out[index] = NA_INTEGER;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  }
  Rf_unprotect(2);
  return out;
}
[[cpp11::register]]
SEXP cpp_sequence_id(SEXP size){
  int size_n = Rf_length(size);
  SEXP size_sexp = Rf_protect(Rf_coerceVector(size, INTSXP));
  if (r_min(size_sexp) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  R_xlen_t N = r_sum(size_sexp, false);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, N));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size_sexp);
  R_xlen_t k = 0;
  int seq_size;
  for (int i = 0; i < size_n; ++i){
    seq_size = p_size[i];
    for (int j = 0; j < seq_size; ++j){
      p_out[k++] = i + 1;
    }
  }
  Rf_unprotect(2);
  return out;
}
