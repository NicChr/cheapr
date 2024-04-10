#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>

// Basic math operations by reference
// All NA and NaN values are ignored

void cpp_check_numeric(SEXP x){
  if (!(Rf_isNumeric(x) && !Rf_isObject(x))){
    Rf_error("x must be a numeric vector");
  }
}

[[cpp11::register]]
SEXP cpp_set_abs(SEXP x){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  switch (TYPEOF(x)){
  case INTSXP: {
    int *p_x = INTEGER(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] != NA_INTEGER) p_x[i] = std::abs(p_x[i]);
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] != NA_INTEGER) p_x[i] = std::abs(p_x[i]);
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::fabs(p_x[i]);
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::fabs(p_x[i]);
      }
    }
    break;
  }
  }
  return x;
}

[[cpp11::register]]
SEXP cpp_set_floor(SEXP x){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  if (Rf_isReal(x)){
    double *p_x = REAL(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::floor(p_x[i]);
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::floor(p_x[i]);
      }
    }
  }
  return x;
}

[[cpp11::register]]
SEXP cpp_set_ceiling(SEXP x){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  if (Rf_isReal(x)){
    double *p_x = REAL(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::ceil(p_x[i]);
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::ceil(p_x[i]);
      }
    }
  }
  return x;
}

[[cpp11::register]]
SEXP cpp_set_trunc(SEXP x){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  if (Rf_isReal(x)){
    double *p_x = REAL(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::trunc(p_x[i]);
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::trunc(p_x[i]);
      }
    }
  }
  return x;
}

[[cpp11::register]]
SEXP cpp_set_round(SEXP x, int digits){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  if (Rf_isReal(x)){
    double *p_x = REAL(x);
    if (n_cores > 1){
      int n_cores = num_cores();
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]){
          double temp = p_x[i];
          temp *= std::pow(10, digits);
          temp = std::round(temp);
          temp *= std::pow(10, -(digits));
          p_x[i] = temp;
        }
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]){
          double temp = p_x[i];
          temp *= std::pow(10, digits);
          temp = std::round(temp);
          // temp -= std::remainder(temp, 1.0);
          temp *= std::pow(10, -(digits));
          // temp = std::nearbyint( temp * 0.5 ) * 2.0 * std::pow(10, -digits);
          p_x[i] = temp;
        }
      }
    }
  }
  return x;
}

[[cpp11::register]]
SEXP cpp_set_change_sign(SEXP x){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  switch (TYPEOF(x)){
  case INTSXP: {
    int *p_x = INTEGER(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] != NA_INTEGER) p_x[i] = -p_x[i];
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] != NA_INTEGER) p_x[i] = -p_x[i];
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    if (n_cores > 1){
      int n_cores = num_cores();
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = -p_x[i];
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = -p_x[i];
      }
    }
    break;
  }
  }
  return x;
}

[[cpp11::register]]
SEXP cpp_set_exp(SEXP x){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  Rf_protect(x = Rf_coerceVector(x, REALSXP));
  double *p_x = REAL(x);
  if (n_cores > 1){
    OMP_PARALLEL_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) {
      if (p_x[i] == p_x[i]) p_x[i] = std::exp(p_x[i]);
    }
  } else {
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) {
      if (p_x[i] == p_x[i]) p_x[i] = std::exp(p_x[i]);
    }
  }
  Rf_unprotect(1);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_sqrt(SEXP x){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  Rf_protect(x = Rf_coerceVector(x, REALSXP));
  double *p_x = REAL(x);
  if (n_cores > 1){
    OMP_PARALLEL_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) {
      if (p_x[i] == p_x[i]) p_x[i] = std::sqrt(p_x[i]);
    }
  } else {
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) {
      if (p_x[i] == p_x[i]) p_x[i] = std::sqrt(p_x[i]);
    }
  }
  Rf_unprotect(1);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_log(SEXP x, double base){
  cpp_check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  Rf_protect(x = Rf_coerceVector(x, REALSXP));
  double *p_x = REAL(x);
  if (n_cores > 1){
    if (base == std::exp(1)){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::log(p_x[i]);
      }
    } else {
      double log_base = std::log(base);
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = (std::log(p_x[i]) / log_base);
      }
    }
  } else {
    if (base == std::exp(1)){
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = std::log(p_x[i]);
      }
    } else {
      double log_base = std::log(base);
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_x[i] == p_x[i]) p_x[i] = (std::log(p_x[i]) / log_base);
      }
    }
  }
  Rf_unprotect(1);
  return x;
}
