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

void copy_warning(){
  Rf_warning("x is not a double vector and has been copied, it will not be replaced by reference.\n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_log(x)`");
}

// Source: https://stackoverflow.com/questions/32746523/ieee-754-compliant-round-half-to-even
double round_nearest_even(double x){
  x -= std::remainder(x, 1.0);
  return x;
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
    }
    else {
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
  if (!Rf_isReal(x)) copy_warning();
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
  if (!Rf_isReal(x)) copy_warning();
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
SEXP cpp_set_add(SEXP x, SEXP y){
  cpp_check_numeric(x);
  cpp_check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= 100000 ? num_cores() : 1;
  if (yn > xn){
    Rf_error("length(y) must be <= length(x)");
  }
  int n_protections = 0;
  switch (TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] == NA_INTEGER || p_y[yi] == NA_INTEGER) ?
      NA_INTEGER : p_x[i] + p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    Rf_protect(x = Rf_coerceVector(x, REALSXP));
    ++n_protections;
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
      NA_REAL : p_x[i] + p_y[yi];
    }
    break;
  }
  }
    break;
  }
  case REALSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    // Rf_protect(y = Rf_coerceVector(y, REALSXP));
    // ++n_protections;
    double *p_x = REAL(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] !=  p_x[i] || p_y[yi] == NA_INTEGER) ?
      NA_REAL : p_x[i] + p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
      NA_REAL : p_x[i] + p_y[yi];
    }
    break;
  }
  }
  }
  }
  Rf_unprotect(n_protections);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_subtract(SEXP x, SEXP y){
  cpp_check_numeric(x);
  cpp_check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= 100000 ? num_cores() : 1;
  if (yn > xn){
    Rf_error("length(y) must be <= length(x)");
  }
  int n_protections = 0;
  switch (TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] == NA_INTEGER || p_y[yi] == NA_INTEGER) ?
      NA_INTEGER : p_x[i] - p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    Rf_protect(x = Rf_coerceVector(x, REALSXP));
    ++n_protections;
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
      NA_REAL : p_x[i] - p_y[yi];
    }
    break;
  }
  }
    break;
  }
  case REALSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    // Rf_protect(y = Rf_coerceVector(y, REALSXP));
    // ++n_protections;
    double *p_x = REAL(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] !=  p_x[i] || p_y[yi] == NA_INTEGER) ?
      NA_REAL : p_x[i] - p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
      NA_REAL : p_x[i] - p_y[yi];
    }
    break;
  }
  }
  }
  }
  Rf_unprotect(n_protections);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_multiply(SEXP x, SEXP y){
  cpp_check_numeric(x);
  cpp_check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= 100000 ? num_cores() : 1;
  if (yn > xn){
    Rf_error("length(y) must be <= length(x)");
  }
  int n_protections = 0;
  switch (TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] == NA_INTEGER || p_y[yi] == NA_INTEGER) ?
      NA_INTEGER : p_x[i] * p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    Rf_protect(x = Rf_coerceVector(x, REALSXP));
    ++n_protections;
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
      NA_REAL : p_x[i] * p_y[yi];
    }
    break;
  }
  }
    break;
  }
  case REALSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    // Rf_protect(y = Rf_coerceVector(y, REALSXP));
    // ++n_protections;
    double *p_x = REAL(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] !=  p_x[i] || p_y[yi] == NA_INTEGER) ?
      NA_REAL : p_x[i] * p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
      NA_REAL : p_x[i] * p_y[yi];
    }
    break;
  }
  }
  }
  }
  Rf_unprotect(n_protections);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_divide(SEXP x, SEXP y){
  cpp_check_numeric(x);
  cpp_check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= 100000 ? num_cores() : 1;
  if (yn > xn){
    Rf_error("length(y) must be <= length(x)");
  }
  if (!Rf_isReal(x)) copy_warning();
  Rf_protect(x = Rf_coerceVector(x, REALSXP));
  switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    double *p_x = REAL(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] !=  p_x[i] || p_y[yi] == NA_INTEGER) ?
      NA_REAL : p_x[i] / p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
      NA_REAL : p_x[i] / p_y[yi];
    }
    break;
  }
  }
  Rf_unprotect(1);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_pow(SEXP x, SEXP y){
  cpp_check_numeric(x);
  cpp_check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  int n_protections = 0;
  int n_cores = xn >= 100000 ? num_cores() : 1;
  if (yn > xn){
    Rf_error("length(y) must be <= length(x)");
  }
  if (!Rf_isReal(x)) copy_warning();
  Rf_protect(x = Rf_coerceVector(x, REALSXP));
  ++n_protections;
  switch (TYPEOF(y)){
  case INTSXP: {
    double *p_x = REAL(x);
    int *p_y = INTEGER(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      if (p_x[i] == 1.0 || p_y[yi] == 0){
        p_x[i] = 1.0;
      } else {
        p_x[i] = (p_x[i] !=  p_x[i] || p_y[yi] == NA_INTEGER) ?
        NA_REAL : std::pow(p_x[i], p_y[yi]);
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    double *p_y = REAL(y);
#pragma omp parallel for simd num_threads(n_cores) if (n_cores > 1)
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t yi = i % yn;
      if (p_x[i] == 1.0 || p_y[yi] == 0.0){
        p_x[i] = 1.0;
      } else {
        p_x[i] = (p_x[i] != p_x[i]  || p_y[yi] != p_y[yi]) ?
        NA_REAL : std::pow(p_x[i], p_y[yi]);
      }
    }
    break;
  }
  }
  Rf_unprotect(n_protections);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_log(SEXP x, SEXP base){
  cpp_check_numeric(x);
  cpp_check_numeric(base);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t basen = Rf_xlength(base);
  if (basen > xn){
    Rf_error("length(y) must be <= length(x)");
  }
  int n_cores = xn >= 100000 ? num_cores() : 1;
  if (!Rf_isReal(x)) copy_warning();
  Rf_protect(x = Rf_coerceVector(x, REALSXP));
  double *p_x = REAL(x);
  double *p_base = REAL(base);
  if (n_cores > 1){
    OMP_PARALLEL_FOR_SIMD
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t basei = i % basen;
      p_x[i] = (p_x[i] != p_x[i]  || p_base[basei] != p_base[basei]) ?
      NA_REAL : std::log(p_x[i]) / std::log(p_base[basei]);
    }
  } else {
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t basei = i % basen;
      p_x[i] = (p_x[i] != p_x[i]  || p_base[basei] != p_base[basei]) ?
      NA_REAL : std::log(p_x[i]) / std::log(p_base[basei]);
    }
  }
  Rf_unprotect(1);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_round(SEXP x, SEXP digits){
  cpp_check_numeric(x);
  cpp_check_numeric(digits);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t digitsn = Rf_xlength(digits);
  if (digitsn > xn){
    Rf_error("length(y) must be <= length(x)");
  }
  int n_cores = xn >= 100000 ? num_cores() : 1;
  // We don't need to round integers.
  if (Rf_isReal(x)){
    switch (TYPEOF(digits)){
    case INTSXP: {
      double *p_x = REAL(x);
      int *p_digits = INTEGER(digits);
      if (n_cores > 1){
        OMP_PARALLEL_FOR_SIMD
        for (R_xlen_t i = 0; i < xn; ++i) {
          R_xlen_t digitsi = i % digitsn;
          if ( (p_x[i] == p_x[i] && p_digits[digitsi] != NA_INTEGER) ){
            double tempx = p_x[i];
            int tempdig = p_digits[digitsi];
            tempx *= std::pow(10.0, tempdig);
            tempx = round_nearest_even(tempx);
            tempx *= std::pow(10.0, -(tempdig));
            p_x[i] = tempx;
          } else {
            p_x[i] = NA_REAL;
          }
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < xn; ++i) {
          R_xlen_t digitsi = i % digitsn;
          if ( (p_x[i] == p_x[i] && p_digits[digitsi] != NA_INTEGER) ){
            double tempx = p_x[i];
            int tempdig = p_digits[digitsi];
            tempx *= std::pow(10.0, tempdig);
            tempx = round_nearest_even(tempx);
            tempx *= std::pow(10.0, -(tempdig));
            p_x[i] = tempx;
          } else {
            p_x[i] = NA_REAL;
          }
        }
      }
      break;
    }
    default: {
      double *p_x = REAL(x);
      double *p_digits = REAL(digits);
      if (n_cores > 1){
        OMP_PARALLEL_FOR_SIMD
        for (R_xlen_t i = 0; i < xn; ++i) {
          R_xlen_t digitsi = i % digitsn;
          if ( (p_x[i] == p_x[i] && p_digits[digitsi] == p_digits[digitsi]) ){
            double tempx = p_x[i];
            int tempdig = p_digits[digitsi];
            tempx *= std::pow(10, tempdig);
            tempx = round_nearest_even(tempx);
            tempx *= std::pow(10, -(tempdig));
            p_x[i] = tempx;
          } else {
            p_x[i] = NA_REAL;
          }
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < xn; ++i) {
          R_xlen_t digitsi = i % digitsn;
          if ( (p_x[i] == p_x[i] && p_digits[digitsi] == p_digits[digitsi]) ){
            double tempx = p_x[i];
            int tempdig = p_digits[digitsi];
            tempx *= std::pow(10, tempdig);
            tempx = round_nearest_even(tempx);
            tempx *= std::pow(10, -(tempdig));
            p_x[i] = tempx;
          } else {
            p_x[i] = NA_REAL;
          }
        }
        break;
      }
    }
    }
  }
  return x;
}
