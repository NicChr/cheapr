#include "cheapr.h"

// Basic math operations by reference
// All NA and NaN values are ignored

void check_numeric(SEXP x){
  if (!(Rf_isNumeric(x) && !Rf_isObject(x))){
    Rf_error("x must be a numeric vector");
  }
}

void copy_warning(){
  Rf_warning("x is not a double vector and has been copied, it will not be replaced by reference.\n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_log(x)`");
}

double round_nearest_even(double x){
  x -= std::remainder(x, 1.0);
  return x;
}

SEXP check_transform_altrep(SEXP x){
  if (ALTREP(x)){
    Rf_warning("Cannot update an ALTREP by reference, a copy has been made. \n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_abs(x)`");
    return altrep_materialise(x);
  } else {
    return x;
  }
}

#define CHEAPR_MATH_INT_LOOP(_fun_)                                         \
for (R_xlen_t i = 0; i < n; ++i) {                                          \
  p_out[i] = p_out[i] == NA_INTEGER ? p_out[i] : _fun_(p_out[i]);           \
}                                                                           \

#define CHEAPR_MATH_REAL_LOOP(_fun_)                                      \
for (R_xlen_t i = 0; i < n; ++i) {                                        \
  p_out[i] = p_out[i] != p_out[i] ? p_out[i] : _fun_(p_out[i]);           \
}                                                                         \


#define CHEAPR_TRUNC(x) (std::trunc(x) + 0.0)

// Convert integer vector to plain double vector

SEXP convert_int_to_real(SEXP x){
  int *p_x = INTEGER(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = SHIELD(new_vec(REALSXP, n));
  double *p_out = REAL(out);
  for (int i = 0; i < n; ++i){
    p_out[i] = p_x[i] != NA_INTEGER ? static_cast<double>(p_x[i]) : NA_REAL;
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_abs(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = INTEGER(out);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_MATH_INT_LOOP(std::abs);
    } else {
      OMP_FOR_SIMD
      CHEAPR_MATH_INT_LOOP(std::abs);
    }
    break;
  }
  case REALSXP: {
    double *p_out = REAL(out);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_MATH_INT_LOOP(std::fabs);
    }
    else {
      OMP_FOR_SIMD
      CHEAPR_MATH_INT_LOOP(std::fabs);
    }
    break;
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_floor(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (Rf_isReal(out)){
    double *p_out = REAL(out);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_MATH_REAL_LOOP(std::floor);
    } else {
      OMP_FOR_SIMD
      CHEAPR_MATH_REAL_LOOP(std::floor);
    }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_ceiling(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (Rf_isReal(out)){
    double *p_out = REAL(out);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_MATH_REAL_LOOP(std::ceil);
    } else {
      OMP_FOR_SIMD
      CHEAPR_MATH_REAL_LOOP(std::ceil);
    }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_trunc(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (Rf_isReal(out)){
    double *p_out = REAL(out);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      CHEAPR_MATH_REAL_LOOP(CHEAPR_TRUNC);
    } else {
      OMP_FOR_SIMD
      CHEAPR_MATH_REAL_LOOP(CHEAPR_TRUNC);
    }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_change_sign(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = INTEGER(out);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        p_out[i] = (p_out[i] == NA_INTEGER) ? p_out[i] : -p_out[i];
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        p_out[i] = (p_out[i] == NA_INTEGER) ? p_out[i] : -p_out[i];
      }
    }
    break;
  }
  case REALSXP: {
    double *p_out = REAL(out);
    if (n_cores > 1){
      int n_cores = num_cores();
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        p_out[i] = (p_out[i] == p_out[i]) ? -p_out[i] : p_out[i];
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) {
        p_out[i] = (p_out[i] == p_out[i]) ? -p_out[i] : p_out[i];
      }
    }
    break;
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_exp(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  SEXP out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = REAL(out);
  if (n_cores > 1){
    OMP_PARALLEL_FOR_SIMD
    CHEAPR_MATH_REAL_LOOP(std::exp);
  } else {
    OMP_FOR_SIMD
    CHEAPR_MATH_REAL_LOOP(std::exp);
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_sqrt(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  SEXP out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = REAL(out);
  if (n_cores > 1){
    OMP_PARALLEL_FOR_SIMD
    CHEAPR_MATH_REAL_LOOP(std::sqrt);
  } else {
    OMP_FOR_SIMD
    CHEAPR_MATH_REAL_LOOP(std::sqrt);
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_add(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  int NP = 0;
  SEXP out = SHIELD(check_transform_altrep(x)); ++NP;
  R_xlen_t xn = Rf_xlength(out);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  if (xn > 0){
    if (yn > xn){
      YIELD(NP);
      Rf_error("length(y) must be <= length(x)");
    }
    if (yn == 0){
      YIELD(NP);
      Rf_error("length(y) must be be non-zero");
    }
  }

  switch (TYPEOF(out)){
  case LGLSXP:
  case INTSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(out);
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
    SHIELD(out = coerce_vec(out, REALSXP));
    ++NP;
    double *p_x = REAL(out);
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
    double *p_x = REAL(out);
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
    double *p_x = REAL(out);
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
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_subtract(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  int NP = 0;
  SEXP out = SHIELD(check_transform_altrep(x));
  ++NP;
  R_xlen_t xn = Rf_xlength(out);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (xn > 0){
    if (yn > xn){
      YIELD(NP);
      Rf_error("length(y) must be <= length(x)");
    }
    if (yn == 0){
      YIELD(NP);
      Rf_error("length(y) must be be non-zero");
    }
  }
  switch (TYPEOF(out)){
  case LGLSXP:
  case INTSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(out);
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
    SHIELD(out = coerce_vec(out, REALSXP));
    ++NP;
    double *p_x = REAL(out);
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
    double *p_x = REAL(out);
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
    double *p_x = REAL(out);
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
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_multiply(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  int NP = 0;
  SEXP out = SHIELD(check_transform_altrep(x));
  ++NP;
  R_xlen_t xn = Rf_xlength(out);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (xn > 0){
    if (yn > xn){
      YIELD(NP);
      Rf_error("length(y) must be <= length(x)");
    }
    if (yn == 0){
      YIELD(NP);
      Rf_error("length(y) must be be non-zero");
    }
  }
  switch (TYPEOF(out)){
  case LGLSXP:
  case INTSXP: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(out);
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
    SHIELD(out = coerce_vec(out, REALSXP));
    ++NP;
    double *p_x = REAL(out);
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
    // SHIELD(y = coerce_vec(y, REALSXP));
    // ++NP;
    double *p_x = REAL(out);
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
    double *p_x = REAL(out);
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
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_divide(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (xn > 0){
    if (yn > xn){
      Rf_error("length(y) must be <= length(x)");
    }
    if (yn == 0){
      Rf_error("length(y) must be be non-zero");
    }
  }
  SEXP out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    double *p_x = REAL(out);
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
    double *p_x = REAL(out);
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
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_pow(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  int n_cores = xn >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (xn > 0){
    if (yn > xn){
      Rf_error("length(y) must be <= length(x)");
    }
    if (yn == 0){
      Rf_error("length(y) must be be non-zero");
    }
  }
  SEXP out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  switch (TYPEOF(y)){
  case INTSXP: {
    double *p_x = REAL(out);
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
    double *p_x = REAL(out);
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
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_log(SEXP x, SEXP base){
  check_numeric(x);
  check_numeric(base);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t basen = Rf_xlength(base);
  if (xn > 0){
    if (basen > xn){
      Rf_error("length(base) must be <= length(x)");
    }
    if (basen == 0){
      Rf_error("length(base) must be be non-zero");
    }
  }
  int n_cores = xn >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  SEXP out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_x = REAL(out);
  double *p_base = REAL(base);
  if (n_cores > 1){
    OMP_PARALLEL_FOR_SIMD
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t basei = i % basen;
      p_x[i] = (p_x[i] != p_x[i] || p_base[basei] != p_base[basei]) ?
      NA_REAL : std::log(p_x[i]) / std::log(p_base[basei]);
    }
  } else {
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t basei = i % basen;
      p_x[i] = (p_x[i] != p_x[i] || p_base[basei] != p_base[basei]) ?
      NA_REAL : std::log(p_x[i]) / std::log(p_base[basei]);
    }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_round(SEXP x, SEXP digits){
  check_numeric(x);
  check_numeric(digits);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t xn = Rf_xlength(out);
  R_xlen_t digitsn = Rf_xlength(digits);
  if (xn > 0){
    if (digitsn > xn){
      Rf_error("length(digits) must be <= length(x)");
    }
    if (digitsn == 0){
      Rf_error("length(digits) must be be non-zero");
    }
  }
  int n_cores = xn >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  // We don't need to round integers.
  if (Rf_isReal(out)){
    switch (TYPEOF(digits)){
    case INTSXP: {
      double *p_x = REAL(out);
      int *p_digits = INTEGER(digits);
      if (n_cores > 1){
        OMP_PARALLEL_FOR_SIMD
        for (R_xlen_t i = 0; i < xn; ++i) {
          R_xlen_t digitsi = i % digitsn;
          if ( (p_x[i] == p_x[i] && p_digits[digitsi] != NA_INTEGER) ){
            double tempx = p_x[i];
            int tempdig = p_digits[digitsi];
            double mfactor = std::pow(10, tempdig);
            tempx *= mfactor;
            tempx = round_nearest_even(tempx);
            tempx /= mfactor;
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
            double mfactor = std::pow(10, tempdig);
            tempx *= mfactor;
            tempx = round_nearest_even(tempx);
            tempx /= mfactor;
            p_x[i] = tempx;
          } else {
            p_x[i] = NA_REAL;
          }
        }
      }
      break;
    }
    default: {
      double *p_x = REAL(out);
      double *p_digits = REAL(digits);
      if (n_cores > 1){
        OMP_PARALLEL_FOR_SIMD
        for (R_xlen_t i = 0; i < xn; ++i) {
          R_xlen_t digitsi = i % digitsn;
          if ( (p_x[i] == p_x[i] && p_digits[digitsi] == p_digits[digitsi]) ){
            double tempx = p_x[i];
            int tempdig = p_digits[digitsi];
            double mfactor = std::pow(10, tempdig);
            tempx *= mfactor;
            tempx = round_nearest_even(tempx);
            tempx /= mfactor;
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
            double mfactor = std::pow(10, tempdig);
            tempx *= mfactor;
            tempx = round_nearest_even(tempx);
            tempx /= mfactor;
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
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_int_sign(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = SHIELD(new_vec(INTSXP, n));
  int* RESTRICT p_out = INTEGER(out);
  int res;
  switch (TYPEOF(x)){
  case LGLSXP: {
    int* RESTRICT p_x = INTEGER(x);
    memmove(&p_out[0], &p_x[0], n * sizeof(int));
    break;
  }
  case INTSXP: {
    int* RESTRICT p_x = INTEGER(x);
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) {
      res = p_x[i] == NA_INTEGER ? NA_INTEGER : (p_x[i] > 0) - (p_x[i] < 0);
      p_out[i] = res;
    }
    break;
  }
  case REALSXP: {
    double* RESTRICT p_x = REAL(x);
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) {
      res = p_x[i] != p_x[i] ? NA_INTEGER : (p_x[i] > 0) - (p_x[i] < 0);
      p_out[i] = res;
    }
    break;
  }
  }
  YIELD(1);
  return out;
}
