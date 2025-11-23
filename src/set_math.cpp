#include "cheapr.h"

using namespace cpp11;

// Basic math operations by reference
// All NA and NaN values are ignored

void check_numeric(SEXP x){
  if (!(Rf_isNumeric(x) && !Rf_isObject(x))){
    stop("x must be a numeric vector");
  }
}

void copy_warning(){
  warning("x is not a double vector and has been copied, it will not be replaced by reference.\n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_log(x)`");
}

SEXP check_transform_altrep(SEXP x){
  if (is_altrep(x)){
    warning("Cannot update an ALTREP by reference, a copy has been made. \n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_abs(x)`");
    return altrep_materialise(x);
  } else {
    return x;
  }
}

#define CHEAPR_MATH_LOOP(FUN)                                               \
for (R_xlen_t i = 0; i < n; ++i) {                                          \
  p_out[i] = is_r_na(p_out[i]) ? p_out[i] : FUN(p_out[i]);                  \
}

#define CHEAPR_PARALLEL_MATH_LOOP(FUN)                         \
if (n_cores > 1){                                              \
  OMP_PARALLEL_FOR_SIMD                                        \
  CHEAPR_MATH_LOOP(FUN);                                       \
} else {                                                       \
  OMP_FOR_SIMD                                                 \
  CHEAPR_MATH_LOOP(FUN);                                       \
}

// Convert integer vector to plain double vector

SEXP convert_int_to_real(SEXP x){
  int *p_x = INTEGER(x);
  R_xlen_t n = Rf_xlength(x);
  writable::doubles out(n);
  double* RESTRICT p_out = REAL(out);
  for (int i = 0; i < n; ++i){
    p_out[i] = as_double(p_x[i]);
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_abs(SEXP x){
  check_numeric(x);
  sexp out = check_transform_altrep(x);
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = INTEGER(out);
    CHEAPR_PARALLEL_MATH_LOOP(std::abs)
    break;
  }
  default: {
    double *p_out = REAL(out);
    CHEAPR_PARALLEL_MATH_LOOP(std::abs)
    break;
  }
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_floor(SEXP x){
  check_numeric(x);
  sexp out = check_transform_altrep(x);
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (Rf_isReal(out)){
    double *p_out = REAL(out);
    CHEAPR_PARALLEL_MATH_LOOP(std::floor)
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_ceiling(SEXP x){
  check_numeric(x);
  sexp out = check_transform_altrep(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (Rf_isReal(out)){
    double *p_out = REAL(out);
    CHEAPR_PARALLEL_MATH_LOOP(std::ceil)
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_trunc(SEXP x){
  check_numeric(x);
  sexp out = check_transform_altrep(x);
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  if (Rf_isReal(out)){
    double *p_out = REAL(out);
    CHEAPR_PARALLEL_MATH_LOOP(cheapr::trunc)
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_change_sign(SEXP x){
  check_numeric(x);
  sexp out = check_transform_altrep(x);
  R_xlen_t n = Rf_xlength(out);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = INTEGER(out);
    CHEAPR_PARALLEL_MATH_LOOP(cheapr::negate)
    break;
  }
  case REALSXP: {
    double *p_out = REAL(out);
    CHEAPR_PARALLEL_MATH_LOOP(cheapr::negate)
    break;
  }
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_exp(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;


  sexp out;

  if (!Rf_isReal(x)){
    copy_warning();
    out = convert_int_to_real(x);
  } else {
    out = check_transform_altrep(x);
  }

  double *p_out = REAL(out);
  CHEAPR_PARALLEL_MATH_LOOP(std::exp)
  return out;
}

[[cpp11::register]]
SEXP cpp_set_sqrt(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  sexp out;

  if (!Rf_isReal(x)){
    copy_warning();
    out = convert_int_to_real(x);
  } else {
    out = check_transform_altrep(x);
  }
  double *p_out = REAL(out);
  CHEAPR_PARALLEL_MATH_LOOP(std::sqrt)
  return out;
}

[[cpp11::register]]
SEXP cpp_set_add(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);

  sexp out = check_transform_altrep(x);
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

  if (xn > 0){
    if (yn > xn){
      stop("length(y) must be <= length(x)");
    }
    if (yn == 0){
      stop("length(y) must be be non-zero");
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
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_INTEGER : p_x[i] + p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    out = coerce_vec(out, REALSXP);
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi]))?
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
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_REAL : p_x[i] + p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] + p_y[yi];
    }
    break;
  }
  }
  }
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_subtract(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);

  sexp out = check_transform_altrep(x);
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

  if (xn > 0){
    if (yn > xn){
      stop("length(y) must be <= length(x)");
    }
    if (yn == 0){
      stop("length(y) must be be non-zero");
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
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_INTEGER : p_x[i] - p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    out = coerce_vec(out, REALSXP);
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi]))?
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
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_REAL : p_x[i] - p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] - p_y[yi];
    }
    break;
  }
  }
  }
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_multiply(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);

  sexp out = check_transform_altrep(x);
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

  if (xn > 0){
    if (yn > xn){
      stop("length(y) must be <= length(x)");
    }
    if (yn == 0){
      stop("length(y) must be be non-zero");
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

    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_INTEGER : p_x[i] * p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    out = coerce_vec(out, REALSXP);
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
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
    double *p_x = REAL(out);
    int *p_y = INTEGER(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_REAL : p_x[i] * p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) ||is_r_na(p_y[yi])) ?
      NA_REAL : p_x[i] * p_y[yi];
    }
    break;
  }
  }
  }
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_divide(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  uint_fast64_t xn = Rf_xlength(x);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

  if (xn > 0){
    if (yn > xn){
      stop("length(y) must be <= length(x)");
    }
    if (yn == 0){
      stop("length(y) must be be non-zero");
    }
  }

  sexp out;

  if (!Rf_isReal(x)){
    copy_warning();
    out = convert_int_to_real(x);
  } else {
    out = check_transform_altrep(x);
  }
  switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    double *p_x = REAL(out);
    int *p_y = INTEGER(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] / p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] / p_y[yi];
    }
    break;
  }
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_pow(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  R_xlen_t yi = 0;
  if (xn > 0){
    if (yn > xn){
      stop("length(y) must be <= length(x)");
    }
    if (yn == 0){
      stop("length(y) must be be non-zero");
    }
  }
  sexp out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = convert_int_to_real(x);
  } else {
    out = check_transform_altrep(x);
  }
  switch (TYPEOF(y)){
  case INTSXP: {
    double *p_x = REAL(out);
    int *p_y = INTEGER(y);
    for (R_xlen_t i = 0; i < xn; recycle_index(yi, yn), ++i) {
      if (p_x[i] == 1.0 || p_y[yi] == 0){
        p_x[i] = 1.0;
      } else {
        p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
        NA_REAL : std::pow(p_x[i], p_y[yi]);
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (R_xlen_t i = 0; i < xn; recycle_index(yi, yn), ++i) {
      if (p_x[i] == 1.0 || p_y[yi] == 0.0){
        p_x[i] = 1.0;
      } else {
        p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
        NA_REAL : std::pow(p_x[i], p_y[yi]);
      }
    }
    break;
  }
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_log(SEXP x, SEXP base){
  check_numeric(x);
  check_numeric(base);
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t basei = 0;
  R_xlen_t basen = Rf_xlength(base);

  if (xn > 0){
    if (basen > xn){
      stop("length(base) must be <= length(x)");
    }
    if (basen == 0){
      stop("length(base) must be be non-zero");
    }
  }
  sexp out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = convert_int_to_real(x);
  } else {
    out = check_transform_altrep(x);
  }
  double *p_x = REAL(out);
  double *p_base = REAL(base);
  for (R_xlen_t i = 0; i < xn; recycle_index(basei, basen), ++i) {
    p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_base[basei])) ?
    NA_REAL : std::log(p_x[i]) / std::log(p_base[basei]);
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_set_round(SEXP x, SEXP digits){
  check_numeric(x);
  check_numeric(digits);
  sexp out = check_transform_altrep(x);
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t digitsn = Rf_xlength(digits);
  uint_fast64_t digitsi = 0;

  double tempx, mfactor;

  if (xn > 0){
    if (digitsn > xn){
      stop("`length(digits)` must be `<= length(x)`");
    }
    if (digitsn == 0){
      stop("`length(digits)` must be be non-zero");
    }
  }

  // We don't need to round integers.

  if (Rf_isReal(out)){
    switch (TYPEOF(digits)){
    case INTSXP: {
      double *p_x = REAL(out);
      const int *p_digits = INTEGER(digits);
        for (uint_fast64_t i = 0; i < xn; recycle_index(digitsi, digitsn), ++i) {
          if ( (!is_r_na(p_x[i]) && !is_r_na(p_digits[digitsi])) ){
            tempx = p_x[i];
            mfactor = std::pow(10, p_digits[digitsi]);
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
    default: {
      double *p_x = REAL(out);
      const double *p_digits = REAL(digits);
      for (uint_fast64_t i = 0; i < xn; recycle_index(digitsi, digitsn), ++i) {
          if ( (!is_r_na(p_x[i]) && !is_r_na(p_digits[digitsi])) ){
            tempx = p_x[i];
            mfactor = std::pow(10, p_digits[digitsi]);
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
  return out;
}

[[cpp11::register]]
SEXP cpp_int_sign(SEXP x){
  check_numeric(x);
  uint_fast64_t n = Rf_xlength(x);
  sexp out = new_vec(INTSXP, n);
  int* RESTRICT p_out = INTEGER(out);
  switch (TYPEOF(x)){
  case LGLSXP: {
    const int *p_x = LOGICAL(x);
    safe_memmove(&p_out[0], &p_x[0], n * sizeof(int));
    break;
  }
  case INTSXP: {
    const int *p_x = INTEGER(x);
    OMP_FOR_SIMD
    for (uint_fast64_t i = 0; i < n; ++i) {
      p_out[i] = is_r_na(p_x[i]) ? NA_INTEGER : sign(p_x[i]);
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    OMP_FOR_SIMD
    for (uint_fast64_t i = 0; i < n; ++i) {
      p_out[i] = is_r_na(p_x[i]) ? NA_INTEGER : sign(p_x[i]);
    }
    break;
  }
  }
  return out;
}
