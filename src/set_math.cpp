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

SEXP check_transform_altrep(SEXP x){
  if (is_altrep(x)){
    Rf_warning("Cannot update an ALTREP by reference, a copy has been made. \n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_abs(x)`");
    return altrep_materialise(x);
  } else {
    return x;
  }
}

#define CHEAPR_MATH_LOOP(FUN)                                             \
for (R_xlen_t i = 0; i < n; ++i) {                                        \
  p_out[i] = is_r_na(p_out[i]) ? p_out[i] : FUN(p_out[i]);                  \
}                                                                         \

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
  SEXP out = SHIELD(new_vec(REALSXP, n));
  double* RESTRICT p_out = REAL(out);
  for (int i = 0; i < n; ++i){
    p_out[i] = as_double(p_x[i]);
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
    CHEAPR_PARALLEL_MATH_LOOP(std::abs)
    break;
  }
  default: {
    double *p_out = REAL(out);
    CHEAPR_PARALLEL_MATH_LOOP(std::abs)
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
    CHEAPR_PARALLEL_MATH_LOOP(std::floor)
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
    CHEAPR_PARALLEL_MATH_LOOP(std::ceil)
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
    CHEAPR_PARALLEL_MATH_LOOP(cheapr::trunc)
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
    CHEAPR_PARALLEL_MATH_LOOP(cheapr::negate)
    break;
  }
  case REALSXP: {
    double *p_out = REAL(out);
    CHEAPR_PARALLEL_MATH_LOOP(cheapr::negate)
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

  SEXP out = R_NilValue;

  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = REAL(out);
  CHEAPR_PARALLEL_MATH_LOOP(std::exp)
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_sqrt(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  SEXP out = R_NilValue;

  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = REAL(out);
  CHEAPR_PARALLEL_MATH_LOOP(std::sqrt)
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_add(SEXP x, SEXP y){
  check_numeric(x);
  check_numeric(y);
  int32_t NP = 0;
  SEXP out = SHIELD(check_transform_altrep(x)); ++NP;
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

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
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_INTEGER : p_x[i] + p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    SHIELD(out = coerce_vec(out, REALSXP)); ++NP;
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
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
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_REAL : p_x[i] + p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] + p_y[yi];
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
  int32_t NP = 0;
  SEXP out = SHIELD(check_transform_altrep(x)); ++NP;
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

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
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_INTEGER : p_x[i] - p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    SHIELD(out = coerce_vec(out, REALSXP)); ++NP;
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
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
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_REAL : p_x[i] - p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] - p_y[yi];
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
  int32_t NP = 0;
  SEXP out = SHIELD(check_transform_altrep(x)); ++NP;
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

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

    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (p_x[i] == NA_INTEGER || p_y[yi] == NA_INTEGER) ?
      NA_INTEGER : p_x[i] * p_y[yi];
    }
    break;
  }
  case REALSXP: {
    copy_warning();
    SHIELD(out = coerce_vec(out, REALSXP)); ++NP;
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
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
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ?
      NA_REAL : p_x[i] * p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) ||is_r_na(p_y[yi])) ?
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
  uint_fast64_t xn = Rf_xlength(x);
  uint_fast64_t yn = Rf_xlength(y);
  uint_fast64_t yi = 0;

  if (xn > 0){
    if (yn > xn){
      Rf_error("length(y) must be <= length(x)");
    }
    if (yn == 0){
      Rf_error("length(y) must be be non-zero");
    }
  }

  SEXP out = R_NilValue;

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
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] / p_y[yi];
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(out);
    double *p_y = REAL(y);
    for (uint_fast64_t i = 0; i < xn; yi = (++yi == yn) ? 0 : yi, ++i){
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_y[yi])) ? NA_REAL : p_x[i] / p_y[yi];
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
  uint_fast64_t xn = Rf_xlength(out);
  uint_fast64_t digitsn = Rf_xlength(digits);
  uint_fast64_t digitsi = 0;

  double tempx, mfactor;

  if (xn > 0){
    if (digitsn > xn){
      YIELD(1);
      Rf_error("`length(digits)` must be `<= length(x)`");
    }
    if (digitsn == 0){
      YIELD(1);
      Rf_error("`length(digits)` must be be non-zero");
    }
  }

  // We don't need to round integers.

  if (Rf_isReal(out)){
    switch (TYPEOF(digits)){
    case INTSXP: {
      double *p_x = REAL(out);
      const int *p_digits = INTEGER(digits);
        for (uint_fast64_t i = 0; i < xn; digitsi = (++digitsi == digitsn) ? 0 : digitsi, ++i) {
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
      for (uint_fast64_t i = 0; i < xn; digitsi = (++digitsi == digitsn) ? 0 : digitsi, ++i) {
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
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_int_sign(SEXP x){
  check_numeric(x);
  uint_fast64_t n = Rf_xlength(x);
  SEXP out = SHIELD(new_vec(INTSXP, n));
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
  YIELD(1);
  return out;
}

// Math functions

// I wanted to export these to C but there is an issue with C allocated R
// vectors and the NAMED mechanism of reference counting
// Vectors created through Rf_allocVector are not named when assigned
// to a C SEXP
// This means we can't use the MAYBE_REFERENCED functions safely in this case
// These are also inefficient as R fns as they are doing an unnecessary loop

// SEXP cpp_abs(SEXP x){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_abs(x)); ++NP;
//   YIELD(NP);
//   return out;
// }
// SEXP cpp_floor(SEXP x){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_floor(x)); ++NP;
//   YIELD(NP);
//   return out;
// }
// SEXP cpp_ceiling(SEXP x){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_ceiling(x)); ++NP;
//   YIELD(NP);
//   return out;
// }
// SEXP cpp_trunc(SEXP x){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_trunc(x)); ++NP;
//   YIELD(NP);
//   return out;
// }
// SEXP cpp_invert_sign(SEXP x){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_change_sign(x)); ++NP;
//   YIELD(NP);
//   return out;
// }
// SEXP cpp_exp(SEXP x){
//   int32_t NP = 0;
//   SHIELD(x = cast<r_numeric_t>(x, R_NilValue)); ++NP;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_exp(x)); ++NP;
//   YIELD(NP);
//   return out;
// }
// SEXP cpp_sqrt(SEXP x){
//   int32_t NP = 0;
//   SHIELD(x = cast<r_numeric_t>(x, R_NilValue)); ++NP;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_sqrt(x)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_log10(SEXP x){
//   int32_t NP = 0;
//   SHIELD(x = cast<r_numeric_t>(x, R_NilValue)); ++NP;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP log10 = SHIELD(as_r_scalar(10.0)); ++NP;
//   SEXP out = SHIELD(cpp_set_log(x, log10)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_log(SEXP x, SEXP base){
//   int32_t NP = 0;
//   SHIELD(x = cast<r_numeric_t>(x, R_NilValue)); ++NP;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_log(x, base)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_round(SEXP x, SEXP digits){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_round(x, digits)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_pow(SEXP x, SEXP y){
//   int32_t NP = 0;
//   SHIELD(x = cast<r_numeric_t>(x, R_NilValue)); ++NP;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_pow(x, y)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_add(SEXP x, SEXP y){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_add(x, y)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_subtract(SEXP x, SEXP y){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_subtract(x, y)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_multiply(SEXP x, SEXP y){
//   int32_t NP = 0;
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_multiply(x, y)); ++NP;
//   YIELD(NP);
//   return out;
// }
//
// SEXP cpp_divide(SEXP x, SEXP y){
//   int32_t NP = 0;
//
//   SHIELD(x = cast<r_numeric_t>(x, R_NilValue)); ++NP;
//
//   if (MAYBE_REFERENCED(x)){
//     SHIELD(x = Rf_duplicate(x)); ++NP;
//   }
//   SEXP out = SHIELD(cpp_set_divide(x, y)); ++NP;
//   YIELD(NP);
//   return out;
// }
