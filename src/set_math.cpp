#include "cheapr.h"

// Basic math operations by reference
// All NA and NaN values are ignored

void check_numeric(SEXP x){
  if (!(Rf_isNumeric(x) && !is_object(x))){
    Rf_error("x must be a numeric vector");
  }
}

void copy_warning(){
  Rf_warning("x is not a double vector and has been copied, it will not be replaced by reference.\n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_log(x)`");
}

SEXP check_transform_altrep(SEXP x){
  if (altrep::is_altrep(x)){
    Rf_warning("Cannot update an ALTREP by reference, a copy has been made. \n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- set_abs(x)`");
    return altrep_materialise(x);
  } else {
    return x;
  }
}

#define CHEAPR_MATH_LOOP(FUN)                                               \
for (R_xlen_t i = 0; i < n; ++i) {                                          \
  p_out[i] = FUN(p_out[i]);                                                 \
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
  int *p_x = integer_ptr(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = SHIELD(new_double(n));
  double* RESTRICT p_out = real_ptr(out);
  for (int i = 0; i < n; ++i){
    p_out[i] = r_cast<double>(p_x[i]);
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_abs(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(out);
  int n_cores = get_cores(n);
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = integer_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_abs)
    break;
  }
  default: {
    double *p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_abs)
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
  int n_cores = get_cores(n);
  if (Rf_isReal(out)){
    double *p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_floor)
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_ceiling(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(x);
  int n_cores = get_cores(n);
  if (Rf_isReal(out)){
    double *p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_ceiling)
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_trunc(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(out);
  int n_cores = get_cores(n);
  if (Rf_isReal(out)){
    double *p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_trunc)
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_change_sign(SEXP x){
  check_numeric(x);
  SEXP out = SHIELD(check_transform_altrep(x));
  R_xlen_t n = Rf_xlength(out);
  int n_cores = get_cores(n);
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = integer_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_negate)
      break;
  }
  case REALSXP: {
    double *p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_negate)
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
  int n_cores = get_cores(n);

  SEXP out = r_null;

  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = real_ptr(out);
  CHEAPR_PARALLEL_MATH_LOOP(std::exp)
    YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_sqrt(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_cores = get_cores(n);

  SEXP out = r_null;

  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = real_ptr(out);
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
    int *p_x = integer_ptr(out);
    int *p_y = integer_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_add(p_x[i], p_y[yi]);
    }
    break;
  }
  default: {
    copy_warning();
    SHIELD(out = vec::coerce_vec(out, REALSXP)); ++NP;
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_add(p_x[i], p_y[yi]);
    }
    break;
  }
  }
    break;
  }
  default: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    double *p_x = real_ptr(out);
    int *p_y = integer_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_add(p_x[i], r_cast<double>(p_y[yi]));
    }
    break;
  }
  default: {
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_add(p_x[i], p_y[yi]);
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
    int *p_x = integer_ptr(out);
    int *p_y = integer_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_subtract(p_x[i], p_y[yi]);
    }
    break;
  }
  default: {
    copy_warning();
    SHIELD(out = vec::coerce_vec(out, REALSXP)); ++NP;
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_subtract(p_x[i], p_y[yi]);
    }
    break;
  }
  }
    break;
  }
  default: {
    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    double *p_x = real_ptr(out);
    int *p_y = integer_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_subtract(p_x[i], r_cast<double>(p_y[yi]));
    }
    break;
  }
  default: {
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_subtract(p_x[i], p_y[yi]);
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
    int *p_x = integer_ptr(out);
    int *p_y = integer_ptr(y);

    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_multiply(p_x[i], p_y[yi]);
    }
    break;
  }
  default: {
    copy_warning();
    SHIELD(out = vec::coerce_vec(out, REALSXP)); ++NP;
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_multiply(p_x[i], p_y[yi]);
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
    double *p_x = real_ptr(out);
    int *p_y = integer_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_multiply(p_x[i], r_cast<double>(p_y[yi]));
    }
    break;
  }
  default: {
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_multiply(p_x[i], p_y[yi]);
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

  SEXP out = r_null;

  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {
    double *p_x = real_ptr(out);
    int *p_y = integer_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_divide(p_x[i], r_cast<double>(p_y[yi]));
    }
    break;
  }
  default: {
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (uint_fast64_t i = 0; i < xn; recycle_index(yi, yn), ++i){
      p_x[i] = r_divide(p_x[i], p_y[yi]);
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
  R_xlen_t yi = 0;
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
    double *p_x = real_ptr(out);
    int *p_y = integer_ptr(y);
    for (R_xlen_t i = 0; i < xn; recycle_index(yi, yn), ++i) {
      p_x[i] = std::pow(r_cast<double>(p_x[i]), r_cast<double>(p_y[yi]));
    }
    break;
  }
  default: {
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (R_xlen_t i = 0; i < xn; recycle_index(yi, yn), ++i) {
      p_x[i] = std::pow(p_x[i], p_y[yi]);
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
  double *p_x = real_ptr(out);
  double *p_base = real_ptr(base);
  if (n_cores > 1){
    OMP_PARALLEL_FOR_SIMD
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t basei = i % basen;
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_base[basei])) ?
      na::real : std::log(p_x[i]) / std::log(p_base[basei]);
    }
  } else {
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < xn; ++i) {
      R_xlen_t basei = i % basen;
      p_x[i] = (is_r_na(p_x[i]) || is_r_na(p_base[basei])) ?
      na::real : std::log(p_x[i]) / std::log(p_base[basei]);
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
      double *p_x = real_ptr(out);
      const int *p_digits = integer_ptr(digits);
      for (uint_fast64_t i = 0; i < xn;
      recycle_index(digitsi, digitsn), ++i) {
        p_x[i] = r_round(p_x[i], r_cast<double>(p_digits[digitsi]));
      }
      break;
    }
    default: {
      double *p_x = real_ptr(out);
      const double *p_digits = real_ptr(digits);
      for (uint_fast64_t i = 0; i < xn;
      recycle_index(digitsi, digitsn), ++i) {
        p_x[i] = r_round(p_x[i], p_digits[digitsi]);
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
  SEXP out = SHIELD(vec::new_integer(n));
  int* RESTRICT p_out = integer_ptr(out);
  switch (TYPEOF(x)){
  case LGLSXP: {
    const int *p_x = integer_ptr(x);
    fast_copy_n(p_x, n, p_out);
    break;
  }
  case INTSXP: {
    const int *p_x = integer_ptr(x);
    OMP_FOR_SIMD
    for (uint_fast64_t i = 0; i < n; ++i) {
      p_out[i] = r_sign(p_x[i]);
    }
    break;
  }
  case REALSXP: {
    double *p_x = real_ptr(x);
    OMP_FOR_SIMD
    for (uint_fast64_t i = 0; i < n; ++i) {
      p_out[i] = r_sign(p_x[i]);
    }
    break;
  }
  }
  YIELD(1);
  return out;
}
