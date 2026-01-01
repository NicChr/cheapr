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
for (R_xlen_t i = 0; i < n; ++i){                                           \
  p_out[i] = FUN(p_x[i]);                                                   \
}

#define CHEAPR_PARALLEL_MATH_LOOP(FUN)                                    \
if (n_threads > 1){                                                       \
  OMP_PARALLEL_FOR_SIMD(n_threads)                                        \
  CHEAPR_MATH_LOOP(FUN);                                                  \
} else {                                                                  \
  OMP_SIMD                                                                \
  CHEAPR_MATH_LOOP(FUN);                                                  \
}

// For fns like log(x, base)
// where length(base) == length(x)

#define CHEAPR_MATH_LOOP2(FUN)                                                                     \
for (R_xlen_t i = 0; i < n; ++i){                                                                  \
  p_out[i] = FUN(p_x[i], p_y[i]);                                                                  \
}

#define CHEAPR_PARALLEL_MATH_LOOP2(FUN)                                             \
if (n_threads > 1){                                                                 \
  OMP_PARALLEL_FOR_SIMD(n_threads)                                                  \
  CHEAPR_MATH_LOOP2(FUN)                                                            \
} else {                                                                            \
  OMP_SIMD                                                                          \
  CHEAPR_MATH_LOOP2(FUN)                                                            \
}


// For fns like log(x, base) and round(x, digits)
// where length(base) == 1

#define CHEAPR_MATH_LOOP3(FUN)                                                                     \
for (R_xlen_t i = 0; i < n; ++i){                                                                  \
  p_out[i] = FUN(p_x[i], y);                                                                       \
}

#define CHEAPR_PARALLEL_MATH_LOOP3(FUN)                                             \
if (n_threads > 1){                                                                 \
  OMP_PARALLEL_FOR_SIMD(n_threads)                                                  \
  CHEAPR_MATH_LOOP3(FUN)                                                            \
} else {                                                                            \
  OMP_SIMD                                                                          \
  CHEAPR_MATH_LOOP3(FUN)                                                            \
}

// General parallelised optimised macro for applying
// cheapr::math binary fns across 2 vectors
#define CHEAPR_VECTORISED_MATH_LOOP(FUN, x, y, p_x, p_y, p_out, _n, n_threads)                                      \
R_xlen_t _xn = vec::length(x), _yn = vec::length(y);                                                                \
if (_xn == 0 || _yn == 0){                                                                                          \
  _n = 0;                                                                                                           \
}                                                                                                                   \
R_xlen_t xi = 0, yi = 0;                                                                                            \
if (_xn == _n && _yn == _n){                                                                                        \
  if (n_threads > 1){                                                                                               \
    OMP_PARALLEL_FOR_SIMD(n_threads)                                                                                \
    for (R_xlen_t i = 0; i < _n; ++i){                                                                              \
      p_out[i] = FUN(p_x[i], p_y[i]);                                                                               \
    }                                                                                                               \
  } else {                                                                                                          \
    OMP_SIMD                                                                                                        \
    for (R_xlen_t i = 0; i < _n; ++i){                                                                              \
      p_out[i] = FUN(p_x[i], p_y[i]);                                                                               \
    }                                                                                                               \
  }                                                                                                                 \
} else if (_xn == 1 && _yn == _n){                                                                                  \
  auto left = p_x[0];                                                                                               \
  if (n_threads > 1){                                                                                               \
    OMP_PARALLEL_FOR_SIMD(n_threads)                                                                                \
    for (R_xlen_t i = 0; i < _n; ++i){                                                                              \
      p_out[i] = FUN(left, p_y[i]);                                                                                 \
    }                                                                                                               \
  } else {                                                                                                          \
    OMP_SIMD                                                                                                        \
    for (R_xlen_t i = 0; i < _n; ++i){                                                                              \
      p_out[i] = FUN(left, p_y[i]);                                                                                 \
    }                                                                                                               \
  }                                                                                                                 \
} else if (_yn == 1){                                                                                               \
  auto right = p_y[0];                                                                                              \
  if (n_threads > 1){                                                                                               \
    OMP_PARALLEL_FOR_SIMD(n_threads)                                                                                \
    for (R_xlen_t i = 0; i < _n; ++i){                                                                              \
      p_out[i] = FUN(p_x[i], right);                                                                                \
    }                                                                                                               \
  } else {                                                                                                          \
    OMP_SIMD                                                                                                        \
    for (R_xlen_t i = 0; i < _n; ++i){                                                                              \
      p_out[i] = FUN(p_x[i], right);                                                                                \
    }                                                                                                               \
  }                                                                                                                 \
} else if (_xn == _n){                                                                                              \
  for (R_xlen_t i = 0; i < _n; recycle_index(yi, _yn), ++i){                                                        \
    p_out[i] = FUN(p_x[i], p_y[yi]);                                                                                \
  }                                                                                                                 \
} else if (_yn == _n){                                                                                              \
  for (R_xlen_t i = 0; i < _n; recycle_index(xi, _xn), ++i){                                                        \
    p_out[i] = FUN(p_x[xi], p_y[i]);                                                                                \
  }                                                                                                                 \
} else {                                                                                                            \
  for (R_xlen_t i = 0; i < _n; recycle_index(xi, _xn), recycle_index(yi, _yn), ++i){                                \
    p_out[i] = FUN(p_x[xi], p_y[yi]);                                                                               \
  }                                                                                                                 \
}


// Convert integer vector to plain double vector

SEXP convert_int_to_real(SEXP x){
  const int *p_x = integer_ptr_ro(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = SHIELD(new_vector<double>(n));
  double* RESTRICT p_out = real_ptr(out);
  for (R_xlen_t i = 0; i < n; ++i){
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
  int n_threads = calc_threads(n);
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = integer_ptr(out);
    int *p_x = p_out;
    CHEAPR_PARALLEL_MATH_LOOP(r_abs)
      break;
  }
  default: {
    double *p_out = real_ptr(out);
    double *p_x = p_out;
    CHEAPR_PARALLEL_MATH_LOOP(r_abs)
      break;
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_abs(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(new_vector<int>(n));
    int *p_x = integer_ptr(x);
    int* RESTRICT p_out = integer_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_abs)
      break;
  }
  default: {
    out = SHIELD(new_vector<double>(n));
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
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
  int n_threads = calc_threads(n);
  if (Rf_isReal(out)){
    double *p_out = real_ptr(out);
    double *p_x = p_out;
    CHEAPR_PARALLEL_MATH_LOOP(r_floor)
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_floor(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(x);
    break;
  }
  default: {
    out = SHIELD(new_vector<double>(n));
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_floor)
      break;
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
  int n_threads = calc_threads(n);
  if (Rf_isReal(out)){
    double *p_out = real_ptr(out);
    double *p_x = p_out;
    CHEAPR_PARALLEL_MATH_LOOP(r_ceiling)
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_ceiling(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(x);
    break;
  }
  default: {
    out = SHIELD(new_vector<double>(n));
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_ceiling)
      break;
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
  int n_threads = calc_threads(n);
  if (Rf_isReal(out)){
    double *p_out = real_ptr(out);
    double *p_x = p_out;
    CHEAPR_PARALLEL_MATH_LOOP(r_trunc)
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_trunc(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(x);
    break;
  }
  default: {
    out = SHIELD(new_vector<double>(n));
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_trunc)
      break;
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
  int n_threads = calc_threads(n);
  switch (TYPEOF(out)){
  case INTSXP: {
    int *p_out = integer_ptr(out);
    int *p_x = p_out;
    CHEAPR_PARALLEL_MATH_LOOP(r_negate)
      break;
  }
  case REALSXP: {
    double *p_out = real_ptr(out);
    double *p_x = p_out;
    CHEAPR_PARALLEL_MATH_LOOP(r_negate)
      break;
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_negate(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(new_vector<int>(n));
    const int *p_x = integer_ptr_ro(x);
    int* RESTRICT p_out = integer_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_negate)
      break;
  }
  default: {
    out = SHIELD(new_vector<double>(n));
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
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
  int n_threads = calc_threads(n);

  SEXP out = r_null;

  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = real_ptr(out);
  double *p_x = p_out;
  CHEAPR_PARALLEL_MATH_LOOP(r_exp)
    YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_exp(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out = SHIELD(new_vector<double>(n));
  switch (TYPEOF(x)){
  case INTSXP: {
    const int *p_x = integer_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_exp)
      break;
  }
  default: {
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_exp)
      break;
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_set_sqrt(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);

  SEXP out = r_null;

  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_out = real_ptr(out);
  double *p_x = p_out;
  CHEAPR_PARALLEL_MATH_LOOP(r_sqrt)
    YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_sqrt(SEXP x){
  check_numeric(x);
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out = SHIELD(new_vector<double>(n));
  switch (TYPEOF(x)){
  case INTSXP: {
    const int *p_x = integer_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_sqrt)
      break;
  }
  default: {
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_PARALLEL_MATH_LOOP(r_sqrt)
      break;
  }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_log(SEXP x, SEXP base){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(base);
  SHIELD(base = coerce_vec(base, REALSXP)); ++NP;
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t basen = Rf_xlength(base);
  int n_threads = calc_threads(n);
  bool is_base_scalar = basen == 1;
  double *p_y = real_ptr(base);
  double y;
  if (is_base_scalar){
    y = p_y[0];
  }
  bool use_natural_log = is_base_scalar && y == std::exp(1);
  bool use_base10 = is_base_scalar && y == 10;
  bool recycle = !is_base_scalar;

  if (recycle){
    SEXP recycled_objs = SHIELD(make_list(x, base)); ++NP;
    n = length_common(recycled_objs);
    recycle_in_place(recycled_objs, n);
    x = VECTOR_ELT(recycled_objs, 0);
    base = VECTOR_ELT(recycled_objs, 1);
    p_y = real_ptr(base);
  }

  SEXP out = SHIELD(new_vector<double>(n)); ++NP;
  switch (TYPEOF(x)){
  case INTSXP: {
    const int *p_x = integer_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    if (use_natural_log){
      CHEAPR_PARALLEL_MATH_LOOP(r_log)
    } else if (use_base10){
      CHEAPR_PARALLEL_MATH_LOOP(r_log10)
    } else if (is_base_scalar){
      CHEAPR_PARALLEL_MATH_LOOP3(r_log)
    } else {
      CHEAPR_PARALLEL_MATH_LOOP2(r_log)
    }
    break;
  }
  default: {
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    if (use_natural_log){
      CHEAPR_PARALLEL_MATH_LOOP(r_log)
    } else if (use_base10){
      CHEAPR_PARALLEL_MATH_LOOP(r_log10)
    } else if (is_base_scalar){
      CHEAPR_PARALLEL_MATH_LOOP3(r_log)
    } else {
      CHEAPR_PARALLEL_MATH_LOOP2(r_log)
    }
    break;
  }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_round(SEXP x, SEXP digits){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(digits);
  SHIELD(digits = coerce_vec(digits, REALSXP)); ++NP;
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t digitsn = Rf_xlength(digits);
  int n_threads = calc_threads(n);
  bool is_digits_scalar = digitsn == 1;
  double *p_y = real_ptr(digits);
  double y;
  if (is_digits_scalar){
    y = p_y[0];
  }
  bool recycle = !is_digits_scalar;

  if (recycle){
    SEXP recycled_objs = SHIELD(make_list(x, digits)); ++NP;
    n = length_common(recycled_objs);
    recycle_in_place(recycled_objs, n);
    x = VECTOR_ELT(recycled_objs, 0);
    digits = VECTOR_ELT(recycled_objs, 1);
    p_y = real_ptr(digits);
  }

  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = x;
    break;
  }
  default: {
    out = SHIELD(new_vector<double>(n)); ++NP;
    const double *p_x = real_ptr_ro(x);
    double* RESTRICT p_out = real_ptr(out);
    if (is_digits_scalar){
      // If digits == 0 we call r_round(x) instead of r_round(x, digits)
      // Which has a more efficient method
      if (y == 0){
        CHEAPR_PARALLEL_MATH_LOOP(r_round)
      } else {
        CHEAPR_PARALLEL_MATH_LOOP3(r_round)
      }
    } else {
      CHEAPR_PARALLEL_MATH_LOOP2(r_round)
    }
    break;
  }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_signif(SEXP x, SEXP digits){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(digits);
  SHIELD(digits = coerce_vec(digits, REALSXP)); ++NP;
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t digitsn = Rf_xlength(digits);
  int n_threads = calc_threads(n);
  bool is_digits_scalar = digitsn == 1;

  if (!is_digits_scalar){
    SEXP recycled_objs = SHIELD(make_list(x, digits)); ++NP;
    n = length_common(recycled_objs);
    recycle_in_place(recycled_objs, n);
    x = VECTOR_ELT(recycled_objs, 0);
    digits = VECTOR_ELT(recycled_objs, 1);
  }

  double *p_digits = real_ptr(digits);

  SEXP out = SHIELD(new_vector<double>(n)); ++NP;
  double* RESTRICT p_out = real_ptr(out);

  switch (TYPEOF(x)){
  case INTSXP: {
    const int *p_x = integer_ptr_ro(x);
    CHEAPR_VECTORISED_MATH_LOOP(r_signif, x, digits, p_x, p_digits, p_out, n, n_threads)
      break;
  }
  default: {
    const double *p_x = real_ptr_ro(x);
    CHEAPR_VECTORISED_MATH_LOOP(r_signif, x, digits, p_x, p_digits, p_out, n, n_threads)
      break;
  }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_add(SEXP x, SEXP y){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(y);

  SEXP objs = SHIELD(make_list(x, y)); ++NP;
  SHIELD(objs = cpp_cast_common(objs)); ++NP;
  x = VECTOR_ELT(objs, 0);
  y = VECTOR_ELT(objs, 1);

  R_xlen_t n = length_common(objs);

  int n_threads = calc_threads(n);

  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(new_vector<int>(n)); ++NP;
    const int *p_x = integer_ptr_ro(x);
    const int *p_y = integer_ptr_ro(y);
    int* RESTRICT p_out = integer_ptr(out);
    CHEAPR_VECTORISED_MATH_LOOP(r_add, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  default: {
    out = SHIELD(new_vector<double>(n)); ++NP;
    const double *p_x = real_ptr_ro(x);
    const double *p_y = real_ptr_ro(y);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_VECTORISED_MATH_LOOP(r_add, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  }
  YIELD(NP);
  return out;
}
[[cpp11::register]]
SEXP cpp_subtract(SEXP x, SEXP y){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(y);

  SEXP objs = SHIELD(make_list(x, y)); ++NP;
  SHIELD(objs = cpp_cast_common(objs)); ++NP;
  x = VECTOR_ELT(objs, 0);
  y = VECTOR_ELT(objs, 1);

  R_xlen_t n = length_common(objs);

  int n_threads = calc_threads(n);

  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(new_vector<int>(n)); ++NP;
    const int *p_x = integer_ptr_ro(x);
    const int *p_y = integer_ptr_ro(y);
    int* RESTRICT p_out = integer_ptr(out);
    CHEAPR_VECTORISED_MATH_LOOP(r_subtract, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  default: {
    out = SHIELD(new_vector<double>(n)); ++NP;
    const double *p_x = real_ptr_ro(x);
    const double *p_y = real_ptr_ro(y);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_VECTORISED_MATH_LOOP(r_subtract, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  }
  YIELD(NP);
  return out;
}
[[cpp11::register]]
SEXP cpp_multiply(SEXP x, SEXP y){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(y);

  SEXP objs = SHIELD(make_list(x, y)); ++NP;
  SHIELD(objs = cpp_cast_common(objs)); ++NP;
  x = VECTOR_ELT(objs, 0);
  y = VECTOR_ELT(objs, 1);

  R_xlen_t n = length_common(objs);

  int n_threads = calc_threads(n);

  SEXP out;
  switch (TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(new_vector<int>(n)); ++NP;
    const int *p_x = integer_ptr_ro(x);
    const int *p_y = integer_ptr_ro(y);
    int* RESTRICT p_out = integer_ptr(out);
    CHEAPR_VECTORISED_MATH_LOOP(r_multiply, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  default: {
    out = SHIELD(new_vector<double>(n)); ++NP;
    const double *p_x = real_ptr_ro(x);
    const double *p_y = real_ptr_ro(y);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_VECTORISED_MATH_LOOP(r_multiply, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  }
  YIELD(NP);
  return out;
}
[[cpp11::register]]
SEXP cpp_divide(SEXP x, SEXP y){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(y);

  SEXP objs = SHIELD(make_list(x, y)); ++NP;
  SHIELD(objs = cpp_cast_common(objs)); ++NP;
  x = VECTOR_ELT(objs, 0);
  y = VECTOR_ELT(objs, 1);

  R_xlen_t n = length_common(objs);

  int n_threads = calc_threads(n);

  SEXP out = SHIELD(new_vector<double>(n)); ++NP;
  double* RESTRICT p_out = real_ptr(out);

  switch (TYPEOF(x)){
  case INTSXP: {
    const int *p_x = integer_ptr_ro(x);
    const int *p_y = integer_ptr_ro(y);
    CHEAPR_VECTORISED_MATH_LOOP(r_divide, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  default: {
    const double *p_x = real_ptr_ro(x);
    const double *p_y = real_ptr_ro(y);
    CHEAPR_VECTORISED_MATH_LOOP(r_divide, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  }
  YIELD(NP);
  return out;
}
[[cpp11::register]]
SEXP cpp_pow(SEXP x, SEXP y){
  int32_t NP = 0;
  check_numeric(x);
  check_numeric(y);

  SEXP objs = SHIELD(make_list(x, y)); ++NP;
  SHIELD(objs = cpp_cast_common(objs)); ++NP;
  x = VECTOR_ELT(objs, 0);
  y = VECTOR_ELT(objs, 1);

  R_xlen_t n = length_common(objs);

  int n_threads = calc_threads(n);

  SEXP out = SHIELD(new_vector<double>(n)); ++NP;
  double* RESTRICT p_out = real_ptr(out);

  switch (TYPEOF(x)){
  case INTSXP: {
    const int *p_x = integer_ptr_ro(x);
    const int *p_y = integer_ptr_ro(y);
    CHEAPR_VECTORISED_MATH_LOOP(r_pow, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  default: {
    const double *p_x = real_ptr_ro(x);
    const double *p_y = real_ptr_ro(y);
    CHEAPR_VECTORISED_MATH_LOOP(r_pow, x, y, p_x, p_y, p_out, n, n_threads)
      break;
  }
  }
  YIELD(NP);
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
      p_x[i] = r_pow(r_cast<double>(p_x[i]), r_cast<double>(p_y[yi]));
    }
    break;
  }
  default: {
    double *p_x = real_ptr(out);
    double *p_y = real_ptr(y);
    for (R_xlen_t i = 0; i < xn; recycle_index(yi, yn), ++i) {
      p_x[i] = r_pow(p_x[i], p_y[yi]);
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
  R_xlen_t basei = 0;
  if (xn > 0){
    if (basen > xn){
      Rf_error("length(base) must be <= length(x)");
    }
    if (basen == 0){
      Rf_error("length(base) must be be non-zero");
    }
  }
  SEXP out;
  if (!Rf_isReal(x)){
    copy_warning();
    out = SHIELD(convert_int_to_real(x));
  } else {
    out = SHIELD(check_transform_altrep(x));
  }
  double *p_x = real_ptr(out);
  double *p_base = real_ptr(base);
  for (R_xlen_t i = 0; i < xn; recycle_index(basei, basen), ++i) {
    p_x[i] = r_log(p_x[i], p_base[basei]);
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
        p_x[i] = r_round(p_x[i], p_digits[digitsi]);
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
  R_xlen_t n = Rf_xlength(x);
  int n_threads = calc_threads(n);
  SEXP out = SHIELD(vec::new_vector<int>(n));
  int* RESTRICT p_out = integer_ptr(out);
  switch (TYPEOF(x)){
  case LGLSXP: {
    const int *p_x = integer_ptr(x);
    std::copy_n(p_x, n, p_out);
    break;
  }
  case INTSXP: {
    const int *p_x = integer_ptr(x);
    CHEAPR_PARALLEL_MATH_LOOP(r_sign)
      break;
  }
  default: {
    double *p_x = real_ptr(x);
    CHEAPR_PARALLEL_MATH_LOOP(r_sign)
      break;
  }
  }
  YIELD(1);
  return out;
}
