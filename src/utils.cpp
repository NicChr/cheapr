#include <vector>
#include "cheapr.h"
#include <R.h> // R_Calloc

// Miscellaneous functions
// Author: Nick Christofides

static SEXP CHEAPR_CORES = NULL;

[[cpp11::register]]
SEXP cpp_is_simple_atomic_vec(SEXP x){
  return scalar_lgl(is_simple_atomic_vec(x));
}

[[cpp11::register]]
SEXP cpp_is_simple_vec(SEXP x){
  return scalar_lgl(is_simple_vec(x));
}

int int_div(int x, int y){
  return x / y;
}

SEXP xlen_to_r(R_xlen_t x){
  return x > integer_max_ ? Rf_ScalarReal(x) : Rf_ScalarInteger(x);
}

R_xlen_t vec_length(SEXP x){
  if (!Rf_isObject(x) || Rf_isVectorAtomic(x)){
    return Rf_xlength(x);
  } else if (is_df(x)){
    return df_nrow(x);
    // Is x a list?
  } else if (TYPEOF(x) == VECSXP){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return Rf_length(x) > 0 ? vec_length(VECTOR_ELT(x, 0)) : 0;
    } else if (Rf_inherits(x, "POSIXlt")){
      const SEXP *p_x = VECTOR_PTR_RO(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
      // return Rf_xlength(VECTOR_ELT(x, 0));
    } else if (Rf_isObject(x)){
      return r_length(x);
    } else {
      return Rf_xlength(x);
    }
    // Catch-all
  } else {
    return r_length(x);
  }
}

[[cpp11::register]]
SEXP cpp_vector_length(SEXP x){
  return xlen_to_r(vec_length(x));
}

int num_cores(){
  if (CHEAPR_CORES == NULL){
    CHEAPR_CORES = install_utf8("cheapr.cores");
  }
  int n_cores = Rf_asInteger(Rf_GetOption1(CHEAPR_CORES));
  return n_cores >= 1 ? n_cores : 1;
}

SEXP r_address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return make_utf8_char(buf);
}

[[cpp11::register]]
SEXP cpp_address(SEXP x){
  return Rf_ScalarString(r_address(x));
}

// Copy atomic elements from source to target

[[cpp11::register]]
void cpp_set_copy_elements(SEXP source, SEXP target){
  if (TYPEOF(source) != TYPEOF(target)){
    Rf_error("`typeof(target)` must match `typeof(source)`");
  }
  R_xlen_t n = Rf_xlength(source);

  if (n != Rf_xlength(target)){
    Rf_error("target and source must have the same length");
  }

  switch (TYPEOF(source)){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_source = INTEGER(source);
    int *p_target = INTEGER(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(int));
    break;
  }
  case REALSXP: {
    double *p_source = REAL(source);
    double *p_target = REAL(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(double));
    break;
  }
  case STRSXP: {
    const SEXP *p_source = STRING_PTR_RO(source);
    for (R_xlen_t i = 0; i < n; ++i){
      SET_STRING_ELT(target, i, p_source[i]);
    }
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_source = COMPLEX(source);
    Rcomplex *p_target = COMPLEX(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(Rcomplex));
    break;
  }
  case RAWSXP: {
    Rbyte *p_source = RAW(source);
    Rbyte *p_target = RAW(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(Rbyte));
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(source)));
  }
  }
}

// Deep copy of data and attributes

[[cpp11::register]]
SEXP r_copy(SEXP x){
  return Rf_duplicate(x);
}


// Shallow copy of data and attributes

[[cpp11::register]]
SEXP cpp_shallow_copy(SEXP x){
  return Rf_shallow_duplicate(x);
}


// Full copy of data and shallow copy of attributes

[[cpp11::register]]
SEXP cpp_semi_copy(SEXP x){

  // If no attributes then we can just full-copy immediately

  if (is_null(ATTRIB(x))){
    return Rf_duplicate(x);
  }

  bool altrep = ALTREP(x);

  if (!altrep && TYPEOF(x) == VECSXP){

    // Lists

    R_xlen_t n = Rf_xlength(x);
    SEXP out = SHIELD(new_vec(VECSXP, n));
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, Rf_duplicate(p_x[i]));
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    YIELD(1);
    return out;
  } else if (!altrep && cpp_is_simple_atomic_vec(x)){

    // Atomic vectors

    SEXP out = SHIELD(new_vec(TYPEOF(x), Rf_xlength(x)));
    cpp_set_copy_elements(x, out);
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    YIELD(1);
    return out;
  } else {

    // All other R objects
    // This method sometimes full copies twice
    // So I don't use it for non-ALTREP atomic vectors

    SEXP out = SHIELD(Rf_shallow_duplicate(x));
    clear_attributes(out);
    SHIELD(out = Rf_duplicate(out));
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    YIELD(2);
    return out;
  }
}

double cpp_sum(SEXP x){

  R_xlen_t n = Rf_xlength(x);

  double sum = 0;

  switch (CHEAPR_TYPEOF(x)){

  case LGLSXP:
  case INTSXP: {

    const int *p_x = INTEGER(x);

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      sum = is_na_dbl(sum) || is_na_int(p_x[i]) ? NA_REAL : sum + p_x[i];
    }
    break;
  }
  case CHEAPR_INT64SXP: {

    const int64_t *p_x = INTEGER64_PTR(x);

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      sum = is_na_dbl(sum) || is_na_int64(p_x[i]) ? NA_REAL : sum + p_x[i];
    }
    break;
  }
  default: {

    const double *p_x = REAL(x);

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) sum += p_x[i];
    break;
  }
  }
  return sum;
}

double cpp_min(SEXP x){

  R_xlen_t n = Rf_xlength(x);
  switch (CHEAPR_TYPEOF(x)){

  case LGLSXP:
  case INTSXP: {

    if (n == 0) return R_PosInf;

    int *p_x = INTEGER(x);
    int out = integer_max_;

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      out = is_na_int(out) || is_na_int(p_x[i]) ? NA_INTEGER : std::min(out, p_x[i]);
    }
    return out == NA_INTEGER ? NA_REAL : out;
  }
  case CHEAPR_INT64SXP: {

    if (n == 0) return R_PosInf;

    int64_t *p_x = INTEGER64_PTR(x);
    int64_t out = integer64_max_;

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      out = is_na_int64(out) || is_na_int64(p_x[i]) ? NA_INTEGER64 : std::min(out, p_x[i]);
    }
    return out == NA_INTEGER64 ? NA_REAL : out;
  }
  default: {

    double *p_x = REAL(x);
    double out = R_PosInf;

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      out = is_na_dbl(out) || is_na_dbl(p_x[i]) ? NA_REAL : std::min(out, p_x[i]);
    }
    return out;
  }
  }
}

// Internal-only function
// Sum of squared-differences
// Used in `cheapr_var()`

[[cpp11::register]]
double var_sum_squared_diff(SEXP x, double mu){
  R_xlen_t n = Rf_xlength(x);

  double out = 0;

  // NA values are always ignored here

  if (!is_na_dbl(mu)){
    switch (TYPEOF(x)){

    case INTSXP: {
      const int *p_x = INTEGER(x);
      for (R_xlen_t i = 0; i < n; ++i){
        if (is_na_int(p_x[i])) continue;
        out += std::pow(p_x[i] - mu, 2);
      }
      break;
    }
    default: {
      const double *p_x = REAL(x);
      for (R_xlen_t i = 0; i < n; ++i){
        if (is_na_dbl(p_x[i])) continue;
        out += std::pow(p_x[i] - mu, 2);
      }
      break;
    }
    }
  } else {
    out = NA_REAL;
  }
  return out;
}

// Credits to R authors
// Re-purposed .bincode
// The main difference is that codes or breaks can be returned efficiently
// Values outside the (right or left) intervals can be included too

#define CHEAPR_BIN_CODES(_IS_NA_, _NA_VAL_)                                                    \
for (R_xlen_t i = 0; i < n; ++i) {                                                             \
  p_out[i] = _NA_VAL_;                                                                         \
  if (!_IS_NA_(p_x[i])) {                                                                      \
    lo = 0;                                                                                    \
    hi = nb1;                                                                                  \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) || \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                     \
      p_out[i] = (left ? hi : lo) + 1;                                                         \
    }                                                                                          \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                         \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                             \
      while (hi - lo >= 2) {                                                                   \
        cutpoint = (hi + lo)/2;                                                                \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                       \
          lo = cutpoint;                                                                       \
        else                                                                                   \
          hi = cutpoint;                                                                       \
      }                                                                                        \
      p_out[i] = lo + 1 + (right && include_oob);                                              \
    }                                                                                          \
  }                                                                                            \
}                                                                                              \

#define CHEAPR_BIN_NCODES(_IS_NA_, _NA_VAL_)                                                                                       \
for (R_xlen_t i = 0; i < n; ++i) {                                                                                                 \
  p_out[i] = _NA_VAL_;                                                                                                             \
  if (!_IS_NA_(p_x[i])) {                                                                                                          \
    lo = 0;                                                                                                                        \
    hi = nb1;                                                                                                                      \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) ||                                     \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                                                         \
      p_out[i] = p_b[(left ? hi : lo)];                                                                                            \
    }                                                                                                                              \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                                                             \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                                                                 \
      while (hi - lo >= 2) {                                                                                                       \
        cutpoint = (hi + lo)/2;                                                                                                    \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                                                           \
          lo = cutpoint;                                                                                                           \
        else                                                                                                                       \
          hi = cutpoint;                                                                                                           \
      }                                                                                                                            \
      p_out[i] = p_b[lo + (right && include_oob)];                                                                                 \
    }                                                                                                                              \
  }                                                                                                                                \
}                                                                                                                                  \

[[cpp11::register]]
SEXP cpp_bin(SEXP x, SEXP breaks, bool codes, bool right,
             bool include_lowest,
             bool include_oob){
  int n = Rf_length(x);
  int lo, hi, cutpoint;
  int nb = Rf_length(breaks);
  int nb1 = nb - 1;
  bool left = !right;
  bool include_border = include_lowest;
  switch(TYPEOF(x)){
  case INTSXP: {
    if (codes){
    SEXP out = SHIELD(new_vec(INTSXP, n));
    SHIELD(breaks = coerce_vec(breaks, REALSXP));
    const int *p_x = INTEGER(x);
    const double *p_b = REAL(breaks);
    int* RESTRICT p_out = INTEGER(out);
    CHEAPR_BIN_CODES(is_na_int, NA_INTEGER);
    YIELD(2);
    return out;
  } else {
    SEXP out = SHIELD(cpp_semi_copy(x));
    SHIELD(breaks = coerce_vec(breaks, REALSXP));
    const int *p_x = INTEGER(x);
    const double *p_b = REAL(breaks);
    int* RESTRICT p_out = INTEGER(out);
    CHEAPR_BIN_NCODES(is_na_int, NA_INTEGER);
    YIELD(2);
    return out;
  }
  }
  default: {
    if (codes){
    SEXP out = SHIELD(new_vec(INTSXP, n));
    SHIELD(breaks = coerce_vec(breaks, REALSXP));
    const double *p_x = REAL(x);
    const double *p_b = REAL(breaks);
    int* RESTRICT p_out = INTEGER(out);
    CHEAPR_BIN_CODES(is_na_dbl, NA_INTEGER);
    YIELD(2);
    return out;
  } else {
    SEXP out = SHIELD(cpp_semi_copy(x));
    SHIELD(breaks = coerce_vec(breaks, REALSXP));
    const double *p_x = REAL(x);
    const double *p_b = REAL(breaks);
    double* RESTRICT p_out = REAL(out);
    CHEAPR_BIN_NCODES(is_na_dbl, NA_REAL);
    YIELD(2);
    return out;
  }
  }
  }
}

[[cpp11::register]]
SEXP cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na){
  int32_t NP = 0; // count num protections
  if (TYPEOF(condition) != LGLSXP){
    Rf_error("condition must be a logical vector");
  }
  if (TYPEOF(yes) != TYPEOF(no)){
    Rf_error("`typeof(yes)` must match `typeof(no)`");
  }
  if (TYPEOF(yes) != TYPEOF(na)){
    Rf_error("`typeof(yes)` must match `typeof(na)`");
  }
  R_xlen_t n = Rf_xlength(condition);
  R_xlen_t yes_size = Rf_xlength(yes);
  R_xlen_t no_size = Rf_xlength(no);
  R_xlen_t na_size = Rf_xlength(na);

  if (yes_size != 1 && yes_size != n){
    Rf_error("`length(yes)` must be 1 or `length(condition)`");
  }
  if (no_size != 1 && no_size != n){
    Rf_error("`length(no)` must be 1 or `length(condition)`");
  }
  if (na_size != 1 && na_size != n){
    Rf_error("`length(na)` must be 1 or `length(condition)`");
  }

  bool yes_scalar = yes_size == 1;
  bool no_scalar = no_size == 1;
  bool na_scalar = na_size == 1;

  const int *p_x = LOGICAL(condition);
  SEXP out = SHIELD(new_vec(TYPEOF(yes), n)); ++NP;

  switch (TYPEOF(yes)){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int* RESTRICT p_out = INTEGER(out);
    const int *p_yes = INTEGER(yes);
    const int *p_no = INTEGER(no);
    const int *p_na = INTEGER(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      p_out[i] = p_yes[yes_scalar ? 0 : i];
      break;
    }
      case false: {
        p_out[i] = p_no[no_scalar ? 0 : i];
        break;
      }
      default: {
        p_out[i] = p_na[na_scalar ? 0 : i];
        break;
      }
      }
    }
    break;
  }
  case REALSXP: {
    double* RESTRICT p_out = REAL(out);
    const double *p_yes = REAL(yes);
    const double *p_no = REAL(no);
    const double *p_na = REAL(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      p_out[i] = p_yes[yes_scalar ? 0 : i];
      break;
    }
      case false: {
        p_out[i] = p_no[no_scalar ? 0 : i];
        break;
      }
      default: {
        p_out[i] = p_na[na_scalar ? 0 : i];
        break;
      }
      }
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_yes = STRING_PTR_RO(yes);
    const SEXP *p_no = STRING_PTR_RO(no);
    const SEXP *p_na = STRING_PTR_RO(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      SET_STRING_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case false: {
        SET_STRING_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_STRING_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  case CPLXSXP: {
    const Rcomplex *p_yes = COMPLEX(yes);
    const Rcomplex *p_no = COMPLEX(no);
    const Rcomplex *p_na = COMPLEX(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      SET_COMPLEX_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case false: {
        SET_COMPLEX_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_COMPLEX_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  case RAWSXP: {
    const Rbyte *p_yes = RAW(yes);
    const Rbyte *p_no = RAW(no);
    const Rbyte *p_na = RAW(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      SET_RAW_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case false: {
        SET_RAW_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_RAW_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  case VECSXP: {
    const SEXP *p_yes = VECTOR_PTR_RO(yes);
    const SEXP *p_no = VECTOR_PTR_RO(no);
    const SEXP *p_na = VECTOR_PTR_RO(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      SET_VECTOR_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case false: {
        SET_VECTOR_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_VECTOR_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  default: {
    YIELD(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(yes)));
  }
  }
  YIELD(NP);
  return out;
}

// Counts number of true, false and NAs in a logical vector in one pass

[[cpp11::register]]
SEXP cpp_lgl_count(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  const int *p_x = LOGICAL(x);

  R_xlen_t i;
  R_xlen_t ntrue = 0, nfalse = 0;

  if (n_cores > 1){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:ntrue, nfalse)
    for (i = 0; i < n; ++i){
      ntrue += p_x[i] == TRUE;
      nfalse += p_x[i] == FALSE;
    }
  } else {
    OMP_FOR_SIMD
    for (i = 0; i < n; ++i){
      ntrue += p_x[i] == TRUE;
      nfalse += p_x[i] == FALSE;
    }
  }
  R_xlen_t nna = n - ntrue - nfalse;

  SEXP out = SHIELD(new_vec(n > integer_max_ ? REALSXP : INTSXP, 3));
  SEXP names = SHIELD(new_vec(STRSXP, 3));
  SET_STRING_ELT(names, 0, make_utf8_char("true"));
  SET_STRING_ELT(names, 1, make_utf8_char("false"));
  SET_STRING_ELT(names, 2, make_utf8_char("na"));

  if (n > integer_max_){
    SET_REAL_ELT(out, 0, static_cast<double>(ntrue));
    SET_REAL_ELT(out, 1, static_cast<double>(nfalse));
    SET_REAL_ELT(out, 2, static_cast<double>(nna));
  } else {
    SET_INTEGER_ELT(out, 0, static_cast<int>(ntrue));
    SET_INTEGER_ELT(out, 1, static_cast<int>(nfalse));
    SET_INTEGER_ELT(out, 2, static_cast<int>(nna));
  }

  set_names(out, names);

  YIELD(2);
  return out;
}

// Essentially `x | y` but updates x by reference

[[cpp11::register]]
SEXP cpp_set_or(SEXP x, SEXP y){
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);

  R_xlen_t n = xn == 0 || yn == 0 ? 0 : xn;

  R_xlen_t i, yi;

  int *p_x = LOGICAL(x);
  const int *p_y = LOGICAL(y);


  for (i = yi = 0; i < n; yi = (++yi == yn) ? 0 : yi, ++i){

    if (p_x[i] != TRUE){
      if (p_y[yi] == TRUE){
        p_x[i] = TRUE;
      } else if ((p_x[i] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
        p_x[i] = NA_LOGICAL;
      } else if (p_x[i] == TRUE || p_y[yi] == TRUE){
        p_x[i] = TRUE;
      }
    }
  }
  return x;
}

// coerceVector() that accounts for int64

SEXP coerce_vector(SEXP source, SEXPTYPE type){
  if (type == CHEAPR_INT64SXP){
    SEXP temp = SHIELD(coerce_vec(source, REALSXP));
    SEXP out = SHIELD(cpp_numeric_to_int64(temp));
    YIELD(2);
    return out;
  } else if (is_int64(source)){
    SEXP temp = SHIELD(cpp_int64_to_numeric(source));
    SEXP out = SHIELD(coerce_vec(temp, type));
    YIELD(2);
    return out;
  } else {
    return coerce_vec(source, type);
  }
}

// Basic growth rate
// i.e the expected percent change per unit of time
// eqn: (new / old ) ^ 1/n
//
// Which is not the same as: new / old
// of which the latter is mistakenly and commonly referred to
// as a 'growth rate'

// o(1) time and space

double growth_rate(double a, double b, double n){
  if (n < 1.0){
    Rf_error("n must be >= 1");
  }
  if (n == R_PosInf){
    Rf_error("n must be finite positive");
  }
  if (n == 1.0){
    return NA_REAL;
  }
  if (a == 0.0 && b == 0.0){
    return 1.0;
  } else {
    return std::pow(b / a, 1.0 / (n - 1.0));
  }
}

#define R_SCALAR_AS_DOUBLE(x, na_val) ((double) (x == na_val ? NA_REAL : x))

[[cpp11::register]]
SEXP cpp_growth_rate(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  if (n == 0){
    return new_vec(REALSXP, 0);
  }
  if (n == 1){
    return Rf_ScalarReal(NA_REAL);
  }
  double a, b;
  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    int x_n = INTEGER(x)[n - 1];
    int x_1 = INTEGER(x)[0];
    a = R_SCALAR_AS_DOUBLE(x_1, NA_INTEGER);
    b = R_SCALAR_AS_DOUBLE(x_n, NA_INTEGER);
    break;
  }
  case CHEAPR_INT64SXP: {
    int64_t x_n = INTEGER64_PTR(x)[n - 1];
    int64_t x_1 = INTEGER64_PTR(x)[0];
    a = R_SCALAR_AS_DOUBLE(x_1, NA_INTEGER64);
    b = R_SCALAR_AS_DOUBLE(x_n, NA_INTEGER64);
    break;
  }
  case REALSXP: {
    b = REAL(x)[n - 1];
    a = REAL(x)[0];
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  return Rf_ScalarReal(growth_rate(a, b, n));
}

SEXP create_df_row_names(int n){
  if (n > 0){
    SEXP out = SHIELD(new_vec(INTSXP, 2));
    INTEGER(out)[0] = NA_INTEGER;
    INTEGER(out)[1] = -n;
    YIELD(1);
    return out;
  } else {
    return new_vec(INTSXP, 0);
  }
}

[[cpp11::register]]
SEXP cpp_name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep){
  int32_t NP = 0;
  if (is_null(names)) return names;

  if (TYPEOF(names) != STRSXP){
    Rf_error("`names` must be a character vector of names in %s", __func__);
  }
  if (TYPEOF(dup_sep) != STRSXP || Rf_length(dup_sep) != 1){
    Rf_error("`dup_sep` must be a character vector of length 1 in %s", __func__);
  }
  if (TYPEOF(empty_sep) != STRSXP || Rf_length(empty_sep) != 1){
    Rf_error("`empty_sep` must be a character vector of length 1 in %s", __func__);
  }
  int n = Rf_length(names);
  SEXP is_dup = SHIELD(Rf_duplicated(names, FALSE)); ++NP;
  SEXP is_dup_from_last = SHIELD(Rf_duplicated(names, TRUE)); ++NP;
  cpp_set_or(is_dup, is_dup_from_last);

  SEXP r_true = SHIELD(scalar_lgl(true)); ++NP;
  SEXP dup_locs = SHIELD(cpp_which_val(is_dup, r_true, false)); ++NP;

  int n_dups = Rf_length(dup_locs);

  SEXP out = SHIELD(new_vec(STRSXP, n)); ++NP;
  cpp_set_copy_elements(names, out);

  SEXP temp, replace;

  if (n_dups > 0){
    temp = SHIELD(sset_vec(names, dup_locs, true)); ++NP;
    replace = SHIELD(base_paste0(temp, dup_sep, dup_locs)); ++NP;
    cpp_loc_set_replace(out, dup_locs, replace);
  }

  SEXP is_empty = SHIELD(new_vec(LGLSXP, n)); ++NP;
  int *p_is_empty = LOGICAL(is_empty);
  bool empty;
  int n_empty = 0;

  for (int i = 0; i < n; ++i){
    empty = (STRING_ELT(names, i) == R_BlankString);
    n_empty += empty;
    p_is_empty[i] = empty;
  }

  SEXP r_n_empty = SHIELD(Rf_ScalarInteger(n_empty)); ++NP;

  if (n_empty > 0){
    SEXP empty_locs = SHIELD(cpp_val_find(is_empty, r_true, false, r_n_empty)); ++NP;
    temp = SHIELD(sset_vec(names, empty_locs, true)); ++NP;
    replace = SHIELD(base_paste0(temp, empty_sep, empty_locs)); ++NP;
    cpp_loc_set_replace(out, empty_locs, replace);
  }

  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_rebuild(SEXP target, SEXP source, SEXP target_attr_names,
                 SEXP source_attr_names, bool shallow_copy){

  int32_t NP = 0;

  if (shallow_copy){
    SHIELD(target = Rf_shallow_duplicate(target)); ++NP;
  }

  if (address_equal(target, source)){
    YIELD(NP);
    return target;
  }

  SEXP target_attrs = ATTRIB(target);
  SEXP source_attrs = ATTRIB(source);

  // Start from clean slate - no attributes
  clear_attributes(target);

  SEXP tag = R_NilValue;
  SEXP current = R_NilValue;

  const SEXP *p_ta = STRING_PTR_RO(target_attr_names);
  const SEXP *p_sa = STRING_PTR_RO(source_attr_names);

  const int n_target = Rf_length(target_attr_names);
  const int n_source = Rf_length(source_attr_names);

  for (int i = 0; i < n_target; ++i){
    current = target_attrs;
    while (!is_null(current)){

      tag = TAG(current);

      if (PRINTNAME(tag) == p_ta[i]){
        Rf_setAttrib(target, tag, CAR(current));
        break;
      }
      current = CDR(current);
    }
  }

  for (int i = 0; i < n_source; ++i){
    current = source_attrs;
    while (!is_null(current)){

      tag = TAG(current);

      if (PRINTNAME(tag) == p_sa[i]){
        Rf_setAttrib(target, tag, CAR(current));
        break;
      }
      current = CDR(current);
    }
  }

  YIELD(NP);
  return target;
}

// Coalesce a list of string vectors

[[cpp11::register]]
SEXP cpp_str_coalesce(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of character vectors in %s", __func__);
  }

  int32_t NP = 0;
  uint_fast64_t n = Rf_xlength(x);
  uint_fast64_t out_size = 0;
  uint_fast64_t m;

  const SEXP *p_x = VECTOR_PTR_RO(x);
  std::vector<const SEXP*> str_ptrs(n);

  SEXP char_vec = R_NilValue;
  uint32_t xtype;

  bool shallow_duplicated = false;

  for (uint_fast64_t i = 0; i < n; ++i){
    char_vec = p_x[i];
    xtype = TYPEOF(char_vec);

    if (xtype != STRSXP){
      if (!shallow_duplicated){
        SHIELD(x = Rf_shallow_duplicate(x)); ++NP;
        p_x = VECTOR_PTR_RO(x);
        shallow_duplicated = true;
      }
      SET_VECTOR_ELT(x, i, base_as_character(char_vec));
      char_vec = p_x[i];
    }

    str_ptrs[i] = STRING_PTR_RO(char_vec);

    if (xtype != NILSXP){
      m = Rf_xlength(char_vec);
      if (m == 0){
        YIELD(NP);
        return Rf_allocVector(STRSXP, 0);
      }
      out_size = std::max(out_size, m);
    }
  }

  SEXP out = SHIELD(Rf_allocVector(STRSXP, out_size)); ++NP;

  SEXP inner_char = R_BlankString;

  uint_fast64_t n_nas;

  for (uint_fast64_t i = 0; i < out_size; ++i){
    n_nas = 0;
    for (uint_fast64_t j = 0; j < n; ++j){
      m = Rf_xlength(p_x[j]);
      if (m == 0) continue;
      inner_char = str_ptrs[j][i % m];
      n_nas += inner_char == NA_STRING;
      if (!(inner_char == R_BlankString || inner_char == NA_STRING)){
        SET_STRING_ELT(out, i, inner_char);
        break;
      }
      // If all ith elements are NA, then return NA
      if (n_nas == n){
        SET_STRING_ELT(out, i, NA_STRING);
      }
    }
  }
  YIELD(NP);
  return out;
}

// Keep this for now, may use in the future

// #include <cpp11.hpp>
// #include "cheapr.h"
// #include <R.h>
// #include <Rinternals.h>
// #include <functional>
// #include <unordered_map>
//
// struct sexp_metadata {
//   void* getter;                                 // Direct pointer to data (e.g., INTEGER(x))
//   std::function<void(R_xlen_t, void*)> setter;    // Setter function (kept for flexibility)
//   void* ptr_to_x;                                 // Pointer to the SEXP object x
//   int type_id;
// };
//
// // Helper function to create metadata for a specific SEXPTYPE (compile-time)
// template<int SEXPTYPE>
// sexp_metadata create_metadata(SEXP x) {
//   sexp_metadata metadata;
//   metadata.ptr_to_x = static_cast<void*>(x);      // Store pointer to x
//   metadata.type_id = SEXPTYPE;
//
//   if constexpr (SEXPTYPE == INTSXP) {
//     void* data = static_cast<void*>(INTEGER(x));  // Single pointer for integer data
//     metadata.getter = data;
//     metadata.setter = [data](R_xlen_t i, void* value) {
//       static_cast<int*>(data)[i] = *static_cast<int*>(value);
//     };
//   } else if constexpr (SEXPTYPE == REALSXP) {
//     void* data = static_cast<void*>(REAL(x));  // Single pointer for integer data
//     metadata.getter = data;
//     metadata.setter = [data](R_xlen_t i, void* value) {
//       static_cast<double*>(data)[i] = *static_cast<double*>(value);
//     };
//   } else if constexpr (SEXPTYPE == LGLSXP) {
//     void* data = static_cast<void*>(LOGICAL(x));  // Single pointer for integer data
//     metadata.getter = data;
//     metadata.setter = [data](R_xlen_t i, void* value) {
//       static_cast<int*>(data)[i] = *static_cast<int*>(value);
//     };
//   } else if constexpr (SEXPTYPE == STRSXP) {
//     void* data = const_cast<void*>(reinterpret_cast<const void*>(STRING_PTR_RO(x)));
//     metadata.getter = data;
//     metadata.setter = [x](R_xlen_t i, void* value) {
//       SET_STRING_ELT(x, i, static_cast<SEXP>(value));
//     };
//   } else {
//     Rf_error("Unknown SEXP type");
//   }
//
//   return metadata;
// }
//
// // Function to get metadata based on runtime TYPEOF(x)
// sexp_metadata get_sexp_metadata(SEXP x) {
//   int type = TYPEOF(x);
//   switch (type) {
//   case INTSXP:
//     return create_metadata<INTSXP>(x);
//   case REALSXP:
//     return create_metadata<REALSXP>(x);
//   case LGLSXP:
//     return create_metadata<LGLSXP>(x);
//   case STRSXP:
//     return create_metadata<STRSXP>(x);
//   default:
//     Rf_error("Unsupported SEXP type: %s", Rf_type2char(TYPEOF(x)));
//   }
// }
// R_xlen_t foo(SEXP x){
//   sexp_metadata metadata = get_sexp_metadata(x);
//
//   R_xlen_t count = 0;
//
//   auto* data = static_cast<int*>(metadata.getter);
//
//   OMP_FOR_SIMD
//   for (R_xlen_t i = 0; i < Rf_xlength(x); ++i){
//     count += is_na_int(data[i]);
//   }
//   return count;
// }
// R_xlen_t bar(SEXP x){
//   sexp_metadata metadata = get_sexp_metadata(x);
//
//   R_xlen_t count = 0;
//
//   // auto* data = static_cast<int*>(metadata.getter);
//   //
//   // int val = 1;
//   //
//   // OMP_FOR_SIMD
//   // for (R_xlen_t i = 0; i < Rf_xlength(x); ++i){
//   //   metadata.setter(i, static_cast<void*>(&val));
//   // }
//
//   auto* data = static_cast<SEXP*>(metadata.getter);
//
//   SEXP val = Rf_protect(Rf_mkChar("ok"));
//
//   for (R_xlen_t i = 0; i < Rf_xlength(x); ++i){
//     metadata.setter(i, static_cast<void*>(val));
//   }
//   YIELD(1);
//   return count;
// }
// R_xlen_t foobar(SEXP x){
//   sexp_metadata metadata = get_sexp_metadata(x);
//
//   R_xlen_t count = 0;
//
//   auto* data = static_cast<SEXP*>(metadata.getter);
//
//   SEXP val = Rf_protect(Rf_mkChar("ok"));
//
//   for (R_xlen_t i = 0; i < Rf_xlength(x); ++i){
//     SET_STRING_ELT(x, i, val);
//   }
//   YIELD(1);
//   return count;
// }

// R's internal tabulate with faster unsigned int check
[[cpp11::register]]
SEXP cpp_tabulate(SEXP x, uint32_t n_bins){

  if (n_bins > integer_max_){
    Rf_error("`n_bins` must be < 2^31 in %s", __func__);
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = SHIELD(new_vec(INTSXP, n_bins));
  const int *p_x = INTEGER_RO(x);
  int* RESTRICT p_out = INTEGER(out);

  // Initialise counts to 0
  memset(p_out, 0, n_bins * sizeof(int));

  uint32_t one = 1;

  OMP_FOR_SIMD
  for (R_xlen_t i = 0 ; i < n; ++i){
    if ((static_cast<uint32_t>(p_x[i]) - one) < n_bins){
      ++p_out[p_x[i] - 1];
    }
  }
  YIELD(1);
  return out;
}
