#include "cheapr.h"
#include <R.h> // R_Calloc

// Miscellaneous functions
// Author: Nick Christofides

[[cpp11::register]]
bool cpp_is_simple_atomic_vec(SEXP x){
  return cheapr_is_simple_atomic_vec(x);
}

[[cpp11::register]]
bool cpp_is_simple_vec(SEXP x){
  return cheapr_is_simple_vec(x);
}

[[cpp11::register]]
SEXP cpp_vector_length(SEXP x){
  return as_r_scalar(vector_length(x));
}

[[cpp11::register]]
SEXP cpp_address(SEXP x){
  return as_r_scalar(address(x));
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
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(int));
    break;
  }
  case REALSXP: {
    double *p_source = REAL(source);
    double *p_target = REAL(target);
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(double));
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
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(Rcomplex));
    break;
  }
  case RAWSXP: {
    Rbyte *p_source = RAW(source);
    Rbyte *p_target = RAW(target);
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(Rbyte));
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
    const SEXP *p_x = LIST_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, Rf_duplicate(p_x[i]));
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    YIELD(1);
    return out;
  } else if (!altrep && cheapr_is_simple_atomic_vec(x)){

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
      sum = is_r_na(sum) || is_r_na(p_x[i]) ? NA_REAL : sum + p_x[i];
    }
    break;
  }
  case CHEAPR_INT64SXP: {

    const int64_t *p_x = INTEGER64_PTR_RO(x);

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      sum = is_r_na(sum) || is_r_na(p_x[i]) ? NA_REAL : sum + p_x[i];
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

// Extremely fast range(x, na.rm = F)
SEXP cpp_range(SEXP x){

  R_xlen_t n = Rf_xlength(x);

  SEXP out = SHIELD(new_vec(REALSXP, 2));
  double lo = R_PosInf;
  double hi = R_NegInf;

  if (n > 0){
    switch (CHEAPR_TYPEOF(x)){

    case LGLSXP:
    case INTSXP: {

      const int *p_x = INTEGER_RO(x);
      int min = std::numeric_limits<int>::max();
      int max = std::numeric_limits<int>::min();

      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        min = std::min(min, p_x[i]);
        max = std::max(max, p_x[i]);
      }
      lo = is_r_na(min) ? NA_REAL : min;
      hi = is_r_na(min) ? NA_REAL : max;
      break;
    }
    case CHEAPR_INT64SXP: {

      const int64_t *p_x = INTEGER64_PTR_RO(x);

      int64_t min = std::numeric_limits<int64_t>::max();
      int64_t max = std::numeric_limits<int64_t>::min();

      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        min = std::min(min, p_x[i]);
        max = std::max(max, p_x[i]);
      }
      lo = is_r_na(min) ? NA_REAL : min;
      hi = is_r_na(min) ? NA_REAL : max;
      break;
    }
    default: {

      const double *p_x = REAL_RO(x);

      double min = R_PosInf;
      double max = R_NegInf;

      for (R_xlen_t i = 0; i < n; ++i){
        if (is_r_na(p_x[i])){
          min = NA_REAL;
          max = NA_REAL;
          break;
        }
        min = std::min(min, p_x[i]);
        max = std::max(max, p_x[i]);
      }
      lo = min;
      hi = max;
      break;
    }
    }
  }
  SET_REAL_ELT(out, 0, lo);
  SET_REAL_ELT(out, 1, hi);
  YIELD(1);
  return out;
}

double cpp_min(SEXP x){
  return REAL_RO(cpp_range(x))[0];
}
double cpp_max(SEXP x){
  return REAL_RO(cpp_range(x))[1];
}

// Internal-only function
// Sum of squared-differences
// Used in `cheapr_var()`

[[cpp11::register]]
double var_sum_squared_diff(SEXP x, double mu){
  R_xlen_t n = Rf_xlength(x);

  double out = 0;

  // NA values are always ignored here

  if (!is_r_na(mu)){
    switch (TYPEOF(x)){

    case INTSXP: {
      const int *p_x = INTEGER(x);
      for (R_xlen_t i = 0; i < n; ++i){
        if (is_r_na(p_x[i])) continue;
        out += std::pow(p_x[i] - mu, 2);
      }
      break;
    }
    default: {
      const double *p_x = REAL(x);
      for (R_xlen_t i = 0; i < n; ++i){
        if (is_r_na(p_x[i])) continue;
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

#define CHEAPR_BIN_CODES                                                                                \
for (R_xlen_t i = 0; i < n; ++i) {                                                                      \
  p_out[i] = na_type(p_out[0]);                                                                         \
  if (!is_r_na(p_x[i])) {                                                                                 \
    lo = 0;                                                                                             \
    hi = nb1;                                                                                           \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) ||          \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                              \
      p_out[i] = (left ? hi : lo) + 1;                                                                  \
    }                                                                                                   \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                                  \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                                      \
      while (hi - lo >= 2) {                                                                            \
        cutpoint = (hi + lo)/2;                                                                         \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                                \
          lo = cutpoint;                                                                                \
        else                                                                                            \
          hi = cutpoint;                                                                                \
      }                                                                                                 \
      p_out[i] = lo + 1 + (right && include_oob);                                                       \
    }                                                                                                   \
  }                                                                                                     \
}

#define CHEAPR_BIN_NCODES                                                                                                            \
for (R_xlen_t i = 0; i < n; ++i) {                                                                                                   \
  p_out[i] = na_type(p_out[0]);                                                                                                      \
  if (!is_r_na(p_x[i])) {                                                                                                              \
    lo = 0;                                                                                                                          \
    hi = nb1;                                                                                                                        \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) ||                                       \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                                                           \
      p_out[i] = p_b[(left ? hi : lo)];                                                                                              \
    }                                                                                                                                \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                                                               \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                                                                   \
      while (hi - lo >= 2) {                                                                                                         \
        cutpoint = (hi + lo)/2;                                                                                                      \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                                                             \
          lo = cutpoint;                                                                                                             \
        else                                                                                                                         \
          hi = cutpoint;                                                                                                             \
      }                                                                                                                              \
      p_out[i] = p_b[lo + (right && include_oob)];                                                                                   \
    }                                                                                                                                \
  }                                                                                                                                  \
}

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
    CHEAPR_BIN_CODES;
    YIELD(2);
    return out;
  } else {
    SEXP out = SHIELD(cpp_semi_copy(x));
    SHIELD(breaks = coerce_vec(breaks, REALSXP));
    const int *p_x = INTEGER(x);
    const double *p_b = REAL(breaks);
    int* RESTRICT p_out = INTEGER(out);
    CHEAPR_BIN_NCODES
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
    CHEAPR_BIN_CODES;
    YIELD(2);
    return out;
  } else {
    SEXP out = SHIELD(cpp_semi_copy(x));
    SHIELD(breaks = coerce_vec(breaks, REALSXP));
    const double *p_x = REAL(x);
    const double *p_b = REAL(breaks);
    double* RESTRICT p_out = REAL(out);
    CHEAPR_BIN_NCODES
    YIELD(2);
    return out;
  }
  }
  }
}

// Counts number of true, false and NAs in a logical vector in one pass

[[cpp11::register]]
SEXP cpp_lgl_count(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  const r_boolean *p_x = BOOLEAN_RO(x);

  R_xlen_t i;
  R_xlen_t ntrue = 0, nfalse = 0;

  if (n_cores > 1){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:ntrue, nfalse)
    for (i = 0; i < n; ++i){
      ntrue += p_x[i] == r_true;
      nfalse += p_x[i] == r_false;
    }
  } else {
    OMP_FOR_SIMD
    for (i = 0; i < n; ++i){
      ntrue += p_x[i] == r_true;
      nfalse += p_x[i] == r_false;
    }
  }

  R_xlen_t nna = n - ntrue - nfalse;

  SEXP names = SHIELD(new_r_chars("true", "false", "na"));
  SEXP out;

  if (n > INTEGER_MAX){
    out = SHIELD(cpp11::writable::doubles(
    {
      static_cast<double>(ntrue),
      static_cast<double>(nfalse),
      static_cast<double>(nna)
    }
    ));
  } else {
    out = SHIELD(cpp11::writable::integers(
    {
      static_cast<int>(ntrue),
      static_cast<int>(nfalse),
      static_cast<int>(nna)
    }
    ));
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

    if (p_x[i] != 1){
      if (p_y[yi] == 1){
        p_x[i] = 1;
      } else if ((p_x[i] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
        p_x[i] = NA_LOGICAL;
      } else if (p_x[i] == 1 || p_y[yi] == 1){
        p_x[i] = 1;
      }
    }
  }
  return x;
}

// SEXP cpp_set_and(SEXP x, SEXP y){
//
//   R_xlen_t xn = Rf_xlength(x);
//   R_xlen_t yn = Rf_xlength(y);
//
//   R_xlen_t n = xn == 0 || yn == 0 ? 0 : xn;
//
//   R_xlen_t i, yi;
//
//   int *p_x = LOGICAL(x);
//   const int *p_y = LOGICAL_RO(y);
//
//   for (i = yi = 0; i < n; yi = (++yi == yn) ? 0 : yi, ++i){
//
//     if (p_x[i] != 0){
//       if (p_y[yi] == 0){
//         p_x[i] = 0;
//       } else if ((p_x[i] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
//         p_x[i] = NA_LOGICAL;
//       } else if (p_x[i] == 1 && p_y[yi] == 1){
//         p_x[i] = 1;
//       }
//     }
//   }
//
//   return x;
// }

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

[[cpp11::register]]
SEXP cpp_growth_rate(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  if (n == 0){
    return new_vec(REALSXP, 0);
  }
  if (n == 1){
    return as_r_scalar(NA_REAL);
  }
  double a, b;
  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    int x_n = INTEGER(x)[n - 1];
    int x_1 = INTEGER(x)[0];
    a = as_double(x_1);
    b = as_double(x_n);
    break;
  }
  case CHEAPR_INT64SXP: {
    int64_t x_n = INTEGER64_PTR(x)[n - 1];
    int64_t x_1 = INTEGER64_PTR(x)[0];
    a = as_double(x_1);
    b = as_double(x_n);
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
  return as_r_scalar(growth_rate(a, b, n));
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

  SEXP r_true = SHIELD(as_r_scalar(true)); ++NP;
  SEXP dup_locs = SHIELD(cpp_which_val(is_dup, r_true, false)); ++NP;

  int n_dups = Rf_length(dup_locs);

  SEXP out = SHIELD(new_vec(STRSXP, n)); ++NP;
  cpp_set_copy_elements(names, out);

  SEXP temp = R_NilValue;
  SEXP replace = R_NilValue;

  if (n_dups > 0){
    temp = SHIELD(sset_vec(names, dup_locs, true)); ++NP;
    replace = SHIELD(r_paste(R_BlankScalarString, R_NilValue, temp, dup_sep, dup_locs)); ++NP;
    static_cast<void>(cpp_replace(out, dup_locs, replace, true, false));
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

  SEXP r_n_empty = SHIELD(as_r_scalar(n_empty)); ++NP;

  if (n_empty > 0){
    SEXP empty_locs = SHIELD(cpp_val_find(is_empty, r_true, false, r_n_empty)); ++NP;
    temp = SHIELD(sset_vec(names, empty_locs, true)); ++NP;
    replace = SHIELD(r_paste(R_BlankScalarString, R_NilValue, temp, empty_sep, empty_locs)); ++NP;
    static_cast<void>(cpp_replace(out, empty_locs, replace, true, false));
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


// R's internal tabulate with faster unsigned int check
[[cpp11::register]]
SEXP cpp_tabulate(SEXP x, uint32_t n_bins){

  if (n_bins > INTEGER_MAX){
    Rf_error("`n_bins` must be < 2^31 in %s", __func__);
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = SHIELD(new_vec(INTSXP, n_bins));
  const int *p_x = INTEGER_RO(x);
  int* RESTRICT p_out = INTEGER(out);

  // Initialise counts to 0
  std::fill(p_out, p_out + n_bins, 0);

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

// Returns true if all numbers are whole numbers
// otherwise false
// Returns NA when na_rm is false and the function can't find any
// non-whole numbers and there is at least 1 NA

[[cpp11::register]]
SEXP cpp_is_whole_number(SEXP x, double tol_, bool na_rm_){
  return as_r_scalar(static_cast<int>(vec_is_whole_number(x, tol_, na_rm_)));
}

SEXP match(SEXP y, SEXP x, int no_match){
  if (Rf_xlength(x) < 100000 && Rf_xlength(y) < 100000){
    return Rf_match(y, x, no_match);
  } else {
    return cheapr_fast_match(x, y, no_match);
  }
}

SEXP get_vec_names(SEXP x){
  if (Rf_isVectorAtomic(x)){
    return get_names(x);
  } else {
    switch(get_r_type(x)){
    case r_null:
    case r_df: {
      return R_NilValue;
    }
    case r_list: {
      return get_names(x);
    }
    case r_unk: {
      SEXP r_names_fn = SHIELD(find_pkg_fun("names", "base", false));
      SEXP expr = SHIELD(Rf_lang2(r_names_fn, x));
      SEXP out = SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
      YIELD(3);
      return out;
    }
    }
    return R_NilValue;
  }
}

SEXP set_vec_names(SEXP x, SEXP names){
  if (is_null(names)){
    return x;
  } else if (Rf_isVectorAtomic(x)){
    SEXP out = SHIELD(shallow_copy(x));
    Rf_namesgets(out, names);
    YIELD(1);
    return out;
  } else {
    switch(get_r_type(x)){
    case r_null:
    case r_df: {
      return x;
    }
    case r_list: {
      SEXP out = SHIELD(shallow_copy(x));
      Rf_namesgets(out, names);
      YIELD(1);
      return out;
    }
    case r_unk: {
      SEXP r_names_fn = SHIELD(find_pkg_fun("names<-", "base", false));
      SEXP expr = SHIELD(Rf_lang3(r_names_fn, x, names));
      SEXP out = SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
      YIELD(3);
      return out;
    }
    }
    return R_NilValue;
  }
}

bool vec_has_names(SEXP x){
  return !is_null(get_vec_names(x));
}
