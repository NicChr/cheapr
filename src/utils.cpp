#include "cheapr.h"
#include <R.h> // R_Calloc

// Miscellaneous functions
// Author: Nick Christofides

[[cpp11::register]]
int cpp_max_threads(){
  return OMP_MAX_THREADS;
}

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
  return as_vector(vec::length(x));
}

[[cpp11::register]]
SEXP cpp_address(SEXP x){
  return as_vector(address(x));
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
    int *p_source = integer_ptr(source);
    int *p_target = integer_ptr(target);
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(int));
    break;
  }
  case REALSXP: {
    double *p_source = real_ptr(source);
    double *p_target = real_ptr(target);
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(double));
    break;
  }
  case STRSXP: {
    const r_string_t *p_source = string_ptr_ro(source);
    for (R_xlen_t i = 0; i < n; ++i){
      set_value<r_string_t>(target, i, p_source[i]);
    }
    break;
  }
  case CPLXSXP: {
    r_complex_t *p_source = complex_ptr(source);
    r_complex_t *p_target = complex_ptr(target);
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(r_complex_t));
    break;
  }
  case RAWSXP: {
    r_byte_t *p_source = raw_ptr(source);
    r_byte_t *p_target = raw_ptr(target);
    safe_memmove(&p_target[0], &p_source[0], n * sizeof(r_byte_t));
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
  return vec::deep_copy(x);
}


// Shallow copy of data and attributes

[[cpp11::register]]
SEXP cpp_shallow_copy(SEXP x){
  return vec::shallow_copy(x);
}


// Full copy of data and shallow copy of attributes

[[cpp11::register]]
SEXP cpp_semi_copy(SEXP x){

  // If no attributes then we can just full-copy immediately

  if (!has_attrs(x)){
    return vec::deep_copy(x);
  }

  bool altrep = ALTREP(x);

  if (!altrep && TYPEOF(x) == VECSXP){

    // Lists

    R_xlen_t n = Rf_xlength(x);
    SEXP out = SHIELD(new_list(n));
    const SEXP *p_x = list_ptr_ro(x);
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, vec::deep_copy(p_x[i]));
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    YIELD(1);
    return out;
  } else if (!altrep && cheapr_is_simple_atomic_vec(x)){

    // Atomic vectors

    SEXP out = SHIELD(internal::new_vec(TYPEOF(x), Rf_xlength(x)));
    cpp_set_copy_elements(x, out);
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    YIELD(1);
    return out;
  } else {

    // All other R objects
    // This method sometimes full copies twice
    // So I don't use it for non-ALTREP atomic vectors

    SEXP out = SHIELD(vec::shallow_copy(x));
    attr::clear_attrs(out);
    SHIELD(out = vec::deep_copy(out));
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

    const int *p_x = integer_ptr(x);

    for (R_xlen_t i = 0; i < n; ++i){
      sum = is_r_na(sum) || is_r_na(p_x[i]) ? na::real : sum + p_x[i];
    }
    break;
  }
  case CHEAPR_INT64SXP: {

    const int64_t *p_x = integer64_ptr_ro(x);

    for (R_xlen_t i = 0; i < n; ++i){
      sum = is_r_na(sum) || is_r_na(p_x[i]) ? na::real : sum + p_x[i];
    }
    break;
  }
  default: {

    const double *p_x = real_ptr(x);

    #pragma omp simd reduction(+:sum)
    for (R_xlen_t i = 0; i < n; ++i) sum += p_x[i];
    break;
  }
  }
  return sum;
}

// Extremely fast range(x, na.rm = F)
SEXP cpp_range(SEXP x){

  R_xlen_t n = Rf_xlength(x);

  double lo = r_limits::r_pos_inf;
  double hi = r_limits::r_neg_inf;

  if (n > 0){
    switch (CHEAPR_TYPEOF(x)){

    case LGLSXP:
    case INTSXP: {

      const int *p_x = integer_ptr_ro(x);
      int min = std::numeric_limits<int>::max();
      int max = std::numeric_limits<int>::min();

      OMP_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        min = std::min(min, p_x[i]);
        max = std::max(max, p_x[i]);
      }
      lo = is_r_na(min) ? na::real : min;
      hi = is_r_na(min) ? na::real : max;
      break;
    }
    case CHEAPR_INT64SXP: {

      const int64_t *p_x = integer64_ptr_ro(x);

      int64_t min = std::numeric_limits<int64_t>::max();
      int64_t max = std::numeric_limits<int64_t>::min();

      OMP_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        min = std::min(min, p_x[i]);
        max = std::max(max, p_x[i]);
      }
      lo = is_r_na(min) ? na::real : min;
      hi = is_r_na(min) ? na::real : max;
      break;
    }
    default: {

      const double *p_x = real_ptr_ro(x);

      double min = r_limits::r_pos_inf;
      double max = r_limits::r_neg_inf;

      for (R_xlen_t i = 0; i < n; ++i){
        if (is_r_na(p_x[i])){
          min = na::real;
          max = na::real;
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
  return combine(lo, hi);
}

double cpp_min(SEXP x){
  return real_ptr_ro(cpp_range(x))[0];
}
double cpp_max(SEXP x){
  return real_ptr_ro(cpp_range(x))[1];
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
      const int *p_x = integer_ptr(x);
      for (R_xlen_t i = 0; i < n; ++i){
        if (is_r_na(p_x[i])) continue;
        out += std::pow(p_x[i] - mu, 2);
      }
      break;
    }
    default: {
      const double *p_x = real_ptr(x);
      for (R_xlen_t i = 0; i < n; ++i){
        if (is_r_na(p_x[i])) continue;
        out += std::pow(p_x[i] - mu, 2);
      }
      break;
    }
    }
  } else {
    out = na::real;
  }
  return out;
}

// Credits to R authors
// Re-purposed .bincode
// The main difference is that codes or breaks can be returned efficiently
// Values outside the (right or left) intervals can be included too

#define CHEAPR_BIN_CODES(na)                                                                                  \
for (R_xlen_t i = 0; i < n; ++i) {                                                                        \
  p_out[i] = na;                                               \
  if (!is_r_na(p_x[i])) {                                                                                 \
    lo = 0;                                                                                               \
    hi = nb1;                                                                                             \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) ||            \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                                \
      p_out[i] = (left ? hi : lo) + 1;                                                                    \
    }                                                                                                     \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                                    \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                                        \
      while (hi - lo >= 2) {                                                                              \
        cutpoint = (hi + lo)/2;                                                                           \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                                  \
          lo = cutpoint;                                                                                  \
        else                                                                                              \
          hi = cutpoint;                                                                                  \
      }                                                                                                   \
      p_out[i] = lo + 1 + (right && include_oob);                                                         \
    }                                                                                                     \
  }                                                                                                       \
}

#define CHEAPR_BIN_NCODES(na)                                                                                                              \
for (R_xlen_t i = 0; i < n; ++i) {                                                                                                     \
  p_out[i] = na;                                               \
  if (!is_r_na(p_x[i])) {                                                                                                              \
    lo = 0;                                                                                                                            \
    hi = nb1;                                                                                                                          \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) ||                                         \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                                                             \
      p_out[i] = p_b[(left ? hi : lo)];                                                                                                \
    }                                                                                                                                  \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                                                                 \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                                                                     \
      while (hi - lo >= 2) {                                                                                                           \
        cutpoint = (hi + lo)/2;                                                                                                        \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                                                               \
          lo = cutpoint;                                                                                                               \
        else                                                                                                                           \
          hi = cutpoint;                                                                                                               \
      }                                                                                                                                \
      p_out[i] = p_b[lo + (right && include_oob)];                                                                                     \
    }                                                                                                                                  \
  }                                                                                                                                    \
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
    SEXP out = SHIELD(vec::new_vector<int>(n));
    SHIELD(breaks = vec::coerce_vec(breaks, REALSXP));
    const int *p_x = integer_ptr(x);
    const double *p_b = real_ptr(breaks);
    int* RESTRICT p_out = vector_ptr<int>(out);
    CHEAPR_BIN_CODES(na_value<int>())
    YIELD(2);
    return out;
  } else {
    SEXP out = SHIELD(cpp_semi_copy(x));
    SHIELD(breaks = vec::coerce_vec(breaks, REALSXP));
    const int *p_x = integer_ptr(x);
    const double *p_b = real_ptr(breaks);
    int* RESTRICT p_out = integer_ptr(out);
    CHEAPR_BIN_NCODES(na_value<int>())
    YIELD(2);
    return out;
  }
  }
  default: {
    if (codes){
    SEXP out = SHIELD(vec::new_vector<int>(n));
    SHIELD(breaks = vec::coerce_vec(breaks, REALSXP));
    const double *p_x = real_ptr(x);
    const double *p_b = real_ptr(breaks);
    int* RESTRICT p_out = integer_ptr(out);
    CHEAPR_BIN_CODES(na_value<int>())
    YIELD(2);
    return out;
  } else {
    SEXP out = SHIELD(cpp_semi_copy(x));
    SHIELD(breaks = vec::coerce_vec(breaks, REALSXP));
    const double *p_x = real_ptr(x);
    const double *p_b = real_ptr(breaks);
    double* RESTRICT p_out = real_ptr(out);
    CHEAPR_BIN_NCODES(na_value<double>())
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
  int n_threads = calc_threads(n);

  const r_bool_t *p_x = logical_ptr_ro(x);

  R_xlen_t i;
  R_xlen_t ntrue = 0, nfalse = 0;

  if (n_threads > 1){
#pragma omp parallel for simd num_threads(n_threads) reduction(+:ntrue, nfalse)
    for (i = 0; i < n; ++i){
      ntrue += is_r_true(p_x[i]);
      nfalse += is_r_false(p_x[i]);
    }
  } else {
#pragma omp simd reduction(+:ntrue, nfalse)
    for (i = 0; i < n; ++i){
      ntrue += is_r_true(p_x[i]);
      nfalse += is_r_false(p_x[i]);
    }
  }

  R_xlen_t nna = n - ntrue - nfalse;

  SEXP out = SHIELD(combine(
    arg("true") = ntrue,
    arg("false") = nfalse,
    arg("na") = nna
  ));

  YIELD(1);
  return out;
}

// Essentially `x | y` but updates x by reference

[[cpp11::register]]
SEXP cpp_set_or(SEXP x, SEXP y){

  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);

  R_xlen_t n = xn == 0 || yn == 0 ? 0 : xn;

  R_xlen_t i, yi;

  r_bool_t *p_x = logical_ptr(x);
  const r_bool_t *p_y = logical_ptr_ro(y);

  for (i = yi = 0; i < n; yi = (++yi == yn) ? 0 : yi, ++i){

    if (!is_r_true(p_x[i])){
      if (is_r_true(p_y[yi])){
        p_x[i] = r_true;
      } else if (is_r_na(p_x[i]) || is_r_na(p_y[yi])){
        p_x[i] = na::logical;
      } else if (is_r_true(p_x[i]) || is_r_true(p_y[yi])){
        p_x[i] = r_true;
      }
    }
  }
  return x;
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
  if (n == r_limits::r_pos_inf){
    Rf_error("n must be finite positive");
  }
  if (n == 1.0){
    return na::real;
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
    return new_vector<double>(0);
  }
  if (n == 1){
    return as_vector(na::real);
  }
  double a, b;
  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    int x_n = integer_ptr(x)[n - 1];
    int x_1 = integer_ptr(x)[0];
    a = r_cast<double>(x_1);
    b = r_cast<double>(x_n);
    break;
  }
  case CHEAPR_INT64SXP: {
    int64_t x_n = integer64_ptr(x)[n - 1];
    int64_t x_1 = integer64_ptr(x)[0];
    a = r_cast<double>(x_1);
    b = r_cast<double>(x_n);
    break;
  }
  case REALSXP: {
    b = real_ptr(x)[n - 1];
    a = real_ptr(x)[0];
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  return as_vector(growth_rate(a, b, n));
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

  SEXP scalar_true = SHIELD(as_vector(r_true)); ++NP;
  SEXP dup_locs = SHIELD(cpp_which_val(is_dup, scalar_true, false)); ++NP;

  int n_dups = Rf_length(dup_locs);

  SEXP out = SHIELD(new_vector<r_string_t>(n)); ++NP;
  cpp_set_copy_elements(names, out);

  SEXP temp = r_null;
  SEXP replace = r_null;

  if (n_dups > 0){
    temp = SHIELD(sset_vec(names, dup_locs, true)); ++NP;
    replace = SHIELD(r_paste(R_BlankScalarString, r_null, temp, dup_sep, dup_locs)); ++NP;
    replace_in_place(out, dup_locs, replace, false);
  }

  SEXP is_empty = SHIELD(new_vector<r_bool_t>(n)); ++NP;
  int *p_is_empty = integer_ptr(is_empty);
  bool empty;
  int n_empty = 0;

  for (int i = 0; i < n; ++i){
    empty = (STRING_ELT(names, i) == blank_r_string);
    n_empty += empty;
    p_is_empty[i] = empty;
  }

  SEXP r_n_empty = SHIELD(as_vector(n_empty)); ++NP;

  if (n_empty > 0){
    SEXP empty_locs = SHIELD(cpp_val_find(is_empty, scalar_true, false, r_n_empty)); ++NP;
    temp = SHIELD(sset_vec(names, empty_locs, true)); ++NP;
    replace = SHIELD(r_paste(R_BlankScalarString, r_null, temp, empty_sep, empty_locs)); ++NP;
    replace_in_place(out, empty_locs, replace, false);
  }

  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_rebuild(SEXP target, SEXP source, SEXP target_attr_names,
                 SEXP source_attr_names, bool shallow_copy){

  int32_t NP = 0;

  if (shallow_copy){
    SHIELD(target = cheapr::vec::shallow_copy(target)); ++NP;
  } else if (address_equal(target, source)){
    YIELD(NP);
    return target;
  }


  SEXP target_attrs = SHIELD(get_attrs(target)); ++NP;
  SEXP source_attrs = SHIELD(get_attrs(source)); ++NP;

  SEXP target_nms = SHIELD(attr::get_old_names(target_attrs)); ++NP;
  SEXP source_nms = SHIELD(attr::get_old_names(source_attrs)); ++NP;


  const r_string_t *p_target_nms = string_ptr_ro(target_nms);
  const r_string_t *p_source_nms = string_ptr_ro(source_nms);

  int n_target_nms = Rf_length(target_nms);
  int n_source_nms = Rf_length(source_nms);

  // Start from clean slate - no attributes
  attr::clear_attrs(target);

  r_string_t curr_tag;

  const r_string_t *p_ta = string_ptr_ro(target_attr_names);
  const r_string_t *p_sa = string_ptr_ro(source_attr_names);

  const int n_target = Rf_length(target_attr_names);
  const int n_source = Rf_length(source_attr_names);


  // Source attributes to keep
  for (int i = 0; i < n_source; ++i){
    for (int j = 0; j < n_source_nms; ++j){

      curr_tag = p_source_nms[j];

      if (curr_tag == p_sa[i]){
        set_attr(target, r_cast<r_symbol_t>(curr_tag), VECTOR_ELT(source_attrs, j));
        break;
      }
    }
  }

  // Target attributes to keep
  for (int i = 0; i < n_target; ++i){
    for (int j = 0; j < n_target_nms; ++j){

      curr_tag = p_target_nms[j];

      if (curr_tag == p_ta[i]){
        set_attr(target, r_cast<r_symbol_t>(curr_tag), VECTOR_ELT(target_attrs, j));
        break;
      }
    }
  }

  YIELD(NP);
  return target;
}


// R's internal tabulate with faster unsigned int check
[[cpp11::register]]
SEXP cpp_tabulate(SEXP x, uint32_t n_bins){

  if (n_bins > r_limits::r_int_max){
    Rf_error("`n_bins` must be < 2^31 in %s", __func__);
  }
  R_xlen_t n = Rf_xlength(x);

  // Initialise counts to 0
  SEXP out = SHIELD(vec::new_vector<int>(n_bins, 0));
  const int *p_x = integer_ptr_ro(x);
  int* RESTRICT p_out = integer_ptr(out);

  uint32_t one = 1;

  OMP_SIMD
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
  return as_vector(vec::all_whole_numbers(x, tol_, na_rm_));
}

SEXP match(SEXP y, SEXP x, int no_match){
  if (Rf_xlength(x) < 100000 && Rf_xlength(y) < 100000){
    return Rf_match(y, x, no_match);
  } else {
    return eval_pkg_fun("fast_match", "cheapr", env::base_env, x, y, no_match);
  }
}

SEXP get_vec_names(SEXP x){
  if (vec::is_atomic(x)){
    return get_old_names(x);
  } else {
    switch(get_r_type(x)){
    case R_null:
    case R_df: {
      return r_null;
    }
    case R_list: {
      return get_old_names(x);
    }
    case R_unk: {
      return eval_pkg_fun("names", "base", env::base_env, x);
    }
    }
    return r_null;
  }
}

void set_vec_names(SEXP x, SEXP names){
  if (is_null(names)){
    return;
  } else if (vec::is_atomic(x)){
    attr::set_old_names(x, names);
    return;
  } else {
    switch(get_r_type(x)){
    case R_null:
    case R_df: {
      return;
    }
    case R_list: {
      attr::set_old_names(x, names);
      return;
    }
    default: {
      SEXP vec_with_names = SHIELD(eval_pkg_fun("names<-", "base", env::base_env, x, names));
      SEXP attrs = SHIELD(get_attrs(vec_with_names));
      attr::set_attrs(x, attrs);
      YIELD(2);
      return;
    }
    }
  }
}

bool vec_has_names(SEXP x){
  return !is_null(get_vec_names(x));
}

[[cpp11::register]]
SEXP cheapr_do_memory_leak_test(){

  // To be run in valgrind

  // Check that 4000 bytes are not lost

  std::vector<int32_t> ints(1000);
  SEXP r_ints = r_safe(Rf_protect)(r_safe(new_vector<int>)(ints.size()));
  SEXP seq = r_safe(Rf_protect)(r_safe(cpp_seq_len)(ints.size()));
  SEXP repl = r_safe(Rf_protect)(r_safe(as_vector)(-1));
  r_safe(replace_in_place)(r_ints, seq, repl, true);
  r_safe(Rf_unprotect)(3);
  r_safe(Rf_error)("%s", "Expected error! This should not cause a C++ memory leak");
  return r_ints; // Never reached
}

// Below will trigger a memory leak, only use for testing purposes
// SEXP cheapr_unsafe_init_memory_leak(){
//   std::vector<int> ints(1000);
//   SEXP r_ints = r_safe(SHIELD)(r_safe(new_vec)(INTSXP, ints.size()));
//   SEXP seq = r_safe(SHIELD)(r_safe(cpp_seq_len)(ints.size()));
//   SEXP repl = r_safe(SHIELD)(r_safe(make_vec)(-1));
//   r_safe(replace_in_place)(r_ints, seq, repl, true);
//   r_safe(YIELD)(3);
//   Rf_error("%s", "Expected error! This will cause a C++ memory leak");
//   return r_ints; // Never reached
// }
