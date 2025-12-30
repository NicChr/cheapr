#include "cheapr.h"

// A more memory-efficient which()
// Author: Nick Christofides
// License: MIT

// Count the number of true values

R_xlen_t count_true(const r_bool_t* RESTRICT px, const uint_fast64_t n){
  uint_fast64_t size = 0;
  int n_threads = calc_threads(n);
  if (n_threads != 1){
#pragma omp parallel for simd num_threads(n_threads) reduction(+:size)
    for (uint_fast64_t j = 0; j != n; ++j) size += is_r_true(px[j]);
    return size;
  } else {
#pragma omp simd reduction(+:size)
    for (uint_fast64_t j = 0; j != n; ++j) size += is_r_true(px[j]);
    return size;
  }
}

#define CHEAPR_WHICH_VAL(VAL)                                        \
if (is_r_na(VAL)){                                                   \
  while (whichi < out_size){                                         \
    p_out[whichi] = i + 1;                                           \
    whichi += is_r_na(p_x[i++]);                                     \
  }                                                                  \
} else {                                                             \
  while (whichi < out_size){                                         \
    p_out[whichi] = i + 1;                                           \
    whichi += (p_x[i++] == VAL);                                     \
  }                                                                  \
}


#define CHEAPR_WHICH_VAL_INVERTED(VAL)                             \
if (is_r_na(VAL)){                                                 \
  while (whichi < out_size){                                       \
    p_out[whichi] = i + 1;                                         \
    whichi += !is_r_na(p_x[i++]);                                  \
  }                                                                \
} else {                                                           \
  while (whichi < out_size){                                       \
    p_out[whichi] = i + 1;                                         \
    whichi += (p_x[i++] != VAL);                                   \
  }                                                                \
}


[[cpp11::register]]
SEXP cpp_which_(SEXP x, bool invert){
  R_xlen_t n = Rf_xlength(x);
  const r_bool_t *p_x = logical_ptr_ro(x);
  bool is_long = (n > r_limits::r_int_max);
  if (invert){
    if (is_long){
      R_xlen_t size = count_true(p_x, n);
      R_xlen_t out_size = n - size;
      SEXP out = SHIELD(new_vector<r_double_t>(out_size));
      r_double_t* RESTRICT p_out = real_ptr(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL_INVERTED(r_true);
      YIELD(1);
      return out;
    } else {
      int size = count_true(p_x, n);
      int out_size = n - size;
      SEXP out = SHIELD(vec::new_vector<r_int_t>(out_size));
      r_int_t* RESTRICT p_out = vector_ptr<r_int_t>(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL_INVERTED(r_true);
      YIELD(1);
      return out;
    }
  } else {
    if (is_long){
      R_xlen_t out_size = count_true(p_x, n);
      SEXP out = SHIELD(new_vector<r_double_t>(out_size));
      r_double_t* RESTRICT p_out = real_ptr(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL(r_true);
      YIELD(1);
      return out;
    } else {
      int out_size = count_true(p_x, n);
      SEXP out = SHIELD(vec::new_vector<r_int_t>(out_size));
      r_int_t* RESTRICT p_out = vector_ptr<r_int_t>(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL(r_true);
      YIELD(1);
      return out;
    }
  }
}

SEXP cpp_val_find(SEXP x, SEXP value, bool invert, SEXP n_values){
  int32_t NP = 0;
  R_xlen_t n = vec::length(x);
  bool is_long = (n > r_limits::r_int_max);
  if (vec::length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }

  if (implicit_na_coercion(value, x)){
    YIELD(NP);
    Rf_error("Value has been implicitly converted to NA, please check");
  }
  SHIELD(n_values = cast<r_doubles_t>(n_values, r_null)); ++NP;
  R_xlen_t n_vals = is_null(n_values) ? scalar_count(x, value, false) : real_ptr(n_values)[0];
  R_xlen_t out_size = invert ? n - n_vals : n_vals;
  R_xlen_t whichi = 0;
  R_xlen_t i = 0;

  switch ( CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    SEXP out = SHIELD(is_long ? new_vector<r_double_t>(out_size) : new_vector<r_int_t>(out_size)); ++NP;
    SHIELD(value = cast<r_integers_t>(value, r_null)); ++NP;
    int val = vector_ptr<r_int_t>(value)[0];
    const r_int_t *p_x = vector_ptr<r_int_t>(x);
    if (is_long){
      r_double_t* RESTRICT p_out = real_ptr(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      r_int_t* RESTRICT p_out = vector_ptr<r_int_t>(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    }
    YIELD(NP);
    return out;
  }
  case REALSXP: {
    SEXP out = SHIELD(is_long ? new_vector<r_double_t>(out_size) : new_vector<r_int_t>(out_size)); ++NP;
    SHIELD(value = cast<r_doubles_t>(value, r_null)); ++NP;
    double val = real_ptr(value)[0];
    const r_double_t *p_x = real_ptr(x);
    if (is_long){
      r_double_t* RESTRICT p_out = real_ptr(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      r_int_t* RESTRICT p_out = vector_ptr<r_int_t>(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    }
    YIELD(NP);
    return out;
  }
  case CHEAPR_INT64SXP: {
    SEXP out = SHIELD(is_long ? new_vector<r_double_t>(out_size) : new_vector<r_int_t>(out_size)); ++NP;
    SHIELD(value = cast<r_integers64_t>(value, r_null)); ++NP;
    int64_t val = integer64_ptr(value)[0];
    const int64_t *p_x = integer64_ptr_ro(x);
    if (is_long){
      r_double_t* RESTRICT p_out = real_ptr(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      r_int_t* RESTRICT p_out = vector_ptr<r_int_t>(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    }
    YIELD(NP);
    return out;
  }
  case STRSXP: {
    SEXP out = SHIELD(is_long ? new_vector<r_double_t>(out_size) : new_vector<r_int_t>(out_size)); ++NP;
    SHIELD(value = cast<r_characters_t>(value, r_null)); ++NP;
    r_string_t val = get_value<r_string_t>(value, 0);
    const r_string_t *p_x = string_ptr_ro(x);
    if (is_long){
      r_double_t *p_out = real_ptr(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      r_int_t *p_out = vector_ptr<r_int_t>(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    }
    YIELD(NP);
    return out;
  }
  case CPLXSXP: {
    SEXP out = SHIELD(is_long ? new_vector<r_double_t>(out_size) : new_vector<r_int_t>(out_size)); ++NP;
    SHIELD(value = cast<r_complexes_t>(value, r_null)); ++NP;
    r_complex_t val = complex_ptr(value)[0];
    const r_complex_t *p_x = complex_ptr_ro(x);
    if (is_long){
      r_double_t* RESTRICT p_out = real_ptr(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      r_int_t* RESTRICT p_out = vector_ptr<r_int_t>(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    }
    YIELD(NP);
    return out;
  }
  default: {
    SEXP is_equal;
    if (cpp_all_na(value, true, false)){
      is_equal = SHIELD(eval_pkg_fun("is_na", "cheapr", env::base_env, x)); ++NP;
    } else {
      is_equal = SHIELD(eval_pkg_fun("==", "base", env::base_env, x, value)); ++NP;
    }
    SEXP out = SHIELD(cpp_which_(is_equal, invert)); ++NP;
    YIELD(NP);
    return out;
  }
  }
}

[[cpp11::register]]
SEXP cpp_which_val(SEXP x, SEXP value, bool invert){
  SEXP n_vals = SHIELD(as_vector(scalar_count(x, value, false)));
  SEXP out = SHIELD(cpp_val_find(x, value, invert, n_vals));
  YIELD(2);
  return out;
}

// Memory-efficient which(is.na(x))

[[cpp11::register]]
SEXP cpp_which_na(SEXP x){
  SEXP na = SHIELD(as_vector(na::integer));
  SEXP out = SHIELD(cpp_which_val(x, na, false));
  YIELD(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_which_not_na(SEXP x){
  SEXP na = SHIELD(as_vector(na::integer));
  SEXP out = SHIELD(cpp_which_val(x, na, true));
  YIELD(2);
  return out;
}

// Return the locations of T, F, and NA in one pass
// Must provide the correct num of T and F as args

[[cpp11::register]]
SEXP cpp_lgl_locs(SEXP x, R_xlen_t n_true, R_xlen_t n_false,
                  bool include_true, bool include_false, bool include_na){
  R_xlen_t n = Rf_xlength(x);
  const r_int_t *p_x = vector_ptr<const r_int_t>(x);

  if (n > r_limits::r_int_max){
    SEXP true_locs = SHIELD(new_vector<r_double_t>(include_true ? n_true : 0));
    SEXP false_locs = SHIELD(new_vector<r_double_t>(include_false ? n_false : 0));
    SEXP na_locs = SHIELD(new_vector<r_double_t>(include_na ? (n - n_true - n_false) : 0));

    r_double_t* RESTRICT p_true = real_ptr(true_locs);
    r_double_t* RESTRICT p_false = real_ptr(false_locs);
    r_double_t* RESTRICT p_na = real_ptr(na_locs);

    R_xlen_t k1 = 0;
    R_xlen_t k2 = 0;
    R_xlen_t k3 = 0;

    for (R_xlen_t i = 0; i < n; ++i){
      if (include_true && p_x[i] == 1){
        p_true[k1++] = i + 1;
      } else if (include_false && p_x[i] == 0){
        p_false[k2++] = i + 1;
      } else if (include_na && is_r_na(p_x[i])){
        p_na[k3++] = i + 1;
      }
    }

    SEXP out = SHIELD(make_list(
      arg("true") = true_locs,
      arg("false") = false_locs,
      arg("na") = na_locs
    ));
    YIELD(4);
    return out;
  } else {
    SEXP true_locs = SHIELD(vec::new_vector<r_int_t>(include_true ? n_true : 0));
    SEXP false_locs = SHIELD(vec::new_vector<r_int_t>(include_false ? n_false : 0));
    SEXP na_locs = SHIELD(vec::new_vector<r_int_t>(include_na ? (n - n_true - n_false) : 0));

    r_int_t* RESTRICT p_true = vector_ptr<r_int_t>(true_locs);
    r_int_t* RESTRICT p_false = vector_ptr<r_int_t>(false_locs);
    r_int_t* RESTRICT p_na = vector_ptr<r_int_t>(na_locs);

    int k1 = 0;
    int k2 = 0;
    int k3 = 0;

    for (int i = 0; i < n; ++i){
      if (include_true && p_x[i] == 1){
        p_true[k1++] = i + 1;
      } else if (include_false && p_x[i] == 0){
        p_false[k2++] = i + 1;
      } else if (include_na && is_r_na(p_x[i])){
        p_na[k3++] = i + 1;
      }
    }
    SEXP out = SHIELD(make_list(
      arg("true") = true_locs,
      arg("false") = false_locs,
      arg("na") = na_locs
    ));

    YIELD(4);
    return out;
  }
}
