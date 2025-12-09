#include "cheapr.h"

// A more memory-efficient which()
// Author: Nick Christofides
// License: MIT

// Count the number of true values

R_xlen_t count_true(const r_bool_t* RESTRICT px, const uint_fast64_t n){
  uint_fast64_t size = 0;
  if (n >= CHEAPR_OMP_THRESHOLD){
#pragma omp parallel for simd num_threads(num_cores()) reduction(+:size)
    for (uint_fast64_t j = 0; j != n; ++j) size += static_cast<uint_fast64_t>(px[j] == r_true);
    return size;
  } else {
    OMP_FOR_SIMD
    for (uint_fast64_t j = 0; j != n; ++j) size += static_cast<uint_fast64_t>(px[j] == r_true);
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
    whichi += eq(p_x[i++], VAL);                                     \
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
    whichi += !eq(p_x[i++], VAL);                                  \
  }                                                                \
}


[[cpp11::register]]
SEXP cpp_which_(SEXP x, bool invert){
  R_xlen_t n = Rf_xlength(x);
  const r_bool_t *p_x = BOOLEAN_RO(x);
  bool is_long = (n > limits::r_int_max);
  if (invert){
    if (is_long){
      R_xlen_t size = count_true(p_x, n);
      R_xlen_t out_size = n - size;
      SEXP out = SHIELD(internal::new_vec(REALSXP, out_size));
      double* RESTRICT p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL_INVERTED(r_true);
      YIELD(1);
      return out;
    } else {
      int size = count_true(p_x, n);
      int out_size = n - size;
      SEXP out = SHIELD(vec::new_integer(out_size));
      int* RESTRICT p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      CHEAPR_WHICH_VAL_INVERTED(r_true);
      YIELD(1);
      return out;
    }
  } else {
    if (is_long){
      R_xlen_t out_size = count_true(p_x, n);
      SEXP out = SHIELD(internal::new_vec(REALSXP, out_size));
      double* RESTRICT p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      CHEAPR_WHICH_VAL(r_true);
      YIELD(1);
      return out;
    } else {
      int out_size = count_true(p_x, n);
      SEXP out = SHIELD(vec::new_integer(out_size));
      int* RESTRICT p_out = INTEGER(out);
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
  bool is_long = (n > limits::r_int_max);
  if (vec::length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }

  if (implicit_na_coercion(value, x)){
    YIELD(NP);
    Rf_error("Value has been implicitly converted to NA, please check");
  }
  SHIELD(n_values = cast<r_double_t>(n_values, r_null)); ++NP;
  R_xlen_t n_vals = is_null(n_values) ? scalar_count(x, value, false) : REAL(n_values)[0];
  R_xlen_t out_size = invert ? n - n_vals : n_vals;
  R_xlen_t whichi = 0;
  R_xlen_t i = 0;

  switch ( CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    SEXP out = SHIELD(internal::new_vec(is_long ? REALSXP : INTSXP, out_size)); ++NP;
    SHIELD(value = cast<r_integer_t>(value, r_null)); ++NP;
    int val = INTEGER(value)[0];
    const int *p_x = INTEGER(x);
    if (is_long){
      double* RESTRICT p_out = REAL(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int* RESTRICT p_out = INTEGER(out);
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
    SEXP out = SHIELD(internal::new_vec(is_long ? REALSXP : INTSXP, out_size)); ++NP;
    SHIELD(value = cast<r_double_t>(value, r_null)); ++NP;
    double val = REAL(value)[0];
    const double *p_x = REAL(x);
    if (is_long){
      double* RESTRICT p_out = REAL(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int* RESTRICT p_out = INTEGER(out);
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
    SEXP out = SHIELD(internal::new_vec(is_long ? REALSXP : INTSXP, out_size)); ++NP;
    SHIELD(value = cast<r_integer64_t>(value, r_null)); ++NP;
    int64_t val = INTEGER64_PTR(value)[0];
    const int64_t *p_x = INTEGER64_PTR_RO(x);
    if (is_long){
      double* RESTRICT p_out = REAL(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int* RESTRICT p_out = INTEGER(out);
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
    SEXP out = SHIELD(internal::new_vec(is_long ? REALSXP : INTSXP, out_size)); ++NP;
    SHIELD(value = cast<r_character_t>(value, r_null)); ++NP;
    SEXP val = SHIELD(STRING_ELT(value, 0)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);
    if (is_long){
      double *p_out = REAL(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int *p_out = INTEGER(out);
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
    SEXP out = SHIELD(internal::new_vec(is_long ? REALSXP : INTSXP, out_size)); ++NP;
    SHIELD(value = cast<r_complex_t>(value, r_null)); ++NP;
    Rcomplex val = as_complex(COMPLEX(value)[0]);
    const Rcomplex *p_x = COMPLEX_RO(x);
    if (is_long){
      double* RESTRICT p_out = REAL(out);
      if (invert){
        CHEAPR_WHICH_VAL_INVERTED(val);
      } else {
        CHEAPR_WHICH_VAL(val);
      }
    } else {
      int* RESTRICT p_out = INTEGER(out);
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
      is_equal = SHIELD(eval_pkg_fun("is_na", "cheapr", R_GetCurrentEnv(), x)); ++NP;
    } else {
      is_equal = SHIELD(eval_pkg_fun("==", "base", R_GetCurrentEnv(), x, value)); ++NP;
    }
    SEXP out = SHIELD(cpp_which_(is_equal, invert)); ++NP;
    YIELD(NP);
    return out;
  }
  }
}

[[cpp11::register]]
SEXP cpp_which_val(SEXP x, SEXP value, bool invert){
  SEXP n_vals = SHIELD(as_vec(scalar_count(x, value, false)));
  SEXP out = SHIELD(cpp_val_find(x, value, invert, n_vals));
  YIELD(2);
  return out;
}

// Memory-efficient which(is.na(x))

[[cpp11::register]]
SEXP cpp_which_na(SEXP x){
  SEXP na = SHIELD(as_vec(na::integer));
  SEXP out = SHIELD(cpp_which_val(x, na, false));
  YIELD(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_which_not_na(SEXP x){
  SEXP na = SHIELD(as_vec(na::integer));
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
  const int *p_x = INTEGER_RO(x);

  if (n > limits::r_int_max){
    SEXP true_locs = SHIELD(internal::new_vec(REALSXP, include_true ? n_true : 0));
    SEXP false_locs = SHIELD(internal::new_vec(REALSXP, include_false ? n_false : 0));
    SEXP na_locs = SHIELD(internal::new_vec(REALSXP, include_na ? (n - n_true - n_false) : 0));

    double* RESTRICT p_true = REAL(true_locs);
    double* RESTRICT p_false = REAL(false_locs);
    double* RESTRICT p_na = REAL(na_locs);

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
    SEXP true_locs = SHIELD(vec::new_integer(include_true ? n_true : 0));
    SEXP false_locs = SHIELD(vec::new_integer(include_false ? n_false : 0));
    SEXP na_locs = SHIELD(vec::new_integer(include_na ? (n - n_true - n_false) : 0));

    int* RESTRICT p_true = INTEGER(true_locs);
    int* RESTRICT p_false = INTEGER(false_locs);
    int* RESTRICT p_na = INTEGER(na_locs);

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
