#include "cheapr.h"

// Scalar-optimised functions for R
// Author: Nick Christofides

void check_atomic(SEXP x){
  if (!vec::is_atomic(x)){
    Rf_error("'cheapr' scalar functions can only accept atomic vectors");
  }
}

bool implicit_na_coercion(SEXP x, SEXP target){
  SEXP coerced = SHIELD(cast_(get_r_type(target), x, target));
  bool out = na_count(x, true) != na_count(coerced, true);
  YIELD(1);
  return out;
}

#define CHEAPR_VAL_COUNT(VAL)                                  \
if (is_r_na(VAL)){                                             \
  _Pragma("omp simd reduction(+:count)")                       \
  for (R_xlen_t i = 0; i < n; ++i){                            \
    count += is_r_na(p_x[i]);                                  \
  }                                                            \
} else {                                                       \
  _Pragma("omp simd reduction(+:count)")                       \
  for (R_xlen_t i = 0; i < n; ++i){                            \
    count += (p_x[i] == VAL);                                  \
  }                                                            \
}


R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive){
  if (vec::length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t count = 0;
  int32_t NP = 0;

  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_integers_t>(value, r_null)); ++NP;
    const int val = integer_ptr(value)[0];
    const int *p_x = integer_ptr(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_doubles_t>(value, r_null)); ++NP;
    const double val = real_ptr(value)[0];
    const double *p_x = real_ptr(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case CHEAPR_INT64SXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_integers64_t>(value, r_null)); ++NP;
    const int64_t val = integer64_ptr(value)[0];
    const int64_t *p_x = integer64_ptr_ro(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_characters_t>(value, r_null)); ++NP;
    const r_string_t val = get_value<r_string_t>(value, 0);
    const r_string_t *p_x = string_ptr_ro(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case CPLXSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_complexes_t>(value, r_null)); ++NP;
    const r_complex_t val = complex_ptr(value)[0];
    const r_complex_t *p_x = complex_ptr_ro(x);
    CHEAPR_VAL_COUNT(val);
    break;
  }
  case VECSXP: {
    if (recursive){
    R_CheckStack(); // Check C Stack size isn't close to the limit
    const SEXP *p_x = list_ptr_ro(x);
    for (R_xlen_t i = 0; i < n; ++i){
      count += scalar_count(p_x[i], value, true);
    }
    break;
  }
  }
  default: {
    if (cpp_all_na(value, true, false)){
    count = na_count(x, false);
  } else {
    SEXP is_equal = SHIELD(eval_pkg_fun("==", "base", env::base_env, x, value)); ++NP;
    count = cpp_sum(is_equal);
  }
  break;
  }
  }
  YIELD(NP);
  return count;
}

[[cpp11::register]]
SEXP cpp_count_val(SEXP x, SEXP value, bool recursive){
  return as_vector(scalar_count(x, value, recursive));
}

[[cpp11::register]]
SEXP cpp_val_replace(SEXP x, SEXP value, SEXP replace, bool recursive){
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);

  if (Rf_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  if (Rf_length(replace) != 1){
    Rf_error("replace must be a vector of length 1");
  }
  bool val_is_na = cpp_any_na(value, true);
  bool any_eq = false;
  bool eql = false;

  SEXP out = x;

  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    int *p_out = integer_ptr(out);
    SHIELD(value = cast<r_integers_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_integers_t>(replace, r_null)); ++NP;
    int val = integer_ptr(value)[0];
    int repl = integer_ptr(replace)[0];
    int *p_x = integer_ptr(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eql = p_x[i] == val;
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(out)); ++NP;
        // Change where pointer is pointing to
        p_out = integer_ptr(out);
      }
      if (eql) p_out[i] = repl;
    }
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    SHIELD(value = cast<r_doubles_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_doubles_t>(replace, r_null)); ++NP;
    double *p_out = real_ptr(out);
    double val = real_ptr(value)[0];
    double repl = real_ptr(replace)[0];
    double *p_x = real_ptr(x);
    if (val_is_na){
      for (R_xlen_t i = 0; i < n; ++i){
        eql = is_r_na(p_x[i]);
        if (!any_eq && eql){
          any_eq = true;
          out = SHIELD(cpp_semi_copy(out)); ++NP;
          // Change where pointer is pointing to
          p_out = real_ptr(out);
        }
        if (eql) p_out[i] = repl;
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i){
        eql = p_x[i] == val;
        if (!any_eq && eql){
          any_eq = true;
          out = SHIELD(cpp_semi_copy(out)); ++NP;
          // Change where pointer is pointing to
          p_out = real_ptr(out);
        }
        if (eql) p_out[i] = repl;
      }
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    SHIELD(value = cast<r_integers64_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_integers64_t>(replace, r_null)); ++NP;
    int64_t *p_out = integer64_ptr(out);
    int64_t val = integer64_ptr(value)[0];
    int64_t repl = integer64_ptr(replace)[0];
    int64_t *p_x = integer64_ptr(x);
    for (R_xlen_t i = 0; i < n; ++i){
      eql = p_x[i] == val;
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
        // Change where pointer is pointing to
        p_out = integer64_ptr(out);
      }
      if (eql) p_out[i] = repl;
    }
    break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    SHIELD(value = cast<r_characters_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_characters_t>(replace, r_null)); ++NP;
    r_string_t val = get_value<r_string_t>(value, 0);
    r_string_t repl = get_value<r_string_t>(replace, 0);
    const r_string_t *p_x = string_ptr_ro(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eql = p_x[i] == val;
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(out)); ++NP;
      }
      if (eql) set_value<r_string_t>(out, i, repl);
    }
    break;
  }
  case CPLXSXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    SHIELD(value = cast<r_complexes_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_complexes_t>(replace, r_null)); ++NP;
    r_complex_t val = get_value<r_complex_t>(value, 0);
    r_complex_t repl = get_value<r_complex_t>(replace, 0);
    const r_complex_t *p_x = complex_ptr_ro(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eql = val_is_na ? is_r_na(p_x[i]) : (p_x[i] == val);
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(out)); ++NP;
      }
      if (eql) set_value<r_complex_t>(out, i, repl);
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    R_CheckStack(); // Check C Stack size isn't close to the limit
    SHIELD(out = vec::shallow_copy(out)); ++NP;
    for (R_xlen_t i = 0; i < n; ++i){
      // Once we extract the vector it maybe needs protecting??
      SET_VECTOR_ELT(out, i, cpp_val_replace(VECTOR_ELT(out, i), value, replace, true));
    }
    break;
  }
  }
  default: {
    YIELD(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  YIELD(NP);
  return out;
}

SEXP cpp_loc_set_replace(SEXP x, SEXP where, SEXP what){
  return cpp_replace(x, where, what, true, true);
}

// Remove elements from a vector very efficiently

[[cpp11::register]]
SEXP cpp_val_remove(SEXP x, SEXP value, bool recursive){
  int32_t NP = 0;
  R_xlen_t n_vals = scalar_count(x, value, true);
  if (n_vals == 0){
    return x;
  } else if (n_vals == Rf_xlength(x)){
    SEXP out = SHIELD(internal::new_vec(TYPEOF(x), 0)); ++NP;
    SEXP attrs = SHIELD(attr::get_attrs(x)); ++NP;
    attr::set_attrs(out, attrs);
    YIELD(NP);
    return out;
  } else {
    R_xlen_t n = Rf_xlength(x);
    R_xlen_t n_keep = n - n_vals;
    R_xlen_t k = 0;

    SEXP out = x;
    switch ( CHEAPR_TYPEOF(x) ){
    case NILSXP: {
      break;
    }
    case LGLSXP:
    case INTSXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(internal::new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = cast<r_integers_t>(value, r_null)); ++NP;
      int val = integer_ptr(value)[0];
      int *p_x = integer_ptr(x);
      int *p_out = integer_ptr(out);

      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] != val){
          p_out[k++] = p_x[i];
        }
      }
      SEXP attrs = SHIELD(attr::get_attrs(x)); ++NP;
      attr::set_attrs(out, attrs);
      break;
    }
    case REALSXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(internal::new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = cast<r_doubles_t>(value, r_null)); ++NP;
      double val = real_ptr(value)[0];
      double *p_x = real_ptr(x);
      double *p_out = real_ptr(out);

      if (cpp_any_na(value, true)){
        for (R_xlen_t i = 0; i < n; ++i){
          if (!is_r_na(p_x[i])){
            p_out[k++] = p_x[i];
          }
        }
      } else {
        for (R_xlen_t i = 0; i < n; ++i){
          if (p_x[i] != val){
            p_out[k++] = p_x[i];
          }
        }
      }
      SEXP attrs = SHIELD(attr::get_attrs(x)); ++NP;
      attr::set_attrs(out, attrs);
      break;
    }
    case CHEAPR_INT64SXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(new_vector<double>(n_keep)); ++NP;
      SHIELD(value = cast<r_integers64_t>(value, r_null)); ++NP;
      int64_t val = integer64_ptr(value)[0];
      int64_t *p_x = integer64_ptr(x);
      int64_t *p_out = integer64_ptr(out);

      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] != val){
          p_out[k++] = p_x[i];
        }
      }
      SEXP attrs = SHIELD(attr::get_attrs(x)); ++NP;
      attr::set_attrs(out, attrs);
      break;
    }
    case STRSXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(new_vector<r_string_t>(n_keep)); ++NP;
      SHIELD(value = cast<r_characters_t>(value, r_null)); ++NP;
      SEXP val = SHIELD(STRING_ELT(value, 0)); ++NP;
      const r_string_t *p_x = string_ptr_ro(x);

      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] != val){
          set_value<r_string_t>(out, k++, p_x[i]);
        }
      }
      SEXP attrs = SHIELD(attr::get_attrs(x)); ++NP;
      attr::set_attrs(out, attrs);
      break;
    }
    case CPLXSXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(new_vector<r_complex_t>(n_keep)); ++NP;
      SHIELD(value = cast<r_complexes_t>(value, r_null)); ++NP;
      r_complex_t val = get_value<r_complex_t>(value, 0);
      const r_complex_t *p_x = complex_ptr_ro(x);

      for (R_xlen_t i = 0; i < n; ++i){
        if ( !(is_r_na(val) ? is_r_na(p_x[i]) : (p_x[i] == val)) ){
          set_value<r_complex_t>(out, k++, p_x[i]);
        }
      }
      SEXP attrs = SHIELD(attr::get_attrs(x)); ++NP;
      attr::set_attrs(out, attrs);
      break;
    }
    case VECSXP: {
      if (recursive){
      R_CheckStack(); // Check C Stack size isn't close to the limit
      SHIELD(out = vec::shallow_copy(out)); ++NP;
      attr::clear_attrs(out);
      for (R_xlen_t i = 0; i < n; ++i){
        SET_VECTOR_ELT(out, i, cpp_val_remove(VECTOR_ELT(out, i), value, true));
      }
      break;
    }
    }
    default: {
      SEXP sexp_n_vals = SHIELD(as_vector(r_cast<double>(n_vals))); ++NP;
      SEXP val_locs = SHIELD(cpp_val_find(x, value, true, sexp_n_vals)); ++NP;
      out = SHIELD(eval_pkg_fun("cheapr_sset", "cheapr", env::base_env, x, val_locs)); ++NP;
      break;
    }
    }
    YIELD(NP);
    return out;
  }
}
