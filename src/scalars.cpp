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
  for (R_xlen_t i = 0; i < n; ++i){                            \
    count += is_r_na(p_x[i]);                                  \
  }                                                            \
} else {                                                       \
  for (R_xlen_t i = 0; i < n; ++i){                            \
    count += eq(p_x[i], VAL);                                  \
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
    SHIELD(value = cast<r_integer_t>(value, r_null)); ++NP;
    const int val = INTEGER(value)[0];
    const int *p_x = INTEGER(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_double_t>(value, r_null)); ++NP;
    const double val = REAL(value)[0];
    const double *p_x = REAL(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case CHEAPR_INT64SXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_integer64_t>(value, r_null)); ++NP;
    const int64_t val = INTEGER64_PTR(value)[0];
    const int64_t *p_x = INTEGER64_PTR_RO(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_character_t>(value, r_null)); ++NP;
    const SEXP val = STRING_ELT(value, 0);
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_VAL_COUNT(val)
      break;
  }
  case CPLXSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_complex_t>(value, r_null)); ++NP;
    const Rcomplex val = as_complex(COMPLEX(value)[0]);
    const Rcomplex *p_x = COMPLEX_RO(x);
    CHEAPR_VAL_COUNT(val);
    break;
  }
  case VECSXP: {
    if (recursive){
    const SEXP *p_x = LIST_PTR_RO(x);
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
    SEXP is_equal = SHIELD(eval_pkg_fun("==", "base", R_GetCurrentEnv(), x, value)); ++NP;
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
  return as_vec(scalar_count(x, value, recursive));
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
    int *p_out = INTEGER(out);
    SHIELD(value = cast<r_integer_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_integer_t>(replace, r_null)); ++NP;
    int val = INTEGER(value)[0];
    int repl = INTEGER(replace)[0];
    int *p_x = INTEGER(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eql = p_x[i] == val;
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(out)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER(out);
      }
      if (eql) p_out[i] = repl;
    }
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    SHIELD(value = cast<r_double_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_double_t>(replace, r_null)); ++NP;
    double *p_out = REAL(out);
    double val = REAL(value)[0];
    double repl = REAL(replace)[0];
    double *p_x = REAL(x);
    if (val_is_na){
      for (R_xlen_t i = 0; i < n; ++i){
        eql = is_r_na(p_x[i]);
        if (!any_eq && eql){
          any_eq = true;
          out = SHIELD(cpp_semi_copy(out)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
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
          p_out = REAL(out);
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
    SHIELD(value = cast<r_integer64_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_integer64_t>(replace, r_null)); ++NP;
    int64_t *p_out = INTEGER64_PTR(out);
    int64_t val = INTEGER64_PTR(value)[0];
    int64_t repl = INTEGER64_PTR(replace)[0];
    int64_t *p_x = INTEGER64_PTR(x);
    for (R_xlen_t i = 0; i < n; ++i){
      eql = p_x[i] == val;
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER64_PTR(out);
      }
      if (eql) p_out[i] = repl;
    }
    break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    SHIELD(value = cast<r_character_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_character_t>(replace, r_null)); ++NP;
    SEXP val = SHIELD(STRING_ELT(value, 0)); ++NP;
    SEXP repl = SHIELD(STRING_ELT(replace, 0)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eql = p_x[i] == val;
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(out)); ++NP;
      }
      if (eql) SET_STRING_ELT(out, i, repl);
    }
    break;
  }
  case CPLXSXP: {
    if (implicit_na_coercion(value, x)){
    break;
  }
    SHIELD(value = cast<r_complex_t>(value, r_null)); ++NP;
    SHIELD(replace = cast<r_complex_t>(replace, r_null)); ++NP;
    Rcomplex val = COMPLEX_ELT(value, 0);
    Rcomplex repl = COMPLEX_ELT(replace, 0);
    const Rcomplex *p_x = COMPLEX_RO(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eql = val_is_na ? is_r_na(p_x[i]) : eq(p_x[i], val);
      if (!any_eq && eql){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(out)); ++NP;
      }
      if (eql) SET_COMPLEX_ELT(out, i, repl);
    }
    break;
  }
  case VECSXP: {
    if (recursive){
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
    cpp_set_add_attributes(out, ATTRIB(x), false);
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
      SHIELD(value = cast<r_integer_t>(value, r_null)); ++NP;
      int val = INTEGER(value)[0];
      int *p_x = INTEGER(x);
      int *p_out = INTEGER(out);

      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] != val){
          p_out[k++] = p_x[i];
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case REALSXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(internal::new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = cast<r_double_t>(value, r_null)); ++NP;
      double val = REAL(value)[0];
      double *p_x = REAL(x);
      double *p_out = REAL(out);

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
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case CHEAPR_INT64SXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(new_double(n_keep)); ++NP;
      SHIELD(value = cast<r_integer64_t>(value, r_null)); ++NP;
      int64_t val = INTEGER64_PTR(value)[0];
      int64_t *p_x = INTEGER64_PTR(x);
      int64_t *p_out = INTEGER64_PTR(out);

      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] != val){
          p_out[k++] = p_x[i];
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case STRSXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(new_character(n_keep)); ++NP;
      SHIELD(value = cast<r_character_t>(value, r_null)); ++NP;
      SEXP val = SHIELD(STRING_ELT(value, 0)); ++NP;
      const SEXP *p_x = STRING_PTR_RO(x);

      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] != val){
          SET_STRING_ELT(out, k++, p_x[i]);
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case CPLXSXP: {
      if (implicit_na_coercion(value, x)){
      break;
    }
      out = SHIELD(new_complex(n_keep)); ++NP;
      SHIELD(value = cast<r_complex_t>(value, r_null)); ++NP;
      Rcomplex val = COMPLEX_ELT(value, 0);
      const Rcomplex *p_x = COMPLEX_RO(x);

      for (R_xlen_t i = 0; i < n; ++i){
        if ( !(is_r_na(val) ? is_r_na(p_x[i]) : eq(p_x[i], val)) ){
          SET_COMPLEX_ELT(out, k++, p_x[i]);
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case VECSXP: {
      if (recursive){
      SHIELD(out = vec::shallow_copy(out)); ++NP;
      attr::clear_attrs(out);
      for (R_xlen_t i = 0; i < n; ++i){
        SET_VECTOR_ELT(out, i, cpp_val_remove(VECTOR_ELT(out, i), value, true));
      }
      break;
    }
    }
    default: {
      SEXP sexp_n_vals = SHIELD(as_vec(r_cast<double>(n_vals))); ++NP;
      SEXP val_locs = SHIELD(cpp_val_find(x, value, true, sexp_n_vals)); ++NP;
      out = SHIELD(eval_pkg_fun("cheapr_sset", "cheapr", R_GetCurrentEnv(), x, val_locs)); ++NP;
      break;
    }
    }
    YIELD(NP);
    return out;
  }
}
