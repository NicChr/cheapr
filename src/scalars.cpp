#include "cheapr.h"

// Scalar-optimised functions for R
// Author: Nick Christofides

void check_atomic(SEXP x){
  if (!Rf_isVectorAtomic(x)){
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
  if (vector_length(value) != 1){
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
    SHIELD(value = cast<r_integer_t>(value, R_NilValue)); ++NP;
    const int val = INTEGER(value)[0];
    const int *p_x = INTEGER(x);
    CHEAPR_VAL_COUNT(val)
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_numeric_t>(value, R_NilValue)); ++NP;
    const double val = REAL(value)[0];
    const double *p_x = REAL(x);
    CHEAPR_VAL_COUNT(val)
    break;
  }
  case CHEAPR_INT64SXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_integer64_t>(value, R_NilValue)); ++NP;
    const int64_t val = INTEGER64_PTR(value)[0];
    const int64_t *p_x = INTEGER64_PTR_RO(x);
    CHEAPR_VAL_COUNT(val)
    break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_character_t>(value, R_NilValue)); ++NP;
    const SEXP val = STRING_ELT(value, 0);
    const SEXP *p_x = STRING_PTR_RO(x);
    CHEAPR_VAL_COUNT(val)
    break;
  }
  case CPLXSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = cast<r_complex_t>(value, R_NilValue)); ++NP;
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
    SEXP expr = SHIELD(Rf_lang3(install_utf8("=="), x, value)); ++NP;
    SEXP is_equal = SHIELD(Rf_eval(expr, R_GetCurrentEnv())); ++NP;
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
  return as_r_scalar(scalar_count(x, value, recursive));
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
  bool eq = false;

  SEXP out = SHIELD(R_NilValue); ++NP;
  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SEXP temp = SHIELD(new_vec(INTSXP, 0)); ++NP;
    int *p_out = INTEGER(temp);
    SHIELD(value = cast<r_integer_t>(value, R_NilValue)); ++NP;
    SHIELD(replace = cast<r_integer_t>(replace, R_NilValue)); ++NP;
    int val = Rf_asInteger(value);
    int repl = Rf_asInteger(replace);
    int *p_x = INTEGER(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER(out);
      }
      if (eq) p_out[i] = repl;
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = SHIELD(x); ++NP;
    }
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SHIELD(value = cast<r_numeric_t>(value, R_NilValue)); ++NP;
    SHIELD(replace = cast<r_numeric_t>(replace, R_NilValue)); ++NP;
    SEXP temp = SHIELD(new_vec(REALSXP, 0)); ++NP;
    double *p_out = REAL(temp);
    double val = Rf_asReal(value);
    double repl = Rf_asReal(replace);
    double *p_x = REAL(x);
    if (val_is_na){
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] != p_x[i];
        if (!any_eq && eq){
          any_eq = true;
          out = SHIELD(cpp_semi_copy(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = SHIELD(x); ++NP;
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!any_eq && eq){
          any_eq = true;
          out = SHIELD(cpp_semi_copy(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = SHIELD(x); ++NP;
      }
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SHIELD(value = cast<r_integer64_t>(value, R_NilValue)); ++NP;
    SHIELD(replace = cast<r_integer64_t>(replace, R_NilValue)); ++NP;
    SEXP temp = SHIELD(new_vec(REALSXP, 0)); ++NP;
    int64_t *p_out = INTEGER64_PTR(temp);
    int64_t val = INTEGER64_PTR(value)[0];
    int64_t repl = INTEGER64_PTR(replace)[0];
    int64_t *p_x = INTEGER64_PTR(x);
    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER64_PTR(out);
      }
      if (eq) p_out[i] = repl;
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = SHIELD(x); ++NP;
    }
    break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SHIELD(value = cast<r_character_t>(value, R_NilValue)); ++NP;
    SHIELD(replace = cast<r_character_t>(replace, R_NilValue)); ++NP;
    SEXP val = SHIELD(Rf_asChar(value)); ++NP;
    SEXP repl = SHIELD(Rf_asChar(replace)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
      }
      if (eq) SET_STRING_ELT(out, i, repl);
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = SHIELD(x); ++NP;
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    out = SHIELD(new_vec(VECSXP, n)); ++NP;
    for (R_xlen_t i = 0; i < n; ++i){
      // Initialise each element
      SET_VECTOR_ELT(out, i, VECTOR_ELT(x, i));

      // Once we extract the vector it maybe needs protecting??
      SET_VECTOR_ELT(out, i, cpp_val_replace(VECTOR_ELT(out, i), value, replace, true));
    }
    // Copy attributes from x to y
    SHALLOW_DUPLICATE_ATTRIB(out, x);
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

// At the moment this doesn't coerce what to type of x
// or handle long vectors

SEXP cpp_loc_set_replace(SEXP x, SEXP where, SEXP what){
  return cpp_replace(x, where, what, true, true);
}

// Remove elements from a vector very efficiently

[[cpp11::register]]
SEXP cpp_val_remove(SEXP x, SEXP value){
  check_atomic(x);
  int32_t NP = 0;
  R_xlen_t n_vals = scalar_count(x, value, true);
  if (n_vals == 0){
    return x;
  } else if (n_vals == Rf_xlength(x)){
    SEXP out = SHIELD(new_vec(TYPEOF(x), 0)); ++NP;
    cpp_set_add_attributes(out, ATTRIB(x), false);
    YIELD(NP);
    return out;
  } else {
    R_xlen_t n = Rf_xlength(x);
    R_xlen_t n_keep = n - n_vals;
    bool eq;
    R_xlen_t k = 0;

    SEXP out;
    switch ( CHEAPR_TYPEOF(x) ){
    case NILSXP: {
      out = SHIELD(R_NilValue); ++NP;
      break;
    }
    case LGLSXP:
    case INTSXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = cast<r_integer_t>(value, R_NilValue)); ++NP;
      int val = Rf_asInteger(value);
      int *p_x = INTEGER(x);
      int *p_out = INTEGER(out);

      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!eq){
          p_out[k++] = p_x[i];
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case REALSXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = cast<r_numeric_t>(value, R_NilValue)); ++NP;
      double val = Rf_asReal(value);
      double *p_x = REAL(x);
      double *p_out = REAL(out);

      if (cpp_any_na(value, true)){
        for (R_xlen_t i = 0; i < n; ++i){
          eq = is_r_na(p_x[i]);
          if (!eq){
            p_out[k++] = p_x[i];
          }
        }
      } else {
        for (R_xlen_t i = 0; i < n; ++i){
          eq = p_x[i] == val;
          if (!eq){
            p_out[k++] = p_x[i];
          }
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case CHEAPR_INT64SXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = cast<r_integer64_t>(value, R_NilValue)); ++NP;
      int64_t val = INTEGER64_PTR(value)[0];
      int64_t *p_x = INTEGER64_PTR(x);
      int64_t *p_out = INTEGER64_PTR(out);

      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!eq){
          p_out[k++] = p_x[i];
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case STRSXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = cast<r_character_t>(value, R_NilValue)); ++NP;
      SEXP val = SHIELD(Rf_asChar(value)); ++NP;
      const SEXP *p_x = STRING_PTR_RO(x);

      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!eq){
          SET_STRING_ELT(out, k++, p_x[i]);
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    default: {
      SEXP sexp_n_vals = SHIELD(as_r_scalar<double>(n_vals)); ++NP;
      SEXP val_locs = SHIELD(cpp_val_find(x, value, true, sexp_n_vals)); ++NP;
      out = SHIELD(cheapr_sset(x, val_locs)); ++NP;
      break;
    }
    }
    YIELD(NP);
    return out;
  }
}
