#include "cheapr_cpp.h"

R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive){
  if (cpp_vec_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t count = 0;
  int NP = 0;
  bool do_parallel = n >= 100000;
  int n_cores = do_parallel ? num_cores() : 1;

  SEXP val_is_na = Rf_protect(cpp_is_na(value));
  ++NP;
  if (Rf_length(val_is_na) == 1 && LOGICAL(val_is_na)[0]){
    Rf_unprotect(NP);
    return na_count(x, recursive);
  }
#define VAL_COUNT(_val_) for (R_xlen_t i = 0; i < n; ++i) count += (p_x[i] == _val_);
  switch ( TYPEOF(x) ){
  case NILSXP: {
    Rf_unprotect(NP);
    return count;
  }
  case LGLSXP:
  case INTSXP: {
    Rf_protect(value = Rf_coerceVector(value, INTSXP));
    ++NP;
    if (Rf_length(value) != 1){
      Rf_unprotect(NP);
      Rf_error("value must be of length 1");
    }
    // Basically if the coercion produces NA the count is 0.
    if (INTEGER(value)[0] == NA_INTEGER) break;
    int val = Rf_asInteger(value);
    int *p_x = INTEGER(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      VAL_COUNT(val)
    } else {
#pragma omp for simd
      VAL_COUNT(val)
    }
    break;
  }
  case REALSXP: {
    Rf_protect(value = Rf_coerceVector(value, REALSXP));
    ++NP;
    if (Rf_length(value) != 1){
      Rf_unprotect(NP);
      Rf_error("value must be of length 1");
    }
    // Basically if the coercion produces NA the count is 0.
    if (REAL(value)[0] != REAL(value)[0]) break;
    double val = Rf_asReal(value);
    double *p_x = REAL(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      VAL_COUNT(val)
    } else {
#pragma omp for simd
      VAL_COUNT(val)
    }
    break;
  }
  case STRSXP: {
    Rf_protect(value = Rf_coerceVector(value, STRSXP));
    ++NP;
    if (Rf_length(value) != 1){
      Rf_unprotect(NP);
      Rf_error("value must be of length 1");
    }
    // Basically if the coercion produces NA the count is 0.
    if (STRING_ELT(value, 0) == NA_STRING) break;
    SEXP val = Rf_protect(Rf_asChar(value));
    ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      VAL_COUNT(val);
    } else {
#pragma omp for simd
      VAL_COUNT(val);
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      count += scalar_count(p_x[i], value, true);
    }
    break;
  }
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
    // SEXP is_equal = Rf_protect(cpp11::package("base")["=="](x, value));
    // ++NP;
    // SEXP r_true = Rf_protect(Rf_ScalarLogical(true));
    // ++NP;
    // count = scalar_count(is_equal, r_true, true);
    // break;
  }
  }
  Rf_unprotect(NP);
  return count;
}

[[cpp11::register]]
SEXP cpp_count_val(SEXP x, SEXP value, bool recursive){
  return xlen_to_r(scalar_count(x, value, recursive));
}

// Quick search to return if x is in y vector

// bool x_in_y(int x, SEXP y){
//   int n = Rf_length(y);
//   bool out;
//   if (n == 1){
//     out = (x == Rf_asInteger(y));
//   } else {
//     int *p_y = INTEGER(y);
//     out = false;
//     for (int i = 0; i < n; ++i){
//       if (x == p_y[i]){
//         out = true;
//         break;
//       }
//     }
//   }
//   return out;
// }


[[cpp11::register]]
SEXP cpp_val_replace(SEXP x, SEXP value, SEXP replace, bool recursive){
  int NP = 0;
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

  SEXP out;
  switch ( TYPEOF(x) ){
  case NILSXP: {
    out = Rf_protect(R_NilValue); ++NP;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    SEXP temp = Rf_protect(Rf_allocVector(INTSXP, 0)); ++NP;
    int *p_out = INTEGER(temp);
    Rf_protect(value = Rf_coerceVector(value, INTSXP)); ++NP;
    Rf_protect(replace = Rf_coerceVector(replace, INTSXP)); ++NP;
    int val = Rf_asInteger(value);
    int repl = Rf_asInteger(replace);
    int *p_x = INTEGER(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = Rf_protect(Rf_duplicate(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER(out);
      }
      if (eq) p_out[i] = repl;
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = Rf_protect(x); ++NP;
    }
    break;
  }
  case REALSXP: {
    SEXP temp = Rf_protect(Rf_allocVector(REALSXP, 0)); ++NP;
    double *p_out = REAL(temp);
    Rf_protect(value = Rf_coerceVector(value, REALSXP)); ++NP;
    Rf_protect(replace = Rf_coerceVector(replace, REALSXP)); ++NP;
    double val = Rf_asReal(value);
    double repl = Rf_asReal(replace);
    double *p_x = REAL(x);
    if (val_is_na){
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] != p_x[i];
        if (!any_eq && eq){
          any_eq = true;
          out = Rf_protect(Rf_duplicate(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = Rf_protect(x); ++NP;
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!any_eq && eq){
          any_eq = true;
          out = Rf_protect(Rf_duplicate(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = Rf_protect(x); ++NP;
      }
    }
    break;
  }
  case STRSXP: {
    Rf_protect(value = Rf_coerceVector(value, STRSXP)); ++NP;
    Rf_protect(replace = Rf_coerceVector(replace, STRSXP)); ++NP;
    SEXP val = Rf_protect(Rf_asChar(value)); ++NP;
    SEXP repl = Rf_protect(Rf_asChar(replace)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = Rf_protect(Rf_duplicate(x)); ++NP;
      }
      if (eq) SET_STRING_ELT(out, i, repl);
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = Rf_protect(x); ++NP;
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    out = Rf_protect(Rf_allocVector(VECSXP, n)); ++NP;
    for (R_xlen_t i = 0; i < n; ++i){
      // Initialise each element
      SET_VECTOR_ELT(out, i, VECTOR_ELT(x, i));

      // Once we extract the vector it maybe needs protecting??
      SET_VECTOR_ELT(out, i, cpp_val_replace(VECTOR_ELT(out, i), value, replace, true));
    }
    // Copy attributes from x to y
    SEXP attrs = Rf_protect(Rf_coerceVector(ATTRIB(x), VECSXP)); ++NP;
    cpp_set_add_attributes(out, attrs, true);
    break;
  }
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(NP);
  return out;
}
// SEXP cpp_val_replace(SEXP x, SEXP value, SEXP replace, bool set){
//
//   //TO-DO: Make sure this can't work on ALTREP
//
//   int NP = 0;
//   R_xlen_t n = Rf_xlength(x);
//   if (Rf_length(value) != 1){
//     Rf_error("value must be a vector of length 1");
//   }
//   if (Rf_length(replace) != 1){
//     Rf_error("replace must be a vector of length 1");
//   }
//   if (Rf_isVectorList(x)){
//     Rf_error("x must not be a list");
//   }
//   bool val_is_na = cpp_any_na(value, true);
//   if (val_is_na && !cpp_any_na(x, true)){
//     Rf_unprotect(NP);
//     return x;
//   }
//   SEXP out;
//   switch ( TYPEOF(x) ){
//   case NILSXP: {
//     out = Rf_protect(R_NilValue);
//     ++NP;
//     break;
//   }
//   case LGLSXP:
//   case INTSXP: {
//     Rf_protect(value = Rf_coerceVector(value, INTSXP));
//     ++NP;
//     Rf_protect(replace = Rf_coerceVector(replace, INTSXP));
//     ++NP;
//     int val = Rf_asInteger(value);
//     int repl = Rf_asInteger(replace);
//     int *p_x = INTEGER(x);
//     out = Rf_protect(set ? x : Rf_duplicate(x));
//     ++NP;
//     int *p_out = INTEGER(out);
//     for (R_xlen_t i = 0; i < n; ++i) if (p_x[i] == val) p_out[i] = repl;
//     break;
//   }
//   case REALSXP: {
//     Rf_protect(value = Rf_coerceVector(value, REALSXP));
//     ++NP;
//     Rf_protect(replace = Rf_coerceVector(replace, REALSXP));
//     ++NP;
//     double val = Rf_asReal(value);
//     double repl = Rf_asReal(replace);
//     double *p_x = REAL(x);
//     out = Rf_protect(set ? x : Rf_duplicate(x));
//     ++NP;
//     double *p_out = REAL(out);
//     if (val_is_na){
//       for (R_xlen_t i = 0; i < n; ++i) if (p_x[i] != p_x[i]) p_out[i] = repl;
//     } else {
//       for (R_xlen_t i = 0; i < n; ++i) if (p_x[i] == val) p_out[i] = repl;
//     }
//     break;
//   }
//   case STRSXP: {
//     Rf_protect(value = Rf_coerceVector(value, STRSXP));
//     ++NP;
//     Rf_protect(replace = Rf_coerceVector(replace, STRSXP));
//     ++NP;
//     SEXP val = Rf_protect(Rf_asChar(value));
//     ++NP;
//     SEXP repl = Rf_protect(Rf_asChar(replace));
//     ++NP;
//     const SEXP *p_x = STRING_PTR_RO(x);
//     out = Rf_protect(set ? x : Rf_duplicate(x));
//     ++NP;
//     for (R_xlen_t i = 0; i < n; ++i) if (p_x[i] == val) SET_STRING_ELT(out, i, repl);
//     break;
//   }
//   default: {
//     Rf_unprotect(NP);
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
//   Rf_unprotect(NP);
//   return out;
// }
