#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>

int int_div(int x, int y){
  return x / y;
}

[[cpp11::register]]
R_xlen_t cpp_vec_length(SEXP x){
  if (Rf_isFrame(x)){
    return cpp_df_nrow(x);
    // Is x a list?
  } else if (Rf_isVectorList(x)){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return cpp_vec_length(VECTOR_ELT(x, 0));
    } else if (Rf_inherits(x, "POSIXlt")){
      return Rf_xlength(VECTOR_ELT(x, 0));
    } else if (Rf_isObject(x)){
      return Rf_asReal(cpp11::package("base")["length"](x));
    } else {
      return Rf_xlength(x);
    }
    // Catch-all
  } else {
    return Rf_xlength(x);
  }
}
int num_cores(){
  int out = Rf_asInteger(Rf_GetOption1(Rf_installChar(Rf_mkChar("cheapr.cores"))));
  if (out >= 1){
    return out;
  } else {
    return 1;
  }
}

SEXP xlen_to_r(R_xlen_t x){
  return x > integer_max_ ? Rf_ScalarReal(x) : Rf_ScalarInteger(x);
}

R_xlen_t cpp_df_nrow(SEXP x){
  return Rf_xlength(Rf_getAttrib(x, R_RowNamesSymbol));
}

R_xlen_t cpp_unnested_length(SEXP x){
  if (!Rf_isVectorList(x)){
    return Rf_xlength(x);
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t out = 0;
  for (R_xlen_t i = 0; i < n; ++i){
   out += Rf_isVectorList(p_x[i]) ? cpp_unnested_length(p_x[i]) : Rf_xlength(p_x[i]);
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_r_unnested_length(SEXP x){
  return xlen_to_r(cpp_unnested_length(x));
}

[[cpp11::register]]
SEXP cpp_lengths(SEXP x) {
  R_xlen_t n = Rf_xlength(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
  int *p_out = INTEGER(out);
  if (!Rf_isVectorList(x)){
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = 1;
    }
    Rf_unprotect(1);
    return out;
  } else {
    const SEXP* p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = cpp_vec_length(p_x[i]);
    }
    Rf_unprotect(1);
    return out;
  }
}

[[cpp11::register]]
SEXP cpp_new_list(R_xlen_t size, SEXP default_value) {
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, size));
  if (!Rf_isNull(default_value)){
    for (R_xlen_t i = 0; i < size; ++i) {
      SET_VECTOR_ELT(out, i, default_value);
    }
  }
  Rf_unprotect(1);
  return out;
}

// Remove NULL elements from list

[[cpp11::register]]
SEXP cpp_list_rm_null(SEXP l) {
  Rf_protect(l = Rf_coerceVector(l, VECSXP));
  const SEXP *p_l = VECTOR_PTR_RO(l);
  int n = Rf_length(l);
  int n_null = 0;
  for (int i = 0; i < n; ++i) {
    n_null += (p_l[i] == R_NilValue);
  }
  int n_keep = n - n_null;
  int whichj = 0;
  int j = 0;
  SEXP keep = Rf_protect(Rf_allocVector(INTSXP, n_keep));
  int *p_keep = INTEGER(keep);
  while (whichj < n_keep){
    p_keep[whichj] = j + 1;
    whichj += (p_l[j] != R_NilValue);
    ++j;
  }
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, n_keep));
  SEXP names = Rf_protect(Rf_duplicate(Rf_getAttrib(l, R_NamesSymbol)));
  bool has_names = !Rf_isNull(names);
  if (has_names){
    SEXP *p_names = STRING_PTR(names);
    SEXP out_names = Rf_protect(Rf_allocVector(STRSXP, n_keep));
    for (int k = 0; k < n_keep; ++k) {
      SET_STRING_ELT(out_names, k, p_names[p_keep[k] - 1]);
      SET_VECTOR_ELT(out, k, p_l[p_keep[k] - 1]);
    }
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    Rf_unprotect(5);
    return out;
  } else {
    for (int k = 0; k < n_keep; ++k) {
      SET_VECTOR_ELT(out, k, p_l[p_keep[k] - 1]);
    }
    Rf_unprotect(4);
    return out;
  }
}

[[cpp11::register]]
SEXP cpp_list_as_df(SEXP x) {
  SEXP out = Rf_protect(cpp_list_rm_null(x));
  int N; // Number of rows
  if (Rf_length(out) == 0){
    N = 0;
  } else {
    N = cpp_vec_length(VECTOR_ELT(out, 0));
  }
  SEXP df_str = Rf_protect(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(df_str, 0, Rf_mkChar("data.frame"));
  if (N > 0){
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 2));
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -N;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    Rf_classgets(out, df_str);
    Rf_unprotect(3);
    return out;
  } else {
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    Rf_classgets(out, df_str);
    Rf_unprotect(3);
    return out;
  }
}

// Remove attributes in-place

[[cpp11::register]]
SEXP cpp_set_rm_attributes(SEXP x){
  SEXP attrs = Rf_protect(cpp11::package("base")["attributes"](x));
  SEXP names = Rf_protect(Rf_getAttrib(attrs, R_NamesSymbol));
  int n = Rf_length(attrs);
  for (int i = 0; i < n; ++i){
    SEXP attrib_nm = Rf_protect(Rf_install(CHAR(STRING_ELT(names, i))));
    Rf_setAttrib(x, attrib_nm, R_NilValue);
  }
  Rf_unprotect(n + 2);
  return x;
}

// Copy specified attributes (character vector of names)
// from source to target (by reference)
// Use with extreme care as it modifies target in-place
// If you use it, make absolutely sure that target is not pointed to by other
// objects as it will modify the attributes of those objects too

[[cpp11::register]]
SEXP cpp_set_copy_attributes(SEXP target, SEXP source, SEXP attrs){
  SEXP *p_attrs = STRING_PTR(attrs);
  int n_attrs = Rf_length(attrs);
  for (int i = 0; i < n_attrs; ++i){
    SEXP attrib_nm = Rf_protect(Rf_install(CHAR(p_attrs[i])));
    Rf_setAttrib(target, attrib_nm, Rf_getAttrib(source, attrib_nm));
  }
  Rf_unprotect(n_attrs);
  return target;
}

// SEXP cpp_unlist(SEXP x, SEXP ptype) {
//   if (!Rf_isVectorList(x)){
//     Rf_error("x must be a list");
//   }
//   int n_protections = 0;
//   R_xlen_t n = Rf_xlength(x);
//   R_xlen_t N = cpp_unnested_length(x);
//   R_xlen_t m;
//   R_xlen_t k = 0;
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//   switch ( TYPEOF(ptype) ){
//   case LGLSXP: {
//     ++n_protections;
//     SEXP out = Rf_protect(Rf_allocVector(LGLSXP, N));
//     int *p_out = LOGICAL(out);
//     for (R_xlen_t i = 0; i < n; ++i) {
//       m = Rf_xlength(p_x[i]);
//       int *p_xj = LOGICAL(p_x[i]);
//       for (R_xlen_t j = 0; j < m; ++j) {
//         p_out[k] = p_xj[j];
//         ++k;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case INTSXP: {
//     ++n_protections;
//     SEXP out = Rf_protect(Rf_allocVector(INTSXP, N));
//     int *p_out = INTEGER(out);
//     for (R_xlen_t i = 0; i < n; ++i) {
//       m = Rf_xlength(p_x[i]);
//       int *p_xj = INTEGER(p_x[i]);
//       for (R_xlen_t j = 0; j < m; ++j) {
//         p_out[k] = p_xj[j];
//         ++k;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case REALSXP: {
//     ++n_protections;
//     SEXP out = Rf_protect(Rf_allocVector(REALSXP, N));
//     double *p_out = REAL(out);
//     for (R_xlen_t i = 0; i < n; ++i) {
//       m = Rf_xlength(p_x[i]);
//       double *p_xj = REAL(p_x[i]);
//       for (R_xlen_t j = 0; j < m; ++j) {
//         p_out[k] = p_xj[j];
//         ++k;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case STRSXP: {
//     ++n_protections;
//     SEXP out = Rf_protect(Rf_allocVector(STRSXP, N));
//     for (R_xlen_t i = 0; i < n; ++i) {
//       m = Rf_xlength(p_x[i]);
//       SEXP *p_xj = STRING_PTR(p_x[i]);
//       for (R_xlen_t j = 0; j < m; ++j) {
//         SET_STRING_ELT(out, k, p_xj[j]);
//         ++k;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(ptype)));
//   }
//   }
// }

// Potentially useful for rolling calculations
// Computes the rolling number of true values in a given
// series of consecutive true values

// SEXP cpp_run_id(SEXP x, bool left_to_right){
//   if (!Rf_isLogical(x)){
//     Rf_error("x must be a logical vector");
//   }
//   R_xlen_t count = 0;
//   R_xlen_t n = Rf_xlength(x);
//   SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
//   int *p_out = INTEGER(out);
//   int *p_x = LOGICAL(x);
//   if (left_to_right){
//     for (R_xlen_t i = 0; i < n; ++i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   } else {
//     for (R_xlen_t i = n - 1; i >= 0; --i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }
