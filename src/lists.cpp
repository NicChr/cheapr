#include "cheapr.h"

R_xlen_t unnested_length(SEXP x){
  if (!Rf_isVectorList(x)){
    return Rf_xlength(x);
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t out = 0;
  for (R_xlen_t i = 0; i < n; ++i){
    out += Rf_isVectorList(p_x[i]) ? unnested_length(p_x[i]) : Rf_xlength(p_x[i]);
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_unnested_length(SEXP x){
  return xlen_to_r(unnested_length(x));
}

[[cpp11::register]]
SEXP cpp_lengths(SEXP x, bool names) {
  R_xlen_t n = Rf_xlength(x);
  SEXP out = SHIELD(new_vec(INTSXP, n));
  int *p_out = INTEGER(out);
  if (!Rf_isVectorList(x)){
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = 1;
    }
  } else {
    const SEXP* p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = vec_length(p_x[i]);
    }
  }
  if (names){
    cpp_copy_names(x, out, false);
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_new_list(R_xlen_t size, SEXP default_value) {
  SEXP out = SHIELD(new_vec(VECSXP, size));
  if (!Rf_isNull(default_value)){
    for (R_xlen_t i = 0; i < size; ++i) {
      SET_VECTOR_ELT(out, i, default_value);
    }
  }
  YIELD(1);
  return out;
}

SEXP shallow_copy(SEXP x){
  if (Rf_isVectorList(x)){
    R_xlen_t n = Rf_xlength(x);
    SEXP out = SHIELD(new_vec(VECSXP,  n));
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, p_x[i]);
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    YIELD(1);
    return out;
  } else {
    return x;
  }
}

// Remove NULL elements from list

[[cpp11::register]]
SEXP cpp_drop_null(SEXP l, bool always_shallow_copy) {
  SHIELD(l = coerce_vec(l, VECSXP));
  const SEXP *p_l = VECTOR_PTR_RO(l);
  int n = Rf_length(l);
  int n_null = 0;
  for (int i = 0; i < n; ++i) {
    n_null += (p_l[i] == R_NilValue);
  }
  if (n_null == 0 && !always_shallow_copy){
    YIELD(1);
    return l;
  }
  int n_keep = n - n_null;
  int whichj = 0;
  int j = 0;

  // Which list elements should we keep?

  SEXP keep = SHIELD(new_vec(INTSXP, n_keep));
  int *p_keep = INTEGER(keep);
  while (whichj < n_keep){
    p_keep[whichj] = j;
    whichj += (p_l[j++] != R_NilValue);
  }

  // Subset on both the list and names of the list

  SEXP out = SHIELD(new_vec(VECSXP, n_keep));
  SEXP names = SHIELD(Rf_getAttrib(l, R_NamesSymbol));
  bool has_names = !Rf_isNull(names);
  if (has_names){
    const SEXP *p_names = STRING_PTR_RO(names);
    SEXP out_names = SHIELD(new_vec(STRSXP, n_keep));
    for (int k = 0; k < n_keep; ++k) {
      SET_STRING_ELT(out_names, k, p_names[p_keep[k]]);
      SET_VECTOR_ELT(out, k, p_l[p_keep[k]]);
    }
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    YIELD(5);
    return out;
  } else {
    for (int k = 0; k < n_keep; ++k) {
      SET_VECTOR_ELT(out, k, p_l[p_keep[k]]);
    }
    YIELD(4);
    return out;
  }
}

// From writing R extensions 5.9.7

SEXP get_list_element(SEXP list, const char *str){
  SEXP out = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);

  for (int i = 0; i < Rf_length(list); i++){
    if (std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
      out = VECTOR_ELT(list, i);
      break;
    }
  }
    return out;
}


[[cpp11::register]]
SEXP cpp_list_as_df(SEXP x) {
  int N; // Number of rows
  int NP = 0; // Number of protects
  SEXP out = SHIELD(cpp_drop_null(x, true)); ++NP;
  int n_items = Rf_length(out);
  if (is_df(x)){
    N = df_nrow(x);
  } else if (n_items == 0){
    N = 0;
  } else {
    N = vec_length(VECTOR_ELT(out, 0));
  }

  SEXP df_str = SHIELD(Rf_mkString("data.frame")); ++NP;
  SEXP row_names = create_df_row_names(N);

  // If no names then add names
  if (Rf_isNull(Rf_getAttrib(out, R_NamesSymbol))){
    SEXP out_names = SHIELD(new_vec(STRSXP, n_items)); ++NP;
    Rf_setAttrib(out, R_NamesSymbol, out_names);
  }
  Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  Rf_classgets(out, df_str);
  YIELD(NP);
  return out;
}

// void cpp_check_nested_lengths(SEXP x, SEXP y){
//   R_xlen_t n1 = Rf_xlength(x);
//   R_xlen_t n2 = Rf_xlength(y);
//   if (n1 != n2){
//     Rf_error("x and y must have the same length");
//   }
//   if (Rf_isVectorList(x) && Rf_isVectorList(y)){
//     R_xlen_t n3, n4;
//     const SEXP *p_x = VECTOR_PTR_RO(x);
//     const SEXP *p_y = VECTOR_PTR_RO(y);
//
//     for (R_xlen_t i = 0; i < n1; ++i){
//       bool xlist = Rf_isVectorList(p_x[i]);
//       bool ylist = Rf_isVectorList(p_y[i]);
//       int both_lists = xlist + ylist;
//       if (both_lists == 1){
//         Rf_error("x and y must have identical nested lengths");
//       } else if (both_lists == 2){
//         // Recurse back through the same function at this point
//         cpp_check_nested_lengths(p_x[i], p_y[i]);
//       } else {
//         n3 = Rf_xlength(p_x[i]);
//         n4 = Rf_xlength(p_y[i]);
//         if (n3 != n4){
//           Rf_error("x and y must have identical nested lengths");
//         }
//       }
//     }
//   } else if (!(!Rf_isVectorList(x) && !Rf_isVectorList(y))){
//     Rf_error("x and y must either be both lists or both not lists");
//   }
// }

// #define cheapr_cast_temp(x, y) cpp11::function cpp11::package("cheapr")["cheapr_cast"];

// SEXP cpp_cast_common(SEXP x, SEXP y){
//   // All length checks will have been done above..
//   // Maybe inefficient but makes things simpler
//   R_xlen_t n = Rf_xlength(x);
//   cpp11::function cheapr_cast = cpp11::package("cheapr")["cheapr_cast"];
//   int n_prot = 0;
//   SEXP out = SHIELD(new_vec(VECSXP, 2));
//   ++n_prot;
//   if (Rf_isVectorList(x) && Rf_isVectorList(y)){
//     // SEXP a = SHIELD(cpp_shallow_copy(x));
//     SEXP a = SHIELD(Rf_shallow_duplicate(x));
//     ++n_prot;
//     SEXP b = SHIELD(Rf_shallow_duplicate(y));
//     // SEXP b = SHIELD(cpp_shallow_copy(y));
//     ++n_prot;
//     const SEXP *p_x = VECTOR_PTR_RO(a);
//     const SEXP *p_y = VECTOR_PTR_RO(b);
//
//     for (R_xlen_t i = 0; i < n; ++i){
//       bool xlist = Rf_isVectorList(p_x[i]);
//       bool ylist = Rf_isVectorList(p_y[i]);
//       if (xlist && ylist){
//         // Recurse back through the same function at this point
//         SEXP temp = SHIELD(cpp_cast_common(p_x[i], p_y[i]));
//         ++n_prot;
//         SET_VECTOR_ELT(a, i, VECTOR_ELT(temp, 0));
//         SET_VECTOR_ELT(b, i, VECTOR_ELT(temp, 1));
//       } else {
//         SET_VECTOR_ELT(a, i, cheapr_cast(p_x[i], p_y[i]));
//         SET_VECTOR_ELT(b, i, cheapr_cast(p_y[i], p_x[i]));
//       }
//     }
//     SET_VECTOR_ELT(out, 0, a);
//     SET_VECTOR_ELT(out, 1, b);
//   } else {
//     SET_VECTOR_ELT(out, 0, cheapr_cast(x, y));
//     SET_VECTOR_ELT(out, 1, cheapr_cast(y, x));
//   }
//   YIELD(n_prot);
//   return out;
// }
