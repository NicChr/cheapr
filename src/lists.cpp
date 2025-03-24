#include "cheapr.h"

R_xlen_t unnested_length(SEXP x){
  if (TYPEOF(x) != VECSXP){
    return Rf_xlength(x);
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t out = 0;
  for (R_xlen_t i = 0; i < n; ++i){
    out += (TYPEOF(p_x[i]) == VECSXP) ? unnested_length(p_x[i]) : Rf_xlength(p_x[i]);
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
  if (TYPEOF(x) != VECSXP){
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
SEXP cpp_new_list(SEXP size, SEXP default_value) {
  if (Rf_length(size) != 1){
   Rf_error("`size` must be a vector of length 1");
  }
  R_xlen_t out_size;
  if (TYPEOF(size) == INTSXP){
    out_size = Rf_asInteger(size);
  } else {
    out_size = Rf_asReal(size);
  }
  SEXP out = SHIELD(new_vec(VECSXP, out_size));
  if (!Rf_isNull(default_value)){
    for (R_xlen_t i = 0; i < out_size; ++i) {
      SET_VECTOR_ELT(out, i, default_value);
    }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP shallow_copy(SEXP x){
  if (TYPEOF(x) == VECSXP){
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

SEXP which_not_null(SEXP x){
  const SEXP *p_x = VECTOR_PTR_RO(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t n_null = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    n_null += (p_x[i] == R_NilValue);
  }
  R_xlen_t n_keep = n - n_null;
  R_xlen_t whichj = 0;
  R_xlen_t j = 0;

  // Which list elements should we keep?

  SEXP keep = SHIELD(new_vec(INTSXP, n_keep));
  int *p_keep = INTEGER(keep);
  while (whichj < n_keep){
    p_keep[whichj] = j + 1;
    whichj += (p_x[j++] != R_NilValue);
  }
  YIELD(1);
  return keep;
}

// From writing R extensions 5.9.7

SEXP get_list_element(SEXP list, const char *str){
  SEXP out = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);

  for (int i = 0; i < Rf_length(list); ++i){
    if (std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
      out = VECTOR_ELT(list, i);
      break;
    }
  }
    return out;
}

// Multi-assign named list elements

[[cpp11::register]]
SEXP cpp_list_assign(SEXP x, SEXP values){
  int NP = 0;

  SEXP names = SHIELD(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
  SEXP col_names = SHIELD(Rf_getAttrib(values, R_NamesSymbol)); ++NP;

  if (TYPEOF(x) != VECSXP){
    YIELD(NP);
    Rf_error("`x` must be a list in %s", __func__);
  }
  if (TYPEOF(values) != VECSXP || Rf_isNull(col_names)){
    YIELD(NP);
    Rf_error("`values` must be a named list in %s", __func__);
  }

  int n = Rf_length(x);
  int n_cols = Rf_length(values);

  if (Rf_isNull(names)){
    SHIELD(names = new_vec(STRSXP, n)); ++NP;
  }

  const SEXP *p_x = VECTOR_PTR_RO(x);
  const SEXP *p_names = STRING_PTR_RO(names);
  const SEXP *p_y = VECTOR_PTR_RO(values);
  const SEXP *p_col_names = STRING_PTR_RO(col_names);

  // We create an int vec to keep track of locations where to add col vecs

  SEXP add_locs = SHIELD(Rf_match(names, col_names, NA_INTEGER)); ++NP;
  int *p_add_locs = INTEGER(add_locs);
  int n_cols_to_add = na_count(add_locs, false);
  int out_size = n + n_cols_to_add;

  int loc;
  SEXP out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
  SEXP out_names = SHIELD(new_vec(STRSXP, out_size)); ++NP;

  // Initialise out
  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, p_x[i]);
    SET_STRING_ELT(out_names, i, p_names[i]);
  }

  bool any_null = false;

  for (int j = 0; j < n_cols; ++j){
    loc = p_add_locs[j];
    any_null = any_null || Rf_isNull(p_y[j]);
    if (is_na_int(loc)){
      SET_VECTOR_ELT(out, n, p_y[j]);
      SET_STRING_ELT(out_names, n, p_col_names[j]);
      ++n;
    } else {
      --loc;
      SET_VECTOR_ELT(out, loc, p_y[j]);
      SET_STRING_ELT(out_names, loc, p_col_names[j]);
    }
  }
  if (any_null){
    SEXP keep = SHIELD(which_not_null(out)); ++NP;
    SHIELD(out = sset_vec(out, keep, false)); ++NP;
    SHIELD(out_names = sset_vec(out_names, keep, false)); ++NP;
  }
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  YIELD(NP);
  return out;
}

// Work-in-progress
// SEXP cpp_list_assign2(SEXP x, SEXP values, SEXP locs){
//   int NP = 0;
//
//   SEXP names = SHIELD(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
//   SEXP col_names = SHIELD(Rf_getAttrib(values, R_NamesSymbol)); ++NP;
//
//   if (TYPEOF(x) != VECSXP){
//     YIELD(NP);
//     Rf_error("`x` must be a list in %s", __func__);
//   }
//   // if (TYPEOF(values) != VECSXP || (Rf_isNull(col_names) && !Rf_isNull(locs))){
//   //   YIELD(NP);
//   //   Rf_error("`values` must be a named list in %s", __func__);
//   // }
//
//   int n = Rf_length(x);
//   int n_cols = Rf_length(values);
//
//   if (Rf_isNull(names)){
//     SHIELD(names = new_vec(STRSXP, n)); ++NP;
//   }
//
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//   const SEXP *p_names = STRING_PTR_RO(names);
//   const SEXP *p_y = VECTOR_PTR_RO(values);
//   const SEXP *p_col_names = STRING_PTR_RO(col_names);
//
//   // We create an int vec to keep track of locations where to add col vecs
//   SEXP add_locs;
//   if (Rf_isNull(locs)){
//     add_locs = SHIELD(Rf_match(names, col_names, NA_INTEGER)); ++NP;
//   } else {
//     SHIELD(locs = coerce_vec(locs, INTSXP)); ++NP;
//     if (Rf_length(locs) != Rf_length(values)){
//       YIELD(NP);
//       Rf_error("`length(locs)` must match `length(values)` when `!is.null(locs)`");
//     }
//     add_locs = SHIELD(Rf_match(names, col_names, NA_INTEGER)); ++NP;
//     SEXP temp1 = SHIELD(cpp_seq_len(n)); ++NP;
//     SEXP temp2 = SHIELD(Rf_match(temp1, add_locs, NA_INTEGER)); ++NP;
//     SHIELD(col_names = sset_vec(names, temp2, true)); ++NP;
//   }
//
//   int *p_add_locs = INTEGER(add_locs);
//   int n_cols_to_add = na_count(add_locs, false);
//   int out_size = n + n_cols_to_add;
//
//   int loc;
//   SEXP out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
//   SEXP out_names = SHIELD(new_vec(STRSXP, out_size)); ++NP;
//
//   // Initialise out
//   for (int i = 0; i < n; ++i){
//     SET_VECTOR_ELT(out, i, p_x[i]);
//     SET_STRING_ELT(out_names, i, p_names[i]);
//   }
//
//   bool any_null = false;
//
//   for (int j = 0; j < n_cols; ++j){
//     loc = p_add_locs[j];
//     any_null = any_null || Rf_isNull(p_y[j]);
//     if (is_na_int(loc)){
//       SET_VECTOR_ELT(out, n, p_y[j]);
//       SET_STRING_ELT(out_names, n, p_col_names[j]);
//       ++n;
//     } else {
//       --loc;
//       SET_VECTOR_ELT(out, loc, p_y[j]);
//       SET_STRING_ELT(out_names, loc, p_col_names[j]);
//     }
//   }
//   if (any_null){
//     SEXP keep = SHIELD(which_not_null(out)); ++NP;
//     SHIELD(out = sset_vec(out, keep, false)); ++NP;
//     SHIELD(out_names = sset_vec(out_names, keep, false)); ++NP;
//   }
//   Rf_setAttrib(out, R_NamesSymbol, out_names);
//   YIELD(NP);
//   return out;
// }

// SEXP df_assign_cols(SEXP x, SEXP cols){
//   if (!is_df(x)){
//     Rf_error("`x` must be a `data.frame` in %s", __func__);
//   }
//
//   int nrows = df_nrow(x);
//   SEXP r_nrows = SHIELD(Rf_ScalarInteger(nrows));
//
//   SEXP out = SHIELD(cpp_list_assign(x, cols));
//   SHIELD(out = cpp_recycle(out, r_nrows));
//
//   Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(nrows));
//   Rf_classgets(out, Rf_mkString("data.frame"));
//   YIELD(3);
//   return out;
// }

// Data-frames

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
  SEXP row_names = SHIELD(create_df_row_names(N)); ++NP;

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

// Can use in the future once I figure out how to re-write named_list() in C++

// SEXP cpp_new_df(SEXP x, SEXP nrows, bool recycle, bool name_repair){
//
//   int NP = 0;
//
//   SEXP out = SHIELD(shallow_copy(x)); ++NP;
//
//   if (recycle){
//     SHIELD(out = cpp_recycle(out, nrows)); ++NP;
//   }
//
//   SEXP row_names;
//   if (Rf_isNull(nrows)){
//     if (Rf_length(out) == 0){
//       row_names = SHIELD(new_vec(INTSXP, 0)); ++NP;
//     } else {
//       row_names = SHIELD(create_df_row_names(vec_length(VECTOR_ELT(out, 0)))); ++NP;
//     }
//   } else {
//     row_names = SHIELD(create_df_row_names(Rf_asInteger(nrows))); ++NP;
//   }
//
//   SEXP out_names = SHIELD(Rf_getAttrib(out, R_NamesSymbol)); ++NP;
//   SHIELD(out_names = coerce_vec(out_names, STRSXP)); ++NP;
//
//   if (name_repair){
//     SEXP sep = SHIELD(Rf_mkString("...")); ++NP;
//     SHIELD(out_names = cpp_name_repair(out_names, sep)); ++NP;
//   }
//   Rf_setAttrib(out, R_NamesSymbol, out_names);
//   Rf_setAttrib(out, R_RowNamesSymbol, row_names);
//   Rf_classgets(out, Rf_mkString("data.frame"));
//   YIELD(NP);
//   return out;
// }

// Multi-assign recycled variables to data frame

[[cpp11::register]]
SEXP cpp_df_assign_cols(SEXP x, SEXP cols){
  int NP = 0;

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame` in %s", __func__);
  }

  SEXP names = SHIELD(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
  SEXP col_names = SHIELD(Rf_getAttrib(cols, R_NamesSymbol)); ++NP;

  if (TYPEOF(cols) != VECSXP || Rf_isNull(col_names)){
    Rf_error("`cols` must be a named list in %s", __func__);
  }

  const SEXP *p_x = VECTOR_PTR_RO(x);
  const SEXP *p_names = STRING_PTR_RO(names);
  const SEXP *p_cols = VECTOR_PTR_RO(cols);
  const SEXP *p_col_names = STRING_PTR_RO(col_names);

  int n = Rf_length(x);
  int n_cols = Rf_length(cols);
  int n_rows = df_nrow(x);

  // We create an int vec to keep track of locations where to add col vecs

  SEXP add_locs = SHIELD(Rf_match(names, col_names, NA_INTEGER)); ++NP;
  int *p_add_locs = INTEGER(add_locs);
  int n_cols_to_add = na_count(add_locs, false);
  int out_size = n + n_cols_to_add;

  int loc;
  SEXP out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
  SEXP out_names = SHIELD(new_vec(STRSXP, out_size)); ++NP;

  // Initialise out
  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, p_x[i]);
    SET_STRING_ELT(out_names, i, p_names[i]);
  }

  bool any_null = false;

  SEXP vec;
  PROTECT_INDEX vec_idx;

  R_ProtectWithIndex(vec = R_NilValue, &vec_idx); ++NP;

  for (int j = 0; j < n_cols; ++j){
    loc = p_add_locs[j];
    R_Reprotect(vec = p_cols[j], vec_idx);
    any_null = any_null || Rf_isNull(vec);
    if (is_na_int(loc)){
      SET_VECTOR_ELT(out, n, cpp_rep_len(vec, n_rows));
      SET_STRING_ELT(out_names, n, p_col_names[j]);
      ++n;
    } else {
      --loc;
      SET_VECTOR_ELT(out, loc, cpp_rep_len(vec, n_rows));
      SET_STRING_ELT(out_names, loc, p_col_names[j]);
    }
  }
  if (any_null){
    SEXP keep = SHIELD(which_not_null(out)); ++NP;
    SHIELD(out = sset_vec(out, keep, false)); ++NP;
    SHIELD(out_names = sset_vec(out_names, keep, false)); ++NP;
  }
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(n_rows));
  Rf_classgets(out, Rf_mkString("data.frame"));
  YIELD(NP);
  return out;
}


[[cpp11::register]]
SEXP cpp_df_reconstruct(SEXP data, SEXP from, bool keep_attrs){
  if (!is_df(data)){
    Rf_error("`data` must be a `data.frame`");
  }
  if (!is_df(from)){
    Rf_error("`from` must be a `data.frame`");
  }

  // Create shallow copies so that we can manipulate attributes freely

  SEXP target = SHIELD(shallow_copy(data));
  SEXP source = SHIELD(shallow_copy(from));

  // The below strips any leftover attributes from `data`,
  // I wonder if these should be kept

  cpp_set_rm_attributes(target);

  if (keep_attrs){
    Rf_setAttrib(source, R_NamesSymbol, R_NilValue);
    Rf_setAttrib(source, R_ClassSymbol, R_NilValue);
    Rf_setAttrib(source, R_RowNamesSymbol, R_NilValue);
    SHALLOW_DUPLICATE_ATTRIB(target, source);
  }

  // Re-add original data attributes as these cannot be changed
  Rf_setAttrib(target, R_NamesSymbol, Rf_getAttrib(data, R_NamesSymbol));
  Rf_setAttrib(target, R_ClassSymbol, Rf_getAttrib(from, R_ClassSymbol));
  Rf_setAttrib(target, R_RowNamesSymbol, create_df_row_names(df_nrow(data)));
  YIELD(2);
  return target;
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
