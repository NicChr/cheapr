#include "cheapr.h"

// From R you pass `list(...)` and `.args` to below function
// which then decides which arguments to use

[[cpp11::register]]
SEXP cpp_list_args(SEXP args1, SEXP args2){
  bool use_dots = Rf_length(args1) != 0;
  bool use_list = args2 != R_NilValue;

  if (use_dots && use_list){
    Rf_error("Please supply either `...` or `.args` in %s", __func__);
  }
  if (use_list && (TYPEOF(args2) != VECSXP || Rf_isObject(args2))){
    Rf_error("`.args` must be a plain list in %s", __func__);
  }
  return use_list ? args2 : args1;
}

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
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = vec_length(p_x[i]);
    }
  }
  if (names){
    set_names(out, get_names(x));
  }
  YIELD(1);
  return out;
}

SEXP new_list(R_xlen_t length, SEXP default_value){
  SEXP out = SHIELD(new_vec(VECSXP, length));
  if (!is_null(default_value)){
    for (R_xlen_t i = 0; i < length; ++i) {
      SET_VECTOR_ELT(out, i, default_value);
    }
  }
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_new_list(SEXP size, SEXP default_value){
  if (Rf_length(size) != 1){
   Rf_error("`size` must be a vector of length 1");
  }
  R_xlen_t out_size;
  if (TYPEOF(size) == INTSXP){
    out_size = INTEGER(size)[0];
  } else {
    out_size = REAL(size)[0];
  }
  return new_list(out_size, default_value);
}

// Remove NULL elements from list

[[cpp11::register]]
SEXP cpp_drop_null(SEXP l, bool always_shallow_copy) {
  const SEXP *p_l = VECTOR_PTR_RO(l);
  int n = Rf_length(l);
  int n_null = 0;
  for (int i = 0; i < n; ++i) {
    n_null += is_null(p_l[i]);
  }
  if (n_null == 0 && !always_shallow_copy){
    return l;
  }
  int n_keep = n - n_null;
  int whichj = 0;
  int j = 0;

  // Which list elements should we keep?

  SEXP keep = SHIELD(new_vec(INTSXP, n_keep));
  int* RESTRICT p_keep = INTEGER(keep);
  while (whichj < n_keep){
    p_keep[whichj] = j;
    whichj += !is_null(p_l[j++]);
  }

  // Subset on both the list and names of the list

  SEXP out = SHIELD(new_vec(VECSXP, n_keep));
  SEXP names = get_names(l);
  bool has_names = !is_null(names);
  if (has_names){
    const SEXP *p_names = STRING_PTR_RO(names);
    SEXP out_names = SHIELD(new_vec(STRSXP, n_keep));
    for (int k = 0; k < n_keep; ++k) {
      SET_STRING_ELT(out_names, k, p_names[p_keep[k]]);
      SET_VECTOR_ELT(out, k, p_l[p_keep[k]]);
    }
    set_names(out, out_names);
    YIELD(3);
    return out;
  } else {
    for (int k = 0; k < n_keep; ++k) {
      SET_VECTOR_ELT(out, k, p_l[p_keep[k]]);
    }
    YIELD(2);
    return out;
  }
}

SEXP which_not_null(SEXP x){
  const SEXP *p_x = VECTOR_PTR_RO(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t n_null = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    n_null += is_null(p_x[i]);
  }
  R_xlen_t n_keep = n - n_null;
  R_xlen_t whichj = 0;
  R_xlen_t j = 0;

  // Which list elements should we keep?

  SEXP keep = SHIELD(new_vec(INTSXP, n_keep));
  int* RESTRICT p_keep = INTEGER(keep);
  while (whichj < n_keep){
    p_keep[whichj] = j + 1;
    whichj += !is_null(p_x[j++]);
  }
  YIELD(1);
  return keep;
}

SEXP get_list_element(SEXP list, SEXP str){
  SEXP out = R_NilValue, names = get_names(list);

  for (int i = 0; i < Rf_length(list); ++i){
    if (STRING_ELT(names, i) == str){
      out = VECTOR_ELT(list, i);
      break;
    }
  }
  return out;
}


// A cpp11 pushback() method for growing a list
// It should in theory be efficient but because of
// cpp11's copy-on-modify, I can't figure out a way to call this function
// in another function repeatedly without copying every time
// To be efficient one must use push.back manually inside a loop in the same
// function

// cpp11::list add_list_element(cpp11::writable::list x, SEXP i, SEXP value){
//
//   cpp11::writable::strings names;
//
//   if (is_null(get_names(x))){
//     names = cpp11::writable::strings(x.size());
//     x.names() = names;
//   } else {
//     names = cpp11::as_sexp(x.names());
//   }
//
//   int k, loc;
//
//   if (TYPEOF(i) == INTSXP){
//     cpp11::integers locs = i;
//     loc = locs[0] - 1;
//     if (loc < x.size()){
//       x[loc] = value;
//     } else {
//       if (loc > x.size()){
//         cpp11::writable::list temp(2);
//         cpp11::writable::list join(loc - x.size());
//         temp[0] = x;
//         temp[1] = join;
//         x = list_c(temp);
//       }
//       x.push_back(value);
//     }
//   }
//   else if (TYPEOF(i) == REALSXP){
//     cpp11::doubles locs = i;
//     loc = locs[0] - 1;
//     if (loc < x.size()){
//       x[loc] = value;
//     } else {
//       if (loc > x.size()){
//         cpp11::writable::list temp(2);
//         cpp11::writable::list join(loc - x.size());
//         temp[0] = x;
//         temp[1] = join;
//         x = list_c(temp);
//       }
//       x.push_back(value);
//     }
//   } else {
//     k = 0;
//     cpp11::strings name = i;
//     for (int i = 0; i < x.size(); ++i, ++k){
//       if (names[i] == name[0]){
//         x[i] = value;
//         break;
//       }
//     }
//     if (k == x.size()){
//       x.push_back(value);
//       names.push_back(name[0]);
//       x.names() = names;
//     }
//   }
//   return x;
// }

// Multi-assign named list elements

[[cpp11::register]]
SEXP cpp_list_assign(SEXP x, SEXP values){
  int NP = 0;
  int n = Rf_length(x);
  int n_cols = Rf_length(values);


  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list in %s", __func__);
  }
  if (TYPEOF(values) != VECSXP){
    Rf_error("`values` must be a named list in %s", __func__);
  }

  SEXP names = get_names(x);
  SEXP col_names = get_names(values);

  if (is_null(names)){
    SHIELD(names = new_vec(STRSXP, n)); ++NP;
  }

  bool empty_value_names = is_null(col_names);

  if (empty_value_names){
    SHIELD(col_names = new_vec(STRSXP, n_cols)); ++NP;
  }

  const SEXP *p_x = VECTOR_PTR_RO(x);
  const SEXP *p_names = STRING_PTR_RO(names);
  const SEXP *p_y = VECTOR_PTR_RO(values);
  const SEXP *p_col_names = STRING_PTR_RO(col_names);

  // We create an int vec to keep track of locations where to add col vecs

  SEXP add_locs;
  int n_cols_to_add;

  // If values is an unnamed list then we can simply append the values
  if (empty_value_names){
    add_locs = SHIELD(new_vec(INTSXP, 0)); ++NP;
    SHIELD(add_locs = cpp_rep_len(add_locs, n_cols)); ++NP;
    n_cols_to_add = n_cols;
  } else {
    add_locs = SHIELD(Rf_match(names, col_names, NA_INTEGER)); ++NP;
    n_cols_to_add = na_count(add_locs, false);
  }
  int* RESTRICT p_add_locs = INTEGER(add_locs);
  int out_size = n + n_cols_to_add;

  int loc;
  SEXP out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
  SEXP out_names = SHIELD(new_vec(STRSXP, out_size)); ++NP;

  // Initialise out
  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, p_x[i]);
    SET_STRING_ELT(out_names, i, p_names[i]);
  }

  cpp11::writable::integers null_locs;

  for (int j = 0; j < n_cols; ++j){
    loc = p_add_locs[j];
    if (is_null(p_y[j])){
      null_locs.push_back(loc);
    }
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
  if (null_locs.size() != 0){
    int n_null = null_locs.size();
    int *p_null_locs = INTEGER(null_locs);
    for (int i = 0; i < n_null; ++i){
      p_null_locs[i] = -p_null_locs[i];
    }
    SEXP keep = SHIELD(exclude_locs(null_locs, out_size)); ++NP;
    SHIELD(out = sset_vec(out, keep, false)); ++NP;
    SHIELD(out_names = sset_vec(out_names, keep, false)); ++NP;
  }
  set_names(out, out_names);
  YIELD(NP);
  return out;
}

// Multi-assign named values using cpp11
// Not as fast as cpp_list_assign though

// cpp11::writable::list cpp_list_assign2(cpp11::writable::list x, cpp11::list values){
//   using namespace cpp11;
//
//   writable::strings names = as_sexp(x.names());
//   strings col_names = as_sexp(values.names());
//
//   if (TYPEOF(x) != VECSXP){
//     stop("`x` must be a list in %s", __func__);
//   }
//   if (TYPEOF(values) != VECSXP || is_null(col_names)){
//     stop("`x` must be a list in %s", __func__);
//   }
//
//   int n = x.size();
//   int n_cols = values.size();
//
//   if (is_null(names)){
//     names = writable::strings(n);
//   }
//
//   integers add_locs = Rf_match(names, col_names, NA_INTEGER);
//   writable::integers null_locs;
//
//   int loc;
//   for (int j = 0; j < n_cols; ++j){
//     loc = add_locs[j];
//     if (is_na_int(loc) && values[j] != R_NilValue){
//       x.push_back(values[j]);
//       names.push_back(col_names[j]);
//     } else {
//       if (values[j] == R_NilValue){
//         null_locs.push_back(loc);
//       } else {
//         --loc;
//         x[loc] = values[j];
//         names[loc] = col_names[j];
//       }
//     }
//   }
//   if (null_locs.size() > 0){
//     cpp_set_change_sign(null_locs);
//     integers not_null_locs = exclude_locs(null_locs, x.size());
//     x = sset_vec(x, not_null_locs, false);
//     names = sset_vec(names, not_null_locs, false);
//   }
//   x.names() = names;
//   return x;
// }


// Work-in-progress
// SEXP cpp_list_assign2(SEXP x, SEXP values, SEXP locs){
//   int NP = 0;
//
//   SEXP names = SHIELD(get_names(x)); ++NP;
//   SEXP col_names = SHIELD(Rf_getAttrib(values, R_NamesSymbol)); ++NP;
//
//   if (TYPEOF(x) != VECSXP){
//     YIELD(NP);
//     Rf_error("`x` must be a list in %s", __func__);
//   }
//   // if (TYPEOF(values) != VECSXP || (is_null(col_names) && !is_null(locs))){
//   //   YIELD(NP);
//   //   Rf_error("`values` must be a named list in %s", __func__);
//   // }
//
//   int n = Rf_length(x);
//   int n_cols = Rf_length(values);
//
//   if (is_null(names)){
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
//   if (is_null(locs)){
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
//     any_null = any_null || is_null(p_y[j]);
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
//   Rf_classgets(out, scalar_utf8_str("data.frame"));
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

  SEXP df_str = SHIELD(make_utf8_str("data.frame")); ++NP;
  SEXP row_names = SHIELD(create_df_row_names(N)); ++NP;

  // If no names then add names
  if (is_null(get_names(out))){
    set_names(out, new_vec(STRSXP, n_items));
  }
  Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  Rf_classgets(out, df_str);
  YIELD(NP);
  return out;
}


// Same as above but always in-place so assuming no NULL elements
void set_list_as_df(SEXP x) {
  int N; // Number of rows
  int NP = 0; // Number of protects
  int n_items = Rf_length(x);
  if (is_df(x)){
    N = df_nrow(x);
  } else if (n_items == 0){
    N = 0;
  } else {
    N = vec_length(VECTOR_ELT(x, 0));
  }

  SEXP df_str = SHIELD(make_utf8_str("data.frame")); ++NP;
  SEXP row_names = SHIELD(create_df_row_names(N)); ++NP;

  // If no names then add names
  if (is_null(get_names(x))){
    set_names(x, new_vec(STRSXP, n_items));
  }
  Rf_setAttrib(x, R_RowNamesSymbol, row_names);
  Rf_classgets(x, df_str);
  YIELD(NP);
}

// Can use in the future once I figure out how to re-write named_list() in C++

[[cpp11::register]]
SEXP cpp_new_df(SEXP x, SEXP nrows, bool recycle, bool name_repair){

  int NP = 0;

  SEXP out = x;

  if (recycle){
    SHIELD(out = cpp_recycle(out, nrows)); ++NP;
  } else {
    SHIELD(out = cpp_drop_null(out, true)); ++NP;
  }

  SEXP row_names;

  if (is_null(nrows)){
    if (Rf_length(out) == 0){
      row_names = SHIELD(new_vec(INTSXP, 0)); ++NP;
    } else {
      row_names = SHIELD(create_df_row_names(vec_length(VECTOR_ELT(out, 0)))); ++NP;
    }
  } else {
    row_names = SHIELD(create_df_row_names(Rf_asInteger(nrows))); ++NP;
  }

  SEXP out_names = SHIELD(get_names(out)); ++NP;
  if (is_null(out_names)){
    SHIELD(out_names = new_vec(STRSXP, Rf_length(out))); ++NP;
  } else {
    SHIELD(out_names = coerce_vec(out_names, STRSXP)); ++NP;
  }

  if (name_repair){
    SEXP dup_sep = SHIELD(make_utf8_str("_")); ++NP;
    SEXP empty_sep = SHIELD(make_utf8_str("col_")); ++NP;
    SHIELD(out_names = cpp_name_repair(out_names, dup_sep, empty_sep)); ++NP;
  }
  set_names(out, out_names);
  Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  Rf_classgets(out, make_utf8_str("data.frame"));
  YIELD(NP);
  return out;
}

// Multi-assign recycled variables to data frame

[[cpp11::register]]
SEXP cpp_df_assign_cols(SEXP x, SEXP cols){
  int NP = 0;

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame` in %s", __func__);
  }

  SEXP names = get_names(x);
  SEXP col_names = get_names(cols);

  if (TYPEOF(cols) != VECSXP || is_null(col_names)){
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
  int* RESTRICT p_add_locs = INTEGER(add_locs);
  int n_cols_to_add = na_count(add_locs, false);
  int out_size = n + n_cols_to_add;

  SEXP out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
  SEXP out_names = SHIELD(new_vec(STRSXP, out_size)); ++NP;

  // Initialise out
  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, p_x[i]);
    SET_STRING_ELT(out_names, i, p_names[i]);
  }

  bool any_null = false;

  int loc;
  SEXP vec = R_NilValue;

  for (int j = 0; j < n_cols; ++j){
    loc = p_add_locs[j];
    vec = p_cols[j];
    any_null = any_null || is_null(vec);
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
  set_names(out, out_names);
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(n_rows));
  Rf_classgets(out, make_utf8_str("data.frame"));
  SHIELD(out = reconstruct(out, x, false)); ++NP;
  YIELD(NP);
  return out;
}
