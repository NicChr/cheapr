#include "cheapr.h"

// From R you pass `list(...)` and `.args` to below function
// which then decides which arguments to use

[[cpp11::register]]
SEXP cpp_list_args(SEXP args1, SEXP args2){
  bool use_dots = Rf_length(args1) != 0;
  bool use_list = !is_null(args2);

  if (use_dots && use_list){
    Rf_error("Please supply either `...` or `.args` in %s", __func__);
  }
  if (use_list && TYPEOF(args2) != VECSXP){
    Rf_error("`.args` must be a list in %s", __func__);
  }
  return use_list ? args2 : args1;
}

R_xlen_t unnested_length(SEXP x){
  R_CheckStack(); // Check C Stack size isn't close to the limit
  if (TYPEOF(x) != VECSXP){
    return Rf_xlength(x);
  }
  const SEXP *p_x = list_ptr_ro(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t out = 0;
  for (R_xlen_t i = 0; i < n; ++i){
    out += (TYPEOF(p_x[i]) == VECSXP) ? unnested_length(p_x[i]) : Rf_xlength(p_x[i]);
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_unnested_length(SEXP x){
  return as_vector(unnested_length(x));
}

[[cpp11::register]]
SEXP cpp_lengths(SEXP x, bool names){
  R_xlen_t n = Rf_xlength(x);
  SEXP out = SHIELD(vec::new_vector<int>(n));
  int* RESTRICT p_out = integer_ptr(out);
  if (TYPEOF(x) != VECSXP){
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = 1;
    }
  } else {
    const SEXP *p_x = list_ptr_ro(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      p_out[i] = vec::length(p_x[i]);
    }
  }
  SEXP x_names = SHIELD(get_old_names(x));
  if (names){
    set_old_names(out, x_names);
  }
  YIELD(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_new_list(SEXP size, SEXP default_value){
  SHIELD(size = cast<r_doubles_t>(size, r_null));
  if (Rf_length(size) != 1){
    Rf_error("`size` must be a vector of length 1");
  }
  R_xlen_t out_size = real_ptr(size)[0];
  SEXP out = SHIELD(new_list(out_size, default_value));
  YIELD(2);
  return out;
}

uint_fast64_t null_count(SEXP x){
  uint_fast64_t n = Rf_xlength(x);
  const SEXP *p_x = list_ptr_ro(x);
  uint_fast64_t n_null = 0;
  for (uint_fast64_t i = 0; i < n; ++i) n_null += is_null(p_x[i]);
  return n_null;
}

// Remove NULL elements from list

[[cpp11::register]]
SEXP cpp_drop_null(SEXP x){
  const SEXP *p_l = list_ptr_ro(x);
  uint_fast64_t n = Rf_xlength(x);
  uint_fast64_t n_null = null_count(x);
  SEXP names = SHIELD(get_old_names(x));

  if (n_null == 0){
    // Always return a plain-list
    SEXP out = SHIELD(new_list(n));
    for (uint_fast64_t i = 0; i < n; ++i) SET_VECTOR_ELT(out, i, p_l[i]);
    set_old_names(out, names);
    YIELD(2);
    return out;
  }


  // Subset on both the list and names of the list
  uint_fast64_t n_keep = n - n_null;
  uint_fast64_t k = 0;
  SEXP out = SHIELD(new_list(n_keep));

  if (is_null(names)){
    for (uint_fast64_t i = 0; i < n; ++i){
      if (is_null(p_l[i])) continue;
      SET_VECTOR_ELT(out, k++, p_l[i]);
    }
    YIELD(2);
    return out;
  } else {
    SEXP out_names = SHIELD(new_vector<r_string_t>(n_keep));
    const r_string_t *p_names = string_ptr_ro(names);
    for (uint_fast64_t i = 0; i < n; ++i){
      if (is_null(p_l[i])) continue;
      SET_VECTOR_ELT(out, k, p_l[i]);
      set_value<r_string_t>(out_names, k, p_names[i]);
      ++k;
    }
    set_old_names(out, out_names);
    YIELD(3);
    return out;
  }
}

SEXP which_not_null(SEXP x){
  const SEXP *p_x = list_ptr_ro(x);
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t n_null = null_count(x);
  R_xlen_t n_keep = n - n_null;
  R_xlen_t whichj = 0;
  R_xlen_t j = 0;

  // Which list elements should we keep?

  SEXP keep = SHIELD(vec::new_vector<int>(n_keep));
  int* RESTRICT p_keep = integer_ptr(keep);
  while (whichj < n_keep){
    p_keep[whichj] = j + 1;
    whichj += !is_null(p_x[j++]);
  }
  YIELD(1);
  return keep;
}

// Multi-assign named list elements

[[cpp11::register]]
SEXP cpp_list_assign(SEXP x, SEXP values){
  int32_t NP = 0;
  int n = Rf_length(x);
  int n_cols = Rf_length(values);

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list in %s", __func__);
  }
  if (TYPEOF(values) != VECSXP){
    Rf_error("`values` must be a named list in %s", __func__);
  }

  SEXP names = SHIELD(get_old_names(x)); ++NP;
  SEXP col_names = SHIELD(get_old_names(values)); ++NP;

  if (is_null(names)){
    SHIELD(names = new_vector<r_string_t>(n)); ++NP;
  }

  bool empty_value_names = is_null(col_names);

  if (empty_value_names){
    SHIELD(col_names = new_vector<r_string_t>(n_cols)); ++NP;
  }

  const SEXP *p_x = list_ptr_ro(x);
  const r_string_t *p_names = string_ptr_ro(names);
  const SEXP *p_y = list_ptr_ro(values);
  const r_string_t *p_col_names = string_ptr_ro(col_names);

  // We create an int vec to keep track of locations where to add col vecs

  SEXP add_locs;
  int n_cols_to_add;

  // If values is an unnamed list then we can simply append the values
  if (empty_value_names){
    add_locs = SHIELD(vec::new_vector<int>(0)); ++NP;
    SHIELD(add_locs = cpp_rep_len(add_locs, n_cols)); ++NP;
    n_cols_to_add = n_cols;
  } else {
    add_locs = SHIELD(match(names, col_names, na::integer)); ++NP;
    n_cols_to_add = na_count(add_locs, false);
  }
  int* RESTRICT p_add_locs = integer_ptr(add_locs);
  int out_size = n + n_cols_to_add;

  int loc;
  SEXP out = SHIELD(new_list(out_size)); ++NP;
  SEXP out_names = SHIELD(new_vector<r_string_t>(out_size)); ++NP;

  // Initialise out
  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, p_x[i]);
    set_value<r_string_t>(out_names, i, p_names[i]);
  }

  int null_count = 0;

  for (int j = 0; j < n_cols; ++j){

    // If loc == NA then we're adding a new element
    // otherwise we're modifying an existing one
    loc = p_add_locs[j];

    if (is_r_na(loc) && is_null(p_y[j])){
      ++null_count;
    } else {
      null_count += is_null(p_y[j]);
    }
  }

  SEXP null_locs = SHIELD(vec::new_vector<int>(null_count)); ++NP;
  int *p_null_locs = integer_ptr(null_locs);
  int nulli = 0;

  for (int j = 0; j < n_cols; ++j){

    // If loc == NA then we're adding a new element
    // otherwise we're modifying an existing one
    loc = p_add_locs[j];

    if (is_r_na(loc)){
      if (is_null(p_y[j])){
        p_null_locs[nulli++] = -(n + 1);
      }
      SET_VECTOR_ELT(out, n, p_y[j]);
      set_value<r_string_t>(out_names, n, p_col_names[j]);
      ++n;
    } else {
      if (is_null(p_y[j])){
        p_null_locs[nulli++] = -loc;
      }
      --loc;
      SET_VECTOR_ELT(out, loc, p_y[j]);
      set_value<r_string_t>(out_names, loc, p_col_names[j]);
    }
  }
  if (null_count != 0){
    SEXP keep = SHIELD(exclude_locs(null_locs, out_size)); ++NP;
    SHIELD(out = sset_vec(out, keep, false)); ++NP;
    SHIELD(out_names = sset_vec(out_names, keep, false)); ++NP;
  }
  set_old_names(out, out_names);
  YIELD(NP);
  return out;
}

// Data-frames

// Remove dimensions from arrays
// because cheapr doesn't work with arrays
SEXP maybe_cast_array(SEXP x){
  if (!Rf_isArray(x)){
    return x;
  } else {
    Rprintf("cheapr pkg cannot handle arrays. Array will be converted to a vector\n");
    SEXP out = SHIELD(vec::shallow_copy(x));
    attr::clear_attrs(out);
    YIELD(1);
    return out;
  }
}


// in-place convert list to df
// No recycling and assumes no elements
// must be a clean data frame like list
void set_list_as_df(SEXP x) {
  int N; // Number of rows
  int32_t NP = 0; // Number of protects
  int n_items = Rf_length(x);
  if (is_df(x)){
    N = df::nrow(x);
  } else if (n_items == 0){
    N = 0;
  } else {
    N = vec::length(VECTOR_ELT(x, 0));
  }

  SEXP df_str = SHIELD(as_vector("data.frame")); ++NP;

  // If no names then add names
  SEXP names = SHIELD(get_old_names(x)); ++NP;
  if (is_null(names)){
    SHIELD(names = new_vector<r_string_t>(n_items)); ++NP;
    set_old_names(x, names);
  }
  df::set_row_names(x, N);
  attr::set_old_class(x, df_str);
  YIELD(NP);
}

[[cpp11::register]]
SEXP cpp_new_df(SEXP x, SEXP nrows, bool recycle, bool name_repair){

  int32_t NP = 0;

  SEXP out = SHIELD(cpp_drop_null(x)); ++NP;

  // Remove array dimensions

  // We can assign in-place because we have shallow-duplicated above
  for (int i = 0; i < Rf_length(out); ++i){
    SET_VECTOR_ELT(out, i, maybe_cast_array(VECTOR_ELT(out, i)));
  }

  if (recycle){
    if (is_null(nrows)){
      recycle_in_place(out, length_common(out));
    } else {
      recycle_in_place(out, Rf_asInteger(nrows));
    }
  }

  int num_row;

  if (is_null(nrows)){
    if (Rf_length(out) == 0){
      num_row = 0;
    } else {
      num_row = vec::length(VECTOR_ELT(out, 0));
    }
  } else {
    SHIELD(nrows = cast<r_integers_t>(nrows, r_null)); ++NP;
    num_row = get_value<int>(nrows, 0);
  }

  SEXP out_names = SHIELD(get_old_names(out)); ++NP;
  if (is_null(out_names)){
    SHIELD(out_names = new_vector<r_string_t>(Rf_length(out))); ++NP;
  } else {
    SHIELD(out_names = vec::coerce_vec(out_names, STRSXP)); ++NP;
  }

  if (name_repair){
    SEXP dup_sep = SHIELD(as_vector("_")); ++NP;
    SEXP empty_sep = SHIELD(as_vector("col_")); ++NP;
    SHIELD(out_names = cpp_name_repair(out_names, dup_sep, empty_sep)); ++NP;
  }
  set_old_names(out, out_names);
  df::set_row_names(out, num_row);
  SEXP df_cls = SHIELD(as_vector("data.frame")); ++NP;
  attr::set_old_class(out, df_cls);
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_list_as_df(SEXP x) {
  return cpp_new_df(x, r_null, false, false);
}

// Multi-assign recycled variables to data frame

[[cpp11::register]]
SEXP cpp_df_assign_cols(SEXP x, SEXP cols){
  int32_t NP = 0;

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame` in %s", __func__);
  }

  SEXP names = SHIELD(get_old_names(x)); ++NP;
  SEXP col_names = SHIELD(get_old_names(cols)); ++NP;

  if (TYPEOF(cols) != VECSXP || is_null(col_names)){
    Rf_error("`cols` must be a named list in %s", __func__);
  }

  const SEXP *p_x = list_ptr_ro(x);
  const r_string_t *p_names = string_ptr_ro(names);
  const SEXP *p_cols = list_ptr_ro(cols);
  const r_string_t *p_col_names = string_ptr_ro(col_names);

  int n = Rf_length(x);
  int n_cols = Rf_length(cols);
  int n_rows = df::nrow(x);

  // We create an int vec to keep track of locations where to add col vecs

  SEXP add_locs = SHIELD(match(names, col_names, na::integer)); ++NP;
  int* RESTRICT p_add_locs = integer_ptr(add_locs);
  int n_cols_to_add = na_count(add_locs, false);
  int out_size = n + n_cols_to_add;

  SEXP out = SHIELD(new_list(out_size)); ++NP;
  SEXP out_names = SHIELD(new_vector<r_string_t>(out_size)); ++NP;

  // Initialise out
  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, p_x[i]);
    set_value<r_string_t>(out_names, i, p_names[i]);
  }

  bool any_null = false;

  int loc;
  SEXP vec = r_null;

  for (int j = 0; j < n_cols; ++j){
    loc = p_add_locs[j];
    vec = p_cols[j];
    any_null = any_null || is_null(vec);
    if (is_r_na(loc)){
      SET_VECTOR_ELT(out, n, cpp_rep_len(vec, n_rows));
      set_value<r_string_t>(out_names, n, p_col_names[j]);
      ++n;
    } else {
      --loc;
      SET_VECTOR_ELT(out, loc, cpp_rep_len(vec, n_rows));
      set_value<r_string_t>(out_names, loc, p_col_names[j]);
    }
  }
  if (any_null){
    SEXP keep = SHIELD(which_not_null(out)); ++NP;
    SHIELD(out = sset_vec(out, keep, false)); ++NP;
    SHIELD(out_names = sset_vec(out_names, keep, false)); ++NP;
  }
  set_old_names(out, out_names);
  df::set_row_names(out, n_rows);
  SEXP df_cls = SHIELD(as_vector("data.frame")); ++NP;
  attr::set_old_class(out, df_cls);
  SHIELD(out = rebuild(out, x, false)); ++NP;
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_as_df(SEXP x){
  if (inherits1(x, "data.frame")){
    SEXP n_rows = SHIELD(as_vector(df::nrow(x)));
    SEXP out = SHIELD(cpp_new_df(x, n_rows, false, false));
    YIELD(2);
    return out;
  } else if (is_null(x)){
    return init<r_data_frame_t>(0, false);
  } else if (Rf_isArray(x)){
    SEXP vec = SHIELD(maybe_cast_array(x));
    SEXP out = SHIELD(cpp_as_df(vec));
    YIELD(2);
    return out;
  } else if (cheapr_is_simple_atomic_vec2(x)){
    SEXP x_names = SHIELD(get_old_names(x));
    SEXP out = SHIELD(make_list(
      arg("name") = x_names,
      arg("value") = x
    ));
    SHIELD(out = cpp_new_df(out, r_null, false, false));
    YIELD(3);
    return out;
  } else if (is_bare_list(x)){
    return cpp_new_df(x, r_null, true, true);
  } else {
    SEXP out = SHIELD(eval_pkg_fun("as.data.frame", "base", env::base_env, x));
    SEXP col_seq = SHIELD(cpp_seq_len(Rf_length(out)));
    SEXP col_str = SHIELD(as_vector("col_"));
    SEXP new_names = SHIELD(r_paste(R_BlankScalarString, r_null, col_str, col_seq));
    set_old_names(out, new_names);
    YIELD(4);
    return out;
  }
}
