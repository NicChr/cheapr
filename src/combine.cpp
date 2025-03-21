#include "cheapr.h"


[[cpp11::register]]
SEXP cpp_rep_len(SEXP x, int length){
  int out_size = length;

  if (Rf_inherits(x, "data.frame")){
    int n_cols = Rf_length(x);
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, n_cols));
    SEXP var;
    for (int i = 0; i < n_cols; ++i){
      var = Rf_protect(VECTOR_ELT(x, i));
      SET_VECTOR_ELT(out, i, cpp_rep_len(var, length));
      Rf_unprotect(1);
    }
    Rf_namesgets(out, Rf_getAttrib(x, R_NamesSymbol));
    Rf_protect(out = cpp_list_as_df(out));
    Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(length));
    Rf_unprotect(2);
    return out;
  } else if (is_simple_atomic_vec(x)){

    int size = Rf_length(x);
    int n_chunks, k, m, chunk_size;

    // Return x if length(x) == length
    if (out_size == size) return x;

    switch (TYPEOF(x)){
    case LGLSXP:
    case INTSXP: {
      int *p_x = INTEGER(x);
      SEXP out = Rf_protect(Rf_allocVector(TYPEOF(x), out_size));
      int *p_out = INTEGER(out);

      if (out_size > 0 && size > 0){
        n_chunks = std::ceil((static_cast<double>(out_size)) / size);
        for (int i = 0; i < n_chunks; ++i){
          k = ( (i + 1) * size) - size;
          m = std::min(k + size - 1, out_size - 1);
          chunk_size = m - k + 1;
          memcpy(&p_out[k], &p_x[0], chunk_size * sizeof(int));
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        OMP_FOR_SIMD
        for (int i = 0; i < out_size; ++i){
          p_out[i] = NA_INTEGER;
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      Rf_unprotect(1);
      return out;
    }
    case REALSXP: {
      double *p_x = REAL(x);
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);

      if (out_size > 0 && size > 0){
        n_chunks = std::ceil((static_cast<double>(out_size)) / size);
        for (int i = 0; i < n_chunks; ++i){
          k = ( (i + 1) * size) - size;
          m = std::min(k + size - 1, out_size - 1);
          chunk_size = m - k + 1;
          memcpy(&p_out[k], &p_x[0], chunk_size * sizeof(double));
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        OMP_FOR_SIMD
        for (int i = 0; i < out_size; ++i){
          p_out[i] = NA_REAL;
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      Rf_unprotect(1);
      return out;
    }
    case STRSXP: {
      const SEXP *p_x = STRING_PTR_RO(x);
      SEXP out = Rf_protect(Rf_allocVector(STRSXP, out_size));

      if (out_size > 0 && size > 0){
        for (int i = 0; i < out_size; ++i){
          SET_STRING_ELT(out, i, p_x[i % size]);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        for (int i = 0; i < out_size; ++i){
          SET_STRING_ELT(out, i, NA_STRING);
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      Rf_unprotect(1);
      return out;
    }
    case CPLXSXP: {
      Rcomplex *p_x = COMPLEX(x);
      SEXP out = Rf_protect(Rf_allocVector(CPLXSXP, out_size));
      Rcomplex *p_out = COMPLEX(out);

      if (out_size > 0 && size > 0){
        for (int i = 0; i < out_size; ++i){
          SET_COMPLEX_ELT(out, i, p_x[i % size]);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        for (int i = 0; i < out_size; ++i){
          p_out[i].r = NA_REAL;
          p_out[i].i = NA_REAL;
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      Rf_unprotect(1);
      return out;
    }
    default: {
      return base_rep(x, cpp11::named_arg("length.out") = length);
    }
    }
  } else {
    return base_rep(x, cpp11::named_arg("length.out") = length);
  }
}

// Recycle elements of a list `x`

[[cpp11::register]]
SEXP cpp_recycle(SEXP x, SEXP length){
  SEXP out = Rf_protect(cpp_drop_null(x, true));
  SEXP sizes = Rf_protect(cpp_lengths(out, false));
  int *p_sizes = INTEGER(sizes);
  bool has_length = !Rf_isNull(length);
  Rf_protect(length = Rf_coerceVector(length, INTSXP));

  int n = 0;
  int n_objs = Rf_length(out);

  if (!has_length){
    if (n_objs > 0){
      // We calculate `max(sizes)`
      // We won't have any NA in sizes so no NA checking is needed
      OMP_FOR_SIMD
      for (int i = 0; i < n_objs; ++i){
        n = p_sizes[i] > n ? p_sizes[i] : n;
      }
    }
  } else {
    n = Rf_asInteger(length);
  }

  SEXP r_zero = Rf_protect(Rf_ScalarInteger(0));
  if (!has_length && scalar_count(sizes, r_zero, false) > 0) n = 0;

  for (int i = 0; i < n_objs; ++i){
    SET_VECTOR_ELT(out, i, cpp_rep_len(VECTOR_ELT(out, i), n));
  }
  Rf_unprotect(4);
  return out;
}


// Fast unique that can be used in C code
// Doesn't return unique df rows
SEXP cpp_unique(SEXP x){

  bool is_simple = is_simple_atomic_vec(x);

  if (is_compact_seq(x)){
    return x;
  } else if (is_simple && Rf_length(x) < 100000){
    SEXP dup = Rf_protect(Rf_duplicated(x, FALSE));
    SEXP unique_locs = Rf_protect(cpp_which_(dup, true));
    SEXP out = Rf_protect(sset_vec(x, unique_locs, false));
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    Rf_unprotect(3);
    return out;
  } else if (is_simple){
    return cheapr_fast_unique(x);
  } else {
    return cpp11::package("base")["unique"](x);
  }
}

SEXP cpp_setdiff(SEXP x, SEXP y){

  bool is_simple = is_simple_atomic_vec(x) && is_simple_atomic_vec(y);

  if (is_simple){
    SEXP ux = Rf_protect(cpp_unique(x));
    SEXP matches;
    if (Rf_length(ux) < 100000){
      matches = Rf_protect(Rf_match(y, ux, NA_INTEGER));
    } else {
      matches = Rf_protect(cheapr_fast_match(ux, y));
    }
    SEXP locs = Rf_protect(cpp_which_na(matches));
    SEXP out = Rf_protect(sset_vec(ux, locs, false));

    SHALLOW_DUPLICATE_ATTRIB(out, x);
    Rf_unprotect(4);
    return out;
  } else {
    return cpp11::package("base")["setdiff"](x, y);
  }

}

SEXP get_ptype(SEXP x){
  SEXP r_zero = Rf_protect(Rf_ScalarInteger(0));
  SEXP out;
  if (Rf_inherits(x, "data.frame")){
    out = Rf_protect(cpp_df_slice(x, r_zero, true));
  } else if (is_simple_atomic_vec(x)){
    out = Rf_protect(cpp_sset(x, r_zero));
  } else {
    out = Rf_protect(base_sset(x, r_zero));
  }
  Rf_unprotect(2);
  return out;
}

SEXP na_init(SEXP x, int n){
  SEXP ptype = Rf_protect(get_ptype(x));
  SEXP out = Rf_protect(cpp_rep_len(ptype, n));
  Rf_unprotect(2);
  return out;
}

// Prototypes of data frame

SEXP get_ptypes(SEXP x){
  int n = Rf_length(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, n));

  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, get_ptype(VECTOR_ELT(x, i)));
  }

  Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));

  Rf_unprotect(1);
  return out;

}

SEXP factor_as_character(SEXP x){
  SEXP levels = Rf_protect(Rf_getAttrib(x, R_LevelsSymbol));
  SEXP out = Rf_protect(sset_vec(levels, x, true));

  Rf_unprotect(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_combine_levels(SEXP x){
  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list of factors in %s", __func__);
  }
  int n = Rf_length(x);
  SEXP levels = Rf_protect(Rf_allocVector(VECSXP, n));
  const SEXP *p_x = VECTOR_PTR_RO(x);

  int out_size = 0;

  SEXP fct_levels;

  PROTECT_INDEX fct_idx;
  R_ProtectWithIndex(fct_levels = R_NilValue, &fct_idx);

  for (int i = 0; i < n; ++i){
    if (Rf_isFactor(p_x[i])){
      R_Reprotect(fct_levels = Rf_getAttrib(p_x[i], R_LevelsSymbol), fct_idx);
    } else {
      R_Reprotect(fct_levels = base_as_character(p_x[i]), fct_idx);
    }
    out_size += Rf_length(fct_levels);
    SET_VECTOR_ELT(levels, i, fct_levels);
  }
  SEXP out = Rf_protect(cpp_c(levels));
  Rf_protect(out = cpp_unique(out));
  // Rf_protect(out = collapse_unique(out));
  Rf_unprotect(4);
  return out;
}

[[cpp11::register]]
SEXP cpp_combine_factors(SEXP x){

  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list of factors in %s", __func__);
  }

  int n = Rf_length(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

  int out_size = 0;

  SEXP levels = Rf_protect(cpp_combine_levels(x));
  SEXP chars = Rf_protect(Rf_allocVector(VECSXP, n));
  SEXP char_vec;

  PROTECT_INDEX char_vec_idx;
  R_ProtectWithIndex(char_vec = R_NilValue, &char_vec_idx);

  for (int i = 0; i < n; ++i){
    if (Rf_isFactor(p_x[i])){
      R_Reprotect(char_vec = factor_as_character(p_x[i]), char_vec_idx);
    } else {
      R_Reprotect(char_vec = base_as_character(p_x[i]), char_vec_idx);
    }
    out_size += Rf_length(char_vec);
    SET_VECTOR_ELT(chars, i, char_vec);
  }

  R_Reprotect(char_vec = cpp_c(chars), char_vec_idx);
  SEXP out = Rf_protect(cheapr_factor(char_vec,
    cpp11::named_arg("levels") = levels,
    cpp11::named_arg("na_exclude") = !cpp_any_na(levels, false)
  ));
  Rf_unprotect(4);
  return out;
}

// Helper to concatenate plain lists

[[cpp11::register]]
SEXP list_c(SEXP x){
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

  R_xlen_t out_size = 0;
  for (R_xlen_t i = 0; i < n; ++i) out_size += Rf_xlength(p_x[i]);

  R_xlen_t k = 0;
  SEXP out;

  out = Rf_protect(Rf_allocVector(VECSXP, out_size)); ++NP;

  SEXP names;
  PROTECT_INDEX nm_idx;
  R_ProtectWithIndex(names = R_NilValue, &nm_idx); ++NP;

  SEXP out_names = Rf_protect(Rf_allocVector(STRSXP, out_size)); ++NP;
  bool any_names = false;

  for (R_xlen_t i = 0; i < n; ++i){
    R_xlen_t m = Rf_xlength(p_x[i]);
    const SEXP *p_temp = VECTOR_PTR_RO(p_x[i]);
    R_Reprotect(names = Rf_getAttrib(p_x[i], R_NamesSymbol), nm_idx);
    any_names = any_names || !Rf_isNull(names);
    if (!Rf_isNull(names)){
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
        SET_STRING_ELT(out_names, k, STRING_ELT(names, j));
      }
    } else {
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
      }
    }
  }
  if (any_names){
    Rf_namesgets(out, out_names);
  }
  Rf_unprotect(NP);
  return out;
}

// Concatenate a list of data frames
SEXP cpp_df_c(SEXP x){

  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list of data frames");
  }

  int n_frames = Rf_length(x);

  if (n_frames == 0) return R_NilValue;

  int NP = 0;
  const SEXP *p_x = VECTOR_PTR_RO(x);

  SEXP df, names;
  PROTECT_INDEX df_idx, names_idx;

  R_ProtectWithIndex(df = p_x[0], &df_idx); ++NP;
  R_ProtectWithIndex(names = Rf_getAttrib(df, R_NamesSymbol), &names_idx); ++NP;

  if (!Rf_inherits(df, "data.frame")){
    Rf_unprotect(NP); Rf_error("Can't combine data frames with non data frames");
  }

  SEXP frames = Rf_protect(Rf_allocVector(VECSXP, n_frames)); ++NP;
  SET_VECTOR_ELT(frames, 0, df);

  SEXP ptypes, new_names, new_ptypes, new_cols,
  temp_list;
  PROTECT_INDEX
  new_names_idx, ptypes_idx, new_ptypes_idx, new_cols_idx,
  temp_list_idx;

  R_ProtectWithIndex(ptypes = get_ptypes(df), &ptypes_idx); ++NP;
  R_ProtectWithIndex(new_names = R_NilValue, &new_names_idx); ++NP;
  R_ProtectWithIndex(new_ptypes = R_NilValue, &new_ptypes_idx); ++NP;
  R_ProtectWithIndex(new_cols = R_NilValue, &new_cols_idx); ++NP;
  R_ProtectWithIndex(temp_list = Rf_allocVector(VECSXP, 2), &temp_list_idx); ++NP;

  // We do 2 passes
  // 1st pass: Check inputs are dfs and construct prototype list
  // 2nd pass: initialise data frame vecs (if need be) and combine simultaneously

  int out_size = df_nrow(df);

  for (int i = 1; i < n_frames; ++i){
    R_Reprotect(df = p_x[i], df_idx);

    if (!Rf_inherits(df, "data.frame")){
      Rf_unprotect(NP); Rf_error("Can't combine data frames with non data frames");
    }

    R_Reprotect(new_names = cpp_setdiff(
      Rf_getAttrib(df, R_NamesSymbol),
      Rf_getAttrib(ptypes, R_NamesSymbol)
    ), new_names_idx);

    // Adjust prototype list
    if (Rf_length(new_names) > 0){
      R_Reprotect(new_cols = cpp_df_select(df, new_names), new_cols_idx);
      R_Reprotect(new_ptypes = get_ptypes(new_cols), new_ptypes_idx);
      SET_VECTOR_ELT(temp_list, 0, ptypes);
      SET_VECTOR_ELT(temp_list, 1, new_ptypes);
      R_Reprotect(ptypes = cpp_c(temp_list), ptypes_idx);
      SET_VECTOR_ELT(temp_list, 0, names);
      SET_VECTOR_ELT(temp_list, 1, new_names);
      R_Reprotect(names = cpp_c(temp_list), names_idx);
      Rf_setAttrib(ptypes, R_NamesSymbol, names);
    }
    out_size += df_nrow(df);
  }

  int n_cols = Rf_length(names);

  SEXP vec, name_locs;
  PROTECT_INDEX vec_idx, name_locs_idx;
  R_ProtectWithIndex(vec = R_NilValue, &vec_idx); ++NP;
  R_ProtectWithIndex(name_locs = R_NilValue, &name_locs_idx); ++NP;

  SEXP out = Rf_protect(Rf_allocVector(VECSXP, n_cols)); ++NP;
  SEXP vectors = Rf_protect(Rf_allocVector(VECSXP, n_frames)); ++NP;

  for (int j = 0; j < n_cols; ++j){
    for (int i = 0; i < n_frames; ++i){
      R_Reprotect(vec = get_list_element(p_x[i], CHAR(STRING_ELT(names, j))), vec_idx);

      if (Rf_isNull(vec)){
        R_Reprotect(vec = VECTOR_ELT(ptypes, j), vec_idx);
        R_Reprotect(vec = na_init(vec, df_nrow(p_x[i])), vec_idx);
      }
      SET_VECTOR_ELT(vectors, i, vec);
    }
    SET_VECTOR_ELT(out, j, cpp_c(vectors));
  }
  Rf_protect(out = cpp_list_as_df(out)); ++NP;
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(out_size));
  Rf_namesgets(out, names);
  Rf_unprotect(NP);
  return out;
}

// `c()` but no concatenation of names
[[cpp11::register]]
SEXP cpp_c(SEXP x){
  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list of vectors");
  }
  int NP = 0;
  int n = Rf_length(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);
  if (n == 1){
    return p_x[0];
  }

  int vector_type = NILSXP;
  R_xlen_t out_size = 0;

  bool is_factor = false, is_df = false, is_classed = false,
    is_simple = false, is_date = false, is_datetime = false;


  // We use the tz info of the first datetime vec to copy to final result
  int first_datetime = integer_max_;

  for (int i = 0; i < n; ++i){
    vector_type = std::max(vector_type, TYPEOF(p_x[i]));
    out_size += Rf_xlength(p_x[i]);
    is_factor = is_factor || Rf_isFactor(p_x[i]);
    is_simple = is_simple || is_simple_atomic_vec(p_x[i]);
    is_date = is_date || Rf_inherits(p_x[i], "Date");
    is_datetime = is_datetime || Rf_inherits(p_x[i], "POSIXct");
    if (is_datetime){
      first_datetime = std::min(i, first_datetime);
    }
    is_df = is_df || Rf_inherits(p_x[i], "data.frame");
    is_classed = is_classed || Rf_isObject(p_x[i]);
  }

  // Date vectors can be ints but datetimes can't
  bool is_date2 = is_date && !is_datetime && (vector_type == INTSXP || vector_type == REALSXP);
  bool is_datetime2 = !is_date && is_datetime && vector_type == REALSXP;

  if (is_df){
    return cpp_df_c(x);
  }

  if (is_factor){
    return cpp_combine_factors(x);
  }

  if (is_classed && !(is_date2 || is_datetime2)){
    SEXP c_char = Rf_protect(Rf_mkString("c"));
    SEXP out = Rf_protect(
      cpp11::package("base")["do.call"](c_char, x)
    );
    Rf_unprotect(2);
    return out;
  }

  R_xlen_t k = 0;
  SEXP out, temp;

  PROTECT_INDEX temp_idx;
  R_ProtectWithIndex(temp = R_NilValue, &temp_idx); ++NP;


  switch(vector_type){
  case NILSXP: {
    out = Rf_protect(R_NilValue); ++NP;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
    int *p_out = INTEGER(out);

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        Rf_coerceVector(p_x[i], vector_type), temp_idx
      );
      int *p_temp = INTEGER(temp);
      memcpy(&p_out[k], &p_temp[0], m * sizeof(int));
      k += m;
    }
    break;
  }
  case REALSXP: {
    out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
    double *p_out = REAL(out);

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        Rf_coerceVector(p_x[i], vector_type), temp_idx
      );
      double *p_temp = REAL(temp);
      memcpy(&p_out[k], &p_temp[0], m * sizeof(double));
      k += m;
    }
    break;
  }
  case STRSXP: {
    out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        Rf_coerceVector(p_x[i], vector_type), temp_idx
      );
      const SEXP *p_temp = STRING_PTR_RO(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_STRING_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }
  case CPLXSXP: {
    out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        Rf_coerceVector(p_x[i], vector_type), temp_idx
      );
      Rcomplex *p_temp = COMPLEX(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_COMPLEX_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }
  case VECSXP: {
    out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        Rf_coerceVector(p_x[i], vector_type), temp_idx
      );
      const SEXP *p_temp = VECTOR_PTR_RO(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }

  default: {
    SEXP c_char = Rf_protect(Rf_mkString("c")); ++NP;
    out = Rf_protect(base_do_call(c_char, x)); ++NP;
    break;
  }
  }
  if (is_date2){
    Rf_classgets(out, Rf_mkChar("Date"));
  }
  if (is_datetime2){
    SHALLOW_DUPLICATE_ATTRIB(out, p_x[first_datetime]);
  }
  Rf_unprotect(NP);
  return out;
}
