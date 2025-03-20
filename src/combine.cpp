#include "cheapr.h"

// Not fast but exists for low overhead C++ calls
SEXP cpp_unique(SEXP x){
  SEXP dup = Rf_protect(Rf_duplicated(x, FALSE));
  SEXP unique_locs = Rf_protect(cpp_which_(dup, true));
  SEXP out = Rf_protect(sset_vec(x, unique_locs, false));
  Rf_unprotect(3);
  return out;
}

SEXP cpp_setdiff(SEXP x, SEXP y){
  Rf_protect(x = cpp_unique(x));

  SEXP matches = Rf_protect(Rf_match(y, x, NA_INTEGER));
  SEXP locs = Rf_protect(cpp_which_na(matches));

  SEXP out = Rf_protect(sset_vec(x, locs, false));

  Rf_unprotect(4);
  return out;
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

[[cpp11::register]]
SEXP na_init(SEXP x, int n){
  SEXP ptype = Rf_protect(get_ptype(x));
  SEXP out = Rf_protect(cpp_rep_len(ptype, n));
  Rf_unprotect(2);
  return out;
}

// Prototypes of data frame

[[cpp11::register]]
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

  // We do 3 passes
  // 1st pass: Check inputs are dfs and construct prototype list
  // 2nd pass: initialise data frames (if need be) and order col names correctly
  // 3rd pass: combine everything

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
  }

  int n_cols = Rf_length(names);

  SEXP new_frame, name_locs;
  PROTECT_INDEX new_frame_idx, name_locs_idx;
  R_ProtectWithIndex(new_frame = R_NilValue, &new_frame_idx); ++NP;
  R_ProtectWithIndex(name_locs = R_NilValue, &name_locs_idx); ++NP;

  for (int i = 0; i < n_frames; ++i){
    R_Reprotect(df = p_x[i], df_idx);

    R_Reprotect(new_names = cpp_setdiff(
      names, Rf_getAttrib(df, R_NamesSymbol)
    ), new_names_idx);

    if (Rf_length(new_names) > 0){
      R_Reprotect(name_locs = Rf_match(names, new_names, NA_INTEGER), name_locs_idx);
      R_Reprotect(new_frame = sset_vec(ptypes, name_locs, false), new_frame_idx);
      Rf_namesgets(new_frame, new_names);
      R_Reprotect(new_frame = cpp_list_as_df(new_frame), new_frame_idx);
      R_Reprotect(
        new_frame = na_init(new_frame, Rf_length(Rf_getAttrib(df, R_RowNamesSymbol))),
        new_frame_idx
      );
      SET_VECTOR_ELT(temp_list, 0, df);
      SET_VECTOR_ELT(temp_list, 1, new_frame);
      R_Reprotect(df = list_c(temp_list), df_idx);
      R_Reprotect(df = cpp_list_as_df(df), df_idx);
    }
    R_Reprotect(df = cpp_df_select(df, names), df_idx);
    SET_VECTOR_ELT(frames, i, df);
  }

  SEXP out = Rf_protect(Rf_allocVector(VECSXP, n_cols)); ++NP;
  SEXP vectors = Rf_protect(Rf_allocVector(VECSXP, n_frames)); ++NP;

  for (int j = 0; j < n_cols; ++j){
    for (int i = 0; i < n_frames; ++i){
      SET_VECTOR_ELT(vectors, i, VECTOR_ELT(VECTOR_ELT(frames, i), j));
    }
    SET_VECTOR_ELT(out, j, cpp_c(vectors));
  }
  Rf_protect(out = cpp_list_as_df(out)); ++NP;
  Rf_namesgets(out, names);
  Rf_unprotect(NP);
  return out;
}

// Working version for data frames with similar col structures
// SEXP cpp_df_c(SEXP x){
//
//   if (!Rf_isVectorList(x)){
//     Rf_error("`x` must be a list of data frames");
//   }
//
//   int n_frames = Rf_length(x);
//
//   if (n_frames == 0) return R_NilValue;
//
//   int NP = 0;
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//
//   SEXP df = p_x[0];
//
//   if (!Rf_inherits(df, "data.frame")){
//     Rf_unprotect(NP); Rf_error("Can't combine data frames with non data frames");
//   }
//
//   SEXP names = Rf_protect(Rf_getAttrib(df, R_NamesSymbol)); ++NP;
//   if (Rf_any_duplicated(names, FALSE) > 0){
//     Rf_unprotect(NP);
//     Rf_error("data frame names aren't unique, please check");
//   }
//
//   int n_cols = Rf_length(names);
//   int out_nrows = 0;
//   SEXP frames = Rf_protect(Rf_allocVector(VECSXP, n_frames)); ++NP;
//   SET_VECTOR_ELT(frames, 0, df);
//
//   for (int i = 1; i < n_frames; ++i){
//     df = Rf_protect(p_x[i]);
//
//     if (!Rf_inherits(df, "data.frame")){
//       Rf_unprotect(NP + 1); Rf_error("Can't combine data frames with non data frames");
//     }
//     if (Rf_length(df) != n_cols){
//       Rf_unprotect(NP + 1); Rf_error("All data frames must contain identically named variables");
//     }
//
//     out_nrows += Rf_length(Rf_getAttrib(df, R_RowNamesSymbol));
//     SET_VECTOR_ELT(frames, i, cpp_df_select(df, names));
//     Rf_unprotect(1);
//   }
//
//   SEXP out = Rf_protect(Rf_allocVector(VECSXP, n_cols)); ++NP;
//   SEXP vectors = Rf_protect(Rf_allocVector(VECSXP, n_frames)); ++NP;
//
//   for (int j = 0; j < n_cols; ++j){
//     for (int i = 0; i < n_frames; ++i){
//       SET_VECTOR_ELT(vectors, i, VECTOR_ELT(VECTOR_ELT(frames, i), j));
//     }
//     SET_VECTOR_ELT(out, j, cpp_c(vectors));
//   }
//   Rf_protect(out = cpp_list_as_df(out)); ++NP;
//   Rf_setAttrib(out, R_NamesSymbol, names);
//   Rf_unprotect(NP);
//   return out;
// }

// `c()` but no concatenation of names
[[cpp11::register]]
SEXP cpp_c(SEXP x){
  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list of vectors");
  }
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);
  if (n == 1){
    return p_x[0];
  }

  int vector_type = NILSXP;
  R_xlen_t out_size = 0;

  bool is_factor = false, is_df = false, is_classed = false;

  for (R_xlen_t i = 0; i < n; ++i){
    vector_type = std::max(vector_type, TYPEOF(p_x[i]));
    out_size += Rf_xlength(p_x[i]);
    is_factor = is_factor || Rf_isFactor(p_x[i]);
    is_df = is_df || Rf_inherits(p_x[i], "data.frame");
    is_classed = is_classed || Rf_isObject(p_x[i]);
  }

  if (is_df){
    return cpp_df_c(x);
  }

  if (is_factor){
    return(cpp11::package("cheapr")["combine_factors"](x));
  }

  if (is_classed){
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

    for (R_xlen_t i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        Rf_coerceVector(p_x[i], vector_type), temp_idx
      );
      int *p_temp = INTEGER(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        p_out[k] = p_temp[j];
      }
    }
    break;
  }
  case REALSXP: {
    out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
    double *p_out = REAL(out);

    for (R_xlen_t i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        Rf_coerceVector(p_x[i], vector_type), temp_idx
      );
      double *p_temp = REAL(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        p_out[k] = p_temp[j];
      }
    }
    break;
  }
  case STRSXP: {
    out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;

    for (R_xlen_t i = 0; i < n; ++i){
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

    for (R_xlen_t i = 0; i < n; ++i){
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

    for (R_xlen_t i = 0; i < n; ++i){
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
  Rf_unprotect(NP);
  return out;
}
