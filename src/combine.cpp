#include "cheapr.h"

// static cpp11::writable::integers CHEAPR_ZERO(1);
// void constants_init(DllInfo* dll){
//   CHEAPR_ZERO[0] = 0;
// }

SEXP rebuild(SEXP x, SEXP source, bool shallow_copy){
  if (is_df(x)){
    if (is_bare_df(source)){
      if (!shallow_copy && is_bare_df(x)){
        return x;
      } else {
        SEXP out = SHIELD(shallow_copy ? Rf_shallow_duplicate(x) : x);
        Rf_classgets(out, Rf_getAttrib(source, R_ClassSymbol));
        YIELD(1);
        return out;
      }
    } else if (is_bare_tbl(source)){

      // Only copy the class if source is a plain tbl

      if (!shallow_copy && is_bare_tbl(x)){
        return x;
      } else {
        SEXP out = SHIELD(shallow_copy ? Rf_shallow_duplicate(x) : x);
        Rf_classgets(out, Rf_getAttrib(source, R_ClassSymbol));
        YIELD(1);
        return out;
      }
    } else {

      // Method dispatch, users can write `rebuild` methods which
      // this will use
      return cheapr_rebuild(x, source, cpp11::named_arg("shallow_copy") = shallow_copy);
    }
  } else {
    return cheapr_rebuild(x, source);
  }
}

[[cpp11::register]]
SEXP cpp_rep_len(SEXP x, int length){
  int out_size = length;

  if (is_null(x)){
    return R_NilValue;
  } else if (is_df(x)){
    if (out_size == df_nrow(x)) return x;
    int n_cols = Rf_length(x);
    SEXP out = SHIELD(new_vec(VECSXP, n_cols));
    const SEXP *p_x = LIST_PTR_RO(x);
    for (int i = 0; i < n_cols; ++i){
      SET_VECTOR_ELT(out, i, cpp_rep_len(p_x[i], length));
    }
    SEXP names = SHIELD(get_names(x));
    set_names(out, names);
    set_list_as_df(out);
    Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(length));
    SHIELD(out = rebuild(out, x, false));
    YIELD(3);
    return out;
  } else if (cheapr_is_simple_vec2(x)){

    int size = Rf_length(x);
    int n_chunks, k, chunk_size;

    // Return x if length(x) == length
    if (out_size == size) return x;

    switch (CHEAPR_TYPEOF(x)){
    case LGLSXP:
    case INTSXP: {
      const int *p_x = INTEGER_RO(x);
      SEXP out = SHIELD(new_vec(TYPEOF(x), out_size));
      int* RESTRICT p_out = INTEGER(out);

      if (size == 1){
        std::fill(p_out, p_out + out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        n_chunks = std::ceil(static_cast<double>(out_size) / size);
        for (int i = 0; i < n_chunks; ++i){
          k = i * size;
          chunk_size = std::min(k + size, out_size) - k;
          std::copy(p_x, p_x + chunk_size, p_out + k);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        std::fill(p_out, p_out + out_size, NA_INTEGER);
      }
      Rf_copyMostAttrib(x, out);
      YIELD(1);
      return out;
    }
    case CHEAPR_INT64SXP: {
      const int64_t *p_x = INTEGER64_PTR_RO(x);
      SEXP out = SHIELD(new_vec(REALSXP, out_size));
      int64_t* RESTRICT p_out = INTEGER64_PTR(out);

      if (size == 1){
        std::fill(p_out, p_out + out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        n_chunks = std::ceil((static_cast<double>(out_size)) / size);
        for (int i = 0; i < n_chunks; ++i){
          k = i * size;
          chunk_size = std::min(k + size, out_size) - k;
          std::copy(p_x, p_x + chunk_size, p_out + k);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        std::fill(p_out, p_out + out_size, NA_INTEGER64);
      }
      Rf_copyMostAttrib(x, out);
      YIELD(1);
      return out;
    }
    case REALSXP: {
      const double *p_x = REAL_RO(x);
      SEXP out = SHIELD(new_vec(REALSXP, out_size));
      double* RESTRICT p_out = REAL(out);

      if (size == 1){
        std::fill(p_out, p_out + out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        n_chunks = std::ceil((static_cast<double>(out_size)) / size);
        for (int i = 0; i < n_chunks; ++i){
          k = i * size;
          chunk_size = std::min(k + size, out_size) - k;
          std::copy(p_x, p_x + chunk_size, p_out + k);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        std::fill(p_out, p_out + out_size, NA_REAL);
      }
      Rf_copyMostAttrib(x, out);
      YIELD(1);
      return out;
    }
    case STRSXP: {
      const SEXP *p_x = STRING_PTR_RO(x);
      SEXP out = SHIELD(new_vec(STRSXP, out_size));

      if (size == 1){
        SEXP val = p_x[0];
        for (int i = 0; i < out_size; ++i){
          SET_STRING_ELT(out, i, val);
        }
      } else if (out_size > 0 && size > 0){
        for (int i = 0, xi = 0; i < out_size; xi = (++xi == size) ? 0 : xi, ++i){
          SET_STRING_ELT(out, i, p_x[xi]);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        for (int i = 0; i < out_size; ++i){
          SET_STRING_ELT(out, i, NA_STRING);
        }
      }
      Rf_copyMostAttrib(x, out);
      YIELD(1);
      return out;
    }
    case CPLXSXP: {
      Rcomplex *p_x = COMPLEX(x);
      SEXP out = SHIELD(new_vec(CPLXSXP, out_size));
      Rcomplex *p_out = COMPLEX(out);

      if (size == 1){
        std::fill(p_out, p_out + out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        for (int i = 0, xi = 0; i < out_size; xi = (++xi == size) ? 0 : xi, ++i){
          SET_COMPLEX_ELT(out, i, p_x[xi]);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        for (int i = 0; i < out_size; ++i){
          p_out[i].r = NA_REAL;
          p_out[i].i = NA_REAL;
        }
      }
      Rf_copyMostAttrib(x, out);
      YIELD(1);
      return out;
    }
    case VECSXP: {
      const SEXP *p_x = LIST_PTR_RO(x);
      SEXP out = SHIELD(new_vec(VECSXP, out_size));

      if (size == 1){
        SEXP val = p_x[0];
        for (int i = 0; i < out_size; ++i){
          SET_VECTOR_ELT(out, i, val);
        }
      } else if (out_size > 0 && size > 0){
        for (int i = 0, xi = 0; i < out_size; xi = (++xi == size) ? 0 : xi, ++i){
          SET_VECTOR_ELT(out, i, p_x[xi]);
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        for (int i = 0; i < out_size; ++i){
          SET_VECTOR_ELT(out, i, R_NilValue);
        }
      }
      Rf_copyMostAttrib(x, out);
      YIELD(1);
      return out;
    }
    default: {
      if (r_length(x) == length){
      return x;
    } else {
      return base_rep(x, cpp11::named_arg("length.out") = length);
    }
    }
    }
  } else {
    if (r_length(x) == length){
      return x;
    } else {
      return base_rep(x, cpp11::named_arg("length.out") = length);
    }
  }
}

[[cpp11::register]]
SEXP cpp_rep(SEXP x, SEXP times){

  R_xlen_t n = vector_length(x);

  R_xlen_t out_size;
  R_xlen_t n_times = Rf_xlength(times);
  SHIELD(times = cast<r_integer_t>(times, R_NilValue));

  if (n_times == 1){
    out_size = n * INTEGER(times)[0];
    SEXP out = SHIELD(cpp_rep_len(x, out_size));
    YIELD(2);
    return out;
  } else {
    int *p_times = INTEGER(times);
    if (n_times != n){
      YIELD(1);
      Rf_error("`times` must be length 1 or `vector_length(x)` in %s", __func__);
    }

    if (is_null(x)){
      YIELD(1);
      return R_NilValue;
    } else if (is_df(x)){
      if (Rf_length(x) == 0){
        SEXP out = SHIELD(Rf_shallow_duplicate(x));
        Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(cpp_sum(times)));
        SHIELD(out = rebuild(out, x, false));
        YIELD(3);
        return out;
      } else {
        SEXP row_names = SHIELD(cpp_seq_len(df_nrow(x)));
        SHIELD(row_names = cpp_rep(row_names, times));
        SEXP out = SHIELD(cpp_sset(x, row_names, true));
        YIELD(4);
        return out;
      }

    } else if (cheapr_is_simple_vec2(x)){
      SEXP out = SHIELD(new_vec(TYPEOF(x), cpp_sum(times)));
      switch (TYPEOF(x)){
      case LGLSXP:
      case INTSXP: {
        int *p_x = INTEGER(x);
        int *p_out = INTEGER(out);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          for (int j = 0; j < p_times[i]; ++j, ++k) p_out[k] = p_x[i];
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case REALSXP: {
        double *p_x = REAL(x);
        double *p_out = REAL(out);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          for (int j = 0; j < p_times[i]; ++j, ++k) p_out[k] = p_x[i];
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case STRSXP: {
        const SEXP *p_x = STRING_PTR_RO(x);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          for (int j = 0; j < p_times[i]; ++j, ++k) SET_STRING_ELT(out, k, p_x[i]);
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case CPLXSXP: {
        const Rcomplex *p_x = COMPLEX(x);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          for (int j = 0; j < p_times[i]; ++j, ++k) SET_COMPLEX_ELT(out, k, p_x[i]);
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case VECSXP: {
        const SEXP *p_x = LIST_PTR_RO(x);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          for (int j = 0; j < p_times[i]; ++j, ++k) SET_VECTOR_ELT(out, k, p_x[i]);
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      default: {
        SEXP out = SHIELD(base_rep(x, times));
        YIELD(3);
        return out;
      }
      }
    } else {
      SEXP out = SHIELD(base_rep(x, times));
      YIELD(2);
      return out;
    }
  }
}

[[cpp11::register]]
SEXP cpp_rep_each(SEXP x, SEXP each){
  int32_t NP = 0;
  SHIELD(each = cast<r_integer_t>(each, R_NilValue)); ++NP;
  if (Rf_length(each) == 1){
    if (INTEGER(each)[0] == 1){
      YIELD(NP);
      return x;
    }
    SHIELD(each = cpp_rep_len(each, vector_length(x))); ++NP;
  }
  SEXP out = SHIELD(cpp_rep(x, each)); ++NP;
  YIELD(NP);
  return out;
}

R_xlen_t length_common(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("x` must be a list");
  }

  R_xlen_t n = Rf_xlength(x);

  const SEXP *p_x = LIST_PTR_RO(x);
  R_xlen_t out = 0;

  for (R_xlen_t i = 0; i < n; ++i){
    if (is_null(p_x[i])) continue;

    if (vector_length(p_x[i]) == 0){
      out = 0;
      break;
    }
    out = std::max(out, vector_length(p_x[i]));
  }
  return out;
}

// Recycle elements of a list `x`

void recycle_in_place(SEXP x, R_xlen_t n){

  int xn = Rf_length(x);
  const SEXP *p_x = LIST_PTR_RO(x);

  for (int i = 0; i < xn; ++i){
    if (!is_null(p_x[i])){
      SET_VECTOR_ELT(x, i, cpp_rep_len(p_x[i], n));
    }
  }
}

[[cpp11::register]]
SEXP cpp_recycle(SEXP x, SEXP length){

  SEXP out = SHIELD(cpp_drop_null(x, true));

  int n = 0;

  if (!is_null(length)){
    n = Rf_asInteger(length);
    if (n < 0){
      Rf_error("Recycled `length` must be >= 0 in %s", __func__);
    }
  } else {
    n = length_common(out);
  }

  recycle_in_place(out, n);

  YIELD(1);
  return out;
}

// Fast unique that can be used in C code
// Doesn't return unique df rows
SEXP cpp_unique(SEXP x, bool names){

  int32_t NP = 0;

  bool is_simple = cheapr_is_simple_atomic_vec(x);

  if (is_compact_seq(x)){
    return x;
  } else if (is_simple && Rf_length(x) < 10000){
    SEXP dup = SHIELD(Rf_duplicated(x, FALSE)); ++NP;
    SEXP unique_locs = SHIELD(cpp_which_(dup, true)); ++NP;
    if (Rf_length(unique_locs) == Rf_length(x)){
      YIELD(NP);
      return x;
    } else {
      SEXP out = SHIELD(sset_vec(x, unique_locs, false)); ++NP;
      Rf_copyMostAttrib(x, out);
      if (names){
        SEXP names = SHIELD(get_names(x)); ++NP;
        SHIELD(names = sset_vec(names, unique_locs, false)); ++NP;
        set_names(out, names);
      }
      YIELD(NP);
      return out;
    }
  } else if (is_simple){
    SEXP out = SHIELD(cheapr_fast_unique(x)); ++NP;
    if (names){
      SEXP names = SHIELD(get_names(x)); ++NP;
      SHIELD(names = cheapr_fast_unique(names)); ++NP;
      set_names(out, names);
    }
    YIELD(NP);
    return out;
  } else {
    SEXP out = SHIELD(cpp11::package("base")["unique"](x)); ++NP;
    if (names){
      SEXP names = SHIELD(cpp11::package("base")["names"](x)); ++NP;
      SHIELD(names = cheapr_fast_unique(names)); ++NP;
      SHIELD(out = cpp11::package("base")["names<-"](out, names)); ++NP;
      YIELD(NP);
      return out;
    } else {
      YIELD(NP);
      return out;
    }
  }
}

[[cpp11::register]]
SEXP cpp_setdiff(SEXP x, SEXP y, bool unique){
  if (unique){
    SHIELD(x = cpp_unique(x, true));
  } else {
    SHIELD(x);
  }
  SEXP matches = SHIELD(match(y, x, NA_INTEGER));
  SEXP locs = SHIELD(cpp_which_na(matches));
  if (Rf_xlength(locs) == Rf_xlength(x)){
    YIELD(3);
    return x;
  } else {
    SEXP out = SHIELD(cpp_sset(x, locs, false));
    Rf_copyMostAttrib(x, out);
    YIELD(4);
    return out;
  }
}

[[cpp11::register]]
SEXP cpp_intersect(SEXP x, SEXP y, bool unique){
  if (unique){
    SHIELD(x = cpp_unique(x, true));
  } else {
    SHIELD(x);
  }
  SEXP matches = SHIELD(match(y, x, NA_INTEGER));
  SEXP locs = SHIELD(cpp_which_not_na(matches));
  if (Rf_xlength(locs) == Rf_xlength(x)){
    YIELD(3);
    return x;
  } else {
    SEXP out = SHIELD(cpp_sset(x, locs, false));
    Rf_copyMostAttrib(x, out);
    YIELD(4);
    return out;
  }
}

[[cpp11::register]]
SEXP cpp_na_init(SEXP x, int n){
  SEXP ptype = SHIELD(get_ptype(x));
  SEXP out = SHIELD(cpp_rep_len(ptype, n));
  YIELD(2);
  return out;
}

SEXP get_ptype(SEXP x){
  return slice_loc(x, 0);
}

// Prototypes of data frame

SEXP get_ptypes(SEXP x){
  int n = Rf_length(x);
  SEXP out = SHIELD(new_vec(VECSXP, n));

  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, get_ptype(VECTOR_ELT(x, i)));
  }

  SEXP names = SHIELD(get_names(x));
  set_names(out, names);

  YIELD(2);
  return out;
}

// Helper to turn a factor into a character vec
SEXP factor_as_character(SEXP x){
  return sset_vec(Rf_getAttrib(x, R_LevelsSymbol), x, true);
}

// Helper to turn a character vec into a factor, given levels
SEXP character_as_factor(SEXP x, SEXP levels){

  if (TYPEOF(x) != STRSXP){
    Rf_error("`x` must be a character vector in %s", __func__);
  }

  int32_t NP = 0;

  SEXP out = SHIELD(match(levels, x, NA_INTEGER)); ++NP;
  SEXP cls = SHIELD(make_utf8_str("factor")); ++NP;
  Rf_setAttrib(out, R_LevelsSymbol, levels);
  Rf_classgets(out, cls);
  YIELD(NP);
  return out;
}

// Helper to concatenate plain lists
// It uses top-level names if and only if x[[i]] is not a list

[[cpp11::register]]
SEXP cpp_list_c(SEXP x){
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = LIST_PTR_RO(x);

  R_xlen_t out_size = 0;
  for (R_xlen_t i = 0; i < n; ++i){
    out_size += (TYPEOF(p_x[i]) == VECSXP ? Rf_xlength(p_x[i]) : 1);
  }

  SEXP x_names = SHIELD(get_names(x)); ++NP;
  bool x_has_names = !is_null(x_names);

  R_xlen_t k = 0;
  SEXP out;

  out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
  SEXP container_list = SHIELD(new_vec(VECSXP, 1)); ++NP;
  set_names(container_list, R_BlankScalarString);

  SEXP names;
  PROTECT_INDEX nm_idx;
  R_ProtectWithIndex(names = R_NilValue, &nm_idx); ++NP;

  SEXP out_names = SHIELD(new_vec(STRSXP, out_size)); ++NP;
  bool any_names = false;

  R_xlen_t m;

  for (R_xlen_t i = 0; i < n; ++i){
    const SEXP *p_temp;

    if (TYPEOF(p_x[i]) == VECSXP){
      p_temp = LIST_PTR_RO(p_x[i]);
      R_Reprotect(names = get_names(p_x[i]), nm_idx);
      m = Rf_xlength(p_x[i]);
    } else {
      SET_VECTOR_ELT(container_list, 0, p_x[i]);
      if (x_has_names){
        R_Reprotect(names = as_r_scalar(STRING_ELT(x_names, i)), nm_idx);
      } else {
        R_Reprotect(names = R_NilValue, nm_idx);
      }
      p_temp = LIST_PTR_RO(container_list);
      m = 1;
    }

    any_names = any_names || !is_null(names);
    if (!is_null(names)){
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
    set_names(out, out_names);
  }
  YIELD(NP);
  return out;
}

// Concatenate a list of data frames

SEXP cpp_df_c(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of data frames");
  }

  int n_frames = Rf_length(x);

  if (n_frames == 0) return R_NilValue;

  int32_t NP = 0;

  // Since we are casting potential non-dfs to dfs we need a new list
  SEXP frames = SHIELD(init<r_list_t>(n_frames, false)); ++NP;

  const SEXP *p_x = LIST_PTR_RO(x);
  const SEXP *p_frames = LIST_PTR_RO(frames);

  SEXP df;
  PROTECT_INDEX df_idx;
  R_ProtectWithIndex(df = cast<r_data_frame_t>(p_x[0], R_NilValue), &df_idx); ++NP;
  SET_VECTOR_ELT(frames, 0, df);

  SEXP new_names, df_names, ptype_names;
  PROTECT_INDEX new_names_idx, df_names_idx, ptype_names_idx;

  R_ProtectWithIndex(new_names = R_NilValue, &new_names_idx); ++NP;
  R_ProtectWithIndex(df_names = R_NilValue, &df_names_idx); ++NP;
  R_ProtectWithIndex(ptype_names = get_names(df), &ptype_names_idx); ++NP;

  int n_cols = Rf_length(ptype_names);

  // We do 2 passes
  // 1st pass: Cast inputs to dfs and construct prototypes
  // 2nd pass: initialise data frame vecs (if need be) and combine simultaneously

  int out_size = df_nrow(df);

  bool na_padding = false;

  for (int i = 1; i < n_frames; ++i){
    R_Reprotect(df = cast<r_data_frame_t>(p_x[i], R_NilValue), df_idx);
    R_Reprotect(df_names = get_names(df), df_names_idx);

    R_Reprotect(new_names = cpp_setdiff(
      df_names, ptype_names, false
    ), new_names_idx);

    // Adjust prototype names
    if (Rf_length(new_names) > 0){
      na_padding = true;
      R_Reprotect(ptype_names = r_combine(ptype_names, new_names), ptype_names_idx);
    }
    na_padding = na_padding || Rf_length(df) != n_cols;
    out_size += df_nrow(df);
    SET_VECTOR_ELT(frames, i, df);
  }

  n_cols = Rf_length(ptype_names);

  SEXP vec;
  PROTECT_INDEX vec_idx;
  R_ProtectWithIndex(vec = R_NilValue, &vec_idx); ++NP;

  SEXP out = SHIELD(new_vec(VECSXP, n_cols)); ++NP;
  SEXP vectors = SHIELD(new_vec(VECSXP, n_frames)); ++NP;

  const SEXP *p_ptype_names = STRING_PTR_RO(ptype_names);

  if (na_padding){
    // Get archetype of each col
    SEXP vec_archetypes = SHIELD(new_vec(VECSXP, n_cols)); ++NP;
    const SEXP *p_vec_archetypes = LIST_PTR_RO(vec_archetypes);
    for (int j = 0; j < n_cols; ++j){
      for (int i = 0; i < n_frames; ++i){
        SET_VECTOR_ELT(vectors, i, get_list_element(p_frames[i], p_ptype_names[j]));
      }
      SET_VECTOR_ELT(vec_archetypes, j, cpp_common_template(vectors));
    }
    // Combine all cols
    for (int j = 0; j < n_cols; ++j){
      for (int i = 0; i < n_frames; ++i){
        vec = get_list_element(p_frames[i], p_ptype_names[j]);
        if (is_null(vec)){
          R_Reprotect(vec = cpp_na_init(p_vec_archetypes[j], df_nrow(p_frames[i])), vec_idx);
        }
        SET_VECTOR_ELT(vectors, i, vec);
      }
      SET_VECTOR_ELT(out, j, combine_internal(vectors, out_size, p_vec_archetypes[j]));
    }
  } else {
    // Simply combine rows with no NA padding
    for (int j = 0; j < n_cols; ++j){
      for (int i = 0; i < n_frames; ++i){
        SET_VECTOR_ELT(vectors, i, get_list_element(p_frames[i], p_ptype_names[j]));
      }
      SET_VECTOR_ELT(out, j, cpp_c(vectors));
    }
  }

  set_list_as_df(out);
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(out_size));
  set_names(out, ptype_names);
  SHIELD(out = rebuild(out, VECTOR_ELT(frames, 0), false)); ++NP;
  YIELD(NP);
  return out;
}

// Helper to concatenate data frames by col

[[cpp11::register]]
SEXP cpp_df_col_c(SEXP x, bool recycle, bool name_repair){
  int32_t NP = 0;

  // Important to recycle first to avoid incorrect size calculations

  if (recycle){
    SHIELD(x = cpp_recycle(x, R_NilValue)); ++NP;
  }

  int n = Rf_length(x);
  const SEXP *p_x = LIST_PTR_RO(x);

  int out_ncols = 0;

  SEXP container_list = SHIELD(new_vec(VECSXP, 1)); ++NP;
  set_names(container_list, R_BlankScalarString);

  std::vector<const SEXP *> df_pointers(n);

  for (int i = 0; i < n; ++i){
    if (is_df(p_x[i])){
      df_pointers[i] = LIST_PTR_RO(p_x[i]);
      out_ncols += Rf_length(p_x[i]);
    } else {
      df_pointers[i] = LIST_PTR_RO(container_list);
      ++out_ncols;
    }
  }

  SEXP x_names = SHIELD(get_names(x)); ++NP;
  bool x_has_names = !is_null(x_names);

  SEXP out = SHIELD(new_vec(VECSXP, out_ncols)); ++NP;

  SEXP names;
  PROTECT_INDEX nm_idx;
  R_ProtectWithIndex(names = R_NilValue, &nm_idx); ++NP;

  SEXP out_names = SHIELD(new_vec(STRSXP, out_ncols)); ++NP;
  bool any_names = false;

  int m;
  int k = 0;

  for (int i = 0; i < n; ++i){

    const SEXP *p_temp = df_pointers[i];

    if (is_df(p_x[i])){
      R_Reprotect(names = get_names(p_x[i]), nm_idx);
      m = Rf_length(p_x[i]);
    } else {
      SET_VECTOR_ELT(container_list, 0, p_x[i]);
      if (x_has_names){
        R_Reprotect(names = as_r_scalar(STRING_ELT(x_names, i)), nm_idx);
      } else {
        R_Reprotect(names = R_NilValue, nm_idx);
      }
      m = 1;
    }

    any_names = any_names || !is_null(names);
    if (!is_null(names)){
      for (int j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
        SET_STRING_ELT(out_names, k, STRING_ELT(names, j));
      }
    } else {
      for (int j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
      }
    }
  }
  if (any_names){
    set_names(out, out_names);
  }

  SEXP r_nrows = SHIELD(R_NilValue); ++NP;
  if (Rf_length(out) == 0 && Rf_length(x) != 0){
    SHIELD(r_nrows = as_r_scalar<int>(vector_length(VECTOR_ELT(x, 0)))); ++NP;
  }

  SHIELD(out = cpp_new_df(out, r_nrows, false, name_repair)); ++NP;

  if (Rf_length(x) != 0 && is_df(VECTOR_ELT(x, 0))){
    SHIELD(out = rebuild(out, VECTOR_ELT(x, 0), false)); ++NP;
  }
  YIELD(NP);
  return out;
}

// define CHEAPR_COMBINE_NAMES
// if (any_names){
//   R_xlen_t ii = 0;
//   SHIELD(combined_names = new_vec(STRSXP, out_size)); ++NP;
//   for (int i = 0; i < n; ++i){
//     R_Reprotect(vec_names = get_vec_names(p_x[i]), vec_names_idx);
//     const SEXP *p_vec_names = STRING_PTR_RO(vec_names);
//     if (is_null(vec_names)){
//       ii += vector_length(p_x[i]);
//     } else {
//       for (R_xlen_t j = 0; j < Rf_xlength(vec_names); ++j, ++ii){
//         SET_STRING_ELT(combined_names, ii, p_vec_names[j]);
//       }
//     }
//   }
// }

// Combine vectors given the following args:
// * A list of vectors
// * The final size of the combined vector
// * A template vector (of length 0) whose type is the common type between all vectors
// This allows us to eliminate a few loops
SEXP combine_internal(SEXP x, const R_xlen_t out_size, SEXP vec_template){

  if (vector_length(vec_template) != 0){
    Rf_error("`vec_template` must be of length 0");
  }

  int n = Rf_length(x);
  const SEXP *p_x = LIST_PTR_RO(x);

  int32_t NP = 0;

  R_xlen_t k = 0;
  R_xlen_t m = 0;
  SEXP out = R_NilValue;

  SEXP vec;
  PROTECT_INDEX vec_idx;
  R_ProtectWithIndex(vec = R_NilValue, &vec_idx); ++NP;

  r_type common = get_r_type(vec_template);
  SEXP combined_names = R_NilValue;
  bool has_top_level_names = vec_has_names(x);

  // If in the future it is decided that inner names are to be kept
  // Use this condition within the inner loop and then use the macro
  // defined above combine_internal
  // any_names = !has_top_level_names && (any_names || has_names(vec));

  switch (common){
  case r_null: {
    break;
  }
  case r_lgl: {

    SHIELD(out = init<r_logical_t>(out_size, false)); ++NP;

    int* RESTRICT p_out = INTEGER(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_logical_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(INTEGER_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_int: {

    SHIELD(out = init<r_integer_t>(out_size, false)); ++NP;

    int* RESTRICT p_out = INTEGER(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_integer_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(INTEGER_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_int64: {

    SHIELD(out = init<r_integer64_t>(out_size, false)); ++NP;

    int64_t* RESTRICT p_out = INTEGER64_PTR(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_integer64_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(INTEGER64_PTR_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_dbl: {

    SHIELD(out = init<r_numeric_t>(out_size, false)); ++NP;

    double* RESTRICT p_out = REAL(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_numeric_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(REAL_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_chr: {

    SHIELD(out = init<r_character_t>(out_size, false)); ++NP;

    for (int i = 0; i < n; ++i){
      R_Reprotect(vec = cast<r_character_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);

      const SEXP *p_vec = STRING_PTR_RO(vec);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_STRING_ELT(out, k, p_vec[j]);
      }
    }
    break;
  }
  case r_cplx: {

    SHIELD(out = init<r_complex_t>(out_size, false)); ++NP;

    Rcomplex *p_out = COMPLEX(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_complex_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(COMPLEX_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_raw: {

    SHIELD(out = init<r_raw_t>(out_size, false)); ++NP;

    Rbyte *p_out = RAW(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_raw_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(RAW_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_list: {

    SHIELD(out = init<r_list_t>(out_size, false)); ++NP;

    for (int i = 0; i < n; ++i){
      R_Reprotect(vec = cast<r_list_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);

      const SEXP *p_vec = LIST_PTR_RO(vec);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_vec[j]);
      }
    }
    break;
  }
  case r_fct: {

    SHIELD(out = cpp_rep_len(vec_template, out_size)); ++NP;

    int* RESTRICT p_out = INTEGER(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_factor_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(INTEGER_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_date: {

    if (TYPEOF(vec_template) == INTSXP){
    SHIELD(out = init<r_integer_t>(out_size, false)); ++NP;
    Rf_classgets(out, make_utf8_str("Date"));

    int* RESTRICT p_out = INTEGER(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_date_t>(p_x[i], vec_template), vec_idx);
      R_Reprotect(vec = coerce_vec(vec, INTSXP), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(INTEGER_RO(vec), m, &p_out[k]);
    }
  } else {
    SHIELD(out = init<r_numeric_t>(out_size, false)); ++NP;
    Rf_classgets(out, make_utf8_str("Date"));

    double* RESTRICT p_out = REAL(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_date_t>(p_x[i], vec_template), vec_idx);
      R_Reprotect(vec = coerce_vec(vec, REALSXP), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(REAL_RO(vec), m, &p_out[k]);
    }
  }
  break;
  }
  case r_pxct: {

    SHIELD(out = init<r_posixt_t>(out_size, false)); ++NP;
    SEXP out_tzone = SHIELD(Rf_getAttrib(vec_template, install_utf8("tzone"))); ++NP;
    Rf_setAttrib(out, install_utf8("tzone"), out_tzone);

    double* RESTRICT p_out = REAL(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_posixt_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(REAL_RO(vec), m, &p_out[k]);
    }
    break;
  }
  case r_df: {
    SHIELD(out = cpp_df_c(x)); ++NP;
    break;
  }
  case r_unk: {
    SEXP call = SHIELD(coerce_vec(x, LISTSXP)); ++NP;
    SHIELD(call = Rf_lcons(install_utf8("c"), call)); ++NP;
    SHIELD(out = Rf_eval(call, R_GetCurrentEnv())); ++NP;
    break;
  }
  default: {
    YIELD(NP);
    Rf_error("Don't know how to combine elements");
  }
  }

  if (has_top_level_names){
    SEXP top_level_names = SHIELD(get_vec_names(x)); ++NP;
    SEXP name_sizes = SHIELD(cpp_lengths(x, false)); ++NP;
    SHIELD(combined_names = cpp_rep(top_level_names, name_sizes)); ++NP;
  }
  SHIELD(out = set_vec_names(out, combined_names)); ++NP;
  YIELD(NP);
  return out;
}

// `c()` but no concatenation of names

[[cpp11::register]]
SEXP cpp_c(SEXP x){

  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list of vectors");
  }

  int n = Rf_length(x);

  // Cast all objects to common type
  const SEXP *p_x = LIST_PTR_RO(x);

  // Figure out final size
  R_xlen_t out_size = 0;
  for (int i = 0; i < n; ++i) out_size += vector_length(p_x[i]);

  // 'vec_template' here acts as a template for the final result
  SEXP vec_template = SHIELD(cpp_common_template(x));
  SEXP out = SHIELD(combine_internal(x, out_size, vec_template));
  YIELD(2);
  return out;
}
