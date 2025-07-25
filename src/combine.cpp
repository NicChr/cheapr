#include "cheapr.h"
#include <vector>

// static cpp11::writable::integers CHEAPR_ZERO(1);
// void constants_init(DllInfo* dll){
//   CHEAPR_ZERO[0] = 0;
// }

[[cpp11::register]]
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
    const SEXP *p_x = VECTOR_PTR_RO(x);
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
  } else if (is_simple_vec2(x)){

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
      const int64_t *p_x = INTEGER64_RO_PTR(x);
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
      const SEXP *p_x = VECTOR_PTR_RO(x);
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

  R_xlen_t n = vec_length(x);

  R_xlen_t out_size;
  R_xlen_t n_times = Rf_xlength(times);
  SHIELD(times = coerce_vector(times, INTSXP));

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

    } else if (is_simple_vec2(x)){
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
        const SEXP *p_x = VECTOR_PTR_RO(x);
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
  SHIELD(each = coerce_vector(each, INTSXP)); ++NP;
  if (Rf_length(each) == 1){
    if (INTEGER(each)[0] == 1){
      Rf_unprotect(NP);
      return x;
    }
    SHIELD(each = cpp_rep_len(each, vec_length(x))); ++NP;
  }
  SEXP out = SHIELD(cpp_rep(x, each)); ++NP;
  YIELD(NP);
  return out;
}

// Recycle elements of a list `x`

[[cpp11::register]]
SEXP cpp_recycle(SEXP x, SEXP length){

  int32_t NP = 0;
  int xn = Rf_length(x);
  int n_objs = xn - static_cast<int>(null_count(x));

  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = SHIELD(new_vec(VECSXP, n_objs)); ++NP;
  SEXP sizes = SHIELD(cpp_lengths(x, false)); ++NP;
  const int* RESTRICT p_sizes = INTEGER_RO(sizes);
  bool has_length = !is_null(length);
  SHIELD(length = coerce_vec(length, INTSXP)); ++NP;

  int n = 0;

  if (!has_length){
    if (n_objs > 0){
      n = -1;
      // We calculate `max(sizes)`
      // We won't have any NA in sizes so no NA checking is needed
      for (int i = 0; i < xn; ++i){
        if (is_null(p_x[i])){
          continue;
        } else if (p_sizes[i] == 0){
          n = 0;
          break;
        } else {
          n = std::max(n, p_sizes[i]);
        }
      }
    }
  } else {
    n = Rf_asInteger(length);
  }

  n = std::max(n, 0);

  int k = 0;

  SEXP names = SHIELD(get_names(x)); ++NP;
  SEXP new_names = R_NilValue;

  if (xn == n_objs || is_null(names)){
    for (int i = 0; i < xn; ++i){
      if (is_null(p_x[i])) continue;
      SET_VECTOR_ELT(out, k++, cpp_rep_len(p_x[i], n));
    }
    set_names(out, names);
  } else {
    new_names = SHIELD(new_vec(STRSXP, n_objs)); ++NP;
    const SEXP *p_names = STRING_PTR_RO(names);
    for (int i = 0; i < xn; ++i){
      if (is_null(p_x[i])) continue;
      SET_VECTOR_ELT(out, k, cpp_rep_len(p_x[i], n));
      SET_STRING_ELT(new_names, k++, p_names[i]);
    }
    set_names(out, new_names);
  }

  YIELD(NP);
  return out;
}

// Fast unique that can be used in C code
// Doesn't return unique df rows
SEXP cpp_unique(SEXP x, bool names){

  int32_t NP = 0;

  bool is_simple = is_simple_atomic_vec(x);

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
  SEXP matches = SHIELD(Rf_match(y, x, NA_INTEGER));
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
  SEXP matches = SHIELD(Rf_match(y, x, NA_INTEGER));
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
// SEXP cpp_setdiff(SEXP x, SEXP y, bool dups){
//
//   bool is_simple = is_simple_atomic_vec(x) && is_simple_atomic_vec(y);
//
//   int32_t NP = 0;
//   if (is_simple){
//     if (!dups){
//       SHIELD(x = cpp_unique(x)); ++NP;
//     }
//     SEXP matches;
//     if (Rf_length(x) < 10000){
//       matches = SHIELD(Rf_match(y, x, NA_INTEGER)); ++NP;
//     } else {
//       matches = SHIELD(cheapr_fast_match(x, y)); ++NP;
//     }
//     SEXP locs = SHIELD(cpp_which_na(matches)); ++NP;
//     SEXP out = SHIELD(sset_vec(x, locs, false)); ++NP;
//
//     Rf_copyMostAttrib(x, out);
//     YIELD(NP);
//     return out;
//   } else {
//     if (!dups){
//       SHIELD(x = cheapr_fast_unique(x)); ++NP;
//     }
//     SEXP matches = SHIELD(cheapr_fast_match(x, y)); ++NP;
//     SEXP locs = SHIELD(cpp_which_na(matches)); ++NP;
//     SEXP out = SHIELD(cheapr_sset(x, locs)); ++NP;
//     YIELD(NP);
//     return out;
//   }
// }

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

SEXP cpp_factor_as_character(SEXP x){
  return sset_vec(Rf_getAttrib(x, R_LevelsSymbol), x, true);
}

[[cpp11::register]]
SEXP cpp_combine_levels(SEXP x){
  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of factors in %s", __func__);
  }
  int n = Rf_length(x);
  SEXP levels = SHIELD(new_vec(VECSXP, n));
  const SEXP *p_x = VECTOR_PTR_RO(x);

  SEXP fct_levels;

  PROTECT_INDEX fct_idx;
  R_ProtectWithIndex(fct_levels = R_NilValue, &fct_idx);

  for (int i = 0; i < n; ++i){
    if (Rf_isFactor(p_x[i])){
      fct_levels = Rf_getAttrib(p_x[i], R_LevelsSymbol);
    } else {
      R_Reprotect(fct_levels = base_as_character(p_x[i]), fct_idx);
    }
    SET_VECTOR_ELT(levels, i, fct_levels);
  }
  SEXP out = SHIELD(cpp_c(levels));
  SHIELD(out = cpp_unique(out, false));
  YIELD(4);
  return out;
}

SEXP cpp_combine_factors(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of factors in %s", __func__);
  }

  int n = Rf_length(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

  SEXP levels = SHIELD(cpp_combine_levels(x));
  SEXP chars = SHIELD(new_vec(VECSXP, n));
  SEXP char_vec;

  PROTECT_INDEX char_vec_idx;
  R_ProtectWithIndex(char_vec = R_NilValue, &char_vec_idx);

  for (int i = 0; i < n; ++i){
    if (Rf_isFactor(p_x[i])){
      R_Reprotect(char_vec = cpp_factor_as_character(p_x[i]), char_vec_idx);
    } else {
      R_Reprotect(char_vec = base_as_character(p_x[i]), char_vec_idx);
    }
    SET_VECTOR_ELT(chars, i, char_vec);
  }

  R_Reprotect(char_vec = cpp_c(chars), char_vec_idx);
  SEXP out = SHIELD(cheapr_factor(char_vec,
    cpp11::named_arg("levels") = levels,
    cpp11::named_arg("na_exclude") = !cpp_any_na(levels, false)
  ));
  YIELD(4);
  return out;
}

// Helper to concatenate plain lists
// It uses top-level names if and only if x[[i]] is not a list

[[cpp11::register]]
SEXP cpp_list_c(SEXP x){
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

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
      p_temp = VECTOR_PTR_RO(p_x[i]);
      R_Reprotect(names = get_names(p_x[i]), nm_idx);
      m = Rf_xlength(p_x[i]);
    } else {
      SET_VECTOR_ELT(container_list, 0, p_x[i]);
      if (x_has_names){
        R_Reprotect(names = Rf_ScalarString(STRING_ELT(x_names, i)), nm_idx);
      } else {
        R_Reprotect(names = R_NilValue, nm_idx);
      }
      p_temp = VECTOR_PTR_RO(container_list);
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

// Helper to concatenate 2 lists

SEXP list_c2(SEXP x, SEXP y){
  SEXP temp = SHIELD(new_vec(VECSXP, 2));

  SET_VECTOR_ELT(temp, 0, x);
  SET_VECTOR_ELT(temp, 1, y);
  SEXP out = SHIELD(cpp_list_c(temp));

  YIELD(2);
  return out;
}

// Helper to concatenate two vectors
SEXP c2(SEXP x, SEXP y){
  SEXP temp = SHIELD(new_vec(VECSXP, 2));

  SET_VECTOR_ELT(temp, 0, x);
  SET_VECTOR_ELT(temp, 1, y);
  SEXP out = SHIELD(cpp_c(temp));

  YIELD(2);
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
  const SEXP *p_x = VECTOR_PTR_RO(x);

  SEXP df = p_x[0];

  SEXP names;
  PROTECT_INDEX names_idx;
  R_ProtectWithIndex(names = get_names(df), &names_idx); ++NP;

  if (!is_df(df)){
    YIELD(NP); Rf_error("Can't combine data frames with non data frames");
  }

  SEXP df_template = df;

  SEXP frames = SHIELD(new_vec(VECSXP, n_frames)); ++NP;
  SET_VECTOR_ELT(frames, 0, df_template);

  SEXP ptypes, new_names, new_ptypes, new_cols,
  temp_list, df_names, ptype_names;
  PROTECT_INDEX
  new_names_idx, ptypes_idx, new_ptypes_idx, new_cols_idx,
  temp_list_idx, df_names_idx, ptype_names_idx;

  R_ProtectWithIndex(ptypes = get_ptypes(df), &ptypes_idx); ++NP;
  R_ProtectWithIndex(new_names = R_NilValue, &new_names_idx); ++NP;
  R_ProtectWithIndex(new_ptypes = R_NilValue, &new_ptypes_idx); ++NP;
  R_ProtectWithIndex(new_cols = R_NilValue, &new_cols_idx); ++NP;
  R_ProtectWithIndex(temp_list = new_vec(VECSXP, 2), &temp_list_idx); ++NP;
  R_ProtectWithIndex(df_names = R_NilValue, &df_names_idx); ++NP;
  R_ProtectWithIndex(ptype_names = R_NilValue, &ptype_names_idx); ++NP;

  // We do 2 passes
  // 1st pass: Check inputs are dfs and construct prototype list
  // 2nd pass: initialise data frame vecs (if need be) and combine simultaneously

  int out_size = df_nrow(df);

  for (int i = 1; i < n_frames; ++i){
    df = p_x[i];

    if (!is_df(df)){
      YIELD(NP); Rf_error("Can't combine data frames with non data frames");
    }

    R_Reprotect(df_names = get_names(df), df_names_idx);
    R_Reprotect(ptype_names = get_names(ptypes), ptype_names_idx);

    R_Reprotect(new_names = cpp_setdiff(
      df_names, ptype_names, false
    ), new_names_idx);

    // Adjust prototype list
    if (Rf_length(new_names) > 0){
      R_Reprotect(new_cols = cpp_df_select(df, new_names), new_cols_idx);
      R_Reprotect(new_ptypes = get_ptypes(new_cols), new_ptypes_idx);
      SET_VECTOR_ELT(temp_list, 0, ptypes);
      SET_VECTOR_ELT(temp_list, 1, new_ptypes);
      R_Reprotect(ptypes = cpp_list_c(temp_list), ptypes_idx);
      SET_VECTOR_ELT(temp_list, 0, names);
      SET_VECTOR_ELT(temp_list, 1, new_names);
      R_Reprotect(names = cpp_c(temp_list), names_idx);
      set_names(ptypes, names);
    }
    out_size += df_nrow(df);
  }

  int n_cols = Rf_length(names);

  SEXP vec, name_locs;
  PROTECT_INDEX vec_idx, name_locs_idx;
  R_ProtectWithIndex(vec = R_NilValue, &vec_idx); ++NP;
  R_ProtectWithIndex(name_locs = R_NilValue, &name_locs_idx); ++NP;

  SEXP out = SHIELD(new_vec(VECSXP, n_cols)); ++NP;
  SEXP vectors = SHIELD(new_vec(VECSXP, n_frames)); ++NP;

  const SEXP *p_ptypes = VECTOR_PTR_RO(ptypes);
  const SEXP *p_names = VECTOR_PTR_RO(names);

  for (int j = 0; j < n_cols; ++j){
    for (int i = 0; i < n_frames; ++i){
      vec = get_list_element(p_x[i], p_names[j]);

      if (is_null(vec)){
        R_Reprotect(vec = cpp_na_init(p_ptypes[j], df_nrow(p_x[i])), vec_idx);
      }
      SET_VECTOR_ELT(vectors, i, vec);
    }
    SET_VECTOR_ELT(out, j, cpp_c(vectors));
  }
  set_list_as_df(out);
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(out_size));
  set_names(out, names);
  SHIELD(out = rebuild(out, df_template, false)); ++NP;
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
  const SEXP *p_x = VECTOR_PTR_RO(x);

  int out_ncols = 0;

  SEXP container_list = SHIELD(new_vec(VECSXP, 1)); ++NP;
  set_names(container_list, R_BlankScalarString);

  std::vector<const SEXP *> df_pointers(n);

  for (int i = 0; i < n; ++i){
    if (is_df(p_x[i])){
      df_pointers[i] = VECTOR_PTR_RO(p_x[i]);
      out_ncols += Rf_length(p_x[i]);
    } else {
      df_pointers[i] = VECTOR_PTR_RO(container_list);
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
        R_Reprotect(names = Rf_ScalarString(STRING_ELT(x_names, i)), nm_idx);
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
    SHIELD(r_nrows = Rf_ScalarInteger(vec_length(VECTOR_ELT(x, 0)))); ++NP;
  }

  SHIELD(out = cpp_new_df(out, r_nrows, false, name_repair)); ++NP;

  if (Rf_length(x) != 0 && is_df(VECTOR_ELT(x, 0))){
    SHIELD(out = rebuild(out, VECTOR_ELT(x, 0), false)); ++NP;
  }
  YIELD(NP);
  return out;
}

// `c()` but no concatenation of names
[[cpp11::register]]
SEXP cpp_c(SEXP x){
  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of vectors");
  }
  int32_t NP = 0;
  int n = Rf_length(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);
  if (n == 1){
    return p_x[0];
  }

  int vector_type = NILSXP;
  R_xlen_t out_size = 0;

  bool is_factor = false, is_frame = false, is_classed = false,
    is_simple = false, is_date = false, is_datetime = false;


  // We use the tz info of the first datetime vec to copy to final result
  int first_datetime = INTEGER_MAX;

  SEXP elem = R_NilValue;

  for (int i = 0; i < n; ++i){
    elem = p_x[i];
    vector_type = std::max(vector_type, TYPEOF(elem));
    out_size += Rf_xlength(elem);
    is_factor = is_factor || Rf_isFactor(elem);
    is_simple = is_simple || is_simple_atomic_vec(elem);
    is_date = is_date || Rf_inherits(elem, "Date");
    is_datetime = is_datetime || Rf_inherits(elem, "POSIXct");
    if (is_datetime){
      first_datetime = std::min(i, first_datetime);
    }
    is_frame = is_frame || is_df(elem);
    is_classed = is_classed || Rf_isObject(elem);
  }

  // Date vectors can be ints but datetimes can't
  bool is_date2 = is_date && !is_datetime && (vector_type == INTSXP || vector_type == REALSXP);
  bool is_datetime2 = !is_date && is_datetime && vector_type == REALSXP;

  if (is_frame){
    return cpp_df_c(x);
  }

  if (is_factor){
    return cpp_combine_factors(x);
  }

  if (is_classed && !(is_date2 || is_datetime2)){
    SEXP call = SHIELD(coerce_vec(x, LISTSXP));
    SHIELD(call = Rf_lcons(install_utf8("c"), call));
    SEXP out = SHIELD(Rf_eval(call, R_GetCurrentEnv()));
    YIELD(3);
    return out;
  }

  R_xlen_t k = 0;
  R_xlen_t m = 0;
  SEXP out, temp;

  PROTECT_INDEX temp_idx;
  R_ProtectWithIndex(temp = R_NilValue, &temp_idx); ++NP;


  switch(vector_type){
  case NILSXP: {
    out = SHIELD(R_NilValue); ++NP;
    break;
  }
  case LGLSXP:
  case INTSXP: {

    out = SHIELD(new_vec(vector_type, out_size)); ++NP;
    int* RESTRICT p_out = INTEGER(out);

    for (int i = 0; i < n; ++i, k += m){
      if (TYPEOF(p_x[i]) == vector_type){
        temp = p_x[i];
      } else {
        R_Reprotect(temp = coerce_vec(p_x[i], vector_type), temp_idx);
      }
      m = Rf_xlength(temp);
      std::copy_n(INTEGER_RO(temp), m, &p_out[k]);
    }
    break;
  }
  case REALSXP: {

    out = SHIELD(new_vec(vector_type, out_size)); ++NP;
    double* RESTRICT p_out = REAL(out);

    for (int i = 0; i < n; ++i, k += m){
      if (TYPEOF(p_x[i]) == vector_type){
        temp = p_x[i];
      } else {
        R_Reprotect(temp = coerce_vec(p_x[i], vector_type), temp_idx);
      }
      m = Rf_xlength(temp);
      std::copy_n(REAL_RO(temp), m, &p_out[k]);
    }
    break;
  }
  case STRSXP: {

    out = SHIELD(new_vec(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      if (TYPEOF(p_x[i]) == vector_type){
        temp = p_x[i];
      } else {
        R_Reprotect(temp = coerce_vec(p_x[i], vector_type), temp_idx);
      }
      m = Rf_xlength(temp);
      const SEXP *p_temp = STRING_PTR_RO(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_STRING_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }
  case CPLXSXP: {

    out = SHIELD(new_vec(vector_type, out_size)); ++NP;
    Rcomplex *p_out = COMPLEX(out);

    for (int i = 0; i < n; ++i, k += m){
      if (TYPEOF(p_x[i]) == vector_type){
        temp = p_x[i];
      } else {
        R_Reprotect(temp = coerce_vec(p_x[i], vector_type), temp_idx);
      }
      m = Rf_xlength(temp);
      std::copy_n(COMPLEX_RO(temp), m, p_out + k);
    }
    break;
  }
  case VECSXP: {

    out = SHIELD(new_vec(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      if (TYPEOF(p_x[i]) == vector_type){
        temp = p_x[i];
      } else {
        R_Reprotect(temp = coerce_vec(p_x[i], vector_type), temp_idx);
      }
      m = Rf_xlength(temp);
      const SEXP *p_temp = VECTOR_PTR_RO(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }

  default: {
    SEXP c_char = SHIELD(make_utf8_str("c")); ++NP;
    out = SHIELD(base_do_call(c_char, x)); ++NP;
    break;
  }
  }
  if (is_date2){
    Rf_classgets(out, make_utf8_str("Date"));
  }
  if (is_datetime2){
    Rf_copyMostAttrib(p_x[first_datetime], out);
  }
  YIELD(NP);
  return out;
}
