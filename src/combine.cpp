#include "cheapr.h"

// static SEXP CHEAPR_ZERO = NULL;
// void constants_init(DllInfo* dll){
//   CHEAPR_ZERO = Rf_ScalarInteger(0);
//   R_PreserveObject(CHEAPR_ZERO);  // Protect from garbage collection
// }

// Helper to avoid repeated calls to R

SEXP fast_df_reconstruct(SEXP x, SEXP source){
  if (is_bare_df(source)){
    return x;
  } else if (is_bare_tbl(source)){

    // Only copy the class if source is a plain tbl

    SEXP out = SHIELD(Rf_shallow_duplicate(x));

    Rf_setAttrib(out, R_ClassSymbol, Rf_getAttrib(source, R_ClassSymbol));
    YIELD(1);
    return out;
  } else {

    // Method dispatch, users can write `reconstruct` methods which
    // this will use
    return cheapr_reconstruct(x, source);
  }
}

[[cpp11::register]]
SEXP cpp_rep_len(SEXP x, int length){
  int out_size = length;

  if (is_null(x)){
    return R_NilValue;
  } else if (is_df(x)){
    int n_cols = Rf_length(x);
    SEXP out = SHIELD(new_vec(VECSXP, n_cols));
    for (int i = 0; i < n_cols; ++i){
      SET_VECTOR_ELT(out, i, cpp_rep_len(VECTOR_ELT(x, i), length));
    }
    Rf_namesgets(out, Rf_getAttrib(x, R_NamesSymbol));
    SHIELD(out = cpp_list_as_df(out));
    Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(length));
    SHIELD(out = fast_df_reconstruct(out, x));
    YIELD(3);
    return out;
  } else if (is_simple_vec(x)){

    int size = Rf_length(x);
    int n_chunks, k, chunk_size;

    // Return x if length(x) == length
    if (out_size == size) return x;

    switch (TYPEOF(x)){
    case LGLSXP:
    case INTSXP: {
      int *p_x = INTEGER(x);
      SEXP out = SHIELD(new_vec(TYPEOF(x), out_size));
      int *p_out = INTEGER(out);

      if (size == 1){
        int val = p_x[0];
        if (val == 0){
          memset(p_out, 0, out_size * sizeof(int));
        } else {
          OMP_FOR_SIMD
          for (int i = 0; i < out_size; ++i) p_out[i] = val;
        }
      } else if (out_size > 0 && size > 0){
        n_chunks = std::ceil((static_cast<double>(out_size)) / size);
        for (int i = 0; i < n_chunks; ++i){
          k = i * size;
          chunk_size = std::min(k + size, out_size) - k;
          memcpy(&p_out[k], &p_x[0], chunk_size * sizeof(int));
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        OMP_FOR_SIMD
        for (int i = 0; i < out_size; ++i){
          p_out[i] = NA_INTEGER;
        }
      }
      Rf_copyMostAttrib(x, out);
      YIELD(1);
      return out;
    }
    case REALSXP: {
      double *p_x = REAL(x);
      SEXP out = SHIELD(new_vec(REALSXP, out_size));
      double *p_out = REAL(out);

      if (size == 1){
       double val = p_x[0];
        if (val == 0.0){
          memset(p_out, 0, out_size * sizeof(double));
        } else {
          OMP_FOR_SIMD
          for (int i = 0; i < out_size; ++i) p_out[i] = val;
        }
      } else if (out_size > 0 && size > 0){
        n_chunks = std::ceil((static_cast<double>(out_size)) / size);
        for (int i = 0; i < n_chunks; ++i){
          k = i * size;
          chunk_size = std::min(k + size, out_size) - k;
          memcpy(&p_out[k], &p_x[0], chunk_size * sizeof(double));
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        OMP_FOR_SIMD
        for (int i = 0; i < out_size; ++i){
          p_out[i] = NA_REAL;
        }
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
        for (int i = 0; i < out_size; ++i){
          SET_STRING_ELT(out, i, p_x[i % size]);
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
        Rcomplex val = p_x[0];
        for (int i = 0; i < out_size; ++i){
          SET_COMPLEX_ELT(out, i, val);
        }
      } else if (out_size > 0 && size > 0){
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
        for (int i = 0; i < out_size; ++i){
          SET_VECTOR_ELT(out, i, p_x[i % size]);
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
      return base_rep(x, cpp11::named_arg("length.out") = length);
    }
    }
  } else {
    return base_rep(x, cpp11::named_arg("length.out") = length);
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
        SHIELD(out = fast_df_reconstruct(out, x));
        YIELD(3);
        return out;
      } else {
        SEXP row_names = SHIELD(cpp_seq_len(df_nrow(x)));
        SHIELD(row_names = cpp_rep(row_names, times));
        SEXP out = SHIELD(cpp_sset(x, row_names, true));
        YIELD(4);
        return out;
      }

    } else if (is_simple_vec(x)){
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
  // SEXP temp_list = SHIELD(new_vec(VECSXP, n_times));
  // SEXP temp;
  // PROTECT_INDEX temp_idx;
  // R_ProtectWithIndex(temp = R_NilValue, &temp_idx);
  //
  // for (R_xlen_t i = 0; i < n_times; ++i){
  //   R_Reprotect(temp = slice_loc(x, i), temp_idx);
  //   SET_VECTOR_ELT(temp_list, i, cpp_rep_len(temp, p_times[i]));
  // }
  //
  // SEXP out = SHIELD(cpp_c(temp_list));
  // YIELD(4);
  // return out;
}

[[cpp11::register]]
SEXP cpp_rep_each(SEXP x, SEXP each){
  int NP = 0;
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
  SEXP out = SHIELD(cpp_drop_null(x, true));
  SEXP sizes = SHIELD(cpp_lengths(out, false));
  int *p_sizes = INTEGER(sizes);
  bool has_length = !is_null(length);
  SHIELD(length = coerce_vec(length, INTSXP));

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

  SEXP r_zero = SHIELD(Rf_ScalarInteger(0));
  if (!has_length && scalar_count(sizes, r_zero, false) > 0) n = 0;

  for (int i = 0; i < n_objs; ++i){
    SET_VECTOR_ELT(out, i, cpp_rep_len(VECTOR_ELT(out, i), n));
  }
  YIELD(4);
  return out;
}

// Fast unique that can be used in C code
// Doesn't return unique df rows
SEXP cpp_unique(SEXP x){

  bool is_simple = is_simple_atomic_vec(x);

  if (is_compact_seq(x)){
    return x;
  } else if (is_simple && Rf_length(x) < 10000){
    SEXP dup = SHIELD(Rf_duplicated(x, FALSE));
    SEXP unique_locs = SHIELD(cpp_which_(dup, true));
    SEXP out = SHIELD(sset_vec(x, unique_locs, false));
    Rf_copyMostAttrib(x, out);
    YIELD(3);
    return out;
  } else if (is_simple){
    return cheapr_fast_unique(x);
  } else {
    return cpp11::package("base")["unique"](x);
  }
}

[[cpp11::register]]
SEXP cpp_setdiff(SEXP x, SEXP y){
  SEXP matches = SHIELD(Rf_match(y, x, NA_INTEGER));
  SEXP locs = SHIELD(cpp_which_na(matches));
  SEXP out = SHIELD(sset_vec(x, locs, false));
  Rf_copyMostAttrib(x, out);
  YIELD(3);
  return out;
}

SEXP cpp_intersect(SEXP x, SEXP y, bool unique){
  if (unique){
    SHIELD(x = cpp_unique(x));
  } else {
    SHIELD(x);
  }
  SEXP matches = SHIELD(Rf_match(y, x, NA_INTEGER));
  SEXP locs = SHIELD(cpp_which_not_na(matches));
  SEXP out = SHIELD(sset_vec(x, locs, false));
  Rf_copyMostAttrib(x, out);
  YIELD(4);
  return out;
}
// SEXP cpp_setdiff(SEXP x, SEXP y, bool dups){
//
//   bool is_simple = is_simple_atomic_vec(x) && is_simple_atomic_vec(y);
//
//   int NP = 0;
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
  SEXP CHEAPR_ZERO = SHIELD(Rf_ScalarInteger(0));
  SEXP out = SHIELD(cpp_sset(x, CHEAPR_ZERO, true));
  YIELD(2);
  return out;
}

// Prototypes of data frame

[[cpp11::register]]
SEXP get_ptypes(SEXP x){
  int n = Rf_length(x);
  SEXP out = SHIELD(new_vec(VECSXP, n));

  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, get_ptype(VECTOR_ELT(x, i)));
  }

  Rf_namesgets(out, Rf_getAttrib(x, R_NamesSymbol));

  YIELD(1);
  return out;
}

// SEXP cpp_cast(SEXP x, SEXP y){
  // int NP = 0;
  //
  // SEXP x_cls = Rf_getAttrib(x, R_ClassSymbol);
  // SEXP y_cls = Rf_getAttrib(y, R_ClassSymbol);
  //
  // SHIELD(x_cls = coerce_vec(x_cls, STRSXP)); ++NP;
  // SHIELD(y_cls = coerce_vec(y_cls, STRSXP)); ++NP;
  //
  // if ( (TYPEOF(x) == TYPEOF(y)) && (Rf_length(x_cls) == Rf_length(y_cls))){
  //
  //   bool class_identical = true;
  //
  //   for (int i = 0; i < Rf_length(x_cls); ++i){
  //     if (std::strcmp(CHAR(STRING_ELT(x_cls, i)), CHAR(STRING_ELT(y_cls, i))) != 0){
  //       class_identical = false; break;
  //     }
  //   }
  //   if (class_identical){
  //     YIELD(NP);
  //     return x;
  //   }
  // }
  //
  // if (is_simple_vec(x) && is_simple_vec(y) &&
  //     !Rf_inherits(x, "factor") && !Rf_inherits(y, "factor")){
  //     if (CHEAPR_TYPEOF(x) >= CHEAPR_TYPEOF(y)){
  //       SEXP out = SHIELD(new_vec(TYPEOF(x), Rf_xlength(x))); ++NP;
  //       Rf_copyMostAttrib(y, out);
  //       YIELD(NP);
  //       return out;
  //     } else {
  //       SEXP out = SHIELD(coerce_vector(x, CHEAPR_TYPEOF(y))); ++NP;
  //       Rf_copyMostAttrib(y, out);
  //       YIELD(NP);
  //       return out;
  //     }
  // } else {
  //   YIELD(NP);
  //   Rf_error("Can't convert `x` based on `y` in %s", __func__);
  // }
// }

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
      R_Reprotect(fct_levels = Rf_getAttrib(p_x[i], R_LevelsSymbol), fct_idx);
    } else {
      R_Reprotect(fct_levels = base_as_character(p_x[i]), fct_idx);
    }
    SET_VECTOR_ELT(levels, i, fct_levels);
  }
  SEXP out = SHIELD(cpp_c(levels));
  SHIELD(out = cpp_unique(out));
  YIELD(4);
  return out;
}

[[cpp11::register]]
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
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

  R_xlen_t out_size = 0;
  for (R_xlen_t i = 0; i < n; ++i){
    out_size += (TYPEOF(p_x[i]) == VECSXP  ? Rf_xlength(p_x[i]) : 1);
  }

  SEXP x_names = SHIELD(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
  bool x_has_names = !is_null(x_names);

  R_xlen_t k = 0;
  SEXP out;

  out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
  SEXP container_list = SHIELD(new_vec(VECSXP, 1)); ++NP;
  Rf_setAttrib(container_list, R_NamesSymbol, R_BlankScalarString);

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
      R_Reprotect(names = Rf_getAttrib(p_x[i], R_NamesSymbol), nm_idx);
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
    Rf_namesgets(out, out_names);
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

  int NP = 0;
  const SEXP *p_x = VECTOR_PTR_RO(x);

  SEXP df, names;
  PROTECT_INDEX df_idx, names_idx;

  R_ProtectWithIndex(df = p_x[0], &df_idx); ++NP;
  R_ProtectWithIndex(names = Rf_getAttrib(df, R_NamesSymbol), &names_idx); ++NP;

  if (!is_df(df)){
    YIELD(NP); Rf_error("Can't combine data frames with non data frames");
  }

  SEXP df_template = df;

  SEXP frames = SHIELD(new_vec(VECSXP, n_frames)); ++NP;
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
  R_ProtectWithIndex(temp_list = new_vec(VECSXP, 2), &temp_list_idx); ++NP;

  // We do 2 passes
  // 1st pass: Check inputs are dfs and construct prototype list
  // 2nd pass: initialise data frame vecs (if need be) and combine simultaneously

  int out_size = df_nrow(df);

  for (int i = 1; i < n_frames; ++i){
    R_Reprotect(df = p_x[i], df_idx);

    if (!is_df(df)){
      YIELD(NP); Rf_error("Can't combine data frames with non data frames");
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
      R_Reprotect(ptypes = cpp_list_c(temp_list), ptypes_idx);
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

  SEXP out = SHIELD(new_vec(VECSXP, n_cols)); ++NP;
  SEXP vectors = SHIELD(new_vec(VECSXP, n_frames)); ++NP;

  for (int j = 0; j < n_cols; ++j){
    for (int i = 0; i < n_frames; ++i){
      R_Reprotect(vec = get_list_element(p_x[i], CHAR(STRING_ELT(names, j))), vec_idx);

      if (is_null(vec)){
        R_Reprotect(vec = VECTOR_ELT(ptypes, j), vec_idx);
        R_Reprotect(vec = cpp_na_init(vec, df_nrow(p_x[i])), vec_idx);
      }
      SET_VECTOR_ELT(vectors, i, vec);
    }
    SET_VECTOR_ELT(out, j, cpp_c(vectors));
  }
  SHIELD(out = cpp_list_as_df(out)); ++NP;
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(out_size));
  Rf_namesgets(out, names);
  SHIELD(out = fast_df_reconstruct(out, df_template)); ++NP;
  YIELD(NP);
  return out;
}

// Helper to concatenate data frames by col

[[cpp11::register]]
SEXP cpp_df_col_c(SEXP x, bool recycle, bool name_repair){
  int NP = 0;

  // Important to recycle first to avoid incorrect size calculations

  if (recycle){
    SHIELD(x = cpp_recycle(x, R_NilValue)); ++NP;
  }

  SEXP out = SHIELD(cpp_list_c(x)); ++NP;

  SEXP r_nrows = SHIELD(R_NilValue); ++NP;
  if (Rf_length(out) == 0 && Rf_length(x) != 0){
    SHIELD(r_nrows = Rf_ScalarInteger(vec_length(VECTOR_ELT(x, 0)))); ++NP;
  }

  SHIELD(out = cpp_new_df(out, r_nrows, false, name_repair)); ++NP;

  if (Rf_length(x) != 0 && is_df(VECTOR_ELT(x, 0))){
    SHIELD(out = fast_df_reconstruct(out, VECTOR_ELT(x, 0))); ++NP;
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
  int NP = 0;
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
    is_frame = is_frame || is_df(p_x[i]);
    is_classed = is_classed || Rf_isObject(p_x[i]);
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
    SEXP c_char = SHIELD(Rf_mkString("c"));
    SEXP out = SHIELD(base_do_call(c_char, x));
    YIELD(2);
    return out;
  }

  R_xlen_t k = 0;
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
    int *p_out = INTEGER(out);

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        coerce_vec(p_x[i], vector_type), temp_idx
      );
      int *p_temp = INTEGER(temp);
      memcpy(&p_out[k], &p_temp[0], m * sizeof(int));
      k += m;
    }
    break;
  }
  case REALSXP: {
    out = SHIELD(new_vec(vector_type, out_size)); ++NP;
    double *p_out = REAL(out);

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        coerce_vec(p_x[i], vector_type), temp_idx
      );
      double *p_temp = REAL(temp);
      memcpy(&p_out[k], &p_temp[0], m * sizeof(double));
      k += m;
    }
    break;
  }
  case STRSXP: {
    out = SHIELD(new_vec(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        coerce_vec(p_x[i], vector_type), temp_idx
      );
      const SEXP *p_temp = STRING_PTR_RO(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_STRING_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }
  case CPLXSXP: {
    out = SHIELD(new_vec(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        coerce_vec(p_x[i], vector_type), temp_idx
      );
      Rcomplex *p_temp = COMPLEX(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_COMPLEX_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }
  case VECSXP: {
    out = SHIELD(new_vec(vector_type, out_size)); ++NP;

    for (int i = 0; i < n; ++i){
      R_xlen_t m = Rf_xlength(p_x[i]);
      R_Reprotect(
        temp = TYPEOF(p_x[i]) == vector_type ?
      p_x[i] :
        coerce_vec(p_x[i], vector_type), temp_idx
      );
      const SEXP *p_temp = VECTOR_PTR_RO(temp);
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
      }
    }
    break;
  }

  default: {
    SEXP c_char = SHIELD(Rf_mkString("c")); ++NP;
    out = SHIELD(base_do_call(c_char, x)); ++NP;
    break;
  }
  }
  if (is_date2){
    Rf_classgets(out, Rf_mkChar("Date"));
  }
  if (is_datetime2){
    Rf_copyMostAttrib(p_x[first_datetime], out);
  }
  YIELD(NP);
  return out;
}
