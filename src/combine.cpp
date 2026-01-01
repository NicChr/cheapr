#include "cheapr.h"
#include "R.h"

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
        SEXP out = SHIELD(shallow_copy ? cheapr::vec::shallow_copy(x) : x);
        attr::set_old_class(out, get_old_class(source));
        YIELD(1);
        return out;
      }
    } else if (is_bare_tbl(source)){

      // Only copy the class if source is a plain tbl

      if (!shallow_copy && is_bare_tbl(x)){
        return x;
      } else {
        SEXP out = SHIELD(shallow_copy ? cheapr::vec::shallow_copy(x) : x);
        attr::set_old_class(out, get_old_class(source));
        YIELD(1);
        return out;
      }
    } else {

      // Method dispatch, users can write `rebuild` methods which
      // this will use
      return eval_pkg_fun("rebuild", "cheapr", env::base_env, x, source, arg("shallow_copy") = shallow_copy);
    }
  } else {
    return eval_pkg_fun("rebuild", "cheapr", env::base_env, x, source, arg("shallow_copy") = shallow_copy);
  }
}
[[cpp11::register]]
SEXP cpp_rep_len(SEXP x, R_xlen_t length){
  R_xlen_t out_size = length;

  if (is_null(x)){
    return r_null;
  } else if (is_df(x)){
    if (out_size == df::nrow(x)) return x;
    int n_cols = Rf_length(x);
    SEXP out = SHIELD(new_list(n_cols));
    const SEXP *p_x = list_ptr_ro(x);
    for (int i = 0; i < n_cols; ++i){
      SET_VECTOR_ELT(out, i, cpp_rep_len(p_x[i], length));
    }
    SEXP names = SHIELD(get_old_names(x));
    set_old_names(out, names);
    set_list_as_df(out);
    df::set_row_names(out, length);
    SHIELD(out = rebuild(out, x, false));
    YIELD(3);
    return out;
  } else if (cheapr_is_simple_vec2(x)){

    R_xlen_t size = Rf_xlength(x);

    // Return x if length(x) == length
    if (out_size == size) return x;

    SEXP out = r_null;

    switch (CHEAPR_TYPEOF(x)){

    case NILSXP: {
      break;
    }
    case LGLSXP: {

      out = SHIELD(internal::new_vec(LGLSXP, out_size));
      auto *p_x = internal::logical_ptr_ro(x);
      auto *p_out = internal::logical_ptr(out);

      if (size == 1){
        r_fill(out, p_out, 0, out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        // Copy first block
        std::copy_n(p_x, size, p_out);

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          std::copy_n(p_out, to_copy, p_out + copied);
          copied += to_copy;
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        r_fill(out, p_out, 0, out_size, na::logical);
      }
      break;
    }
    case INTSXP: {
      out = SHIELD(internal::new_vec(INTSXP, out_size));
      auto *p_x = internal::integer_ptr_ro(x);
      auto *p_out = internal::integer_ptr(out);

      if (size == 1){
        r_fill(out, p_out, 0, out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        // Copy first block
        std::copy_n(p_x, size, p_out);

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          std::copy_n(p_out, to_copy, p_out + copied);
          copied += to_copy;
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        r_fill(out, p_out, 0, out_size, na::integer);
      }
      break;
    }
    case REALSXP: {
      out = SHIELD(internal::new_vec(REALSXP, out_size));
      auto *p_x = internal::real_ptr_ro(x);
      auto *p_out = internal::real_ptr(out);

      if (size == 1){
        r_fill(out, p_out, 0, out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        // Copy first block
        std::copy_n(p_x, size, p_out);

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          std::copy_n(p_out, to_copy, p_out + copied);
          copied += to_copy;
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        r_fill(out, p_out, 0, out_size, na::real);
      }
      break;
    }
    case CHEAPR_INT64SXP: {
      out = SHIELD(new_vector<int64_t>(out_size));
      auto *p_x = internal::integer64_ptr_ro(x);
      auto *p_out = internal::integer64_ptr(out);

      if (size == 1){
        r_fill(out, p_out, 0, out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        // Copy first block
        std::copy_n(p_x, size, p_out);

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          std::copy_n(p_out, to_copy, p_out + copied);
          copied += to_copy;
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        r_fill(out, p_out, 0, out_size, na::integer64);
      }
      break;
    }
    case STRSXP: {
      out = SHIELD(internal::new_vec(STRSXP, out_size));
      auto *p_x = internal::string_ptr_ro(x);
      auto *p_out = internal::string_ptr_ro(out);

      if (size == 1){
        for (R_xlen_t i = 0; i < out_size; ++i){
          SET_STRING_ELT(out, i, p_x[0]);
        }
      } else if (out_size > 0 && size > 0){

        for (R_xlen_t i = 0; i < size; ++i){
          SET_STRING_ELT(out, i, p_x[i]);
        }

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          for (R_xlen_t i = 0; i < to_copy; ++i){
            SET_STRING_ELT(out, copied + i, p_out[i]);
          }
          copied += to_copy;
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        for (R_xlen_t i = 0; i < out_size; ++i){
          SET_STRING_ELT(out, i, na::string);
        }
      }
      break;
    }
    case CPLXSXP: {
      out = SHIELD(internal::new_vec(CPLXSXP, out_size));
      auto *p_x = internal::complex_ptr_ro(x);
      auto *p_out = internal::complex_ptr(out);

      if (size == 1){
        r_fill(out, p_out, 0, out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        // Copy first block
        std::copy_n(p_x, size, p_out);

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          std::copy_n(p_out, to_copy, p_out + copied);
          copied += to_copy;
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        r_fill(out, p_out, 0, out_size, na::complex);
      }
      break;
    }
    case RAWSXP: {
      out = SHIELD(internal::new_vec(RAWSXP, out_size));
      auto *p_x = internal::raw_ptr_ro(x);
      auto *p_out = internal::raw_ptr(out);

      if (size == 1){
        r_fill(out, p_out, 0, out_size, p_x[0]);
      } else if (out_size > 0 && size > 0){
        // Copy first block
        std::copy_n(p_x, size, p_out);

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          std::copy_n(p_out, to_copy, p_out + copied);
          copied += to_copy;
        }
        // If length > 0 but length(x) == 0 then fill with NA
      } else if (size == 0 && out_size > 0){
        r_fill(out, p_out, 0, out_size, na::raw);
      }
      break;
    }
    case VECSXP: {
      out = SHIELD(internal::new_vec(VECSXP, out_size));
      auto *p_x = VECTOR_PTR_RO(x);
      auto *p_out = VECTOR_PTR_RO(out);

      if (size == 1){
        for (R_xlen_t i = 0; i < out_size; ++i){
          SET_VECTOR_ELT(out, i, p_x[0]);
        }
      } else if (out_size > 0 && size > 0){

        for (R_xlen_t i = 0; i < size; ++i){
          SET_VECTOR_ELT(out, i, p_x[i]);
        }

        // copy result to itself, doubling each iteration
        R_xlen_t copied = size;
        while (copied < out_size) {
          R_xlen_t to_copy = std::min(copied, out_size - copied);
          for (R_xlen_t i = 0; i < to_copy; ++i){
            SET_VECTOR_ELT(out, copied + i, p_out[i]);
          }
          copied += to_copy;
        }
      }
      break;
    }
    default: {
      if (r_length(x) == length){
      return x;
    } else {
      return eval_pkg_fun("rep", "base", env::base_env, x, arg("length.out") = length);
    }
    }
    }
    Rf_copyMostAttrib(x, out);
    YIELD(1);
    return out;
  } else {
    if (r_length(x) == length){
      return x;
    } else {
      return eval_pkg_fun("rep", "base", env::base_env, x, arg("length.out") = length);
    }
  }
}

[[cpp11::register]]
SEXP cpp_rep(SEXP x, SEXP times){

  R_xlen_t n = vec::length(x);

  R_xlen_t out_size;
  R_xlen_t n_times = Rf_xlength(times);
  SHIELD(times = cast<r_integers_t>(times, r_null));

  if (n_times == 1){
    out_size = n * integer_ptr(times)[0];
    SEXP out = SHIELD(cpp_rep_len(x, out_size));
    YIELD(2);
    return out;
  } else {
    int *p_times = integer_ptr(times);
    if (n_times != n){
      YIELD(1);
      Rf_error("`times` must be length 1 or `vec::length(x)` in %s", __func__);
    }

    if (is_null(x)){
      YIELD(1);
      return r_null;
    } else if (is_df(x)){
      if (Rf_length(x) == 0){
        SEXP out = SHIELD(vec::shallow_copy(x));
        df::set_row_names(out, cpp_sum(times));
        SHIELD(out = rebuild(out, x, false));
        YIELD(3);
        return out;
      } else {
        SEXP row_names = SHIELD(cpp_seq_len(df::nrow(x)));
        SHIELD(row_names = cpp_rep(row_names, times));
        SEXP out = SHIELD(cpp_sset(x, row_names, true));
        YIELD(4);
        return out;
      }

    } else if (cheapr_is_simple_vec2(x)){
      SEXP out = SHIELD(internal::new_vec(TYPEOF(x), cpp_sum(times)));
      switch (TYPEOF(x)){
      case LGLSXP:
      case INTSXP: {
        int *p_x = integer_ptr(x);
        int *p_out = integer_ptr(out);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          std::fill_n(p_out + k, p_times[i], p_x[i]);
          k += p_times[i];
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case REALSXP: {
        double *p_x = real_ptr(x);
        double *p_out = real_ptr(out);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          std::fill_n(p_out + k, p_times[i], p_x[i]);
          k += p_times[i];
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case STRSXP: {
        const r_string_t *p_x = string_ptr_ro(x);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          for (int j = 0; j < p_times[i]; ++j, ++k) set_value<r_string_t>(out, k, p_x[i]);
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case CPLXSXP: {
        const r_complex_t *p_x = complex_ptr(x);
        r_complex_t *p_out = complex_ptr(x);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          std::fill_n(p_out + k, p_times[i], p_x[i]);
          k += p_times[i];
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      case VECSXP: {
        const SEXP *p_x = list_ptr_ro(x);
        R_xlen_t k = 0;
        for (R_xlen_t i = 0; i < n; ++i){
          for (int j = 0; j < p_times[i]; ++j, ++k) SET_VECTOR_ELT(out, k, p_x[i]);
        }
        Rf_copyMostAttrib(x, out);
        YIELD(2);
        return out;
      }
      default: {
        SEXP out = SHIELD(eval_pkg_fun("rep", "base", env::base_env, x, times));
        YIELD(3);
        return out;
      }
      }
    } else {
      SEXP out = SHIELD(eval_pkg_fun("rep", "base", env::base_env, x, times));
      YIELD(2);
      return out;
    }
  }
}

[[cpp11::register]]
SEXP cpp_rep_each(SEXP x, SEXP each){
  int32_t NP = 0;
  SHIELD(each = cast<r_integers_t>(each, r_null)); ++NP;
  if (Rf_length(each) == 1){
    if (integer_ptr(each)[0] == 1){
      YIELD(NP);
      return x;
    }
    SHIELD(each = cpp_rep_len(each, vec::length(x))); ++NP;
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

  const SEXP *p_x = list_ptr_ro(x);
  R_xlen_t out = 0;

  for (R_xlen_t i = 0; i < n; ++i){
    if (is_null(p_x[i])) continue;

    if (vec::length(p_x[i]) == 0){
      out = 0;
      break;
    }
    out = std::max(out, vec::length(p_x[i]));
  }
  return out;
}

// Recycle elements of a list `x`

void recycle_in_place(SEXP x, R_xlen_t n){

  int xn = Rf_length(x);
  const SEXP *p_x = list_ptr_ro(x);

  for (int i = 0; i < xn; ++i){
    if (!is_null(p_x[i])){
      SET_VECTOR_ELT(x, i, cpp_rep_len(p_x[i], n));
    }
  }
}

[[cpp11::register]]
SEXP cpp_recycle(SEXP x, SEXP length){

  SEXP out = SHIELD(cpp_drop_null(x));

  int n = 0;
  int32_t NP = 1;

  if (!is_null(length)){
    if (TYPEOF(length) == INTSXP){
      n = integer_ptr(length)[0];
    } else {
      SHIELD(length = cast<r_doubles_t>(length, r_null));
      NP = 2;
      n = real_ptr(length)[0];
    }
    if (n < 0){
      YIELD(NP);
      Rf_error("Recycled `length` must be >= 0 in %s", __func__);
    }
  } else {
    n = length_common(out);
  }

  recycle_in_place(out, n);

  YIELD(NP);
  return out;
}

// Fast unique that can be used in C code
// Doesn't return unique df rows
SEXP cpp_unique(SEXP x, bool names){

  int32_t NP = 0;

  bool is_simple = cheapr_is_simple_atomic_vec(x);

  if (is_compact_seq(x)){
    YIELD(NP);
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
        SEXP names = SHIELD(get_old_names(x)); ++NP;
        SHIELD(names = sset_vec(names, unique_locs, false)); ++NP;
        set_old_names(out, names);
      }
      YIELD(NP);
      return out;
    }
  } else if (is_simple){
    SEXP out = SHIELD(eval_pkg_fun("unique_", "cheapr", env::base_env, x)); ++NP;
    YIELD(NP);
    return out;
  } else {
    SEXP out = SHIELD(eval_pkg_fun("unique", "base", env::base_env, x)); ++NP;
    YIELD(NP);
    return out;
  }
}

[[cpp11::register]]
SEXP cpp_setdiff(SEXP x, SEXP y, bool unique){
  if (unique){
    SHIELD(x = cpp_unique(x, true));
  } else {
    SHIELD(x);
  }
  SEXP matches = SHIELD(match(y, x, na::integer));
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
  SEXP matches = SHIELD(match(y, x, na::integer));
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
  SEXP out = SHIELD(new_list(n));

  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, get_ptype(VECTOR_ELT(x, i)));
  }

  SEXP names = SHIELD(get_old_names(x));
  set_old_names(out, names);

  YIELD(2);
  return out;
}

// Helper to turn a factor into a character vec
SEXP factor_as_character(SEXP x){
  return sset_vec(get_attr(x, symbol::levels_sym), x, true);
}

// Helper to turn a character vec into a factor, given levels
SEXP character_as_factor(SEXP x, SEXP levels){

  if (TYPEOF(x) != STRSXP){
    Rf_error("`x` must be a character vector in %s", __func__);
  }

  int32_t NP = 0;

  SEXP out = SHIELD(match(levels, x, na::integer)); ++NP;
  SEXP cls = SHIELD(as_vector("factor")); ++NP;
  set_attr(out, symbol::levels_sym, levels);
  attr::set_old_class(out, cls);
  YIELD(NP);
  return out;
}

// Helper to concatenate plain lists
// It uses top-level names if and only if x[[i]] is not a list

[[cpp11::register]]
SEXP cpp_list_c(SEXP x){
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = list_ptr_ro(x);

  R_xlen_t out_size = 0;
  for (R_xlen_t i = 0; i < n; ++i){
    out_size += (TYPEOF(p_x[i]) == VECSXP ? Rf_xlength(p_x[i]) : 1);
  }

  SEXP x_names = SHIELD(get_old_names(x)); ++NP;
  bool x_has_names = !is_null(x_names);

  R_xlen_t k = 0;
  SEXP out = SHIELD(new_list(out_size)); ++NP;
  SEXP container_list = SHIELD(make_list(arg("") = r_null)); ++NP;

  SEXP names;
  PROTECT_INDEX nm_idx;
  R_ProtectWithIndex(names = r_null, &nm_idx); ++NP;

  SEXP out_names = SHIELD(new_vector<r_string_t>(out_size)); ++NP;
  bool any_names = false;

  R_xlen_t m;

  for (R_xlen_t i = 0; i < n; ++i){
    const SEXP *p_temp;

    if (TYPEOF(p_x[i]) == VECSXP){
      p_temp = list_ptr_ro(p_x[i]);
      R_Reprotect(names = get_old_names(p_x[i]), nm_idx);
      m = Rf_xlength(p_x[i]);
    } else {
      SET_VECTOR_ELT(container_list, 0, p_x[i]);
      if (x_has_names){
        R_Reprotect(names = as_vector(get_value<r_string_t>(x_names, i)), nm_idx);
      } else {
        names = r_null;
      }
      p_temp = list_ptr_ro(container_list);
      m = 1;
    }

    any_names = any_names || !is_null(names);
    if (!is_null(names)){
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
        set_value<r_string_t>(out_names, k, get_value<r_string_t>(names, j));
      }
    } else {
      for (R_xlen_t j = 0; j < m; ++k, ++j){
        SET_VECTOR_ELT(out, k, p_temp[j]);
      }
    }
  }
  if (any_names){
    set_old_names(out, out_names);
  }
  YIELD(NP);
  return out;
}

SEXP get_list_element(SEXP list, r_string_t col){
  const r_string_t *names = vector_ptr<const r_string_t>(get_old_names(list));
  for (int i = 0; i < Rf_length(list); ++i){
    if (names[i] == col){
      return VECTOR_ELT(list, i);
    }
  }
  return r_null;
}

// Concatenate a list of data frames

SEXP cpp_df_c(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of data frames");
  }

  int n_frames = Rf_length(x);

  if (n_frames == 0) return r_null;

  int32_t NP = 0;

  // Since we are casting potential non-dfs to dfs we need a new list
  SEXP frames = SHIELD(init<r_list_t>(n_frames, false)); ++NP;

  const SEXP *p_x = list_ptr_ro(x);
  const SEXP *p_frames = list_ptr_ro(frames);

  SEXP df;
  PROTECT_INDEX df_idx;
  R_ProtectWithIndex(df = cast<r_data_frame_t>(p_x[0], r_null), &df_idx); ++NP;
  SET_VECTOR_ELT(frames, 0, df);

  SEXP new_names, df_names, ptype_names;
  PROTECT_INDEX new_names_idx, df_names_idx, ptype_names_idx;

  R_ProtectWithIndex(new_names = r_null, &new_names_idx); ++NP;
  R_ProtectWithIndex(df_names = r_null, &df_names_idx); ++NP;
  R_ProtectWithIndex(ptype_names = get_old_names(df), &ptype_names_idx); ++NP;

  int n_cols = Rf_length(ptype_names);

  // We do 2 passes
  // 1st pass: Cast inputs to dfs and construct prototypes
  // 2nd pass: initialise data frame vecs (if need be) and combine simultaneously

  int out_size = df::nrow(df);

  bool na_padding = false;

  for (int i = 1; i < n_frames; ++i){
    R_Reprotect(df = cast<r_data_frame_t>(p_x[i], r_null), df_idx);
    R_Reprotect(df_names = get_old_names(df), df_names_idx);

    R_Reprotect(new_names = cpp_setdiff(
      df_names, ptype_names, false
    ), new_names_idx);

    // Adjust prototype names
    if (Rf_length(new_names) > 0){
      na_padding = true;
      R_Reprotect(ptype_names = combine(ptype_names, new_names), ptype_names_idx);
    }
    na_padding = na_padding || Rf_length(df) != n_cols;
    out_size += df::nrow(df);
    SET_VECTOR_ELT(frames, i, df);
  }

  n_cols = Rf_length(ptype_names);

  SEXP vec;
  PROTECT_INDEX vec_idx;
  R_ProtectWithIndex(vec = r_null, &vec_idx); ++NP;

  SEXP out = SHIELD(new_list(n_cols)); ++NP;
  SEXP vectors = SHIELD(new_list(n_frames)); ++NP;

  const r_string_t *p_ptype_names = string_ptr_ro(ptype_names);

  if (na_padding){
    // Get archetype of each col
    SEXP vec_archetypes = SHIELD(new_list(n_cols)); ++NP;
    const SEXP *p_vec_archetypes = list_ptr_ro(vec_archetypes);
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
          R_Reprotect(vec = cpp_na_init(p_vec_archetypes[j], df::nrow(p_frames[i])), vec_idx);
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
  df::set_row_names(out, out_size);
  set_old_names(out, ptype_names);
  SHIELD(out = rebuild(out, VECTOR_ELT(frames, 0), false)); ++NP;
  YIELD(NP);
  return out;
}

// Helper to concatenate data frames by col

[[cpp11::register]]
SEXP cpp_df_col_c(SEXP x, bool recycle, bool name_repair){

  int32_t NP = 0;
  R_xlen_t common_size = length_common(x);
  SEXP out = SHIELD(cpp_list_c(x)); ++NP;
  SEXP df_nrows = SHIELD(as_vector(r_cast<int>(common_size))); ++NP;
  SHIELD(out = cpp_new_df(out, df_nrows, recycle, name_repair)); ++NP;

  if (Rf_length(x) != 0 && is_df(VECTOR_ELT(x, 0))){
    SHIELD(out = rebuild(out, VECTOR_ELT(x, 0), false)); ++NP;
  }

  YIELD(NP);
  return out;
}

// Combine vectors given the following args:
// * A list of vectors
// * The final size of the combined vector
// * A template vector (of length 0) whose type is the common type between all vectors
// This allows us to eliminate a few loops
SEXP combine_internal(SEXP x, const R_xlen_t out_size, SEXP vec_template){

  if (vec::length(vec_template) != 0){
    Rf_error("`vec_template` must be of length 0");
  }

  int n = Rf_length(x);
  const SEXP *p_x = list_ptr_ro(x);

  int32_t NP = 0;

  R_xlen_t k = 0;
  R_xlen_t m = 0;
  SEXP out = r_null;

  SEXP vec;
  PROTECT_INDEX vec_idx;
  R_ProtectWithIndex(vec = r_null, &vec_idx); ++NP;

  r_type common = get_r_type(vec_template);
  SEXP combined_names = r_null;
  bool has_top_level_names = vec_has_names(x);

  // If in the future it is decided that inner names are to be kept
  // Use this condition within the inner loop and then use the macro
  // defined above combine_internal
  // any_names = !has_top_level_names && (any_names || has_names(vec));

  switch (common){
  case R_null: {
    break;
  }
  case R_lgl: {

    SHIELD(out = init<r_logicals_t>(out_size, false)); ++NP;

    r_bool_t* RESTRICT p_out = logical_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_logicals_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(vector_ptr<const r_bool_t>(vec), m, p_out + k);
    }
    break;
  }
  case R_int: {

    SHIELD(out = init<r_integers_t>(out_size, false)); ++NP;

    int* RESTRICT p_out = integer_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_integers_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::integer_ptr_ro(vec), m, p_out + k);
    }
    break;
  }
  case R_int64: {

    SHIELD(out = init<r_integers64_t>(out_size, false)); ++NP;

    int64_t* RESTRICT p_out = integer64_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_integers64_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::integer64_ptr_ro(vec), m, p_out + k);
    }
    break;
  }
  case R_dbl: {

    SHIELD(out = init<r_doubles_t>(out_size, false)); ++NP;

    double* RESTRICT p_out = real_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_doubles_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::real_ptr_ro(vec), m, p_out + k);
    }
    break;
  }
  case R_chr: {

    SHIELD(out = init<r_characters_t>(out_size, false)); ++NP;

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_characters_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      for (R_xlen_t i = 0; i < m; ++i) {
        SET_STRING_ELT(out, k + i, STRING_ELT(vec, i));
      }
    }
    break;
  }
  case R_cplx: {

    SHIELD(out = init<r_complexes_t>(out_size, false)); ++NP;

    r_complex_t *p_out = complex_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_complexes_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::complex_ptr_ro(vec), m, p_out + k);
    }
    break;
  }
  case R_raw: {

    SHIELD(out = init<r_raws_t>(out_size, false)); ++NP;

    r_byte_t *p_out = raw_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_raws_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::raw_ptr_ro(vec), m, p_out + k);
    }
    break;
  }
  case R_list: {

    SHIELD(out = init<r_list_t>(out_size, false)); ++NP;

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_list_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      for (R_xlen_t i = 0; i < m; ++i) {
        SET_VECTOR_ELT(out, k + i, VECTOR_ELT(vec, i));
      }
    }
    break;
  }
  case R_fct: {

    SHIELD(out = cpp_rep_len(vec_template, out_size)); ++NP;

    int* RESTRICT p_out = integer_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_factors_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::integer_ptr_ro(vec), m, p_out + k);
    }
    break;
  }
  case R_date: {

    if (TYPEOF(vec_template) == INTSXP){
    SHIELD(out = init<r_integers_t>(out_size, false)); ++NP;
    SEXP date_cls = SHIELD(as_vector("Date")); ++NP;
    attr::set_old_class(out, date_cls);

    int* RESTRICT p_out = integer_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_dates_t>(p_x[i], vec_template), vec_idx);
      R_Reprotect(vec = vec::coerce_vec(vec, INTSXP), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::integer_ptr_ro(vec), m, p_out + k);
    }
  } else {
    SHIELD(out = init<r_doubles_t>(out_size, false)); ++NP;
    SEXP date_cls = SHIELD(as_vector("Date")); ++NP;
    attr::set_old_class(out, date_cls);

    double* RESTRICT p_out = real_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_dates_t>(p_x[i], vec_template), vec_idx);
      R_Reprotect(vec = vec::coerce_vec(vec, REALSXP), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::real_ptr_ro(vec), m, p_out + k);
    }
  }
  break;
  }
  case R_pxt: {

    SHIELD(out = init<r_posixts_t>(out_size, false)); ++NP;
    SEXP out_tzone = SHIELD(get_attr(vec_template, r_cast<r_symbol_t>("tzone"))); ++NP;
    set_attr(out, r_cast<r_symbol_t>("tzone"), out_tzone);

    double* RESTRICT p_out = real_ptr(out);

    for (int i = 0; i < n; ++i, k += m){
      R_Reprotect(vec = cast<r_posixts_t>(p_x[i], vec_template), vec_idx);
      m = Rf_xlength(vec);
      std::copy_n(internal::real_ptr_ro(vec), m, p_out + k);
    }
    break;
  }
  case R_df: {
    SHIELD(out = cpp_df_c(x)); ++NP;
    break;
  }
  case R_unk: {
    SEXP call = SHIELD(vec::coerce_vec(x, LISTSXP)); ++NP;
    SHIELD(call = Rf_lcons(r_cast<r_symbol_t>("c"), call)); ++NP;
    SHIELD(out = eval(call, env::base_env)); ++NP;
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
  set_vec_names(out, combined_names);
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
  const SEXP *p_x = list_ptr_ro(x);

  // Figure out final size
  R_xlen_t out_size = 0;
  for (int i = 0; i < n; ++i) out_size += vec::length(p_x[i]);

  // 'vec_template' here acts as a template for the final result
  SEXP vec_template = SHIELD(cpp_common_template(x));
  SEXP out = SHIELD(combine_internal(x, out_size, vec_template));
  YIELD(2);
  return out;
}
