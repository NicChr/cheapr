#include "cheapr.h"
#include <R.h>

// Coalesce a list of string vectors

[[cpp11::register]]
SEXP cpp_str_coalesce(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of character vectors in %s", __func__);
  }

  int32_t NP = 0;

  SEXP chars = SHIELD(cpp_recycle(x, r_null)); ++NP;
  if (MAYBE_REFERENCED(chars)){
    SHIELD(chars = vec::shallow_copy(chars)); ++NP;
  }

  R_xlen_t n_chars = Rf_xlength(chars);

  if (n_chars == 0){
    SEXP out = SHIELD(new_vector<r_string_t>(0)); ++NP;
    YIELD(NP);
    return out;
  }

  const SEXP *p_chars = list_ptr_ro(chars);

  // First cast all to character vectors

  for (R_xlen_t i = 0; i < n_chars; ++i){
    SET_VECTOR_ELT(chars, i, cast<r_characters_t>(p_chars[i], r_null));
  }

  // The reason for not assigning these ptrs in the previous loop is because
  // cast<> may fail and the R_Calloc'd memory will leak if this happens
  // At this point onwards the code is assumed to be safe

  const r_string_t **char_ptrs = (const r_string_t **) R_Calloc(n_chars, const r_string_t *);

  for (R_xlen_t i = 0; i < n_chars; ++i){
    char_ptrs[i] = string_ptr_ro(p_chars[i]);
  }

  R_xlen_t n_strings = Rf_xlength(VECTOR_ELT(chars, 0));

  if (n_strings == 0){
    SEXP out = SHIELD(new_vector<r_string_t>(0)); ++NP;
    R_Free(char_ptrs);
    YIELD(NP);
    return out;
  }

  SEXP out = SHIELD(new_vector<r_string_t>(n_strings)); ++NP;
  r_string_t inner_char = blank_r_string;

  R_xlen_t n_nas;

  for (R_xlen_t i = 0; i < n_strings; ++i){
    n_nas = 0;
    for (R_xlen_t j = 0; j < n_chars; ++j){
      inner_char = char_ptrs[j][i];
      n_nas += is_r_na(inner_char);
      if (!(inner_char == blank_r_string || is_r_na(inner_char))){
        set_value(out, i, inner_char);
        break;
      }
      // If all ith elements are NA, then return NA
      if (n_nas == n_chars){
        set_value(out, i, na::string);
      }
    }
  }
  R_Free(char_ptrs);
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_paste(SEXP x, SEXP sep, SEXP collapse){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of character vectors in %s", __func__);
  }

  int32_t NP = 0;

  SHIELD(sep = cast<r_characters_t>(sep, r_null)); ++NP;

  if (Rf_length(sep) != 1){
    YIELD(NP);
    Rf_error("`sep` must be a `<character>` vector of length 1 in %s", __func__);
  }

  SEXP chars = SHIELD(cpp_recycle(x, r_null)); ++NP;
  if (MAYBE_REFERENCED(chars)){
    SHIELD(chars = vec::shallow_copy(chars)); ++NP;
  }

  R_xlen_t n_chars = Rf_xlength(chars);

  if (n_chars == 0){
    SEXP out = SHIELD(new_vector<r_string_t>(0)); ++NP;
    YIELD(NP);
    return out;
  }

  const SEXP *p_chars = list_ptr_ro(chars);

  // First cast all to character vectors

  for (R_xlen_t i = 0; i < n_chars; ++i){
    SET_VECTOR_ELT(chars, i, cast<r_characters_t>(p_chars[i], r_null));
  }

  if (!is_null(collapse)){
    SHIELD(collapse = cast<r_characters_t>(collapse, r_null)); ++NP;
    if (Rf_length(collapse) != 1){
      YIELD(NP);
      Rf_error("`.collapse` must be a length 1 character vector in %s", __func__);
    }
  }

  const r_string_t **char_ptrs = (const r_string_t **) R_Calloc(n_chars, const r_string_t *);

  for (R_xlen_t i = 0; i < n_chars; ++i){
    char_ptrs[i] = string_ptr_ro(p_chars[i]);
  }

  R_xlen_t n_strings = Rf_xlength(VECTOR_ELT(chars, 0));

  if (n_strings == 0){
    SEXP out = SHIELD(new_vector<r_string_t>(0)); ++NP;
    R_Free(char_ptrs);
    YIELD(NP);
    return out;
  }

  std::string sep1 = utf8_char(get_value<r_string_t>(sep, 0));
  std::string sep2 = utf8_char(blank_r_string);
  std::string strng = utf8_char(blank_r_string);

  if (!is_null(collapse)){

    sep2 = utf8_char(get_value<r_string_t>(collapse, 0));

    for (R_xlen_t j = 0; j < n_strings; ++j){
      if (j != 0) strng += sep2;
      for (R_xlen_t i = 0; i < n_chars; ++i){
        strng += utf8_char(char_ptrs[i][j]);
        if (i != (n_chars - 1)) strng += sep1;
      }
    }

    SEXP out = SHIELD(as_vector(strng.c_str())); ++NP;
    R_Free(char_ptrs);
    YIELD(NP);
    return out;

  } else {

    SEXP out = SHIELD(new_vector<r_string_t>(n_strings)); ++NP;

    // Along the length of each char vec

    for (R_xlen_t j = 0; j < n_strings; ++j){

      strng = utf8_char(char_ptrs[0][j]);

      for (R_xlen_t i = 1; i < n_chars; ++i){
        str_paste(strng, sep1, utf8_char(char_ptrs[i][j]));
      }
      set_value(out, j, r_cast<r_string_t>(strng.c_str()));
    }

    R_Free(char_ptrs);
    YIELD(NP);
    return out;
  }
}
