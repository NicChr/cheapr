#include "cheapr.h"

// Coalesce a list of string vectors

[[cpp11::register]]
SEXP cpp_str_coalesce(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of character vectors in %s", __func__);
  }

  int32_t NP = 0;

  SEXP chars = SHIELD(cpp_recycle(x, R_NilValue)); ++NP;
  if (MAYBE_REFERENCED(chars)){
    SHIELD(chars = Rf_shallow_duplicate(chars)); ++NP;
  }

  R_xlen_t n_chars = Rf_xlength(chars);

  if (n_chars == 0){
    SEXP out = SHIELD(new_vec(STRSXP, 0)); ++NP;
    YIELD(NP);
    return out;
  }

  const SEXP *p_chars = LIST_PTR_RO(chars);
  std::vector<const SEXP *> char_ptrs(n_chars);

  // First cast all to character vectors

  for (R_xlen_t i = 0; i < n_chars; ++i){
    SET_VECTOR_ELT(chars, i, cast<r_character_t>(p_chars[i], R_NilValue));
    char_ptrs[i] = STRING_PTR_RO(p_chars[i]);
  }

  R_xlen_t n_strings = Rf_xlength(VECTOR_ELT(chars, 0));

  if (n_strings == 0){
    SEXP out = SHIELD(new_vec(STRSXP, 0)); ++NP;
    YIELD(NP);
    return out;
  }

  SEXP out = SHIELD(Rf_allocVector(STRSXP, n_strings)); ++NP;
  SEXP inner_char = R_BlankString;

  R_xlen_t n_nas;

  for (R_xlen_t i = 0; i < n_strings; ++i){
    n_nas = 0;
    for (R_xlen_t j = 0; j < n_chars; ++j){
      inner_char = char_ptrs[j][i];
      n_nas += is_r_na(inner_char);
      if (!(inner_char == R_BlankString || is_r_na(inner_char))){
        SET_STRING_ELT(out, i, inner_char);
        break;
      }
      // If all ith elements are NA, then return NA
      if (n_nas == n_chars){
        SET_STRING_ELT(out, i, NA_STRING);
      }
    }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_paste(SEXP x, SEXP sep, SEXP collapse){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of character vectors in %s", __func__);
  }

  int32_t NP = 0;

  SHIELD(sep = cast<r_character_t>(sep, R_NilValue)); ++NP;

  if (Rf_length(sep) != 1){
    YIELD(NP);
    Rf_error("`sep` must be a `<character>` vector of length 1 in %s", __func__);
  }

  SEXP chars = SHIELD(cpp_recycle(x, R_NilValue)); ++NP;
  if (MAYBE_REFERENCED(chars)){
    SHIELD(chars = Rf_shallow_duplicate(chars)); ++NP;
  }

  R_xlen_t n_chars = Rf_xlength(chars);

  if (n_chars == 0){
    SEXP out = SHIELD(new_vec(STRSXP, 0)); ++NP;
    YIELD(NP);
    return out;
  }

  const SEXP *p_chars = LIST_PTR_RO(chars);
  std::vector<const SEXP *> char_ptrs(n_chars);

  // First cast all to character vectors

  for (R_xlen_t i = 0; i < n_chars; ++i){
    SET_VECTOR_ELT(chars, i, cast<r_character_t>(p_chars[i], R_NilValue));
    char_ptrs[i] = STRING_PTR_RO(p_chars[i]);
  }

  R_xlen_t n_strings = Rf_xlength(VECTOR_ELT(chars, 0));

  if (n_strings == 0){
    SEXP out = SHIELD(new_vec(STRSXP, 0)); ++NP;
    YIELD(NP);
    return out;
  }

  std::string sep1 = utf8_char(STRING_ELT(sep, 0));
  std::string sep2 = utf8_char(R_BlankString);
  std::string strng = utf8_char(R_BlankString);

  if (!is_null(collapse)){

    SHIELD(collapse = cast<r_character_t>(collapse, R_NilValue)); ++NP;
    if (Rf_length(collapse) != 1){
      YIELD(NP);
      Rf_error("`.collapse` must be a length 1 character vector in %s", __func__);
    }

    sep2 = utf8_char(STRING_ELT(collapse, 0));

    for (R_xlen_t j = 0; j < n_strings; ++j){
      if (j != 0) strng += sep2;
      for (R_xlen_t i = 0; i < n_chars; ++i){
        strng += utf8_char(char_ptrs[i][j]);
        if (i != (n_chars - 1)) strng += sep1;
      }
    }

    SEXP out = SHIELD(Rf_mkString(strng.c_str())); ++NP;

    YIELD(NP);
    return out;

  } else {

    SEXP out = SHIELD(new_vec(STRSXP, n_strings)); ++NP;

    // Along the length of each char vec

    for (R_xlen_t j = 0; j < n_strings; ++j){

      strng = utf8_char(char_ptrs[0][j]);

      for (R_xlen_t i = 1; i < n_chars; ++i){
        str_paste(strng, sep1, utf8_char(char_ptrs[i][j]));
      }
      SET_STRING_ELT(out, j, Rf_mkChar(strng.c_str()));
    }

    YIELD(NP);
    return out;
  }
}
