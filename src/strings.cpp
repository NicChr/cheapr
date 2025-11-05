#include "cheapr.h"

// Coalesce a list of string vectors

[[cpp11::register]]
SEXP cpp_str_coalesce(SEXP x){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of character vectors in %s", __func__);
  }

  int32_t NP = 0;
  uint_fast64_t n = Rf_xlength(x);
  uint_fast64_t out_size = 0;
  uint_fast64_t m;

  const SEXP *p_x = VECTOR_PTR_RO(x);
  std::vector<const SEXP*> str_ptrs(n);

  SEXP char_vec = R_NilValue;
  uint32_t xtype;

  bool shallow_duplicated = false;

  for (uint_fast64_t i = 0; i < n; ++i){
    char_vec = p_x[i];
    xtype = TYPEOF(char_vec);

    if (xtype != STRSXP){
      if (!shallow_duplicated){
        SHIELD(x = Rf_shallow_duplicate(x)); ++NP;
        p_x = VECTOR_PTR_RO(x);
        shallow_duplicated = true;
      }
      SET_VECTOR_ELT(x, i, base_as_character(char_vec));
      char_vec = p_x[i];
    }

    str_ptrs[i] = STRING_PTR_RO(char_vec);

    if (xtype != NILSXP){
      m = Rf_xlength(char_vec);
      if (m == 0){
        YIELD(NP);
        return Rf_allocVector(STRSXP, 0);
      }
      out_size = std::max(out_size, m);
    }
  }

  SEXP out = SHIELD(Rf_allocVector(STRSXP, out_size)); ++NP;

  SEXP inner_char = R_BlankString;

  uint_fast64_t n_nas;

  for (uint_fast64_t i = 0; i < out_size; ++i){
    n_nas = 0;
    for (uint_fast64_t j = 0; j < n; ++j){
      m = Rf_xlength(p_x[j]);
      if (m == 0) continue;
      inner_char = str_ptrs[j][i % m];
      n_nas += inner_char == NA_STRING;
      if (!(inner_char == R_BlankString || inner_char == NA_STRING)){
        SET_STRING_ELT(out, i, inner_char);
        break;
      }
      // If all ith elements are NA, then return NA
      if (n_nas == n){
        SET_STRING_ELT(out, i, NA_STRING);
      }
    }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_paste(SEXP x, std::string sep, SEXP collapse){

  if (TYPEOF(x) != VECSXP){
    Rf_error("`x` must be a list of character vectors in %s", __func__);
  }

  R_xlen_t n_chars = Rf_xlength(x);

  if (n_chars == 0){
    return new_vec(STRSXP, 0);
  }

  int32_t NP = 0;

  SEXP chars = SHIELD(cpp_recycle(x, R_NilValue)); ++NP;
  if (MAYBE_REFERENCED(chars)){
    SHIELD(chars = Rf_shallow_duplicate(chars)); ++NP;
  }

  const SEXP *p_chars = VECTOR_PTR_RO(chars);
  std::vector<const SEXP *> char_ptrs(n_chars);

  // First cast all to character vectors

  for (R_xlen_t i = 0; i < n_chars; ++i){
    SET_VECTOR_ELT(chars, i, cast<r_character_t>(p_chars[i], R_NilValue));
    char_ptrs[i] = STRING_PTR_RO(p_chars[i]);
  }

  R_xlen_t n_strings = Rf_xlength(VECTOR_ELT(chars, 0));

  std::string strng;

  if (!is_null(collapse)){

    SHIELD(collapse = cast<r_character_t>(collapse, R_NilValue)); ++NP;
    if (Rf_length(collapse) != 1){
      YIELD(NP);
      Rf_error("`.collapse` must be a length 1 character vector in %s", __func__);
    }

    sep += static_cast<std::string>(utf8_char(STRING_ELT(collapse, 0)));
    SEXP out = SHIELD(new_vec(STRSXP, 1)); ++NP;

    // Initialise string to first collapsed string

    strng = utf8_char(char_ptrs[0][0]);
    for (R_xlen_t i = 1; i < n_chars; ++i){
      str_paste(strng, sep, utf8_char(char_ptrs[i][0]));
    }

    for (R_xlen_t j = 1; j < n_strings; ++j){
      for (R_xlen_t i = 0; i < n_chars; ++i){
        str_paste(strng, sep, utf8_char(char_ptrs[i][j]));
      }
    }

    SET_STRING_ELT(out, 0, Rf_mkChar(strng.c_str()));

    YIELD(NP);
    return out;

  } else {

    SEXP out = SHIELD(new_vec(STRSXP, n_strings)); ++NP;

    // Along the length of each char vec

    for (R_xlen_t j = 0; j < n_strings; ++j){

      strng = utf8_char(char_ptrs[0][j]);

      for (R_xlen_t i = 1; i < n_chars; ++i){
        str_paste(strng, sep, utf8_char(char_ptrs[i][j]));
      }
      SET_STRING_ELT(out, j, Rf_mkChar(strng.c_str()));
    }

    YIELD(NP);
    return out;
  }

}
// SEXP cpp_paste(SEXP x, const std::string sep, SEXP collapse){
//
//   if (TYPEOF(x) != VECSXP){
//     Rf_error("`x` must be a list of character vectors in %s", __func__);
//   }
//
//   R_xlen_t n_chars = Rf_xlength(x);
//
//   if (n_chars == 0){
//     return new_vec(STRSXP, 0);
//   }
//
//   int32_t NP = 0;
//
//   SEXP chars = SHIELD(cpp_recycle(x, R_NilValue)); ++NP;
//   if (MAYBE_REFERENCED(chars)){
//     SHIELD(chars = Rf_shallow_duplicate(chars)); ++NP;
//   }
//
//   const SEXP *p_chars = VECTOR_PTR_RO(chars);
//   std::vector<const SEXP *> char_ptrs(n_chars);
//
//   // First cast all to character vectors
//
//   for (R_xlen_t i = 0; i < n_chars; ++i){
//     SET_VECTOR_ELT(chars, i, cast<r_character_t>(p_chars[i], R_NilValue));
//     char_ptrs[i] = STRING_PTR_RO(p_chars[i]);
//   }
//
//   R_xlen_t n_strings = Rf_xlength(VECTOR_ELT(chars, 0));
//
//   std::string strng;
//
//   if (collapse){
//
//     SEXP out = SHIELD(new_vec(STRSXP, 1)); ++NP;
//
//     // Initialise string to first collapsed string
//
//     strng = utf8_char(char_ptrs[0][0]);
//     for (R_xlen_t i = 1; i < n_chars; ++i){
//       str_paste(strng, sep, utf8_char(char_ptrs[i][0]));
//     }
//
//     for (R_xlen_t j = 1; j < n_strings; ++j){
//       for (R_xlen_t i = 0; i < n_chars; ++i){
//         str_paste(strng, sep, utf8_char(char_ptrs[i][j]));
//       }
//     }
//
//     SET_STRING_ELT(out, 0, Rf_mkChar(strng.c_str()));
//
//     YIELD(NP);
//     return out;
//
//   } else {
//
//     SEXP out = SHIELD(new_vec(STRSXP, n_strings)); ++NP;
//
//     // Along the length of each char vec
//
//     for (R_xlen_t j = 0; j < n_strings; ++j){
//
//       strng = utf8_char(char_ptrs[0][j]);
//
//       for (R_xlen_t i = 1; i < n_chars; ++i){
//         str_paste(strng, sep, utf8_char(char_ptrs[i][j]));
//       }
//       SET_STRING_ELT(out, j, Rf_mkChar(strng.c_str()));
//     }
//
//     YIELD(NP);
//     return out;
//   }
//
// }
