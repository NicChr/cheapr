#include "cheapr.h"

SEXP as_lgl = NULL;
SEXP as_int = NULL;
SEXP as_dbl = NULL;
SEXP as_char = NULL;
SEXP as_cplx = NULL;
SEXP as_raw = NULL;
SEXP as_date = NULL;
SEXP as_posixct = NULL;
SEXP as_list = NULL;

// The order of functions here MUST MATCH the order of defined r types
const cast_fn CAST_FNS[15] = {
  cast_null, cast_logical, cast_integer, cast_integer64,
  cast_numeric, cast_complex, cast_raw, cast_date, cast_posixt,
  cast_vctrs_rcrd, cast_character, cast_factor, cast_list,
  cast_data_frame, cast_unknown
};

const init_fn INIT_FNS[15] = {
  init_null, init_logical, init_integer, init_integer64,
  init_numeric, init_complex, init_raw, init_date, init_posixt,
  init_vctrs_rcrd, init_character, init_factor, init_list,
  init_data_frame, init_unknown
};

r_type r_common_type(SEXP x){

  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list");
  }

  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

  // Initialise to null
  r_type common = 0;

  for (R_xlen_t i = 0; i < n; ++i){
    common = common_type(common, get_r_type(p_x[i]));
  }
  return common;
}

// Return common template from a list of vectors
[[cpp11::register]]
SEXP cpp_common_template(SEXP x){

  R_xlen_t n = Rf_xlength(x);
  int32_t NP = 0;

  SEXP out = R_NilValue;

  r_type common = r_common_type(x);

  if (common == r_unk){

    const SEXP *p_x = VECTOR_PTR_RO(x);

    SEXP vec_template;
    PROTECT_INDEX vec_template_idx;
    R_ProtectWithIndex(vec_template = R_NilValue, &vec_template_idx); ++NP;

    for (R_xlen_t i = 0; i < n; ++i){
      R_Reprotect(vec_template = cast<r_unknown_t>(vec_template, p_x[i]), vec_template_idx);
    }

    out = vec_template;

  } else {
    SHIELD(out = init_(common, 0, false)); ++NP;

    switch (common){
    case r_fct: {

      SEXP all_lvls, new_lvls;

      PROTECT_INDEX all_lvls_idx, new_lvls_idx;
      R_ProtectWithIndex(all_lvls = new_vec(STRSXP, 0), &all_lvls_idx); ++NP;
      R_ProtectWithIndex(new_lvls = new_vec(STRSXP, 0), &new_lvls_idx); ++NP;

      const SEXP *p_x = VECTOR_PTR_RO(x);
      int n = Rf_length(x);

      for (int i = 0; i < n; ++i){
        if (Rf_isFactor(p_x[i])){
          R_Reprotect(new_lvls = Rf_getAttrib(p_x[i], R_LevelsSymbol), new_lvls_idx);
        } else {
          R_Reprotect(new_lvls = cast<r_character_t>(p_x[i], R_NilValue), new_lvls_idx);
        }
        if (!R_compute_identical(all_lvls, new_lvls, 0)){
          R_Reprotect(new_lvls = cpp_setdiff(new_lvls, all_lvls, false), new_lvls_idx);
          if (Rf_length(new_lvls) != 0){
            R_Reprotect(all_lvls = c2(all_lvls, new_lvls), all_lvls_idx);
          }
        }
      }
      Rf_setAttrib(out, R_LevelsSymbol, all_lvls);
      break;
    }
    case r_pxct: {
      if (Rf_xlength(x) > 0){
      SHIELD(out = cast<r_posixt_t>(out, VECTOR_ELT(x, 0))); ++NP;
    }
      break;
    }
    default: {
      break;
    }
    }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_cast_common(SEXP x){

  int32_t NP = 0;

  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

  if (n <= 1){
    return x;
  }

  SEXP out = SHIELD(new_vec(VECSXP, n)); ++NP;
  SEXP vec_template = SHIELD(cpp_common_template(x)); ++NP;
  r_type common = get_r_type(vec_template);

  switch (common){
  case r_null: {
    break;
  }
  case r_lgl: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_logical_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_int: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_integer_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_int64: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_integer64_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_dbl: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_numeric_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_chr: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_character_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_cplx: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_complex_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_raw: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_raw_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_list: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_list_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_fct: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_factor_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_date: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_date_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_pxct: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_posixt_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_rcrd: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_vctrs_rcrd_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_df: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_data_frame_t>(p_x[i], vec_template));
  }
    break;
  }
  case r_unk: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_unknown_t>(p_x[i], vec_template));
  }
    break;
  }
    // This should never be reached because of the r_unk (unknown) case above
  default: {
    YIELD(NP);
    Rf_error("Unimplemented cast type");
  }
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_cast(SEXP x, SEXP y){
  return cast_(get_r_type(y), x, y);
}

[[cpp11::register]]
SEXP cpp_type(SEXP x){
  return make_utf8_str(r_type_char(x));
}
