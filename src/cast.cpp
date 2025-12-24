#include "cheapr.h"

[[cpp11::init]]
void init_assert_symmetric_r_types(DllInfo* dll){
  for (int i = 0; i < n_types; ++i){
    for (int j = 0; j < n_types; ++j){
      if (r_type_pairs[i][j] != r_type_pairs[j][i]){
        Rf_error(
          "Internal error, `r_type_pairs[%d][%d]` returns %s while `r_type_pairs[%d][%d]` returns %s",
          i, j, r_type_names[r_type_pairs[i][j]], j, i, r_type_names[r_type_pairs[j][i]]
        );
      }
    }
  }
}

r_type r_common_type(SEXP x){

  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list");
  }

  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = list_ptr_ro(x);

  // Initialise to null
  r_type common = R_null;

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

  SEXP out = r_null;

  r_type common = r_common_type(x);

  const SEXP *p_x = list_ptr_ro(x);

  if (common == R_unk){

    SEXP vec_template;
    PROTECT_INDEX vec_template_idx;
    R_ProtectWithIndex(vec_template = r_null, &vec_template_idx); ++NP;

    for (R_xlen_t i = 0; i < n; ++i){
      R_Reprotect(vec_template = cast<r_unknown_t>(vec_template, p_x[i]), vec_template_idx);
    }

    out = vec_template;

  } else {

    // Initialise based on R type
    SHIELD(out = init_(common, 0, false)); ++NP;

    switch (common){

    // Combine levels
    case R_fct: {

      SEXP all_lvls, new_lvls;

      PROTECT_INDEX all_lvls_idx, new_lvls_idx;
      R_ProtectWithIndex(all_lvls = new_vector<r_string_t>(0), &all_lvls_idx); ++NP;
      R_ProtectWithIndex(new_lvls = new_vector<r_string_t>(0), &new_lvls_idx); ++NP;

      int n = Rf_length(x);

      for (int i = 0; i < n; ++i){
        if (Rf_isFactor(p_x[i])){
          R_Reprotect(new_lvls = get_attr(p_x[i], symbol::levels_sym), new_lvls_idx);
        } else {
          R_Reprotect(new_lvls = cast<r_characters_t>(p_x[i], r_null), new_lvls_idx);
        }
        if (!R_compute_identical(all_lvls, new_lvls, 0)){
          R_Reprotect(new_lvls = cpp_setdiff(new_lvls, all_lvls, false), new_lvls_idx);
          if (Rf_length(new_lvls) != 0){
            R_Reprotect(all_lvls = combine(all_lvls, new_lvls), all_lvls_idx);
          }
        }
      }
      set_attr(out, symbol::levels_sym, all_lvls);
      break;
    }
    // Figure out if date is integer or double
    case R_date: {

      if (Rf_xlength(x) > 0){
      int date_type = INTSXP;

      for (int i = 0; i < n; ++i){
        if (TYPEOF(p_x[i]) != INTSXP){
          date_type = REALSXP;
          break;
        }
      }
      SHIELD(out = vec::coerce_vec(out, date_type)); ++NP;
    }

      break;
    }
    case R_pxt: {
      // Initialised date-time gets timezone of first date-time
      if (Rf_xlength(x) > 0){
      SHIELD(out = cast<r_posixts_t>(out, p_x[0])); ++NP;
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
  const SEXP *p_x = list_ptr_ro(x);

  if (n <= 1){
    return x;
  }

  SEXP out = SHIELD(new_list(n)); ++NP;
  SEXP vec_template = SHIELD(cpp_common_template(x)); ++NP;
  r_type common = get_r_type(vec_template);

  switch (common){
  case R_null: {
    break;
  }
  case R_lgl: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_logicals_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_int: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_integers_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_int64: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_integers64_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_dbl: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_doubles_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_chr: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_characters_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_cplx: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_complexes_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_raw: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_raws_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_list: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_list_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_fct: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_factors_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_date: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_dates_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_pxt: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_posixts_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_df: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_data_frame_t>(p_x[i], vec_template));
  }
    break;
  }
  case R_unk: {
    for (R_xlen_t i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, cast<r_unknown_t>(p_x[i], vec_template));
  }
    break;
  }
    // This should never be reached because of the R_unk (unknown) case above
  default: {
    YIELD(NP);
    Rf_error("Unimplemented cast type");
  }
  }
  SEXP names = SHIELD(get_old_names(x)); ++NP;
  set_old_names(out, names);

  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_cast(SEXP x, SEXP y){
  return cast_(get_r_type(y), x, y);
}

[[cpp11::register]]
SEXP cpp_type(SEXP x){
  return as_vector(r_type_char(x));
}
