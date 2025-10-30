#include "cheapr.h"
#include <vector>

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

// Fast casting of objects to common type (via typeof)
SEXP fast_cast(SEXP x){

  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list");
  }

  R_xlen_t n = Rf_xlength(x);

  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = SHIELD(new_vec(VECSXP, n));

  std::vector<int16_t> types(n);

  int16_t out_type = 0;
  int16_t type;

  for (R_xlen_t i = 0; i < n; ++i){
    type = TYPEOF(p_x[i]);
    out_type = std::max(out_type, type);
    types[i] = type;
  }

  // Cast all objects
  for (R_xlen_t i = 0; i < n; ++i){
    if (types[i] != out_type){
      SET_VECTOR_ELT(out, i, coerce_vec(p_x[i], out_type));
    } else {
      SET_VECTOR_ELT(out, i, p_x[i]);
    }
  }

  YIELD(1);
  return out;
}

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

[[cpp11::register]]
SEXP cpp_cast_common(SEXP x){

  int32_t NP = 0;

  r_type common = r_common_type(x);

  R_xlen_t n = Rf_xlength(x);
  const SEXP *p_x = VECTOR_PTR_RO(x);

  if (n <= 1){
    return x;
  }

  SEXP out = SHIELD(new_vec(VECSXP, n)); ++NP;

  SEXP temp;
  PROTECT_INDEX temp_idx;
  R_ProtectWithIndex(temp = R_NilValue, &temp_idx); ++NP;

#define CAST_LOOP(cast_fn)                                     \
  for (R_xlen_t i = 0; i < (n - 1); ++i){                      \
    R_Reprotect(temp = cast_fn(p_x[i], p_x[i + 1]), temp_idx); \
    SET_VECTOR_ELT(out, i, temp);                              \
  }                                                            \
  SET_VECTOR_ELT(out, n - 1, cast_fn(p_x[n - 1], temp));


switch (common){
case r_null: {
  break;
}
case r_lgl: {
  CAST_LOOP(cast<r_logical_t>)
  break;
}
case r_int: {
  CAST_LOOP(cast<r_integer_t>)
  break;
}
case r_int64: {
  CAST_LOOP(cast<r_integer64_t>)
  break;
}
case r_dbl: {
  CAST_LOOP(cast<r_numeric_t>)
  break;
}
case r_chr: {
  CAST_LOOP(cast<r_character_t>)
  break;
}
case r_cplx: {
  CAST_LOOP(cast<r_complex_t>)
  break;
}
case r_raw: {
  CAST_LOOP(cast<r_raw_t>)
  break;
}
case r_list: {
  CAST_LOOP(cast<r_list_t>)
  break;
}
case r_fct: {

  SEXP lvls = SHIELD(combine_levels(x)); ++NP;

  for (R_xlen_t i = 0; i < n; ++i){
    R_Reprotect(temp = cast<r_character_t>(p_x[i], R_NilValue), temp_idx);
    SET_VECTOR_ELT(out, i, character_as_factor(temp, lvls));
  }
  break;
}
case r_date: {
  CAST_LOOP(cast<r_date_t>)
  break;
}
case r_pxct: {
  CAST_LOOP(cast<r_posixt_t>)
  break;
}
case r_rcrd: {
  CAST_LOOP(cast<r_vctrs_rcrd_t>)
  break;
}
case r_df: {
  CAST_LOOP(cast<r_data_frame_t>)
  break;
}
case r_unk: {
  CAST_LOOP(cast<r_unknown_t>)
  break;
}
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
