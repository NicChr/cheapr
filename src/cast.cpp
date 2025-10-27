#include "cheapr.h"
#include <string>
#include <unordered_map>
#include <functional>

static SEXP as_lgl = NULL;
static SEXP as_int = NULL;
static SEXP as_dbl = NULL;
static SEXP as_char = NULL;
static SEXP as_cplx = NULL;
static SEXP as_raw = NULL;
static SEXP as_date = NULL;
static SEXP as_posixct = NULL;
static SEXP as_list = NULL;

// r types

struct r_null_t {};
struct r_logical_t {};
struct r_integer_t {};
struct r_integer64_t {};
struct r_numeric_t {};
struct r_character_t {};
struct r_complex_t {};
struct r_raw_t {};
struct r_list_t {};
struct r_factor_t {};
struct r_date_t {};
struct r_posixt_t {};
struct r_data_frame_t {};
struct r_unknown_t {};


// r type constants
using r_type = uint8_t;
enum : r_type {
  r_null = (uint8_t) 0,
  r_lgl = (uint8_t) 1,
  r_int = (uint8_t) 2,
  r_int64 = (uint8_t) 3,
  r_dbl = (uint8_t) 4,
  r_chr = (uint8_t) 5,
  r_cplx = (uint8_t) 6,
  r_raw = (uint8_t) 7,
  r_list = (uint8_t) 8,
  r_fct = (uint8_t) 9,
  r_date = (uint8_t) 10,
  r_pxct = (uint8_t) 11,
  r_df = (uint8_t) 12,
  r_unk = (uint8_t) 13,
};

// Symmetric lattice with columns/rows in the new order:
// NULL, LGL, INT, I64, DBL, CHR, CPLX, RAW, LIST, FCT, DATE, PXCT, DF
constexpr uint8_t r_type_pairs[14][14] = {
  /*            NULL    LGL     INT     I64     DBL     CHR     CPLX    RAW     LIST    FCT     DATE    PXCT    DF   Unknown */
  /* NULL */  { r_null, r_lgl,  r_int,  r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* LGL  */  { r_lgl,  r_lgl,  r_int,  r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* INT  */  { r_int,  r_int,  r_int,  r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* I64  */  { r_int64,r_int64,r_int64,r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* DBL  */  { r_dbl,  r_dbl,  r_dbl,  r_dbl,  r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* CHR  */  { r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_list, r_fct,  r_chr,  r_chr,  r_df, r_unk },
  /* CPLX */  { r_cplx, r_cplx, r_cplx, r_cplx, r_cplx, r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_fct,  r_df,   r_df, r_unk },
  /* RAW  */  { r_raw,  r_raw,  r_raw,  r_raw,  r_raw,  r_chr,  r_raw,  r_raw,  r_list, r_df,   r_df,   r_df,   r_df, r_unk },
  /* LIST */  { r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_df, r_unk },
  /* FCT  */  { r_fct,  r_fct,  r_fct,  r_fct,  r_fct,  r_fct,  r_fct,  r_df,   r_list, r_fct,  r_fct,  r_fct,  r_df, r_unk },
  /* DATE */  { r_date, r_date, r_date, r_date, r_date, r_chr,  r_fct,  r_df,   r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* PXCT */  { r_pxct, r_pxct, r_pxct, r_pxct, r_pxct, r_chr,  r_df,   r_df,   r_list, r_fct,  r_pxct, r_pxct, r_df, r_unk },
  /* DF   */  { r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df, r_unk },
  /* Unknown   */  { r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk, r_unk }
};

using namespace cpp11;

// `class()`
const char* get_class(SEXP obj){
  if (Rf_isObject(obj)){
    SEXP klass = Rf_getAttrib(obj, R_ClassSymbol);

    int n = Rf_length(klass);

    return CHAR(STRING_ELT(klass, n - 1));
  } else {
    switch(TYPEOF(obj)) {
    case LGLSXP: {
      return "logical";
    }
    case INTSXP: {
      return "integer";
    }
    case REALSXP: {
      return "numeric";
    }
    case STRSXP: {
      return "character";
    }
    case CPLXSXP: {
      return "complex";
    }
    case RAWSXP: {
      return "raw";
    }
    case VECSXP: {
      return "list";
    }
    default: {
      return "Unknown";
    }
    }
  }
}

inline r_type common_type(const r_type &a, const r_type &b) {
  // return (a <= b) ? r_type_pairs[a][b] : r_type_pairs[b][a];
  return r_type_pairs[a][b];
}

// Convert single SEXP into r_* code.
inline const r_type get_r_type(SEXP x) {

  if (!Rf_isObject(x)){
    switch (TYPEOF(x)) {
    case NILSXP:  return r_null;
    case LGLSXP:  return r_lgl;
    case INTSXP:  return r_int;
    case REALSXP: return r_dbl;
    case STRSXP:  return r_chr;
    case CPLXSXP: return r_cplx;
    case RAWSXP:  return r_raw;
    case VECSXP:  return r_list;
    default: return r_unk;
    }
  } else {

    if (Rf_inherits(x, "data.frame")) return r_df;
    if (Rf_inherits(x, "POSIXct"))    return r_pxct;
    if (Rf_inherits(x, "Date"))       return r_date;
    if (Rf_inherits(x, "factor"))     return r_fct;
    if (Rf_inherits(x, "integer64"))  return r_int64;
    return r_unk;
  }
}

// cast template with specialisations
template<typename T>
SEXP cast(SEXP x, SEXP y) {
  stop("Unimplemented type conversion");
}

template<>
inline SEXP cast<r_null_t>(SEXP x, SEXP y) {
  return R_NilValue;
}

template<>
inline SEXP cast<r_logical_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "logical")){
    return x;
  } else if (Rf_isObject(x)){
    as_lgl = as_lgl != NULL ? as_lgl : Rf_install("as.logical");
    return Rf_eval(Rf_lang2(as_lgl, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, LGLSXP);
  }
}

template<>
inline SEXP cast<r_integer_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "integer")){
    return x;
  } else if (Rf_isObject(x)){
    as_int = as_int != NULL ? as_int : Rf_install("as.integer");
    return Rf_eval(Rf_lang2(as_int, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, INTSXP);
  }
}

template<>
inline SEXP cast<r_integer64_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "integer64")){
    return x;
  } else {
    return coerce_vector(x, CHEAPR_INT64SXP);
  }
}

template<>
inline SEXP cast<r_numeric_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "numeric")){
    return x;
  } else if (Rf_isObject(x)){
    as_dbl = as_dbl != NULL ? as_dbl : Rf_install("as.numeric");
    return Rf_eval(Rf_lang2(as_dbl, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, REALSXP);
  }
}

template<>
inline SEXP cast<r_character_t>(SEXP x, SEXP y) {

  if (Rf_inherits(x, "character")){
    return x;
  } else if (Rf_isFactor(x)){
    return factor_as_character(x);
  } else if (Rf_isObject(x)){
    as_char = as_char != NULL ? as_char : Rf_install("as.character");
    return Rf_eval(Rf_lang2(as_char, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, STRSXP);
  }
}

template<>
inline SEXP cast<r_complex_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "complex")){
    return x;
  } else if (Rf_isObject(x)){
    as_cplx = as_cplx != NULL ? as_cplx : Rf_install("as.complex");
    return Rf_eval(Rf_lang2(as_cplx, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, CPLXSXP);
  }
}

template<>
inline SEXP cast<r_raw_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "raw")){
    return x;
  } else if (Rf_isObject(x)){
    as_raw = as_raw != NULL ? as_raw : Rf_install("as.raw");
    return Rf_eval(Rf_lang2(as_raw, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, RAWSXP);
  }
}

template<>
inline SEXP cast<r_list_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "list")){
    return x;
  } else if (Rf_isObject(x)){
    as_list = as_list != NULL ? as_list : Rf_install("as.list");
    return Rf_eval(Rf_lang2(as_list, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, VECSXP);
  }
}

template<>
inline SEXP cast<r_factor_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "factor") && !Rf_inherits(y, "factor")){
    return x;
  } else if (Rf_inherits(y, "factor")){
    SEXP fctrs = SHIELD(new_vec(VECSXP, 2));
    SET_VECTOR_ELT(fctrs, 0, x);
    SET_VECTOR_ELT(fctrs, 1, y);
    SEXP all_levels = SHIELD(cpp_combine_levels(fctrs));
    SEXP out = SHIELD(cheapr_factor(x, cpp11::named_arg("levels") = all_levels));
    YIELD(3);
    return out;
  } else {
    return cheapr_factor(x);
  }
}

template<>
inline SEXP cast<r_date_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "Date")){
    return x;
  } else if (Rf_isObject(x)){
    as_date = as_date != NULL ? as_date : Rf_install("as.date");
    return Rf_eval(Rf_lang2(as_date, x), R_GetCurrentEnv());
  } else {

    int32_t NP = 0;

    SEXP out = SHIELD(Rf_shallow_duplicate(x)); ++NP;
    if (TYPEOF(x) != INTSXP){
      SHIELD(out = coerce_vec(x, REALSXP)); ++NP;
    }
    Rf_classgets(out, make_utf8_str("Date"));
    YIELD(NP);
    return out;
  }
}

template<>
inline SEXP cast<r_posixt_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "POSIXct")){
    return x;
  } else if (Rf_inherits(x, "Date") && Rf_inherits(y, "POSIXct")){
    R_xlen_t n = Rf_xlength(x);
    SEXP out = SHIELD(new_vec(REALSXP, n));
    SEXP out_class = SHIELD(new_vec(STRSXP, 2));
    SEXP out_tzone = SHIELD(Rf_getAttrib(y, install_utf8("tzone")));

    SET_STRING_ELT(out_class, 0, make_utf8_char("POSIXct"));
    SET_STRING_ELT(out_class, 1, make_utf8_char("POSIXt"));
    Rf_classgets(out, out_class);
    Rf_setAttrib(out, install_utf8("tzone"), out_tzone);

    double* RESTRICT p_out = REAL(out);

    if (TYPEOF(x) == INTSXP){
      const int *p_x = INTEGER_RO(x);
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = is_na_int(p_x[i]) ? NA_REAL : static_cast<double>(p_x[i]) * 86400;

    } else {
      const double *p_x = REAL_RO(x);
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = p_x[i] * 86400;
    }

    YIELD(3);
    return out;
  } else {
    as_posixct = as_posixct != NULL ? as_posixct : Rf_install("as.POSIXct");
    return Rf_eval(Rf_lang2(as_posixct, x), R_GetCurrentEnv());
  }
}

template<>
inline SEXP cast<r_data_frame_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "data.frame") && Rf_inherits(y, "data.frame")){
    return rebuild(x, y, true);
  } else if (is_simple_atomic_vec2(x)){
    SEXP out = SHIELD(new_vec(VECSXP, 1));
    SET_VECTOR_ELT(out, 0, x);
    SEXP names = SHIELD(make_utf8_str("x"));
    set_names(out, names);
    set_list_as_df(out); // as data_frame in-place
    YIELD(2);
    return out;
  } else {
    int32_t NP = 0;
    SEXP out = SHIELD(cheapr_as_df(x)); ++NP;

    if (Rf_inherits(y, "data.frame")){
      SHIELD(out = rebuild(out, y, true)); ++NP;
    }

    YIELD(NP);
    return out;
  }
}

template<>
inline SEXP cast<r_unknown_t>(SEXP x, SEXP y) {
  return cheapr_cast(x, y);
}


// Wrapper functions for cast fns map
using cast_fn = SEXP(*)(SEXP, SEXP);
inline SEXP cast_null(SEXP x, SEXP y) { return cast<r_null_t>(x, y); }
inline SEXP cast_logical(SEXP x, SEXP y) { return cast<r_logical_t>(x, y); }
inline SEXP cast_integer(SEXP x, SEXP y) { return cast<r_integer_t>(x, y); }
inline SEXP cast_integer64(SEXP x, SEXP y) { return cast<r_integer64_t>(x, y); }
inline SEXP cast_numeric(SEXP x, SEXP y) { return cast<r_numeric_t>(x, y); }
inline SEXP cast_character(SEXP x, SEXP y) { return cast<r_character_t>(x, y); }
inline SEXP cast_complex(SEXP x, SEXP y) { return cast<r_complex_t>(x, y); }
inline SEXP cast_raw(SEXP x, SEXP y) { return cast<r_raw_t>(x, y); }
inline SEXP cast_list(SEXP x, SEXP y) { return cast<r_list_t>(x, y); }
inline SEXP cast_factor(SEXP x, SEXP y) { return cast<r_factor_t>(x, y); }
inline SEXP cast_date(SEXP x, SEXP y) { return cast<r_date_t>(x, y); }
inline SEXP cast_posixt(SEXP x, SEXP y) { return cast<r_posixt_t>(x, y); }
inline SEXP cast_data_frame(SEXP x, SEXP y) { return cast<r_data_frame_t>(x, y); }

static const cast_fn CAST_FNS[13] = {
  cast_null, cast_logical, cast_integer, cast_integer64,
  cast_numeric, cast_character, cast_complex, cast_raw, cast_list,
  cast_factor, cast_date, cast_posixt, cast_data_frame
};

// Dispatcher function
inline SEXP cast_(r_type cast_type, SEXP x, SEXP y) {
  return CAST_FNS[cast_type](x, y);
}

[[cpp11::register]]
SEXP cpp_cast(SEXP x, SEXP y) {
  r_type common = common_type(get_r_type(x), get_r_type(y));
  return cast_(common, x, y);
}

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
    if (common == r_unk) break;
  }
  return common;
}

[[cpp11::register]]
SEXP cpp_cast_all(SEXP x){

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

    SEXP lvls = SHIELD(cpp_combine_levels(x)); ++NP;
    SEXP fctr_cls = SHIELD(make_utf8_str("factor")); ++NP;
    // R_Reprotect(temp = cheapr_factor(cpp11::named_arg("levels") = lvls), temp_idx);


    for (R_xlen_t i = 0; i < n; ++i){
      R_Reprotect(temp = cast<r_character_t>(p_x[i], R_NilValue), temp_idx);
      R_Reprotect(temp = Rf_match(lvls, temp, NA_INTEGER), temp_idx);
      Rf_setAttrib(temp, R_LevelsSymbol, lvls);
      Rf_classgets(temp, fctr_cls);
      SET_VECTOR_ELT(out, i, temp);
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
