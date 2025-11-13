#ifndef CHEAPR_CAST_H
#define CHEAPR_CAST_H

#include "cheapr_core.h"
#include "types.h"
#include "variadic.h"

// Defining custom R types
// for casting, initialising, combining and assigning

// Symbols for R conversion fns

inline SEXP as_lgl = R_NilValue;
inline SEXP as_int = R_NilValue;
inline SEXP as_dbl = R_NilValue;
inline SEXP as_char = R_NilValue;
inline SEXP as_cplx = R_NilValue;
inline SEXP as_raw = R_NilValue;
inline SEXP as_date = R_NilValue;
inline SEXP as_posixct = R_NilValue;
inline SEXP as_list = R_NilValue;

// Custom r types

struct r_null_t {};
struct r_logical_t {};
struct r_integer_t {};
struct r_integer64_t {};
struct r_numeric_t {};
struct r_complex_t {};
struct r_raw_t {};
struct r_date_t {};
struct r_posixt_t {};
struct r_character_t {};
struct r_factor_t {};
struct r_list_t {};
struct r_data_frame_t {};
struct r_unknown_t {};

// r type constants
enum : cheapr::r_type {
  r_null = 0,
    r_lgl = 1,
    r_int = 2,
    r_int64 = 3,
    r_dbl = 4,
    r_cplx = 5,
    r_raw = 6,
    r_date = 7,
    r_pxct = 8,
    r_chr = 9,
    r_fct = 10,
    r_list = 11,
    r_df = 12,
    r_unk = 13,
    r_err = 14 // Special type to signal incompatible cast (Currently unused)
};

// R type chars
inline constexpr const char* r_type_names[15] = {
  "NULL",        // 0
  "logical",     // 1
  "integer",     // 2
  "integer64",   // 3
  "numeric",     // 4
  "complex",     // 5
  "raw",         // 6
  "Date",        // 7
  "POSIXct",     // 8
  "character",  // 10
  "factor",     // 11
  "list",       // 12
  "data.frame", // 13
  "unknown"     // 14
};

// An n x n matrix of r types and their common cast type

inline constexpr cheapr::r_type r_type_pairs[14][14] = {
  /*            NULL    LGL     INT     I64     DBL     CPLX    RAW     DATE    PXCT    RCRD    CHR     FCT     LIST    DF      Unknown */
  /* NULL */  { r_null, r_lgl,  r_int,  r_int64, r_dbl,  r_cplx, r_raw,  r_date, r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* LGL  */  { r_lgl,  r_lgl,  r_int,  r_int64, r_dbl,  r_cplx, r_raw,  r_date, r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* INT  */  { r_int,  r_int,  r_int,  r_int64, r_dbl,  r_cplx, r_raw,  r_date, r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* I64  */  { r_int64,r_int64, r_int64, r_int64, r_dbl, r_cplx, r_raw, r_date, r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* DBL  */  { r_dbl,  r_dbl,  r_dbl,  r_dbl,   r_dbl,  r_cplx, r_raw, r_date, r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* CPLX */  { r_cplx, r_cplx, r_cplx, r_cplx,  r_cplx, r_cplx, r_raw, r_date,  r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* RAW  */  { r_raw,  r_raw,  r_raw,  r_raw,   r_raw,  r_raw,  r_raw, r_unk,   r_unk, r_chr,  r_fct,   r_list, r_df,   r_unk },
  /* DATE */  { r_date, r_date, r_date, r_date,  r_date, r_date,  r_unk, r_date, r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* PXCT */  { r_pxct, r_pxct, r_pxct, r_pxct,  r_pxct, r_pxct, r_unk, r_pxct, r_pxct, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* CHR  */  { r_chr,  r_chr,  r_chr,  r_chr,   r_chr,  r_chr,  r_chr,  r_chr,  r_chr, r_chr,  r_fct,  r_list, r_df,   r_unk },
  /* FCT  */  { r_fct,  r_fct,  r_fct,  r_fct,   r_fct,  r_fct,  r_df,   r_fct,  r_fct, r_fct,  r_fct,  r_list, r_df,   r_unk },
  /* LIST */  { r_list, r_list, r_list, r_list,  r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_df,   r_unk },
  /* DF   */  { r_df,   r_df,   r_df,   r_df,    r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_unk },
  /* Unknown */ { r_unk,  r_unk,  r_unk,  r_unk,   r_unk,  r_unk,  r_unk,  r_unk,  r_unk,  r_unk,  r_unk,  r_unk,  r_unk,  r_unk }
};

inline cheapr::r_type common_type(const cheapr::r_type &a, const cheapr::r_type &b) {
  return r_type_pairs[a][b];
}

// Convert single SEXP into r_* code.
inline cheapr::r_type get_r_type(SEXP x) {

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

    if (Rf_inherits(x, "factor"))     return r_fct;
    if (Rf_inherits(x, "Date"))       return r_date;
    if (Rf_inherits(x, "POSIXct"))    return r_pxct;
    if (Rf_inherits(x, "data.frame")) return r_df;
    if (Rf_inherits(x, "integer64"))  return r_int64;
    return r_unk;
  }
}

inline const char *r_type_char(SEXP x){

  cheapr::r_type type = get_r_type(x);

  // If unknown type
  if (type == r_unk){
    return cheapr::r_class(x);
  } else {
    return r_type_names[type];
  }
}

// initialise template with specialisations
template<typename T>
inline SEXP init(R_xlen_t n, bool with_na) {
  Rf_error("Unimplemented initialisation");
  return R_NilValue;
}

template<>
inline SEXP init<r_null_t>(R_xlen_t n, bool with_na) {
  return R_NilValue;
}

template<>
inline SEXP init<r_logical_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    SEXP out = cheapr::SHIELD(cheapr::new_vec(LGLSXP, n));
    int* RESTRICT p_out = INTEGER(out);
    std::fill(p_out, p_out + n, NA_LOGICAL);
    cheapr::YIELD(1);
    return out;
  } else {
    return cheapr::new_vec(LGLSXP, n);
  }

}

template<>
inline SEXP init<r_integer_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    SEXP out = cheapr::SHIELD(cheapr::new_vec(INTSXP, n));
    int* RESTRICT p_out = INTEGER(out);
    std::fill(p_out, p_out + n, NA_INTEGER);
    cheapr::YIELD(1);
    return out;
  } else {
    return cheapr::new_vec(INTSXP, n);
  }
}

template<>
inline SEXP init<r_integer64_t>(R_xlen_t n, bool with_na) {
  SEXP out = cheapr::SHIELD(cheapr::new_vec(REALSXP, n));
  if (with_na){
    int64_t* RESTRICT p_out = cheapr::INTEGER64_PTR(out);
    std::fill(p_out, p_out + n, cheapr::NA_INTEGER64);
  }
  Rf_classgets(out, cheapr::make_utf8_str("integer64"));
  cheapr::YIELD(1);
  return out;
}

template<>
inline SEXP init<r_numeric_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    SEXP out = cheapr::SHIELD(cheapr::new_vec(REALSXP, n));
    double* RESTRICT p_out = REAL(out);
    std::fill(p_out, p_out + n, NA_REAL);
    cheapr::YIELD(1);
    return out;
  } else {
    return cheapr::new_vec(REALSXP, n);
  }

}

template<>
inline SEXP init<r_character_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    SEXP out = cheapr::SHIELD(cheapr::new_vec(STRSXP, 0));
    cheapr::SHIELD(out = cpp_na_init(out, n));
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::new_vec(STRSXP, n);
  }
}

template<>
inline SEXP init<r_complex_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    SEXP out = cheapr::SHIELD(cheapr::new_vec(CPLXSXP, 0));
    cheapr::SHIELD(out = cpp_na_init(out, n));
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::new_vec(CPLXSXP, n);
  }

}

template<>
inline SEXP init<r_raw_t>(R_xlen_t n, bool with_na) {
  return cheapr::new_vec(RAWSXP, n);
}

template<>
inline SEXP init<r_list_t>(R_xlen_t n, bool with_na) {
  return cheapr::new_vec(VECSXP, n);
}

template<>
inline SEXP init<r_factor_t>(R_xlen_t n, bool with_na) {
  SEXP out = cheapr::SHIELD(init<r_integer_t>(n, with_na));
  SEXP lvls = cheapr::SHIELD(cheapr::new_vec(STRSXP, 0));
  SEXP cls = cheapr::SHIELD(cheapr::make_utf8_str("factor"));
  Rf_setAttrib(out, R_LevelsSymbol, lvls);
  Rf_classgets(out, cls);
  cheapr::YIELD(3);
  return out;
}

template<>
inline SEXP init<r_date_t>(R_xlen_t n, bool with_na){
  SEXP out = cheapr::SHIELD(init<r_numeric_t>(n, with_na));
  SEXP cls = cheapr::SHIELD(cheapr::make_utf8_str("Date"));
  Rf_classgets(out, cls);
  cheapr::YIELD(2);
  return out;
}

template<>
inline SEXP init<r_posixt_t>(R_xlen_t n, bool with_na) {
  SEXP out = cheapr::SHIELD(init<r_numeric_t>(n, with_na));
  SEXP tz = cheapr::SHIELD(cheapr::new_vec(STRSXP, 1));
  SEXP cls = cheapr::SHIELD(cheapr::make_r_chars("POSIXct", "POSIXt"));
  Rf_classgets(out, cls);
  Rf_setAttrib(out, cheapr::install_utf8("tzone"), tz);
  cheapr::YIELD(3);
  return out;
}

template<>
inline SEXP init<r_data_frame_t>(R_xlen_t n, bool with_na) {
  SEXP out = cheapr::SHIELD(cheapr::new_vec(VECSXP, 0));
  cheapr::SHIELD(out = cpp_new_df(out, R_NilValue, false, false));
  cheapr::SHIELD(out = cpp_na_init(out, n));
  cheapr::YIELD(3);
  return out;
}

template<>
inline SEXP init<r_unknown_t>(R_xlen_t n, bool with_na) {
  Rf_error("Don't know how to initialise unknown type");
  return R_NilValue;
}

inline void check_casted_length(SEXP x, SEXP out){

  if (cheapr::vector_length(x) != cheapr::vector_length(out)){
    Rf_error(
      "Bad cast from type %s to type %s. \n`vector_length(x)` %lld doesn't match returned length of %lld",
      r_type_char(x),
      r_type_char(out),
      static_cast<long long int>(cheapr::vector_length(x)),
      static_cast<long long int>(cheapr::vector_length(out))
    );
  }
}

inline void signal_bad_cast(SEXP from, SEXP to){
  Rf_error(
    "Don't know how to cast from type %s to type %s",
    r_type_char(from), r_type_char(to)
  );
}

// cast template with specialisations
template<typename T>
inline SEXP cast(SEXP x, SEXP y) {
  Rf_error("Unimplemented cast specialisation");
  return R_NilValue;
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
    as_lgl = !cheapr::is_null(as_lgl) ? as_lgl : Rf_install("as.logical");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_lgl, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::coerce_vec(x, LGLSXP);
  }
}

template<>
inline SEXP cast<r_integer_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "integer")){
    return x;
  } else if (Rf_isObject(x)){
    as_int = !cheapr::is_null(as_int) ? as_int : Rf_install("as.integer");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_int, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::coerce_vec(x, INTSXP);
  }
}

template<>
inline SEXP cast<r_numeric_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "numeric")){
    return x;
  } else if (Rf_isObject(x)){
    as_dbl = !cheapr::is_null(as_dbl) ? as_dbl : Rf_install("as.double");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_dbl, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::coerce_vec(x, REALSXP);
  }
}

template<>
inline SEXP cast<r_integer64_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "integer64")){
    return x;
  } else {
    SEXP out = cheapr::SHIELD(cast<r_numeric_t>(x, R_NilValue));
    cheapr::SHIELD(out = cpp_numeric_to_int64(x));
    cheapr::YIELD(2);
    return out;
  }
}

template<>
inline SEXP cast<r_character_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "character")){
    return x;
  } else if (Rf_inherits(x, "factor")){
    return factor_as_character(x);
  } else if (Rf_isObject(x)){
    as_char = !cheapr::is_null(as_char) ? as_char : Rf_install("as.character");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_char, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::coerce_vec(x, STRSXP);
  }
}

template<>
inline SEXP cast<r_complex_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "complex")){
    return x;
  } else if (Rf_isObject(x)){
    as_cplx = !cheapr::is_null(as_cplx) ? as_cplx : Rf_install("as.complex");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_cplx, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::coerce_vec(x, CPLXSXP);
  }
}

template<>
inline SEXP cast<r_raw_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "raw")){
    return x;
  } else if (Rf_isObject(x)){
    as_raw = !cheapr::is_null(as_raw) ? as_raw : Rf_install("as.raw");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_raw, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::coerce_vec(x, RAWSXP);
  }
}

template<>
inline SEXP cast<r_list_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "list")){
    return x;
  } else if (Rf_isObject(x)){
    as_list = !cheapr::is_null(as_list) ? as_list : Rf_install("as.list");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_list, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {
    return cheapr::coerce_vec(x, VECSXP);
  }
}

template<>
inline SEXP cast<r_factor_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "factor") && !Rf_inherits(y, "factor")){
    return x;
  } else if (Rf_inherits(y, "factor")){
    SEXP x_lvls = cheapr::SHIELD(Rf_getAttrib(x, R_LevelsSymbol));
    SEXP out_lvls = cheapr::SHIELD(Rf_getAttrib(y, R_LevelsSymbol));
    if (R_compute_identical(x_lvls, out_lvls, 0)){
      cheapr::YIELD(2);
      return x;
    }
    SEXP out = cheapr::SHIELD(cast<r_character_t>(x, R_NilValue));
    cheapr::SHIELD(out = character_as_factor(out, out_lvls));
    cheapr::YIELD(4);
    return out;
  } else if (cheapr::is_null(x)){
    return init<r_factor_t>(0, true);
  } else {
    SEXP out = cheapr::SHIELD(cheapr::cheapr_factor(x));
    check_casted_length(x, out);
    cheapr::YIELD(1);
    return out;
  }
}

template<>
inline SEXP cast<r_date_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "Date") && !Rf_inherits(y, "Date")){
    return x;
  } else if  (Rf_inherits(x, "Date") && Rf_inherits(y, "Date")){
    if (TYPEOF(x) == TYPEOF(y)){
      return x;
    } else {
      return cheapr::coerce_vec(x, TYPEOF(y));
    }
  } else if (cheapr::is_null(x) && cheapr::is_null(y)){
    return init<r_date_t>(cheapr::vector_length(x), true);
  } else if (Rf_isObject(x)){
    as_date = !cheapr::is_null(as_date) ? as_date : Rf_install("as.Date");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_date, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::YIELD(2);
    return out;
  } else {

    int32_t NP = 0;

    SEXP out = cheapr::SHIELD(Rf_shallow_duplicate(x)); ++NP;
    if (TYPEOF(x) != INTSXP){
      cheapr::SHIELD(out = cheapr::coerce_vec(x, REALSXP)); ++NP;
    }
    Rf_classgets(out, cheapr::make_utf8_str("Date"));
    cheapr::YIELD(NP);
    return out;
  }
}

template<>
inline SEXP cast<r_posixt_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "POSIXct") && !Rf_inherits(y, "POSIXct")){
    return x;
    // Copy timezone information
  } else if (Rf_inherits(x, "POSIXct") && Rf_inherits(y, "POSIXct")){
    SEXP x_tzone = cheapr::SHIELD(Rf_getAttrib(x, cheapr::install_utf8("tzone")));
    SEXP out_tzone = cheapr::SHIELD(Rf_getAttrib(y, cheapr::install_utf8("tzone")));

    if (R_compute_identical(x_tzone, out_tzone, 0)){
      cheapr::YIELD(2);
      return x;
    }
    SEXP out = cheapr::SHIELD(Rf_shallow_duplicate(x));
    Rf_setAttrib(out, cheapr::install_utf8("tzone"), out_tzone);
    cheapr::YIELD(3);
    return out;
  } else if (cheapr::is_null(x) && cheapr::is_null(y)){
    return init<r_posixt_t>(0, true);
    // Fast method for converting into a date into a date-time
  } else if (Rf_inherits(x, "Date") && Rf_inherits(y, "POSIXct")){
    R_xlen_t n = Rf_xlength(x);
    SEXP out = cheapr::SHIELD(cheapr::new_vec(REALSXP, n));
    SEXP out_class = cheapr::SHIELD(cheapr::make_r_chars("POSIXct", "POSIXt"));
    SEXP out_tzone = cheapr::SHIELD(Rf_getAttrib(y, cheapr::install_utf8("tzone")));
    Rf_classgets(out, out_class);
    Rf_setAttrib(out, cheapr::install_utf8("tzone"), out_tzone);

    double* RESTRICT p_out = REAL(out);

    if (TYPEOF(x) == INTSXP){
      const int *p_x = INTEGER_RO(x);
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = cheapr::is_na<int>(p_x[i]) ? NA_REAL : static_cast<double>(p_x[i]) * 86400;

    } else {
      const double *p_x = REAL_RO(x);
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = p_x[i] * 86400;
    }

    cheapr::YIELD(3);
    return out;
  } else if (!Rf_isObject(x)){
    SEXP out = cheapr::SHIELD(cheapr::coerce_vec(x, REALSXP));
    SEXP out_class = cheapr::SHIELD(cheapr::make_r_chars("POSIXct", "POSIXt"));
    SEXP out_tzone = cheapr::SHIELD(cheapr::new_vec(STRSXP, 1));
    Rf_classgets(out, out_class);
    Rf_setAttrib(out, cheapr::install_utf8("tzone"), out_tzone);
    cheapr::SHIELD(out = cast<r_posixt_t>(out, y)); // To set the correct attributes
    cheapr::YIELD(4);
    return out;
  } else {
    as_posixct = !cheapr::is_null(as_posixct) ? as_posixct : Rf_install("as.POSIXct");
    SEXP expr = cheapr::SHIELD(Rf_lang2(as_posixct, x));
    SEXP out = cheapr::SHIELD(Rf_eval(expr, R_GetCurrentEnv()));
    check_casted_length(x, out);
    cheapr::SHIELD(out = cast<r_posixt_t>(out, y)); // To set the correct attributes
    cheapr::YIELD(3);
    return out;
  }
}

template<>
inline SEXP cast<r_data_frame_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "data.frame")){
    if (!Rf_inherits(y, "data.frame")){
      return x;
    } else {
      return rebuild(x, y, static_cast<bool>(MAYBE_REFERENCED(x)));
    }
  }

  int32_t NP = 0;
  SEXP out = cheapr::SHIELD(cpp_as_df(x)); ++NP;

  if (Rf_inherits(y, "data.frame")){
    cheapr::SHIELD(out = rebuild(out, y, static_cast<bool>(MAYBE_REFERENCED(out)))); ++NP;
  }
  cheapr::YIELD(NP);
  return out;
}

template<>
inline SEXP cast<r_unknown_t>(SEXP x, SEXP y) {
  if (std::strcmp(cheapr::r_class(x), cheapr::r_class(y)) == 0){
    return x;
  } else {
    SEXP out = cheapr::SHIELD(cheapr::base_cast(x, y));
    check_casted_length(x, out);
    cheapr::YIELD(1);
    return out;
  }
}


// Wrapper functions for cast fns
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
inline SEXP cast_unknown(SEXP x, SEXP y) { return cast<r_unknown_t>(x, y); }

// Wrapper functions for init fns
using init_fn = SEXP(*)(R_xlen_t, bool);
inline SEXP init_null(R_xlen_t n, bool with_na) { return init<r_null_t>(n, with_na); }
inline SEXP init_logical(R_xlen_t n, bool with_na) { return init<r_logical_t>(n, with_na); }
inline SEXP init_integer(R_xlen_t n, bool with_na) { return init<r_integer_t>(n, with_na); }
inline SEXP init_integer64(R_xlen_t n, bool with_na) { return init<r_integer64_t>(n, with_na); }
inline SEXP init_numeric(R_xlen_t n, bool with_na) { return init<r_numeric_t>(n, with_na); }
inline SEXP init_character(R_xlen_t n, bool with_na) { return init<r_character_t>(n, with_na); }
inline SEXP init_complex(R_xlen_t n, bool with_na) { return init<r_complex_t>(n, with_na); }
inline SEXP init_raw(R_xlen_t n, bool with_na) { return init<r_raw_t>(n, with_na); }
inline SEXP init_list(R_xlen_t n, bool with_na) { return init<r_list_t>(n, with_na); }
inline SEXP init_factor(R_xlen_t n, bool with_na) { return init<r_factor_t>(n, with_na); }
inline SEXP init_date(R_xlen_t n, bool with_na) { return init<r_date_t>(n, with_na); }
inline SEXP init_posixt(R_xlen_t n, bool with_na) { return init<r_posixt_t>(n, with_na); }
inline SEXP init_data_frame(R_xlen_t n, bool with_na) { return init<r_data_frame_t>(n, with_na); }
inline SEXP init_unknown(R_xlen_t n, bool with_na) { return init<r_unknown_t>(n, with_na); }

// The order of functions here MUST MATCH the order of defined r types
inline const cast_fn CAST_FNS[14] = {
  cast_null, cast_logical, cast_integer, cast_integer64,
  cast_numeric, cast_complex, cast_raw, cast_date, cast_posixt,
  cast_character, cast_factor, cast_list,
  cast_data_frame, cast_unknown
};

inline const init_fn INIT_FNS[14] = {
  init_null, init_logical, init_integer, init_integer64,
  init_numeric, init_complex, init_raw, init_date, init_posixt,
  init_character, init_factor, init_list,
  init_data_frame, init_unknown
};

// Dispatcher functions
inline SEXP cast_(cheapr::r_type cast_type, SEXP x, SEXP y) {
  return CAST_FNS[cast_type](x, y);
}
inline SEXP init_(cheapr::r_type cast_type, R_xlen_t n, bool with_na) {
  return INIT_FNS[cast_type](n, with_na);
}

#endif
