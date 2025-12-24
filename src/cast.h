#ifndef CHEAPR_CAST_H
#define CHEAPR_CAST_H

#include <cheapr/internal/c_core.h>
#include "types.h"
#include "variadic.h"

// Defining custom R types
// for casting, initialising, combining and assigning

// Symbols for R conversion fns

inline SEXP r_as_lgl = cheapr::r_null;
inline SEXP r_as_int = cheapr::r_null;
inline SEXP r_as_dbl = cheapr::r_null;
inline SEXP r_as_char = cheapr::r_null;
inline SEXP r_as_cplx = cheapr::r_null;
inline SEXP r_as_raw = cheapr::r_null;
inline SEXP r_as_date = cheapr::r_null;
inline SEXP r_as_posixct = cheapr::r_null;
inline SEXP r_as_list = cheapr::r_null;

using cheapr::internal::r_type;

// Custom r types

struct r_null_t {};
struct r_logicals_t {};
struct r_integers_t {};
struct r_integers64_t {};
struct r_doubles_t {};
struct r_complexes_t {};
struct r_raws_t {};
struct r_dates_t {};
struct r_posixts_t {};
struct r_characters_t {};
struct r_factors_t {};
struct r_list_t {};
struct r_data_frame_t {};
struct r_unknown_t {};

// r type constants
enum : r_type {
  R_null = 0,
    R_lgl = 1,
    R_int = 2,
    R_int64 = 3,
    R_dbl = 4,
    R_cplx = 5,
    R_raw = 6,
    R_date = 7,
    R_pxt = 8,
    R_chr = 9,
    R_fct = 10,
    R_list = 11,
    R_df = 12,
    R_unk = 13
};

inline constexpr int n_types = 14;

// R type chars
inline constexpr const char* r_type_names[14] = {
  "NULL",        // 0
  "logical",     // 1
  "integer",     // 2
  "integer64",   // 3
  "numeric",     // 4
  "complex",     // 5
  "raw",         // 6
  "Date",        // 7
  "POSIXct",     // 8
  "character",  // 9
  "factor",     // 10
  "list",       // 11
  "data.frame", // 12
  "unknown"     // 13
};

// An n x n matrix of r types and their common cast type

inline constexpr r_type r_type_pairs[14][14] = {
  /*                0-NULL  1-LGL   2-INT   3-I64   4-DBL    5-CPLX  6-RAW   7-DATE  8-PXCT  9-CHR   10-FCT  11-LIST 12-DF 13-UNK */
  /* 0 - NULL */  { R_null, R_lgl,  R_int,  R_int64, R_dbl,  R_cplx, R_raw,  R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 1 - LGL  */  { R_lgl,  R_lgl,  R_int,  R_int64, R_dbl,  R_cplx, R_raw,  R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 2 - INT  */  { R_int,  R_int,  R_int,  R_int64, R_dbl,  R_cplx, R_raw,  R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 3 - I64  */  { R_int64,R_int64, R_int64, R_int64, R_dbl, R_cplx, R_raw, R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 4 - DBL  */  { R_dbl,  R_dbl,  R_dbl,  R_dbl,   R_dbl,  R_cplx, R_raw, R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 5 - CPLX */  { R_cplx, R_cplx, R_cplx, R_cplx,  R_cplx, R_cplx, R_raw, R_date,  R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 6 - RAW  */  { R_raw,  R_raw,  R_raw,  R_raw, R_raw,  R_raw,  R_raw, R_unk, R_unk, R_chr,  R_fct, R_list, R_df, R_unk },
  /* 7 - DATE */  { R_date, R_date, R_date, R_date,  R_date, R_date,  R_unk, R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 8 - PXCT */  { R_pxt, R_pxt, R_pxt, R_pxt,  R_pxt, R_pxt, R_unk, R_pxt, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 9 - CHR  */  { R_chr,  R_chr,  R_chr,  R_chr, R_chr,  R_chr,  R_chr,  R_chr,  R_chr, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 10 - FCT  */  { R_fct,  R_fct,  R_fct,  R_fct, R_fct,  R_fct,  R_fct, R_fct,  R_fct, R_fct,  R_fct,  R_list, R_df, R_unk },
  /* 11 - LIST */  { R_list, R_list, R_list, R_list,  R_list, R_list, R_list, R_list, R_list, R_list, R_list, R_list, R_df, R_unk },
  /* 12 - DF */    { R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_unk },
  /* 13 - Unknown */ { R_unk,  R_unk,  R_unk,  R_unk, R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk }
};

inline r_type common_type(const r_type a, const r_type b) {
  return r_type_pairs[a][b];
}

// Convert single SEXP into r_* code.
inline r_type get_r_type(SEXP x) {

  if (!cheapr::vec::is_object(x)){
    switch (TYPEOF(x)) {
    case NILSXP:  return R_null;
    case LGLSXP:  return R_lgl;
    case INTSXP:  return R_int;
    case REALSXP: return R_dbl;
    case STRSXP:  return R_chr;
    case CPLXSXP: return R_cplx;
    case RAWSXP:  return R_raw;
    case VECSXP:  return R_list;
    default: return R_unk;
    }
  } else {

    if (cheapr::internal::inherits1(x, "factor"))     return R_fct;
    if (cheapr::internal::inherits1(x, "Date"))       return R_date;
    if (cheapr::internal::inherits1(x, "POSIXct"))    return R_pxt;
    if (cheapr::internal::inherits1(x, "data.frame")) return R_df;
    if (cheapr::internal::inherits1(x, "integer64"))  return R_int64;
    return R_unk;
  }
}

inline const char* r_class1(SEXP obj){
  if (cheapr::vec::is_object(obj)){
    return CHAR(STRING_ELT(cheapr::attr::get_old_class(obj), 0));
  } else {
    switch(TYPEOF(obj)) {
    case CLOSXP:
    case SPECIALSXP:
    case BUILTINSXP: {
      return "function";
    }
    case SYMSXP: {
      return "name";
    }
    case OBJSXP: {
      return Rf_isS4(obj) ? "S4" : "object";
    }
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
      return Rf_type2char(TYPEOF(obj));
    }
    }
  }
}

inline const char *r_type_char(SEXP x){

  r_type type = get_r_type(x);

  // If unknown type
  if (type == R_unk){
    return r_class1(x);
  } else {
    return r_type_names[type];
  }
}

// initialise template with specialisations
template<typename T>
inline SEXP init(R_xlen_t n, bool with_na) {
  Rf_error("Unimplemented initialisation");
  return cheapr::r_null;
}

template<>
inline SEXP init<r_null_t>(R_xlen_t n, bool with_na) {
  return cheapr::r_null;
}

template<>
inline SEXP init<r_logicals_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    return cheapr::vec::new_vector<cheapr::r_bool_t>(n, cheapr::na::logical);
  } else {
    return cheapr::vec::new_vector<cheapr::r_bool_t>(n);
  }

}

template<>
inline SEXP init<r_integers_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    return cheapr::vec::new_vector<int>(n, cheapr::na::integer);
  } else {
    return cheapr::vec::new_vector<int>(n);
  }
}

template<>
inline SEXP init<r_integers64_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    return cheapr::vec::new_vector<int64_t>(n, cheapr::na::integer64);
  } else {
    return cheapr::vec::new_vector<int64_t>(n);
  }
}

template<>
inline SEXP init<r_doubles_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    return cheapr::vec::new_vector<double>(n, cheapr::na::real);
  } else {
    return cheapr::vec::new_vector<double>(n);
  }
}

template<>
inline SEXP init<r_characters_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    return cheapr::vec::new_vector<cheapr::r_string_t>(n, cheapr::na::string);
  } else {
    return cheapr::vec::new_vector<cheapr::r_string_t>(n);
  }
}

template<>
inline SEXP init<r_complexes_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    return cheapr::vec::new_vector<cheapr::r_complex_t>(0, cheapr::na::complex);
  } else {
    return cheapr::vec::new_vector<cheapr::r_complex_t>(n);
  }

}

template<>
inline SEXP init<r_raws_t>(R_xlen_t n, bool with_na) {
  if (with_na){
    return cheapr::vec::new_vector<cheapr::r_byte_t>(n, cheapr::r_byte_t{0});
  } else {
    return cheapr::vec::new_vector<cheapr::r_byte_t>(n);
  }
}

template<>
inline SEXP init<r_list_t>(R_xlen_t n, bool with_na) {
  return cheapr::vec::new_list(n);
}

template<>
inline SEXP init<r_factors_t>(R_xlen_t n, bool with_na) {
  SEXP out = SHIELD(init<r_integers_t>(n, with_na));
  SEXP lvls = SHIELD(cheapr::vec::new_vector<cheapr::r_string_t>(0));
  SEXP cls = SHIELD(cheapr::vec::as_vector("factor"));
  cheapr::attr::set_attr(out, cheapr::symbol::levels_sym, lvls);
  cheapr::attr::set_old_class(out, cls);
  YIELD(3);
  return out;
}

template<>
inline SEXP init<r_dates_t>(R_xlen_t n, bool with_na){
  SEXP out = SHIELD(init<r_doubles_t>(n, with_na));
  SEXP cls = SHIELD(cheapr::vec::as_vector("Date"));
  cheapr::attr::set_old_class(out, cls);
  YIELD(2);
  return out;
}

template<>
inline SEXP init<r_posixts_t>(R_xlen_t n, bool with_na) {
  SEXP out = SHIELD(init<r_doubles_t>(n, with_na));
  SEXP tz = SHIELD(cheapr::vec::new_vector<cheapr::r_string_t>(1));
  SEXP cls = SHIELD(cheapr::vec::combine("POSIXct", "POSIXt"));
  cheapr::attr::set_old_class(out, cls);
  cheapr::attr::set_attr(out, cheapr::r_cast<cheapr::r_symbol_t>("tzone"), tz);
  YIELD(3);
  return out;
}

template<>
inline SEXP init<r_data_frame_t>(R_xlen_t n, bool with_na) {
  SEXP out = SHIELD(cheapr::vec::new_list(0));
  SHIELD(out = cpp_new_df(out, cheapr::r_null, false, false));
  SHIELD(out = cpp_na_init(out, n));
  YIELD(3);
  return out;
}

template<>
inline SEXP init<r_unknown_t>(R_xlen_t n, bool with_na) {
  Rf_error("Don't know how to initialise unknown type");
  return cheapr::r_null;
}

inline void check_casted_length(SEXP x, SEXP out){

  if (cheapr::vec::length(x) != cheapr::vec::length(out)){
    Rf_error(
      "Bad cast from type %s to type %s. \n`vec::length(x)` %lld doesn't match returned length of %lld",
      r_type_char(x),
      r_type_char(out),
      static_cast<long long int>(cheapr::vec::length(x)),
      static_cast<long long int>(cheapr::vec::length(out))
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
  return cheapr::r_null;
}

template<>
inline SEXP cast<r_null_t>(SEXP x, SEXP y) {
  return cheapr::r_null;
}

template<>
inline SEXP cast<r_logicals_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "logical")){
    return x;
  } else if (cheapr::vec::is_object(x)){
    r_as_lgl = !cheapr::is_null(r_as_lgl) ? r_as_lgl : cheapr::r_cast<cheapr::r_symbol_t>("as.logical");
    SEXP expr = SHIELD(Rf_lang2(r_as_lgl, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {
    return cheapr::vec::coerce_vec(x, LGLSXP);
  }
}

template<>
inline SEXP cast<r_integers_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "integer")){
    return x;
  } else if (cheapr::vec::is_object(x)){
    r_as_int = !cheapr::is_null(r_as_int) ? r_as_int : cheapr::r_cast<cheapr::r_symbol_t>("as.integer");
    SEXP expr = SHIELD(Rf_lang2(r_as_int, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {
    return cheapr::vec::coerce_vec(x, INTSXP);
  }
}

template<>
inline SEXP cast<r_doubles_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "numeric")){
    return x;
  } else if (cheapr::vec::is_object(x)){
    r_as_dbl = !cheapr::is_null(r_as_dbl) ? r_as_dbl : cheapr::r_cast<cheapr::r_symbol_t>("as.double");
    SEXP expr = SHIELD(Rf_lang2(r_as_dbl, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {
    return cheapr::vec::coerce_vec(x, REALSXP);
  }
}

template<>
inline SEXP cast<r_integers64_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "integer64")){
    return x;
  } else {
    SEXP out = SHIELD(cast<r_doubles_t>(x, cheapr::r_null));
    SHIELD(out = cpp_numeric_to_int64(x));
    YIELD(2);
    return out;
  }
}

template<>
inline SEXP cast<r_characters_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "character")){
    return x;
  } else if (cheapr::internal::inherits1(x, "factor")){
    return factor_as_character(x);
  } else if (cheapr::vec::is_object(x)){
    r_as_char = !cheapr::is_null(r_as_char) ? r_as_char : cheapr::r_cast<cheapr::r_symbol_t>("as.character");
    SEXP expr = SHIELD(Rf_lang2(r_as_char, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {
    return cheapr::vec::coerce_vec(x, STRSXP);
  }
}

template<>
inline SEXP cast<r_complexes_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "complex")){
    return x;
  } else if (cheapr::vec::is_object(x)){
    r_as_cplx = !cheapr::is_null(r_as_cplx) ? r_as_cplx : cheapr::r_cast<cheapr::r_symbol_t>("as.complex");
    SEXP expr = SHIELD(Rf_lang2(r_as_cplx, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {
    return cheapr::vec::coerce_vec(x, CPLXSXP);
  }
}

template<>
inline SEXP cast<r_raws_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "raw")){
    return x;
  } else if (cheapr::vec::is_object(x)){
    r_as_raw = !cheapr::is_null(r_as_raw) ? r_as_raw : cheapr::r_cast<cheapr::r_symbol_t>("as.raw");
    SEXP expr = SHIELD(Rf_lang2(r_as_raw, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {
    return cheapr::vec::coerce_vec(x, RAWSXP);
  }
}

template<>
inline SEXP cast<r_list_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "list")){
    return x;
  } else if (cheapr::vec::is_object(x)){
    r_as_list = !cheapr::is_null(r_as_list) ? r_as_list : cheapr::r_cast<cheapr::r_symbol_t>("as.list");
    SEXP expr = SHIELD(Rf_lang2(r_as_list, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {
    return cheapr::vec::coerce_vec(x, VECSXP);
  }
}

template<>
inline SEXP cast<r_factors_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "factor") && !cheapr::internal::inherits1(y, "factor")){
    return x;
  } else if (cheapr::internal::inherits1(y, "factor")){
    SEXP x_lvls = SHIELD(cheapr::attr::get_attr(x, cheapr::symbol::levels_sym));
    SEXP out_lvls = SHIELD(cheapr::attr::get_attr(y, cheapr::symbol::levels_sym));
    if (R_compute_identical(x_lvls, out_lvls, 0)){
      YIELD(2);
      return x;
    }
    SEXP out = SHIELD(cast<r_characters_t>(x, cheapr::r_null));
    SHIELD(out = character_as_factor(out, out_lvls));
    YIELD(4);
    return out;
  } else if (cheapr::is_null(x)){
    return init<r_factors_t>(0, true);
  } else {
    SEXP cheapr_factor = SHIELD(cheapr::fn::find_pkg_fun("factor_", "cheapr", true));
    SEXP expr = SHIELD(Rf_lang2(cheapr_factor, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(3);
    return out;
  }
}

template<>
inline SEXP cast<r_dates_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "Date") && !cheapr::internal::inherits1(y, "Date")){
    return x;
  } else if  (cheapr::internal::inherits1(x, "Date") && cheapr::internal::inherits1(y, "Date")){
    if (TYPEOF(x) == TYPEOF(y)){
      return x;
    } else {
      return cheapr::vec::coerce_vec(x, TYPEOF(y));
    }
  } else if (cheapr::is_null(x) && cheapr::is_null(y)){
    return init<r_dates_t>(cheapr::vec::length(x), true);
  } else if (cheapr::vec::is_object(x)){
    r_as_date = !cheapr::is_null(r_as_date) ? r_as_date : cheapr::r_cast<cheapr::r_symbol_t>("as.Date");
    SEXP expr = SHIELD(Rf_lang2(r_as_date, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  } else {

    int32_t NP = 0;

    SEXP out = SHIELD(cheapr::vec::shallow_copy(x)); ++NP;
    if (TYPEOF(x) != INTSXP){
      SHIELD(out = cheapr::vec::coerce_vec(x, REALSXP)); ++NP;
    }
    SEXP date_cls = SHIELD(cheapr::vec::as_vector("Date")); ++NP;
    cheapr::attr::set_old_class(out, date_cls);
    YIELD(NP);
    return out;
  }
}

template<>
inline SEXP cast<r_posixts_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "POSIXct") && !cheapr::internal::inherits1(y, "POSIXct")){
    return x;
    // Copy timezone information
  } else if (cheapr::internal::inherits1(x, "POSIXct") && cheapr::internal::inherits1(y, "POSIXct")){
    SEXP x_tzone = SHIELD(cheapr::attr::get_attr(x, cheapr::r_cast<cheapr::r_symbol_t>("tzone")));
    SEXP out_tzone = SHIELD(cheapr::attr::get_attr(y, cheapr::r_cast<cheapr::r_symbol_t>("tzone")));

    if (R_compute_identical(x_tzone, out_tzone, 0)){
      YIELD(2);
      return x;
    }
    SEXP out = SHIELD(cheapr::vec::shallow_copy(x));
    cheapr::attr::set_attr(out, cheapr::r_cast<cheapr::r_symbol_t>("tzone"), out_tzone);
    YIELD(3);
    return out;
  } else if (cheapr::is_null(x) && cheapr::is_null(y)){
    return init<r_posixts_t>(0, true);
    // Fast method for converting into a date into a date-time
  } else if (cheapr::internal::inherits1(x, "Date") && cheapr::internal::inherits1(y, "POSIXct")){
    R_xlen_t n = Rf_xlength(x);
    SEXP out = SHIELD(cheapr::vec::new_vector<double>(n));
    SEXP out_class = SHIELD(cheapr::vec::combine("POSIXct", "POSIXt"));
    SEXP out_tzone = SHIELD(cheapr::attr::get_attr(y, cheapr::r_cast<cheapr::r_symbol_t>("tzone")));
    cheapr::attr::set_old_class(out, out_class);
    cheapr::attr::set_attr(out, cheapr::r_cast<cheapr::r_symbol_t>("tzone"), out_tzone);

    double* RESTRICT p_out = cheapr::internal::real_ptr(out);

    if (TYPEOF(x) == INTSXP){
      const int *p_x = cheapr::internal::integer_ptr_ro(x);
      OMP_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = cheapr::r_cast<double>(p_x[i]) * 86400.0;

    } else {
      const double *p_x = cheapr::internal::real_ptr_ro(x);
      OMP_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = p_x[i] * 86400;
    }

    YIELD(3);
    return out;
  } else if (!cheapr::vec::is_object(x)){
    SEXP out = SHIELD(cheapr::vec::coerce_vec(x, REALSXP));
    SEXP out_class = SHIELD(cheapr::vec::combine("POSIXct", "POSIXt"));
    SEXP out_tzone = SHIELD(cheapr::vec::new_vector<cheapr::r_string_t>(1));
    cheapr::attr::set_old_class(out, out_class);
    cheapr::attr::set_attr(out, cheapr::r_cast<cheapr::r_symbol_t>("tzone"), out_tzone);
    SHIELD(out = cast<r_posixts_t>(out, y)); // To set the correct attributes
    YIELD(4);
    return out;
  } else {
    r_as_posixct = !cheapr::is_null(r_as_posixct) ? r_as_posixct : cheapr::r_cast<cheapr::r_symbol_t>("as.POSIXct");
    SEXP expr = SHIELD(Rf_lang2(r_as_posixct, x));
    SEXP out = SHIELD(cheapr::eval(expr, cheapr::env::base_env));
    check_casted_length(x, out);
    SHIELD(out = cast<r_posixts_t>(out, y)); // To set the correct attributes
    YIELD(3);
    return out;
  }
}

template<>
inline SEXP cast<r_data_frame_t>(SEXP x, SEXP y) {
  if (cheapr::internal::inherits1(x, "data.frame")){
    if (!cheapr::internal::inherits1(y, "data.frame")){
      return x;
    } else {
      return rebuild(x, y, static_cast<bool>(MAYBE_REFERENCED(x)));
    }
  }

  int32_t NP = 0;
  SEXP out = SHIELD(cpp_as_df(x)); ++NP;

  if (cheapr::internal::inherits1(y, "data.frame")){
    SHIELD(out = rebuild(out, y, static_cast<bool>(MAYBE_REFERENCED(out)))); ++NP;
  }
  YIELD(NP);
  return out;
}

template<>
inline SEXP cast<r_unknown_t>(SEXP x, SEXP y) {
  if (std::strcmp(r_class1(x), r_class1(y)) == 0){
    return x;
  } else {
    SEXP base_cast_fn = SHIELD(cheapr::fn::find_pkg_fun("base_cast", "cheapr", true));
    SEXP out = SHIELD(cheapr::fn::eval_fn(base_cast_fn, cheapr::env::base_env, x, y));
    check_casted_length(x, out);
    YIELD(2);
    return out;
  }
}


// Wrapper functions for cast fns
using cast_fn = SEXP(*)(SEXP, SEXP);
inline SEXP cast_null(SEXP x, SEXP y) { return cast<r_null_t>(x, y); }
inline SEXP cast_logical(SEXP x, SEXP y) { return cast<r_logicals_t>(x, y); }
inline SEXP cast_integer(SEXP x, SEXP y) { return cast<r_integers_t>(x, y); }
inline SEXP cast_integer64(SEXP x, SEXP y) { return cast<r_integers64_t>(x, y); }
inline SEXP cast_numeric(SEXP x, SEXP y) { return cast<r_doubles_t>(x, y); }
inline SEXP cast_character(SEXP x, SEXP y) { return cast<r_characters_t>(x, y); }
inline SEXP cast_complex(SEXP x, SEXP y) { return cast<r_complexes_t>(x, y); }
inline SEXP cast_raw(SEXP x, SEXP y) { return cast<r_raws_t>(x, y); }
inline SEXP cast_list(SEXP x, SEXP y) { return cast<r_list_t>(x, y); }
inline SEXP cast_factor(SEXP x, SEXP y) { return cast<r_factors_t>(x, y); }
inline SEXP cast_date(SEXP x, SEXP y) { return cast<r_dates_t>(x, y); }
inline SEXP cast_posixt(SEXP x, SEXP y) { return cast<r_posixts_t>(x, y); }
inline SEXP cast_data_frame(SEXP x, SEXP y) { return cast<r_data_frame_t>(x, y); }
inline SEXP cast_unknown(SEXP x, SEXP y) { return cast<r_unknown_t>(x, y); }

// Wrapper functions for init fns
using init_fn = SEXP(*)(R_xlen_t, bool);
inline SEXP init_null(R_xlen_t n, bool with_na) { return init<r_null_t>(n, with_na); }
inline SEXP init_logical(R_xlen_t n, bool with_na) { return init<r_logicals_t>(n, with_na); }
inline SEXP init_integer(R_xlen_t n, bool with_na) { return init<r_integers_t>(n, with_na); }
inline SEXP init_integer64(R_xlen_t n, bool with_na) { return init<r_integers64_t>(n, with_na); }
inline SEXP init_numeric(R_xlen_t n, bool with_na) { return init<r_doubles_t>(n, with_na); }
inline SEXP init_character(R_xlen_t n, bool with_na) { return init<r_characters_t>(n, with_na); }
inline SEXP init_complex(R_xlen_t n, bool with_na) { return init<r_complexes_t>(n, with_na); }
inline SEXP init_raw(R_xlen_t n, bool with_na) { return init<r_raws_t>(n, with_na); }
inline SEXP init_list(R_xlen_t n, bool with_na) { return init<r_list_t>(n, with_na); }
inline SEXP init_factor(R_xlen_t n, bool with_na) { return init<r_factors_t>(n, with_na); }
inline SEXP init_date(R_xlen_t n, bool with_na) { return init<r_dates_t>(n, with_na); }
inline SEXP init_posixt(R_xlen_t n, bool with_na) { return init<r_posixts_t>(n, with_na); }
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
inline SEXP cast_(r_type cast_type, SEXP x, SEXP y) {
  return CAST_FNS[cast_type](x, y);
}
inline SEXP init_(r_type cast_type, R_xlen_t n, bool with_na) {
  return INIT_FNS[cast_type](n, with_na);
}

#endif
