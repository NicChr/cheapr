#ifndef CHEAPR_CAST_H
#define CHEAPR_CAST_H

#include <cheapr/internal/c_core.h>

namespace cheapr {

// Defining custom R types
// for casting, initialising, combining and assigning
// Here R type must be an object that directly interfaces with R
// Things like `NULL`, vectors, factors, data frames, etc
// Which is different from the lower-level RVal concept for C/C++ we have built

namespace internal {
using r_type = uint8_t;
}

using internal::r_type;

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

template <typename T>
inline r_type get_r_type(T x) {
  static_assert(always_false<T>, "Unsupported type for `as_r_string`");
}



template <RVal U>
inline r_type get_r_type<r_vec<r_lgl>>(r_vec<U> x) {

}





// Convert single SEXP into r_* code.
inline r_type get_r_type(SEXP x) {

  if (!vec::is_object(x)){
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

    if (attr::inherits1(x, "factor"))     return R_fct;
    if (attr::inherits1(x, "Date"))       return R_date;
    if (attr::inherits1(x, "POSIXct"))    return R_pxt;
    if (attr::inherits1(x, "data.frame")) return R_df;
    if (attr::inherits1(x, "integer64"))  return R_int64;
    return R_unk;
  }
}

inline const char* r_class1(SEXP obj){
  if (vec::is_object(obj)){
    return CHAR(STRING_ELT(attr::get_old_class(obj), 0));
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
inline SEXP init(r_size_t n, bool with_na) {
  cpp11::stop("Unimplemented initialisation");
  return r_null;
}

template<>
inline SEXP init<r_null_t>(r_size_t n, bool with_na) {
  return r_null;
}

template<>
inline SEXP init<r_logicals_t>(r_size_t n, bool with_na) {
  if (with_na){
    return r_vec<r_lgl>(n, na::logical);
  } else {
    return r_vec<r_lgl>(n);
  }

}

template<>
inline SEXP init<r_integers_t>(r_size_t n, bool with_na) {
  if (with_na){
    return r_vec<r_int>(n, na::integer);
  } else {
    return r_vec<r_int>(n);
  }
}

template<>
inline SEXP init<r_integers64_t>(r_size_t n, bool with_na) {
  if (with_na){
    return r_vec<r_int64>(n, na::integer64);
  } else {
    return r_vec<r_int64>(n);
  }
}

template<>
inline SEXP init<r_doubles_t>(r_size_t n, bool with_na) {
  if (with_na){
    return r_vec<r_dbl>(n, na::real);
  } else {
    return r_vec<r_dbl>(n);
  }
}

template<>
inline SEXP init<r_characters_t>(r_size_t n, bool with_na) {
  if (with_na){
    return r_vec<r_str>(n, na::string);
  } else {
    return r_vec<r_str>(n);
  }
}

template<>
inline SEXP init<r_complexes_t>(r_size_t n, bool with_na) {
  if (with_na){
    return r_vec<r_cplx>(0, na::complex);
  } else {
    return r_vec<r_cplx>(n);
  }

}

template<>
inline SEXP init<r_raws_t>(r_size_t n, bool with_na) {
  if (with_na){
    return r_vec<r_raw>(n, r_raw{0});
  } else {
    return r_vec<r_raw>(n);
  }
}

template<>
inline SEXP init<r_list_t>(r_size_t n, bool with_na) {
  return r_vec<r_sexp>(n);
}

template<>
inline SEXP init<r_factors_t>(r_size_t n, bool with_na) {
  SEXP out = init<r_integers_t>(n, with_na);
  SEXP lvls = r_vec<r_str>(0);
  SEXP cls = as_vector("factor");
  attr::set_attr(out, symbol::levels_sym, lvls);
  attr::set_old_class(out, cls);
  return out;
}

template<>
inline SEXP init<r_dates_t>(r_size_t n, bool with_na){
  SEXP out = init<r_doubles_t>(n, with_na);
  SEXP cls = as_vector("Date");
  attr::set_old_class(out, cls);
  return out;
}

template<>
inline SEXP init<r_posixts_t>(r_size_t n, bool with_na) {
  SEXP out = init<r_doubles_t>(n, with_na);
  SEXP tz = r_vec<r_str>(1);
  SEXP cls = vec::combine("POSIXct", "POSIXt");
  attr::set_old_class(out, cls);
  attr::set_attr(out, as<r_sym>("tzone"), tz);
  return out;
}

template<>
inline SEXP init<r_data_frame_t>(r_size_t n, bool with_na) {
  SEXP out = r_vec<r_sexp>(0);
  out = cpp_new_df(out, r_null, false, false);
  out = cpp_na_init(out, n);
  return out;
}

template<>
inline SEXP init<r_unknown_t>(r_size_t n, bool with_na) {
  cpp11::stop("Don't know how to initialise unknown type");
  return r_null;
}

inline void check_casted_length(SEXP x, SEXP out){

  if (vec::length(x) != vec::length(out)){
    cpp11::stop(
      "Bad cast from type %s to type %s. \n`vec::length(x)` %lld doesn't match returned length of %lld",
      r_type_char(x),
      r_type_char(out),
      static_cast<long long int>(vec::length(x)),
      static_cast<long long int>(vec::length(out))
    );
  }
}

inline void signal_bad_cast(SEXP from, SEXP to){
  cpp11::stop(
    "Don't know how to cast from type %s to type %s",
    r_type_char(from), r_type_char(to)
  );
}

// cast template with specialisations
template<typename T>
inline SEXP cast(SEXP x, SEXP y) {
  cpp11::stop("Unimplemented cast specialisation");
  return r_null;
}

template<>
inline SEXP cast<r_null_t>(SEXP x, SEXP y) {
  return r_null;
}

template<>
inline SEXP cast<r_logicals_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "logical")){
    return x;
  } else if (vec::is_object(x)){
    r_as_lgl = r_as_lgl != NULL ? r_as_lgl : as<r_sym>("as.logical");
    SEXP expr = Rf_lang2(r_as_lgl, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {
    return internal::coerce_vec(x, LGLSXP);
  }
}

template<>
inline SEXP cast<r_integers_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "integer")){
    return x;
  } else if (vec::is_object(x)){
    r_as_int = r_as_int != NULL ? r_as_int : as<r_sym>("as.integer");
    SEXP expr = Rf_lang2(r_as_int, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {
    return internal::coerce_vec(x, INTSXP);
  }
}

template<>
inline SEXP cast<r_doubles_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "numeric")){
    return x;
  } else if (vec::is_object(x)){
    r_as_dbl = r_as_dbl != NULL  ? r_as_dbl : as<r_sym>("as.double");
    SEXP expr = Rf_lang2(r_as_dbl, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {
    return internal::coerce_vec(x, REALSXP);
  }
}

template<>
inline SEXP cast<r_integers64_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "integer64")){
    return x;
  } else {
    SEXP out = cast<r_doubles_t>(x, r_null);
    out = cpp_numeric_to_int64(x);
    return out;
  }
}

template<>
inline SEXP cast<r_characters_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "character")){
    return x;
  } else if (internal::inherits1(x, "factor")){
    return factor_as_character(x);
  } else if (vec::is_object(x)){
    r_as_char = r_as_char != NULL  ? r_as_char : as<r_sym>("as.character");
    SEXP expr = Rf_lang2(r_as_char, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {
    return internal::coerce_vec(x, STRSXP);
  }
}

template<>
inline SEXP cast<r_complexes_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "complex")){
    return x;
  } else if (vec::is_object(x)){
    r_as_cplx = r_as_cplx != NULL  ? r_as_cplx : as<r_sym>("as.complex");
    SEXP expr = Rf_lang2(r_as_cplx, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {
    return internal::coerce_vec(x, CPLXSXP);
  }
}

template<>
inline SEXP cast<r_raws_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "raw")){
    return x;
  } else if (vec::is_object(x)){
    r_as_raw = r_as_raw != NULL  ? r_as_raw : as<r_sym>("as.raw");
    SEXP expr = Rf_lang2(r_as_raw, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {
    return internal::coerce_vec(x, RAWSXP);
  }
}

template<>
inline SEXP cast<r_list_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "list")){
    return x;
  } else if (vec::is_object(x)){
    r_as_list = r_as_list != NULL  ? r_as_list : as<r_sym>("as.list");
    SEXP expr = Rf_lang2(r_as_list, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {
    return internal::coerce_vec(x, VECSXP);
  }
}

template<>
inline SEXP cast<r_factors_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "factor") && !internal::inherits1(y, "factor")){
    return x;
  } else if (internal::inherits1(y, "factor")){
    SEXP x_lvls = attr::get_attr(x, symbol::levels_sym);
    SEXP out_lvls = attr::get_attr(y, symbol::levels_sym);
    if (R_compute_identical(x_lvls, out_lvls, 0)){
      return x;
    }
    SEXP out = cast<r_characters_t>(x, r_null);
    out = character_as_factor(out, out_lvls);
    return out;
  } else if (x.is_null()){
    return init<r_factors_t>(0, true);
  } else {
    SEXP cheapr_factor = fn::find_pkg_fun("factor_", "cheapr", true);
    SEXP expr = Rf_lang2(cheapr_factor, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  }
}

template<>
inline SEXP cast<r_dates_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "Date") && !internal::inherits1(y, "Date")){
    return x;
  } else if  (internal::inherits1(x, "Date") && internal::inherits1(y, "Date")){
    if (TYPEOF(x) == TYPEOF(y)){
      return x;
    } else {
      return internal::coerce_vec(x, TYPEOF(y));
    }
  } else if (x.is_null() && y.is_null()){
    return init<r_dates_t>(vec::length(x), true);
  } else if (vec::is_object(x)){
    r_as_date = r_as_date != NULL  ? r_as_date : as<r_sym>("as.Date");
    SEXP expr = Rf_lang2(r_as_date, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    return out;
  } else {


    SEXP out = vec::shallow_copy(x);
    if (TYPEOF(x) != INTSXP){
      out = internal::coerce_vec(x, REALSXP);
    }
    SEXP date_cls = as_vector("Date");
    attr::set_old_class(out, date_cls);
    return out;
  }
}

template<>
inline SEXP cast<r_posixts_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "POSIXct") && !internal::inherits1(y, "POSIXct")){
    return x;
    // Copy timezone information
  } else if (internal::inherits1(x, "POSIXct") && internal::inherits1(y, "POSIXct")){
    SEXP x_tzone = attr::get_attr(x, as<r_sym>("tzone"));
    SEXP out_tzone = attr::get_attr(y, as<r_sym>("tzone"));

    if (R_compute_identical(x_tzone, out_tzone, 0)){
      return x;
    }
    SEXP out = vec::shallow_copy(x);
    attr::set_attr(out, as<r_sym>("tzone"), out_tzone);
    return out;
  } else if (x.is_null() && y.is_null()){
    return init<r_posixts_t>(0, true);
    // Fast method for converting into a date into a date-time
  } else if (internal::inherits1(x, "Date") && internal::inherits1(y, "POSIXct")){
    r_size_t n = Rf_xlength(x);
    SEXP out = r_vec<r_dbl>(n);
    SEXP out_class = vec::combine("POSIXct", "POSIXt");
    SEXP out_tzone = attr::get_attr(y, as<r_sym>("tzone"));
    attr::set_old_class(out, out_class);
    attr::set_attr(out, as<r_sym>("tzone"), out_tzone);

    auto* RESTRICT p_out = internal::real_ptr(out);

    if (TYPEOF(x) == INTSXP){
      auto *p_x = vector_ptr<const r_int>(x);
      OMP_SIMD
      for (r_size_t i = 0; i < n; ++i) p_out[i] = as<r_dbl>(p_x[i]) * 86400.0;

    } else {
      auto *p_x = vector_ptr<const r_dbl>(x);
      OMP_SIMD
      for (r_size_t i = 0; i < n; ++i) p_out[i] = p_x[i] * 86400;
    }

    return out;
  } else if (!vec::is_object(x)){
    SEXP out = internal::coerce_vec(x, REALSXP);
    SEXP out_class = vec::combine("POSIXct", "POSIXt");
    SEXP out_tzone = r_vec<r_str>(1);
    attr::set_old_class(out, out_class);
    attr::set_attr(out, as<r_sym>("tzone"), out_tzone);
    out = cast<r_posixts_t>(out, y); // To set the correct attributes
    return out;
  } else {
    r_as_posixct = r_as_posixct != NULL  ? r_as_posixct : as<r_sym>("as.POSIXct");
    SEXP expr = Rf_lang2(r_as_posixct, x);
    SEXP out = eval(expr, env::base_env);
    check_casted_length(x, out);
    out = cast<r_posixts_t>(out, y); // To set the correct attributes
    return out;
  }
}

template<>
inline SEXP cast<r_data_frame_t>(SEXP x, SEXP y) {
  if (internal::inherits1(x, "data.frame")){
    if (!internal::inherits1(y, "data.frame")){
      return x;
    } else {
      return rebuild(x, y, static_cast<bool>(MAYBE_REFERENCED(x)));
    }
  }

  SEXP out = cpp_as_df(x);

  if (internal::inherits1(y, "data.frame")){
    out = rebuild(out, y, static_cast<bool>(MAYBE_REFERENCED(out)));
  }
  return out;
}

template<>
inline SEXP cast<r_unknown_t>(SEXP x, SEXP y) {
  if (std::strcmp(r_class1(x), r_class1(y)) == 0){
    return x;
  } else {
    SEXP base_cast_fn = fn::find_pkg_fun("base_cast", "cheapr", true);
    SEXP out = fn::eval_fn(base_cast_fn, env::base_env, x, y);
    check_casted_length(x, out);
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
using init_fn = SEXP(*)(r_size_t, bool);
inline SEXP init_null(r_size_t n, bool with_na) { return init<r_null_t>(n, with_na); }
inline SEXP init_logical(r_size_t n, bool with_na) { return init<r_logicals_t>(n, with_na); }
inline SEXP init_integer(r_size_t n, bool with_na) { return init<r_integers_t>(n, with_na); }
inline SEXP init_integer64(r_size_t n, bool with_na) { return init<r_integers64_t>(n, with_na); }
inline SEXP init_numeric(r_size_t n, bool with_na) { return init<r_doubles_t>(n, with_na); }
inline SEXP init_character(r_size_t n, bool with_na) { return init<r_characters_t>(n, with_na); }
inline SEXP init_complex(r_size_t n, bool with_na) { return init<r_complexes_t>(n, with_na); }
inline SEXP init_raw(r_size_t n, bool with_na) { return init<r_raws_t>(n, with_na); }
inline SEXP init_list(r_size_t n, bool with_na) { return init<r_list_t>(n, with_na); }
inline SEXP init_factor(r_size_t n, bool with_na) { return init<r_factors_t>(n, with_na); }
inline SEXP init_date(r_size_t n, bool with_na) { return init<r_dates_t>(n, with_na); }
inline SEXP init_posixt(r_size_t n, bool with_na) { return init<r_posixts_t>(n, with_na); }
inline SEXP init_data_frame(r_size_t n, bool with_na) { return init<r_data_frame_t>(n, with_na); }
inline SEXP init_unknown(r_size_t n, bool with_na) { return init<r_unknown_t>(n, with_na); }

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
inline SEXP init_(r_type cast_type, r_size_t n, bool with_na) {
  return INIT_FNS[cast_type](n, with_na);
}

}

#endif
