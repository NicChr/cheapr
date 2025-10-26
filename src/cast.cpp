#include "cheapr.h"
#include <string.h>
#include <unordered_map>
#include <functional>

static SEXP as_lgl = NULL;
static SEXP as_int = NULL;
static SEXP as_dbl = NULL;
static SEXP as_cplx = NULL;
static SEXP as_raw = NULL;
static SEXP as_char = NULL;
static SEXP as_date = NULL;
static SEXP as_posixct = NULL;
static SEXP as_list = NULL;

std::string combine_types(const std::string &a, const std::string &b) {
  if (a < b) {
    return a + "_" + b;
  } else {
    return b + "_" + a;
  }
}

const std::unordered_map<std::string, std::string> type_pairs = {
  {"NULL_NULL", "NULL"},
  {"integer_integer", "integer"},
  {"numeric_numeric", "numeric"},
  {"logical_logical", "logical"},
  {"complex_complex", "complex"},
  {"raw_raw", "raw"},
  {"list_list", "list"},
  {"character_character", "character"},
  {"Date_Date", "Date"},
  {"POSIXt_POSIXt", "POSIXt"},
  {"factor_factor", "factor"},
  {"integer64_integer64", "integer64"},
  {"data.frame_data.frame", "data.frame"},
  {"NULL_integer", "integer"},
  {"NULL_numeric", "numeric"},
  {"NULL_logical", "logical"},
  {"NULL_complex", "complex"},
  {"NULL_raw", "raw"},
  {"NULL_list", "list"},
  {"NULL_character", "character"},
  {"Date_NULL", "Date"},
  {"NULL_POSIXt", "POSIXt"},
  {"NULL_factor", "factor"},
  {"NULL_integer64", "integer64"},
  {"NULL_data.frame", "data.frame"},
  {"integer_numeric", "numeric"},
  {"integer_logical", "integer"},
  {"complex_integer", "complex"},
  {"integer_raw", "raw"},
  {"integer_list", "list"},
  {"character_integer", "character"},
  {"Date_integer", "Date"},
  {"POSIXt_integer", "POSIXt"},
  {"factor_integer", "factor"},
  {"integer_integer64", "integer64"},
  {"data.frame_integer", "data.frame"},
  {"logical_numeric", "numeric"},
  {"complex_numeric", "complex"},
  {"numeric_raw", "raw"},
  {"list_numeric", "list"},
  {"character_numeric", "character"},
  {"Date_numeric", "Date"},
  {"POSIXt_numeric", "POSIXt"},
  {"factor_numeric", "factor"},
  {"integer64_numeric", "numeric"},
  {"data.frame_numeric", "data.frame"},
  {"complex_logical", "complex"},
  {"logical_raw", "raw"},
  {"list_logical", "list"},
  {"character_logical", "character"},
  {"Date_logical", "Date"},
  {"POSIXt_logical", "POSIXt"},
  {"factor_logical", "factor"},
  {"integer64_logical", "integer64"},
  {"data.frame_logical", "data.frame"},
  {"complex_list", "list"},
  {"character_complex", "character"},
  {"complex_factor", "factor"},
  {"complex_data.frame", "data.frame"},
  {"list_raw", "list"},
  {"character_raw", "character"},
  {"data.frame_raw", "data.frame"},
  {"character_list", "list"},
  {"Date_list", "list"},
  {"POSIXt_list", "list"},
  {"factor_list", "list"},
  {"integer64_list", "list"},
  {"data.frame_list", "data.frame"},
  {"Date_character", "character"},
  {"POSIXt_character", "character"},
  {"character_factor", "factor"},
  {"character_integer64", "character"},
  {"character_data.frame", "data.frame"},
  {"POSIXt_Date", "POSIXt"},
  {"Date_factor", "factor"},
  {"Date_integer64", "Date"},
  {"Date_data.frame", "data.frame"},
  {"POSIXt_factor", "factor"},
  {"POSIXt_integer64", "POSIXt"},
  {"POSIXt_data.frame", "data.frame"},
  {"factor_integer64", "factor"},
  {"data.frame_factor", "data.frame"},
  {"data.frame_integer64", "data.frame"},
};

using namespace cpp11;

// lang2str() and R_data_class() taken directly from R, all thanks go to R team
// R_data_class() has been slightly modified

SEXP lang2str(SEXP obj){
  SEXP symb = CAR(obj);
  static SEXP if_sym = 0, while_sym, for_sym, eq_sym, gets_sym,
    lpar_sym, lbrace_sym, call_sym;
  if(!if_sym) {
    if_sym = Rf_install("if");
    while_sym = Rf_install("while");
    for_sym = Rf_install("for");
    eq_sym = Rf_install("=");
    gets_sym = Rf_install("<-");
    lpar_sym = Rf_install("(");
    lbrace_sym = Rf_install("{");
    call_sym = Rf_install("call");
  }
  if(Rf_isSymbol(symb)) {
    if(symb == if_sym || symb == for_sym || symb == while_sym ||
       symb == lpar_sym || symb == lbrace_sym ||
       symb == eq_sym || symb == gets_sym)
      return PRINTNAME(symb);
  }
  return PRINTNAME(call_sym);
}

// `class()`
SEXP r_data_class(SEXP obj){
  SEXP value, klass = Rf_getAttrib(obj, R_ClassSymbol);
  int n = Rf_length(klass);
  if(n > 0){
    return(klass);
  }
  SEXP dim = Rf_getAttrib(obj, R_DimSymbol);
  int nd = Rf_length(dim);

  if(nd > 0) {
    if(nd == 2) {
      SHIELD(klass = new_vec(STRSXP, 2));
      SET_STRING_ELT(klass, 0, Rf_mkChar("matrix"));
      SET_STRING_ELT(klass, 1, Rf_mkChar("array"));
      YIELD(1);
      return klass;
    }
    else {
      klass = Rf_mkChar("array");
    }
  } else {
    SEXPTYPE t = TYPEOF(obj);
    switch(t) {
    case CLOSXP:
    case SPECIALSXP:
    case BUILTINSXP: {
      klass = Rf_mkChar("function");
      break;
    }
    case REALSXP: {
      klass = Rf_mkChar("numeric");
      break;
    }
    case SYMSXP: {
      klass = Rf_mkChar("name");
      break;
    }
    case LANGSXP: {
      klass = lang2str(obj);
      break;
    }
    case OBJSXP: {
      klass = Rf_mkChar(IS_S4_OBJECT(obj) ? "S4" : "object");
      break;
    }
    default:{
      klass = Rf_type2str(t);
      break;
    }
    }
  }
  SHIELD(klass);
  value = Rf_ScalarString(klass);
  YIELD(1);
  return value;
}

SEXP get_classes(SEXP x){
  return r_data_class(x);
}

const char* get_class(SEXP x){
  SEXP classes = SHIELD(get_classes(x));
  int n = Rf_length(classes);
  const char *out = CHAR(STRING_ELT(classes, n - 1));
  YIELD(1);
  return out;
}

std::string common_type(const std::string &a, const std::string &b) {

  std::string data_pair_type = combine_types(a, b);

  auto it = type_pairs.find(data_pair_type);

  if (it == type_pairs.end()) {
    stop("Can't find suitable cast between <%s> and <%s>", a, b);
  }

  return it->second;
}

struct r_null {};
struct r_logical {};
struct r_integer {};
struct r_integer64 {};
struct r_numeric {};
struct r_complex {};
struct r_raw {};
struct r_character {};
struct r_factor {};
struct r_list {};
struct r_date {};
struct r_POSIXt {};
struct r_data_frame {};

// cast template with specialisations
template<typename T>
SEXP cast(SEXP x, SEXP y) {
  stop("Unimplemented type conversion");
}

template<>
inline SEXP cast<r_null>(SEXP x, SEXP y) {
  return R_NilValue;
}

template<>
inline SEXP cast<r_logical>(SEXP x, SEXP y) {
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
inline SEXP cast<r_integer>(SEXP x, SEXP y) {
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
inline SEXP cast<r_integer64>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "integer64")){
    return x;
  } else {
    return coerce_vector(x, CHEAPR_INT64SXP);
  }
}

template<>
inline SEXP cast<r_numeric>(SEXP x, SEXP y) {
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
inline SEXP cast<r_character>(SEXP x, SEXP y) {

  if (Rf_inherits(x, "character")){
    return x;
  } else if (Rf_isObject(x)){
    as_char = as_char != NULL ? as_char : Rf_install("as.character");
    return Rf_eval(Rf_lang2(as_char, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, STRSXP);
  }
}

template<>
inline SEXP cast<r_complex>(SEXP x, SEXP y) {
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
inline SEXP cast<r_raw>(SEXP x, SEXP y) {
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
inline SEXP cast<r_list>(SEXP x, SEXP y) {
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
inline SEXP cast<r_factor>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "factor")){
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

// TO-DO: FIX THIS METHOD
template<>
inline SEXP cast<r_date>(SEXP x, SEXP y) {
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
inline SEXP cast<r_POSIXt>(SEXP x, SEXP y) {
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
inline SEXP cast<r_data_frame>(SEXP x, SEXP y) {
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


// Wrapper functions for cast fns map
inline SEXP cast_null(SEXP x, SEXP y) { return cast<r_null>(x, y); }
inline SEXP cast_logical(SEXP x, SEXP y) { return cast<r_logical>(x, y); }
inline SEXP cast_integer(SEXP x, SEXP y) { return cast<r_integer>(x, y); }
inline SEXP cast_integer64(SEXP x, SEXP y) { return cast<r_integer64>(x, y); }
inline SEXP cast_numeric(SEXP x, SEXP y) { return cast<r_numeric>(x, y); }
inline SEXP cast_character(SEXP x, SEXP y) { return cast<r_character>(x, y); }
inline SEXP cast_complex(SEXP x, SEXP y) { return cast<r_complex>(x, y); }
inline SEXP cast_raw(SEXP x, SEXP y) { return cast<r_raw>(x, y); }
inline SEXP cast_list(SEXP x, SEXP y) { return cast<r_list>(x, y); }
inline SEXP cast_factor(SEXP x, SEXP y) { return cast<r_factor>(x, y); }
inline SEXP cast_date(SEXP x, SEXP y) { return cast<r_date>(x, y); }
inline SEXP cast_posixt(SEXP x, SEXP y) { return cast<r_POSIXt>(x, y); }
inline SEXP cast_data_frame(SEXP x, SEXP y) { return cast<r_data_frame>(x, y); }

// Cast functions
const std::unordered_map<std::string, std::function<SEXP(SEXP, SEXP)>> cast_fns = {
  {"NULL", cast_null},
  {"logical", cast_logical},
  {"integer", cast_integer},
  {"integer64", cast_integer64},
  {"numeric", cast_numeric},
  {"character", cast_character},
  {"complex", cast_complex},
  {"raw", cast_raw},
  {"list", cast_list},
  {"factor", cast_factor},
  {"Date", cast_date},
  {"POSIXt", cast_posixt},
  {"data.frame", cast_data_frame}
};

// Dispatcher function
inline SEXP cast_(const std::string& cast_type, SEXP x, SEXP y) {
  auto it = cast_fns.find(cast_type);
  if (it != cast_fns.end()) {
    return it->second(x, y);
  } else {
    stop("Unknown cast type: %s", cast_type.c_str());
  }
}

[[cpp11::register]]
SEXP cpp_cast(SEXP x, SEXP y) {
  const char *a = get_class(x);
  const char *b = get_class(y);
  std::string cast_type = common_type(a, b);
  return cast_(cast_type, x, y);
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

// hash common type for switch statement in cpp_cast_all() to work
constexpr unsigned int hash_type(const char *s, int off = 0) {
  return !s[off] ? 5381 : (hash_type(s, off + 1) * 33) ^ s[off];
}

[[cpp11::register]]
SEXP cpp_cast_all(SEXP x){

  if (!Rf_isVectorList(x)){
    Rf_error("`x` must be a list");
  }

  R_xlen_t n = Rf_xlength(x);
  if (n <= 1){
    return x;
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);

  // n is guaranteed to be >=2 because of earlier check

  const char *a = get_class(p_x[0]);
  const char *b = get_class(p_x[1]);

  std::string data_pair_type = combine_types(a, b);
  auto it = type_pairs.find(data_pair_type);

  if (it == type_pairs.end()){
    Rf_error("Can't find suitable cast between <%s> and <%s>", a, b);
  }

  std::string common_type = it->second;

  for (R_xlen_t i = 1; i < n; ++i){

    a = common_type.c_str();
    b = get_class(p_x[i]);

    std::string data_pair_type = combine_types(a, b);
    it = type_pairs.find(data_pair_type);

    if (it == type_pairs.end()){
      Rf_error("Can't find suitable cast between <%s> and <%s>", a, b);
    }
    common_type = it->second;
  }

  SEXP out = SHIELD(new_vec(VECSXP, n));

  SEXP temp;
  PROTECT_INDEX temp_idx;
  R_ProtectWithIndex(temp = R_NilValue, &temp_idx);

#define CAST_LOOP(cast_fn)                                     \
  for (R_xlen_t i = 0; i < (n - 1); ++i){                      \
    R_Reprotect(temp = cast_fn(p_x[i], p_x[i + 1]), temp_idx); \
    SET_VECTOR_ELT(out, i, temp);                              \
  }                                                            \
  SET_VECTOR_ELT(out, n - 1, cast_fn(p_x[n - 1], temp));


  switch (hash_type(common_type.c_str())){
  case hash_type("NULL"): {
    break;
  }
  case hash_type("logical"): {
    CAST_LOOP(cast<r_logical>)
    break;
  }
  case hash_type("integer"): {
    CAST_LOOP(cast<r_integer>)
    break;
  }
  case hash_type("integer64"): {
    CAST_LOOP(cast<r_integer64>)
    break;
  }
  case hash_type("numeric"): {
    CAST_LOOP(cast<r_numeric>)
    break;
  }
  case hash_type("character"): {
    CAST_LOOP(cast<r_character>)
    break;
  }
  case hash_type("complex"): {
    CAST_LOOP(cast<r_complex>)
    break;
  }
  case hash_type("raw"): {
    CAST_LOOP(cast<r_raw>)
    break;
  }
  case hash_type("list"): {
    CAST_LOOP(cast<r_list>)
    break;
  }
  case hash_type("factor"): {
    CAST_LOOP(cast<r_factor>)
    break;
  }
  case hash_type("Date"): {
    CAST_LOOP(cast<r_date>)
    break;
  }
  case hash_type("POSIXt"): {
    CAST_LOOP(cast<r_POSIXt>)
    break;
  }
  case hash_type("data.frame"): {
    CAST_LOOP(cast<r_data_frame>)
    break;
  }
  default: {
    YIELD(2);
    Rf_error("Unimplemented cast type");
  }
  }
  YIELD(2);
  return out;
}
