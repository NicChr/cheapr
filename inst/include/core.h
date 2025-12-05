#ifndef CHEAPR_CORE_H
#define CHEAPR_CORE_H

// cheapr Core definitions and templates
// Author: Nick Christofides
// License: MIT

#include <cpp11.hpp>

#ifdef _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#ifdef _OPENMP
#include <omp.h>
#define OMP_NUM_PROCS omp_get_num_procs()
#define OMP_THREAD_LIMIT omp_get_thread_limit()
#define OMP_MAX_THREADS omp_get_max_threads()
#define OMP_PARALLEL _Pragma("omp parallel num_threads(n_cores) ")
#define OMP_FOR_SIMD _Pragma("omp for simd ")
#define OMP_PARALLEL_FOR_SIMD	_Pragma("omp parallel for simd num_threads(n_cores) ")
#else
#define OMP_NUM_PROCS 1
#define OMP_THREAD_LIMIT 1
#define OMP_MAX_THREADS 1
#define OMP_PARALLEL
#define OMP_FOR_SIMD
#define OMP_PARALLEL_FOR_SIMD
#endif

// Kept for dependency reasons

#ifndef VECTOR_PTR_RO
#define VECTOR_PTR_RO(x) ((const SEXP *) DATAPTR_RO(x))
#endif

#ifndef SHIELD
#define SHIELD(x) (Rf_protect(x))
#endif

#ifndef YIELD
#define YIELD(n) (Rf_unprotect(n))
#endif

namespace cheapr {

// type-safe bool type, similar to Rboolean
// I would use cpp11::r_bool but it doesn't seem to work in
// many cases, throwing ambiguous conversion errors

enum class r_boolean : int {
  r_true = 1,
  r_false = 0,
  r_na = INT_MIN
};

inline constexpr r_boolean r_true = r_boolean::r_true;
inline constexpr r_boolean r_false = r_boolean::r_false;
inline constexpr r_boolean r_na = r_boolean::r_na;

// Constants

inline constexpr int INTEGER_MIN = std::numeric_limits<int>::min() + 1;
inline constexpr int INTEGER_MAX = std::numeric_limits<int>::max();

inline constexpr int64_t INTEGER64_MIN = std::numeric_limits<int64_t>::min() + 1;
inline constexpr int64_t INTEGER64_MAX = std::numeric_limits<int64_t>::max();
inline constexpr int64_t NA_INTEGER64 = std::numeric_limits<int64_t>::min();


inline constexpr int CHEAPR_OMP_THRESHOLD = 100000;
inline constexpr SEXPTYPE CHEAPR_INT64SXP = 64;

inline const SEXP r_null = R_NilValue;
inline const Rcomplex NA_COMPLEX = {{NA_REAL, NA_REAL}};

inline const SEXP empty_env = R_EmptyEnv;
inline const SEXP base_env = R_BaseEnv;

// NAs

namespace na {
  inline const r_boolean logical = r_na;
  inline const int integer = NA_INTEGER;
  inline const double numeric = NA_REAL;
  inline const Rcomplex complex = NA_COMPLEX;
  inline const SEXP string = NA_STRING;
  inline const SEXP list = r_null;
}

// Functions

inline const SEXP* LIST_PTR_RO(SEXP x) {
  return static_cast<const SEXP*>(DATAPTR_RO(x));
}
inline int64_t* INTEGER64_PTR(SEXP x) {
  return reinterpret_cast<int64_t*>(REAL(x));
}
inline const int64_t* INTEGER64_PTR_RO(SEXP x) {
  return reinterpret_cast<const int64_t*>(REAL_RO(x));
}
inline r_boolean* BOOLEAN(SEXP x) {
  return reinterpret_cast<r_boolean*>(INTEGER(x));
}
inline const r_boolean* BOOLEAN_RO(SEXP x) {
  return reinterpret_cast<const r_boolean*>(INTEGER_RO(x));
}

inline SEXPTYPE CHEAPR_TYPEOF(SEXP x){
  return Rf_inherits(x, "integer64") ? CHEAPR_INT64SXP : TYPEOF(x);
}

// Check that n = 0 to avoid R CMD warnings
inline void *safe_memmove(void *dst, const void *src, size_t n){
  return n ? memmove(dst, src, n) : dst;
}

template<typename T>
inline constexpr bool between(const T x, const T lo, const T hi) {
  return x >= lo && x <= hi;
}

// Recycle loop indices
template<typename T>
inline constexpr void recycle_index(T& v, const T size) {
  v = (++v == size) ? static_cast<T>(0) : v;
}

// UTF-8 helpers

inline const char* utf8_char(SEXP x){
  return Rf_translateCharUTF8(x);
}

inline SEXP make_utf8_char(const char *x){
  return Rf_mkCharCE(x, CE_UTF8);
}

inline SEXP make_utf8_str(const char *x){
  return Rf_ScalarString(make_utf8_char(x));
}

inline SEXP install_utf8(const char *x){
  return Rf_installChar(make_utf8_char(x));
}

inline const char* char_as_utf8(const char *x){
  return CHAR(make_utf8_char(x));
}

inline bool is_null(SEXP x){
  return x == r_null;
}

inline bool is_altrep(SEXP x){
  return ALTREP(x);
}

inline bool is_object(SEXP x){
  return Rf_isObject(x);
}

inline SEXP get_attrib(SEXP x, SEXP which){
  return Rf_getAttrib(x, which);
}

inline void set_attrib(SEXP x, SEXP which, SEXP value){
  Rf_setAttrib(x, which, value);
}

inline void set_class(SEXP x, SEXP cls){
  Rf_classgets(x, cls);
}


inline SEXP new_vec(SEXPTYPE type, R_xlen_t n){
  return Rf_allocVector(type, n);
}

inline SEXP new_immutable_vec(SEXPTYPE type, R_xlen_t n){
  SEXP out = SHIELD(new_vec(type, n));
  MARK_NOT_MUTABLE(out);
  YIELD(1);
  return out;
}

inline SEXP coerce_vec(SEXP x, SEXPTYPE type){
  return Rf_coerceVector(x, type);
}

inline int df_nrow(SEXP x){
  return Rf_length(get_attrib(x, R_RowNamesSymbol));
}

inline bool is_r_inf(const double x){
  return x == R_PosInf || x == R_NegInf;
}

// C++ templates

template<typename T>
inline bool is_r_na(const T x) {
  Rf_error("Unimplemented `is_r_na` specialisation");
  return false;
}

template<>
inline bool is_r_na<r_boolean>(const r_boolean x){
  return x == r_na;
}

template<>
inline bool is_r_na<Rboolean>(const Rboolean x){
  return x == NA_LOGICAL;
}

template<>
inline bool is_r_na<cpp11::r_bool>(const cpp11::r_bool x){
  return x == cpp11::na<cpp11::r_bool>();
}

template<>
inline bool is_r_na<int>(const int x){
  return x == NA_INTEGER;
}

template<>
inline bool is_r_na<double>(const double x){
  return x != x;
}

template<>
inline bool is_r_na<int64_t>(const int64_t x){
  return x == NA_INTEGER64;
}

template<>
inline bool is_r_na<Rcomplex>(const Rcomplex x){
  return is_r_na<double>(x.r) || is_r_na<double>(x.i);
}

template<>
inline bool is_r_na<Rbyte>(const Rbyte x){
  return false;
}

template<>
inline bool is_r_na<cpp11::r_string>(const cpp11::r_string x){
  return x == NA_STRING;
}

// Works for CHARSXP and NULL
template<>
inline bool is_r_na<SEXP>(const SEXP x){
  return is_null(x) || x == NA_STRING;
}

// NA type

template<typename T>
inline T na_value(const T x) {
  Rf_error("Unimplemented `na_value` specialisation");
}

template<>
inline r_boolean na_value<r_boolean>(const r_boolean x){
  return r_na;
}

template<>
inline cpp11::r_bool na_value<cpp11::r_bool>(const cpp11::r_bool x){
  return cpp11::na<cpp11::r_bool>();
}

template<>
inline int na_value<int>(const int x){
  return NA_INTEGER;
}

template<>
inline double na_value<double>(const double x){
  return NA_REAL;
}

template<>
inline int64_t na_value<int64_t>(const int64_t x){
  return NA_INTEGER64;
}

template<>
inline Rcomplex na_value<Rcomplex>(const Rcomplex x){
  return NA_COMPLEX;
}

template<>
inline Rbyte na_value<Rbyte>(const Rbyte x){
  return 0;
}

template<>
inline cpp11::r_string na_value<cpp11::r_string>(const cpp11::r_string x){
  return cpp11::na<cpp11::r_string>();
}

template<>
inline SEXP na_value<SEXP>(const SEXP x){
  switch (TYPEOF(x)){
  case CHARSXP: {
    return NA_STRING;
  }
  case VECSXP: {
    return r_null;
  }
  default: {
   Rf_error("No `na_value` specialisation for R type %s", Rf_type2char(TYPEOF(x)));
  }
  }
}

// equals template that doesn't support NA values
// use is_r_na template functions
template<typename T>
inline bool eq(const T x, const T y) {
  return x == y;
}
template<>
inline bool eq<Rcomplex>(const Rcomplex x, const Rcomplex y) {
  return eq(x.r, y.r) && eq(x.i, y.i);
}

template<typename T>
inline SEXP as_r_scalar(const T x){
  if constexpr (std::is_integral<T>::value){
    return as_r_scalar<int64_t>(x);
  } else if constexpr (std::is_convertible_v<T, SEXP>){
    return as_r_scalar<SEXP>(x);
  } else {
    Rf_error("Unimplemented scalar constructor");
  }
}
template<>
inline SEXP as_r_scalar<bool>(const bool x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_r_scalar<r_boolean>(const r_boolean x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_r_scalar<Rboolean>(const Rboolean x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_r_scalar<int>(const int x){
  return Rf_ScalarInteger(x);
}
template<>
inline SEXP as_r_scalar<int64_t>(const int64_t x){
  if (is_r_na(x)){
    return Rf_ScalarInteger(NA_INTEGER);
  } else if (between<int64_t>(x, INTEGER_MIN, INTEGER_MAX)){
    return Rf_ScalarInteger(static_cast<int>(x));
  } else {
    return Rf_ScalarReal(static_cast<double>(x));
  }
}
template<>
inline SEXP as_r_scalar<double>(const double x){
  return Rf_ScalarReal(x);
}
template<>
inline SEXP as_r_scalar<Rcomplex>(const Rcomplex x){
  return Rf_ScalarComplex(x);
}
template<>
inline SEXP as_r_scalar<Rbyte>(const Rbyte x){
  return Rf_ScalarRaw(x);
}
template<>
inline SEXP as_r_scalar<const char *>(const char * const x){
  return make_utf8_str(x);
}
template<>
inline SEXP as_r_scalar<std::string>(const std::string x){
  return as_r_scalar<const char *>(x.c_str());
}
template<>
inline SEXP as_r_scalar<cpp11::r_string>(const cpp11::r_string x){
  return Rf_ScalarString(x);
}

// Scalar string
template<>
inline SEXP as_r_scalar<SEXP>(const SEXP x){
  switch (TYPEOF(x)){
  case CHARSXP: {
    return Rf_ScalarString(x);
  }
  default: {
    SEXP out = SHIELD(new_vec(VECSXP, 1));
    SET_VECTOR_ELT(out, 0, x);
    YIELD(1);
    return out;
  }
  }
}

template<typename T>
inline SEXP as_r_vec(const T x){
  if constexpr (std::is_convertible_v<T, SEXP>){
    return as_r_vec<SEXP>(x);
  } else {
    return as_r_scalar(x);
  }
}
template<>
inline SEXP as_r_vec<SEXP>(const SEXP x){
  switch (TYPEOF(x)){
  case CHARSXP: {
    return Rf_ScalarString(x);
  }
  default: {
    return x;
  }
  }
}

// inline cpp11::sexp as_sexp<r_boolean>(SEXP x){
//   cpp11::as_
// }

// Coerce functions that account for NA
template<typename T>
inline r_boolean as_r_boolean(T x){
  return is_r_na(x) ? r_na : static_cast<r_boolean>(x);
}
template<typename T>
inline int as_int(T x){
  return is_r_na(x) ? NA_INTEGER : static_cast<int>(x);
}
template<typename T>
inline int64_t as_int64(T x){
  return is_r_na(x) ? NA_INTEGER64 : static_cast<int64_t>(x);
}
template<typename T>
inline double as_double(T x){
  return is_r_na(x) ? NA_REAL : static_cast<double>(x);
}
template<typename T>
inline Rcomplex as_complex(T x){
  return is_r_na(x) ? NA_COMPLEX : static_cast<Rcomplex>(x);
}
template<typename T>
inline Rbyte as_raw(T x){
  return is_r_na(x) ? static_cast<Rbyte>(0) : static_cast<Rbyte>(x);
}
// As CHARSXP
template<typename T>
inline SEXP as_char(T x){
  if (is_r_na(x)){
   return NA_STRING;
  } else {
    SEXP scalar = SHIELD(as_r_vec(x));
    SEXP str = SHIELD(coerce_vec(scalar, STRSXP));
    SEXP out = STRING_ELT(str, 0);
    YIELD(2);
    return out;
  }
}

// R fns

// Memory address of R obj
inline SEXP address(SEXP x) {
  char buf[1000];
  snprintf(buf, 1000, "%p", static_cast<void*>(x));
  return make_utf8_char(buf);
}

inline SEXP eval(SEXP expr, SEXP env){
  return Rf_eval(expr, env);
}

// Return R function from a specified package
inline SEXP find_pkg_fun(const char *name, const char *pkg, bool all_fns){

  SEXP expr = r_null;

  if (all_fns){
    expr = SHIELD(Rf_lang3(R_TripleColonSymbol, Rf_install(pkg), Rf_install(name)));
  } else {
    expr = SHIELD(Rf_lang3(R_DoubleColonSymbol, Rf_install(pkg), Rf_install(name)));
  }
  SEXP out = SHIELD(eval(expr, base_env));
  YIELD(2);
  return out;
}

inline SEXP r_length_sym = r_null;

inline R_xlen_t vector_length(SEXP x){
  if (!is_object(x) || Rf_isVectorAtomic(x)){
    return Rf_xlength(x);
  } else if (Rf_inherits(x, "data.frame")){
    return df_nrow(x);
    // Is x a list?
  } else if (TYPEOF(x) == VECSXP){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return Rf_length(x) > 0 ? vector_length(VECTOR_ELT(x, 0)) : 0;
    } else if (Rf_inherits(x, "POSIXlt")){
      const SEXP *p_x = LIST_PTR_RO(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
    } else {
      if (is_null(r_length_sym)){
        SHIELD(r_length_sym = install_utf8("length"));
      } else {
        SHIELD(r_length_sym);
      }
      SEXP expr = SHIELD(Rf_lang2(r_length_sym, x));
      SEXP r_len = SHIELD(eval(expr, R_GetCurrentEnv()));
      R_xlen_t out = TYPEOF(r_len) == INTSXP ? INTEGER_ELT(r_len, 0) : REAL_ELT(r_len, 0);
      YIELD(3);
      return out;
    }
    // Catch-all
  } else {
    if (is_null(r_length_sym)){
      SHIELD(r_length_sym = install_utf8("length"));
    } else {
      SHIELD(r_length_sym);
    }
    SEXP expr = SHIELD(Rf_lang2(r_length_sym, x));
    SEXP r_len = SHIELD(eval(expr, R_GetCurrentEnv()));
    R_xlen_t out = TYPEOF(r_len) == INTSXP ? INTEGER_ELT(r_len, 0) : REAL_ELT(r_len, 0);
    YIELD(3);
    return out;
  }
}

// Can't use `Rf_namesgets(x, r_null)`
// as it adds empty names instead of NULL
inline void set_names(SEXP x, SEXP names){
  is_null(names) ? set_attrib(x, R_NamesSymbol, r_null) : static_cast<void>(Rf_namesgets(x, names));
}
inline SEXP get_names(SEXP x){
  return get_attrib(x, R_NamesSymbol);
}
inline bool has_names(SEXP x){
  SEXP names = SHIELD(get_names(x));
  bool out = !is_null(names);
  YIELD(1);
  return out;
}

// Attributes of x as a list
inline SEXP attributes(SEXP x){
  SEXP a = ATTRIB(x);
  int n = Rf_length(a);

  SEXP out = SHIELD(new_vec(VECSXP, n));
  SEXP names = SHIELD(new_vec(STRSXP, n));
  SEXP current = a;

  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, CAR(current));
    SEXP nm = PRINTNAME(TAG(current));
    if (!is_null(nm)){
      SET_STRING_ELT(names, i, nm);
    }
    current = CDR(current);
  }
  set_names(out, names);
  YIELD(2);
  return out;
}

inline double r_round(double x){
  return is_r_na(x) ? na_value(x) : x - std::remainder(x, 1.0);
}

inline r_boolean is_whole_number(const double x, const double tolerance){
  return is_r_na(x) || is_r_na(tolerance) ? na::logical : static_cast<r_boolean>(std::fabs(x - std::round(x)) < tolerance);
}

inline SEXP deep_copy(SEXP x){
  return Rf_duplicate(x);
}

inline SEXP shallow_copy(SEXP x){
  return Rf_shallow_duplicate(x);
}

inline double gcd2(double x, double y, double tol, bool na_rm){

  if (!na_rm && ( is_r_na(x) || is_r_na(y))){
    return NA_REAL;
  }
  // GCD(0,0)=0
  if (x == 0.0 && y == 0.0){
    return 0.0;
  }
  // GCD(a,0)=a
  if (x == 0.0){
    return y;
  }
  // GCD(a,0)=a
  if (y == 0.0){
    return x;
  }
  double r;
  // Taken from number theory lecture notes
  while(std::fabs(y) > tol){
    r = std::fmod(x, y);
    x = y;
    y = r;
  }
  return x;
}

template<typename T>
inline void set_val(SEXP x, const R_xlen_t i, T val, T* p_x = nullptr){
  Rf_error("Unimplemented `set_val` specialisation");
}

inline void set_val(SEXP x, const R_xlen_t i, bool val, int* p_x = nullptr){
  if (p_x != nullptr){
    p_x[i] = static_cast<int>(val);
  } else {
    SET_LOGICAL_ELT(x, i, static_cast<int>(val));
  }
}
inline void set_val(SEXP x, const R_xlen_t i, r_boolean val, r_boolean* p_x = nullptr){
  if (p_x != nullptr){
    p_x[i] = val;
  } else {
    SET_LOGICAL_ELT(x, i, static_cast<int>(val));
  }
}
inline void set_val(SEXP x, const R_xlen_t i, cpp11::r_bool val, int* p_x = nullptr){
  if (p_x != nullptr){
    p_x[i] = static_cast<int>(val);
  } else {
    SET_LOGICAL_ELT(x, i, static_cast<int>(val));
  }
}
inline void set_val(SEXP x, const R_xlen_t i, int val, int* p_x = nullptr){
  if (p_x != nullptr){
    p_x[i] = val;
  } else {
    SET_INTEGER_ELT(x, i, val);
  }
}
inline void set_val(SEXP x, const R_xlen_t i, double val, double* p_x = nullptr){
  if (p_x != nullptr){
    p_x[i] = val;
  } else {
    SET_REAL_ELT(x, i, val);
  }
}
inline void set_val(SEXP x, const R_xlen_t i, Rcomplex val, Rcomplex* p_x = nullptr){
  if (p_x != nullptr){
    p_x[i] = val;
  } else {
    SET_COMPLEX_ELT(x, i, val);
  }
}
inline void set_val(SEXP x, const R_xlen_t i, Rbyte val, Rbyte* p_x = nullptr){
  if (p_x != nullptr){
    p_x[i] = val;
  } else {
    SET_RAW_ELT(x, i, val);
  }
}
inline void set_val(SEXP x, const R_xlen_t i, const char* val, const SEXP* p_x = nullptr){
  SET_STRING_ELT(x, i, make_utf8_char(val));
}
inline void set_val(SEXP x, const R_xlen_t i, std::string val, const SEXP* p_x = nullptr){
  SET_STRING_ELT(x, i, make_utf8_char(val.c_str()));
}
inline void set_val(SEXP x, const R_xlen_t i, cpp11::r_string val, const SEXP* p_x = nullptr){
  SET_STRING_ELT(x, i, val);
}
// Never use the pointer here to assign
inline void set_val(SEXP x, const R_xlen_t i, SEXP val, const SEXP *p_x = nullptr){
  switch (TYPEOF(val)){
  case CHARSXP: {
    SET_STRING_ELT(x, i, val);
    break;
  }
  default: {
    SET_VECTOR_ELT(x, i, val);
    break;
  }
  }
}

// Named argument

struct arg {
  const char* name;
  SEXP value;

  explicit arg(const char* n) : name(n), value(r_null) {}
  explicit arg(const char* n, SEXP v) : name(n), value(v) {}

  template<typename T>
  arg operator=(T v) const {
    return arg(name, as_r_vec(v));
  }
};

// Variadic list constructor
template<typename... Args>
inline SEXP new_r_list(Args... args) {
  constexpr int n = sizeof...(args);

  if (n == 0){
    return new_vec(VECSXP, 0);
  } else {
    SEXP out = SHIELD(new_vec(VECSXP, n));

    // Are any args named?
    constexpr bool any_named = (std::is_same_v<std::decay_t<Args>, arg> || ...);

    SEXP nms;

    if (any_named){
      nms = SHIELD(new_vec(STRSXP, n));
    } else {
      nms = SHIELD(r_null);
    }

    int i = 0;
    (([&]() {
      if constexpr (std::is_same_v<std::decay_t<Args>, arg>) {
        SET_VECTOR_ELT(out, i, args.value);
        SET_STRING_ELT(nms, i, make_utf8_char(args.name));
      } else {
        SET_VECTOR_ELT(out, i, as_r_vec(args));
      }
      ++i;
    }()), ...);

    set_names(out, nms);
    YIELD(2);
    return out;
  }
}

template<typename... Args>
inline SEXP new_r_pairlist(Args... args) {
  constexpr int n = sizeof...(args);

  if (n == 0){
    return Rf_allocList(0);
  } else {
    SEXP out = SHIELD(Rf_allocList(n));

    SEXP current = out;

    (([&]() {
      if constexpr (std::is_same_v<std::decay_t<Args>, arg>) {
        SETCAR(current, args.value);
        SET_TAG(current, install_utf8(args.name));
      } else {
        SETCAR(current, as_r_vec(args));
      }
      current = CDR(current);
    }()), ...);

    YIELD(1);
    return out;
  }
}

template<typename... Args>
inline SEXP eval_fun(SEXP r_fn, SEXP envir, Args... args){
  // Expression
  SEXP call = SHIELD(Rf_lcons(r_fn, new_r_pairlist(args...)));
  // Evaluate expression
  SEXP out = SHIELD(eval(call, envir));

  YIELD(2);
  return out;
}

// Compact seq generator as ALTREP, same as `seq_len()`
inline SEXP compact_seq_len(R_xlen_t n){
  if (n < 0){
    Rf_error("`n` must be >= 0");
  }
  if (n == 0){
    return new_vec(INTSXP, 0);
  }
  SEXP colon_fn = SHIELD(find_pkg_fun(":", "base", base_env));
  SEXP out = SHIELD(eval_fun(colon_fn, base_env, 1, n));
  YIELD(2);
  return out;
}

namespace vec {
// int not bool because bool can't be NA
inline r_boolean is_whole_number(SEXP x, double tol_, bool na_rm_){

  R_xlen_t n = Rf_xlength(x);

  // Use int instead of bool as int can hold NA
  r_boolean out = r_true;
  bool any_na = false;

  switch ( CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP:
  case CHEAPR_INT64SXP: {
    break;
  }
  case REALSXP: {
    const double *p_x = REAL_RO(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      if (is_r_na(p_x[i])){
        any_na = true;
        continue;
      }
      out = static_cast<r_boolean>(cheapr::is_whole_number(p_x[i], tol_));
      if (out == r_false){
        break;
      }
    }
    if (!na_rm_ && any_na){
      out = r_na;
    }
    break;
  }
  default: {
    out = r_false;
    break;
  }
  }
  return out;
}
}

// Wrap any callable f, and return a new callable that:
//   - takes (auto&&... args)
//   - calls f(args...) inside cpp11::unwind_protect

// Like cpp11::safe but works also  for variadic fns
template <typename F>
auto r_safe_impl(F f) {
  return [f](auto&&... args)
    -> decltype(f(std::forward<decltype(args)>(args)...)) {

      using result_t = decltype(f(std::forward<decltype(args)>(args)...));

      if constexpr (std::is_void_v<result_t>) {
        cpp11::unwind_protect([&] {
          f(std::forward<decltype(args)>(args)...);
        });
        // no return; result_t is void
      } else {
        return cpp11::unwind_protect([&]() -> result_t {
          return f(std::forward<decltype(args)>(args)...);
        });
      }
    };
}

#define r_safe(F)                                                            \
r_safe_impl(                                                                 \
  [&](auto&&... args)                                                        \
    -> decltype(F(std::forward<decltype(args)>(args)...)) {                  \
      return F(std::forward<decltype(args)>(args)...);                       \
    }                                                                        \
)

} // End of cheapr namespace

#endif
