#ifndef CHEAPR_C_CORE_H
#define CHEAPR_C_CORE_H

// cheapr Core definitions and templates
// Author: Nick Christofides
// License: MIT

#include <cpp11.hpp>
#include <optional>
#include <type_traits>

#ifdef _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#if !defined(OBJSXP) && defined(S4SXP)
#define OBJSXP S4SXP
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

// These should be functions in cheapr namespace but
// rchk produces errors with that method
#ifndef SHIELD
#define SHIELD(x) (Rf_protect(x))
#endif

#ifndef YIELD
#define YIELD(n) (Rf_unprotect(n))
#endif

namespace cheapr {

// Constants

namespace internal {

inline constexpr int CHEAPR_OMP_THRESHOLD = 100000;
inline constexpr SEXPTYPE CHEAPR_INT64SXP = 64;
}

namespace r_limits {
inline constexpr int r_int_min = -std::numeric_limits<int>::max();
inline constexpr int r_int_max = std::numeric_limits<int>::max();
inline constexpr int64_t r_int64_min = -std::numeric_limits<int64_t>::max();
inline constexpr int64_t r_int64_max = std::numeric_limits<int64_t>::max();
inline constexpr double r_pos_inf = std::numeric_limits<double>::infinity();
inline constexpr double r_neg_inf = -std::numeric_limits<double>::infinity();
}

// bool type, similar to Rboolean

enum r_bool_t : int {
  r_true = 1,
  r_false = 0,
  r_na = std::numeric_limits<int>::min()
};

inline const SEXP r_null = R_NilValue;

template <typename T>
inline constexpr bool is_r_integral_v =
std::is_integral_v<T> ||
std::is_same_v<T, r_bool_t> ||
std::is_same_v<T, Rboolean> ||
std::is_same_v<T, cpp11::r_bool>;

template <typename T>
inline constexpr bool is_r_arithmetic_v =
std::is_arithmetic_v<T> ||
is_r_integral_v<T>;

namespace env {
inline const SEXP empty_env = R_EmptyEnv;
inline const SEXP base_env = R_BaseEnv;
}

// NAs

namespace na {
  inline constexpr r_bool_t logical = r_na;
  inline constexpr int integer = std::numeric_limits<int>::min();
  inline constexpr int64_t integer64 = std::numeric_limits<int64_t>::min();
  inline const double numeric = NA_REAL;
  inline const Rcomplex complex = {{NA_REAL, NA_REAL}};
  inline constexpr Rbyte raw = static_cast<Rbyte>(0);
  inline const SEXP string = NA_STRING;
  inline const SEXP list = r_null;
}

// Functions

namespace internal {

inline bool inherits1(SEXP x, const char *r_cls){
  return Rf_inherits(x, r_cls);
}

inline SEXPTYPE CHEAPR_TYPEOF(SEXP x){
  return inherits1(x, "integer64") ? internal::CHEAPR_INT64SXP : TYPEOF(x);
}

inline const SEXP* LIST_PTR_RO(SEXP x) {
  return static_cast<const SEXP*>(DATAPTR_RO(x));
}
inline int64_t* INTEGER64_PTR(SEXP x) {
  return reinterpret_cast<int64_t*>(REAL(x));
}
inline const int64_t* INTEGER64_PTR_RO(SEXP x) {
  return reinterpret_cast<const int64_t*>(REAL_RO(x));
}
// Check that n = 0 to avoid R CMD warnings
inline void *safe_memmove(void *dst, const void *src, size_t n){
  return n ? std::memmove(dst, src, n) : dst;
}

inline SEXP new_vec(SEXPTYPE type, R_xlen_t n){
  return Rf_allocVector(type, n);
}

// UTF-8 helpers

inline const char* utf8_char(SEXP x){
  return Rf_translateCharUTF8(x);
}

inline SEXP make_utf8_charsxp(const char *x){
  return Rf_mkCharCE(x, CE_UTF8);
}

inline SEXP make_utf8_strsxp(const char *x){
  return Rf_ScalarString(make_utf8_charsxp(x));
}

inline const char* char_as_utf8(const char *x){
  return CHAR(make_utf8_charsxp(x));
}
}

inline SEXP make_symbol(const char *x){
  return Rf_installChar(internal::make_utf8_charsxp(x));
}

// Memory address
inline SEXP address(SEXP x) {
  char buf[1000];
  std::snprintf(buf, 1000, "%p", static_cast<void*>(x));
  return internal::make_utf8_charsxp(buf);
}

inline r_bool_t* BOOLEAN(SEXP x) {
  return reinterpret_cast<r_bool_t*>(INTEGER(x));
}
inline const r_bool_t* BOOLEAN_RO(SEXP x) {
  return reinterpret_cast<const r_bool_t*>(INTEGER_RO(x));
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

inline bool is_null(SEXP x){
  return x == r_null;
}

namespace altrep {
inline bool is_altrep(SEXP x){
  return ALTREP(x);
}
}

namespace attr {

inline SEXP get_attr(SEXP x, SEXP sym){
  return Rf_getAttrib(x, sym);
}

inline void set_attr(SEXP x, SEXP sym, SEXP value){
  Rf_setAttrib(x, sym, value);
}

}

namespace vec {

inline bool is_object(SEXP x){
  return Rf_isObject(x);
}

inline bool is_atomic(SEXP x){
  return Rf_isVectorAtomic(x);
}

inline bool is_vec(SEXP x){
  return Rf_isVector(x);
}

inline bool is_bare(SEXP x){
  return Rf_length(ATTRIB(x)) == 0;
}

inline bool is_logical(SEXP x){
  return TYPEOF(x) == LGLSXP;
}
inline bool is_integer(SEXP x){
  return TYPEOF(x) == INTSXP;
}
inline bool is_integer64(SEXP x){
  return internal::inherits1(x, "integer64");
}
inline bool is_double(SEXP x){
  return TYPEOF(x) == REALSXP;
}
inline bool is_character(SEXP x){
  return TYPEOF(x) == STRSXP;
}
inline bool is_list(SEXP x){
  return TYPEOF(x) == VECSXP;
}
inline bool is_complex(SEXP x){
  return TYPEOF(x) == CPLXSXP;
}
inline bool is_raw(SEXP x){
  return TYPEOF(x) == RAWSXP;
}
inline bool is_date(SEXP x){
  return internal::inherits1(x, "Date");
}
inline bool is_datetime(SEXP x){
  return internal::inherits1(x, "POSIXt");
}
inline bool is_factor(SEXP x){
  return is_integer(x) && internal::inherits1(x, "factor");
}
inline bool is_df(SEXP x){
  return internal::inherits1(x, "data.frame");
}

}

namespace internal {

inline SEXP r_length_sym = r_null;

inline void set_r_names(SEXP x, SEXP names){
  is_null(names) ? attr::set_attr(x, R_NamesSymbol, r_null) : static_cast<void>(Rf_namesgets(x, names));
}
inline SEXP get_r_names(SEXP x){
  return attr::get_attr(x, R_NamesSymbol);
}
inline bool has_r_names(SEXP x){
  SEXP names = SHIELD(get_r_names(x));
  bool out = !is_null(names);
  YIELD(1);
  return out;
}

}

namespace attr {

inline void set_class(SEXP x, SEXP cls){
  Rf_classgets(x, cls);
}

inline void clear_attrs(SEXP x){
  SEXP current = ATTRIB(x);
  while (!is_null(current)){
    set_attr(x, TAG(current), r_null);
    current = CDR(current);
  }
}

inline bool inherits(SEXP x, SEXP classes){
  R_xlen_t n = Rf_xlength(classes);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (internal::inherits1(x, CHAR(STRING_ELT(classes, i)))){
      return true;
    }
  }
  return false;
}
}

namespace vec {

inline SEXP coerce_vec(SEXP x, SEXPTYPE type){
  return Rf_coerceVector(x, type);
}

inline SEXP new_logical(R_xlen_t n, std::optional<r_bool_t> default_value = std::nullopt) {
  if (default_value.has_value()) {
    r_bool_t val = *default_value;
    SEXP out = SHIELD(internal::new_vec(LGLSXP, n));
    r_bool_t *p_out = BOOLEAN(out);
    std::fill(p_out, p_out + n, val);
    YIELD(1);
    return out;
  } else {
    return internal::new_vec(LGLSXP, n);
  }
}
inline SEXP new_integer(R_xlen_t n, std::optional<int> default_value = std::nullopt){
  if (default_value.has_value()) {
    int val = *default_value;
    SEXP out = SHIELD(internal::new_vec(INTSXP, n));
    int *p_out = INTEGER(out);
    std::fill(p_out, p_out + n, val);
    YIELD(1);
    return out;
  } else {
    return internal::new_vec(INTSXP, n);
  }
}
inline SEXP new_integer64(R_xlen_t n, std::optional<int64_t> default_value = std::nullopt){
  SEXP out = SHIELD(internal::new_vec(REALSXP, n));
  if (default_value.has_value()) {
    int64_t val = *default_value;
    int64_t *p_out = internal::INTEGER64_PTR(out);
    std::fill(p_out, p_out + n, val);
  }
  attr::set_class(out, SHIELD(internal::make_utf8_strsxp("integer64")));
  YIELD(2);
  return out;
}
inline SEXP new_double(R_xlen_t n, std::optional<double> default_value = std::nullopt){
  if (default_value.has_value()){
    double val = *default_value;
    SEXP out = SHIELD(internal::new_vec(REALSXP, n));
    double *p_out = REAL(out);
    std::fill(p_out, p_out + n, val);
    YIELD(1);
    return out;
  } else {
    return internal::new_vec(REALSXP, n);
  }
}
inline SEXP new_character(R_xlen_t n, std::optional<SEXP> default_value = std::nullopt){
  if (default_value.has_value()){
    SEXP val = *default_value;
    SEXP out = SHIELD(internal::new_vec(STRSXP, n));
    if (val != R_BlankString){
      for (R_xlen_t i = 0; i < n; ++i){
        SET_STRING_ELT(out, i, val);
      }
    }
    YIELD(1);
    return out;
  } else {
    return internal::new_vec(STRSXP, n);
  }
}

inline SEXP new_complex(R_xlen_t n, std::optional<Rcomplex> default_value = std::nullopt){
  if (default_value.has_value()){
    Rcomplex val = *default_value;
    SEXP out = SHIELD(internal::new_vec(CPLXSXP, n));
    Rcomplex *p_out = COMPLEX(out);
    std::fill(p_out, p_out + n, val);
    YIELD(1);
    return out;
  } else {
    return internal::new_vec(CPLXSXP, n);
  }
}
inline SEXP new_raw(R_xlen_t n, std::optional<Rbyte> default_value = std::nullopt){
  if (default_value.has_value()){
    Rbyte val = *default_value;
    SEXP out = SHIELD(internal::new_vec(RAWSXP, n));
    Rbyte *p_out = RAW(out);
    std::fill(p_out, p_out + n, val);
    YIELD(1);
    return out;
  } else {
    return internal::new_vec(RAWSXP, n);
  }
}
inline SEXP new_list(R_xlen_t n, std::optional<SEXP> default_value = std::nullopt){
  if (default_value.has_value()){
    SEXP val = *default_value;
    SEXP out = SHIELD(internal::new_vec(VECSXP, n));
    if (!is_null(val)){
      for (R_xlen_t i = 0; i < n; ++i){
        SET_VECTOR_ELT(out, i, val);
      }
    }
    YIELD(1);
    return out;
  } else {
    return internal::new_vec(VECSXP, n);
  }
}
}

namespace df {

inline bool is_df(SEXP x){
  return vec::is_df(x);
}

inline int nrow(SEXP x){
  return Rf_length(attr::get_attr(x, R_RowNamesSymbol));
}
inline int ncol(SEXP x){
  return Rf_length(x);
}
inline SEXP new_row_names(int n){
  if (n > 0){
    SEXP out = SHIELD(vec::new_integer(2));
    SET_INTEGER_ELT(out, 0, na::integer);
    SET_INTEGER_ELT(out, 1, -n);
    YIELD(1);
    return out;
  } else {
    return vec::new_integer(0);
  }
}
inline void set_row_names(SEXP x, int n){
  SEXP row_names = SHIELD(new_row_names(n));
  attr::set_attr(x, R_RowNamesSymbol, row_names);
  YIELD(1);
}
}

inline constexpr bool is_r_inf(const double x){
  return x == r_limits::r_pos_inf || x == r_limits::r_neg_inf;
}

// C++ templates

template<typename T>
inline bool is_r_na(const T x) {
  return false;
}

template<>
inline constexpr bool is_r_na<r_bool_t>(const r_bool_t x){
  return x == na::logical;
}

template<>
inline constexpr bool is_r_na<Rboolean>(const Rboolean x){
  return static_cast<r_bool_t>(x) == na::logical;
}

template<>
inline bool is_r_na<cpp11::r_bool>(const cpp11::r_bool x){
  return x == na::logical;
}

template<>
inline constexpr bool is_r_na<int>(const int x){
  return x == na::integer;
}

template<>
inline constexpr bool is_r_na<double>(const double x){
  return x != x;
}

template<>
inline constexpr bool is_r_na<int64_t>(const int64_t x){
  return x == na::integer64;
}

template<>
inline constexpr bool is_r_na<Rcomplex>(const Rcomplex x){
  return is_r_na<double>(x.r) || is_r_na<double>(x.i);
}

template<>
inline constexpr bool is_r_na<Rbyte>(const Rbyte x){
  return false;
}

template<>
inline bool is_r_na<cpp11::r_string>(const cpp11::r_string x){
  return x == na::string;
}

// Works for CHARSXP and NULL
template<>
inline bool is_r_na<SEXP>(const SEXP x){
  return is_null(x) || x == na::string;
}

// NA type

template<typename T>
inline T na_value(const T x) {
  Rf_error("Unimplemented `na_value` specialisation");
}

template<>
inline constexpr r_bool_t na_value<r_bool_t>(const r_bool_t x){
  return na::logical;
}

template<>
inline cpp11::r_bool na_value<cpp11::r_bool>(const cpp11::r_bool x){
  return cpp11::na<cpp11::r_bool>();
}

template<>
inline constexpr int na_value<int>(const int x){
  return na::integer;
}

template<>
inline double na_value<double>(const double x){
  return na::numeric;
}

template<>
inline constexpr int64_t na_value<int64_t>(const int64_t x){
  return na::integer64;
}

template<>
inline Rcomplex na_value<Rcomplex>(const Rcomplex x){
  return na::complex;
}

template<>
inline constexpr Rbyte na_value<Rbyte>(const Rbyte x){
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
    return na::string;
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
inline constexpr bool eq(const T x, const T y) {
  return x == y;
}
template<>
inline constexpr bool eq<Rcomplex>(const Rcomplex x, const Rcomplex y) {
  return eq(x.r, y.r) && eq(x.i, y.i);
}


namespace vec {

template<typename T>
inline SEXP as_vec(const T x){
  if constexpr (is_r_integral_v<T>){
    return as_vec<int64_t>(x);
  } else if constexpr (std::is_convertible_v<T, SEXP>){
    return as_vec<SEXP>(x);
  } else {
    Rf_error("Unimplemented scalar constructor");
  }
}
template<>
inline SEXP as_vec<bool>(const bool x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vec<r_bool_t>(const r_bool_t x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vec<Rboolean>(const Rboolean x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vec<int>(const int x){
  return Rf_ScalarInteger(x);
}
template<>
inline SEXP as_vec<int64_t>(const int64_t x){
  if (is_r_na(x)){
    return Rf_ScalarInteger(na::integer);
  } else if (between<int64_t>(x, r_limits::r_int_min, r_limits::r_int_max)){
    return Rf_ScalarInteger(static_cast<int>(x));
  } else {
    return Rf_ScalarReal(static_cast<double>(x));
  }
}
template<>
inline SEXP as_vec<double>(const double x){
  return Rf_ScalarReal(x);
}
template<>
inline SEXP as_vec<Rcomplex>(const Rcomplex x){
  return Rf_ScalarComplex(x);
}
template<>
inline SEXP as_vec<Rbyte>(const Rbyte x){
  return Rf_ScalarRaw(x);
}
template<>
inline SEXP as_vec<const char *>(const char * const x){
  return internal::make_utf8_strsxp(x);
}
template<>
inline SEXP as_vec<std::string>(const std::string x){
  return as_vec<const char *>(x.c_str());
}
template<>
inline SEXP as_vec<cpp11::r_string>(const cpp11::r_string x){
  return Rf_ScalarString(x);
}
template<>
inline SEXP as_vec<cpp11::r_bool>(const cpp11::r_bool x){
  return Rf_ScalarLogical(static_cast<int>(x));
}

// Scalar string
template<>
inline SEXP as_vec<SEXP>(const SEXP x){
  switch (TYPEOF(x)){
  case CHARSXP: {
    return Rf_ScalarString(x);
  }
  case LGLSXP:
  case INTSXP:
  case REALSXP:
  case STRSXP:
  case VECSXP:
  case CPLXSXP:
  case RAWSXP: {
    return x;
  }
  default: {
    return new_list(1, x);
  }
  }
}

}

template<typename T>
inline SEXP as_r_obj(const T x){
  if constexpr (std::is_convertible_v<T, SEXP>){
    switch (TYPEOF(x)){
    case CHARSXP: {
      return Rf_ScalarString(x);
    }
    default: {
      return x;
    }
    }
  } else {
    return vec::as_vec(x);
  }
}

namespace internal {
// Assumes no NAs at all
template<typename T>
inline constexpr bool can_be_int(T x){
  // If x is an int type whose size is <= int OR
  // an arithmetic type (e.g. double)
  if constexpr (std::is_integral_v<T> && sizeof(T) <= sizeof(int)){
    return true;
  } else if constexpr (is_r_arithmetic_v<T>){
    return between<T>(x, r_limits::r_int_min, r_limits::r_int_max);
  } else {
    return false;
  }
}
template<typename T>
inline constexpr bool can_be_int64(T x){
  // If x is an int64 type whose size is <= int64 OR
  // an arithmetic type (e.g. double)
  if constexpr (is_r_integral_v<T> && sizeof(T) <= sizeof(int64_t)){
    return true;
  } else if constexpr (is_r_arithmetic_v<T>){
    return between<T>(x, r_limits::r_int64_min, r_limits::r_int64_max);
  } else {
    return false;
  }
}
}

// Coerce functions that account for NA
template<typename T>
inline constexpr r_bool_t as_bool(T x){
  if constexpr (std::is_same_v<T, int>){
    return static_cast<r_bool_t>(x);
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) ? na::logical : static_cast<r_bool_t>(static_cast<bool>(x));
  } else {
    return na::logical;
  }
}
template<typename T>
inline constexpr int as_int(T x){
  if constexpr (std::is_same_v<T, int>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) || !internal::can_be_int(x) ? na::integer : static_cast<int>(x);
  } else {
    return na::integer;
  }
}
template<typename T>
inline constexpr int64_t as_int64(T x){
  if constexpr (std::is_same_v<T, int64_t>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) || !internal::can_be_int64(x) ? na::integer64 : static_cast<int64_t>(x);
  } else {
    return na::integer64;
  }
}
template<typename T>
inline constexpr double as_double(T x){
  if constexpr (std::is_same_v<T, double>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) ? na::numeric : static_cast<double>(x);
  } else {
    return na::numeric;
  }
}
template<typename T>
inline constexpr Rcomplex as_complex(T x){
  if constexpr (std::is_same_v<T, Rcomplex>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return {{as_double(x), 0.0}};
  } else {
    return na::complex;
  }
}
template<typename T>
inline constexpr Rbyte as_raw(T x){
  if constexpr (std::is_same_v<T, Rbyte>){
    return x;
  } else if constexpr (is_r_integral_v<T> && sizeof(T) <= sizeof(int8_t)){
    return is_r_na(x) || x < 0 ? na::raw : static_cast<Rbyte>(x);
  } else if constexpr (std::is_convertible_v<T, Rbyte>){
    return is_r_na(x) || !between(x, static_cast<T>(0), static_cast<T>(255)) ? na::raw : static_cast<Rbyte>(x);
  } else {
    return na::raw;
  }
}
// As CHARSXP
template<typename T>
inline SEXP as_r_string(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, const char *>){
    return internal::make_utf8_charsxp(x);
  } else {
    if (is_r_na(x)){
      return na::string;
    } else {
      SEXP scalar = SHIELD(vec::as_vec(x));
      SEXP str = SHIELD(vec::coerce_vec(scalar, STRSXP));
      SEXP out = STRING_ELT(str, 0);
      YIELD(2);
      return out;
    }
  }
}
namespace internal {
// R version of static_cast
template<typename T, typename U>
struct r_cast_impl {
  static T cast(U x) {
    Rf_error("Can't `r_cast` this type, use `static_cast`");
  }
};

// Specializations for each target type
template<typename U>
struct r_cast_impl<r_bool_t, U> {
  static r_bool_t cast(U x) {
    return as_bool(x);
  }
};

template<typename U>
struct r_cast_impl<int, U> {
  static int cast(U x) {
    return as_int(x);
  }
};

template<typename U>
struct r_cast_impl<int64_t, U> {
  static int64_t cast(U x) {
    return as_int64(x);
  }
};

template<typename U>
struct r_cast_impl<double, U> {
  static double cast(U x) {
    return as_double(x);
  }
};

template<typename U>
struct r_cast_impl<Rcomplex, U> {
  static Rcomplex cast(U x) {
    return as_complex(x);
  }
};

template<typename U>
struct r_cast_impl<Rbyte, U> {
  static Rbyte cast(U x) {
    return as_raw(x);
  }
};

template<typename U>
struct r_cast_impl<SEXP, U> {
  static SEXP cast(U x) {
    return as_r_obj(x);
  }
};
}

template<typename T, typename U>
inline T r_cast(U x) {
  return internal::r_cast_impl<T, U>::cast(x);
}


// R fns

template<typename T>
inline T r_abs(T x){
  if constexpr (is_r_integral_v<T>){
    return is_r_na(x) ? x : static_cast<T>(std::abs(x));
  } else {
    return static_cast<T>(std::abs(x));
  }
}

inline double r_abs(Rcomplex x){
  if (is_r_na(x)){
    return na::numeric;
  } else {
    return std::sqrt(x.r * x.r + x.i * x.i);
  }
}


inline double r_round(double x){
  return is_r_na(x) ? na_value(x) : x - std::remainder(x, 1.0);
}

inline double abs_diff(const double x, const double y){
  return std::fabs(x - y);
}

inline r_bool_t is_whole_number(const double x, const double tolerance){
  return is_r_na(x) || is_r_na(tolerance) ? na::logical : static_cast<r_bool_t>(abs_diff(x, std::round(x)) < tolerance);
}

template<
  typename T,
  typename = typename std::enable_if<is_r_arithmetic_v<T>>::type
>
inline T gcd2(T x, T y, bool na_rm = true, T tol = std::sqrt(std::numeric_limits<T>::epsilon())){

  if (is_r_na(x) || is_r_na(y)){
   if (na_rm){
     if (is_r_na(x)){
       return r_abs(y);
     } else {
       return r_abs(x);
     }
   } else {
     return na_value(x);
   }
  }

  T ax = std::abs(x);
  T ay = std::abs(y);

  if constexpr (is_r_integral_v<T>){

    // Taken from number theory lecture notes

    // GCD(0,0)=0
    if (ax == 0 && ay == 0){
      return 0;
    }
    // GCD(a,0)=a
    if (ax == 0){
      return ay;
    }
    // GCD(a,0)=a
    if (ay == 0){
      return ax;
    }

    T r;
    while(ay != 0){
      r = ax % ay;
      ax = ay;
      ay = r;
    }
    return ax;
  } else {

    // GCD(0,0)=0
    if (ax <= tol && ay <= tol){
      return 0.0;
    }
    // GCD(a,0)=a
    if (ax <= tol){
      return ay;
    }
    // GCD(a,0)=a
    if (ay <= tol){
      return ax;
    }

    T r;
    while(ay > tol){
      r = std::fmod(ax, ay);
      ax = ay;
      ay = r;
    }
    return ax;
  }
}


// Overloaded lowest-common-multiple fn
template<typename T,
         typename = typename std::enable_if<is_r_arithmetic_v<T>>::type>
inline T lcm2(
    T x, T y, bool na_rm = true, T tol = std::sqrt(std::numeric_limits<T>::epsilon())
){
  if (is_r_na(x) || is_r_na(y)){
    if (na_rm){
      if (is_r_na(x)){
        return y;
      } else {
        return x;
      }
    } else {
      return na_value(x);
    }
  }

  if constexpr (is_r_integral_v<T>){
    if (x == 0 && y == 0){
      return 0;
    }
    T res = std::abs(x) / gcd2(x, y, na_rm);
    if (y != 0 && (std::abs(res) > (std::numeric_limits<T>::max() / std::abs(y)))){
      return na_value(x);
    }
    return res * std::abs(y);
  } else {
    if (std::fabs(x) <= tol && std::fabs(y) <= tol){
      return 0.0;
    }
    return ( std::fabs(x) / gcd2(x, y, na_rm, tol) ) * std::fabs(y);
  }
}

inline SEXP eval(SEXP expr, SEXP env){
  return Rf_eval(expr, env);
}

namespace fn {
// Return R function from a specified package
inline SEXP find_pkg_fun(const char *name, const char *pkg, bool all_fns){

  SEXP expr = r_null;

  if (all_fns){
    expr = SHIELD(Rf_lang3(R_TripleColonSymbol, Rf_install(pkg), Rf_install(name)));
  } else {
    expr = SHIELD(Rf_lang3(R_DoubleColonSymbol, Rf_install(pkg), Rf_install(name)));
  }
  SEXP out = SHIELD(eval(expr, env::base_env));
  YIELD(2);
  return out;
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
    return arg(name, as_r_obj(v));
  }
};

namespace vec {

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
inline void set_val(SEXP x, const R_xlen_t i, r_bool_t val, r_bool_t* p_x = nullptr){
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
  SET_STRING_ELT(x, i, internal::make_utf8_charsxp(val));
}
inline void set_val(SEXP x, const R_xlen_t i, std::string val, const SEXP* p_x = nullptr){
  SET_STRING_ELT(x, i, internal::make_utf8_charsxp(val.c_str()));
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

// Variadic list constructor
template<typename... Args>
inline SEXP make_list(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return new_list(0);
  } else {
    SEXP out = SHIELD(vec::new_list(n));

    // Are any args named?
    constexpr bool any_named = (std::is_same_v<std::decay_t<Args>, arg> || ...);

    SEXP nms = r_null;

    if (any_named){
      nms = vec::new_character(n);
    }
    SHIELD(nms);

    int i = 0;
    (([&]() {
      if constexpr (std::is_same_v<std::decay_t<Args>, arg>) {
        SET_VECTOR_ELT(out, i, args.value);
        SET_STRING_ELT(nms, i, internal::make_utf8_charsxp(args.name));
      } else {
        SET_VECTOR_ELT(out, i, as_r_obj(args));
      }
      ++i;
    }()), ...);

    internal::set_r_names(out, nms);
    YIELD(2);
    return out;
  }
}

}

namespace internal {

template<typename... Args>
inline SEXP make_pairlist(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return Rf_allocList(0);
  } else {
    SEXP out = SHIELD(Rf_allocList(n));

    SEXP current = out;

    (([&]() {
      if constexpr (std::is_same_v<std::decay_t<Args>, arg>) {
        SETCAR(current, args.value);
        SET_TAG(current, make_symbol(args.name));
      } else {
        SETCAR(current, as_r_obj(args));
      }
      current = CDR(current);
    }()), ...);

    YIELD(1);
    return out;
  }
}

}

namespace fn {
template<typename... Args>
inline SEXP eval_fn(SEXP r_fn, SEXP envir, Args... args){
  // Expression
  SEXP call = SHIELD(Rf_lcons(r_fn, internal::make_pairlist(args...)));
  // Evaluate expression
  SEXP out = SHIELD(eval(call, envir));

  YIELD(2);
  return out;
}
}


namespace vec {

inline R_xlen_t length(SEXP x){
  if (!vec::is_object(x) || vec::is_atomic(x)){
    return Rf_xlength(x);
  } else if (internal::inherits1(x, "data.frame")){
    return df::nrow(x);
    // Is x a list?
  } else if (TYPEOF(x) == VECSXP){
    if (internal::inherits1(x, "vctrs_rcrd")){
      return Rf_length(x) > 0 ? vec::length(VECTOR_ELT(x, 0)) : 0;
    } else if (internal::inherits1(x, "POSIXlt")){
      const SEXP *p_x = internal::LIST_PTR_RO(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
    } else {
      if (is_null(internal::r_length_sym)){
        internal::r_length_sym = make_symbol("length");
      }
      SEXP expr = SHIELD(Rf_lang2(internal::r_length_sym, x));
      SEXP r_len = SHIELD(eval(expr, R_GetCurrentEnv()));
      R_xlen_t out = TYPEOF(r_len) == INTSXP ? INTEGER_ELT(r_len, 0) : REAL_ELT(r_len, 0);
      YIELD(2);
      return out;
    }
    // Catch-all
  } else {
    if (is_null(internal::r_length_sym)){
      internal::r_length_sym = make_symbol("length");
    }
    SEXP expr = SHIELD(Rf_lang2(internal::r_length_sym, x));
    SEXP r_len = SHIELD(eval(expr, R_GetCurrentEnv()));
    R_xlen_t out = TYPEOF(r_len) == INTSXP ? INTEGER_ELT(r_len, 0) : REAL_ELT(r_len, 0);
    YIELD(2);
    return out;
  }
}

inline SEXP deep_copy(SEXP x){
  return Rf_duplicate(x);
}

inline SEXP shallow_copy(SEXP x){
  return Rf_shallow_duplicate(x);
}

// Compact seq generator as ALTREP, same as `seq_len()`
inline SEXP compact_seq_len(R_xlen_t n){
  if (n < 0){
    Rf_error("`n` must be >= 0");
  }
  if (n == 0){
    return vec::new_integer(0);
  }
  SEXP colon_fn = SHIELD(fn::find_pkg_fun(":", "base", false));
  SEXP out = SHIELD(fn::eval_fn(colon_fn, env::base_env, 1, n));
  YIELD(2);
  return out;
}

// r_bool_t not bool because bool can't be NA
inline r_bool_t all_whole_numbers(SEXP x, double tol_, bool na_rm_){

  R_xlen_t n = Rf_xlength(x);

  // Use int instead of bool as int can hold NA
  r_bool_t out = r_true;
  bool any_na = false;

  switch ( internal::CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP:
  case internal::CHEAPR_INT64SXP: {
    break;
  }
  case REALSXP: {
    const double *p_x = REAL_RO(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      if (is_r_na(p_x[i])){
        any_na = true;
        continue;
      }
      out = static_cast<r_bool_t>(cheapr::is_whole_number(p_x[i], tol_));
      if (out == r_false){
        break;
      }
    }
    if (!na_rm_ && any_na){
      out = na::logical;
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

namespace internal {
// template<typename... Args>
inline void add_attrs(SEXP x, SEXP attrs) {

  int32_t NP = 0;

  switch (TYPEOF(attrs)){
  case NILSXP: {
    break;
  }
  case VECSXP: {
    SEXP names = SHIELD(internal::get_r_names(attrs)); ++NP;
    if (is_null(names)){
      YIELD(NP);
      Rf_error("attributes must be a named list");
    }
    const SEXP *p_attributes = internal::LIST_PTR_RO(attrs);
    const SEXP *p_names = STRING_PTR_RO(names);

    SEXP attr_nm = r_null;

    for (int i = 0; i < Rf_length(names); ++i){
      if (p_names[i] != R_BlankString){
        attr_nm = make_symbol(CHAR(p_names[i]));
        if (address(x) == address(p_attributes[i])){
          SEXP dup_attr = SHIELD(vec::deep_copy(p_attributes[i])); ++NP;
          attr::set_attr(x, attr_nm, dup_attr);
        } else {
          attr::set_attr(x, attr_nm, p_attributes[i]);
        }
      }
    }
    break;
  }
  case LISTSXP: {
    SEXP addr_x = SHIELD(address(x)); ++NP;

    SEXP current = attrs;

    while (!is_null(current)){
      if (is_null(TAG(current)) || PRINTNAME(TAG(current)) == R_BlankString){
        YIELD(NP);
        Rf_error("Please only supply named attributes in %s", __func__);
      }
      if (addr_x == address(CAR(current))){
        SEXP dup_attr = SHIELD(vec::deep_copy(CAR(current))); ++NP;
        attr::set_attr(x, TAG(current), dup_attr);
      } else {
        attr::set_attr(x, TAG(current), CAR(current));
      }
      // Next node
      current = CDR(current);
    }
    break;
  }
  default: {
    Rf_error("`attrs` must be a named list");
  }
  }
    YIELD(NP);
  }

}

namespace attr {

// Attributes of x as a list
inline SEXP get_attrs(SEXP x){
  SEXP a = ATTRIB(x);
  int n = Rf_length(a);

  SEXP out = SHIELD(internal::new_vec(VECSXP, n));
  SEXP names = SHIELD(internal::new_vec(STRSXP, n));
  SEXP current = a;

  for (int i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, CAR(current));
    if (!is_null(TAG(current))){
      SET_STRING_ELT(names, i, PRINTNAME(TAG(current)));
    }
    current = CDR(current);
  }
  internal::set_r_names(out, names);
  YIELD(2);
  return out;
}

template<typename... Args>
inline void modify_attrs(SEXP x, Args... args) {
  SEXP attrs = SHIELD(internal::make_pairlist(args...));;
  internal::add_attrs(x, attrs);
  YIELD(1);
}

inline void set_attrs(SEXP x, SEXP attrs){
  clear_attrs(x);
  internal::add_attrs(x, attrs);
}

}

namespace env {
inline SEXP get(SEXP sym, SEXP env, bool inherits = true){

  int32_t NP = 0;

  if (TYPEOF(sym) != SYMSXP){
    SHIELD(sym = vec::coerce_vec(sym, SYMSXP)); ++NP;
  }

  if (TYPEOF(env) != ENVSXP){
    Rf_error("second argument to '%s' must be an environment", __func__);
  }

  SEXP val = inherits ? Rf_findVar(sym, env) : Rf_findVarInFrame(env, sym);

  if (val == R_MissingArg){
    YIELD(NP);
    Rf_error("arg `sym` cannot be missing");
  } else if (val == R_UnboundValue){
    YIELD(NP);
    return r_null;
  } else if (TYPEOF(val) == PROMSXP){
    SHIELD(val);
    val = eval(val, env);
    YIELD(1);
  }
  YIELD(NP);
  return val;
}
}


namespace internal {
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
}

#define r_safe(F)                                                                      \
internal::r_safe_impl(                                                                 \
  [&](auto&&... args)                                                                  \
    -> decltype(F(std::forward<decltype(args)>(args)...)) {                            \
      return F(std::forward<decltype(args)>(args)...);                                 \
    }                                                                                  \
)

} // End of cheapr namespace

#endif
