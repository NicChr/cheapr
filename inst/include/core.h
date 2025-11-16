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

namespace cheapr {

// type-safe bool type, similar to Rboolean

enum class r_bool : int {
  r_true = 1,
  r_false = 0,
  r_na = INT_MIN
};

inline constexpr r_bool r_true = r_bool::r_true;
inline constexpr r_bool r_false = r_bool::r_false;
inline constexpr r_bool r_na = r_bool::r_na;

// Constants

inline constexpr int INTEGER_MIN = std::numeric_limits<int>::min() + 1;
inline constexpr int INTEGER_MAX = std::numeric_limits<int>::max();

inline constexpr int64_t INTEGER64_MIN = std::numeric_limits<int64_t>::min() + 1;
inline constexpr int64_t INTEGER64_MAX = std::numeric_limits<int64_t>::max();
inline constexpr int64_t NA_INTEGER64 = std::numeric_limits<int64_t>::min();


inline constexpr int CHEAPR_OMP_THRESHOLD = 100000;
inline constexpr SEXPTYPE CHEAPR_INT64SXP = 64;

// inline const std::string NA_CPP_STRING = CHAR(NA_STRING);

inline const Rcomplex NA_COMPLEX = {NA_REAL, NA_REAL};

// Functions

inline SEXP* UNSAFE_VECTOR_PTR(SEXP x) {
  return static_cast<SEXP*>(DATAPTR(x));
}
inline SEXP* UNSAFE_STRING_PTR(SEXP x) {
  return static_cast<SEXP*>(DATAPTR(x));
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
inline r_bool* BOOLEAN(SEXP x) {
  return reinterpret_cast<r_bool*>(INTEGER(x));
}
inline const r_bool* BOOLEAN_RO(SEXP x) {
  return reinterpret_cast<const r_bool*>(INTEGER_RO(x));
}

inline SEXP SHIELD(SEXP x){
  return Rf_protect(x);
}

inline void YIELD(int n){
  Rf_unprotect(n);
}

inline SEXPTYPE CHEAPR_TYPEOF(SEXP x){
  return Rf_inherits(x, "integer64") ? CHEAPR_INT64SXP : TYPEOF(x);
}

inline void assert_charsxp(SEXP x){
  if (TYPEOF(x) != CHARSXP){
    Rf_error("`x` must be a CHARSXP, not a %s", Rf_type2char(TYPEOF(x)));
  }
}


// is_na C++ template

template<typename T>
inline bool is_na(T x) {
  Rf_error("Unimplemented `is_na` specialisation");
  return false;
}

template<>
inline bool is_na<r_bool>(r_bool x){
  return x == r_na;
}

template<>
inline bool is_na<int>(int x){
  return x == NA_INTEGER;
}

template<>
inline bool is_na<double>(double x){
  return x != x;
}

template<>
inline bool is_na<int64_t>(int64_t x){
  return x == NA_INTEGER64;
}

template<>
inline bool is_na<Rcomplex>(Rcomplex x){
  return is_na<double>(x.r) || is_na<double>(x.i);
}

template<>
inline bool is_na<Rbyte>(Rbyte x){
  return false;
}

// is_na<SEXP> will always be for CHARSXP
template<>
inline bool is_na<SEXP>(SEXP x){
  return x == NA_STRING;
}

// is_na<SEXP> will always be for CHARSXP
// template<>
// inline bool is_na<std::string>(std::string& x){
//   return x == NA_CPP_STRING;
// }

// NA type

template<typename T>
inline T na_type(T x) {
  Rf_error("Unimplemented `na_type` specialisation");
}

template<>
inline r_bool na_type<r_bool>(r_bool x){
  return r_na;
}

template<>
inline int na_type<int>(int x){
  return NA_INTEGER;
}

template<>
inline double na_type<double>(double x){
  return NA_REAL;
}

template<>
inline int64_t na_type<int64_t>(int64_t x){
  return NA_INTEGER64;
}

template<>
inline Rcomplex na_type<Rcomplex>(Rcomplex x){
  return NA_COMPLEX;
}

template<>
inline SEXP na_type<SEXP>(SEXP x){
  return NA_STRING;
}
// template<>
// inline std::string na_type<std::string>(std::string& x){
//   return NA_CPP_STRING;
// }

// equals template that doesn't support NA values
// use is_na template functions
template<typename T>
inline bool eq(const T& x, const T& y) {
  return x == y;
}
template<>
inline bool eq<Rcomplex>(const Rcomplex& x, const Rcomplex& y) {
  return eq(x.r, y.r) && eq(x.i, y.i);
}

template<typename T>
inline SEXP as_vec_scalar(T x){
  Rf_error("Unimplemented scalar constructor");
}
template<>
inline SEXP as_vec_scalar<bool>(bool x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vec_scalar<r_bool>(r_bool x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vec_scalar<int>(int x){
  return Rf_ScalarInteger(x);
}
template<>
inline SEXP as_vec_scalar<R_xlen_t>(R_xlen_t x){
  if (x <= INTEGER_MAX){
    return Rf_ScalarInteger(static_cast<int>(x));
  } else {
    return Rf_ScalarReal(static_cast<double>(x));
  }
}
// template<>
// inline SEXP as_vec_scalar<int64_t>(int64_t x){
//   SEXP out = SHIELD(Rf_allocVector(REALSXP, 1));
//   Rf_classgets(out, Rf_ScalarString(Rf_mkCharCE("integer64", CE_UTF8)));
//   INTEGER64_PTR(out)[0] = x;
//   YIELD(1);
//   return out;
// }
template<>
inline SEXP as_vec_scalar<double>(double x){
  return Rf_ScalarReal(x);
}
template<>
inline SEXP as_vec_scalar<Rcomplex>(Rcomplex x){
  return Rf_ScalarComplex(x);
}
template<>
inline SEXP as_vec_scalar<Rbyte>(Rbyte x){
  return Rf_ScalarRaw(x);
}
// Scalar string
template<>
inline SEXP as_vec_scalar<SEXP>(SEXP x){
  return Rf_ScalarString(x);
}
// template<>
// inline SEXP as_vec_scalar<std::string>(std::string x){
//   return Rf_mkString(x.c_str());
// }

// Coerce functions that account for NA
template<typename T>
inline r_bool as_r_bool(T x){
  return is_na(x) ? r_na : static_cast<r_bool>(x);
}
template<typename T>
inline int as_int(T x){
  return is_na(x) ? NA_INTEGER : static_cast<int>(x);
}
template<typename T>
inline int64_t as_int64(T x){
  return is_na(x) ? NA_INTEGER64 : static_cast<int64_t>(x);
}
template<typename T>
inline double as_double(T x){
  return is_na(x) ? NA_REAL : static_cast<double>(x);
}
template<typename T>
inline Rcomplex as_complex(T x){
  return is_na(x) ? NA_COMPLEX : static_cast<Rcomplex>(x);
}
template<typename T>
inline Rbyte as_raw(T x){
  return is_na(x) ? static_cast<Rbyte>(0) : static_cast<Rbyte>(x);
}
// As CHARSXP
template<typename T>
inline SEXP as_char(T x){
  if (is_na(x)){
   return NA_STRING;
  } else {
    SEXP scalar = SHIELD(as_vec_scalar(x));
    SEXP str = SHIELD(Rf_coerceVector(scalar, STRSXP));
    SEXP out = STRING_ELT(str, 0);
    YIELD(2);
    return out;
  }
}

// template<typename T>
// inline std::string as_cpp_string(T x){
//   if (is_na(x)){
//     return NA_CPP_STRING;
//   } else {
//     SEXP scalar = SHIELD(as_vec_scalar(x));
//     SEXP str = SHIELD(Rf_coerceVector(scalar, STRSXP));
//     SEXP out = STRING_ELT(str, 0);
//     YIELD(2);
//     return out;
//   }
// }

// R fns

inline cpp11::function cheapr_sset = cpp11::package("cheapr")["cheapr_sset"];
inline cpp11::function cheapr_is_na = cpp11::package("cheapr")["is_na"];
inline cpp11::function cheapr_factor = cpp11::package("cheapr")["factor_"];
inline cpp11::function base_rep = cpp11::package("base")["rep"];
inline cpp11::function cheapr_fast_match = cpp11::package("cheapr")["fast_match"];
inline cpp11::function cheapr_fast_unique = cpp11::package("cheapr")["fast_unique"];
inline cpp11::function cheapr_rebuild = cpp11::package("cheapr")["rebuild"];
inline cpp11::function base_cast = cpp11::package("cheapr")["base_cast"];
inline cpp11::function base_assign = cpp11::package("cheapr")["base_assign_at"];
inline cpp11::function base_length = cpp11::package("base")["length"];

// Check that n = 0 to avoid R CMD warnings
inline void *safe_memmove(void *dst, const void *src, size_t n){
  return n ? memmove(dst, src, n) : dst;
}

template<typename T>
inline constexpr bool between(T x, T lo, T hi) {
  return x >= lo && x <= hi;
}

template<typename T>
inline constexpr bool is_integerable(T x){
  return between<T>(x, INTEGER_MIN, INTEGER_MAX);
}

// Recycle loop indices
template<typename T>
inline constexpr void recycle_index(T& v, T size) {
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
  return Rf_ScalarString(Rf_mkCharCE(x, CE_UTF8));
}

inline SEXP install_utf8(const char *x){
  return Rf_installChar(Rf_mkCharCE(x, CE_UTF8));
}

inline const char* char_as_utf8(const char *x){
  return CHAR(Rf_mkCharCE(x, CE_UTF8));
}

// string paste helper
inline void str_paste(std::string &x, const std::string &sep, const std::string &y){
  x += sep;
  x += y;
}

inline bool is_null(SEXP x){
  return x == R_NilValue;
}

inline bool is_altrep(SEXP x){
  return ALTREP(x);
}

inline SEXP new_vec(SEXPTYPE type, R_xlen_t n){
  return Rf_allocVector(type, n);
}

inline SEXP new_immutable_vec(SEXPTYPE type, R_xlen_t n){
  SEXP out = SHIELD(Rf_allocVector(type, n));
  MARK_NOT_MUTABLE(out);
  YIELD(1);
  return out;
}


inline SEXP coerce_vec(SEXP x, SEXPTYPE type){
  return Rf_coerceVector(x, type);
}

inline bool is_int64(SEXP x){
  return Rf_isReal(x) && Rf_inherits(x, "integer64");
}

inline bool is_df(SEXP x){
  return Rf_inherits(x, "data.frame");
}

inline int df_nrow(SEXP x){
  return Rf_length(Rf_getAttrib(x, R_RowNamesSymbol));
}

inline SEXP CHEAPR_CORES = R_NilValue;

inline int num_cores(){
  if (is_null(CHEAPR_CORES)){
    CHEAPR_CORES = install_utf8("cheapr.cores");
  }
  int n_cores = Rf_asInteger(Rf_GetOption1(CHEAPR_CORES));
  return n_cores >= 1 ? n_cores : 1;
}

// Memory address of R obj
inline SEXP address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return make_utf8_char(buf);
}

inline bool address_equal(SEXP x, SEXP y){
  return address(x) == address(y);
}

// Return R function from a specified package
inline SEXP find_pkg_fun(const char *name, const char *pkg, bool all_fns){

  SEXP expr = R_NilValue;

  if (all_fns){
    expr = SHIELD(Rf_lang3(R_TripleColonSymbol, Rf_install(pkg), Rf_install(name)));
  } else {
    expr = SHIELD(Rf_lang3(R_DoubleColonSymbol, Rf_install(pkg), Rf_install(name)));
  }
  SEXP out = SHIELD(Rf_eval(expr, R_BaseEnv));
  YIELD(2);
  return out;
}

inline R_xlen_t r_length(SEXP x){
  SEXP length_fn = SHIELD(find_pkg_fun("length", "base", false));
  SEXP expr = SHIELD(Rf_lang2(length_fn, x));
  SEXP r_length = SHIELD(Rf_eval(expr, R_BaseEnv));
  R_xlen_t out = TYPEOF(r_length) == INTSXP ? INTEGER(r_length)[0] : REAL(r_length)[0];
  YIELD(3);
  return out;
}

inline R_xlen_t vector_length(SEXP x){
  if (!Rf_isObject(x) || Rf_isVectorAtomic(x)){
    return Rf_xlength(x);
  } else if (is_df(x)){
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
      return r_length(x);
    }
    // Catch-all
  } else {
    return r_length(x);
  }
}

// Definition of simple vector is one in which
// - It is a vector or it is a list with no class
// - Attributes are data-independent
//
// Care must be taken when combining different simple vectors
// as attributes may only be applicable within a single vector
// e.g. for factors (different levels) and POSIXct (different timezones)
//

inline bool cheapr_is_simple_atomic_vec(SEXP x){
  return (
      Rf_isVectorAtomic(x) && (
          !Rf_isObject(x) || (
              Rf_inherits(x, "Date") || Rf_inherits(x, "factor") ||
              Rf_inherits(x, "POSIXct")
          )
      )
  );
}

inline bool is_bare_list(SEXP x){
  return (!Rf_isObject(x) && TYPEOF(x) == VECSXP);
}

inline bool is_bare_atomic(SEXP x){
  return !Rf_isObject(x) && Rf_isVectorAtomic(x);
}

// Sometimes bare lists can be easily handled
inline bool cheapr_is_simple_vec(SEXP x){
  return (cheapr_is_simple_atomic_vec(x) || is_bare_list(x));
}

inline bool cheapr_is_simple_atomic_vec2(SEXP x){
  return cheapr_is_simple_atomic_vec(x) || is_int64(x);
}

inline bool cheapr_is_simple_vec2(SEXP x){
  return cheapr_is_simple_vec(x) || is_int64(x);
}

inline bool is_bare_df(SEXP x){
  SEXP cls = Rf_getAttrib(x, R_ClassSymbol);
  return Rf_length(cls) == 1 &&
    std::strcmp(CHAR(STRING_ELT(cls, 0)), "data.frame") == 0;
}

inline bool is_bare_tbl(SEXP x){
  SEXP xclass = Rf_getAttrib(x, R_ClassSymbol);

  return Rf_length(xclass) == 3 &&
    std::strcmp(CHAR(STRING_ELT(xclass, 0)), "tbl_df") == 0 &&
    std::strcmp(CHAR(STRING_ELT(xclass, 1)), "tbl") == 0 &&
    std::strcmp(CHAR(STRING_ELT(xclass, 2)), "data.frame") == 0;
}

inline const char* r_class(SEXP obj){
  if (Rf_isObject(obj)){
    return CHAR(STRING_ELT(Rf_getAttrib(obj, R_ClassSymbol), 0));
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

// Can't use `Rf_namesgets(x, R_NilValue)`
// as it adds empty names instead of NULL
inline void set_names(SEXP x, SEXP names){
  names == R_NilValue ? Rf_setAttrib(x, R_NamesSymbol, R_NilValue) : Rf_namesgets(x, names);
}
inline SEXP get_names(SEXP x){
  return Rf_getAttrib(x, R_NamesSymbol);
}

template <typename T>
inline int sign(T x) {
  return (T(0) < x) - (x < T(0));
}

template<typename T>
inline T negate(T x){
  return -x;
}

// trunc but eliminates negative zeroes
template<typename T>
inline double trunc(T x){
  return std::trunc(x) + 0.0;
}

inline double round_nearest_even(double x){
  return x - std::remainder(x, 1.0);
}

inline bool is_whole_number(double x, double tolerance){
  return (std::fabs(x - std::round(x)) < tolerance);
}

inline SEXP deep_copy(SEXP x){
  return Rf_duplicate(x);
}

inline SEXP shallow_copy(SEXP x){
  return Rf_shallow_duplicate(x);
}

// int not bool because bool can't be NA
inline r_bool vec_is_whole_number(SEXP x, double tol_, bool na_rm_){

  R_xlen_t n = Rf_xlength(x);

  // Use int instead of bool as int can hold NA
  r_bool out = r_true;
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
      if (is_na(p_x[i])){
        any_na = true;
        continue;
      }
      out = static_cast<r_bool>(is_whole_number(p_x[i], tol_));
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

inline double gcd2(double x, double y, double tol, bool na_rm){

  if (!na_rm && ( is_na(x) || is_na(y))){
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

// Variadic function to create R list
template<typename... Args>
inline SEXP make_r_list(Args... args){
  constexpr int n = sizeof...(args);
  if (n == 0){
    return new_vec(VECSXP, 0);
  } else {
    SEXP out = SHIELD(new_vec(VECSXP, n));
    int i = 0;
    ((SET_VECTOR_ELT(out, i++, args), void()), ...);
    YIELD(1);
    return out;
  }
}

// Make a character vec from const char ptrs
template<typename... Args>
inline SEXP make_r_chars(Args... args){
  constexpr int n = sizeof...(args);
  if (n == 0){
    return new_vec(STRSXP, 0);
  } else {
    SEXP out = SHIELD(new_vec(STRSXP, n));
    int i = 0;
    ((SET_STRING_ELT(out, i++, make_utf8_char(args)), void()), ...);
    YIELD(1);
    return out;
  }
}

} // End of cheapr namespace

#endif
