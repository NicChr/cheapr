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

#ifndef UNSAFE_VECTOR_PTR
#define UNSAFE_VECTOR_PTR(x) ((SEXP*) DATAPTR(x))
#endif
#ifndef UNSAFE_STRING_PTR
#define UNSAFE_STRING_PTR(x) ((SEXP*) DATAPTR(x))
#endif
#ifndef VECTOR_PTR_RO
#define VECTOR_PTR_RO(x) ((const SEXP*) DATAPTR_RO(x))
#endif
#ifndef INTEGER64_PTR
#define INTEGER64_PTR(x) ((int64_t*) REAL(x))
#endif
#ifndef INTEGER64_PTR_RO
#define INTEGER64_PTR_RO(x) ((const int64_t*) REAL_RO(x))
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

// Not always guaranteed to be 32-bits

#ifndef INTEGER_MIN
#define INTEGER_MIN std::numeric_limits<int>::min() + 1
#endif

#ifndef INTEGER_MAX
#define INTEGER_MAX std::numeric_limits<int>::max()
#endif

// Guaranteed to be 64-bits

#ifndef INTEGER64_MIN
#define INTEGER64_MIN std::numeric_limits<int64_t>::min() + 1
#endif

#ifndef INTEGER64_MAX
#define INTEGER64_MAX std::numeric_limits<int64_t>::max()
#endif

#ifndef NA_INTEGER64
#define NA_INTEGER64 std::numeric_limits<int64_t>::min()
#endif


// is_na macro functions

#ifndef is_na_lgl
#define is_na_lgl(x) ((bool) (x == NA_LOGICAL))
#endif

#ifndef is_na_int
#define is_na_int(x) ((bool) (x == NA_INTEGER))
#endif

#ifndef is_na_dbl
#define is_na_dbl(x) ((bool) (x != x))
#endif

#ifndef is_na_str
#define is_na_str(x) ((bool) (x == NA_STRING))
#endif

#ifndef is_na_cplx
#define is_na_cplx(x) ((bool) (x.r != x.r) || (x.i != x.i))
#endif

#ifndef is_na_raw
#define is_na_raw(x) ((bool) false)
#endif

#ifndef is_na_int64
#define is_na_int64(x) ((bool) (x == NA_INTEGER64))
#endif


#ifndef CHEAPR_INT_TO_INT64
#define CHEAPR_INT_TO_INT64(x) ((int64_t) (x == NA_INTEGER ? NA_INTEGER64 : x))
#endif
#ifndef CHEAPR_DBL_TO_INT64
#define CHEAPR_DBL_TO_INT64(x) ((int64_t) (x != x ? NA_INTEGER64 : x))
#endif
#ifndef CHEAPR_INT64_TO_INT
#define CHEAPR_INT64_TO_INT(x) ((int) (x == NA_INTEGER64 ? NA_INTEGER : x))
#endif
#ifndef CHEAPR_INT64_TO_DBL
#define CHEAPR_INT64_TO_DBL(x) ((double) (x == NA_INTEGER64 ? NA_REAL : x))
#endif

#ifndef CHEAPR_OMP_THRESHOLD
#define CHEAPR_OMP_THRESHOLD 100000
#endif

#ifndef CHEAPR_INT64SXP
#define CHEAPR_INT64SXP 64
#endif

#ifndef CHEAPR_TYPEOF
#define CHEAPR_TYPEOF(x) ( (SEXPTYPE) (Rf_inherits(x, "integer64") ? CHEAPR_INT64SXP : TYPEOF(x)) )
#endif

#ifndef SHIELD
#define SHIELD(x) (Rf_protect(x))
#endif

#ifndef YIELD
#define YIELD(n) (Rf_unprotect(n))
#endif

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
inline constexpr void recycle(T& v, T size) {
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

inline SEXP CHEAPR_CORES = NULL;

inline int num_cores(){
  if (CHEAPR_CORES == NULL){
    CHEAPR_CORES = install_utf8("cheapr.cores");
  }
  int n_cores = Rf_asInteger(Rf_GetOption1(CHEAPR_CORES));
  return n_cores >= 1 ? n_cores : 1;
}

// Memory address of R obj
inline SEXP r_address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return make_utf8_char(buf);
}

inline bool address_equal(SEXP x, SEXP y){
  return r_address(x) == r_address(y);
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

inline R_xlen_t vec_length(SEXP x){
  if (!Rf_isObject(x) || Rf_isVectorAtomic(x)){
    return Rf_xlength(x);
  } else if (is_df(x)){
    return df_nrow(x);
    // Is x a list?
  } else if (TYPEOF(x) == VECSXP){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return Rf_length(x) > 0 ? vec_length(VECTOR_ELT(x, 0)) : 0;
    } else if (Rf_inherits(x, "POSIXlt")){
      const SEXP *p_x = VECTOR_PTR_RO(x);
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

inline bool is_simple_atomic_vec(SEXP x){
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
inline bool is_simple_vec(SEXP x){
  return (is_simple_atomic_vec(x) || is_bare_list(x));
}

inline bool is_simple_atomic_vec2(SEXP x){
  return is_simple_atomic_vec(x) || is_int64(x);
}

inline bool is_simple_vec2(SEXP x){
  return is_simple_vec(x) || is_int64(x);
}

// Because Rf_ScalarLogical sometimes crashes R?.. Need to look into this
inline SEXP scalar_lgl(bool x){
  SEXP out = SHIELD(Rf_allocVector(LGLSXP, 1));
  LOGICAL(out)[0] = x;
  YIELD(1);
  return out;
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

// Can't use `Rf_namesgets(x, R_NilValue)`
// as it adds empty names instead of NULL
inline void set_names(SEXP x, SEXP names){
  names == R_NilValue ? Rf_setAttrib(x, R_NamesSymbol, R_NilValue) : Rf_namesgets(x, names);
}
inline SEXP get_names(SEXP x){
  return Rf_getAttrib(x, R_NamesSymbol);
}

inline double round_nearest_even(double x){
  return x - std::remainder(x, 1.0);
}

inline bool is_whole_number(double x, double tolerance){
  return (std::fabs(x - std::round(x)) < tolerance);
}

// Variadic function to create R list
template<typename... Args>
inline SEXP make_r_list(Args... args){
  constexpr int n = sizeof...(args);
  SEXP out = SHIELD(new_vec(VECSXP, n));
  int i = 0;
  int dummy[] = {(SET_VECTOR_ELT(out, i++, args), 0)...};
  static_cast<void>(dummy);
  YIELD(1);
  return out;
}

// Make a character vec from const char ptrs
template<typename... Args>
inline SEXP make_r_chars(Args... args){
  constexpr int n = sizeof...(args);
  SEXP out = SHIELD(new_vec(STRSXP, n));
  int i = 0;
  int dummy[] = {(SET_STRING_ELT(out, i++, make_utf8_char(args)), 0)...};
  static_cast<void>(dummy);
  YIELD(1);
  return out;
}

#endif
