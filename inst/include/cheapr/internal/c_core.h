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
#define OMP_PRAGMA(x) _Pragma(#x)
#define OMP_NUM_PROCS omp_get_num_procs()
#define OMP_THREAD_LIMIT omp_get_thread_limit()
#define OMP_MAX_THREADS omp_get_max_threads()
#define OMP_PARALLEL(n_threads) OMP_PRAGMA(omp parallel if ((n_threads) > 1) num_threads((n_threads)))
#define OMP_FOR_SIMD OMP_PRAGMA(omp for simd)
#define OMP_SIMD OMP_PRAGMA(omp simd)
#define OMP_PARALLEL_FOR_SIMD(n_threads) OMP_PRAGMA(omp parallel for simd if ((n_threads) > 1) num_threads((n_threads)))

#define OMP_DO_NOTHING
#else
#define OMP_PRAGMA(x)
#define OMP_NUM_PROCS 1
#define OMP_THREAD_LIMIT 1
#define OMP_MAX_THREADS 1
#define OMP_PARALLEL(n_threads)
#define OMP_SIMD
#define OMP_FOR_SIMD
#define OMP_PARALLEL_FOR_SIMD(n_threads)
#define OMP_DO_NOTHING
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

inline constexpr int64_t CHEAPR_OMP_THRESHOLD = 100000;
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

// Alias type for CHARSXP
struct r_string_t {
  SEXP value;
  r_string_t() : value(R_BlankString) {}
  explicit constexpr r_string_t(SEXP x) : value(x) {}
  constexpr operator SEXP() const { return value; }
};

inline const r_string_t blank_r_string = r_string_t();

// Alias type for SYMSXP
struct r_symbol_t {
  SEXP value;
  r_symbol_t() : value(R_MissingArg) {}
  explicit constexpr r_symbol_t(SEXP x) : value(x) {}
  constexpr operator SEXP() const { return value; }
};

namespace symbol {
inline r_symbol_t class_sym = static_cast<r_symbol_t>(R_ClassSymbol);
inline r_symbol_t names_sym = static_cast<r_symbol_t>(R_NamesSymbol);
inline r_symbol_t dim_sym = static_cast<r_symbol_t>(R_DimSymbol);
inline r_symbol_t dim_names_sym = static_cast<r_symbol_t>(R_DimNamesSymbol);
inline r_symbol_t row_names_sym = static_cast<r_symbol_t>(R_RowNamesSymbol);
inline r_symbol_t levels_sym = static_cast<r_symbol_t>(R_LevelsSymbol);
inline r_symbol_t double_colon_sym = static_cast<r_symbol_t>(R_DoubleColonSymbol);
inline r_symbol_t triple_colon_sym = static_cast<r_symbol_t>(R_TripleColonSymbol);
inline r_symbol_t dollar_sym = static_cast<r_symbol_t>(R_DollarSymbol);
inline r_symbol_t bracket_sym = static_cast<r_symbol_t>(R_BracketSymbol);
inline r_symbol_t double_brackets_sym = static_cast<r_symbol_t>(R_Bracket2Symbol);
inline r_symbol_t brace_sym = static_cast<r_symbol_t>(R_BraceSymbol);
inline r_symbol_t dots_sym = static_cast<r_symbol_t>(R_DotsSymbol);
inline r_symbol_t tsp_sym = static_cast<r_symbol_t>(R_TspSymbol);
inline r_symbol_t name_sym = static_cast<r_symbol_t>(R_NameSymbol);
inline r_symbol_t base_sym = static_cast<r_symbol_t>(R_BaseSymbol);
inline r_symbol_t quote_sym = static_cast<r_symbol_t>(R_QuoteSymbol);
inline r_symbol_t function_sym = static_cast<r_symbol_t>(R_FunctionSymbol);
inline r_symbol_t namespace_env_sym = static_cast<r_symbol_t>(R_NamespaceEnvSymbol);
inline r_symbol_t package_sym = static_cast<r_symbol_t>(R_PackageSymbol);
inline r_symbol_t seeds_sym = static_cast<r_symbol_t>(R_SeedsSymbol);
inline r_symbol_t na_rm_sym = static_cast<r_symbol_t>(R_NaRmSymbol);
inline r_symbol_t source_sym = static_cast<r_symbol_t>(R_SourceSymbol);
inline r_symbol_t mode_sym = static_cast<r_symbol_t>(R_ModeSymbol);
inline r_symbol_t device_sym = static_cast<r_symbol_t>(R_DeviceSymbol);
inline r_symbol_t last_value_sym = static_cast<r_symbol_t>(R_LastvalueSymbol);
inline r_symbol_t spec_sym = static_cast<r_symbol_t>(R_SpecSymbol);
inline r_symbol_t previous_sym = static_cast<r_symbol_t>(R_PreviousSymbol);
inline r_symbol_t sort_list_sym = static_cast<r_symbol_t>(R_SortListSymbol);
inline r_symbol_t eval_sym = static_cast<r_symbol_t>(R_EvalSymbol);
inline r_symbol_t drop_sym = static_cast<r_symbol_t>(R_DropSymbol);
inline r_symbol_t missing_arg = static_cast<r_symbol_t>(R_MissingArg);

}

template <typename T>
inline constexpr bool is_r_integral_v =
std::is_integral_v<T> ||
std::is_same_v<std::decay_t<T>, r_bool_t> ||
std::is_same_v<std::decay_t<T>, Rboolean> ||
std::is_same_v<std::decay_t<T>, cpp11::r_bool>;

template <typename T>
inline constexpr bool is_r_arithmetic_v =
is_r_integral_v<T> ||
std::is_arithmetic_v<T>;

namespace env {
inline const SEXP empty_env = R_EmptyEnv;
inline const SEXP base_env = R_BaseEnv;
}

// NAs

namespace na {
  inline constexpr r_bool_t logical = r_na;
  inline constexpr int integer = std::numeric_limits<int>::min();
  inline constexpr int64_t integer64 = std::numeric_limits<int64_t>::min();
  inline const double real = NA_REAL;
  inline const Rcomplex complex = {{NA_REAL, NA_REAL}};
  inline constexpr Rbyte raw = static_cast<Rbyte>(0);
  inline const r_string_t string = static_cast<r_string_t>(NA_STRING);
}

// Pointers

namespace r_ptr {
inline int* integer_ptr(SEXP x){
  return INTEGER(x);
}
inline const int* integer_ptr_ro(SEXP x){
  return INTEGER_RO(x);
}
inline r_bool_t* logical_ptr(SEXP x){
  return reinterpret_cast<r_bool_t*>(integer_ptr(x));
}
inline const r_bool_t* logical_ptr_ro(SEXP x){
  return reinterpret_cast<const r_bool_t*>(integer_ptr_ro(x));
}
inline double* real_ptr(SEXP x){
  return REAL(x);
}
inline const double* real_ptr_ro(SEXP x){
  return REAL_RO(x);
}
inline int64_t* integer64_ptr(SEXP x){
  return reinterpret_cast<int64_t*>(real_ptr(x));
}
inline const int64_t* integer64_ptr_ro(SEXP x){
  return reinterpret_cast<const int64_t*>(real_ptr_ro(x));
}
inline Rcomplex* complex_ptr(SEXP x){
  return COMPLEX(x);
}
inline const Rcomplex* complex_ptr_ro(SEXP x){
  return COMPLEX_RO(x);
}
inline Rbyte* raw_ptr(SEXP x){
  return RAW(x);
}
inline const Rbyte* raw_ptr_ro(SEXP x){
  return RAW_RO(x);
}
inline const r_string_t* string_ptr_ro(SEXP x){
  return reinterpret_cast<const r_string_t*>(STRING_PTR_RO(x));
}
inline const SEXP* list_ptr_ro(SEXP x){
  return VECTOR_PTR_RO(x);
}
}

// Functions

namespace internal {

inline bool inherits1(SEXP x, const char *r_cls){
  return Rf_inherits(x, r_cls);
}

inline SEXPTYPE CHEAPR_TYPEOF(SEXP x){
  return inherits1(x, "integer64") ? internal::CHEAPR_INT64SXP : TYPEOF(x);
}

inline int64_t* INTEGER64_PTR(SEXP x) {
  return r_ptr::integer64_ptr(x);
}
inline const int64_t* INTEGER64_PTR_RO(SEXP x) {
  return r_ptr::integer64_ptr_ro(x);
}
// Check that n = 0 to avoid R CMD warnings
inline void *safe_memmove(void *dst, const void *src, size_t n){
  return n ? std::memmove(dst, src, n) : dst;
}

inline SEXP new_vec(SEXPTYPE type, R_xlen_t n){
  return Rf_allocVector(type, n);
}

// UTF-8 helpers

inline const char* utf8_char(r_string_t x){
  return Rf_translateCharUTF8(static_cast<SEXP>(x));
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

namespace symbol {

inline r_symbol_t tag(SEXP x){
  return static_cast<r_symbol_t>(TAG(x));
}
}

// Memory address
inline r_string_t address(SEXP x) {
  char buf[1000];
  std::snprintf(buf, 1000, "%p", static_cast<void*>(x));
  return static_cast<r_string_t>(internal::make_utf8_charsxp(buf));
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

inline SEXP get_attr(SEXP x, r_symbol_t sym){
  return Rf_getAttrib(x, static_cast<SEXP>(sym));
}

inline void set_attr(SEXP x, r_symbol_t sym, SEXP value){
  Rf_setAttrib(x, static_cast<SEXP>(sym), value);
}

}

namespace internal {
inline SEXP CHEAPR_CORES = r_null;

inline int get_threads(){
  if (is_null(CHEAPR_CORES)){
    CHEAPR_CORES = Rf_installChar(make_utf8_charsxp("cheapr.cores"));
  }
  int n_threads = Rf_asInteger(Rf_GetOption1(CHEAPR_CORES));
  n_threads = std::min(n_threads, OMP_MAX_THREADS);
  return n_threads > 1 ? n_threads : 1;
}

inline int calc_threads(R_xlen_t data_size){
  return data_size >= CHEAPR_OMP_THRESHOLD ? get_threads() : 1;
}

// inline void set_omp_threshold(int64_t n){
//   CHEAPR_OMP_THRESHOLD = n;
// }

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
inline bool is_double(SEXP x){
  return TYPEOF(x) == REALSXP;
}
inline bool is_integer64(SEXP x){
  return is_double(x) && internal::inherits1(x, "integer64");
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

// R vector getters + setters

inline r_bool_t get_r_bool(const r_bool_t* p_x, const R_xlen_t i){
  return static_cast<r_bool_t>(p_x[i]);
}
inline r_bool_t get_r_bool(SEXP x, const R_xlen_t i){
  return static_cast<r_bool_t>(LOGICAL_ELT(x, i));
}
inline int get_int(const int* p_x, const R_xlen_t i){
  return p_x[i];
}
inline int get_int(SEXP x, const R_xlen_t i){
  return INTEGER_ELT(x, i);
}
inline double get_double(const double* p_x, const R_xlen_t i){
  return p_x[i];
}
inline double get_double(SEXP x, const R_xlen_t i){
  return REAL_ELT(x, i);
}
inline double get_int64(const int64_t *p_x, const R_xlen_t i){
  return static_cast<int64_t>(p_x[i]);
}
inline double get_int64(SEXP x, const R_xlen_t i){
  return static_cast<int64_t>(REAL_ELT(x, i));
}
inline Rcomplex get_complex(const Rcomplex* p_x, const R_xlen_t i){
  return p_x[i];
}
inline Rcomplex get_complex(SEXP x, const R_xlen_t i){
  return COMPLEX_ELT(x, i);
}
inline Rbyte get_raw(const Rbyte* p_x, const R_xlen_t i){
  return p_x[i];
}
inline Rbyte get_raw(SEXP x, const R_xlen_t i){
  return RAW_ELT(x, i);
}
inline r_string_t get_r_string(const r_string_t* p_x, const R_xlen_t i){
  return static_cast<r_string_t>(p_x[i]);
}
inline r_string_t get_r_string(SEXP x, const R_xlen_t i){
  return static_cast<r_string_t>(STRING_ELT(x, i));
}
inline r_symbol_t get_r_string(const r_symbol_t* p_x, const R_xlen_t i){
  return static_cast<r_symbol_t>(p_x[i]);
}
inline r_symbol_t get_r_sym(SEXP x, const R_xlen_t i){
  return static_cast<r_symbol_t>(VECTOR_ELT(x, i));
}
inline SEXP get_r_obj(const SEXP* p_x, const R_xlen_t i){
  return p_x[i];
}
inline SEXP get_r_obj(SEXP x, const R_xlen_t i){
  return VECTOR_ELT(x, i);
}

template<typename T>
inline void set_val(SEXP x, const R_xlen_t i, T val){
  static_assert(
    sizeof(T) == 0,
    "Unimplemented `set_val` specialisation"
  );
  return T{};
}

inline void set_val(bool* p_x, const R_xlen_t i, bool val){
  p_x[i] = val;
}
inline void set_val(int* p_x, const R_xlen_t i, bool val){
  p_x[i] = static_cast<int>(val);
}
inline void set_val(SEXP x, const R_xlen_t i, bool val){
  SET_LOGICAL_ELT(x, i, static_cast<int>(val));
}
inline void set_val(r_bool_t* p_x, const R_xlen_t i, r_bool_t val){
  p_x[i] = val;
}
inline void set_val(SEXP x, const R_xlen_t i, r_bool_t val){
  SET_LOGICAL_ELT(x, i, static_cast<int>(val));
}
inline void set_val(int* p_x, const R_xlen_t i, cpp11::r_bool val){
  p_x[i] = static_cast<int>(val);
}
inline void set_val(SEXP x, const R_xlen_t i, cpp11::r_bool val){
  SET_LOGICAL_ELT(x, i, static_cast<int>(val));
}
inline void set_val(int* p_x, const R_xlen_t i, int val){
  p_x[i] = val;
}
inline void set_val(SEXP x, const R_xlen_t i, int val){
  SET_INTEGER_ELT(x, i, val);
}
inline void set_val(int64_t* p_x, const R_xlen_t i, int64_t val){
  p_x[i] = val;
}
inline void set_val(SEXP x, const R_xlen_t i, int64_t val){
  r_ptr::integer64_ptr(x)[i] = val;
}
inline void set_val(double* p_x, const R_xlen_t i, double val){
  p_x[i] = val;
}
inline void set_val(SEXP x, const R_xlen_t i, double val){
  SET_REAL_ELT(x, i, val);
}

inline void set_val(Rcomplex* p_x, const R_xlen_t i, Rcomplex val){
  p_x[i].r = val.r;
  p_x[i].i = val.i;
}
inline void set_val(SEXP x, const R_xlen_t i, Rcomplex val){
  SET_COMPLEX_ELT(x, i, val);
}
inline void set_val(Rbyte* p_x, const R_xlen_t i, Rbyte val){
  p_x[i] = val;
}
inline void set_val(SEXP x, const R_xlen_t i, Rbyte val){
  SET_RAW_ELT(x, i, val);
}
inline void set_val(SEXP x, const R_xlen_t i, const char* val){
  SET_STRING_ELT(x, i, internal::make_utf8_charsxp(val));
}
inline void set_val(SEXP x, const R_xlen_t i, std::string val){
  SET_STRING_ELT(x, i, internal::make_utf8_charsxp(val.c_str()));
}
inline void set_val(SEXP x, const R_xlen_t i, cpp11::r_string val){
  SET_STRING_ELT(x, i, val);
}
inline void set_val(SEXP x, const R_xlen_t i, r_string_t val){
  SET_STRING_ELT(x, i, static_cast<SEXP>(val));
}
inline void set_val(SEXP x, const R_xlen_t i, r_symbol_t val){
  SET_VECTOR_ELT(x, i, static_cast<SEXP>(val));
}
// Never use the pointer here to assign
inline void set_val(SEXP x, const R_xlen_t i, SEXP val){
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

template <typename T>
inline void fast_fill(T *first, T *last, const T val) {

  if constexpr (is_r_arithmetic_v<T> || std::is_same_v<std::decay_t<T>, Rcomplex>){
    R_xlen_t size = last - first;
    int n_threads = internal::calc_threads(size);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < size; ++i) {
        vec::set_val(first, i, val);
      }
    } else {
      std::fill(first, last, val);
    }
  } else {
    for (T *it = first; it != last; ++it) {
      vec::set_val(first, it - first, val);
    }
  }
}

template <typename T>
inline void fast_replace(T *first, T *last, const T old_val, const T new_val) {

  if constexpr (is_r_arithmetic_v<T> || std::is_same_v<std::decay_t<T>, Rcomplex>){
    R_xlen_t size = last - first;
    int n_threads = internal::calc_threads(size);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < size; ++i) {
        if (eq(first[i], old_val)){
          vec::set_val(first, i, new_val);
        }
      }
    } else {
      std::replace(first, last, old_val, new_val);
    }
  } else {
    for (T *it = first; it != last; ++it) {
      R_xlen_t i = it - first;
      if (eq(first[i], old_val)){
        vec::set_val(first, i, new_val);
      }
    }
  }
}
template <typename T>
inline void fast_copy_n(const T *source, R_xlen_t n, T *target){

  if constexpr (is_r_arithmetic_v<T> || std::is_same_v<std::decay_t<T>, Rcomplex>){
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        vec::set_val(target, i, source[i]);
      }
    } else {
      std::copy_n(source, n, target);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      vec::set_val(target, i, source[i]);
    }
  }
}

namespace internal {

inline SEXP r_length_sym = r_null;

inline void set_r_names(SEXP x, SEXP names){
  is_null(names) ? attr::set_attr(x, symbol::names_sym, r_null) : static_cast<void>(Rf_namesgets(x, names));
}
inline SEXP get_r_names(SEXP x){
  return attr::get_attr(x, symbol::names_sym);
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
    set_attr(x, symbol::tag(current), r_null);
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

template <typename T>
inline SEXP new_vector(R_xlen_t n, const T default_value) {
  static_assert(
    sizeof(T) == 0,
    "Unimplemented `new_vector` specialisation"
  );
  return T{};
}
// One-parameter template version
template <typename T>
inline SEXP new_vector(R_xlen_t n) {
  static_assert(
    sizeof(T) == 0,
    "Unimplemented `new_vector` specialisation"
  );
  return T{};
}

template <>
inline SEXP new_vector<r_bool_t>(R_xlen_t n) {
  return internal::new_vec(LGLSXP, n);
}
template <>
inline SEXP new_vector<r_bool_t>(R_xlen_t n, const r_bool_t default_value) {
  SEXP out = SHIELD(new_vector<r_bool_t>(n));
  r_bool_t* RESTRICT p_out = r_ptr::logical_ptr(out);
  fast_fill(p_out, p_out + n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<int>(R_xlen_t n){
  return internal::new_vec(INTSXP, n);
}
template <>
inline SEXP new_vector<int>(R_xlen_t n, const int default_value){
  SEXP out = SHIELD(new_vector<int>(n));
  int* RESTRICT p_out = r_ptr::integer_ptr(out);
  fast_fill(p_out, p_out + n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<double>(R_xlen_t n){
  return internal::new_vec(REALSXP, n);
}
template <>
inline SEXP new_vector<double>(R_xlen_t n, const double default_value){
  SEXP out = SHIELD(new_vector<double>(n));
  double* RESTRICT p_out = r_ptr::real_ptr(out);
  fast_fill(p_out, p_out + n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<int64_t>(R_xlen_t n){
  SEXP out = SHIELD(new_vector<double>(n));
  attr::set_class(out, SHIELD(internal::make_utf8_strsxp("integer64")));
  YIELD(2);
  return out;
}
template <>
inline SEXP new_vector<int64_t>(R_xlen_t n, const int64_t default_value){
  SEXP out = SHIELD(new_vector<int64_t>(n));
  int64_t* RESTRICT p_out = r_ptr::integer64_ptr(out);
  fast_fill(p_out, p_out + n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<r_string_t>(R_xlen_t n){
  return internal::new_vec(STRSXP, n);
}
template <>
inline SEXP new_vector<r_string_t>(R_xlen_t n, const r_string_t default_value){
  SEXP out = SHIELD(new_vector<r_string_t>(n));
  if (default_value != blank_r_string){
    for (R_xlen_t i = 0; i < n; ++i){
      SET_STRING_ELT(out, i, default_value);
    }
  }
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<Rcomplex>(R_xlen_t n){
  return internal::new_vec(CPLXSXP, n);
}
template <>
inline SEXP new_vector<Rcomplex>(R_xlen_t n, const Rcomplex default_value){
  SEXP out = SHIELD(new_vector<Rcomplex>(n));
  Rcomplex* RESTRICT p_out = r_ptr::complex_ptr(out);
  fast_fill(p_out, p_out + n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<Rbyte>(R_xlen_t n){
  return internal::new_vec(RAWSXP, n);
}
template <>
inline SEXP new_vector<Rbyte>(R_xlen_t n,const Rbyte default_value){
  SEXP out = SHIELD(new_vector<Rbyte>(n));
  Rbyte *p_out = r_ptr::raw_ptr(out);
  fast_fill(p_out, p_out + n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<SEXP>(R_xlen_t n){
  return internal::new_vec(VECSXP, n);
}
template <>
inline SEXP new_vector<SEXP>(R_xlen_t n, const SEXP default_value){
  SEXP out = SHIELD(internal::new_vec(VECSXP, n));
  if (!is_null(default_value)){
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, default_value);
    }
  }
  YIELD(1);
  return out;
}

inline SEXP new_list(R_xlen_t n){
  return new_vector<SEXP>(n);
}
inline SEXP new_list(R_xlen_t n, const SEXP default_value){
  return new_vector<SEXP>(n, default_value);
}

}

namespace df {

inline bool is_df(SEXP x){
  return vec::is_df(x);
}

inline int nrow(SEXP x){
  return Rf_length(attr::get_attr(x, symbol::row_names_sym));
}
inline int ncol(SEXP x){
  return Rf_length(x);
}
inline SEXP new_row_names(int n){
  if (n > 0){
    SEXP out = SHIELD(vec::new_vector<int>(2));
    vec::set_val(out, 0, na::integer);
    vec::set_val(out, 1, -n);
    YIELD(1);
    return out;
  } else {
    return vec::new_vector<int>(0);
  }
}
inline void set_row_names(SEXP x, int n){
  SEXP row_names = SHIELD(new_row_names(n));
  attr::set_attr(x, symbol::row_names_sym, row_names);
  YIELD(1);
}
}

template <typename T>
inline constexpr bool is_r_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_inf<double>(const double x){
  return x == r_limits::r_pos_inf || x == r_limits::r_neg_inf;
}

template <typename T>
inline constexpr bool is_r_pos_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_pos_inf<double>(const double x){
  return x == r_limits::r_pos_inf;
}

template <typename T>
inline constexpr bool is_r_neg_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_neg_inf<double>(const double x){
  return x == r_limits::r_neg_inf;
}

// C++ templates

template<typename T>
inline constexpr bool is_r_na(const T x) {
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

template<>
inline constexpr bool is_r_na<r_symbol_t>(const r_symbol_t x){
  return false;
}

// Only CHARSXP has NA
template<>
inline bool is_r_na<SEXP>(const SEXP x){
  return x == na::string;
}

template<>
inline bool is_r_na<r_string_t>(const r_string_t x){
  return x == na::string;
}

// NA type

template<typename T>
inline constexpr T na_value(const T x) {
  static_assert(
    sizeof(T) == 0,
    "Unimplemented `na_value` specialisation"
  );
  return T{};
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
  return na::real;
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
inline r_string_t na_value<r_string_t>(const r_string_t x){
  return na::string;
}

template<>
inline SEXP na_value<SEXP>(const SEXP x){
  switch (TYPEOF(x)){
  case CHARSXP: {
    return na::string;
  }
  default: {
   Rf_error("No `na_value` specialisation for R type %s", Rf_type2char(TYPEOF(x)));
  }
  }
}


namespace vec {

template<typename T>
inline SEXP as_vector(const T x){
  if constexpr (is_r_integral_v<T>){
    return as_vector<int64_t>(x);
  } else if constexpr (std::is_convertible_v<T, SEXP>){
    return as_vector<SEXP>(x);
  } else {
    static_assert(
      sizeof(T) == 0,
      "Unimplemented `as_vector` specialisation"
    );
    return T{};
  }
}
template<>
inline SEXP as_vector<bool>(const bool x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vector<r_bool_t>(const r_bool_t x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vector<Rboolean>(const Rboolean x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vector<int>(const int x){
  return Rf_ScalarInteger(x);
}
template<>
inline SEXP as_vector<int64_t>(const int64_t x){
  if (is_r_na(x)){
    return Rf_ScalarInteger(na::integer);
  } else if (between<int64_t>(x, r_limits::r_int_min, r_limits::r_int_max)){
    return Rf_ScalarInteger(static_cast<int>(x));
  } else {
    return Rf_ScalarReal(static_cast<double>(x));
  }
}
template<>
inline SEXP as_vector<double>(const double x){
  return Rf_ScalarReal(x);
}
template<>
inline SEXP as_vector<Rcomplex>(const Rcomplex x){
  return Rf_ScalarComplex(x);
}
template<>
inline SEXP as_vector<Rbyte>(const Rbyte x){
  return Rf_ScalarRaw(x);
}
template<>
inline SEXP as_vector<const char *>(const char * const x){
  return internal::make_utf8_strsxp(x);
}
template<>
inline SEXP as_vector<std::string>(const std::string x){
  return as_vector<const char *>(x.c_str());
}
template<>
inline SEXP as_vector<cpp11::r_string>(const cpp11::r_string x){
  return Rf_ScalarString(x);
}
template<>
inline SEXP as_vector<cpp11::r_bool>(const cpp11::r_bool x){
  return Rf_ScalarLogical(static_cast<int>(x));
}
template<>
inline SEXP as_vector<r_symbol_t>(const r_symbol_t x){
  return new_vector<SEXP>(1, x);
}

// Scalar string
template<>
inline SEXP as_vector<r_string_t>(const r_string_t x){
  return Rf_ScalarString(x);
}
template<>
inline SEXP as_vector<SEXP>(const SEXP x){
  switch (TYPEOF(x)){
  case CHARSXP: {
    return as_vector<r_string_t>(static_cast<r_string_t>(x));
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
    return new_vector<SEXP>(1, x);
  }
  }
}

}

template<typename T>
inline SEXP as_r_obj(const T x){
  if constexpr (std::is_convertible_v<T, SEXP>){
    switch (TYPEOF(x)){
    case CHARSXP: {
      return vec::as_vector<r_string_t>(static_cast<r_string_t>(x));
    }
    default: {
      return x;
    }
    }
  } else {
    return vec::as_vector(x);
  }
}

template<>
inline SEXP as_r_obj<r_string_t>(const r_string_t x){
  return vec::as_vector<r_string_t>(x);
}
template<>
inline SEXP as_r_obj<r_symbol_t>(const r_symbol_t x){
  return static_cast<SEXP>(x);
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
  if constexpr (std::is_same_v<std::decay_t<T>, int>){
    return static_cast<r_bool_t>(x);
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) ? na::logical : static_cast<r_bool_t>(static_cast<bool>(x));
  } else {
    return na::logical;
  }
}
template<typename T>
inline constexpr int as_int(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, int>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) || !internal::can_be_int(x) ? na::integer : static_cast<int>(x);
  } else {
    return na::integer;
  }
}
template<typename T>
inline constexpr int64_t as_int64(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, int64_t>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) || !internal::can_be_int64(x) ? na::integer64 : static_cast<int64_t>(x);
  } else {
    return na::integer64;
  }
}
template<typename T>
inline constexpr double as_double(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, double>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return is_r_na(x) ? na::real : static_cast<double>(x);
  } else {
    return na::real;
  }
}
template<typename T>
inline constexpr Rcomplex as_complex(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, Rcomplex>){
    return x;
  } else if constexpr (is_r_arithmetic_v<T>){
    return {{as_double(x), 0.0}};
  } else {
    return na::complex;
  }
}
template<typename T>
inline constexpr Rbyte as_raw(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, Rbyte>){
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
inline r_string_t as_r_string(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, r_string_t>){
    return x;
  } else if constexpr (std::is_same_v<std::decay_t<T>, const char *>){
    return static_cast<r_string_t>(internal::make_utf8_charsxp(x));
  } else if constexpr (std::is_same_v<std::decay_t<T>, r_symbol_t>){
    return static_cast<r_string_t>(PRINTNAME(static_cast<SEXP>(x)));
  } else {
    SEXP scalar = SHIELD(as_r_obj(x));
    switch (TYPEOF(scalar)){
    case CHARSXP: {
      YIELD(1);
      return static_cast<r_string_t>(scalar);
    }
    case SYMSXP: {
      r_string_t out = static_cast<r_string_t>(SHIELD(PRINTNAME(scalar)));
      YIELD(2);
      return out;
    }
    default: {
      if (Rf_length(scalar) != 1){
      YIELD(1);
      Rf_error("`x` is a non-scalar vector and cannot be convered to an `r_string_t` in %s", __func__);
    }
      if (is_r_na(x)){
        YIELD(1);
        return na::string;
      }
      SEXP str = SHIELD(vec::coerce_vec(scalar, STRSXP));
      r_string_t out = vec::get_r_string(str, 0);
      YIELD(2);
      return out;
    }
    }

  }
}

// As SYMSXP
template<typename T>
inline r_symbol_t as_r_sym(T x){
  if constexpr (std::is_same_v<std::decay_t<T>, r_symbol_t>){
    return x;
  } else if constexpr (std::is_same_v<std::decay_t<T>, const char *>){
    return static_cast<r_symbol_t>(Rf_installChar(internal::make_utf8_charsxp(x)));
  } else if constexpr (std::is_same_v<std::decay_t<T>, r_string_t>){
    return as_r_sym(CHAR(static_cast<SEXP>(x)));
  } else {
    SEXP scalar = SHIELD(as_r_obj(x));
    switch (TYPEOF(scalar)){
    case SYMSXP: {
      YIELD(1);
      return static_cast<r_symbol_t>(scalar);
    }
    default: {
      if (Rf_length(scalar) != 1){
      YIELD(1);
      Rf_error("`x` is a non-scalar vector and cannot be convered to an `r_symbol_t` in %s", __func__);
    }
      SEXP str = SHIELD(vec::coerce_vec(scalar, STRSXP));
      r_symbol_t out = as_r_sym(vec::get_r_string(str, 0));
      YIELD(2);
      return out;
    }
    }

  }
}

namespace internal {
// R version of static_cast
template<typename T, typename U>
struct r_cast_impl {
  static T cast(U x) {
    static_assert(
      sizeof(T) == 0,
      "Can't `r_cast` this type, use `static_cast`"
    );
    return T{};
  }
};

// Specializations for each target type
template<typename U>
struct r_cast_impl<r_bool_t, U> {
  static constexpr r_bool_t cast(U x) {
    return as_bool(x);
  }
};

template<typename U>
struct r_cast_impl<int, U> {
  static constexpr int cast(U x) {
    return as_int(x);
  }
};

template<typename U>
struct r_cast_impl<int64_t, U> {
  static constexpr int64_t cast(U x) {
    return as_int64(x);
  }
};

template<typename U>
struct r_cast_impl<double, U> {
  static constexpr double cast(U x) {
    return as_double(x);
  }
};

template<typename U>
struct r_cast_impl<Rcomplex, U> {
  static constexpr Rcomplex cast(U x) {
    return as_complex(x);
  }
};

template<typename U>
struct r_cast_impl<Rbyte, U> {
  static constexpr Rbyte cast(U x) {
    return as_raw(x);
  }
};

template<typename U>
struct r_cast_impl<r_string_t, U> {
  static r_string_t cast(U x) {
    return as_r_string(x);
  }
};

template<typename U>
struct r_cast_impl<r_symbol_t, U> {
  static r_symbol_t cast(U x) {
    return as_r_sym(x);
  }
};

template<typename U>
struct r_cast_impl<SEXP, U> {
  static constexpr SEXP cast(U x) {
    return as_r_obj(x);
  }
};
}

template<typename T, typename U>
inline constexpr T r_cast(U x) {
  return internal::r_cast_impl<T, U>::cast(x);
}


// R math fns
namespace internal {

inline double round_to_even(double x){
  return x - std::remainder(x, 1.0);
}

}

namespace math {

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
    return na::real;
  } else {
    return std::sqrt(x.r * x.r + x.i * x.i);
  }
}

template<typename T>
inline T r_floor(T x){
  return is_r_na(x) ? x : std::floor(x);
}
template<>
inline double r_floor(double x){
  return std::floor(x);
}

template<typename T>
inline T r_ceiling(T x){
  return is_r_na(x) ? x : std::ceil(x);
}
template<>
inline double r_ceiling(double x){
  return std::ceil(x);
}

template<typename T>
inline T r_trunc(T x){
  return is_r_na(x) ? x : std::trunc(x);
}

template <>
inline double r_trunc(double x){
  return std::trunc(x) + 0.0;
}

template <typename T>
inline int r_sign(T x) {
  return is_r_na(x) ? na::integer : (T(0) < x) - (x < T(0));
}

template<typename T>
inline T r_negate(T x){
  return is_r_na(x) ? x : -x;
}
template<>
inline double r_negate(double x){
  return -x;
}

template<typename T>
inline double r_sqrt(T x){
  return std::sqrt(r_cast<double>(x));
}

template<typename T>
inline double r_pow(T x, T y){
  if (y == 2){
    double left = r_cast<double>(x);
    return left * left;
  } else {
    return std::pow(r_cast<double>(x), r_cast<double>(y));
  }
}

template<typename T>
inline double r_log10(T x){
  return std::log10(r_cast<double>(x));
}

template<typename T>
inline double r_exp(T x){
  return std::exp(r_cast<double>(x));
}

template<typename T, typename U>
inline double r_log(T x, U base){
  return std::log(r_cast<double>(x)) / std::log(r_cast<double>(base));
}
template<typename T>
inline double r_log(T x){
  return std::log(r_cast<double>(x));
}
inline Rcomplex r_log(Rcomplex x){
  if (is_r_na(x)){
    return x;
  }
  Rcomplex out;
  out.r = 0.5 * (r_log(r_pow(x.r, 2.0) + r_pow(x.i, 2.0)));
  out.i = std::atan2(r_cast<double>(x.i), r_cast<double>(x.r));
  return out;
}


template<typename T, typename U>
inline double r_round(T x, const U digits){
  if (is_r_na(x)){
    return r_cast<double>(x);
  } else if (is_r_na(digits)){
    return na::real;
  } else if (is_r_inf(x)){
    return x;
  } else if (is_r_neg_inf(digits)){
    return 0.0;
  } else if (is_r_pos_inf(digits)){
    return x;
  } else {
    double scale = std::pow(10.0, r_cast<double>(digits));
    return internal::round_to_even(r_cast<double>(x) * scale) / scale;
  }
}

template<typename T>
inline double r_round(T x){
  if (is_r_na(x)){
    return r_cast<double>(x);
  } else if (is_r_inf(x)){
    return x;
  } else {
    return internal::round_to_even(r_cast<double>(x));
  }
}

template<typename T, typename U>
inline double r_signif(T x, const U digits){
  U new_digits = std::max(U(1), digits);
  if (is_r_na(x)){
    return r_cast<double>(x);
  } else if (is_r_na(new_digits)){
    return na::real;
  } else if (is_r_pos_inf(digits)){
    return x;
  } else {
    new_digits -= std::ceil(std::log10(std::abs(x)));
    double scale = std::pow(10, new_digits);
    return internal::round_to_even(scale * x) / scale;
  }
}

template<typename T>
inline T r_add(T x, T y){
  return is_r_na(x) || is_r_na(y) ? na_value(x) : x + y;
}
template<>
inline double r_add(double x, double y){
  return x + y;
}
template<typename T>
inline T r_subtract(T x, T y){
  return is_r_na(x) || is_r_na(y) ? na_value(x) : x - y;
}
template<>
inline double r_subtract(double x, double y){
  return x - y;
}
template<typename T>
inline T r_multiply(T x, T y){
  return is_r_na(x) || is_r_na(y) ? na_value(x) : x * y;
}
template<>
inline double r_multiply(double x, double y){
  return x * y;
}
template<typename T, typename U>
inline double r_divide(T x, U y){
  return r_cast<double>(x) / r_cast<double>(y);
}

template<typename T>
inline T r_abs_diff(const T x, const T y){
  return std::abs(r_subtract(x, y));
}

inline r_bool_t is_whole_number(const double x, const double tolerance){
  return is_r_na(x) || is_r_na(tolerance) ? na::logical : static_cast<r_bool_t>(r_abs_diff(x, std::round(x)) <= tolerance);
}

template<
  typename T,
  typename = typename std::enable_if<is_r_arithmetic_v<T>>::type
>
inline T r_gcd(T x, T y, bool na_rm = true, T tol = std::sqrt(std::numeric_limits<T>::epsilon())){

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
inline T r_lcm(
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
    T res = std::abs(x) / r_gcd(x, y, na_rm);
    if (y != 0 && (std::abs(res) > (std::numeric_limits<T>::max() / std::abs(y)))){
      return na_value(x);
    }
    return res * std::abs(y);
  } else {
    if (std::fabs(x) <= tol && std::fabs(y) <= tol){
      return 0.0;
    }
    return ( std::fabs(x) / r_gcd(x, y, na_rm, tol) ) * std::fabs(y);
  }
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
    expr = SHIELD(Rf_lang3(symbol::triple_colon_sym, r_cast<r_symbol_t>(pkg), r_cast<r_symbol_t>(name)));
  } else {
    expr = SHIELD(Rf_lang3(symbol::double_colon_sym, r_cast<r_symbol_t>(pkg), r_cast<r_symbol_t>(name)));
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

// Variadic list constructor
template<typename... Args>
inline SEXP make_list(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return new_vector<SEXP>(0);
  } else {
    SEXP out = SHIELD(vec::new_vector<SEXP>(n));

    // Are any args named?
    constexpr bool any_named = (std::is_same_v<std::decay_t<Args>, arg> || ...);

    SEXP nms = r_null;

    if (any_named){
      nms = vec::new_vector<r_string_t>(n);
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
        SET_TAG(current, r_cast<r_symbol_t>(args.name));
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
      const SEXP *p_x = r_ptr::list_ptr_ro(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
    } else {
      if (is_null(internal::r_length_sym)){
        internal::r_length_sym = r_cast<r_symbol_t>("length");
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
      internal::r_length_sym = r_cast<r_symbol_t>("length");
    }
    SEXP expr = SHIELD(Rf_lang2(internal::r_length_sym, x));
    SEXP r_len = SHIELD(eval(expr, R_GetCurrentEnv()));
    R_xlen_t out = TYPEOF(r_len) == INTSXP ? INTEGER_ELT(r_len, 0) : REAL_ELT(r_len, 0);
    YIELD(2);
    return out;
  }
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
    return vec::new_vector<int>(0);
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
  R_xlen_t na_count = 0;

  switch ( internal::CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP:
  case internal::CHEAPR_INT64SXP: {
    break;
  }
  case REALSXP: {
    const double *p_x = r_ptr::real_ptr_ro(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      out = static_cast<r_bool_t>(math::is_whole_number(p_x[i], tol_));
      na_count += is_r_na(out);
      if (out == r_false){
        break;
      }
    }
    if (out && !na_rm_ && na_count > 0){
      out = r_na;
    } else if (na_rm_ && na_count == n){
      out = r_true;
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
    const SEXP *p_attributes = r_ptr::list_ptr_ro(attrs);
    const r_string_t *p_names = r_ptr::string_ptr_ro(names);

    r_symbol_t attr_nm;

    for (int i = 0; i < Rf_length(names); ++i){
      if (p_names[i] != blank_r_string){
        attr_nm = r_cast<r_symbol_t>(p_names[i]);
        if (address(x) == address(p_attributes[i])){
          SEXP dup_attr = SHIELD(Rf_duplicate(p_attributes[i])); ++NP;
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
      if (is_null(symbol::tag(current)) || r_cast<r_string_t>(symbol::tag(current)) == blank_r_string){
        YIELD(NP);
        Rf_error("Please only supply named attributes in %s", __func__);
      }
      if (addr_x == address(CAR(current))){
        SEXP dup_attr = SHIELD(Rf_duplicate(CAR(current))); ++NP;
        attr::set_attr(x, symbol::tag(current), dup_attr);
      } else {
        attr::set_attr(x, symbol::tag(current), CAR(current));
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

  if (is_null(a)){
    return r_null;
  }

  int n = Rf_length(a);

  SEXP out = SHIELD(vec::new_list(n));
  SEXP names = SHIELD(vec::new_vector<r_string_t>(n));
  SEXP current = a;

  for (int i = 0; i < n; ++i){
    vec::set_val(out, i, CAR(current));
    if (!is_null(symbol::tag(current))){
      vec::set_val(names, i, r_cast<r_string_t>(symbol::tag(current)));
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
  if (!is_null(x)){
    clear_attrs(x);
    internal::add_attrs(x, attrs);
  }
}

}

namespace vec {
inline SEXP deep_copy(SEXP x){
  return Rf_duplicate(x);
  // int32_t NP = 0;
  // SEXP out = r_null;
  // R_xlen_t n = Rf_xlength(x);
  // SEXP attrs = r_null;
  //
  // switch (TYPEOF(x)){
  // case NILSXP: {
  //   break;
  // }
  // case LGLSXP: {
  //   out = SHIELD(new_vector<r_bool_t>(n)); ++NP;
  //   fast_copy_n(r_ptr::logical_ptr_ro(x), n, r_ptr::logical_ptr(out));
  //   break;
  // }
  // case INTSXP: {
  //   out = SHIELD(new_vector<int>(n)); ++NP;
  //   fast_copy_n(r_ptr::integer_ptr_ro(x), n, r_ptr::integer_ptr(out));
  //   break;
  // }
  // case REALSXP: {
  //   out = SHIELD(new_vector<double>(n)); ++NP;
  //   fast_copy_n(r_ptr::real_ptr_ro(x), n, r_ptr::real_ptr(out));
  //   break;
  // }
  // case STRSXP: {
  //   out = SHIELD(new_vector<r_string_t>(n)); ++NP;
  //   const r_string_t *p_x = r_ptr::string_ptr_ro(x);
  //   for (R_xlen_t i = 0; i < n; ++i){
  //     SET_STRING_ELT(out, i, p_x[i]);
  //   }
  //   break;
  // }
  // case CPLXSXP: {
  //   out = SHIELD(new_vector<Rcomplex>(n)); ++NP;
  //   fast_copy_n(r_ptr::complex_ptr_ro(x), n, r_ptr::complex_ptr(out));
  //   break;
  // }
  // case RAWSXP: {
  //   out = SHIELD(new_vector<Rbyte>(n)); ++NP;
  //   fast_copy_n(r_ptr::raw_ptr_ro(x), n, r_ptr::raw_ptr(out));
  //   break;
  // }
  // case VECSXP: {
  //   out = SHIELD(new_vector<SEXP>(n)); ++NP;
  //   const SEXP *p_x = r_ptr::list_ptr_ro(x);
  //   for (R_xlen_t i = 0; i < n; ++i){
  //     SET_VECTOR_ELT(out, i, deep_copy(p_x[i]));
  //   }
  //   break;
  // }
  // default: {
  //   out = SHIELD(Rf_duplicate(x)); ++NP;
  //   YIELD(NP);
  //   return out;
  // }
  // }
  //
  // if (!is_null(x)){
  //   SHIELD(attrs = attr::get_attrs(x)); ++NP;
  //   int n_attrs = Rf_length(attrs);
  //   for (R_xlen_t i = 0; i < n_attrs; ++i){
  //     SET_VECTOR_ELT(attrs, i, deep_copy(VECTOR_ELT(attrs, i)));
  //   }
  //   attr::set_attrs(out, attrs);
  // }
  //
  // YIELD(NP);
  // return out;
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

// We call R fn`cheapr::set_threads` to make sure the R option is set
inline void set_threads(uint16_t n){
  uint16_t max_threads = OMP_MAX_THREADS;
  uint16_t threads = std::min(n, max_threads);
  SEXP cheapr_set_threads = SHIELD(fn::find_pkg_fun("set_threads", "cheapr", true));
  SEXP r_threads = SHIELD(vec::as_vector(r_cast<int>(threads)));
  SHIELD(fn::eval_fn(cheapr_set_threads, R_BaseEnv, r_threads));
  YIELD(3);
}


namespace internal {

// A type that signals the `SEXP` type is unsupported for the
// current calculation
struct unsupported_sexp_t {
  SEXP value;
  explicit unsupported_sexp_t(SEXP x) : value(x) {}
};

// Retrieve the pointer of x
// To be used in a lambda
// E.g. with_read_only_data(x, [&](auto p_x) {})
template <class F>
decltype(auto) with_read_only_data(SEXP x, F&& f) {
  switch (CHEAPR_TYPEOF(x)) {
  case NILSXP:          return f(static_cast<const int*>(nullptr));
  case LGLSXP:          return f(r_ptr::logical_ptr_ro(x));
  case INTSXP:          return f(r_ptr::integer_ptr_ro(x));
  case CHEAPR_INT64SXP: return f(r_ptr::integer64_ptr_ro(x));
  case REALSXP:         return f(r_ptr::real_ptr_ro(x));
  case STRSXP:          return f(r_ptr::string_ptr_ro(x));
  case VECSXP:          return f(r_ptr::list_ptr_ro(x));
  case CPLXSXP:         return f(r_ptr::complex_ptr_ro(x));
  case RAWSXP:          return f(r_ptr::raw_ptr_ro(x));
  default:              return f(unsupported_sexp_t{x});
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
