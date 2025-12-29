#ifndef CHEAPR_C_CORE_H
#define CHEAPR_C_CORE_H

// cheapr Core definitions and templates
// Author: Nick Christofides
// License: MIT

// template<typename T>
// struct r_vector_t {
//   SEXP sexp;
  
//   // Constructor that wraps new_vector<T>
//   explicit r_vector_t(R_xlen_t size) : sexp(new_vector<T>(size)) {}
  
//   // Constructor from existing SEXP
//   explicit r_vector_t(SEXP s) : sexp(s) {}
  
//   // Get/set using your existing templates
//   void set(R_xlen_t index, T value) {
//     set_value<T>(sexp, index, value);
//   }
  
//   T get(R_xlen_t index) const {
//     return get_value<T>(sexp, index);
//   }
  
//   // Implicit conversion to SEXP
//   operator SEXP() const {
//     return sexp;
//   }
  
//   // Get size
//   R_xlen_t size() const {
//     return Rf_xlength(sexp);
//   }
// };

// Can also play around with automatic protection (unlikely with this method..)
// template<typename T>
// struct r_vector_t {
//   SEXP sexp;
//   bool owned;
  
//   explicit r_vector_t(R_xlen_t size) 
//     : sexp(PROTECT(new_vector<T>(size))), owned(true) {}
  
//   explicit r_vector_t(SEXP s, bool protect = false) 
//     : sexp(protect ? PROTECT(s) : s), owned(protect) {}
  
//   ~r_vector_t() {
//     if (owned) UNPROTECT(1);
//   }
  
//   // Disable copying (Rule of Three)
//   r_vector_t(const r_vector_t&) = delete;
//   r_vector_t& operator=(const r_vector_t&) = delete;
  
//   // Methods as before...
// };

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
#define SHIELD(x)                                              \
static_cast<std::remove_cvref_t<decltype(x)>>(Rf_protect(static_cast<SEXP>(x)))
#endif

#ifndef YIELD
#define YIELD(n) (Rf_unprotect(n))
#endif

namespace cheapr {

// Constants

// R C NULL
inline const SEXP r_null = R_NilValue;

namespace internal {

inline constexpr int64_t CHEAPR_OMP_THRESHOLD = 100000;
inline constexpr SEXPTYPE CHEAPR_INT64SXP = 64;
}

// An idea for unified wrapper type
// template<typename T>
// struct r_wrapper_t {
// private:
//   T val;
  
// public:
//   // Construction from underlying type
//   explicit r_wrapper_t(T v) : val(v) {}
  
//   // Implicit conversion to underlying type
//   operator T&() { return val; }
//   operator const T&() const { return val; }
  
//   // Explicit extraction if needed
//   T& get() { return val; }
//   const T& get() const { return val; }
// };

// using r_int_t = r_wrapper_t<int>;
// using r_double_t = r_wrapper_t<double>;

// bool type, similar to Rboolean

struct r_bool_t {
  int value;
  r_bool_t() : value{0} {}
  explicit constexpr r_bool_t(bool x) : value{x} {}
  explicit constexpr r_bool_t(int x) : value{x} {}
  constexpr operator int() const { return value; }
  constexpr operator bool() const { return value; }
};

inline constexpr r_bool_t r_true{1};
inline constexpr r_bool_t r_false{0};
inline constexpr r_bool_t r_na{std::numeric_limits<int>::min()};

struct r_int_t {
  int value;
  r_int_t() : value{0} {}
  constexpr r_int_t(int x) : value{x} {}
  constexpr operator int() const { return value; }
};

struct r_double_t {
  double value;
  r_double_t() : value{0} {}
  constexpr r_double_t(double x) : value{x} {}
  constexpr operator double() const { return value; }
};

// Alias type for CHARSXP
struct r_string_t {
  SEXP value;
  r_string_t() : value{R_BlankString} {}
  // Explicit SEXP -> r_string_t
  explicit constexpr r_string_t(SEXP x) : value{x} {}
  // Implicit r_string_t -> SEXP
  constexpr operator SEXP() const { return value; }
};

inline const r_string_t blank_r_string = r_string_t();

// Alias type for SYMSXP
struct r_symbol_t {
  SEXP value;
  r_symbol_t() : value{R_MissingArg} {}
  explicit constexpr r_symbol_t(SEXP x) : value{x} {}
  constexpr operator SEXP() const { return value; }
};


// Alias type for Rcomplex
struct r_complex_t {
  Rcomplex value;

  // Constructors
  constexpr r_complex_t() : value{0.0, 0.0} {}
  constexpr r_complex_t(r_double_t r, r_double_t i) : value{r, i} {}

  // Conversion handling
  explicit constexpr r_complex_t(Rcomplex x) : value{x} {}
  constexpr operator Rcomplex() const { return value; }

  // Get real and imaginary parts
    constexpr r_double_t re() const { return r_double_t{value.r}; }
    constexpr r_double_t im() const { return r_double_t{value.i}; }
  // constexpr r_double_t re() { return r_double_t{value.r}; }
  // constexpr r_double_t im() { return r_double_t{value.i}; }
  // constexpr const r_double_t& re() const { return value.r; }
  // constexpr const r_double_t& im() const { return value.i; }

  // Equality operators
  // constexpr r_bool_t operator==(const r_complex_t& other) const {
  //   if (re() != re() || im() != im() || other.re() != other.re() || other.im() != other.im()){
  //     return r_na;
  //   } else {
  //     return r_bool_t{re() == other.re() && im() == other.im()};
  //   }
  // }
  //
  // constexpr r_bool_t operator!=(const r_complex_t& other) const {
  //   if (re() != re() || im() != im() || other.re() != other.re() || other.im() != other.im()){
  //     return r_na;
  //   } else {
  //     return r_bool_t{!(re() == other.re() && im() == other.im())};
  //   }
  // }

};

// Alias type for r_byte_t
struct r_byte_t {
  Rbyte value;

  // Constructors
  constexpr r_byte_t() : value{static_cast<Rbyte>(0)} {}

  // Conversion handling
  explicit constexpr r_byte_t(Rbyte x) : value{x} {}
  constexpr operator Rbyte() const { return value; }
};

namespace r_limits {
inline constexpr r_int_t r_int_min = r_int_t{-std::numeric_limits<int>::max()};
inline constexpr r_int_t r_int_max = r_int_t{std::numeric_limits<int>::max()};
inline constexpr int64_t r_int64_min = -std::numeric_limits<int64_t>::max();
inline constexpr int64_t r_int64_max = std::numeric_limits<int64_t>::max();
inline constexpr r_double_t r_pos_inf = r_double_t{std::numeric_limits<double>::infinity()};
inline constexpr r_double_t r_neg_inf = r_double_t{-std::numeric_limits<double>::infinity()};
}

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
inline r_symbol_t unbound_value = static_cast<r_symbol_t>(R_UnboundValue);

}

// Compile-time type check `is<>`
template<typename T, typename U>
    inline constexpr bool is = std::is_same_v<std::remove_cvref_t<T>, std::remove_cvref_t<U>>;

// Concepts to enable R type templates
template<typename T>
concept RMathType = std::same_as<std::remove_cvref_t<T>, r_bool_t> ||
std::same_as<std::remove_cvref_t<T>, r_int_t> ||
std::same_as<std::remove_cvref_t<T>, r_double_t>;

template<typename T>
concept RType = RMathType<T> ||
std::same_as<std::remove_cvref_t<T>, r_complex_t> ||
std::same_as<std::remove_cvref_t<T>, r_string_t> ||
std::same_as<std::remove_cvref_t<T>, r_symbol_t> ||
std::same_as<std::remove_cvref_t<T>, SEXP>;


template <class... T>
inline constexpr bool always_false = false;

template <typename T>
inline constexpr bool is_r_or_cpp_integral_v =
std::is_integral_v<T> ||
is<T, r_bool_t> ||
is<T, r_int_t>;

template <typename T>
inline constexpr bool is_r_or_cpp_arithmetic_v =
is_r_or_cpp_integral_v<T> ||
std::is_arithmetic_v<T> ||
is<T, r_double_t>;


template <typename T>
inline constexpr bool is_r_ptr_writable_v =
is_r_or_cpp_arithmetic_v<T> ||
is<T, r_complex_t> ||
is<T, r_byte_t>;


template<typename T>
concept MathType = is_r_or_cpp_arithmetic_v<T>;

namespace env {
inline const SEXP empty_env = R_EmptyEnv;
inline const SEXP base_env = R_BaseEnv;
}

// Get the R type from the C type (useful for mixed type operators)

namespace internal {
template<typename T>
inline auto as_r_type(T x) {
  if constexpr (is<T, bool>) {
    return r_bool_t{x};
  } else if constexpr (is<T, int>) {
    return r_int_t{x};
  } else if constexpr (is<T, double>) {
    return r_double_t{x};    
  } else if constexpr (is<T, r_string_t>) {
    return r_string_t{x};
  } else if constexpr (is<T, Rcomplex>) {
    return r_complex_t{x};
  }  else if constexpr (is<T, Rbyte>) {
    return r_byte_t{x};
  } else {
  static_assert(
    always_false<T>,
    "Unsupported type for `as_r_type`"
  );
  }
}

}

// NAs

namespace na {
inline constexpr r_bool_t logical = r_na;
inline constexpr r_int_t integer = r_int_t{std::numeric_limits<int>::min()};
inline constexpr int64_t integer64 = std::numeric_limits<int64_t>::min();
inline const r_double_t real = r_double_t{NA_REAL};
inline const r_complex_t complex = r_complex_t{NA_REAL, NA_REAL};
inline constexpr r_byte_t raw = r_byte_t{0};
inline const r_string_t string = r_string_t{NA_STRING};
inline const SEXP nil = r_null;
}

namespace internal {

inline r_int_t* integer_ptr(SEXP x){
  return reinterpret_cast<r_int_t*>(INTEGER(x));
}
inline const r_int_t* integer_ptr_ro(SEXP x){
  return reinterpret_cast<const r_int_t*>(INTEGER_RO(x));
}
inline r_bool_t* logical_ptr(SEXP x){
  return reinterpret_cast<r_bool_t*>(integer_ptr(x));
}
inline const r_bool_t* logical_ptr_ro(SEXP x){
  return reinterpret_cast<const r_bool_t*>(integer_ptr_ro(x));
}
inline r_double_t* real_ptr(SEXP x){
  return reinterpret_cast<r_double_t*>(REAL(x));
}
inline const r_double_t* real_ptr_ro(SEXP x){
  return reinterpret_cast<const r_double_t*>(REAL_RO(x));
}
inline int64_t* integer64_ptr(SEXP x){
  return reinterpret_cast<int64_t*>(real_ptr(x));
}
inline const int64_t* integer64_ptr_ro(SEXP x){
  return reinterpret_cast<const int64_t*>(real_ptr_ro(x));
}
inline r_complex_t* complex_ptr(SEXP x){
  return reinterpret_cast<r_complex_t*>(COMPLEX(x));
}
inline const r_complex_t* complex_ptr_ro(SEXP x){
  return reinterpret_cast<const r_complex_t*>(COMPLEX_RO(x));
}
inline r_byte_t* raw_ptr(SEXP x){
  return reinterpret_cast<r_byte_t*>(RAW(x));
}
inline const r_byte_t* raw_ptr_ro(SEXP x){
  return reinterpret_cast<const r_byte_t*>(RAW_RO(x));
}
inline const r_string_t* string_ptr_ro(SEXP x){
  return reinterpret_cast<const r_string_t*>(STRING_PTR_RO(x));
}

}

inline const SEXP* list_ptr_ro(SEXP x){
  return VECTOR_PTR_RO(x);
}

template<typename T>
inline T* vector_ptr(SEXP x) {
  static_assert(
    always_false<T>,
    "Unsupported type for vector_ptr"
  );
  return nullptr;
}


template<> 
inline r_bool_t* vector_ptr<r_bool_t>(SEXP x) {
  return internal::logical_ptr(x);
}

template<>
inline r_int_t* vector_ptr<r_int_t>(SEXP x) {
  return internal::integer_ptr(x);
}

template<>
inline r_double_t* vector_ptr<r_double_t>(SEXP x) {
  return internal::real_ptr(x);
}

template<>
inline int64_t* vector_ptr<int64_t>(SEXP x) {
  return internal::integer64_ptr(x);
}

template<>
inline r_complex_t* vector_ptr<r_complex_t>(SEXP x) {
  return internal::complex_ptr(x);
}

template<>
inline r_byte_t* vector_ptr<r_byte_t>(SEXP x) {
  return internal::raw_ptr(x);
}

template<>
inline const r_bool_t* vector_ptr<const r_bool_t>(SEXP x) {
  return internal::logical_ptr_ro(x);
}

template<>
inline const r_int_t* vector_ptr<const r_int_t>(SEXP x) {
  return internal::integer_ptr_ro(x);
}

template<>
inline const r_double_t* vector_ptr<const r_double_t>(SEXP x) {
  return internal::real_ptr_ro(x);
}

template<>
inline const int64_t* vector_ptr<const int64_t>(SEXP x) {
  return internal::integer64_ptr_ro(x);
}

template<>
inline const r_complex_t* vector_ptr<const r_complex_t>(SEXP x) {
  return internal::complex_ptr_ro(x);
}

template<>
inline const r_string_t* vector_ptr<const r_string_t>(SEXP x) {
  return internal::string_ptr_ro(x);
}

// template<>
// inline const r_byte_t* vector_ptr<const r_byte_t>(SEXP x) {
//   return internal::raw_ptr_ro(x);
// }

template<>
inline const SEXP* vector_ptr<const SEXP>(SEXP x) {
  return list_ptr_ro(x);
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
  return internal::integer64_ptr(x);
}
inline const int64_t* INTEGER64_PTR_RO(SEXP x) {
  return internal::integer64_ptr_ro(x);
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


namespace internal {
inline SEXP BASE_ATTRIBUTES = r_null;
inline SEXP CHEAPR_CORES = r_null;
inline SEXP BASE_LENGTH = r_null;
}

namespace attr {

// Attributes of x as a list
inline SEXP get_attrs(SEXP x){
  if (is_null(internal::BASE_ATTRIBUTES)){
    internal::BASE_ATTRIBUTES = Rf_install("attributes");
  }
  SEXP expr = SHIELD(Rf_lang2(internal::BASE_ATTRIBUTES, x));
  SEXP out = SHIELD(Rf_eval(expr, R_BaseEnv));
  YIELD(2);
  return out;
}

inline bool has_attrs(SEXP x){
  return !is_null(get_attrs(x));
}

inline SEXP get_attr(SEXP x, r_symbol_t sym){
  return Rf_getAttrib(x, static_cast<SEXP>(sym));
}

inline void set_attr(SEXP x, r_symbol_t sym, SEXP value){
  Rf_setAttrib(x, static_cast<SEXP>(sym), value);
}

}

namespace internal {

inline int get_threads(){
  if (is_null(CHEAPR_CORES)){
    CHEAPR_CORES = Rf_install("cheapr.cores");
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
  return !is_object(x);
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

}

namespace internal {

// R vector getters + setters

template<typename T>
inline T get_value_impl(T *p_x, const R_xlen_t i){
  static_assert(
    always_false<T>,
    "Unimplemented `get_value_impl` specialisation"
  );
  return T{};
}

template<typename T>
inline T get_value_impl(SEXP x, const R_xlen_t i){
  static_assert(
    always_false<T>,
    "Unimplemented `get_value_impl` specialisation"
  );
  return T{};
}

template<>
inline r_bool_t get_value_impl<r_bool_t>(r_bool_t* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline r_bool_t get_value_impl<r_bool_t>(SEXP x, const R_xlen_t i){
  return internal::logical_ptr_ro(x)[i];
}
template<>
inline r_int_t get_value_impl<r_int_t>(r_int_t* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline r_int_t get_value_impl<r_int_t>(SEXP x, const R_xlen_t i){
  return internal::integer_ptr_ro(x)[i];
}
template<>
inline int64_t get_value_impl<int64_t>(int64_t* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline int64_t get_value_impl<int64_t>(SEXP x, const R_xlen_t i){
  return internal::integer64_ptr_ro(x)[i];
}
template<>
inline r_double_t get_value_impl<r_double_t>(r_double_t* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline r_double_t get_value_impl<r_double_t>(SEXP x, const R_xlen_t i){
  return internal::real_ptr_ro(x)[i];
}
template<>
inline r_complex_t get_value_impl<r_complex_t>(r_complex_t* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline r_complex_t get_value_impl<r_complex_t>(SEXP x, const R_xlen_t i){
  return internal::complex_ptr_ro(x)[i];
}
template<>
inline r_byte_t get_value_impl<r_byte_t>(r_byte_t* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline r_byte_t get_value_impl<r_byte_t>(SEXP x, const R_xlen_t i){
  return internal::raw_ptr_ro(x)[i];
}
template<>
inline r_string_t get_value_impl<r_string_t>(r_string_t* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline r_string_t get_value_impl<r_string_t>(SEXP x, const R_xlen_t i){
  return internal::string_ptr_ro(x)[i];
}
template<>
inline SEXP get_value_impl<SEXP>(SEXP* p_x, const R_xlen_t i){
  return p_x[i];
}
template<>
inline SEXP get_value_impl<SEXP>(SEXP x, const R_xlen_t i){
  return list_ptr_ro(x)[i];
}

}
namespace vec {

template<typename T>
inline auto get_value(T *p_x, const R_xlen_t i) {
  return internal::get_value_impl<std::remove_cvref_t<T>>(p_x, i);
}

template<typename T>
inline auto get_value(SEXP x, const R_xlen_t i) {
  return internal::get_value_impl<std::remove_cvref_t<T>>(x, i);
}

template<typename T>
inline void set_value(T *p_x, const R_xlen_t i, T val){
  static_assert(
    always_false<T>,
    "Unimplemented `set_value` specialisation"
  );
  return T{};
}

template<typename T>
inline void set_value(SEXP x, const R_xlen_t i, T val){
  static_assert(
    always_false<T>,
    "Unimplemented `set_value` specialisation"
  );
  return T{};
}

template<>
inline void set_value<r_bool_t>(r_bool_t* p_x, const R_xlen_t i, r_bool_t val){
  p_x[i] = val;
}
template<>
inline void set_value<r_bool_t>(SEXP x, const R_xlen_t i, r_bool_t val){
  set_value(internal::logical_ptr(x), i, val);
}
template<>
inline void set_value<r_int_t>(r_int_t* p_x, const R_xlen_t i, r_int_t val){
  p_x[i] = val;
}
template<>
inline void set_value<r_int_t>(SEXP x, const R_xlen_t i, r_int_t val){
  set_value(internal::integer_ptr(x), i, val);
}
template<>
inline void set_value<int64_t>(int64_t* p_x, const R_xlen_t i, int64_t val){
  p_x[i] = val;
}
template<>
inline void set_value<int64_t>(SEXP x, const R_xlen_t i, int64_t val){
  set_value(internal::integer64_ptr(x), i, val);
}
template<>
inline void set_value<r_double_t>(r_double_t* p_x, const R_xlen_t i, r_double_t val){
  p_x[i] = val;
}
template<>
inline void set_value<r_double_t>(SEXP x, const R_xlen_t i, r_double_t val){
  set_value(internal::real_ptr(x), i, val);
}
template<>
inline void set_value<r_complex_t>(r_complex_t* p_x, const R_xlen_t i, r_complex_t val){
  p_x[i].re() = val.re();
  p_x[i].im() = val.im();
}
template<>
inline void set_value<r_complex_t>(SEXP x, const R_xlen_t i, r_complex_t val){
  set_value(internal::complex_ptr(x), i, val);
}
template<>
inline void set_value<r_byte_t>(r_byte_t* p_x, const R_xlen_t i, r_byte_t val){
  p_x[i] = val;
}
template<>
inline void set_value<r_byte_t>(SEXP x, const R_xlen_t i, r_byte_t val){
  set_value(internal::raw_ptr(x), i, val);
}
template<>
inline void set_value<r_string_t>(SEXP x, const R_xlen_t i, r_string_t val){
  SET_STRING_ELT(x, i, static_cast<SEXP>(val));
}
template<>
inline void set_value<const char *>(SEXP x, const R_xlen_t i, const char* val){
  set_value<r_string_t>(x, i, static_cast<r_string_t>(internal::make_utf8_charsxp(val)));
}
template<>
inline void set_value<r_symbol_t>(SEXP x, const R_xlen_t i, r_symbol_t val){
  SET_VECTOR_ELT(x, i, static_cast<SEXP>(val));
}

// Never use the pointer here to assign
template<>
inline void set_value<SEXP>(SEXP x, const R_xlen_t i, SEXP val){
  SET_VECTOR_ELT(x, i, val);
}

}

template <typename T>
inline void r_fill(
    [[maybe_unused]] SEXP target,
    [[maybe_unused]] T *p_target,
    R_xlen_t start, R_xlen_t n, const T val
){

  if constexpr (is_r_ptr_writable_v<T>){
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        vec::set_value(p_target, start + i, val);
      }
    } else {
      std::fill_n(p_target + start, n, val);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      vec::set_value(target, start + i, val);
    }
  }
}

template <typename T>
inline void r_fill(
    [[maybe_unused]] SEXP target,
    [[maybe_unused]] const T *p_target,
    R_xlen_t start, R_xlen_t n, const T val
){
  using r_t = std::remove_const_t<T>;
  r_fill(target, const_cast<r_t*>(p_target), start, n, val);
}

template <typename T>
inline void r_replace(
    [[maybe_unused]] SEXP target,
    [[maybe_unused]] T *p_target,
    R_xlen_t start, R_xlen_t n, const T old_val, const T new_val
){
  if constexpr (is_r_ptr_writable_v<T>){
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        if (p_target[start + i] == old_val){
          vec::set_value(p_target, start + i, new_val);
        }
      }
    } else {
      std::replace(p_target + start, p_target + start + n, old_val, new_val);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      R_xlen_t idx = start + i;
      if (p_target[idx] == old_val){
        vec::set_value(target, idx, new_val);
      }
    }
  }
}

template <typename T>
inline void r_replace(
    [[maybe_unused]] SEXP target,
    [[maybe_unused]] const T *p_target,
    R_xlen_t start, R_xlen_t n, const T old_val, const T new_val
){
  using r_t = std::remove_const_t<T>;
  r_replace(target, const_cast<r_t*>(p_target), start, n, old_val, new_val);
}

template <typename T>
inline void r_copy_n(
    [[maybe_unused]] SEXP target,
    [[maybe_unused]] T *p_target,
    const T *p_source, R_xlen_t target_offset, R_xlen_t n
){

  if constexpr (is_r_ptr_writable_v<T>){
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        vec::set_value(p_target, target_offset + i, p_source[i]);
      }
    } else {
      std::copy_n(p_source, n, p_target + target_offset);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      vec::set_value(target, target_offset + i, p_source[i]);
    }
  }
}
template <typename T>
inline void r_copy_n(
    [[maybe_unused]] SEXP target,
    [[maybe_unused]] const T *p_target,
    const T *p_source, R_xlen_t target_offset, R_xlen_t n
){
  using r_t = std::remove_const_t<T>;
  r_copy_n(target, const_cast<r_t*>(p_target), p_source, target_offset, n);
}


namespace attr {

inline void set_old_names(SEXP x, SEXP names){
  if (is_null(names)){
    attr::set_attr(x, symbol::names_sym, r_null);
  } else {
    Rf_namesgets(x, names);
  }
}
inline SEXP get_old_names(SEXP x){
  return attr::get_attr(x, symbol::names_sym);
}
inline bool has_r_names(SEXP x){
  SEXP names = SHIELD(get_old_names(x));
  bool out = !is_null(names);
  YIELD(1);
  return out;
}

inline SEXP get_old_class(SEXP x){
  return get_attr(x, symbol::class_sym);
}
inline void set_old_class(SEXP x, SEXP cls){
  Rf_classgets(x, cls);
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
    always_false<T>,
    "Unimplemented `new_vector` specialisation"
  );
  return T{};
}
// One-parameter template version
template <typename T>
inline SEXP new_vector(R_xlen_t n) {
  static_assert(
    always_false<T>,
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
  r_bool_t* RESTRICT p_out = internal::logical_ptr(out);
  r_fill(out, p_out, 0, n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<r_int_t>(R_xlen_t n){
  return internal::new_vec(INTSXP, n);
}
template <>
inline SEXP new_vector<r_int_t>(R_xlen_t n, const r_int_t default_value){
  SEXP out = SHIELD(new_vector<r_int_t>(n));
  r_int_t* RESTRICT p_out = internal::integer_ptr(out);
  r_fill(out, p_out, 0, n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<r_double_t>(R_xlen_t n){
  return internal::new_vec(REALSXP, n);
}
template <>
inline SEXP new_vector<r_double_t>(R_xlen_t n, const r_double_t default_value){
  SEXP out = SHIELD(new_vector<r_double_t>(n));
  r_double_t* RESTRICT p_out = internal::real_ptr(out);
  r_fill(out, p_out, 0, n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<int64_t>(R_xlen_t n){
  SEXP out = SHIELD(new_vector<r_double_t>(n));
  attr::set_old_class(out, SHIELD(internal::make_utf8_strsxp("integer64")));
  YIELD(2);
  return out;
}
template <>
inline SEXP new_vector<int64_t>(R_xlen_t n, const int64_t default_value){
  SEXP out = SHIELD(new_vector<int64_t>(n));
  int64_t* RESTRICT p_out = internal::integer64_ptr(out);
  r_fill(out, p_out, 0, n, default_value);
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
      set_value<r_string_t>(out, i, default_value);
    }
  }
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<r_complex_t>(R_xlen_t n){
  return internal::new_vec(CPLXSXP, n);
}
template <>
inline SEXP new_vector<r_complex_t>(R_xlen_t n, const r_complex_t default_value){
  SEXP out = SHIELD(new_vector<r_complex_t>(n));
  r_complex_t* RESTRICT p_out = internal::complex_ptr(out);
  r_fill(out, p_out, 0, n, default_value);
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<r_byte_t>(R_xlen_t n){
  return internal::new_vec(RAWSXP, n);
}
template <>
inline SEXP new_vector<r_byte_t>(R_xlen_t n, const r_byte_t default_value){
  SEXP out = SHIELD(new_vector<r_byte_t>(n));
  r_byte_t *p_out = internal::raw_ptr(out);
  r_fill(out, p_out, 0, n, default_value);
  YIELD(1);
  return out;
}
inline SEXP new_list(R_xlen_t n){
  return internal::new_vec(VECSXP, n);
}
inline SEXP new_list(R_xlen_t n, const SEXP default_value){
  SEXP out = SHIELD(internal::new_vec(VECSXP, n));
  if (!is_null(default_value)){
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, default_value);
    }
  }
  YIELD(1);
  return out;
}
template <>
inline SEXP new_vector<SEXP>(R_xlen_t n){
  return new_list(n);
}
template <>
inline SEXP new_vector<SEXP>(R_xlen_t n, const SEXP default_value){
  return new_list(n, default_value);
}

}

namespace df {

inline bool is_df(SEXP x){
  return vec::is_df(x);
}

inline r_int_t nrow(SEXP x){
  return Rf_length(attr::get_attr(x, symbol::row_names_sym));
}
inline r_int_t ncol(SEXP x){
  return Rf_length(x);
}
inline SEXP new_row_names(int n){
  if (n > 0){
    SEXP out = SHIELD(vec::new_vector<r_int_t>(2));
    vec::set_value(out, 0, na::integer);
    vec::set_value(out, 1, -n);
    YIELD(1);
    return out;
  } else {
    return vec::new_vector<r_int_t>(0);
  }
}
inline void set_row_names(SEXP x, int n){
  SEXP row_names = SHIELD(new_row_names(n));
  attr::set_attr(x, symbol::row_names_sym, row_names);
  YIELD(1);
}
}

namespace math {
template <typename T>
inline constexpr bool is_r_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_inf<r_double_t>(const r_double_t x){
  return x == r_limits::r_pos_inf || x == r_limits::r_neg_inf;
}

template <typename T>
inline constexpr bool is_r_pos_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_pos_inf<r_double_t>(const r_double_t x){
  return x == r_limits::r_pos_inf;
}

template <typename T>
inline constexpr bool is_r_neg_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_neg_inf<r_double_t>(const r_double_t x){
  return x == r_limits::r_neg_inf;
}

}

// C++ templates


template<typename T>
inline constexpr bool is_r_true(const T x) {
  return false;
}

template<>
inline constexpr bool is_r_true<r_bool_t>(const r_bool_t x){
  return x == r_true;
}

template<typename T>
inline constexpr bool is_r_false(const T x) {
  return false;
}

template<>
inline constexpr bool is_r_false<r_bool_t>(const r_bool_t x){
  return x == r_false;
}

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
inline constexpr bool is_r_na<r_int_t>(const r_int_t x){
  return x == na::integer;
}

template<>
inline constexpr bool is_r_na<r_double_t>(const r_double_t x){
  return x.value != x.value;
}

template<>
inline constexpr bool is_r_na<int64_t>(const int64_t x){
  return x == na::integer64;
}

template<>
inline constexpr bool is_r_na<r_complex_t>(const r_complex_t x){
  return is_r_na<r_double_t>(x.re()) || is_r_na<r_double_t>(x.im());
}

template<>
inline constexpr bool is_r_na<r_byte_t>(const r_byte_t x){
  return false;
}

template<>
inline constexpr bool is_r_na<r_symbol_t>(const r_symbol_t x){
  return false;
}

template<>
inline bool is_r_na<r_string_t>(const r_string_t x){
  return x == na::string;
}

// NULL is treated as NA of general R objects
template<>
inline bool is_r_na<SEXP>(const SEXP x){
  return is_null(x);
}

// NA type
namespace internal {

template<typename T>
inline constexpr T na_value_impl() {
  static_assert(
    always_false<T>,
    "Unimplemented `na_value` specialisation"
  );
  return T{};
}

template<>
inline constexpr r_bool_t na_value_impl<r_bool_t>(){
  return na::logical;
}

template<>
inline constexpr r_int_t na_value_impl<r_int_t>(){
  return na::integer;
}

template<>
inline r_double_t na_value_impl<r_double_t>(){
  return na::real;
}

template<>
inline constexpr int64_t na_value_impl<int64_t>(){
  return na::integer64;
}

template<>
inline r_complex_t na_value_impl<r_complex_t>(){
  return na::complex;
}

template<>
inline constexpr r_byte_t na_value_impl<r_byte_t>(){
  return r_byte_t{0};
}

template<>
inline r_string_t na_value_impl<r_string_t>(){
  return na::string;
}

template<>
inline SEXP na_value_impl<SEXP>(){
  return r_null;
}

}

template<typename T>
inline constexpr auto na_value() {
  return internal::na_value_impl<std::remove_cvref_t<T>>();
}


namespace vec {

template<typename T>
inline SEXP as_vector(const T x){
  if constexpr (is_r_or_cpp_integral_v<T>){
    return as_vector<int64_t>(x);
  } else if constexpr (std::is_convertible_v<T, SEXP>){
    return as_vector<SEXP>(x);
  } else {
    static_assert(
      always_false<T>,
      "Unimplemented `as_vector` specialisation"
    );
    return T{};
  }
}
template<>
inline SEXP as_vector<bool>(const bool x){
  return new_vector<r_bool_t>(1, static_cast<r_bool_t>(x));
}
template<>
inline SEXP as_vector<r_bool_t>(const r_bool_t x){
  return new_vector<r_bool_t>(1, x);
}
template<>
inline SEXP as_vector<Rboolean>(const Rboolean x){
  return new_vector<r_bool_t>(1, static_cast<r_bool_t>(x));
}
template<>
inline SEXP as_vector<r_int_t>(const r_int_t x){
  return new_vector<r_int_t>(1, x);
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
inline SEXP as_vector<r_double_t>(const r_double_t x){
  return new_vector<r_double_t>(1, x);
}
template<>
inline SEXP as_vector<r_complex_t>(const r_complex_t x){
  return new_vector<r_complex_t>(1, x);
}
template<>
inline SEXP as_vector<r_byte_t>(const r_byte_t x){
  return new_vector<r_byte_t>(1, x);
}
template<>
inline SEXP as_vector<const char *>(const char *x){
  return internal::make_utf8_strsxp(x);
}
template<>
inline SEXP as_vector<r_symbol_t>(const r_symbol_t x){
  return new_list(1, x);
}

// Scalar string
template<>
inline SEXP as_vector<r_string_t>(const r_string_t x){
  return new_vector<r_string_t>(1, x);
}
template<>
inline SEXP as_vector<SEXP>(const SEXP x){
  switch (TYPEOF(x)){
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

namespace internal {
template<typename T>
inline SEXP as_r_obj(const T x){
  if constexpr (std::is_convertible_v<T, SEXP>){
    return x;
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

// Assumes no NAs at all
template<typename T>
inline constexpr bool can_be_int(T x){
  // If x is an int type whose size is <= int OR
  // an arithmetic type (e.g. double)
  if constexpr (std::is_integral_v<T> && sizeof(T) <= sizeof(int)){
    return true;
  } else if constexpr (is_r_or_cpp_arithmetic_v<T>){
    return between<T>(x, r_limits::r_int_min, r_limits::r_int_max);
  } else {
    return false;
  }
}
template<typename T>
inline constexpr bool can_be_int64(T x){
  // If x is an int64 type whose size is <= int64 OR
  // an arithmetic type (e.g. double)
  if constexpr (is_r_or_cpp_integral_v<T> && sizeof(T) <= sizeof(int64_t)){
    return true;
  } else if constexpr (is_r_or_cpp_arithmetic_v<T>){
    return between<T>(x, r_limits::r_int64_min, r_limits::r_int64_max);
  } else {
    return false;
  }
}

// Coerce functions that account for NA
template<typename T>
inline constexpr r_bool_t as_bool(T x){
  if constexpr (is<T, int> || is<T, r_bool_t>){
    return static_cast<r_bool_t>(x);
  } else if constexpr (is_r_or_cpp_arithmetic_v<T>){
    return is_r_na(x) ? na::logical : static_cast<r_bool_t>(static_cast<bool>(x));
  } else {
    return na::logical;
  }
}
template<typename T>
inline constexpr r_int_t as_int(T x){
  if constexpr (is<T, int> || is<T, r_int_t>){
    return static_cast<r_int_t>(x);
  } else if constexpr (is_r_or_cpp_arithmetic_v<T>){
    return is_r_na(x) || !internal::can_be_int(x) ? na::integer : static_cast<r_int_t>(x);
  } else {
    return na::integer;
  }
}
template<typename T>
inline constexpr int64_t as_int64(T x){
  if constexpr (is<T, int64_t>){
    return x;
  } else if constexpr (is_r_or_cpp_arithmetic_v<T>){
    return is_r_na(x) || !internal::can_be_int64(x) ? na::integer64 : static_cast<int64_t>(x);
  } else {
    return na::integer64;
  }
}
template<typename T>
inline constexpr r_double_t as_double(T x){
  if constexpr (is<T, double> || is<T, r_double_t>){
    return static_cast<r_double_t>(x);
  } else if constexpr (is_r_or_cpp_arithmetic_v<T>){
    return is_r_na(x) ? na::real : static_cast<r_double_t>(x);
  } else {
    return na::real;
  }
}
template<typename T>
inline constexpr r_complex_t as_complex(T x){
  if constexpr (is<T, r_complex_t>){
    return x;
  } else if constexpr (is_r_or_cpp_arithmetic_v<T>){
    return r_complex_t{as_double(x), 0.0};
  } else {
    return na::complex;
  }
}
template<typename T>
inline constexpr r_byte_t as_raw(T x){
  if constexpr (is<T, r_byte_t>){
    return x;
  } else if constexpr (is_r_or_cpp_integral_v<T> && sizeof(T) <= sizeof(int8_t)){
    return is_r_na(x) || x < 0 ? na::raw : static_cast<r_byte_t>(x);
  } else if constexpr (std::is_convertible_v<T, r_byte_t>){
    return is_r_na(x) || !between(x, static_cast<T>(0), static_cast<T>(255)) ? na::raw : static_cast<r_byte_t>(x);
  } else {
    return na::raw;
  }
}
// As CHARSXP
template<typename T>
inline r_string_t as_r_string(T x){
  if constexpr (is<T, r_string_t>){
    return x;
  } else if constexpr (is<T, const char *>){
    return static_cast<r_string_t>(internal::make_utf8_charsxp(x));
  } else if constexpr (is<T, r_symbol_t>){
    return static_cast<r_string_t>(PRINTNAME(static_cast<SEXP>(x)));
  } else {
    SEXP scalar = SHIELD(as_r_obj(x));
    if (Rf_length(scalar) != 1){
      YIELD(1);
      Rf_error("`x` is a non-scalar vector and cannot be convered to an `r_string_t` in %s", __func__);
    }
    if (is_r_na(x)){
      YIELD(1);
      return na::string;
    }
    SEXP str = SHIELD(vec::coerce_vec(scalar, STRSXP));
    r_string_t out = vec::get_value<r_string_t>(str, 0);
    YIELD(2);
    return out;
  }
}

// As SYMSXP
template<typename T>
inline r_symbol_t as_r_sym(T x){
  if constexpr (is<T, r_symbol_t>){
    return x;
  } else if constexpr (is<T, const char *>){
    return static_cast<r_symbol_t>(Rf_install(x));
  } else if constexpr (is<T, r_string_t>){
    return as_r_sym(CHAR(static_cast<SEXP>(x)));
  } else {
    SEXP scalar = SHIELD(as_r_obj(x));
    if (Rf_length(scalar) != 1){
      YIELD(1);
      Rf_error("`x` is a non-scalar vector and cannot be convered to an `r_symbol_t` in %s", __func__);
    }
    SEXP str = SHIELD(vec::coerce_vec(scalar, STRSXP));
    r_symbol_t out = as_r_sym(vec::get_value<r_string_t>(str, 0));
    YIELD(2);
    return out;
  }
}

// R version of static_cast
template<typename T, typename U>
struct as_impl {
  static T cast(U x) {
    static_assert(
      always_false<T>,
      "Can't `as` this type, use `static_cast`"
    );
    return T{};
  }
};

// Specializations for each target type

template<typename U>
struct as_impl<r_bool_t, U> {
  static constexpr r_bool_t cast(U x) {
    return as_bool(x);
  }
};

template<typename U>
struct as_impl<r_int_t, U> {
  static constexpr r_int_t cast(U x) {
    return as_int(x);
  }
};

template<typename U>
struct as_impl<int64_t, U> {
  static constexpr int64_t cast(U x) {
    return as_int64(x);
  }
};

template<typename U>
struct as_impl<r_double_t, U> {
  static constexpr r_double_t cast(U x) {
    return as_double(x);
  }
};

template<typename U>
struct as_impl<r_complex_t, U> {
  static constexpr r_complex_t cast(U x) {
    return as_complex(x);
  }
};

template<typename U>
struct as_impl<r_byte_t, U> {
  static constexpr r_byte_t cast(U x) {
    return as_raw(x);
  }
};

template<typename U>
struct as_impl<r_string_t, U> {
  static r_string_t cast(U x) {
    return as_r_string(x);
  }
};

template<typename U>
struct as_impl<r_symbol_t, U> {
  static r_symbol_t cast(U x) {
    return as_r_sym(x);
  }
};

template<typename U>
struct as_impl<SEXP, U> {
  static constexpr SEXP cast(U x) {
    return as_r_obj(x);
  }
};
}

template<typename T, typename U>
inline constexpr T as(U x) {
  return internal::as_impl<std::remove_cvref_t<T>, U>::cast(x);
}

// Methods for custom R types

// r_complex_t methods

// Unary minus
inline constexpr r_complex_t operator-(const r_complex_t& x) {
  return r_complex_t{-x.re(), -x.im()};
}

// Binary arithmetic operators
inline constexpr r_complex_t operator+(const r_complex_t& lhs, const r_complex_t& rhs) {
  return r_complex_t{lhs.re() + rhs.re(), lhs.im() + rhs.im()};
}

inline constexpr r_complex_t operator-(const r_complex_t& lhs, const r_complex_t& rhs) {
  return r_complex_t{lhs.re() - rhs.re(), lhs.im() - rhs.im()};
}

inline constexpr r_complex_t operator*(const r_complex_t& lhs, const r_complex_t& rhs) {
  // (a+bi) * (c+di) = (ac-bd) + (ad+bc)i
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  return r_complex_t{a*c - b*d, a*d + b*c};
}

inline constexpr r_complex_t operator/(const r_complex_t& lhs, const r_complex_t& rhs) {
  // (a+bi) / (c+di) = [(ac+bd)/(c²+d²) + (bc-ad)/(c²+d²)i]
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  r_double_t denom = c*c + d*d;
  return r_complex_t{(a*c + b*d) / denom, (b*c - a*d) / denom};
}

// Compound assignment operators
inline constexpr r_complex_t& operator+=(r_complex_t& lhs, const r_complex_t& rhs) {
  lhs.value.r += rhs.value.r;
  lhs.value.i += rhs.value.i;
  return lhs;
}

inline constexpr r_complex_t& operator-=(r_complex_t& lhs, const r_complex_t& rhs) {
  lhs.value.r -= rhs.value.r;
  lhs.value.i -= rhs.value.i;
  return lhs;
}

inline constexpr r_complex_t& operator*=(r_complex_t& lhs, const r_complex_t& rhs) {
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  lhs.re() = a*c - b*d;
  lhs.im() = a*d + b*c;
  return lhs;
}

inline constexpr r_complex_t& operator/=(r_complex_t& lhs, const r_complex_t& rhs) {
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  r_double_t denom = c*c + d*d;
  lhs.re() = (a*c + b*d) / denom;
  lhs.im() = (b*c - a*d) / denom;
  return lhs;
}

template<RType T, RType U>
inline constexpr r_bool_t operator==(const T lhs, const U rhs) {

  // Check for NA in either operand
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return r_na;
  }

  if constexpr (is<T, r_complex_t> && is<U, r_complex_t>){
    // Compare complex types by components
    return r_bool_t{lhs.re() == rhs.re() && lhs.im() == rhs.im()};
  } else {
    return r_bool_t{lhs.value == rhs.value};
  }
}

// Other comparison operators
template<RType T, RType U>
inline constexpr r_bool_t operator!=(const T lhs, const U rhs) {

  // Check for NA in either operand
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return r_na;
  }

  if constexpr (is<T, r_complex_t> && is<U, r_complex_t>){
    return r_bool_t{lhs.re() != rhs.re() || lhs.im() != rhs.im()};
  } else {
    return r_bool_t{lhs.value != rhs.value};
  }
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator<(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return r_na;
  }
  return r_bool_t{lhs.value < rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator<=(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return r_na;
  }
  return r_bool_t{lhs.value <= rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator>(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return r_na;
  }
  return r_bool_t{lhs.value > rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator>=(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return r_na;
  }
  return r_bool_t{lhs.value >= rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr T operator+=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value += rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_double_t operator+=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value += rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
inline constexpr T operator+(const T lhs, const U rhs) {
  auto res = lhs;
  res += rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T operator-=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value -= rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_double_t operator-=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value -= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
inline constexpr T operator-(const T lhs, const U rhs) {
  auto res = lhs;
  res -= rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T operator*=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value *= rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_double_t operator*=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value *= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
inline constexpr T operator*(const T lhs, const U rhs) {
  auto res = lhs;
  res *= rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T operator/=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value /= rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_double_t operator/=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value /= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
inline constexpr T operator/(const T lhs, const U rhs) {
  auto res = lhs;
  res /= rhs;
  return res;
}

template<RMathType T>
inline constexpr T operator-(const T x) {
  if (is_r_na(x)){
    return x;
  }
  return T{-x.value};
}
template<>
inline constexpr r_double_t operator-(const r_double_t x) {
  return r_double_t{-x.value};
}

// R math fns
namespace internal {

inline r_double_t round_to_even(r_double_t x){
  return x - r_double_t{std::remainder(x.value, 1.0)};
}

}

namespace math {

template<RMathType T>
inline constexpr T r_abs(T x){
  if constexpr (is_r_or_cpp_integral_v<T>){
    return is_r_na(x) ? x : static_cast<T>(std::abs(x));
  } else {
    return static_cast<T>(std::abs(x));
  }
}

inline r_double_t r_abs(r_complex_t x){
  return std::sqrt(x.re() * x.re() + x.im() * x.im());
}

template<RMathType T>
inline T r_floor(T x){
  return is_r_na(x) ? x : std::floor(x);
}
template<>
inline r_double_t r_floor(r_double_t x){
  return std::floor(x);
}

template<typename T>
inline T r_ceiling(T x){
  return is_r_na(x) ? x : std::ceil(x);
}
template<>
inline r_double_t r_ceiling(r_double_t x){
  return std::ceil(x);
}

template<typename T>
inline T r_trunc(T x){
  return is_r_na(x) ? x : std::trunc(x);
}

template <>
inline r_double_t r_trunc(r_double_t x){
  return std::trunc(x) + 0.0;
}

template <typename T>
inline r_int_t r_sign(T x) {
  return is_r_na(x) ? na::integer : (T(0) < x) - (x < T(0));
}

template<MathType T>
inline T r_negate(T x){
  return -x;
}

template<MathType T>
inline r_double_t r_sqrt(T x){
  return std::sqrt(as<r_double_t>(x));
}

template<typename T, typename U>
  requires (is_r_or_cpp_arithmetic_v<T> || is_r_or_cpp_arithmetic_v<U>)
inline r_double_t r_pow(T x, U y){
  if (y == 0) return 1.0;
  if (x == 1) return 1.0;
  if (y == 2){
    r_double_t left = as<r_double_t>(x);
    return left * left;
  }
  return std::pow(as<r_double_t>(x), as<r_double_t>(y));
}

template<typename T>
inline r_double_t r_log10(T x){
  return std::log10(as<r_double_t>(x));
}

template<typename T>
inline r_double_t r_exp(T x){
  return std::exp(as<r_double_t>(x));
}

template<typename T, typename U>
inline r_double_t r_log(T x, U base){
  return std::log(as<r_double_t>(x)) / std::log(as<r_double_t>(base));
}
template<typename T>
inline r_double_t r_log(T x){
  return std::log(as<r_double_t>(x));
}
inline r_complex_t r_log(r_complex_t x){
  if (is_r_na(x)){ 
    return x;
  }
  r_double_t real = 0.5 * (r_log(r_pow(x.re(), 2.0) + r_pow(x.im(), 2.0)));
  r_double_t imag = std::atan2(as<r_double_t>(x.im()), as<r_double_t>(x.re()));
  return r_complex_t{real, imag};
}


template<typename T, typename U>
inline r_double_t r_round(T x, const U digits){
  if (is_r_na(x)){
    return as<r_double_t>(x);
  } else if (is_r_na(digits)){
    return na::real;
  } else if (is_r_inf(x)){
    return x;
  } else if (is_r_neg_inf(digits)){
    return 0.0;
  } else if (is_r_pos_inf(digits)){
    return x;
  } else {
    r_double_t scale = std::pow(10.0, as<r_double_t>(digits));
    return internal::round_to_even(as<r_double_t>(x) * scale) / scale;
  }
}

template<typename T>
inline r_double_t r_round(T x){
  if (is_r_na(x)){
    return as<r_double_t>(x);
  } else if (is_r_inf(x)){
    return x;
  } else {
    return internal::round_to_even(as<r_double_t>(x));
  }
}

template<typename T, typename U>
inline r_double_t r_signif(T x, const U digits){
  U new_digits = std::max(U(1), digits);
  if (is_r_na(x)){
    return as<r_double_t>(x);
  } else if (is_r_na(new_digits)){
    return na::real;
  } else if (is_r_pos_inf(digits)){
    return x;
  } else {
    new_digits -= std::ceil(std::log10(std::abs(x)));
    r_double_t scale = std::pow(10, new_digits);
    return internal::round_to_even(scale * x) / scale;
  }
}

template<typename T>
inline T r_abs_diff(const T x, const T y){
  return r_abs(x - y);
}

inline r_bool_t is_whole_number(const r_double_t x, const r_double_t tolerance){
  return is_r_na(x) || is_r_na(tolerance) ? na::logical : static_cast<r_bool_t>(r_abs_diff(x, std::round(x)) <= tolerance);
}


// Greatest common divisor
template<
  typename T,
  typename = typename std::enable_if<is_r_or_cpp_arithmetic_v<T>>::type
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
        return na_value<decltype(x)>();
      }
    }

    T ax = std::abs(x);
    T ay = std::abs(y);

    if constexpr (is_r_or_cpp_integral_v<T>){

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


// Lowest common multiple
template<typename T,
         typename = typename std::enable_if<is_r_or_cpp_arithmetic_v<T>>::type>
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
        return na_value<decltype(x)>();
      }
    }

    if constexpr (is_r_or_cpp_integral_v<T>){
      if (x == 0 && y == 0){
        return 0;
      }
      T res = std::abs(x) / r_gcd(x, y, na_rm);
      if (y != 0 && (std::abs(res) > (std::numeric_limits<T>::max() / std::abs(y)))){
        return na_value<decltype(x)>();
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
    expr = SHIELD(Rf_lang3(symbol::triple_colon_sym, as<r_symbol_t>(pkg), as<r_symbol_t>(name)));
  } else {
    expr = SHIELD(Rf_lang3(symbol::double_colon_sym, as<r_symbol_t>(pkg), as<r_symbol_t>(name)));
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
    return arg(name, as<SEXP>(v));
  }
};


namespace vec {

// Variadic list constructor
template<typename... Args>
inline SEXP make_list(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return new_list(0);
  } else {
    SEXP out = SHIELD(vec::new_list(n));

    // Are any args named?
    constexpr bool any_named = (is<Args, arg> || ...);

    SEXP nms = r_null;

    if (any_named){
      nms = vec::new_vector<r_string_t>(n);
    }
    SHIELD(nms);

    int i = 0;
    (([&]() {
      if constexpr (is<Args, arg>) {
        SET_VECTOR_ELT(out, i, args.value);
        set_value<r_string_t>(nms, i, as<r_string_t>(args.name));
      } else {
        SET_VECTOR_ELT(out, i, as<SEXP>(args));
      }
      ++i;
    }()), ...);

    attr::set_old_names(out, nms);
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
      if constexpr (is<Args, arg>) {
        SETCAR(current, args.value);
        SET_TAG(current, as<r_symbol_t>(args.name));
      } else {
        SETCAR(current, as<SEXP>(args));
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

inline R_xlen_t old_length(SEXP x){
  return Rf_xlength(x);
}

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
      const SEXP *p_x = list_ptr_ro(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
    } else {
      if (is_null(internal::BASE_LENGTH)){
        internal::BASE_LENGTH = as<r_symbol_t>("length");
      }
      SEXP expr = SHIELD(Rf_lang2(internal::BASE_LENGTH, x));
      SEXP r_len = SHIELD(eval(expr, env::base_env));
      R_xlen_t out = TYPEOF(r_len) == INTSXP ? INTEGER_ELT(r_len, 0) : REAL_ELT(r_len, 0);
      YIELD(2);
      return out;
    }
    // Catch-all
  } else {
    if (is_null(internal::BASE_LENGTH)){
      internal::BASE_LENGTH = as<r_symbol_t>("length");
    }
    SEXP expr = SHIELD(Rf_lang2(internal::BASE_LENGTH, x));
    SEXP r_len = SHIELD(eval(expr, env::base_env));
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
    return vec::new_vector<r_int_t>(0);
  }
  SEXP colon_fn = SHIELD(fn::find_pkg_fun(":", "base", false));
  SEXP out = SHIELD(fn::eval_fn(colon_fn, env::base_env, 1, n));
  YIELD(2);
  return out;
}

// r_bool_t not bool because bool can't be NA
inline r_bool_t all_whole_numbers(SEXP x, r_double_t tol_, bool na_rm_){

  R_xlen_t n = Rf_xlength(x);

  // Use r_bool_t instead of bool as r_bool_t can hold NA
  r_bool_t out = r_true;
  R_xlen_t na_count = 0;

  switch ( internal::CHEAPR_TYPEOF(x) ){
  case LGLSXP:
  case INTSXP:
  case internal::CHEAPR_INT64SXP: {
    break;
  }
  case REALSXP: {
    const r_double_t *p_x = internal::real_ptr_ro(x);
    for (R_xlen_t i = 0; i < n; ++i) {
      out = static_cast<r_bool_t>(math::is_whole_number(p_x[i], tol_));
      na_count += is_r_na(out);
      if (is_r_false(out)){
        break;
      }
    }
    if (is_r_true(out) && !na_rm_ && na_count > 0){
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

inline void add_attrs(SEXP x, SEXP attrs) {

  if (is_null(x)){
    Rf_error("Cannot add attributes to `NULL`");
  }

  int32_t NP = 0;

  switch (TYPEOF(attrs)){
  case NILSXP: {
    break;
  }
  case VECSXP: {
    SEXP names = SHIELD(attr::get_old_names(attrs)); ++NP;
    if (is_null(names)){
      YIELD(NP);
      Rf_error("attributes must be a named list");
    }
    const SEXP *p_attributes = list_ptr_ro(attrs);
    const r_string_t *p_names = vector_ptr<const r_string_t>(names);

    r_symbol_t attr_nm;

    for (int i = 0; i < Rf_length(names); ++i){
      if (p_names[i] != blank_r_string){
        attr_nm = as<r_symbol_t>(p_names[i]);
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
    r_string_t addr_x = SHIELD(address(x)); ++NP;

    SEXP current = attrs;

    while (!is_null(current)){
      if (is_null(symbol::tag(current)) || as<r_string_t>(symbol::tag(current)) == blank_r_string){
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

inline void clear_attrs(SEXP x){
  SEXP attrs = SHIELD(get_attrs(x));
  if (is_null(attrs)){
    YIELD(1);
    return;
  }
  SEXP names = SHIELD(attr::get_old_names(attrs));
  const r_string_t *p_names = vector_ptr<const r_string_t>(names);

  int n = Rf_length(attrs);
  for (R_xlen_t i = 0; i < n; ++i){
    r_symbol_t target_sym = as<r_symbol_t>(p_names[i]);
    set_attr(x, target_sym, r_null);
  }
  YIELD(2);
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
  int32_t NP = 0;
  SEXP out = r_null;
  R_xlen_t n = Rf_xlength(x);
  SEXP attrs = r_null;

  switch (TYPEOF(x)){
  case NILSXP: {
    break;
  }
  case LGLSXP: {
    using r_t = r_bool_t;
    out = SHIELD(new_vector<r_t>(n)); ++NP;
    r_copy_n(out, vector_ptr<r_t>(out), vector_ptr<const r_t>(x), 0, n);
    break;
  }
  case INTSXP: {
    using r_t = r_int_t;
    out = SHIELD(new_vector<r_t>(n)); ++NP;
    r_copy_n(out, vector_ptr<r_t>(out), vector_ptr<const r_t>(x), 0, n);
    break;
  }
  case REALSXP: {
    using r_t = r_double_t;
    out = SHIELD(new_vector<r_t>(n)); ++NP;
    r_copy_n(out, vector_ptr<r_t>(out), vector_ptr<const r_t>(x), 0, n);
    break;
  }
  case STRSXP: {
    using r_t = r_string_t;
    out = SHIELD(new_vector<r_t>(n)); ++NP;
    r_copy_n(out, vector_ptr<const r_t>(out), vector_ptr<const r_t>(x), 0, n);
    break;
  }
  case CPLXSXP: {
    using r_t = r_complex_t;
    out = SHIELD(new_vector<r_t>(n)); ++NP;
    r_copy_n(out, vector_ptr<r_t>(out), vector_ptr<const r_t>(x), 0, n);
    break;
  }
  case RAWSXP: {
    using r_t = r_byte_t;
    out = SHIELD(new_vector<r_t>(n)); ++NP;
    r_copy_n(out, vector_ptr<r_t>(out), vector_ptr<const r_t>(x), 0, n);
    break;
  }
  case VECSXP: {
    using r_t = SEXP;
    out = SHIELD(new_vector<r_t>(n)); ++NP;
    const r_t *p_x = vector_ptr<const r_t>(x);
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, deep_copy(p_x[i]));
    }
    break;
  }
  default: {
    out = SHIELD(Rf_duplicate(x)); ++NP;
    YIELD(NP);
    return out;
  }
  }

  if (!is_null(x)){
    SHIELD(attrs = attr::get_attrs(x)); ++NP;
    int n_attrs = Rf_length(attrs);
    for (R_xlen_t i = 0; i < n_attrs; ++i){
      SET_VECTOR_ELT(attrs, i, deep_copy(VECTOR_ELT(attrs, i)));
    }
    attr::set_attrs(out, attrs);
  }

  YIELD(NP);
  return out;
}

}


namespace env {
inline SEXP get(r_symbol_t sym, SEXP env, bool inherits = true){

  if (TYPEOF(env) != ENVSXP){
    Rf_error("second argument to '%s' must be an environment", __func__);
  }

  SEXP val = inherits ? Rf_findVar(sym, env) : Rf_findVarInFrame(env, sym);

  if (val == static_cast<SEXP>(symbol::missing_arg)){
    Rf_error("arg `sym` cannot be missing");
  } else if (val == static_cast<SEXP>(symbol::unbound_value)){
    return r_null;
  } else if (TYPEOF(val) == PROMSXP){
    SHIELD(val);
    val = eval(val, env);
    YIELD(1);
  }
  return val;
}
}

// We call R fn`cheapr::set_threads` to make sure the R option is set
inline void set_threads(uint16_t n){
  uint16_t max_threads = OMP_MAX_THREADS;
  uint16_t threads = std::min(n, max_threads);
  SEXP cheapr_set_threads = SHIELD(fn::find_pkg_fun("set_threads", "cheapr", true));
  SEXP r_threads = SHIELD(vec::as_vector(as<r_int_t>(threads)));
  SHIELD(fn::eval_fn(cheapr_set_threads, R_BaseEnv, r_threads));
  YIELD(3);
}


namespace internal {

// A cleaner lambda-based alternative to
// using the canonical switch(TYPEOF(x))
//
// Pass both the SEXP and an auto variable inside a lambda
// and visit_vector() will assign the auto variable to the
// correct pointer
// Then simply deduce its type (via decltype) for further manipulation
// To be used in a lambda
// E.g. visit_r_ptr(x, [&](auto p_x) {})

// One must account for the default case via
// if constexpr (std::is_same_v<T, std::nullptr_t>)
// Since `NULL` is included in the default case, if you want
// separate logic to handle this case, just do the below inside the default case
// if (is_null(x)){
// ...
// } else {
// ...
// }
template <class F>
decltype(auto) visit_vector(SEXP x, F&& f) {
  switch (CHEAPR_TYPEOF(x)) {
  case LGLSXP:          return f(vector_ptr<const r_bool_t>(x));
  case INTSXP:          return f(vector_ptr<const r_int_t>(x));
  case CHEAPR_INT64SXP: return f(vector_ptr<const int64_t>(x));
  case REALSXP:         return f(vector_ptr<const r_double_t>(x));
  case STRSXP:          return f(vector_ptr<const r_string_t>(x));
  case VECSXP:          return f(vector_ptr<const SEXP>(x));
  case CPLXSXP:         return f(vector_ptr<const r_complex_t>(x));
  case RAWSXP:          return f(vector_ptr<const r_byte_t>(x));
  default:              return f(nullptr);
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
