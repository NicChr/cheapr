#ifndef cheapr_h
#define cheapr_h

#include <cpp11.hpp>
#include <Rinternals.h>

#ifdef _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT __restrict__
#endif

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#ifndef VECTOR_PTR_RO
#define VECTOR_PTR_RO(x) ((const SEXP*) DATAPTR_RO(x))
#endif
#ifndef INTEGER64_PTR
#define INTEGER64_PTR(x) ((int64_t*) REAL(x))
#endif
#ifndef INTEGER64_RO_PTR
#define INTEGER64_RO_PTR(x) ((int64_t*) REAL_RO(x))
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

template<typename T>
inline constexpr bool between(T x, T lo, T hi) {
  return x >= lo && x <= hi;
}

template<typename T>
inline constexpr bool is_integerable(T x){
  return between<T>(x, INTEGER_MIN, INTEGER_MAX);
}


int num_cores();
SEXP cpp_which_(SEXP x, bool invert);
SEXP cpp_missing_row(SEXP x, double threshold, bool threshold_is_prop);
int int_div(int x, int y);
SEXP xlen_to_r(R_xlen_t x);
R_xlen_t vec_length(SEXP x);
SEXP r_address(SEXP x);
R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive);
SEXP cpp_list_as_df(SEXP x);
SEXP cpp_new_df(SEXP x, SEXP nrows, bool recycle, bool name_repair);
SEXP cpp_is_na(SEXP x);
SEXP cpp_which_na(SEXP x);
SEXP cpp_which_not_na(SEXP x);
SEXP check_transform_altrep(SEXP x);
SEXP altrep_materialise(SEXP x);
SEXP compact_seq_data(SEXP x);
bool is_compact_seq(SEXP x);
R_xlen_t na_count(SEXP x, bool recursive);
bool cpp_any_na(SEXP x, bool recursive);
SEXP cpp_int64_to_double(SEXP x);
SEXP cpp_numeric_to_int64(SEXP x);
SEXP cpp_int64_to_numeric(SEXP x);
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add);
SEXP cpp_set_rm_attributes(SEXP x);
SEXP coerce_vector(SEXP source, SEXPTYPE type);
bool implicit_na_coercion(SEXP x, SEXP target);
SEXP cpp_val_find(SEXP x, SEXP value, bool invert, SEXP n_values);
SEXP cpp_set_divide(SEXP x, SEXP y);
SEXP cpp_val_remove(SEXP x, SEXP value);
SEXP cpp_seq_len(R_xlen_t n);
SEXP create_df_row_names(int n);
SEXP cpp_shallow_copy(SEXP x);
SEXP exclude_locs(SEXP exclude, R_xlen_t xn);
R_xlen_t unnested_length(SEXP x);
SEXP cpp_drop_null(SEXP l, bool always_shallow_copy);
SEXP cpp_lengths(SEXP x, bool names);
SEXP sset_vec(SEXP x, SEXP indices, bool check);
SEXP cpp_sset(SEXP x, SEXP indices, bool check);
SEXP cpp_df_slice(SEXP x, SEXP indices, bool check);
SEXP cpp_df_select(SEXP x, SEXP locs);
SEXP cpp_df_subset(SEXP x, SEXP i, SEXP j, bool check);
SEXP cpp_which_val(SEXP x, SEXP value, bool invert);
SEXP cpp_sequence(SEXP size, SEXP from, SEXP by);
SEXP cpp_rep_len(SEXP x, int length);
SEXP cpp_rep(SEXP x, SEXP times);
SEXP cpp_recycle(SEXP x, SEXP length);
SEXP cpp_c(SEXP x);
SEXP cpp_list_c(SEXP x);
SEXP cpp_loc_set_replace(SEXP x, SEXP where, SEXP what);
SEXP cpp_name_repair(SEXP names, SEXP dup_sep, SEXP empty_sep);
SEXP cpp_unique(SEXP x, bool names);
SEXP cpp_setdiff(SEXP x, SEXP y, bool unique);
SEXP cpp_intersect(SEXP x, SEXP y, bool unique);
SEXP get_ptype(SEXP x);
SEXP get_list_element(SEXP list, SEXP str);
SEXP list_c2(SEXP x, SEXP y);
SEXP c2(SEXP x, SEXP y);
SEXP rebuild(SEXP x, SEXP source, bool shallow_copy);
SEXP cpp_df_assign_cols(SEXP x, SEXP cols);
SEXP cpp_df_col_c(SEXP x, bool recycle, bool name_repair);
SEXP cpp_list_assign(SEXP x, SEXP values);
SEXP slice_loc(SEXP x, R_xlen_t i);
double cpp_sum(SEXP x);
double cpp_min(SEXP x);
SEXP cpp_str_coalesce(SEXP x);
SEXP cpp_na_init(SEXP x, int n);
SEXP new_list(R_xlen_t length, SEXP default_value);
void set_list_as_df(SEXP x);
SEXP cpp_semi_copy(SEXP x);
void clear_attributes(SEXP x);
uint_fast64_t null_count(SEXP x);

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

inline bool is_null(SEXP x){
  return x == R_NilValue;
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

inline R_xlen_t r_length(SEXP x){
  return Rf_asReal(cpp11::package("base")["length"](x));
}

inline int df_nrow(SEXP x){
  return Rf_length(Rf_getAttrib(x, R_RowNamesSymbol));
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
  SEXP out = SHIELD(new_vec(LGLSXP, 1));
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
  const SEXP *p_x = STRING_PTR_RO(xclass);

  return Rf_length(xclass) == 3 &&
    std::strcmp(CHAR(p_x[0]), "tbl_df") == 0 &&
    std::strcmp(CHAR(p_x[1]), "tbl") == 0 &&
    std::strcmp(CHAR(p_x[2]), "data.frame") == 0;
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

inline cpp11::function cheapr_sset = cpp11::package("cheapr")["sset"];
inline cpp11::function base_sset = cpp11::package("base")["["];
inline cpp11::function cheapr_is_na = cpp11::package("cheapr")["is_na"];
inline cpp11::function cheapr_factor = cpp11::package("cheapr")["factor_"];
inline cpp11::function base_colon = cpp11::package("base")[":"];
inline cpp11::function base_rep = cpp11::package("base")["rep"];
inline cpp11::function base_do_call = cpp11::package("base")["do.call"];
inline cpp11::function base_as_character = cpp11::package("base")["as.character"];
inline cpp11::function base_paste0 = cpp11::package("base")["paste0"];
inline cpp11::function cheapr_fast_match = cpp11::package("cheapr")["fast_match"];
inline cpp11::function cheapr_fast_unique = cpp11::package("cheapr")["fast_unique"];
inline cpp11::function cheapr_rebuild = cpp11::package("cheapr")["rebuild"];

inline bool address_equal(SEXP x, SEXP y){
  return r_address(x) == r_address(y);
}

#endif
