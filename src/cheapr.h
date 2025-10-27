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

// Make function definitions visible to all C++ files

int num_cores();
SEXP cpp_which_(SEXP x, bool invert);
SEXP cpp_missing_row(SEXP x, double threshold, bool threshold_is_prop);
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
SEXP cpp_sequence(SEXP size, SEXP from, SEXP by, bool as_list, bool add_id);
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
SEXP compact_seq_len(R_xlen_t n);
SEXP clean_indices(SEXP indices, SEXP x, bool count);
SEXP cpp_combine_levels(SEXP x);
SEXP fast_cast(SEXP x);
SEXP cpp_lgl_count(SEXP x);
SEXP cpp_lgl_locs(SEXP x, R_xlen_t n_true, R_xlen_t n_false,
                  bool include_true, bool include_false, bool include_na);
SEXP cpp_cast_all(SEXP x);
SEXP factor_as_character(SEXP x);

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

inline cpp11::function cheapr_sset = cpp11::package("cheapr")["cheapr_sset"];
inline cpp11::function cheapr_is_na = cpp11::package("cheapr")["is_na"];
inline cpp11::function cheapr_factor = cpp11::package("cheapr")["factor_"];
inline cpp11::function base_rep = cpp11::package("base")["rep"];
inline cpp11::function base_do_call = cpp11::package("base")["do.call"];
inline cpp11::function base_as_character = cpp11::package("base")["as.character"];
inline cpp11::function base_paste0 = cpp11::package("base")["paste0"];
inline cpp11::function cheapr_fast_match = cpp11::package("cheapr")["fast_match"];
inline cpp11::function cheapr_fast_unique = cpp11::package("cheapr")["fast_unique"];
inline cpp11::function cheapr_rebuild = cpp11::package("cheapr")["rebuild"];
inline cpp11::function cheapr_as_df = cpp11::package("cheapr")["as_df"];
inline cpp11::function cheapr_cast = cpp11::package("cheapr")["cast"];

inline bool address_equal(SEXP x, SEXP y){
  return r_address(x) == r_address(y);
}

// Types

static SEXP as_lgl = NULL;
static SEXP as_int = NULL;
static SEXP as_dbl = NULL;
static SEXP as_char = NULL;
static SEXP as_cplx = NULL;
static SEXP as_raw = NULL;
static SEXP as_date = NULL;
static SEXP as_posixct = NULL;
static SEXP as_list = NULL;

// r types

struct r_null_t {};
struct r_logical_t {};
struct r_integer_t {};
struct r_integer64_t {};
struct r_numeric_t {};
struct r_character_t {};
struct r_complex_t {};
struct r_raw_t {};
struct r_list_t {};
struct r_factor_t {};
struct r_date_t {};
struct r_posixt_t {};
struct r_data_frame_t {};
struct r_unknown_t {};


// r type constants
using r_type = uint8_t;
enum : r_type {
  r_null = (uint8_t) 0,
    r_lgl = (uint8_t) 1,
    r_int = (uint8_t) 2,
    r_int64 = (uint8_t) 3,
    r_dbl = (uint8_t) 4,
    r_chr = (uint8_t) 5,
    r_cplx = (uint8_t) 6,
    r_raw = (uint8_t) 7,
    r_list = (uint8_t) 8,
    r_fct = (uint8_t) 9,
    r_date = (uint8_t) 10,
    r_pxct = (uint8_t) 11,
    r_df = (uint8_t) 12,
    r_unk = (uint8_t) 13,
};

// Symmetric lattice with columns/rows in the new order:
// NULL, LGL, INT, I64, DBL, CHR, CPLX, RAW, LIST, FCT, DATE, PXCT, DF
inline constexpr uint8_t r_type_pairs[14][14] = {
  /*            NULL    LGL     INT     I64     DBL     CHR     CPLX    RAW     LIST    FCT     DATE    PXCT    DF   Unknown */
  /* NULL */  { r_null, r_lgl,  r_int,  r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* LGL  */  { r_lgl,  r_lgl,  r_int,  r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* INT  */  { r_int,  r_int,  r_int,  r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* I64  */  { r_int64,r_int64,r_int64,r_int64, r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* DBL  */  { r_dbl,  r_dbl,  r_dbl,  r_dbl,  r_dbl,  r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* CHR  */  { r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_chr,  r_list, r_fct,  r_chr,  r_chr,  r_df, r_unk },
  /* CPLX */  { r_cplx, r_cplx, r_cplx, r_cplx, r_cplx, r_chr,  r_cplx, r_raw,  r_list, r_fct,  r_fct,  r_df,   r_df, r_unk },
  /* RAW  */  { r_raw,  r_raw,  r_raw,  r_raw,  r_raw,  r_chr,  r_raw,  r_raw,  r_list, r_df,   r_df,   r_df,   r_df, r_unk },
  /* LIST */  { r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_list, r_df, r_unk },
  /* FCT  */  { r_fct,  r_fct,  r_fct,  r_fct,  r_fct,  r_fct,  r_fct,  r_df,   r_list, r_fct,  r_fct,  r_fct,  r_df, r_unk },
  /* DATE */  { r_date, r_date, r_date, r_date, r_date, r_chr,  r_fct,  r_df,   r_list, r_fct,  r_date, r_pxct, r_df, r_unk },
  /* PXCT */  { r_pxct, r_pxct, r_pxct, r_pxct, r_pxct, r_chr,  r_df,   r_df,   r_list, r_fct,  r_pxct, r_pxct, r_df, r_unk },
  /* DF   */  { r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df,   r_df, r_unk },
  /* Unknown   */  { r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk,   r_unk, r_unk }
};

// cast template with specialisations
template<typename T>
inline SEXP cast(SEXP x, SEXP y) {
  Rf_error("Unimplemented type conversion");
}

template<>
inline SEXP cast<r_null_t>(SEXP x, SEXP y) {
  return R_NilValue;
}

template<>
inline SEXP cast<r_logical_t>(SEXP x, SEXP y) {
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
inline SEXP cast<r_integer_t>(SEXP x, SEXP y) {
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
inline SEXP cast<r_integer64_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "integer64")){
    return x;
  } else {
    return coerce_vector(x, CHEAPR_INT64SXP);
  }
}

template<>
inline SEXP cast<r_numeric_t>(SEXP x, SEXP y) {
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
inline SEXP cast<r_character_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "character")){
    return x;
  } else if (Rf_isFactor(x)){
    return factor_as_character(x);
  } else if (Rf_isObject(x)){
    as_char = as_char != NULL ? as_char : Rf_install("as.character");
    return Rf_eval(Rf_lang2(as_char, x), R_GetCurrentEnv());
  } else {
    return coerce_vec(x, STRSXP);
  }
}

template<>
inline SEXP cast<r_complex_t>(SEXP x, SEXP y) {
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
inline SEXP cast<r_raw_t>(SEXP x, SEXP y) {
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
inline SEXP cast<r_list_t>(SEXP x, SEXP y) {
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
inline SEXP cast<r_factor_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "factor") && !Rf_inherits(y, "factor")){
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

template<>
inline SEXP cast<r_date_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "Date")){
    return x;
  } else if (Rf_isObject(x)){
    as_date = as_date != NULL ? as_date : Rf_install("as.Date");
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
inline SEXP cast<r_posixt_t>(SEXP x, SEXP y) {
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
inline SEXP cast<r_data_frame_t>(SEXP x, SEXP y) {
  if (Rf_inherits(x, "data.frame") && !Rf_inherits(y, "data.frame")){
    return x;
  } else if (Rf_inherits(x, "data.frame") && Rf_inherits(y, "data.frame")){
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

template<>
inline SEXP cast<r_unknown_t>(SEXP x, SEXP y) {
  return cheapr_cast(x, y);
}

inline r_type common_type(const r_type &a, const r_type &b) {
  return r_type_pairs[a][b];
}

// Convert single SEXP into r_* code.
inline const r_type get_r_type(SEXP x) {

  if (!Rf_isObject(x)){
    switch (TYPEOF(x)) {
    case NILSXP:  return r_null;
    case LGLSXP:  return r_lgl;
    case INTSXP:  return r_int;
    case REALSXP: return r_dbl;
    case STRSXP:  return r_chr;
    case CPLXSXP: return r_cplx;
    case RAWSXP:  return r_raw;
    case VECSXP:  return r_list;
    default: return r_unk;
    }
  } else {

    if (Rf_inherits(x, "data.frame")) return r_df;
    if (Rf_inherits(x, "POSIXct"))    return r_pxct;
    if (Rf_inherits(x, "Date"))       return r_date;
    if (Rf_inherits(x, "factor"))     return r_fct;
    if (Rf_inherits(x, "integer64"))  return r_int64;
    return r_unk;
  }
}

// Wrapper functions for cast fns map
using cast_fn = SEXP(*)(SEXP, SEXP);
inline SEXP cast_null(SEXP x, SEXP y) { return cast<r_null_t>(x, y); }
inline SEXP cast_logical(SEXP x, SEXP y) { return cast<r_logical_t>(x, y); }
inline SEXP cast_integer(SEXP x, SEXP y) { return cast<r_integer_t>(x, y); }
inline SEXP cast_integer64(SEXP x, SEXP y) { return cast<r_integer64_t>(x, y); }
inline SEXP cast_numeric(SEXP x, SEXP y) { return cast<r_numeric_t>(x, y); }
inline SEXP cast_character(SEXP x, SEXP y) { return cast<r_character_t>(x, y); }
inline SEXP cast_complex(SEXP x, SEXP y) { return cast<r_complex_t>(x, y); }
inline SEXP cast_raw(SEXP x, SEXP y) { return cast<r_raw_t>(x, y); }
inline SEXP cast_list(SEXP x, SEXP y) { return cast<r_list_t>(x, y); }
inline SEXP cast_factor(SEXP x, SEXP y) { return cast<r_factor_t>(x, y); }
inline SEXP cast_date(SEXP x, SEXP y) { return cast<r_date_t>(x, y); }
inline SEXP cast_posixt(SEXP x, SEXP y) { return cast<r_posixt_t>(x, y); }
inline SEXP cast_data_frame(SEXP x, SEXP y) { return cast<r_data_frame_t>(x, y); }
inline SEXP cast_unknown(SEXP x, SEXP y) { return cast<r_unknown_t>(x, y); }

static const cast_fn CAST_FNS[14] = {
  cast_null, cast_logical, cast_integer, cast_integer64,
  cast_numeric, cast_character, cast_complex, cast_raw, cast_list,
  cast_factor, cast_date, cast_posixt, cast_data_frame, cast_unknown
};

// Dispatcher function
inline SEXP cast_(r_type cast_type, SEXP x, SEXP y) {
  return CAST_FNS[cast_type](x, y);
}

r_type r_common_type(SEXP x);

#endif
