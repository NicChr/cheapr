#ifndef cheapr_cpp_funs
#define cheapr_cpp_funs

#include <cpp11.hpp>
#include <Rinternals.h>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#ifndef VECTOR_PTR_RO
#define VECTOR_PTR_RO(x) ((const SEXP*) DATAPTR_RO(x))
#endif
#ifndef INTEGER64_PTR
#define INTEGER64_PTR(x) ((long long*) REAL(x))
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

#ifndef integer_max_
#define integer_max_ std::numeric_limits<int>::max()
#endif

#ifndef cheapr_is_na_int
#define cheapr_is_na_int(x) ((bool) (x == NA_INTEGER))
#endif

#ifndef cheapr_is_na_str
#define cheapr_is_na_str(x) ((bool) (x == NA_STRING))
#endif

#ifndef cheapr_is_na_cplx
#define cheapr_is_na_cplx(x) ((bool) (x.r != x.r) || (x.i != x.i))
#endif

#ifndef cheapr_is_na_dbl
#define cheapr_is_na_dbl(x) ((bool) (x != x))
#endif

#ifndef NA_INTEGER64
#define NA_INTEGER64 LLONG_MIN
#endif

#ifndef cheapr_is_na_int64
#define cheapr_is_na_int64(x) ((bool) (x == NA_INTEGER64))
#endif

#ifndef CHEAPR_INT_TO_INT64
#define CHEAPR_INT_TO_INT64(x) ((long long int) (x == NA_INTEGER ? NA_INTEGER64 : x))
#endif
#ifndef CHEAPR_DBL_TO_INT64
#define CHEAPR_DBL_TO_INT64(x) ((long long int) (x != x ? NA_INTEGER64 : x))
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
#define CHEAPR_TYPEOF(x)  ( (SEXPTYPE) (Rf_inherits(x, "integer64") ? CHEAPR_INT64SXP : TYPEOF(x)) )
#endif

inline bool is_int64(SEXP x){
  return Rf_isReal(x) && Rf_inherits(x, "integer64");
}

inline bool is_df(SEXP x){
  return Rf_inherits(x, "data.frame");
}

inline R_xlen_t r_length(SEXP x){
  return Rf_asReal(cpp11::package("base")["length"](x));
}

inline cpp11::function base_match = cpp11::package("base")["match"];
inline cpp11::function cheapr_sset = cpp11::package("cheapr")["sset"];
inline cpp11::function cheapr_is_na = cpp11::package("cheapr")["is_na"];
inline cpp11::function base_colon = cpp11::package("base")[":"];


// Definition of simple atomic vector is one in which
// it is both atomic and all attributes data-independent
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

int num_cores();
SEXP cpp_which_(SEXP x, bool invert);
SEXP cpp_missing_row(SEXP x, double threshold, bool threshold_is_prop);
int int_div(int x, int y);
R_xlen_t cpp_df_nrow(SEXP x);
SEXP xlen_to_r(R_xlen_t x);
R_xlen_t vec_length(SEXP x);
SEXP r_address(SEXP x);
R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive);
SEXP cpp_list_as_df(SEXP x);
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
void cpp_copy_names(SEXP source, SEXP target, bool deep_copy);
void cpp_copy_attributes(SEXP source, SEXP target, bool deep_copy);
SEXP coerce_vector(SEXP source, SEXPTYPE type);
bool implicit_na_coercion(SEXP x, SEXP target);
SEXP cpp_val_find(SEXP x, SEXP value, bool invert, SEXP n_values);
double round_nearest_even(double x);
SEXP cpp_set_divide(SEXP x, SEXP y);
SEXP cpp_val_remove(SEXP x, SEXP value);
SEXP cpp_seq_len(R_xlen_t n);
SEXP create_df_row_names(int n);
SEXP shallow_copy(SEXP x);
SEXP exclude_locs(SEXP exclude, R_xlen_t xn);

#endif
