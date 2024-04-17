#include <cpp11.hpp>
#include <Rinternals.h>

#ifndef cheapr_cpp_funs
#define cheapr_cpp_funs

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#define VECTOR_PTR(x) ((SEXP *) DATAPTR(x))
#define VECTOR_PTR_RO(x) ((const SEXP*) DATAPTR_RO(x))

#ifdef _OPENMP
#include <omp.h>
#define OMP_NUM_PROCS omp_get_num_procs()
#define OMP_THREAD_LIMIT omp_get_thread_limit()
#define OMP_MAX_THREADS omp_get_max_threads()
#define OMP_PARALLEL _Pragma("omp parallel num_threads(n_cores) ")
#define OMP_FOR_SIMD _Pragma("omp for simd ")
#define OMP_PARALLEL_FOR_SIMD	_Pragma("omp parallel for simd num_threads(n_cores) ")
// #if _OPENMP >= 201307
//   #define OMP_VER_4
// #endif
#else
#define OMP_NUM_PROCS 1
#define OMP_THREAD_LIMIT 1
#define OMP_MAX_THREADS 1
#define OMP_PARALLEL
#define OMP_FOR_SIMD
#define OMP_PARALLEL_FOR_SIMD
#endif

#define integer_max_ std::numeric_limits<int>::max()

// #ifdef OMP_VER_4
// #define SAFE_SIMD _Pragma("omp simd")
// #define SAFE_FOR_SIMD _Pragma("omp for simd")
// #define SAFE_PARALLEL_FOR_SIMD _Pragma("omp parallel for simd")
// #else
// #define SAFE_SIMD
// #define SAFE_FOR_SIMD
// #define SAFE_PARALLEL_FOR_SIMD
// #endif

int num_cores();
SEXP cpp_which_(SEXP x, bool invert);
SEXP cpp_missing_row(SEXP x, double threshold, bool threshold_is_prop);
int int_div(int x, int y);
R_xlen_t cpp_df_nrow(SEXP x);
R_xlen_t cpp_unnested_length(SEXP x);
SEXP xlen_to_r(R_xlen_t x);
R_xlen_t cpp_vec_length(SEXP x);
SEXP r_address(SEXP x);
R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive);
SEXP cpp_list_as_df(SEXP x);
SEXP cpp_is_na(SEXP x);
SEXP cpp_which_na(SEXP x);
SEXP cpp_which_not_na(SEXP x);

#endif
