#include <cpp11.hpp>
#include <Rinternals.h>

#ifndef cheapr_cpp_funs
#define cheapr_cpp_funs

#define R_NO_REMAP
#define VECTOR_PTR_RO(x) ((const SEXP*) DATAPTR_RO(x))

#ifdef _OPENMP
#include <omp.h>
#define OMP_NUM_PROCS omp_get_num_procs()
#define OMP_THREAD_LIMIT omp_get_thread_limit()
#define OMP_MAX_THREADS omp_get_max_threads()
// #if _OPENMP >= 201307
//   #define OMP_VER_4
// #endif
#else
#define OMP_NUM_PROCS 1
#define OMP_THREAD_LIMIT 1
#define OMP_MAX_THREADS 1
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
R_xlen_t count_true(int *px, R_xlen_t n);
SEXP cpp_missing_row(SEXP x, double threshold, bool threshold_is_prop);
int int_div(int x, int y);
R_xlen_t cpp_df_nrow(SEXP x);
R_xlen_t cpp_unnested_length(SEXP x);
SEXP xlen_to_r(R_xlen_t x);
R_xlen_t cpp_vec_length(SEXP x);

#endif
