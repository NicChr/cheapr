#include "cheapr_cpp.h"

// The below works but probably not a good idea to use internal R code
// Direct inclusion of internal structure definitions
// typedef struct sxpinfo_struct {
//   SEXPTYPE type      : 5;
//   unsigned int obj   : 1;
//   unsigned int named : 2;
//   unsigned int gp    : 16;
//   unsigned int mark  : 1;
//   unsigned int debug : 1;
//   unsigned int trace : 1;
//   unsigned int spare : 1;
//   unsigned int gcgen : 1;
//   unsigned int gccls : 3;
// } sxpinfo_struct;
//
// typedef struct SEXPREC {
//   sxpinfo_struct sxpinfo;
//   struct SEXPREC *attrib;
//   struct SEXPREC *gengc_next_node, *gengc_prev_node;
//   union {
//     struct primsxp_struct {int offset; } primsxp;
//     struct symsxp_struct {
//       struct SEXPREC *pname;
//       struct SEXPREC *value;
//       struct SEXPREC *internal;
//     } symsxp;
//     struct listsxp_struct {
//       struct SEXPREC *carval;
//       struct SEXPREC *cdrval;
//       struct SEXPREC *tagval;
//     } listsxp;
//     struct envsxp_struct {
//       struct SEXPREC *frame;
//       struct SEXPREC *enclos;
//       struct SEXPREC *hashtab;
//     } envsxp;
//     struct closxp_struct {
//       struct SEXPREC *formals;
//       struct SEXPREC *body;
//       struct SEXPREC *env;
//     } closxp;
//     struct promsxp_struct {
//       struct SEXPREC *value;
//       struct SEXPREC *expr;
//       struct SEXPREC *env;
//     } promsxp;
//   } u;
// } SEXPREC, *SEXP;
//
// #define SET_TYPEOF2(x, v) ((x)->sxpinfo.type = (v))


R_xlen_t na_count(SEXP x, bool recursive){
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t count = 0;
  int n_prot = 0;
  int n_cores = n >= 100000 ? num_cores() : 1;
  bool do_parallel = n_cores > 1;
  switch ( TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      for (R_xlen_t i = 0; i < n; ++i) count += (p_x[i] == NA_INTEGER);
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) count += (p_x[i] == NA_INTEGER);
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      for (R_xlen_t i = 0; i < n; ++i) count += (p_x[i] != p_x[i]);
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) count += (p_x[i] != p_x[i]);
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      for (R_xlen_t i = 0; i < n; ++i) count += (p_x[i] == NA_STRING);
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) count += (p_x[i] == NA_STRING);
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      for (R_xlen_t i = 0; i < n; ++i){
        count += is_na_cplx(p_x[i]);
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        count += is_na_cplx(p_x[i]);
      }
    }
    break;
  }
  case VECSXP: {
    // We use a recursive method if recursive is true
    // Otherwise we skip to the default section below
    if (recursive){
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      count += na_count(p_x[i], true);
    }
    break;
  }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x));
    ++n_prot;
    SEXP r_true = Rf_protect(Rf_ScalarLogical(true));
    ++n_prot;
    count = scalar_count(is_missing, r_true, true);
    break;
  }
  }
  Rf_unprotect(n_prot);
  return count;
}

[[cpp11::register]]
SEXP cpp_num_na(SEXP x, bool recursive){
  return xlen_to_r(na_count(x, recursive));
}


[[cpp11::register]]
bool cpp_any_na(SEXP x, bool recursive){
  int n_prot = 0;
  R_xlen_t n = Rf_xlength(x);
  bool out = false;
  switch ( TYPEOF(x) ){
  case NILSXP: {
    return out;
  }
  case LGLSXP:
  case INTSXP: {
    // The commented-out code is a way to
    // do a pseudo while-loop in openMP
    // volatile bool flag = false;
    int *p_x = INTEGER(x);
    // #pragma omp parallel for num_threads(num_cores()) shared(flag)
    for (R_xlen_t i = 0; i < n; ++i){
      // if (flag) continue;
      if (p_x[i] == NA_INTEGER){
        out = true;
        // flag = true;
        break;
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (p_x[i] != p_x[i]){
        out = true;
        break;
      }
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (p_x[i] == NA_STRING){
        out = true;
        break;
      }
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (( ((p_x[i]).r != (p_x[i]).r) || ((p_x[i]).i != (p_x[i]).i) )){
        out = true;
        break;
      }
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    for (int i = 0; i < n; ++i){
      out = cpp_any_na(VECTOR_ELT(x, i), true);
      if (out) break;
    }
    break;
  }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x));
    ++n_prot;
    SEXP any_missing = Rf_protect(cpp11::package("base")["any"](is_missing));
    ++n_prot;
    out = Rf_asLogical(any_missing);
    break;
  }
  }
  Rf_unprotect(n_prot);
  return out;
}

[[cpp11::register]]
bool cpp_all_na(SEXP x, bool return_true_on_empty, bool recursive){
  int n_prot = 0;
  R_xlen_t n = Rf_xlength(x);
  bool out = true;
  if (n == 0){
    if (return_true_on_empty){
      return true;
    } else {
      return false;
    }
  }
  switch ( TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (p_x[i] != NA_INTEGER){
        out = false;
        break;
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (p_x[i] == p_x[i]){
        out = false;
        break;
      }
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (p_x[i] != NA_STRING){
        out = false;
        break;
      }
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (( ((p_x[i]).r == (p_x[i]).r) && ((p_x[i]).i == (p_x[i]).i) )){
        out = false;
        break;
      }
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    for (int i = 0; i < n; ++i){
      out = cpp_all_na(VECTOR_ELT(x, i), return_true_on_empty, true);
      if (!out) break;
    }
    break;
  }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x));
    ++n_prot;
    SEXP all_missing = Rf_protect(cpp11::package("base")["all"](is_missing));
    ++n_prot;
    out = Rf_asLogical(all_missing);
    break;
  }
  }
  Rf_unprotect(n_prot);
  return out;
}

// A multi-threaded version of `is.na()`
// lists are handled differently in that each element
// must contain only NA in all nested elements to be regarded as NA

[[cpp11::register]]
SEXP cpp_is_na(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= 100000 ? num_cores() : 1;
  SEXP out;
  switch ( TYPEOF(x) ){
  case NILSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, 0));
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    int *p_x = INTEGER(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = (p_x[i] == NA_INTEGER);
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = (p_x[i] == NA_INTEGER);
    }

    break;
  }
  case REALSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    double *p_x = REAL(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = (p_x[i] != p_x[i]);
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = (p_x[i] != p_x[i]);
    }
    break;
  }
  case STRSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    const SEXP *p_x = STRING_PTR_RO(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = (p_x[i] == NA_STRING);
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i) p_out[i] = (p_x[i] == NA_STRING);
    }
    break;
  }
  case RAWSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    memset(p_out, 0, n * sizeof(int));
    break;
  }
  case CPLXSXP: {
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    Rcomplex *p_x = COMPLEX(x);
    if (n_cores > 1){
      OMP_PARALLEL_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        p_out[i] = ( ((p_x[i]).r != (p_x[i]).r) || ((p_x[i]).i != (p_x[i]).i) );
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        p_out[i] = ( ((p_x[i]).r != (p_x[i]).r) || ((p_x[i]).i != (p_x[i]).i) );
      }
    }
    break;
  }
  case VECSXP: {
    if (!Rf_isObject(x)){
    out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = cpp_all_na(p_x[i], false, true);
    }
    break;
  }
  }
  default: {
    out = Rf_protect(cpp11::package("base")["is.na"](x));
    break;
  }
  }
  Rf_unprotect(1);
  return out;
}

// Memory-efficient which(is.na(x))

[[cpp11::register]]
SEXP cpp_which_na(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  bool is_short = (n <= integer_max_);
  switch ( TYPEOF(x) ){
  case NILSXP: {
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_unprotect(1);
    return out;
  }
  case LGLSXP:
  case INTSXP: {
    R_xlen_t count = na_count(x, true);
    int *p_x = INTEGER(x);
    if (is_short){
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, count));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] == NA_INTEGER);
      }
      Rf_unprotect(1);
      return out;
    } else {
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, count));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] == NA_INTEGER);
      }
      Rf_unprotect(1);
      return out;
    }
  }
  case REALSXP: {
    R_xlen_t count = na_count(x, true);
    double *p_x = REAL(x);
    if (is_short){
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, count));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += (p_x[i] != p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    } else {
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, count));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += (p_x[i] != p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    }
  }
  case STRSXP: {
    R_xlen_t count = na_count(x, true);
    const SEXP *p_x = STRING_PTR_RO(x);
    if (is_short){
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, count));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] == NA_STRING);
      }
      Rf_unprotect(1);
      return out;
    } else {
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, count));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] == NA_STRING);
      }
      Rf_unprotect(1);
      return out;
    }
  }
  case RAWSXP: {
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_unprotect(1);
    return out;
  }
  case CPLXSXP: {
    R_xlen_t count = na_count(x, true);
    Rcomplex *p_x = COMPLEX(x);
    if (is_short){
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, count));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += ( ( !((p_x[i]).r == (p_x[i]).r) ) || ( !((p_x[i]).i == (p_x[i]).i) ) );
        ++i;
      }
      Rf_unprotect(1);
      return out;
    } else {
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, count));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < count){
        p_out[whichi] = i + 1;
        whichi += ( ( !((p_x[i]).r == (p_x[i]).r) ) || ( !((p_x[i]).i == (p_x[i]).i) ) );
        ++i;
      }
      Rf_unprotect(1);
      return out;
    }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x));
    SEXP out = Rf_protect(cpp_which_(is_missing, false));
    Rf_unprotect(2);
    return out;
    break;
  }
  }
}

[[cpp11::register]]
SEXP cpp_which_not_na(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  bool is_short = (n <= integer_max_);
  switch ( TYPEOF(x) ){
  case NILSXP: {
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, 0));
    Rf_unprotect(1);
    return out;
  }
  case LGLSXP:
  case INTSXP: {
    R_xlen_t count = na_count(x, true);
    int *p_x = INTEGER(x);
    if (is_short){
      int out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] != NA_INTEGER);
      }
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] != NA_INTEGER);
      }
      Rf_unprotect(1);
      return out;
    }
  }
  case REALSXP: {
    R_xlen_t count = na_count(x, true);
    double *p_x = REAL(x);
    if (is_short){
      int out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i] == p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i] == p_x[i]);
        ++i;
      }
      Rf_unprotect(1);
      return out;
    }
  }
  case STRSXP: {
    R_xlen_t count = na_count(x, true);
    const SEXP *p_x = STRING_PTR_RO(x);
    if (is_short){
      int out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] != NA_STRING);
      }
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += (p_x[i++] != NA_STRING);
      }
      Rf_unprotect(1);
      return out;
    }
  }
  case RAWSXP: {
    if (is_short){
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
    int *p_out = INTEGER(out);
    for (int i = 0; i < n; ++i){
      p_out[i] = i + 1;
    }
    Rf_unprotect(1);
    return out;
  } else {
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
    double *p_out = REAL(out);
    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = i + 1;
    }
    Rf_unprotect(1);
    return out;
  }
  }
  case CPLXSXP: {
    R_xlen_t count = na_count(x, true);
    Rcomplex *p_x = COMPLEX(x);
    if (is_short){
      int out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
      int *p_out = INTEGER(out);
      int whichi = 0;
      int i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += ( ( ((p_x[i]).r == (p_x[i]).r) ) && ( ((p_x[i]).i == (p_x[i]).i) ) );
        ++i;
      }
      Rf_unprotect(1);
      return out;
    } else {
      R_xlen_t out_size = n - count;
      SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
      double *p_out = REAL(out);
      R_xlen_t whichi = 0;
      R_xlen_t i = 0;
      while (whichi < out_size){
        p_out[whichi] = i + 1;
        whichi += ( ( ((p_x[i]).r == (p_x[i]).r) ) && ( ((p_x[i]).i == (p_x[i]).i) ) );
        ++i;
      }
      Rf_unprotect(1);
      return out;
    }
  }
  default: {
    SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](x));
    SEXP out = Rf_protect(cpp_which_(is_missing, true));
    Rf_unprotect(2);
    return out;
    break;
  }
  }
}

[[cpp11::register]]
SEXP cpp_row_na_counts(SEXP x){
  if (!Rf_isFrame(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  int num_col = Rf_length(x);
  int n_prot = 0;
  R_xlen_t num_row = cpp_df_nrow(x);
  SEXP n_empty = Rf_protect(Rf_allocVector(INTSXP, num_row));
  ++n_prot;
  int *p_n_empty = INTEGER(n_empty);
  memset(p_n_empty, 0, num_row * sizeof(int));
  int do_parallel = num_row >= 100000;
  int n_cores = do_parallel ? num_cores() : 1;
  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case LGLSXP:
    case INTSXP: {
      int *p_xj = INTEGER(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_n_empty[i] += (p_xj[i] == NA_INTEGER);
      }
      break;
    }
    case REALSXP: {
      double *p_xj = REAL(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_n_empty[i] += (p_xj[i] != p_xj[i]);
      }
      break;
    }
    case STRSXP: {
      const SEXP *p_xj = STRING_PTR_RO(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_n_empty[i] += (p_xj[i] == NA_STRING);
      }
      break;
    }
    case RAWSXP: {
      break;
    }
    case CPLXSXP: {
      Rcomplex *p_xj = COMPLEX(p_x[j]);
#pragma omp parallel num_threads(n_cores) if(do_parallel)
#pragma omp for simd
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_n_empty[i] += (p_xj[i]).r != (p_xj[i]).r || (p_xj[i]).i != (p_xj[i]).i;
      }
      break;
    }
    case VECSXP: {
      if (Rf_isObject(p_x[j])){
      SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](p_x[j]));
      ++n_prot;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing);
        ++n_prot;
        SEXP names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
        Rf_unprotect(n_prot);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 CHAR(STRING_ELT(names, j)), element_length, int_nrows);
      }
      int *p_is_missing = LOGICAL(is_missing);
      for (R_xlen_t k = 0; k < num_row; ++k){
        p_n_empty[k] += p_is_missing[k];
      }
    } else {
      const SEXP *p_xj = VECTOR_PTR_RO(p_x[j]);
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_n_empty[i] += cpp_all_na(p_xj[i], false, true);
      }
    }
    break;
    }
    default: {
      Rf_unprotect(n_prot);
      Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(p_x[j])));
    }
    }
  }
  Rf_unprotect(n_prot);
  return n_empty;
}

[[cpp11::register]]
SEXP cpp_col_na_counts(SEXP x){
  if (!Rf_isFrame(x)){
    Rf_error("x must be a data frame");
  }
  const SEXP *p_x = VECTOR_PTR_RO(x);
  int num_col = Rf_length(x);
  int n_prot = 0;
  R_xlen_t num_row = cpp_df_nrow(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, num_col));
  ++n_prot;
  int *p_out = INTEGER(out);
  memset(p_out, 0, num_col * sizeof(int));
  for (int j = 0; j < num_col; ++j){
    switch ( TYPEOF(p_x[j]) ){
    case VECSXP: {
      if (Rf_isObject(p_x[j])){
      SEXP is_missing = Rf_protect(cpp11::package("cheapr")["is_na"](p_x[j]));
      ++n_prot;
      if (Rf_xlength(is_missing) != num_row){
        int int_nrows = num_row;
        int element_length = Rf_xlength(is_missing);
        ++n_prot;
        SEXP names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
        Rf_unprotect(n_prot);
        Rf_error("is.na method for list variable %s produces a length (%d) vector which does not equal the number of rows (%d)",
                 CHAR(STRING_ELT(names, j)), element_length, int_nrows);
      }
      int *p_is_missing = LOGICAL(is_missing);
      for (R_xlen_t k = 0; k < num_row; ++k){
        p_out[j] += p_is_missing[k];
      }
    } else {
      for (R_xlen_t i = 0; i < num_row; ++i){
        p_out[j] += cpp_all_na(VECTOR_ELT(p_x[j], i), false, true);
      }
    }
    break;
    }
    default: {
      p_out[j] = na_count(p_x[j], false);
      break;
    }
    }
  }
  Rf_unprotect(n_prot);
  return out;
}

R_xlen_t cpp_clean_threshold(double threshold, bool threshold_is_prop, R_xlen_t n){
  if (threshold != threshold){
    Rf_error("threshold cannot be NA");
  }
  R_xlen_t out = threshold;
  if (threshold_is_prop){
    if (threshold < 0){
      out = 0;
    } else if (threshold == R_PosInf){
      out = n + 1;
    } else {
      out = std::floor( (threshold * n) + 0.0000000001);
    }
  } else {
    if (threshold < 0){
      out = 0;
    }
    if (threshold == R_PosInf){
      out = n + 1;
    }
  }
  return out;
}

// Matrix methods
// The methods for matrices are substantially different

[[cpp11::register]]
SEXP cpp_matrix_row_na_counts(SEXP x){
  if (!Rf_isMatrix(x)){
    Rf_error("x must be a matrix");
  }
  int num_row = Rf_nrows(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, num_row));
  int *p_out = INTEGER(out);
  memset(p_out, 0, num_row * sizeof(int));
  int do_parallel = n >= 100000;
  int n_cores = do_parallel ? num_cores() : 1;
#pragma omp parallel num_threads(n_cores) if(do_parallel)
  switch ( TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += (p_x[i] == NA_INTEGER);
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += (p_x[i] != p_x[i]);
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += (p_x[i] == NA_STRING);
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[i % num_row] += (p_x[i]).r != (p_x[i]).r || (p_x[i]).i != (p_x[i]).i;
    }
    break;
  }
  default: {
    Rf_unprotect(1);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(1);
  return out;
}
[[cpp11::register]]
SEXP cpp_matrix_col_na_counts(SEXP x){
  if (!Rf_isMatrix(x)){
    Rf_error("x must be a matrix");
  }
  int num_row = Rf_nrows(x);
  int num_col = Rf_ncols(x);
  R_xlen_t n = Rf_xlength(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, num_col));
  int *p_out = INTEGER(out);
  memset(p_out, 0, num_col * sizeof(int));
  int do_parallel = n >= 100000;
  int n_cores = do_parallel ? num_cores() : 1;
#pragma omp parallel num_threads(n_cores) if(do_parallel)
  switch ( TYPEOF(x) ){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += (p_x[i] == NA_INTEGER);
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += (p_x[i] != p_x[i]);
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += (p_x[i] == NA_STRING);
    }
    break;
  }
  case RAWSXP: {
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
#pragma omp for
    for (R_xlen_t i = 0; i < n; ++i){
#pragma omp atomic
      p_out[int_div(i, num_row)] += (p_x[i]).r != (p_x[i]).r || (p_x[i]).i != (p_x[i]).i;
    }
    break;
  }
  default: {
    Rf_unprotect(1);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(1);
  return out;
}
