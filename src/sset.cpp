#include "cheapr_cpp.h"
// #include <vector>
// using namespace cpp11;

// Subsetting vectors and data frames
// Includes a unique optimisation on range subsetting

// Altrep utils

bool is_altrep(SEXP x){
  return ALTREP(x);
}
SEXP alt_class(SEXP x){
  if (is_altrep(x)){
    SEXP alt_attribs = Rf_protect(Rf_coerceVector(ATTRIB(ALTREP_CLASS(x)), VECSXP));
    SEXP out = Rf_protect(Rf_coerceVector(VECTOR_ELT(alt_attribs, 0), STRSXP));
    Rf_unprotect(2);
    return out;
  } else {
    return R_NilValue;
  }
}
SEXP alt_pkg(SEXP x){
  if (is_altrep(x)){
    SEXP alt_attribs = Rf_protect(Rf_coerceVector(ATTRIB(ALTREP_CLASS(x)), VECSXP));
    SEXP out = Rf_protect(Rf_coerceVector(VECTOR_ELT(alt_attribs, 1), STRSXP));
    Rf_unprotect(2);
    return out;
  } else {
    return R_NilValue;
  }
}

SEXP alt_data1(SEXP x){
  if (is_altrep(x)){
    return R_altrep_data1(x);
  } else {
    return R_NilValue;
  }
}
[[cpp11::register]]
bool is_alt_compact_seq(SEXP x){
  if (!is_altrep(x)) return false;
  SEXP alt_class_nm = Rf_protect(alt_class(x));
  SEXP alt_pkg_nm = Rf_protect(alt_pkg(x));
  SEXP intseq_char = Rf_protect(Rf_mkChar("compact_intseq"));
  SEXP realseq_char = Rf_protect(Rf_mkChar("compact_realseq"));
  SEXP base_char = Rf_protect(Rf_mkChar("base"));
  bool out = (STRING_ELT(alt_class_nm, 0) == intseq_char ||
              STRING_ELT(alt_class_nm, 0) == realseq_char) &&
              STRING_ELT(alt_pkg_nm, 0) == base_char;
  Rf_unprotect(5);
  return out;
}

[[cpp11::register]]
SEXP alt_compact_seq_data(SEXP x){
  if (!is_alt_compact_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = Rf_protect(Rf_coerceVector(alt_data1(x), REALSXP));
  double alt_size = REAL(alt_data)[0];
  double alt_from = REAL(alt_data)[1];
  double alt_by = REAL(alt_data)[2];
  double alt_to = (std::fmax(alt_size - 1.0, 0.0) * alt_by) + alt_from;
  SEXP out = Rf_protect(Rf_allocVector(REALSXP, 4));
  double *p_out = REAL(out);
  p_out[0] = alt_from;
  p_out[1] = alt_to;
  p_out[2] = alt_by;
  p_out[3] = alt_size;
  Rf_unprotect(2);
  return out;
}

// Subset with no checks, indices vector must be pre-curated

SEXP cpp_sset_unsafe(SEXP x, int *pind, int n, int n_cores){
  bool do_parallel = n_cores > 1;
  switch ( TYPEOF(x) ){

  case NILSXP: {
    return R_NilValue;
  }
  case LGLSXP: {
    int *p_x = LOGICAL(x);
    SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
    int *p_out = LOGICAL(out);
    if (do_parallel){
      OMP_PARALLEL_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    } else {
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case INTSXP: {
    int *p_x = INTEGER(x);
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
    int *p_out = INTEGER(out);
    if (do_parallel){
      OMP_PARALLEL_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    } else {
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
    double *p_out = REAL(out);
    if (do_parallel){
      OMP_PARALLEL_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    } else {
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case STRSXP: {
    SEXP *p_x = STRING_PTR(x);
    SEXP out = Rf_protect(Rf_allocVector(STRSXP, n));
    OMP_FOR_SIMD
    for (int i = 0; i < n; ++i){
      SET_STRING_ELT(out, i, p_x[pind[i] - 1]);
    }
    Rf_unprotect(1);
    return out;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    SEXP out = Rf_protect(Rf_allocVector(CPLXSXP, n));
    OMP_FOR_SIMD
    for (int i = 0; i < n; ++i){
      SET_COMPLEX_ELT(out, i, p_x[pind[i] - 1]);
    }
    Rf_unprotect(1);
    return out;
  }
  case RAWSXP: {
    Rbyte *p_x = RAW(x);
    SEXP out = Rf_protect(Rf_allocVector(RAWSXP, n));
    OMP_FOR_SIMD
    for (int i = 0; i < n; ++i){
      SET_RAW_ELT(out, i, p_x[pind[i] - 1]);
    }
    Rf_unprotect(1);
    return out;
  }
  case VECSXP: {
    const SEXP *p_x = VECTOR_PTR_RO(x);
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, n));
    OMP_FOR_SIMD
    for (int i = 0; i < n; ++i){
      SET_VECTOR_ELT(out, i, p_x[pind[i] - 1]);
    }
    Rf_unprotect(1);
    return out;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}

// A range-based subset method
// Can be readily used when indices are an altrep compact integer sequence
// Also works with negative-indexing

[[cpp11::register]]
SEXP cpp_sset_range(SEXP x, R_xlen_t from, R_xlen_t to, R_xlen_t by){
  R_xlen_t n = Rf_xlength(x);
  if (by != 1 && by != -1){
    Rf_error("by increment must be 1 or -1");
  }
  if ( (from > 0 && to < 0) ||
       (from < 0 && to > 0) ){
    Rf_error("Cannot mix positive and negative indices");
  }
  if ( (from > to && by > 0) || (from < to && by < 0)){
    Rf_error("Wrong increment sign in by arg");
  }
  R_xlen_t istart = from;
  R_xlen_t iend = to;
  R_xlen_t out_size, istart1, istart2, iend1, iend2, diff;
  bool double_loop = false;
  // Negative indexing is complicated

  // Assuming N = length(x)
  // Because there are 6 scenarios...

  // Scenario 1: abs(m) <= N but abs(n) > N
  // Scenario 2: both abs(m) & abs(n) > N
  // Scenario 3: -1:-N (Empty vector result)
  // Scenario 4: -1:mE{m:-m < N} - We subset from n+1:N
  // Scenario 5: k:NE{k:-k < N}
  // Scenario 6: m:kE{-m > 1 & -k < N} - Everything but middle chunk of vector

  if (istart == 0 && iend == 0){
    istart = 1;
    iend = 0;
    by = 1;
    out_size = 0;
  } else if (istart < 0 || iend < 0){
    if (istart == 0){
      istart = -1;
    }
    if (iend == 0){
      iend = -1;
    }
    // We first switch them
    if (istart < iend){
      R_xlen_t iend_temp = iend;
      iend = istart;
      istart = iend_temp;
    }
    // Out-of-bounds adjustments

    // Scenario 1
    if (std::abs(istart) <= n && std::abs(iend) > n){
      iend = std::abs(istart) - 1;
      istart = 1;
      by = 1;
      out_size = (iend - istart) + 1;
      // Scenario 2
    } else if (std::abs(istart) > n && std::abs(iend) > n){
      istart = 1;
      iend  = n;
      by = 1;
      out_size = n;
      // Scenario 3
    } else if (istart == -1 && iend == -n){
      istart = n + 1;
      iend = n;
      by = 1;
      out_size = 0;
      // Scenario 4
    } else if (istart == -1 && std::abs(iend) < n){
      istart = std::abs(iend) + 1;
      iend = n;
      by = 1;
      out_size = (iend - istart) + 1;
      // Scenario 5
    } else if (std::abs(istart) < n && std::abs(iend) == n){
      iend = std::abs(istart) - 1;
      istart = 1;
      by = 1;
      out_size = (iend - istart) + 1;
    } else {
      // Scenario 6
      double_loop = true;
      istart1 = 1;
      iend1 = std::abs(istart) - 1;
      istart2 = std::abs(iend) + 1;
      iend2 = n;
      by = 1;
      out_size = (iend1 - istart1) + (iend2 - istart2) + 2;
    }
  } else {
    if (istart == 0){
      istart = 1;
    }
    if (iend == 0){
      iend = 1;
    }
    out_size = ((iend - istart) / by) + 1;
  }
  bool do_parallel = out_size >= 100000;
  int n_cores = do_parallel ? num_cores() : 1;
  unsigned int k = 0;

  switch ( TYPEOF(x) ){
  case NILSXP: {
    return R_NilValue;
  }
  case LGLSXP: {
    int *p_x = LOGICAL(x);
    SEXP out = Rf_protect(Rf_allocVector(LGLSXP, out_size));
    int *p_out = LOGICAL(out);
    if (double_loop){
      diff = iend1 - istart1 + 2;
#pragma omp parallel num_threads(n_cores) if (do_parallel)
      OMP_FOR_SIMD
      for (R_xlen_t i = istart1 - 1; i < iend1; ++i){
        p_out[i - istart1 + 1] = i < n ? p_x[i] : NA_LOGICAL;
      }
      OMP_FOR_SIMD
      for (R_xlen_t i = istart2 - 1; i < iend2; ++i){
        p_out[i - istart2 + diff] = i < n ? p_x[i] : NA_LOGICAL;
      }
    } else {
      if (by > 0){
        if (do_parallel){
          OMP_PARALLEL_FOR_SIMD
          for (R_xlen_t i = istart - 1; i < iend; ++i){
            p_out[i - istart + 1] = i < n ? p_x[i] : NA_LOGICAL;
          }
        } else {
          OMP_FOR_SIMD
          for (R_xlen_t i = istart - 1; i < iend; ++i){
            p_out[i - istart + 1] = i < n ? p_x[i] : NA_LOGICAL;
          }
        }
      } else {
        if (do_parallel){
          OMP_PARALLEL_FOR_SIMD
          for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
            p_out[istart - i - 1] = i < n ? p_x[i] : NA_LOGICAL;
          }
        } else {
          OMP_FOR_SIMD
          for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
            p_out[istart - i - 1] = i < n ? p_x[i] : NA_LOGICAL;
          }
        }
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case INTSXP: {
    int *p_x = INTEGER(x);
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
    int *p_out = INTEGER(out);
    if (double_loop){
      diff = iend1 - istart1 + 2;
#pragma omp parallel num_threads(n_cores) if (do_parallel)
      OMP_FOR_SIMD
      for (R_xlen_t i = istart1 - 1; i < iend1; ++i){
        p_out[i - istart1 + 1] = (i < n ? p_x[i] : NA_INTEGER);
      }
      OMP_FOR_SIMD
      for (R_xlen_t i = istart2 - 1; i < iend2; ++i){
        p_out[i - istart2 + diff] = (i < n ? p_x[i] : NA_INTEGER);
      }
    } else {
      if (by > 0){
        if (do_parallel){
          OMP_PARALLEL_FOR_SIMD
          for (R_xlen_t i = istart - 1; i < iend; ++i){
            p_out[i - istart + 1] = i < n ? p_x[i] : NA_INTEGER;
          }
        } else {
          OMP_FOR_SIMD
          for (R_xlen_t i = istart - 1; i < iend; ++i){
            p_out[i - istart + 1] = i < n ? p_x[i] : NA_INTEGER;
          }
        }

      } else {
        if (do_parallel){
          OMP_PARALLEL_FOR_SIMD
          for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
            p_out[istart - i - 1] = i < n ? p_x[i] : NA_INTEGER;
          }
        } else {
          OMP_FOR_SIMD
          for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
            p_out[istart - i - 1] = i < n ? p_x[i] : NA_INTEGER;
          }
        }
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
    double *p_out = REAL(out);
    if (double_loop){
      diff = iend1 - istart1 + 2;
#pragma omp parallel num_threads(n_cores) if (do_parallel)
      OMP_FOR_SIMD
      for (R_xlen_t i = (istart1 - 1); i < iend1; ++i){
        p_out[i - istart1 + 1] = (i < n ? p_x[i] : NA_REAL);
      }
      OMP_FOR_SIMD
      for (R_xlen_t i = istart2 - 1; i < iend2; ++i){
        p_out[i - istart2 + diff] = i < n ? p_x[i] : NA_REAL;
      }
    } else {
      if (by > 0){
        if (do_parallel){
          OMP_PARALLEL_FOR_SIMD
          for (R_xlen_t i = istart - 1; i < iend; ++i){
            p_out[i - istart + 1] = i < n ? p_x[i] : NA_REAL;
          }
        } else {
          OMP_FOR_SIMD
          for (R_xlen_t i = istart - 1; i < iend; ++i){
            p_out[i - istart + 1] = i < n ? p_x[i] : NA_REAL;
          }
        }

      } else {
        if (do_parallel){
          OMP_PARALLEL_FOR_SIMD
          for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
            p_out[istart - i - 1] = i < n ? p_x[i] : NA_REAL;
          }
        } else {
          OMP_FOR_SIMD
          for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
            p_out[istart - i - 1] = i < n ? p_x[i] : NA_REAL;
          }
        }

      }
    }
    Rf_unprotect(1);
    return out;
  }
  case STRSXP: {
    SEXP *p_x = STRING_PTR(x);
    SEXP out = Rf_protect(Rf_allocVector(STRSXP, out_size));
    if (double_loop){
      OMP_FOR_SIMD
      for (R_xlen_t i = istart1 - 1; i < iend1; ++i){
        SET_STRING_ELT(out, k++, i < n ? p_x[i] : NA_STRING);
      }
      OMP_FOR_SIMD
      for (R_xlen_t j = istart2 - 1; j < iend2; ++j){
        SET_STRING_ELT(out, k++, j < n ? p_x[j] : NA_STRING);
      }
    } else {
      if (by > 0){
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i < iend; ++i){
          SET_STRING_ELT(out, k++, i < n ? p_x[i] : NA_STRING);
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
          SET_STRING_ELT(out, k++, i < n ? p_x[i] : NA_STRING);
        }
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    SEXP out = Rf_protect(Rf_allocVector(CPLXSXP, out_size));
    SEXP na_complex_sexp = Rf_protect(Rf_allocVector(CPLXSXP, 1));
    Rcomplex *p_na_complex = COMPLEX(na_complex_sexp);
    p_na_complex[0].i = NA_REAL;
    p_na_complex[0].r = NA_REAL;
    Rcomplex na_complex = Rf_asComplex(na_complex_sexp);
    if (double_loop){
      OMP_FOR_SIMD
      for (R_xlen_t i = istart1 - 1; i < iend1; ++i){
        SET_COMPLEX_ELT(out, k++, i < n ? p_x[i] : na_complex);
      }
      OMP_FOR_SIMD
      for (R_xlen_t j = istart2 - 1; j < iend2; ++j){
        SET_COMPLEX_ELT(out, k++, j < n ? p_x[j] : na_complex);
      }
    } else {
      if (by > 0){
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i < iend; ++i){
          SET_COMPLEX_ELT(out, k++, i < n ? p_x[i] : na_complex);
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
          SET_COMPLEX_ELT(out, k++, i < n ? p_x[i] : na_complex);
        }
      }
    }
    Rf_unprotect(2);
    return out;
  }
  case RAWSXP: {
    Rbyte *p_x = RAW(x);
    SEXP out = Rf_protect(Rf_allocVector(RAWSXP, out_size));
    if (double_loop){
      OMP_FOR_SIMD
      for (R_xlen_t i = istart1 - 1; i < iend1; ++i){
        SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
      }
      OMP_FOR_SIMD
      for (R_xlen_t j = istart2 - 1; j < iend2; ++j){
        SET_RAW_ELT(out, k++, j < n ? p_x[j] : 0);
      }
    } else {
      if (by > 0){
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i < iend; ++i){
          SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
          SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
        }
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case VECSXP: {
    const SEXP *p_x = VECTOR_PTR_RO(x);
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, out_size));
    if (double_loop){
      OMP_FOR_SIMD
      for (R_xlen_t i = istart1 - 1; i < iend1; ++i){
        if (i < n){
          SET_VECTOR_ELT(out, k, p_x[i]);
        }
        ++k;
      }
      OMP_FOR_SIMD
      for (R_xlen_t j = istart2 - 1; j < iend2; ++j){
        if (j < n){
          SET_VECTOR_ELT(out, k, p_x[j]);
        }
        ++k;
      }
    } else {
      if (by > 0){
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i < iend; ++i){
          if (i < n){
            SET_VECTOR_ELT(out, k, p_x[i]);
          }
          ++k;
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
          if (i < n){
            SET_VECTOR_ELT(out, k, p_x[i]);
          }
          ++k;
        }
      }
    }
    Rf_unprotect(1);
    return out;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}

// Helper to convert altrep sequences into the final subsetted length

R_xlen_t get_alt_final_sset_size(R_xlen_t n, R_xlen_t from, R_xlen_t to, R_xlen_t by){
  R_xlen_t istart = from;
  R_xlen_t iend = to;
  R_xlen_t out;
  if (istart == 0 && iend == 0){
    out = 0;
  } else if (istart < 0 || iend < 0){
    if (istart == 0){
      istart = -1;
    }
    if (iend == 0){
      iend = -1;
    }
    // We first switch them
    if (istart < iend){
      R_xlen_t iend_temp = iend;
      iend = istart;
      istart = iend_temp;
    }
    // Out-of-bounds adjustments

    // Scenario 1
    if (std::abs(istart) <= n && std::abs(iend) > n){
      iend = std::abs(istart) - 1;
      istart = 1;
      out = (iend - istart) + 1;
      // Scenario 2
    } else if (std::abs(istart) > n && std::abs(iend) > n){
      out = n;
      // Scenario 3
    } else if (istart == -1 && iend == -n){
      out = 0;
      // Scenario 4
    } else if (istart == -1 && std::abs(iend) < n){
      istart = std::abs(iend) + 1;
      iend = n;
      out = (iend - istart) + 1;
      // Scenario 5
    } else if (std::abs(istart) < n && std::abs(iend) == n){
      iend = std::abs(istart) - 1;
      istart = 1;
      out = (iend - istart) + 1;
    } else {
      // Scenario 6
      R_xlen_t istart1 = 1;
      R_xlen_t iend1 = std::abs(istart) - 1;
      R_xlen_t istart2 = std::abs(iend) + 1;
      R_xlen_t iend2 = n;
      out = (iend1 - istart1) + (iend2 - istart2) + 2;
    }
  } else {
    if (istart == 0){
      istart = 1;
    }
    if (iend == 0){
      iend = 1;
    }
    out = ((iend - istart) / by) + 1;
  }
  return out;
}

[[cpp11::register]]
SEXP cpp_sset_df(SEXP x, SEXP indices){
  int xn = cpp_df_nrow(x);
  int ncols = Rf_length(x);
  int n = Rf_length(indices);
  int n_protections = 0;
  int zero_count = 0;
  int pos_count = 0;
  int oob_count = 0;
  int na_count = 0;
  int out_size;
  bool do_parallel = n >= 10000;
  int n_cores = do_parallel ? num_cores() : 1;
  cpp11::function cheapr_sset = cpp11::package("cheapr")["sset"];
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, ncols));
  ++n_protections;
  // SEXP *p_out = VECTOR_PTR(out);

  // If indices is a special type of ALTREP compact int sequence, we can
  // Use a range-based subset instead

  if (is_alt_compact_seq(indices)){

    // ALTREP integer sequence method

    SEXP seq_data = Rf_protect(alt_compact_seq_data(indices));
    ++n_protections;
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    for (int j = 0; j < ncols; ++j){
      SEXP df_var = Rf_protect(p_x[j]);
      if (!Rf_isObject(df_var) ||
          Rf_inherits(df_var, "Date") ||
          Rf_inherits(df_var, "POSIXct") ||
          Rf_inherits(df_var, "factor")){
        SEXP list_var = Rf_protect(cpp_sset_range(df_var, from, to, by));
        Rf_copyMostAttrib(df_var, list_var);
        int has_names = Rf_getAttrib(df_var, R_NamesSymbol) != R_NilValue;
        if (has_names){
          SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
          SEXP new_names = Rf_protect(cpp_sset_range(
            old_names, from, to, by)
          );
          Rf_setAttrib(list_var, R_NamesSymbol, new_names);
        }
        // R_PreserveObject(list_var);
        // p_out[j] = list_var; // Unsafe?  (test this by running sset(iris, 1:10^7))

        SET_VECTOR_ELT(out, j, list_var);

        // We un-protect below the original df variables, as well as
        // old and new names
        // Once they are added to the data frame, they are protected
        // If we didn't do this we would easily reach the protection stack limit
        // of 10,000
        Rf_unprotect(1 + (has_names * 2));
      } else {
        SET_VECTOR_ELT(out, j, cheapr_sset(df_var, indices));
      }
      // Unprotecting new data frame variable
      Rf_unprotect(1);
    }
    out_size = get_alt_final_sset_size(xn, from, to, by);
  } else {
    int *pi = INTEGER(indices);

    // Counting the number of:
    // Zeroes
    // Out-of-bounds indices
    // Positive indices
    // NA indices
    // From this we can also work out the number of negatives

    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:zero_count,pos_count,oob_count,na_count)
      for (int j = 0; j < n; ++j){
        zero_count += (pi[j] == 0);
        pos_count += (pi[j] > 0);
        oob_count += (std::fabs(pi[j]) > xn);
        na_count += (pi[j] == NA_INTEGER);
      }
    } else {
      OMP_FOR_SIMD
      for (int j = 0; j < n; ++j){
        zero_count += (pi[j] == 0);
        pos_count += (pi[j] > 0);
        oob_count += (std::fabs(pi[j]) > xn);
        na_count += (pi[j] == NA_INTEGER);
      }
    }
    int neg_count = n - pos_count - zero_count - na_count;
    if ( (pos_count > 0 && neg_count > 0) ||
         (neg_count > 0 && na_count > 0)){
      Rf_error("Cannot mix positive and negative indices");
    }
    // Should a simplified sset method be used?

    bool simple_sset = zero_count == 0 && oob_count == 0 && na_count == 0 && pos_count == n;

    // Final length of output

    out_size = na_count + pos_count;

    // If Index vector is clean we can use fast subset

    if (simple_sset){
      for (int j = 0; j < ncols; ++j){
        SEXP df_var = Rf_protect(p_x[j]);
        if (!Rf_isObject(df_var) ||
            Rf_inherits(df_var, "Date") ||
            Rf_inherits(df_var, "POSIXct") ||
            Rf_inherits(df_var, "factor")){
          SEXP list_var = Rf_protect(cpp_sset_unsafe(df_var, pi, out_size, n_cores));
          Rf_copyMostAttrib(df_var, list_var);
          int has_names = Rf_getAttrib(df_var, R_NamesSymbol) != R_NilValue;
          if (has_names){
            SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
            SEXP new_names = Rf_protect(cpp_sset_unsafe(
              old_names, pi, out_size, n_cores
            ));
            Rf_setAttrib(list_var, R_NamesSymbol, new_names);
          }
          SET_VECTOR_ELT(out, j, list_var);
          Rf_unprotect(1 + (has_names * 2));
        } else {
          SET_VECTOR_ELT(out, j, cheapr_sset(df_var, indices));
        }
        Rf_unprotect(1);
      }

      // Negative indexing

    } else if (neg_count > 0){
      SEXP indices2 = Rf_protect(cpp11::package("cheapr")["neg_indices_to_pos"](indices, xn));
      ++n_protections;
      out_size = Rf_length(indices2);
      int *pi2 = INTEGER(indices2);
      for (int j = 0; j < ncols; ++j){
        SEXP df_var = Rf_protect(p_x[j]);
        if (!Rf_isObject(df_var) ||
            Rf_inherits(df_var, "Date") ||
            Rf_inherits(df_var, "POSIXct") ||
            Rf_inherits(df_var, "factor")){
          SEXP list_var = Rf_protect(cpp_sset_unsafe(df_var, pi2, out_size, n_cores));
          Rf_copyMostAttrib(df_var, list_var);
          int has_names = Rf_getAttrib(df_var, R_NamesSymbol) != R_NilValue;
          if (has_names){
            SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
            SEXP new_names = Rf_protect(cpp_sset_unsafe(
              old_names, pi2, out_size, n_cores
            ));
            Rf_setAttrib(list_var, R_NamesSymbol, new_names);
          }
          SET_VECTOR_ELT(out, j, list_var);
          Rf_unprotect(1 + (has_names * 2));
        } else {
          SET_VECTOR_ELT(out, j, cheapr_sset(df_var, indices2));
        }
        Rf_unprotect(1);
      }
      // If index vector is clean except for existence of zeroes
    } else if (zero_count > 0 && oob_count == 0 && na_count == 0){
      SEXP r_zero = Rf_protect(Rf_ScalarInteger(0));
      ++n_protections;
      SEXP indices2 = Rf_protect(cpp11::package("cheapr")["val_rm"](indices, r_zero));
      ++n_protections;
      int *pi2 = INTEGER(indices2);
      for (int j = 0; j < ncols; ++j){
        SEXP df_var = Rf_protect(p_x[j]);
        if (!Rf_isObject(df_var) ||
            Rf_inherits(df_var, "Date") ||
            Rf_inherits(df_var, "POSIXct") ||
            Rf_inherits(df_var, "factor")){
          SEXP list_var = Rf_protect(cpp_sset_unsafe(df_var, pi2, out_size, n_cores));
          Rf_copyMostAttrib(df_var, list_var);
          int has_names = Rf_getAttrib(df_var, R_NamesSymbol) != R_NilValue;
          if (has_names){
            SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
            SEXP new_names = Rf_protect(cpp_sset_unsafe(
              old_names, pi2, out_size, n_cores
            ));
            Rf_setAttrib(list_var, R_NamesSymbol, new_names);
          }
          SET_VECTOR_ELT(out, j, list_var);
          Rf_unprotect(1 + (has_names * 2));
        } else {
          SET_VECTOR_ELT(out, j, cheapr_sset(df_var, indices2));
        }
        Rf_unprotect(1);
      }
    } else {
      for (int j = 0; j < ncols; ++j){
        SEXP df_var = Rf_protect(p_x[j]);
        SET_VECTOR_ELT(out, j, cheapr_sset(df_var, indices));
        Rf_unprotect(1);
      }
    }
  }
  SEXP names = Rf_protect(Rf_duplicate(Rf_getAttrib(x, R_NamesSymbol)));
  ++n_protections;
  Rf_setAttrib(out, R_NamesSymbol, names);

  // list to data frame object
  SEXP df_str = Rf_protect(Rf_ScalarString(Rf_mkChar("data.frame")));
  ++n_protections;
  if (out_size > 0){
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 2));
    ++n_protections;
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -out_size;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  } else {
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 0));
    ++n_protections;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  }
  Rf_classgets(out, df_str);
  Rf_unprotect(n_protections);
  return out;
}

// SEXP cpp_sset(SEXP x, SEXP indices){
//   if (!Rf_isObject(x) && Rf_isNull(Rf_getAttrib(x, R_NamesSymbol)) && is_alt_compact_seq(indices)){
//     SEXP int_seq_data = Rf_protect(Rf_coerceVector(alt_data1(indices), INTSXP));
//     int size = INTEGER(int_seq_data)[0];
//     int from = INTEGER(int_seq_data)[1];
//     int by = INTEGER(int_seq_data)[2];
//     int to = from + (std::max(size - 1, 0) * by);
//     Rf_unprotect(1);
//     return cpp_sset_range(x, from, to, by);
//   }
//   int *pi = INTEGER(indices);
//   int xn = Rf_length(x);
//   int n = Rf_length(indices);
//   int n_protections = 0;
//   int zero_count = 0;
//   int pos_count = 0;
//   int oob_count = 0;
//   int na_count = 0;
//   int k = 0;
//   int out_size;
//   // bool do_parallel = n >= 10000;
//   // int n_cores = do_parallel ? num_cores() : 1;
//
//   // Counting the number of:
//   // Zeroes
//   // Out-of-bounds indices
//   // Positive indices
//   // From this we can also work out the number of negatives
//
// //   if (do_parallel){
// // #pragma omp parallel for simd num_threads(n_cores) reduction(+:zero_count,pos_count,oob_count,na_count)
//   //   for (int j = 0; j < n; ++j){
//   //     zero_count += (pi[j] == 0);
//   //     pos_count += (pi[j] > 0);
//   //     na_count += (pi[j] == NA_INTEGER);
//   //     oob_count += (std::fabs(pi[j]) > xn);
//   //   }
//   // } else {
// // #pragma omp for simd
//     for (int j = 0; j < n; ++j){
//       zero_count += (pi[j] == 0);
//       pos_count += (pi[j] > 0);
//       na_count += pi[j] == NA_INTEGER;
//       oob_count += (std::abs(pi[j]) > xn);
//     }
//   // }
//   int neg_count = n - pos_count - zero_count - na_count;
//   if ( pos_count > 0 && neg_count > 0){
//     Rf_error("Cannot mix positive and negative indices");
//   }
//   bool simple_sset = zero_count == 0 && oob_count == 0 && na_count == 0 && pos_count == n;
//
//   // Convert negative index vector to positive
//
//   if (neg_count > 0){
//     SEXP indices2 = Rf_protect(cpp11::package("cheapr")["neg_indices_to_pos"](indices, xn));
//     // ++n_protections;
//     // int *pi2 = INTEGER(indices2);
//     // pi = pi2;
//     // out_size = Rf_xlength(indices2);
//     // n = out_size;
//     // ++n_protections;
//     Rf_unprotect(n_protections + 1);
//     return cpp_sset_unsafe(x, indices2);
//   } else {
//     out_size = n - zero_count;
//   }
//   if (simple_sset){
//     Rf_unprotect(n_protections);
//     return cpp_sset_unsafe(x, indices);
//   }
//   switch ( TYPEOF(x) ){
//   case NILSXP: {
//     return R_NilValue;
//   }
//   case INTSXP: {
//     int *p_x = INTEGER(x);
//     SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
//     ++n_protections;
//     int *p_out = INTEGER(out);
//     for (R_xlen_t i = 0; i < n; ++i){
//       int xi = pi[i];
//       if (xi != 0){
//         p_out[k++] = (xi <= xn && xi != NA_INTEGER) ? p_x[xi - 1] : NA_INTEGER;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
// }


// A subset method using c++ vectors

// list cpp_sset(SEXP x, integers i){
//   int xn = Rf_xlength(x);
//   int n = i.size();
//   switch ( TYPEOF(x) ){
//   case INTSXP: {
//     std::vector<int> out;
//     int *p_x = INTEGER(x);
//     out.reserve(n);
//     for (int j = 0; j < n; ++j){
//       if (i[j] > 0 && i[j] <= xn){
//         int val = p_x[i[j] - 1];
//         out.push_back(val);
//       } else {
//         out.push_back(NA_INTEGER);
//       }
//     }
//       return writable::list({
//         "out"_nm = out
//       });
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
// }
