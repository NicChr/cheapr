#include "cheapr.h"
#include <R.h>
// #include <R_ext/Memory.h>

// Subsetting vectors and data frames
// Includes a unique optimisation on range subsetting

// OOB and NA is checked when `check = T` but zeros are never checked

SEXP cpp_sset_unsafe(SEXP x, SEXP indices, bool check){
  int NP = 0;
  int *pind = INTEGER(indices);
  int n = Rf_length(indices), xn = Rf_length(x);

  SEXP out;

  switch ( TYPEOF(x) ){

  case NILSXP: {
    out = R_NilValue;
    break;
  }
  case LGLSXP: {
    int *p_x = LOGICAL(x);
    out = Rf_protect(Rf_allocVector(LGLSXP, n)); ++NP;
    int *p_out = LOGICAL(out);

    if (check){
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = (pind[i] == NA_INTEGER || pind[i] > xn) ? NA_LOGICAL : p_x[pind[i] - 1];
      }
    } else {
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    break;
  }
  case INTSXP: {
    int *p_x = INTEGER(x);
    out = Rf_protect(Rf_allocVector(INTSXP, n)); ++NP;
    int *p_out = INTEGER(out);
    if (check){
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = (pind[i] == NA_INTEGER || pind[i] > xn) ? NA_INTEGER : p_x[pind[i] - 1];
      }
    } else {
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    break;
  }
  // case CHEAPR_INT64SXP: {
  //   long long int *p_x = INTEGER64_PTR(x);
  //   out = Rf_protect(Rf_allocVector(REALSXP, n)); ++NP;
  //   Rf_classgets(out, Rf_getAttrib(x, R_ClassSymbol));
  //   long long int *p_out = INTEGER64_PTR(out);
  //   if (check){
  //     OMP_FOR_SIMD
  //     for (int i = 0; i < n; ++i){
  //       p_out[i] = (pind[i] == NA_INTEGER || pind[i] > xn) ? NA_INTEGER64 : p_x[pind[i] - 1];
  //     }
  //   } else {
  //     OMP_FOR_SIMD
  //     for (int i = 0; i < n; ++i){
  //       p_out[i] = p_x[pind[i] - 1];
  //     }
  //   }
  //   break;
  // }
  case REALSXP: {
    double *p_x = REAL(x);
    out = Rf_protect(Rf_allocVector(REALSXP, n)); ++NP;
    double *p_out = REAL(out);
    if (check){
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = (pind[i] == NA_INTEGER || pind[i] > xn) ? NA_REAL : p_x[pind[i] - 1];
      }
    } else {
      OMP_FOR_SIMD
      for (int i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(STRSXP, n)); ++NP;

    if (check){
      for (int i = 0; i < n; ++i){
        SET_STRING_ELT(out, i, (pind[i] == NA_INTEGER || pind[i] > xn) ? NA_STRING : p_x[pind[i] - 1]);
      }
    } else {
      for (int i = 0; i < n; ++i) SET_STRING_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    out = Rf_protect(Rf_allocVector(CPLXSXP, n)); ++NP;
    Rcomplex *p_out = COMPLEX(out);
    if (check){
      for (int i = 0; i < n; ++i){
        if (pind[i] == NA_INTEGER || pind[i] > xn){
          p_out[i].r = NA_REAL;
          p_out[i].i = NA_REAL;
        } else {
          SET_COMPLEX_ELT(out, i, p_x[pind[i] - 1]);
        }
      }
    } else {
      for (int i = 0; i < n; ++i) SET_COMPLEX_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  case RAWSXP: {
    Rbyte *p_x = RAW(x);
    out = Rf_protect(Rf_allocVector(RAWSXP, n)); ++NP;
    if (check){
      for (int i = 0; i < n; ++i){
        SET_RAW_ELT(out, i, (pind[i] == NA_INTEGER || pind[i] > xn) ? 0 : p_x[pind[i] - 1]);
      }
    } else {
      for (int i = 0; i < n; ++i) SET_RAW_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  case VECSXP: {
    const SEXP *p_x = VECTOR_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(VECSXP, n)); ++NP;
    if (check){
      for (int i = 0; i < n; ++i){
        SET_VECTOR_ELT(out, i, (pind[i] == NA_INTEGER || pind[i] > xn) ? R_NilValue : p_x[pind[i] - 1]);
      }
    } else {
      for (int i = 0; i < n; ++i) SET_VECTOR_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  // if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){
  //   SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
  //   SEXP new_names = Rf_protect(cpp_sset_unsafe(old_names, indices, check)); ++NP;
  //   Rf_setAttrib(out, R_NamesSymbol, new_names);
  // }
  Rf_unprotect(NP);
  return out;

}

// A range-based subset method
// Can be readily used when indices are an altrep compact integer sequence
// Also works with negative-indexing


[[cpp11::register]]
SEXP cpp_sset_range(SEXP x, R_xlen_t from, R_xlen_t to, R_xlen_t by){
  int NP = 0;
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
  R_xlen_t out_size, istart1, istart2;
  R_xlen_t iend1 = 0, iend2 = 0;
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

  R_xlen_t k = 0;

  // Out-of-bounds
  R_xlen_t n_oob = std::max(( by > 0) ? iend - n : istart - n, (R_xlen_t) 0);
  // Adjustment for when all values are oob
  if ( ( by > 0 && istart > n ) || (by < 0 && iend > n)){
    n_oob = out_size;
  }
  R_xlen_t in_bounds_size = std::max(out_size - n_oob, (R_xlen_t) 0);

  SEXP out;

  switch ( TYPEOF(x) ){
  case NILSXP: {
    out = R_NilValue;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    out = Rf_protect(Rf_allocVector(TYPEOF(x), out_size)); ++NP;
    int *p_out = INTEGER(out);
    if (double_loop){
      memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * sizeof(int));
      memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * sizeof(int));
    } else {
      if (by > 0){
        memmove(p_out, p_x + (istart - 1), in_bounds_size * sizeof(int));
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n_oob; ++i) p_out[in_bounds_size + i] = NA_INTEGER;
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n_oob; ++i) p_out[i] = NA_INTEGER;
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i) p_out[istart - i - 1] = p_x[i];
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    out = Rf_protect(Rf_allocVector(REALSXP, out_size)); ++NP;
    double *p_out = REAL(out);
    if (double_loop){
      memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * sizeof(double));
      memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * sizeof(double));
    } else {
      if (by > 0){
        memmove(p_out, p_x + (istart - 1), in_bounds_size * sizeof(double));
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n_oob; ++i) p_out[in_bounds_size + i] = NA_REAL;
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n_oob; ++i) p_out[i] = NA_REAL;
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i) p_out[istart - i - 1] = p_x[i];
      }
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(STRSXP, out_size)); ++NP;
    if (double_loop){
      for (R_xlen_t i = istart1 - 1, k = 0; i < iend1; ++i, ++k){
        SET_STRING_ELT(out, k, p_x[i]);
      }
      for (R_xlen_t j = istart2 - 1, k = iend1; j < iend2; ++j, ++k){
        SET_STRING_ELT(out, k, p_x[j]);
      }
    } else {
      if (by > 0){
        for (R_xlen_t i = istart - 1, k = 0; i < (iend - n_oob); ++i, ++k){
          SET_STRING_ELT(out, k, p_x[i]);
        }
        for (R_xlen_t i = 0; i < n_oob; ++i){
          SET_STRING_ELT(out, in_bounds_size + i, NA_STRING);
        }
      } else {
        for (R_xlen_t i = 0; i < n_oob; ++i){
          SET_STRING_ELT(out, i, NA_STRING);
        }
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i){
          SET_STRING_ELT(out, istart - i - 1, p_x[i]);
        }
      }
    }
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    out = Rf_protect(Rf_allocVector(CPLXSXP, out_size)); ++NP;
    Rcomplex *p_out = COMPLEX(out);
    if (double_loop){
      memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * sizeof(Rcomplex));
      memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * sizeof(Rcomplex));
    } else {
      if (by > 0){
        memmove(p_out, &p_x[istart - 1], in_bounds_size * sizeof(Rcomplex));
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n_oob; ++i){
          R_xlen_t tempi = in_bounds_size + i;
          p_out[tempi].r = NA_REAL;
          p_out[tempi].i = NA_REAL;
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n_oob; ++i){
          p_out[i].r = NA_REAL;
          p_out[i].i = NA_REAL;
        }
        OMP_FOR_SIMD
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i){
          R_xlen_t tempi = istart - i - 1;
          p_out[tempi].r = p_x[i].r;
          p_out[tempi].i = p_x[i].i;
        }
      }
    }
    break;
  }
  case RAWSXP: {
    Rbyte *p_x = RAW(x);
    out = Rf_protect(Rf_allocVector(RAWSXP, out_size)); ++NP;
    if (double_loop){
      for (R_xlen_t i = istart1 - 1, k = 0; i < iend1; ++i, ++k){
        SET_RAW_ELT(out, k, p_x[i]);
      }
      for (R_xlen_t j = istart2 - 1, k = iend1; j < iend2; ++j, ++k){
        SET_RAW_ELT(out, k, p_x[j]);
      }
    } else {
      if (by > 0){
        for (R_xlen_t i = istart - 1; i < iend; ++i){
          SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
        }
      } else {
        for (R_xlen_t i = istart - 1; i >= iend - 1; --i){
          SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
        }
      }
    }
    break;
  }
  case VECSXP: {
    const SEXP *p_x = VECTOR_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(VECSXP, out_size)); ++NP;
    if (double_loop){
      for (R_xlen_t i = istart1 - 1, k = 0; i < iend1; ++i, ++k){
        SET_VECTOR_ELT(out, k, p_x[i]);
      }
      for (R_xlen_t j = istart2 - 1, k = iend1; j < iend2; ++j, ++k){
        SET_VECTOR_ELT(out, k, p_x[j]);
      }
    } else {
      if (by > 0){
        for (R_xlen_t i = istart - 1, k = 0; i < iend; ++i, ++k){
          if (i < n){
            SET_VECTOR_ELT(out, k, p_x[i]);
          }
        }
      } else {
        for (R_xlen_t i = istart - 1, k = 0; i >= iend - 1; --i, ++k){
          if (i < n){
            SET_VECTOR_ELT(out, k, p_x[i]);
          }
        }
      }
    }
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(NP);
  return out;
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

// expects only zero and negative elements in `exclude`

SEXP exclude_elements(SEXP x, SEXP exclude) {
  int n = Rf_length(x);
  int m = Rf_length(exclude);
  int *p_excl = INTEGER(exclude);
  int idx;

  // Which elements do we keep?
  int *p_keep = R_Calloc(n, int);

  OMP_FOR_SIMD
  for (int i = 0; i < n; ++i) p_keep[i] = TRUE;

  int exclude_count = 0;
  for (int j = 0; j < m; ++j) {
    if (p_excl[j] == NA_INTEGER) continue;
    if (p_excl[j] > 0){
      R_Free(p_keep);
      Rf_error("Cannot mix positive and negative subscripts");
    }
    idx = -p_excl[j];
    // Check keep array for already assigned FALSE to avoid double counting
    if (idx > 0 && idx <= n && p_keep[idx - 1] == TRUE){
      p_keep[idx - 1] = FALSE;
      ++exclude_count;
    }
  }

  SEXP out = Rf_protect(Rf_allocVector(INTSXP, n - exclude_count));
  int *p_out = INTEGER(out);
  int *p_x = INTEGER(x);
  int k = 0;
  for (int i = 0; i < n; ++i){
    if (p_keep[i]) p_out[k++] = p_x[i];
  }
  R_Free(p_keep);
  Rf_unprotect(1);
  return out;
}

// Cleans indices for subsetting
// Removes zeros unless indices is a compact seq because that
// is handled nicely by cpp_sset_range
// Because of this niche difference, do not export it to users

SEXP clean_indices(SEXP indices, int xn){
  long long int llxn = xn;
  int n = Rf_length(indices);
  int NP = 0;
  int zero_count = 0,
    pos_count = 0,
    oob_count = 0,
    na_count = 0;
  bool do_parallel = n >= CHEAPR_OMP_THRESHOLD;
  int n_cores = do_parallel ? num_cores() : 1;

  // If indices is a special type of ALTREP compact int sequence, we can
  // Use a range-based subset instead

  Rf_protect(indices = Rf_coerceVector(indices, INTSXP)); ++NP;

  int out_size;
  bool check_indices;
  SEXP clean_indices;

  if (is_compact_seq(indices)){
    clean_indices = indices;

    SEXP seq_data = Rf_protect(compact_seq_data(indices)); ++NP;
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    out_size = get_alt_final_sset_size(xn, from, to, by);
    check_indices = true;

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
        // oob_count counts NA as true so adjust after the fact
        oob_count += (std::llabs(pi[j]) > llxn);
        na_count += (pi[j] == NA_INTEGER);
      }
    } else {
      OMP_FOR_SIMD
      for (int j = 0; j < n; ++j){
        zero_count += (pi[j] == 0);
        pos_count += (pi[j] > 0);
        // oob_count counts NA as true so adjust after the fact
        oob_count += (std::llabs(pi[j]) > llxn);
        na_count += (pi[j] == NA_INTEGER);
      }
    }
    // adjust oob_count
    oob_count = oob_count - na_count;
    int neg_count = n - pos_count - zero_count - na_count;
    if ( (pos_count > 0 && neg_count > 0) ||
         (neg_count > 0 && na_count > 0)){
      Rf_error("Cannot mix positive and negative indices");
    }

    // Should a simplified sset method be used?

    check_indices = !(oob_count == 0 && na_count == 0);

    if (neg_count > 0){
      SEXP row_seq = Rf_protect(cpp_seq_len(xn)); ++NP;
      clean_indices = Rf_protect(exclude_elements(row_seq, indices)); ++NP;
      check_indices = false;
    } else if (zero_count > 0){
      SEXP zero = Rf_protect(Rf_ScalarInteger(0)); ++NP;
      clean_indices = Rf_protect(cpp_val_remove(indices, zero)); ++NP;
    } else {
      clean_indices = indices;
    }
    out_size = Rf_length(clean_indices);
  }

  SEXP out = Rf_protect(Rf_allocVector(VECSXP, 3)); ++NP;

  // There are the `Rf_Scalar` shortcuts BUT R crashes sometimes when
  // using the scalar logical shortcuts so I avoid it

  SEXP out_size_sexp = Rf_protect(Rf_allocVector(INTSXP, 1)); ++NP;
  SEXP check_indices_sexp = Rf_protect(Rf_allocVector(LGLSXP, 1)); ++NP;
  INTEGER(out_size_sexp)[0] = out_size;
  LOGICAL(check_indices_sexp)[0] = check_indices;
  SET_VECTOR_ELT(out, 0, clean_indices);
  SET_VECTOR_ELT(out, 1, out_size_sexp);
  SET_VECTOR_ELT(out, 2, check_indices_sexp);

  Rf_unprotect(NP);
  return out;
}


// Data frame subsetting

// Fast col select
// Supports
//  integer locations
//  character vectors
//  negative subscripting
//  NULL to signify all locs (shallow copy)

[[cpp11::register]]
SEXP cpp_df_select(SEXP x, SEXP locs){
  int NP = 0,
    n_cols = Rf_length(x),
    n_rows = Rf_length(Rf_getAttrib(x, R_RowNamesSymbol)),
    n_locs = Rf_length(locs);

  // Flag to check indices
  bool check = true;

  SEXP names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol)); ++NP;

  SEXP cols;

  if (Rf_isNull(locs)){
    // If NULL then select all cols
    cols = Rf_protect(cpp_seq_len(n_cols)); ++NP;
    n_locs = n_cols;
    check = false;
  } else if (Rf_isString(locs)){
    cols = Rf_protect(base_match(locs, names)); ++NP;
  } else if (Rf_isLogical(locs)){
    // If logical then find locs using `which_()`
    if (Rf_length(locs) != n_cols){
      Rf_unprotect(NP);
      Rf_error("`length(j)` must match `ncol(x)` when `j` is a logical vector");
    }
    cols = Rf_protect(cpp_which_(locs, false)); ++NP;
    n_locs = Rf_length(cols);
    check = false;
  } else {
    // Catch-all make sure cols is an int vector
    cols = Rf_protect(Rf_coerceVector(locs, INTSXP)); ++NP;
  }

  // Negative subscripting
  if (n_locs > 0 && INTEGER(cols)[0] != NA_INTEGER && INTEGER(cols)[0] < 0){
    Rf_protect(cols = exclude_elements(cpp_seq_len(n_cols), cols)); ++NP;
    n_locs = Rf_length(cols);
    check = false;
  }

  SEXP out = Rf_protect(Rf_allocVector(VECSXP, n_locs)); ++NP;
  SEXP out_names = Rf_protect(Rf_allocVector(STRSXP, n_locs)); ++NP;

  int *p_cols = INTEGER(cols);
  const SEXP *p_x = VECTOR_PTR_RO(x);
  const SEXP *p_names = STRING_PTR_RO(names);
  int k = 0;
  int col;

  if (check){
    for (int i = 0; i < n_locs; ++i) {
      col = p_cols[i];
      if (col != NA_INTEGER && col >= 1 && col <= n_cols){
        SET_VECTOR_ELT(out, k, p_x[col - 1]);
        SET_STRING_ELT(out_names, k, p_names[col - 1]);
        ++k;
      } else if (col != NA_INTEGER && col < 0){
        // This can only happen when there is a mix of pos & neg
        // but wasn't captured by the negative subscripting section
        // because that only looks to see if the 1st element is neg
        Rf_error("Cannot mix positive and negative subscripts");
      }
    }
  } else {
    for (int i = 0; i < n_locs; ++i) {
      col = p_cols[i];
      SET_VECTOR_ELT(out, i, p_x[col - 1]);
      SET_STRING_ELT(out_names, i, p_names[col - 1]);
    }
  }

  if (check && k != n_locs){
    Rf_protect(out = Rf_xlengthgets(out, k)); ++NP;
    Rf_protect(out_names = Rf_xlengthgets(out_names, k)); ++NP;
  }

  // Make a plain data frame
  SEXP row_names;
  if (n_rows > 0){
    row_names = Rf_protect(Rf_allocVector(INTSXP, 2)); ++NP;
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -n_rows;
  } else {
    row_names = Rf_protect(Rf_allocVector(INTSXP, 0)); ++NP;
  }
  Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  Rf_classgets(out, Rf_mkString("data.frame"));
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  Rf_unprotect(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_df_slice(SEXP x, SEXP indices){

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame`, not a %s", Rf_type2char(TYPEOF(x)));
  }

  if (Rf_isNull(indices)){
    return x;
  }

  int xn = cpp_df_nrow(x);
  int ncols = Rf_length(x);
  int NP = 0;
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, ncols)); ++NP;

  const SEXP clean_indices_info = Rf_protect(clean_indices(indices, xn)); ++NP;
  Rf_protect(indices = VECTOR_ELT(clean_indices_info, 0)); ++NP;
  int out_size = INTEGER(VECTOR_ELT(clean_indices_info, 1))[0];
  bool check_indices = LOGICAL(VECTOR_ELT(clean_indices_info, 2))[0];

  // If indices is a special type of ALTREP compact int sequence, we can
  // Use a range-based subset instead

  if (is_compact_seq(indices)){

    // ALTREP integer sequence method

    SEXP seq_data = Rf_protect(compact_seq_data(indices)); ++NP;
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    for (int j = 0; j < ncols; ++j){
      SEXP df_var = Rf_protect(p_x[j]);
      if (is_base_atomic_vec(df_var)){
        SEXP list_var = Rf_protect(cpp_sset_range(df_var, from, to, by));
        Rf_copyMostAttrib(df_var, list_var);
        int has_names = !Rf_isNull(Rf_getAttrib(df_var, R_NamesSymbol));
        if (has_names){
          SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
          SEXP new_names = Rf_protect(cpp_sset_range(
            old_names, from, to, by)
          );
          Rf_setAttrib(list_var, R_NamesSymbol, new_names);
        }
        SET_VECTOR_ELT(out, j, list_var);
        Rf_unprotect(1 + (has_names * 2));
      } else {
        SET_VECTOR_ELT(out, j, cheapr_sset(df_var, indices));
      }
      // Unprotecting new data frame variable
      Rf_unprotect(1);
    }
  } else {
    // If Index vector is clean we can use fast subset

      for (int j = 0; j < ncols; ++j){
        SEXP df_var = Rf_protect(p_x[j]);
        if (is_base_atomic_vec(df_var)){
          SEXP list_var = Rf_protect(cpp_sset_unsafe(df_var, indices, check_indices));
          Rf_copyMostAttrib(df_var, list_var);
          int has_names = !Rf_isNull(Rf_getAttrib(df_var, R_NamesSymbol));
          if (has_names){
            SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
            SEXP new_names = Rf_protect(cpp_sset_unsafe(
              old_names, indices, check_indices
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
    }

  cpp_copy_names(x, out, true);

  // list to data frame object
  if (out_size > 0){
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 2)); ++NP;
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -out_size;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  } else {
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 0)); ++NP;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  }
  Rf_classgets(out, Rf_mkString("data.frame"));
  Rf_unprotect(NP);
  return out;
}

// Subset that does both selecting and slicing

[[cpp11::register]]
SEXP cpp_df_subset(SEXP x, SEXP i, SEXP j, bool keep_attrs){

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame`, not a %s", Rf_type2char(TYPEOF(x)));
  }

  int NP = 0, n_rows = cpp_df_nrow(x);

  if (Rf_isLogical(i)){
    if (Rf_length(i) != n_rows){
      Rf_error("`length(i)` must match `nrow(x)` when `i` is a logical vector");
    }
    Rf_protect(i = cpp_which_(i, false)); ++NP;
  }

  // Subset columns

  // `cpp_df_select()` always creates a shallow copy

  SEXP out = Rf_protect(cpp_df_select(x, j)); ++NP;
  Rf_protect(out = cpp_df_slice(out, i)); ++NP;

  if (keep_attrs){
    SEXP names = Rf_protect(Rf_getAttrib(out, R_NamesSymbol)); ++NP;
    SEXP row_names = Rf_protect(Rf_getAttrib(out, R_RowNamesSymbol)); ++NP;
    Rf_copyMostAttrib(x, out);
    Rf_setAttrib(out, R_NamesSymbol, names);
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  }
  Rf_unprotect(NP);
  return out;
}

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


// This is kept to not break dependencies, will remove later

[[cpp11::register]]
SEXP cpp_sset_df(SEXP x, SEXP indices){
  int xn = cpp_df_nrow(x);
  int ncols = Rf_length(x);
  int NP = 0;
  // cpp11::function cheapr_sset = cpp11::package("cheapr")["sset"];
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, ncols)); ++NP;

  const SEXP clean_indices_info = Rf_protect(clean_indices(indices, xn)); ++NP;
  Rf_protect(indices = VECTOR_ELT(clean_indices_info, 0)); ++NP;
  int out_size = INTEGER(VECTOR_ELT(clean_indices_info, 1))[0];
  bool check_indices = LOGICAL(VECTOR_ELT(clean_indices_info, 2))[0];

  // If indices is a special type of ALTREP compact int sequence, we can
  // Use a range-based subset instead

  if (is_compact_seq(indices)){

    // ALTREP integer sequence method

    SEXP seq_data = Rf_protect(compact_seq_data(indices)); ++NP;
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    for (int j = 0; j < ncols; ++j){
      SEXP df_var = Rf_protect(p_x[j]);
      if (is_base_atomic_vec(df_var)){
        SEXP list_var = Rf_protect(cpp_sset_range(df_var, from, to, by));
        Rf_copyMostAttrib(df_var, list_var);
        int has_names = !Rf_isNull(Rf_getAttrib(df_var, R_NamesSymbol));
        if (has_names){
          SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
          SEXP new_names = Rf_protect(cpp_sset_range(
            old_names, from, to, by)
          );
          Rf_setAttrib(list_var, R_NamesSymbol, new_names);
        }
        SET_VECTOR_ELT(out, j, list_var);
        Rf_unprotect(1 + (has_names * 2));
      } else {
        SET_VECTOR_ELT(out, j, cheapr_sset(df_var, indices));
      }
      // Unprotecting new data frame variable
      Rf_unprotect(1);
    }
  } else {
    // If Index vector is clean we can use fast subset

    for (int j = 0; j < ncols; ++j){
      SEXP df_var = Rf_protect(p_x[j]);
      if (is_base_atomic_vec(df_var)){
        SEXP list_var = Rf_protect(cpp_sset_unsafe(df_var, indices, check_indices));
        Rf_copyMostAttrib(df_var, list_var);
        int has_names = !Rf_isNull(Rf_getAttrib(df_var, R_NamesSymbol));
        if (has_names){
          SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
          SEXP new_names = Rf_protect(cpp_sset_unsafe(
            old_names, indices, check_indices
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
  }

  cpp_copy_names(x, out, false);

  // list to data frame object
  if (out_size > 0){
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 2)); ++NP;
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -out_size;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  } else {
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 0)); ++NP;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  }
  Rf_classgets(out, Rf_mkString("data.frame"));
  Rf_unprotect(NP);
  return out;
}
