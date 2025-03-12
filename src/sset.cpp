#include "cheapr.h"
#include <R.h> // R_Calloc

// Subsetting vectors and data frames
// Includes a unique optimisation on range subsetting

// Author: Nick Christofides

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

// A cheaper negative subscript for integer vectors
// expects only zero-value and negative elements in `exclude`

SEXP exclude_locs(SEXP exclude, R_xlen_t xn) {

  R_xlen_t n = xn;
  R_xlen_t m = Rf_length(exclude);
  R_xlen_t out_size, idx;
  R_xlen_t exclude_count = 0;
  R_xlen_t i = 0, k = 0;

  // Which elements do we keep?
  int *p_keep = R_Calloc(n, int);

  OMP_FOR_SIMD
  for (int i = 0; i < n; ++i) p_keep[i] = TRUE;

  SEXP x_seq = Rf_protect(cpp_seq_len(xn));

  if (xn > integer_max_){
    Rf_protect(exclude = Rf_coerceVector(exclude, REALSXP));
    double *p_x = REAL(x_seq);
    double *p_excl = REAL(exclude);

    for (int j = 0; j < m; ++j) {
      if (p_excl[j] != p_excl[j]) continue;
      if (p_excl[j] > 0){
        R_Free(p_keep);
        Rf_unprotect(2);
        Rf_error("Cannot mix positive and negative subscripts");
      }
      idx = -p_excl[j];
      // Check keep array for already assigned FALSE to avoid double counting
      if (idx > 0 && idx <= n && p_keep[idx - 1] == TRUE){
        p_keep[idx - 1] = FALSE;
        ++exclude_count;
      }
    }
    out_size = n - exclude_count;
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
    double *p_out = REAL(out);
    while(k != out_size){
      if (p_keep[i]){
        p_out[k++] = p_x[i];
      }
      ++i;
    }
    R_Free(p_keep);
    Rf_unprotect(3);
    return out;
  } else {
    int *p_x = INTEGER(x_seq);
    int *p_excl = INTEGER(exclude);

    for (int j = 0; j < m; ++j) {
      if (p_excl[j] == NA_INTEGER) continue;
      if (p_excl[j] > 0){
        R_Free(p_keep);
        Rf_unprotect(1);
        Rf_error("Cannot mix positive and negative subscripts");
      }
      idx = -p_excl[j];
      // Check keep array for already assigned FALSE to avoid double counting
      if (idx > 0 && idx <= n && p_keep[idx - 1] == TRUE){
        p_keep[idx - 1] = FALSE;
        ++exclude_count;
      }
    }
    out_size = n - exclude_count;
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
    int *p_out = INTEGER(out);
    while(k != out_size){
      if (p_keep[i]){
        p_out[k++] = p_x[i];
      }
      ++i;
    }
    R_Free(p_keep);
    Rf_unprotect(2);
    return out;
  }
}

// Cleans indices for subsetting
// Also returns metadata regarding final vec size and if indices should be
// checked (internal flag)

[[cpp11::register]]
SEXP clean_indices(SEXP indices, R_xlen_t xn){
  int NP = 0;
  R_xlen_t zero_count = 0,
    pos_count = 0,
    oob_count = 0,
    na_count = 0,
    neg_count = 0;

  R_xlen_t n = Rf_xlength(indices);
  bool do_parallel = n >= CHEAPR_OMP_THRESHOLD;
  int n_cores = do_parallel ? num_cores() : 1;

  R_xlen_t out_size;
  bool check_indices = true;
  SEXP clean_indices;

  if (is_compact_seq(indices)){
    clean_indices = indices;

    SEXP seq_data = Rf_protect(compact_seq_data(indices)); ++NP;
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    out_size = get_alt_final_sset_size(xn, from, to, by);
    check_indices = true;

  } else if (Rf_isLogical(indices)){
    if (Rf_length(indices) != xn){
      Rf_error("`length(i)` must match `length(x)` when `i` is a logical vector");
    }
    clean_indices = Rf_protect(cpp_which_(indices, false)); ++NP;
    n = Rf_xlength(clean_indices);
    out_size = n;
    check_indices = false;
  } else if (xn > integer_max_){
    Rf_protect(indices = Rf_coerceVector(indices, REALSXP)); ++NP;

    double *pi = REAL(indices);

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
        na_count += (pi[j] != pi[j]);
      }
    } else {
      OMP_FOR_SIMD
      for (int j = 0; j < n; ++j){
        zero_count += (pi[j] == 0);
        pos_count += (pi[j] > 0);
        oob_count += (std::fabs(pi[j]) > xn);
        na_count += (pi[j] != pi[j]);
      }
    }
    neg_count = n - pos_count - zero_count - na_count;
    if ( (pos_count > 0 && neg_count > 0) ||
         (neg_count > 0 && na_count > 0)){
      Rf_error("Cannot mix positive and negative indices");
    }

    // Should a simplified sset method be used?

    check_indices = !(oob_count == 0 && na_count == 0 && zero_count == 0);

    if (neg_count > 0){
      clean_indices = Rf_protect(exclude_locs(indices, xn)); ++NP;
      check_indices = false;
      out_size = Rf_length(clean_indices);
    } else {
      clean_indices = indices;
      out_size = pos_count + na_count;
    }
  } else {

    Rf_protect(indices = Rf_coerceVector(indices, INTSXP)); ++NP;

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
        oob_count += (std::llabs(pi[j]) > xn);
        na_count += (pi[j] == NA_INTEGER);
      }
    } else {
      OMP_FOR_SIMD
      for (int j = 0; j < n; ++j){
        zero_count += (pi[j] == 0);
        pos_count += (pi[j] > 0);
        // oob_count counts NA as true so adjust after the fact
        oob_count += (std::llabs(pi[j]) > xn);
        na_count += (pi[j] == NA_INTEGER);
      }
    }
    // adjust oob_count
    oob_count = oob_count - na_count;
    neg_count = n - pos_count - zero_count - na_count;
    if ( (pos_count > 0 && neg_count > 0) ||
         (neg_count > 0 && na_count > 0)){
      Rf_error("Cannot mix positive and negative indices");
    }

    // Should a simplified sset method be used?

    check_indices = !(oob_count == 0 && na_count == 0 && zero_count == 0);

    if (neg_count > 0){
      clean_indices = Rf_protect(exclude_locs(indices, xn)); ++NP;
      check_indices = false;
      out_size = Rf_length(clean_indices);
    } else {
      clean_indices = indices;
      out_size = pos_count + na_count;
    }
  }

  SEXP out = Rf_protect(Rf_allocVector(VECSXP, 3)); ++NP;

  // There are the `Rf_Scalar` shortcuts BUT R crashes sometimes when
  // using the scalar logical shortcuts so I avoid it

  SEXP out_size_sexp = Rf_protect(Rf_allocVector(REALSXP, 1)); ++NP;
  SEXP check_indices_sexp = Rf_protect(Rf_allocVector(LGLSXP, 1)); ++NP;
  REAL(out_size_sexp)[0] = out_size;
  LOGICAL(check_indices_sexp)[0] = check_indices;
  SET_VECTOR_ELT(out, 0, clean_indices);
  SET_VECTOR_ELT(out, 1, out_size_sexp);
  SET_VECTOR_ELT(out, 2, check_indices_sexp);

  Rf_unprotect(NP);
  return out;
}


// A range-based subset method
// Can be readily used when indices are an altrep compact integer sequence
// Also works with negative-indexing

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

// Vector subset
// OOB, zeros and NA are checked when `check = T`

SEXP sset_vec(SEXP x, SEXP indices, bool check){

  if (is_compact_seq(indices)){
    SEXP seq_data = Rf_protect(compact_seq_data(indices));
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    SEXP out = Rf_protect(cpp_sset_range(x, from, to, by));
    Rf_unprotect(2);
    return out;
  }

  int NP = 0;
  R_xlen_t
  n = Rf_xlength(indices),
    out_size = n,
    k = 0,
    j,
    xn = Rf_xlength(x),
    na_int = NA_INTEGER;

  SEXP out;

  if (xn > integer_max_){

    double *pind = REAL(indices);

    switch ( TYPEOF(x) ){

    case NILSXP: {
      out = R_NilValue;
      break;
    }
    case LGLSXP:
    case INTSXP: {
      int *p_x = INTEGER(x);
      out = Rf_protect(Rf_allocVector(TYPEOF(x), n)); ++NP;
      int *p_out = INTEGER(out);

      if (check){
        for (R_xlen_t i = 0; i < n; ++i){
          j = pind[i];
          if (j != 0){
            p_out[k++] = (pind[i] != pind[i] || j > xn) ? NA_INTEGER : p_x[j - 1];
          } else {
            --out_size;
          }
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n; ++i){
          p_out[i] = p_x[(R_xlen_t)pind[i] - 1];
        }
      }
      break;
    }
    case REALSXP: {
      double *p_x = REAL(x);
      out = Rf_protect(Rf_allocVector(REALSXP, n)); ++NP;
      double *p_out = REAL(out);
      if (check){
        for (R_xlen_t i = 0; i < n; ++i){
          j = pind[i];
          if (j != 0){
            p_out[k++] = (pind[i] != pind[i] || j > xn) ? NA_REAL : p_x[j - 1];
          } else {
            --out_size;
          }
        }
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < n; ++i){
          p_out[i] = p_x[(R_xlen_t)pind[i] - 1];
        }
      }
      break;
    }
    case STRSXP: {
      const SEXP *p_x = STRING_PTR_RO(x);
      out = Rf_protect(Rf_allocVector(STRSXP, n)); ++NP;

      if (check){
        for (R_xlen_t i = 0; i < n; ++i){
          j = pind[i];
          if (j != 0){
            SET_STRING_ELT(out, k++, (pind[i] != pind[i] || j > xn) ? NA_STRING : p_x[j - 1]);
          } else {
            --out_size;
          }
        }
      } else {
        for (R_xlen_t i = 0; i < n; ++i) SET_STRING_ELT(out, i, p_x[(R_xlen_t)pind[i] - 1]);
      }
      break;
    }
    case CPLXSXP: {
      Rcomplex *p_x = COMPLEX(x);
      out = Rf_protect(Rf_allocVector(CPLXSXP, n)); ++NP;
      Rcomplex *p_out = COMPLEX(out);
      if (check){
        for (R_xlen_t i = 0; i < n; ++i){

          j = pind[i];

          if (j != 0){
            if (pind[i] != pind[i] || j > xn){
              p_out[k].r = NA_REAL;
              p_out[k].i = NA_REAL;
            } else {
              SET_COMPLEX_ELT(out, k, p_x[j - 1]);
            }
            ++k;
          } else {
            --out_size;
          }

        }
      } else {
        for (R_xlen_t i = 0; i < n; ++i) SET_COMPLEX_ELT(out, i, p_x[(R_xlen_t)pind[i] - 1]);
      }
      break;
    }
    case RAWSXP: {
      Rbyte *p_x = RAW(x);
      out = Rf_protect(Rf_allocVector(RAWSXP, n)); ++NP;
      if (check){
        for (R_xlen_t i = 0; i < n; ++i){
          j = pind[i];
          if (j != 0){
            SET_RAW_ELT(out, k++, (pind[i] != pind[i] || j > xn) ? 0 : p_x[j - 1]);
          } else {
            --out_size;
          }
        }
      } else {
        for (R_xlen_t i = 0; i < n; ++i) SET_RAW_ELT(out, i, p_x[(R_xlen_t)pind[i] - 1]);
      }
      break;
    }
    case VECSXP: {
      const SEXP *p_x = VECTOR_PTR_RO(x);
      out = Rf_protect(Rf_allocVector(VECSXP, n)); ++NP;
      if (check){
        for (R_xlen_t i = 0; i < n; ++i){
          j = pind[i];
          if (j != 0){
            SET_VECTOR_ELT(out, k++, (pind[i] != pind[i] || j > xn) ? R_NilValue : p_x[j - 1]);
          } else {
            --out_size;
          }
        }
      } else {
        for (R_xlen_t i = 0; i < n; ++i) SET_VECTOR_ELT(out, i, p_x[(R_xlen_t)pind[i] - 1]);
      }
      break;
    }
    default: {
      Rf_unprotect(NP);
      Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
    }
    }
  } else {

  int *pind = INTEGER(indices);
  switch ( TYPEOF(x) ){

  case NILSXP: {
    out = R_NilValue;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);
    out = Rf_protect(Rf_allocVector(TYPEOF(x), n)); ++NP;
    int *p_out = INTEGER(out);

    if (check){
      for (R_xlen_t i = 0; i < n; ++i){
        j = pind[i];
        if (j != 0){
          p_out[k++] = (j == na_int || j > xn) ? NA_INTEGER : p_x[j - 1];
        } else {
          --out_size;
        }
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    out = Rf_protect(Rf_allocVector(REALSXP, n)); ++NP;
    double *p_out = REAL(out);
    if (check){
      for (R_xlen_t i = 0; i < n; ++i){
        j = pind[i];
        if (j != 0){
          p_out[k++] = (j == na_int || j > xn) ? NA_REAL : p_x[j - 1];
        } else {
          --out_size;
        }
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        p_out[i] = p_x[pind[i] - 1];
      }
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_x = STRING_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(STRSXP, n)); ++NP;

    if (check){
      for (R_xlen_t i = 0; i < n; ++i){
        j = pind[i];
        if (j != 0){
          SET_STRING_ELT(out, k++, (j == na_int || j > xn) ? NA_STRING : p_x[j - 1]);
        } else {
          --out_size;
        }
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i) SET_STRING_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    out = Rf_protect(Rf_allocVector(CPLXSXP, n)); ++NP;
    Rcomplex *p_out = COMPLEX(out);
    if (check){
      for (R_xlen_t i = 0; i < n; ++i){

        j = pind[i];

        if (j != 0){
          if (j == na_int || j > xn){
            p_out[k].r = NA_REAL;
            p_out[k].i = NA_REAL;
          } else {
            SET_COMPLEX_ELT(out, k, p_x[j - 1]);
          }
          ++k;
        } else {
          --out_size;
        }

      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i) SET_COMPLEX_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  case RAWSXP: {
    Rbyte *p_x = RAW(x);
    out = Rf_protect(Rf_allocVector(RAWSXP, n)); ++NP;
    if (check){
      for (R_xlen_t i = 0; i < n; ++i){
        j = pind[i];
        if (j != 0){
          SET_RAW_ELT(out, k++, (j == na_int || j > xn) ? 0 : p_x[j - 1]);
        } else {
          --out_size;
        }
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i) SET_RAW_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  case VECSXP: {
    const SEXP *p_x = VECTOR_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(VECSXP, n)); ++NP;
    if (check){
      for (R_xlen_t i = 0; i < n; ++i){
        j = pind[i];
        if (j != 0){
          SET_VECTOR_ELT(out, k++, (j == na_int || j > xn) ? R_NilValue : p_x[j - 1]);
        } else {
          --out_size;
        }
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i) SET_VECTOR_ELT(out, i, p_x[pind[i] - 1]);
    }
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  }
  if (!Rf_isNull(out) && out_size != n){
    Rf_protect(out = Rf_xlengthgets(out, out_size)); ++NP;
  }
  Rf_unprotect(NP);
  return out;
}

// Basic and fast vector subset, exported to R

[[cpp11::register]]
SEXP cpp_sset(SEXP x, SEXP indices){
  int NP = 0;

  SEXP indices_metadata = Rf_protect(clean_indices(indices, Rf_xlength(x))); ++NP;
  SEXP clean_indices = Rf_protect(VECTOR_ELT(indices_metadata, 0)); ++NP;
  bool check_indices = LOGICAL(VECTOR_ELT(indices_metadata, 2))[0];
  SEXP out = Rf_protect(sset_vec(x, clean_indices, check_indices)); ++NP;

  Rf_copyMostAttrib(x, out);

  // Subset names

  SEXP xnames = Rf_protect(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
  if (!Rf_isNull(xnames)){
    SEXP outnames = Rf_protect(sset_vec(xnames, clean_indices, check_indices)); ++NP;
    Rf_setAttrib(out, R_NamesSymbol, outnames);
  }
  Rf_unprotect(NP);
  return out;
}

// Safe and very efficient reverse in-place (or with copy)

// matrix structure is preserved with cpp_rev

// For `set = F` the same can be accomplished (slightly faster) with
// the range-based subset
// Keeping this implementation as it works seamlessly for in-place and
// normal reverse

[[cpp11::register]]
SEXP cpp_rev(SEXP x, bool set){
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t half = int_div(n, 2);
  R_xlen_t n2 = n - 1; // Offset n for 0-indexing
  R_xlen_t k;
  int NP = 0;
  bool rev_names = true;
  bool is_altrep = ALTREP(x);
  if (set && is_altrep){
    Rf_warning("Cannot update an ALTREP by reference, a copy has been made.\n\tEnsure the result is assigned to an object if used in further calculations");
  }
  Rf_protect(x = altrep_materialise(x)); ++NP;

  // altrep will have already been materialised so this should be safe
  // and avoids a second copy
  if (is_altrep){
    set = true;
  }
  SEXP out;
  switch (TYPEOF(x)){
  case NILSXP: {
    out = R_NilValue;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    int *p_out = INTEGER(out);
    int left;
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      left = p_out[i];
      p_out[i] = p_out[k];
      p_out[k] = left;
    }
    break;
  }
  case REALSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    double *p_out = REAL(out);
    double left;
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      left = p_out[i];
      p_out[i] = p_out[k];
      p_out[k] = left;
    }
    break;
  }
  case STRSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    const SEXP *p_out = STRING_PTR_RO(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      SEXP left = Rf_protect(p_out[i]);
      SET_STRING_ELT(out, i, p_out[k]);
      SET_STRING_ELT(out, k, left);
      Rf_unprotect(1);
    }
    break;
  }
  case CPLXSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    Rcomplex *p_out = COMPLEX(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      Rcomplex left = p_out[i];
      SET_COMPLEX_ELT(out, i, p_out[k]);
      SET_COMPLEX_ELT(out, k, left);
    }
    break;
  }
  case RAWSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    Rbyte *p_out = RAW(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      Rbyte left = p_out[i];
      SET_RAW_ELT(out, i, p_out[k]);
      SET_RAW_ELT(out, k, left);
    }
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
    // We can reverse data frames in-place with the below commented-out code
    // It will not work properly though for data.tables

    // case VECSXP: {
    //   if (recursive){
    //   rev_names = false;
    //   out = Rf_protect(Rf_allocVector(VECSXP, n)); ++NP;
    //   SHALLOW_DUPLICATE_ATTRIB(out, x);
    //   for (R_xlen_t i = 0; i < n; ++i){
    //     SET_VECTOR_ELT(out, i, cpp_rev(VECTOR_ELT(x, i), true, set));
    //   }
    //   break;
    // } else if (!recursive && (!Rf_isObject(x) || is_df(x))){
    //   out = Rf_protect(set ? x : list_shallow_copy(x, false)); ++NP;
    //   if (!set){
    //     // SHALLOW_DUPLICATE_ATTRIB(out, x);
    //     cpp_copy_names(x, out, true);
    //   }
    //   const SEXP *p_out = VECTOR_PTR_RO(out);
    //   for (R_xlen_t i = 0; i < half; ++i) {
    //     k = n2 - i;
    //     SEXP left = Rf_protect(p_out[i]);
    //     SET_VECTOR_ELT(out, i, p_out[k]);
    //     SET_VECTOR_ELT(out, k, left);
    //     Rf_unprotect(1);
    //   }
    //     break;
    //   }
    // }
    // default: {
    //   // set rev_names to false because catch-all r_rev() will handle that
    //   rev_names = false;
    //   cpp11::function r_rev = cpp11::package("base")["rev"];
    //   if (set){
    //     Rf_unprotect(NP);
    //     Rf_error("Can't reverse in place here");
    //   }
    //   out = Rf_protect(r_rev(x)); ++NP;
    //   break;
    // }
  }
  // // If x has names, reverse them too
  if (rev_names && !Rf_isNull(Rf_getAttrib(out, R_NamesSymbol))){
    SEXP old_names = Rf_protect(Rf_getAttrib(out, R_NamesSymbol)); ++NP;
    // should be okay to rev in-place here because we already copied the names
    Rf_setAttrib(out, R_NamesSymbol, cpp_rev(old_names, true));
  }
  Rf_unprotect(NP);
  return out;
}

// SEXP cpp_cheapr_rev(SEXP x){
//   SEXP indices = Rf_protect(base_colon(vec_length(x), 0));
//   SEXP out;
//   if (is_simple_atomic_vec(x)){
//     out = Rf_protect(cpp_sset(x, indices));
//   } else {
//     out = Rf_protect(cheapr_sset(x, indices));
//   }
//   Rf_unprotect(2);
//   return out;
// }


// Data frame subsetting

// Fast col select
// Supports
//  integer locations
//  character vectors
//  negative subscripting
//  NULL to signify all locs (shallow copy)

[[cpp11::register]]
SEXP cpp_df_select(SEXP x, SEXP locs){

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame`, not a %s", Rf_type2char(TYPEOF(x)));
  }

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
    Rf_protect(cols = exclude_locs(cols, n_cols)); ++NP;
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
    Rf_protect(out = Rf_lengthgets(out, k)); ++NP;
    Rf_protect(out_names = Rf_lengthgets(out_names, k)); ++NP;
  }

  // Make a plain data frame
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(n_rows));
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


  // Clean indices and get metadata

  const SEXP clean_indices_info = Rf_protect(clean_indices(indices, xn)); ++NP;
  Rf_protect(indices = VECTOR_ELT(clean_indices_info, 0)); ++NP;
  int out_size = REAL(VECTOR_ELT(clean_indices_info, 1))[0];
  bool check_indices = LOGICAL(VECTOR_ELT(clean_indices_info, 2))[0];


  // Subset columns

  for (int j = 0; j < ncols; ++j){
    SEXP df_var = Rf_protect(p_x[j]);
    SEXP names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
    SEXP list_var;
    if (is_simple_atomic_vec(df_var)){
      list_var = Rf_protect(sset_vec(df_var, indices, check_indices));
      Rf_copyMostAttrib(df_var, list_var);
      Rf_setAttrib(list_var, R_NamesSymbol, sset_vec(names, indices, check_indices));
    } else {
      list_var = Rf_protect(cheapr_sset(df_var, indices));
    }
    SET_VECTOR_ELT(out, j, list_var);
    Rf_unprotect(3); // Unprotect var, names and new var
  }

  cpp_copy_names(x, out, false);

  // list to data frame object
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(out_size));
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

  int NP = 0;

  // Subset columns
  // `cpp_df_select()` always creates a shallow copy
  SEXP out = Rf_protect(cpp_df_select(x, j)); ++NP;
  // Subset rows
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


// The below function is deprecated and only kept to not break dependencies

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
  int out_size = REAL(VECTOR_ELT(clean_indices_info, 1))[0];
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
      if (is_simple_atomic_vec(df_var)){
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
      if (is_simple_atomic_vec(df_var)){
        SEXP list_var = Rf_protect(sset_vec(df_var, indices, check_indices));
        Rf_copyMostAttrib(df_var, list_var);
        int has_names = !Rf_isNull(Rf_getAttrib(df_var, R_NamesSymbol));
        if (has_names){
          SEXP old_names = Rf_protect(Rf_getAttrib(df_var, R_NamesSymbol));
          SEXP new_names = Rf_protect(sset_vec(
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
