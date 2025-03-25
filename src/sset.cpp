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

  SEXP x_seq = SHIELD(cpp_seq_len(xn));

  if (xn > integer_max_){
    SHIELD(exclude = coerce_vec(exclude, REALSXP));
    double *p_x = REAL(x_seq);
    double *p_excl = REAL(exclude);

    for (int j = 0; j < m; ++j) {
      if (p_excl[j] != p_excl[j]) continue;
      if (p_excl[j] > 0){
        R_Free(p_keep);
        YIELD(2);
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
    SEXP out = SHIELD(new_vec(REALSXP, out_size));
    double *p_out = REAL(out);
    while(k != out_size){
      if (p_keep[i]){
        p_out[k++] = p_x[i];
      }
      ++i;
    }
    R_Free(p_keep);
    YIELD(3);
    return out;
  } else {
    int *p_x = INTEGER(x_seq);
    int *p_excl = INTEGER(exclude);

    for (int j = 0; j < m; ++j) {
      if (p_excl[j] == NA_INTEGER) continue;
      if (p_excl[j] > 0){
        R_Free(p_keep);
        YIELD(1);
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
    SEXP out = SHIELD(new_vec(INTSXP, out_size));
    int *p_out = INTEGER(out);
    while(k != out_size){
      if (p_keep[i]){
        p_out[k++] = p_x[i];
      }
      ++i;
    }
    R_Free(p_keep);
    YIELD(2);
    return out;
  }
}

// Cleans indices for subsetting
// Also returns metadata regarding final vec size and if indices should be
// checked (internal flag)

[[cpp11::register]]
SEXP clean_indices(SEXP indices, SEXP x){
  R_xlen_t xn = vec_length(x);
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

  if (TYPEOF(indices) == STRSXP){
    SEXP names = SHIELD(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
    if (Rf_isNull(names)){
      YIELD(NP);
      Rf_error("Cannot subset on the names of an unnamed vector");
    }
    if (n < 10000){
      SHIELD(indices = Rf_match(names, indices, NA_INTEGER)); ++NP;
    } else {
      SHIELD(indices = cheapr_fast_match(indices, names)); ++NP;
    }
  }

  if (is_compact_seq(indices)){
    clean_indices = indices;

    SEXP seq_data = SHIELD(compact_seq_data(indices)); ++NP;
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    out_size = get_alt_final_sset_size(xn, from, to, by);
    check_indices = true;

  } else if (TYPEOF(indices) == LGLSXP){
    if (Rf_length(indices) != xn){
      Rf_error("`length(i)` must match `length(x)` when `i` is a logical vector");
    }
    clean_indices = SHIELD(cpp_which_(indices, false)); ++NP;
    n = Rf_xlength(clean_indices);
    out_size = n;
    check_indices = false;
  } else if (xn > integer_max_){
    SHIELD(indices = coerce_vec(indices, REALSXP)); ++NP;

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
      clean_indices = SHIELD(exclude_locs(indices, xn)); ++NP;
      check_indices = false;
      out_size = Rf_length(clean_indices);
    } else {
      clean_indices = indices;
      out_size = pos_count + na_count;
    }
  } else {

    SHIELD(indices = coerce_vec(indices, INTSXP)); ++NP;

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
      clean_indices = SHIELD(exclude_locs(indices, xn)); ++NP;
      check_indices = false;
      out_size = Rf_length(clean_indices);
    } else {
      clean_indices = indices;
      out_size = pos_count + na_count;
    }
  }

  SEXP out = SHIELD(new_vec(VECSXP, 3)); ++NP;

  // There are the `Rf_Scalar` shortcuts BUT R crashes sometimes when
  // using the scalar logical shortcuts so I avoid it

  SEXP out_size_sexp = SHIELD(Rf_ScalarReal(out_size)); ++NP;
  SEXP check_indices_sexp = SHIELD(scalar_lgl(check_indices)); ++NP;
  SET_VECTOR_ELT(out, 0, clean_indices);
  SET_VECTOR_ELT(out, 1, out_size_sexp);
  SET_VECTOR_ELT(out, 2, check_indices_sexp);

  YIELD(NP);
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
    out = SHIELD(new_vec(TYPEOF(x), out_size)); ++NP;
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
    out = SHIELD(new_vec(REALSXP, out_size)); ++NP;
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
    out = SHIELD(new_vec(STRSXP, out_size)); ++NP;
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
    out = SHIELD(new_vec(CPLXSXP, out_size)); ++NP;
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
    out = SHIELD(new_vec(RAWSXP, out_size)); ++NP;
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
    out = SHIELD(new_vec(VECSXP, out_size)); ++NP;
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
  YIELD(NP);
  return out;
}

// Vector subset
// OOB, zeros and NA are checked when `check = T`

SEXP sset_vec(SEXP x, SEXP indices, bool check){

  if (is_compact_seq(indices)){
    SEXP seq_data = SHIELD(compact_seq_data(indices));
    R_xlen_t from = REAL(seq_data)[0];
    R_xlen_t to = REAL(seq_data)[1];
    R_xlen_t by = REAL(seq_data)[2];
    SEXP out = SHIELD(cpp_sset_range(x, from, to, by));
    YIELD(2);
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
      out = SHIELD(new_vec(TYPEOF(x), n)); ++NP;
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
      out = SHIELD(new_vec(REALSXP, n)); ++NP;
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
      out = SHIELD(new_vec(STRSXP, n)); ++NP;

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
      out = SHIELD(new_vec(CPLXSXP, n)); ++NP;
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
      out = SHIELD(new_vec(RAWSXP, n)); ++NP;
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
      out = SHIELD(new_vec(VECSXP, n)); ++NP;
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
      YIELD(NP);
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
    out = SHIELD(new_vec(TYPEOF(x), n)); ++NP;
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
    out = SHIELD(new_vec(REALSXP, n)); ++NP;
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
    out = SHIELD(new_vec(STRSXP, n)); ++NP;

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
    out = SHIELD(new_vec(CPLXSXP, n)); ++NP;
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
    out = SHIELD(new_vec(RAWSXP, n)); ++NP;
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
    out = SHIELD(new_vec(VECSXP, n)); ++NP;
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
    YIELD(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  }
  if (!Rf_isNull(out) && out_size != n){
    SHIELD(out = Rf_xlengthgets(out, out_size)); ++NP;
  }
  YIELD(NP);
  return out;
}

// Basic and fast vector subset, exported to R

[[cpp11::register]]
SEXP cpp_sset(SEXP x, SEXP indices){
  int NP = 0;

  SEXP indices_metadata = SHIELD(clean_indices(indices, x)); ++NP;
  SEXP clean_indices = SHIELD(VECTOR_ELT(indices_metadata, 0)); ++NP;
  bool check_indices = LOGICAL(VECTOR_ELT(indices_metadata, 2))[0];
  SEXP out = SHIELD(sset_vec(x, clean_indices, check_indices)); ++NP;

  Rf_copyMostAttrib(x, out);

  // Subset names

  SEXP xnames = SHIELD(Rf_getAttrib(x, R_NamesSymbol)); ++NP;
  if (!Rf_isNull(xnames)){
    Rf_setAttrib(out, R_NamesSymbol, sset_vec(xnames, clean_indices, check_indices));
  }
  YIELD(NP);
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
  SHIELD(x = altrep_materialise(x)); ++NP;

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
    out = SHIELD(set ? x : Rf_duplicate(x)); ++NP;
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
    out = SHIELD(set ? x : Rf_duplicate(x)); ++NP;
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
    out = SHIELD(set ? x : Rf_duplicate(x)); ++NP;
    const SEXP *p_out = STRING_PTR_RO(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      SEXP left = SHIELD(p_out[i]);
      SET_STRING_ELT(out, i, p_out[k]);
      SET_STRING_ELT(out, k, left);
      YIELD(1);
    }
    break;
  }
  case CPLXSXP: {
    out = SHIELD(set ? x : Rf_duplicate(x)); ++NP;
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
    out = SHIELD(set ? x : Rf_duplicate(x)); ++NP;
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
    YIELD(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
    // We can reverse data frames in-place with the below commented-out code
    // It will not work properly though for data.tables

    // case VECSXP: {
    //   if (recursive){
    //   rev_names = false;
    //   out = SHIELD(new_vec(VECSXP, n)); ++NP;
    //   SHALLOW_DUPLICATE_ATTRIB(out, x);
    //   for (R_xlen_t i = 0; i < n; ++i){
    //     SET_VECTOR_ELT(out, i, cpp_rev(VECTOR_ELT(x, i), true, set));
    //   }
    //   break;
    // } else if (!recursive && (!Rf_isObject(x) || is_df(x))){
    //   out = SHIELD(set ? x : list_shallow_copy(x, false)); ++NP;
    //   if (!set){
    //     // SHALLOW_DUPLICATE_ATTRIB(out, x);
    //     cpp_copy_names(x, out, true);
    //   }
    //   const SEXP *p_out = VECTOR_PTR_RO(out);
    //   for (R_xlen_t i = 0; i < half; ++i) {
    //     k = n2 - i;
    //     SEXP left = SHIELD(p_out[i]);
    //     SET_VECTOR_ELT(out, i, p_out[k]);
    //     SET_VECTOR_ELT(out, k, left);
    //     YIELD(1);
    //   }
    //     break;
    //   }
    // }
    // default: {
    //   // set rev_names to false because catch-all r_rev() will handle that
    //   rev_names = false;
    //   cpp11::function r_rev = cpp11::package("base")["rev"];
    //   if (set){
    //     YIELD(NP);
    //     Rf_error("Can't reverse in place here");
    //   }
    //   out = SHIELD(r_rev(x)); ++NP;
    //   break;
    // }
  }
  // // If x has names, reverse them too
  if (rev_names && !Rf_isNull(Rf_getAttrib(out, R_NamesSymbol))){
    SEXP old_names = SHIELD(Rf_getAttrib(out, R_NamesSymbol)); ++NP;
    // should be okay to rev in-place here because we already copied the names
    Rf_setAttrib(out, R_NamesSymbol, cpp_rev(old_names, true));
  }
  YIELD(NP);
  return out;
}

// SEXP cpp_cheapr_rev(SEXP x){
//   SEXP indices = SHIELD(base_colon(vec_length(x), 0));
//   SEXP out;
//   if (is_simple_atomic_vec(x)){
//     out = SHIELD(cpp_sset(x, indices));
//   } else {
//     out = SHIELD(cheapr_sset(x, indices));
//   }
//   YIELD(2);
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
    n_rows = df_nrow(x),
    n_locs = Rf_length(locs);

  // Flag to check indices
  bool check = true;

  SEXP names = SHIELD(Rf_getAttrib(x, R_NamesSymbol)); ++NP;

  SEXP cols;

  if (Rf_isNull(locs)){
    // If NULL then select all cols
    cols = SHIELD(cpp_seq_len(n_cols)); ++NP;
    n_locs = n_cols;
    check = false;
  } else if (Rf_isString(locs)){
    cols = SHIELD(Rf_match(names, locs, NA_INTEGER)); ++NP;
  } else if (Rf_isLogical(locs)){
    // If logical then find locs using `which_()`
    if (Rf_length(locs) != n_cols){
      YIELD(NP);
      Rf_error("`length(j)` must match `ncol(x)` when `j` is a logical vector");
    }
    cols = SHIELD(cpp_which_(locs, false)); ++NP;
    n_locs = Rf_length(cols);
    check = false;
  } else {
    // Catch-all make sure cols is an int vector
    cols = SHIELD(coerce_vec(locs, INTSXP)); ++NP;
  }

  // Negative subscripting
  if (n_locs > 0 && INTEGER(cols)[0] != NA_INTEGER && INTEGER(cols)[0] < 0){
    SHIELD(cols = exclude_locs(cols, n_cols)); ++NP;
    n_locs = Rf_length(cols);
    check = false;
  }

  SEXP out = SHIELD(new_vec(VECSXP, n_locs)); ++NP;
  SEXP out_names = SHIELD(new_vec(STRSXP, n_locs)); ++NP;

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
    SHIELD(out = Rf_lengthgets(out, k)); ++NP;
    SHIELD(out_names = Rf_lengthgets(out_names, k)); ++NP;
  }

  // Make a plain data frame
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(n_rows));
  Rf_classgets(out, Rf_mkString("data.frame"));
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_df_slice(SEXP x, SEXP indices, bool check){

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame`, not a %s", Rf_type2char(TYPEOF(x)));
  }

  if (Rf_isNull(indices)){
    return x;
  }

  int xn = df_nrow(x);
  int ncols = Rf_length(x);
  int NP = 0;
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = SHIELD(new_vec(VECSXP, ncols)); ++NP;


  // Clean indices and get metadata

  int out_size;
  bool check_indices;

  if (check){
    const SEXP clean_indices_info = SHIELD(clean_indices(indices, x)); ++NP;
    SHIELD(indices = VECTOR_ELT(clean_indices_info, 0)); ++NP;
    out_size = REAL(VECTOR_ELT(clean_indices_info, 1))[0];
    check_indices = LOGICAL(VECTOR_ELT(clean_indices_info, 2))[0];
  } else {
    out_size = Rf_length(indices);
    check_indices = false;
  }

  // Subset columns

  PROTECT_INDEX index1, index2, index3;
  SEXP df_var, names, list_var;
  R_ProtectWithIndex(df_var = R_NilValue, &index1); ++NP;
  R_ProtectWithIndex(names = R_NilValue, &index2); ++NP;
  R_ProtectWithIndex(list_var = R_NilValue, &index3); ++NP;

  for (int j = 0; j < ncols; ++j){
    R_Reprotect(df_var = p_x[j], index1);
    R_Reprotect(names = Rf_getAttrib(df_var, R_NamesSymbol), index2);
    if (is_simple_vec(df_var)){
      R_Reprotect(list_var = sset_vec(df_var, indices, check_indices), index3);
      Rf_copyMostAttrib(df_var, list_var);
      Rf_setAttrib(list_var, R_NamesSymbol, sset_vec(names, indices, check_indices));
    } else {
      R_Reprotect(list_var = cheapr_sset(df_var, indices), index3);
    }
    SET_VECTOR_ELT(out, j, list_var);
  }

  cpp_copy_names(x, out, false);

  // list to data frame object
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(out_size));
  Rf_classgets(out, Rf_mkString("data.frame"));
  YIELD(NP);
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
  SEXP out = SHIELD(cpp_df_select(x, j)); ++NP;
  // Subset rows
  SHIELD(out = cpp_df_slice(out, i, true)); ++NP;

  if (keep_attrs){
    SEXP names = SHIELD(Rf_getAttrib(out, R_NamesSymbol)); ++NP;
    SEXP row_names = SHIELD(Rf_getAttrib(out, R_RowNamesSymbol)); ++NP;
    Rf_copyMostAttrib(x, out);
    Rf_setAttrib(out, R_NamesSymbol, names);
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  }
  YIELD(NP);
  return out;
}
