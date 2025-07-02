#include <vector>
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

  int32_t NP = 0;
  R_xlen_t n = xn;
  R_xlen_t m = Rf_length(exclude);
  R_xlen_t out_size, idx;
  R_xlen_t exclude_count = 0;
  R_xlen_t i = 0, k = 0;

  // Which elements do we keep?
  std::vector<uint8_t> keep(n);
  std::fill(keep.begin(), keep.end(), 1);

  if (xn > INTEGER_MAX){
    SHIELD(exclude = coerce_vec(exclude, REALSXP)); ++NP;
    double *p_excl = REAL(exclude);

    for (int j = 0; j < m; ++j) {
      if (is_na_dbl(p_excl[j])) continue;
      if (p_excl[j] > 0){
        YIELD(NP);
        Rf_error("Cannot mix positive and negative subscripts");
      }
      idx = -p_excl[j];
      // Check keep array for already assigned FALSE to avoid double counting
      if (idx > 0 && idx <= n && keep[idx - 1] == 1){
        keep[idx - 1] = 0;
        ++exclude_count;
      }
    }
    out_size = n - exclude_count;
    SEXP out = SHIELD(new_vec(REALSXP, out_size)); ++NP;
    double* RESTRICT p_out = REAL(out);
    while(k != out_size){
      if (keep[i] == 1){
        p_out[k++] = i + 1;
      }
      ++i;
    }
    YIELD(NP);
    return out;
  } else {
    int *p_excl = INTEGER(exclude);

    for (int j = 0; j < m; ++j) {
      if (is_na_int(p_excl[j])) continue;
      if (p_excl[j] > 0){
        YIELD(NP);
        Rf_error("Cannot mix positive and negative subscripts");
      }
      idx = -p_excl[j];
      // Check keep array for already assigned FALSE to avoid double counting
      if (idx > 0 && idx <= n && keep[idx - 1] == 1){
        keep[idx - 1] = 0;
        ++exclude_count;
      }
    }
    out_size = n - exclude_count;
    SEXP out = SHIELD(new_vec(INTSXP, out_size)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
    while(k != out_size){
      if (keep[i++] == 1){
        p_out[k++] = i;
      }
    }
    YIELD(NP);
    return out;
  }
}

// Cleans indices for subsetting
// Also returns metadata regarding final vec size and if indices should be
// checked (internal flag)

SEXP clean_indices(SEXP indices, SEXP x, bool count){
  R_xlen_t xn = vec_length(x);
  int32_t NP = 0;
  R_xlen_t zero_count = 0,
    pos_count = 0,
    oob_count = 0,
    na_count = 0,
    neg_count = 0;

  R_xlen_t n = Rf_xlength(indices);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  int_fast64_t out_size = NA_INTEGER64;
  bool check_indices = true;
  SEXP clean_indices = R_NilValue;

  if (TYPEOF(indices) == STRSXP){
    SEXP names = SHIELD(get_names(x)); ++NP;
    if (is_null(names)){
      YIELD(NP);
      Rf_error("Cannot subset on the names of an unnamed vector");
    }
    if (is_df(x)){
      YIELD(NP);
      Rf_error("Cannot subset rows of a data frame using a character vector");
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
      YIELD(NP);
      Rf_error("`length(i)` must match `length(x)` when `i` is a logical vector");
    }
    clean_indices = SHIELD(cpp_which_(indices, false)); ++NP;
    n = Rf_xlength(clean_indices);
    out_size = n;
    check_indices = false;
  } else if (xn > INTEGER_MAX){

    SHIELD(clean_indices = coerce_vec(indices, REALSXP)); ++NP;

    if (count){

      const double *pi = REAL_RO(clean_indices);

      // Counting the number of:
      // Zeroes
      // Out-of-bounds indices
      // Positive indices
      // NA indices
      // From this we can also work out the number of negatives

      if (n_cores > 1){
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
        YIELD(NP);
        Rf_error("Cannot mix positive and negative indices");
      }

      // Should a simplified sset method be used?

      check_indices = !(oob_count == 0 && na_count == 0 && zero_count == 0);

      if (neg_count > 0){
        clean_indices = SHIELD(exclude_locs(clean_indices, xn)); ++NP;
        check_indices = false;
        out_size = Rf_length(clean_indices);
      } else {
        out_size = pos_count + na_count;
      }
    }
  } else {

    SHIELD(clean_indices = coerce_vec(indices, INTSXP)); ++NP;


    if (count){

      const int *pi = INTEGER_RO(clean_indices);

      // Counting the number of:
      // Zeroes
      // Out-of-bounds indices
      // Positive indices
      // NA indices
      // From this we can also work out the number of negatives

      if (n_cores > 1){
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
        YIELD(NP);
        Rf_error("Cannot mix positive and negative indices");
      }

      // Should a simplified sset method be used?

      check_indices = !(oob_count == 0 && na_count == 0 && zero_count == 0);

      if (neg_count > 0){
        clean_indices = SHIELD(exclude_locs(clean_indices, xn)); ++NP;
        check_indices = false;
        out_size = Rf_length(clean_indices);
      } else {
        out_size = pos_count + na_count;
      }
    }
  }

  SEXP out = SHIELD(new_vec(VECSXP, 3)); ++NP;

  // There are the `Rf_Scalar` shortcuts BUT R crashes sometimes when
  // using the scalar logical shortcuts so I avoid it
  SET_VECTOR_ELT(out, 0, clean_indices);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(is_na_int64(out_size) ? NA_REAL : static_cast<double>(out_size)));
  SET_VECTOR_ELT(out, 2, scalar_lgl(check_indices));

  YIELD(NP);
  return out;
}


// A range-based subset method
// Can be readily used when indices are an altrep compact integer sequence
// Also works with negative-indexing

SEXP cpp_sset_range(SEXP x, R_xlen_t from, R_xlen_t to, R_xlen_t by){
  int32_t NP = 0;
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
    const int *p_x = INTEGER_RO(x);
    out = SHIELD(new_vec(TYPEOF(x), out_size)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
    if (double_loop){
      safe_memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * sizeof(int));
      safe_memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * sizeof(int));
    } else {
      if (by > 0){
        safe_memmove(&p_out[0], &p_x[istart - 1], in_bounds_size * sizeof(int));
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
    const double *p_x = REAL_RO(x);
    out = SHIELD(new_vec(REALSXP, out_size)); ++NP;
    double* RESTRICT p_out = REAL(out);
    if (double_loop){
      safe_memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * sizeof(double));
      safe_memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * sizeof(double));
    } else {
      if (by > 0){
        if (in_bounds_size != 0){
          safe_memmove(&p_out[0], &p_x[istart - 1], in_bounds_size * sizeof(double));
        }
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
    Rcomplex* RESTRICT p_out = COMPLEX(out);
    if (double_loop){
      safe_memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * sizeof(Rcomplex));
      safe_memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * sizeof(Rcomplex));
    } else {
      if (by > 0){
        if (in_bounds_size != 0){
          safe_memmove(&p_out[0], &p_x[istart - 1], in_bounds_size * sizeof(Rcomplex));
        }
        for (R_xlen_t i = 0; i < n_oob; ++i){
          R_xlen_t tempi = in_bounds_size + i;
          p_out[tempi].r = NA_REAL;
          p_out[tempi].i = NA_REAL;
        }
      } else {
        for (R_xlen_t i = 0; i < n_oob; ++i){
          p_out[i].r = NA_REAL;
          p_out[i].i = NA_REAL;
        }
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
// OOB, zeros, NA and negative values are checked when `check = T`

SEXP sset_vec(SEXP x, SEXP indices, bool check){

  SEXP out = R_NilValue;
  int xtype = TYPEOF(x);

  if (check){

    if (Rf_xlength(x) > INTEGER_MAX){

      int_fast64_t xn = Rf_xlength(x);

      int_fast64_t
      n = Rf_xlength(indices), k = 0, j;

      const double* pind = REAL_RO(indices);

      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(R_NilValue);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int* p_x = INTEGER_RO(x);
        out = SHIELD(new_vec(xtype, n));
        int* RESTRICT p_out = INTEGER(out);

          for (int_fast64_t i = 0; i < n; ++i){
            j = pind[i];
            if (j < 0){
              SEXP new_indices = SHIELD(exclude_locs(indices, xn));
              SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
              YIELD(3);
              return out2;
            } else if (j != 0){
              p_out[k++] = (pind[i] != pind[i] || j > xn) ? NA_INTEGER : p_x[j - 1];
            }
          }
        break;
      }
      case REALSXP: {
        const double* p_x = REAL_RO(x);
        out = SHIELD(new_vec(REALSXP, n));
        double* RESTRICT p_out = REAL(out);
          for (int_fast64_t i = 0; i < n; ++i){
            j = pind[i];
            if (j < 0){
              SEXP new_indices = SHIELD(exclude_locs(indices, xn));
              SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
              YIELD(3);
              return out2;
            } else if (j != 0){
              p_out[k++] = (pind[i] != pind[i] || j > xn) ? NA_REAL : p_x[j - 1];
            }
          }
        break;
      }
      case STRSXP: {
        const SEXP *p_x = STRING_PTR_RO(x);
        out = SHIELD(new_vec(STRSXP, n));

          for (int_fast64_t i = 0; i < n; ++i){
            j = pind[i];
            if (j < 0){
              SEXP new_indices = SHIELD(exclude_locs(indices, xn));
              SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
              YIELD(3);
              return out2;
            } else if (j != 0){
              SET_STRING_ELT(out, k++, (pind[i] != pind[i] || j > xn) ? NA_STRING : p_x[j - 1]);
            }
          }
        break;
      }
      case CPLXSXP: {
        const Rcomplex* p_x = COMPLEX_RO(x);
        out = SHIELD(new_vec(CPLXSXP, n));
        Rcomplex* RESTRICT p_out = COMPLEX(out);
          for (int_fast64_t i = 0; i < n; ++i){

            j = pind[i];

            if (j < 0){
              SEXP new_indices = SHIELD(exclude_locs(indices, xn));
              SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
              YIELD(3);
              return out2;
            } else if (j != 0){
              if (pind[i] != pind[i] || j > xn){
                p_out[k].r = NA_REAL;
                p_out[k].i = NA_REAL;
              } else {
                p_out[k].r = p_x[j - 1].r;
                p_out[k].i = p_x[j - 1].i;
              }
              ++k;
            }

          }
        break;
      }
      case RAWSXP: {
        const Rbyte *p_x = RAW_RO(x);
        out = SHIELD(new_vec(RAWSXP, n));
          for (int_fast64_t i = 0; i < n; ++i){
            j = pind[i];
            if (j < 0){
              SEXP new_indices = SHIELD(exclude_locs(indices, xn));
              SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
              YIELD(3);
              return out2;
            } else if (j != 0){
              SET_RAW_ELT(out, k++, (pind[i] != pind[i] || j > xn) ? 0 : p_x[j - 1]);
            }
          }
        break;
      }
      case VECSXP: {
        const SEXP *p_x = VECTOR_PTR_RO(x);
        out = SHIELD(new_vec(VECSXP, n));
          for (int_fast64_t i = 0; i < n; ++i){
            j = pind[i];
            if (j < 0){
              SEXP new_indices = SHIELD(exclude_locs(indices, xn));
              SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
              YIELD(3);
              return out2;
            } else if (j != 0){
              SET_VECTOR_ELT(out, k++, (pind[i] != pind[i] || j > xn) ? R_NilValue : p_x[j - 1]);
            }
          }
        break;
      }
      default: {
        YIELD(1);
        Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(xtype));
      }
      }
      if (!is_null(out) && k != n){
        SHIELD(out = Rf_xlengthgets(out, k));
        YIELD(2);
        return out;
      } else {
        YIELD(1);
        return out;
      }
    } else {

      unsigned int
      xn = Rf_length(x),
        n = Rf_xlength(indices),
        k = 0,
        na_val = NA_INTEGER,
        j;

      const int *pind = INTEGER_RO(indices);
      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(R_NilValue);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int* p_x = INTEGER_RO(x);
        out = SHIELD(new_vec(xtype, n));
        int* RESTRICT p_out = INTEGER(out);
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            p_out[k++] = p_x[--j];
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            p_out[k++] = NA_INTEGER;
          }
        }
        break;
      }
      case REALSXP: {
        const double* p_x = REAL_RO(x);
        out = SHIELD(new_vec(REALSXP, n));
        double* RESTRICT p_out = REAL(out);
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            p_out[k++] = p_x[--j];
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            p_out[k++] = NA_REAL;
          }
        }
        break;
      }
      case STRSXP: {
        const SEXP *p_x = STRING_PTR_RO(x);
        out = SHIELD(new_vec(STRSXP, n));
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            SET_STRING_ELT(out, k++, p_x[--j]);
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            SET_STRING_ELT(out, k++, NA_STRING);
          }
        }
        break;
      }
      case CPLXSXP: {
        const Rcomplex* p_x = COMPLEX_RO(x);
        out = SHIELD(new_vec(CPLXSXP, n));
        Rcomplex* RESTRICT p_out = COMPLEX(out);
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            p_out[k++] = p_x[--j];
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            p_out[k].r = NA_REAL;
            p_out[k++].i = NA_REAL;
          }
        }
        break;
      }
      case RAWSXP: {
        const Rbyte *p_x = RAW_RO(x);
        out = SHIELD(new_vec(RAWSXP, n));
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
           SET_RAW_ELT(out, k++, p_x[--j]);
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            SET_RAW_ELT(out, k++, 0);
          }
        }
        break;
      }
      case VECSXP: {
        const SEXP *p_x = VECTOR_PTR_RO(x);
        out = SHIELD(new_vec(VECSXP, n));
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            SET_VECTOR_ELT(out, k++, p_x[--j]);
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            SET_VECTOR_ELT(out, k++, R_NilValue);
          }
        }
        break;
      }
      default: {
        YIELD(1);
        Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(xtype));
      }
      }
      if (!is_null(out) && k != n){
        SHIELD(out = Rf_lengthgets(out, k));
        YIELD(2);
        return out;
      } else {
        YIELD(1);
        return out;
      }
    }
  } else {

    if (Rf_xlength(x) > INTEGER_MAX){

      int_fast64_t n = Rf_xlength(indices);

      const double *pind = REAL_RO(indices);

      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(R_NilValue);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int *p_x = INTEGER_RO(x);
        out = SHIELD(new_vec(xtype, n));
        int* RESTRICT p_out = INTEGER(out);
        OMP_FOR_SIMD
        for (int_fast64_t i = 0; i < n; ++i){
          p_out[i] = p_x[static_cast<int_fast64_t>(pind[i] - 1.0)];
        }
        break;
      }
      case REALSXP: {
        const double *p_x = REAL_RO(x);
        out = SHIELD(new_vec(REALSXP, n));
        double* RESTRICT p_out = REAL(out);
        OMP_FOR_SIMD
        for (int_fast64_t i = 0; i < n; ++i){
          p_out[i] = p_x[static_cast<int_fast64_t>(pind[i] - 1.0)];
        }
        break;
      }
      case STRSXP: {
        const SEXP *p_x = STRING_PTR_RO(x);
        out = SHIELD(new_vec(STRSXP, n));
        for (int_fast64_t i = 0; i < n; ++i) SET_STRING_ELT(out, i, p_x[static_cast<int_fast64_t>(pind[i] - 1.0)]);
        break;
      }
      case CPLXSXP: {
        const Rcomplex *p_x = COMPLEX_RO(x);
        out = SHIELD(new_vec(CPLXSXP, n));
        Rcomplex* RESTRICT p_out = COMPLEX(out);
        for (int_fast64_t i = 0; i < n; ++i){
          p_out[i].r = p_x[static_cast<int_fast64_t>(pind[i] - 1.0)].r;
          p_out[i].i = p_x[static_cast<int_fast64_t>(pind[i] - 1.0)].i;
        }
        break;
      }
      case RAWSXP: {
        const Rbyte *p_x = RAW_RO(x);
        out = SHIELD(new_vec(RAWSXP, n));
        for (int_fast64_t i = 0; i < n; ++i) SET_RAW_ELT(out, i, p_x[static_cast<int_fast64_t>(pind[i] - 1.0)]);
        break;
      }
      case VECSXP: {
        const SEXP *p_x = VECTOR_PTR_RO(x);
        out = SHIELD(new_vec(VECSXP, n));
        for (int_fast64_t i = 0; i < n; ++i) SET_VECTOR_ELT(out, i, p_x[static_cast<int_fast64_t>(pind[i] - 1.0)]);
        break;
      }
      default: {
        YIELD(1);
        Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(xtype));
      }
      }
      YIELD(1);
      return out;
    } else {

      int n = Rf_length(indices);

      const int *pind = INTEGER_RO(indices);
      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(R_NilValue);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int *p_x = INTEGER_RO(x);
        out = SHIELD(new_vec(xtype, n));
        int* RESTRICT p_out = INTEGER(out);
        OMP_FOR_SIMD
        for (int i = 0; i != n; ++i){
          p_out[i] = p_x[pind[i] - 1];
        }
        break;
      }
      case REALSXP: {
        const double *p_x = REAL_RO(x);
        out = SHIELD(new_vec(REALSXP, n));
        double* RESTRICT p_out = REAL(out);
        OMP_FOR_SIMD
        for (int i = 0; i != n; ++i){
          p_out[i] = p_x[pind[i] - 1];
        }
        break;
      }
      case STRSXP: {
        const SEXP *p_x = STRING_PTR_RO(x);
        out = SHIELD(new_vec(STRSXP, n));
        for (int i = 0; i != n; ++i) SET_STRING_ELT(out, i, p_x[pind[i] - 1]);
        break;
      }
      case CPLXSXP: {
        const Rcomplex *p_x = COMPLEX_RO(x);
        out = SHIELD(new_vec(CPLXSXP, n));
        Rcomplex* RESTRICT p_out = COMPLEX(out);
        for (int i = 0; i != n; ++i){
          p_out[i].r = p_x[pind[i] - 1].r;
          p_out[i].i = p_x[pind[i] - 1].i;
        }
        break;
      }
      case RAWSXP: {
        const Rbyte *p_x = RAW_RO(x);
        out = SHIELD(new_vec(RAWSXP, n));
        for (int i = 0; i != n; ++i) SET_RAW_ELT(out, i, p_x[pind[i] - 1]);
        break;
      }
      case VECSXP: {
        const SEXP *p_x = VECTOR_PTR_RO(x);
        out = SHIELD(new_vec(VECSXP, n));
        for (int i = 0; i != n; ++i) SET_VECTOR_ELT(out, i, p_x[pind[i] - 1]);
        break;
      }
      default: {
        YIELD(1);
        Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(xtype));
      }
      }
      YIELD(1);
      return out;
    }
  }
}

// Safe and very efficient reverse in-place (or with copy)

// matrix structure is preserved with cpp_rev

// For `set = F` the same can be accomplished (slightly faster) with
// the range-based subset
// Keeping this implementation as it works seamlessly for in-place and
// normal reverse

SEXP rev(SEXP x, bool set){
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t half = n / 2;
  R_xlen_t n2 = n - 1; // Offset n for 0-indexing
  R_xlen_t k;
  int32_t NP = 0;
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
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
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
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
    double* RESTRICT p_out = REAL(out);
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
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
    const SEXP *p_out = STRING_PTR_RO(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      SEXP left = p_out[i];
      SET_STRING_ELT(out, i, p_out[k]);
      SET_STRING_ELT(out, k, left);
    }
    break;
  }
  case CPLXSXP: {
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
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
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
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
  }
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_rev(SEXP x, bool set){
  SEXP out = SHIELD(rev(x, set));
  SEXP names = SHIELD(get_names(x));
  SHIELD(names = rev(names, set && !ALTREP(x)));
  set_names(out, names);
  YIELD(3);
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

  int32_t NP = 0,
    n_cols = Rf_length(x),
    n_rows = df_nrow(x),
    n_locs = Rf_length(locs);

  // Flag to check indices
  bool check = true;

  SEXP names = SHIELD(get_names(x)); ++NP;

  SEXP cols;
  int loc_type = TYPEOF(locs);

  if (loc_type == NILSXP){
    // If NULL then select all cols
    cols = SHIELD(cpp_seq_len(n_cols)); ++NP;
    n_locs = n_cols;
    check = false;
  } else if (loc_type == STRSXP){
    cols = SHIELD(Rf_match(names, locs, NA_INTEGER)); ++NP;
    int *match_locs = INTEGER(cols);
    if (cpp_any_na(cols, false)){
      for (int i = 0; i < Rf_length(cols); ++i){
        if (is_na_int(match_locs[i])){
          const char *bad_loc = utf8_char(STRING_ELT(locs, i));
          YIELD(NP);
          Rf_error("Column %s does not exist", bad_loc);
        }
      }
    }
    check = false;
  } else if (loc_type == LGLSXP){
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

  const int *p_cols = INTEGER_RO(cols);

  SEXP out = SHIELD(new_vec(VECSXP, n_locs)); ++NP;
  SEXP out_names = SHIELD(new_vec(STRSXP, n_locs)); ++NP;

  const SEXP *p_x = VECTOR_PTR_RO(x);
  const SEXP *p_names = STRING_PTR_RO(names);
  int k = 0;
  int col;

  if (check){
    for (int i = 0; i < n_locs; ++i) {
      col = p_cols[i];
      if (col == NA_INTEGER){
        YIELD(NP);
        Rf_error("Cannot select `NA` column locations in %s", __func__);
      } else if (col >= 1 && col <= n_cols){
        --col;
        SET_VECTOR_ELT(out, k, p_x[col]);
        SET_STRING_ELT(out_names, k, p_names[col]);
        ++k;
      } else if (col < 0){
        // This can only happen when there is a mix of pos & neg
        // but wasn't captured by the negative subscripting section
        // because that only looks to see if the 1st element is neg
        YIELD(NP);
        Rf_error("Cannot mix positive and negative subscripts in %s", __func__);
      } else if (col != 0){
        YIELD(NP);
        Rf_error("There is no column location %d: ", col);
      }
    }
  } else {
    for (int i = 0; i < n_locs; ++i){
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
  Rf_classgets(out, make_utf8_str("data.frame"));
  set_names(out, out_names);
  YIELD(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_df_slice(SEXP x, SEXP indices, bool check){

  if (is_null(indices)){
    return x;
  }
  int ncols = Rf_length(x);
  int32_t NP = 0;
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = SHIELD(new_vec(VECSXP, ncols)); ++NP;

  // Clean indices and get metadata

  int out_size;
  bool check_indices;

  if (check){
    SEXP clean_indices_info = SHIELD(clean_indices(indices, x, true)); ++NP;
    SHIELD(indices = VECTOR_ELT(clean_indices_info, 0)); ++NP;
    out_size = REAL(VECTOR_ELT(clean_indices_info, 1))[0];
    check_indices = LOGICAL(VECTOR_ELT(clean_indices_info, 2))[0];
  } else {
    out_size = Rf_length(indices);
    check_indices = false;
  }

  // Subset columns

  for (int j = 0; j < ncols; ++j){
    SET_VECTOR_ELT(out, j, cpp_sset(p_x[j], indices, check_indices));
  }

  SEXP names = SHIELD(get_names(x)); ++NP;
  set_names(out, names);

  // list to data frame object
  Rf_setAttrib(out, R_RowNamesSymbol, create_df_row_names(out_size));
  Rf_classgets(out, make_utf8_str("data.frame"));
  YIELD(NP);
  return out;
}


// Subset that does both selecting and slicing

[[cpp11::register]]
SEXP cpp_df_subset(SEXP x, SEXP i, SEXP j, bool check){

  if (!is_df(x)){
    Rf_error("`x` must be a `data.frame`, not a %s", Rf_type2char(TYPEOF(x)));
  }

  int32_t NP = 0;

  // Subset columns
  // `cpp_df_select()` always creates a shallow copy
  SEXP out = SHIELD(cpp_df_select(x, j)); ++NP;
  // Subset rows
  SHIELD(out = cpp_df_slice(out, i, check)); ++NP;
  SHIELD(out = rebuild(out, x, false)); ++NP;
  YIELD(NP);
  return out;
}


// Fast vector/data frame subset, exported to R

[[cpp11::register]]
SEXP cpp_sset(SEXP x, SEXP indices, bool check){

  if (is_simple_vec(x)){

    int32_t NP = 0;

    bool check2 = check;

    SEXP out = R_NilValue;
    SEXP names = R_NilValue;

    if (is_compact_seq(indices)){
      SEXP seq_data = SHIELD(compact_seq_data(indices)); ++NP;
      const double *p_data = REAL_RO(seq_data);
      out = SHIELD(cpp_sset_range(x, p_data[0], p_data[1], p_data[2])); ++NP;

      // Subset names
      names = SHIELD(get_names(x)); ++NP;
      SHIELD(names = cpp_sset_range(names, p_data[0], p_data[1], p_data[2])); ++NP;
    } else {
      if (check){
        SEXP indices_metadata = SHIELD(clean_indices(indices, x, false)); ++NP;
        SHIELD(indices = VECTOR_ELT(indices_metadata, 0)); ++NP;
        check2 = LOGICAL(VECTOR_ELT(indices_metadata, 2))[0];
      }
      out = SHIELD(sset_vec(x, indices, check2)); ++NP;

      // Subset names
      names = SHIELD(get_names(x)); ++NP;
      SHIELD(names = sset_vec(names, indices, check2)); ++NP;
    }
    Rf_copyMostAttrib(x, out);
    set_names(out, names);
    YIELD(NP);
    return out;
  } else if (is_df(x)){
    return cpp_df_subset(x, indices, R_NilValue, check);
  } else {
    // Normally we would use base_sset here BUT
    // we want to dispatch on some of sset's methods
    // This can all be re-worked and simplified by passing `...` to cpp_sset
    // from sset.default

    ////// IMPORTANT //////
    // The reason this doesn't result in infinite recursion is because
    // both this function and sset.default check that `x` is a simple vec

    // This can be more safely avoided by
    // writing an internal dispatch method in R
    // e.g. sset can be a non-generic function that calls an
    // internal generic function, e.g. `cheapr_sset`
    // I don't like this approach as much because the user can't see
    // all the available arguments like `j` as its hidden by the dots `...`
    // So as previously stated, fastest is to handle the args in C but
    // I don't currently know how to! :)

    return cheapr_sset(x, indices);
  }
}

// scalar subset
SEXP slice_loc(SEXP x, R_xlen_t i){

  if (i < 0){
    Rf_error("`i` must be >= 0");
  }

  if (Rf_isObject(x)){
    SEXP loc;
    if (i <= INTEGER_MAX){
      loc = SHIELD(Rf_ScalarInteger(i));
    } else {
      loc = SHIELD(Rf_ScalarReal(i));
    }
    SEXP out = SHIELD(cpp_sset(x, loc, true));
    YIELD(2);
    return out;
  }

  if (i == 0){
    return new_vec(TYPEOF(x), 0);
  }

  if (i > Rf_xlength(x)){
    return cpp_na_init(x, 1);
  }

  --i;

  switch ( TYPEOF(x) ){

  case NILSXP: {
    return R_NilValue;
  }
  case LGLSXP: {
    return scalar_lgl(LOGICAL(x)[i]);
  }
  case INTSXP: {
    return Rf_ScalarInteger(INTEGER(x)[i]);
  }
  case REALSXP: {
    return Rf_ScalarReal(REAL(x)[i]);
  }
  case STRSXP: {
    return Rf_ScalarString(STRING_ELT(x, i));
  }
  case CPLXSXP: {
    return Rf_ScalarComplex(COMPLEX(x)[i]);
  }
  case RAWSXP: {
    return Rf_ScalarRaw(RAW(x)[i]);
  }
  case VECSXP: {
    return VECTOR_ELT(x, i);
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}

// template<int SEXPTYPE, typename T>
// SEXP foo(T px, T pi, int n) {
//   if constexpr (SEXPTYPE == INTSXP) {
//     SEXP out = SHIELD(new_vec(INTSXP, n));
//     int* RESTRICT p_out = INTEGER(out);
//     OMP_FOR_SIMD
//     for (int i = 0; i < n; ++i){
//       p_out[i] = px[pi[i] - 1];
//     }
//     YIELD(1);
//     return out;
//   } else {
//     Rf_error("error");
//   }
// }
