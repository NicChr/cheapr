#include "cheapr.h"

// Subsetting vectors and data frames
// Includes a unique optimisation on range subsetting

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
    const SEXP *p_x = STRING_PTR_RO(x);
    SEXP out = Rf_protect(Rf_allocVector(STRSXP, n));
    for (int i = 0; i < n; ++i){
      SET_STRING_ELT(out, i, p_x[pind[i] - 1]);
    }
    Rf_unprotect(1);
    return out;
  }
  case CPLXSXP: {
    Rcomplex *p_x = COMPLEX(x);
    SEXP out = Rf_protect(Rf_allocVector(CPLXSXP, n));
    for (int i = 0; i < n; ++i){
      SET_COMPLEX_ELT(out, i, p_x[pind[i] - 1]);
    }
    Rf_unprotect(1);
    return out;
  }
  case RAWSXP: {
    Rbyte *p_x = RAW(x);
    SEXP out = Rf_protect(Rf_allocVector(RAWSXP, n));
    for (int i = 0; i < n; ++i){
      SET_RAW_ELT(out, i, p_x[pind[i] - 1]);
    }
    Rf_unprotect(1);
    return out;
  }
  case VECSXP: {
    const SEXP *p_x = VECTOR_PTR_RO(x);
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, n));
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
    out = Rf_protect(Rf_allocVector(TYPEOF(x), out_size));
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
    out = Rf_protect(Rf_allocVector(REALSXP, out_size));
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
    out = Rf_protect(Rf_allocVector(STRSXP, out_size));
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
    out = Rf_protect(Rf_allocVector(CPLXSXP, out_size));
    Rcomplex *p_out = COMPLEX(out);
    if (double_loop){
      memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * sizeof(Rcomplex));
      memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * sizeof(Rcomplex));
      // memmove(&p_out[0], &p_x[istart1 - 1], (iend1 - istart1 + 1) * 2 * sizeof(double));
      // memmove(&p_out[iend1 - istart1 + 1], &p_x[istart2 - 1], (iend2 - istart2 + 1) * 2 * sizeof(double));
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
    out = Rf_protect(Rf_allocVector(RAWSXP, out_size));
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
    out = Rf_protect(Rf_allocVector(VECSXP, out_size));
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
  Rf_unprotect(1);
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

[[cpp11::register]]
SEXP cpp_sset_df(SEXP x, SEXP indices){
  int xn = cpp_df_nrow(x);
  long long int llxn = xn;
  int ncols = Rf_length(x);
  int n = Rf_length(indices);
  int NP = 0;
  int zero_count = 0;
  int pos_count = 0;
  int oob_count = 0;
  int na_count = 0;
  int out_size;
  bool do_parallel = n >= CHEAPR_OMP_THRESHOLD;
  int n_cores = do_parallel ? num_cores() : 1;
  cpp11::function cheapr_sset = cpp11::package("cheapr")["sset"];
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, ncols));
  ++NP;
  // SEXP *p_out = VECTOR_PTR(out);

  // If indices is a special type of ALTREP compact int sequence, we can
  // Use a range-based subset instead

  if (is_compact_seq(indices)){

    // ALTREP integer sequence method

    SEXP seq_data = Rf_protect(compact_seq_data(indices));
    ++NP;
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
      ++NP;
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
      ++NP;
      SEXP indices2 = Rf_protect(cpp11::package("cheapr")["val_rm"](indices, r_zero));
      ++NP;
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
  ++NP;
  Rf_setAttrib(out, R_NamesSymbol, names);

  // list to data frame object
  SEXP df_str = Rf_protect(Rf_ScalarString(Rf_mkChar("data.frame")));
  ++NP;
  if (out_size > 0){
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 2));
    ++NP;
    INTEGER(row_names)[0] = NA_INTEGER;
    INTEGER(row_names)[1] = -out_size;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  } else {
    SEXP row_names = Rf_protect(Rf_allocVector(INTSXP, 0));
    ++NP;
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
  }
  Rf_classgets(out, df_str);
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
