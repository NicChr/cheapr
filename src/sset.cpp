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
  R_xlen_t istart1 = 0, istart2 = 0;
  R_xlen_t iend1 = 0, iend2 = 0;
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
      istart1 = 1;
      iend1 = std::abs(istart) - 1;
      istart2 = std::abs(iend) + 1;
      iend2 = n;
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
  bool *keep = (bool *) R_Calloc(n, bool);
  std::fill(keep, keep + n, true);

  if (xn > r_limits::r_int_max){
    SHIELD(exclude = vec::coerce_vec(exclude, REALSXP)); ++NP;
    double *p_excl = real_ptr(exclude);

    for (int j = 0; j < m; ++j) {
      if (is_r_na(p_excl[j])) continue;
      if (p_excl[j] > 0){
        R_Free(keep);
        YIELD(NP);
        Rf_error("Cannot mix positive and negative subscripts");
      }
      idx = -p_excl[j];
      // Check keep array for already assigned FALSE to avoid double counting
      if (idx > 0 && idx <= n && keep[idx - 1] == 1){
        keep[idx - 1] = false;
        ++exclude_count;
      }
    }
    out_size = n - exclude_count;
    SEXP out = SHIELD(new_vector<double>(out_size)); ++NP;
    double* RESTRICT p_out = real_ptr(out);
    while(k != out_size){
      if (keep[i] == true){
        p_out[k++] = i + 1;
      }
      ++i;
    }
    R_Free(keep);
    YIELD(NP);
    return out;
  } else {
    int *p_excl = integer_ptr(exclude);

    for (int j = 0; j < m; ++j) {
      if (is_r_na(p_excl[j])) continue;
      if (p_excl[j] > 0){
        R_Free(keep);
        YIELD(NP);
        Rf_error("Cannot mix positive and negative subscripts");
      }
      idx = -p_excl[j];
      // Check keep array for already assigned FALSE to avoid double counting
      if (idx > 0 && idx <= n && keep[idx - 1] == 1){
        keep[idx - 1] = false;
        ++exclude_count;
      }
    }
    out_size = n - exclude_count;
    SEXP out = SHIELD(vec::new_vector<int>(out_size)); ++NP;
    int* RESTRICT p_out = integer_ptr(out);
    while(k != out_size){
      if (keep[i++] == true){
        p_out[k++] = i;
      }
    }
    R_Free(keep);
    YIELD(NP);
    return out;
  }
}

// Cleans indices for subsetting
// Also returns metadata regarding final vec size and if indices should be
// checked (internal flag)

SEXP clean_indices(SEXP indices, SEXP x, bool count){
  R_xlen_t xn = vec::length(x);
  int32_t NP = 0;
  R_xlen_t zero_count = 0,
    pos_count = 0,
    oob_count = 0,
    na_count = 0,
    neg_count = 0;

  if (is_null(indices)){
    SHIELD(indices = compact_seq_len(xn)); ++NP;
  }

  R_xlen_t n = Rf_xlength(indices);
  int n_threads = calc_threads(n);

  int_fast64_t out_size = na::integer64;
  bool check_indices = true;
  SEXP clean_indices = r_null;

  if (TYPEOF(indices) == STRSXP){
    SEXP names = SHIELD(get_old_names(x)); ++NP;
    if (is_null(names)){
      YIELD(NP);
      Rf_error("Cannot subset on the names of an unnamed vector");
    }
    if (is_df(x)){
      YIELD(NP);
      Rf_error("Cannot subset rows of a data frame using a character vector");
    }
    SHIELD(indices = match(names, indices, na::integer)); ++NP;
  }

  if (is_compact_seq(indices)){
    clean_indices = indices;

    SEXP seq_data = SHIELD(compact_seq_data(indices)); ++NP;
    R_xlen_t from = real_ptr(seq_data)[0];
    R_xlen_t to = real_ptr(seq_data)[1];
    R_xlen_t by = real_ptr(seq_data)[2];
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
  } else if (xn > r_limits::r_int_max){

    SHIELD(clean_indices = vec::coerce_vec(indices, REALSXP)); ++NP;

    if (count){

      const double *pi = real_ptr_ro(clean_indices);

      // Counting the number of:
      // Zeroes
      // Out-of-bounds indices
      // Positive indices
      // NA indices
      // From this we can also work out the number of negatives

      if (n_threads > 1){
#pragma omp parallel for simd num_threads(n_threads) reduction(+:zero_count,pos_count,oob_count,na_count)
        for (int j = 0; j < n; ++j){
          zero_count += (pi[j] == 0);
          pos_count += (pi[j] > 0);
          oob_count += (std::fabs(pi[j]) > xn);
          na_count += is_r_na(pi[j]);
        }
      } else {
#pragma omp simd reduction(+:zero_count,pos_count,oob_count,na_count)
        for (int j = 0; j < n; ++j){
          zero_count += (pi[j] == 0);
          pos_count += (pi[j] > 0);
          oob_count += (std::fabs(pi[j]) > xn);
          na_count += is_r_na(pi[j]);
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

    SHIELD(clean_indices = vec::coerce_vec(indices, INTSXP)); ++NP;


    if (count){

      const int *pi = integer_ptr_ro(clean_indices);

      // Counting the number of:
      // Zeroes
      // Out-of-bounds indices
      // Positive indices
      // NA indices
      // From this we can also work out the number of negatives

      if (n_threads > 1){
#pragma omp parallel for simd num_threads(n_threads) reduction(+:zero_count,pos_count,oob_count,na_count)
        for (int j = 0; j < n; ++j){
          zero_count += (pi[j] == 0);
          pos_count += (pi[j] > 0);
          // oob_count counts NA as true so adjust after the fact
          oob_count += (std::llabs(pi[j]) > xn);
          na_count += is_r_na(pi[j]);
        }
      } else {
#pragma omp simd reduction(+:zero_count,pos_count,oob_count,na_count)
        for (int j = 0; j < n; ++j){
          zero_count += (pi[j] == 0);
          pos_count += (pi[j] > 0);
          // oob_count counts NA as true so adjust after the fact
          oob_count += (std::llabs(pi[j]) > xn);
          na_count += is_r_na(pi[j]);
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

  SEXP r_out_size = SHIELD(as_vector(r_cast<double>(out_size))); ++NP;
  SEXP r_check_indices = SHIELD(as_vector(check_indices)); ++NP;

  SEXP out = SHIELD(make_list(
    clean_indices,
    r_out_size,
    r_check_indices
  )); ++NP;

  YIELD(NP);
  return out;
}

// Clean indices
// Removes NAs, zeros and out-of-bounds
// Converts logical to integer locations
// Converts character to integer locations

SEXP clean_locs(SEXP locs, SEXP x){

  int_fast64_t xn = vec::length(x);

  int32_t NP = 0;

  if (is_null(locs)){
    SHIELD(locs = vec::new_vector<int>(0)); ++NP;
  }

  int_fast64_t n = Rf_xlength(locs);

  if (get_r_type(locs) == R_chr){
    SEXP names = SHIELD(get_old_names(x)); ++NP;
    if (is_null(names)){
      YIELD(NP);
      Rf_error("Cannot subset on the names of an unnamed vector");
    }
    if (is_df(x)){
      YIELD(NP);
      Rf_error("Cannot subset rows of a data frame using a character vector");
    }
    SHIELD(locs = match(names, locs, na::integer)); ++NP;
    SEXP na_int = SHIELD(as_vector(na::integer)); ++NP;
    SHIELD(locs = cpp_val_remove(locs, na_int, false)); ++NP;
    YIELD(NP);
    return locs;
  }

  if (get_r_type(locs) == R_lgl){

    if (Rf_length(locs) != xn){
      YIELD(NP);
      Rf_error("`length(i)` must match `length(x)` when `i` is a logical vector");
    }
    SHIELD(locs = cpp_which_(locs, false)); ++NP;
    YIELD(NP);
    return locs;
  }

  SHIELD(locs = cast<r_integers_t>(locs, r_null)); ++NP;

  const int *p_locs = integer_ptr_ro(locs);

  int zero_count = 0,
    pos_count = 0,
    oob_count = 0,
    na_count = 0,
    neg_count = 0;

  int loc;
  for (int i = 0; i < n; ++i){
    loc = p_locs[i];
    zero_count += (loc == 0);
    pos_count += (loc > 0);
    oob_count += std::abs(static_cast<int_fast64_t>(loc)) > xn;
    na_count += is_r_na(loc);
  }

  oob_count = oob_count - na_count;
  neg_count = n - pos_count - zero_count - na_count;

  if ( (pos_count > 0 && neg_count > 0) ||
       (neg_count > 0 && na_count > 0)){
    YIELD(NP);
    Rf_error("Cannot mix positive and negative indices");
  }

  if (neg_count > 0){
    SHIELD(locs = exclude_locs(locs, xn)); ++NP;
    YIELD(NP);
    return locs;
  }
  if (zero_count > 0 || oob_count > 0 || na_count > 0){
    int out_size = pos_count - oob_count;
    SEXP out = SHIELD(vec::new_vector<int>(out_size)); ++NP;
    int* RESTRICT p_out = integer_ptr(out);

    int k = 0;
    for (int i = 0; i < n; ++i){
      if (between<int>(p_locs[i], 1, xn)){
        p_out[k++] = p_locs[i];
      }
    }
    YIELD(NP);
    return out;
  }

  YIELD(NP);
  return locs;
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
  R_xlen_t out_size, istart1 = 0, istart2 = 0;
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

  // Out-of-bounds
  R_xlen_t n_oob = std::max(( by > 0) ? iend - n : istart - n, (R_xlen_t) 0);
  // Adjustment for when all values are oob
  if ( ( by > 0 && istart > n ) || (by < 0 && iend > n)){
    n_oob = out_size;
  }
  R_xlen_t in_bounds_size = std::max(out_size - n_oob, (R_xlen_t) 0);

  SEXP out = r_null;

  switch ( TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    const int *p_x = integer_ptr_ro(x);
    out = SHIELD(internal::new_vec(TYPEOF(x), out_size)); ++NP;
    int* RESTRICT p_out = integer_ptr(out);
    if (double_loop){
      std::copy_n(&p_x[istart1 - 1], iend1 - istart1 + 1, &p_out[0]);
      std::copy_n(&p_x[istart2 - 1], iend2 - istart2 + 1, &p_out[iend1 - istart1 + 1]);
    } else {
      if (by > 0){
        std::copy_n(&p_x[istart - 1], in_bounds_size, &p_out[0]);
        std::fill(p_out + in_bounds_size, p_out + in_bounds_size + n_oob, na::integer);
      } else {
        std::fill(p_out, p_out + n_oob, na::integer);
        OMP_SIMD
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i) p_out[istart - i - 1] = p_x[i];
      }
    }
    break;
  }
  case REALSXP: {
    const double *p_x = real_ptr_ro(x);
    out = SHIELD(new_vector<double>(out_size)); ++NP;
    double* RESTRICT p_out = real_ptr(out);
    if (double_loop){
      std::copy_n(&p_x[istart1 - 1], iend1 - istart1 + 1, &p_out[0]);
      std::copy_n(&p_x[istart2 - 1], iend2 - istart2 + 1, &p_out[iend1 - istart1 + 1]);
    } else {
      if (by > 0){
        std::copy_n(&p_x[istart - 1], in_bounds_size, &p_out[0]);
        std::fill(p_out + in_bounds_size, p_out + in_bounds_size + n_oob, na::real);
      } else {
        std::fill(p_out, p_out + n_oob, na::real);
        OMP_SIMD
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i) p_out[istart - i - 1] = p_x[i];
      }
    }
    break;
  }
  case STRSXP: {
    const r_string_t *p_x = string_ptr_ro(x);
    out = SHIELD(new_vector<r_string_t>(out_size)); ++NP;
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
          SET_STRING_ELT(out, in_bounds_size + i, na::string);
        }
      } else {
        for (R_xlen_t i = 0; i < n_oob; ++i){
          SET_STRING_ELT(out, i, na::string);
        }
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i){
          SET_STRING_ELT(out, istart - i - 1, p_x[i]);
        }
      }
    }
    break;
  }
  case CPLXSXP: {
    const r_complex_t *p_x = complex_ptr_ro(x);
    out = SHIELD(new_vector<r_complex_t>(out_size)); ++NP;
    r_complex_t* RESTRICT p_out = complex_ptr(out);
    if (double_loop){
      std::copy_n(&p_x[istart1 - 1], iend1 - istart1 + 1, &p_out[0]);
      std::copy_n(&p_x[istart2 - 1], iend2 - istart2 + 1, &p_out[iend1 - istart1 + 1]);
    } else {
      if (by > 0){
        std::copy_n(&p_x[istart - 1], in_bounds_size, &p_out[0]);
        std::fill(p_out + in_bounds_size, p_out + in_bounds_size + n_oob, na::complex);
      } else {
        std::fill(p_out, p_out + n_oob, na::complex);
        OMP_SIMD
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i) set_value(p_out, istart - i - 1, p_x[i]);
      }
    }
    break;
  }
  case RAWSXP: {
    const r_byte_t *p_x = raw_ptr_ro(x);
    out = SHIELD(new_vector<r_byte_t>(out_size)); ++NP;
    r_byte_t* RESTRICT p_out = raw_ptr(out);
    if (double_loop){
      std::copy_n(&p_x[istart1 - 1], iend1 - istart1 + 1, &p_out[0]);
      std::copy_n(&p_x[istart2 - 1], iend2 - istart2 + 1, &p_out[iend1 - istart1 + 1]);
    } else {
      if (by > 0){
        std::copy_n(&p_x[istart - 1], in_bounds_size, &p_out[0]);
        std::fill(p_out + in_bounds_size, p_out + in_bounds_size + n_oob, na::raw);
      } else {
        std::fill(p_out, p_out + n_oob, na::raw);
        OMP_SIMD
        for (R_xlen_t i = istart - 1 - n_oob; i >= iend - 1; --i) set_value(p_out, istart - i - 1, p_x[i]);
      }
    }
    break;
  }
  case VECSXP: {
    const SEXP *p_x = list_ptr_ro(x);
    out = SHIELD(new_list(out_size)); ++NP;
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

  SEXP out = r_null;
  int xtype = TYPEOF(x);

  if (check){

    if (Rf_xlength(x) > r_limits::r_int_max){

      int_fast64_t xn = Rf_xlength(x);

      int_fast64_t
      n = Rf_xlength(indices), k = 0, j;

      const double* pind = real_ptr_ro(indices);

      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(r_null);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int* p_x = integer_ptr_ro(x);
        out = SHIELD(internal::new_vec(xtype, n));
        int* RESTRICT p_out = integer_ptr(out);

        for (int_fast64_t i = 0; i < n; ++i){
          j = pind[i];
          if (j < 0){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0){
            p_out[k++] = (is_r_na(pind[i]) || j > xn) ? na::integer : p_x[j - 1];
          }
        }
        break;
      }
      case REALSXP: {
        const double* p_x = real_ptr_ro(x);
        out = SHIELD(new_vector<double>(n));
        double* RESTRICT p_out = real_ptr(out);
        for (int_fast64_t i = 0; i < n; ++i){
          j = pind[i];
          if (j < 0){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0){
            p_out[k++] = (is_r_na(pind[i]) || j > xn) ? na::real : p_x[j - 1];
          }
        }
        break;
      }
      case STRSXP: {
        const r_string_t *p_x = string_ptr_ro(x);
        out = SHIELD(new_vector<r_string_t>(n));
        for (int_fast64_t i = 0; i < n; ++i){
          j = pind[i];
          if (j < 0){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0){
            set_value(out, k++, (is_r_na(pind[i]) || j > xn) ? na::string : p_x[j - 1]);
          }
        }
        break;
      }
      case CPLXSXP: {
        const r_complex_t* p_x = complex_ptr_ro(x);
        out = SHIELD(new_vector<r_complex_t>(n));
        r_complex_t* RESTRICT p_out = complex_ptr(out);
        for (int_fast64_t i = 0; i < n; ++i){

          j = pind[i];

          if (j < 0){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0){
            if (is_r_na(pind[i]) || j > xn){
              set_value(p_out, k++, na::complex);
            } else {
              set_value(p_out, k++, p_x[j - 1]);
            }
          }

        }
        break;
      }
      case RAWSXP: {
        const r_byte_t *p_x = raw_ptr_ro(x);
        out = SHIELD(new_vector<r_byte_t>(n));
        for (int_fast64_t i = 0; i < n; ++i){
          j = pind[i];
          if (j < 0){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0){
            set_value<r_byte_t>(out, k++, (is_r_na(pind[i]) || j > xn) ? na::raw : p_x[j - 1]);
          }
        }
        break;
      }
      case VECSXP: {
        const SEXP *p_x = list_ptr_ro(x);
        out = SHIELD(new_list(n));
        for (int_fast64_t i = 0; i < n; ++i){
          j = pind[i];
          if (j < 0){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0){
            SET_VECTOR_ELT(out, k++, (is_r_na(pind[i]) || j > xn) ? r_null : p_x[j - 1]);
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
        na_val = na::integer,
        j;

      const int *pind = integer_ptr_ro(indices);
      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(r_null);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int* p_x = integer_ptr_ro(x);
        out = SHIELD(new_vec(xtype, n));
        int* RESTRICT p_out = integer_ptr(out);
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
            p_out[k++] = na::integer;
          }
        }
        break;
      }
      case REALSXP: {
        const double* p_x = real_ptr_ro(x);
        out = SHIELD(new_vector<double>(n));
        double* RESTRICT p_out = real_ptr(out);
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
            p_out[k++] = na::real;
          }
        }
        break;
      }
      case STRSXP: {
        const r_string_t *p_x = string_ptr_ro(x);
        out = SHIELD(new_vector<r_string_t>(n));
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            set_value(out, k++, p_x[--j]);
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            set_value(out, k++, na::string);
          }
        }
        break;
      }
      case CPLXSXP: {
        const r_complex_t* p_x = complex_ptr_ro(x);
        out = SHIELD(new_vector<r_complex_t>(n));
        r_complex_t* RESTRICT p_out = complex_ptr(out);
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            set_value(p_out, k++, p_x[--j]);
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            set_value(p_out, k++, na::complex);
          }
        }
        break;
      }
      case RAWSXP: {
        const r_byte_t *p_x = raw_ptr_ro(x);
        out = SHIELD(new_vector<r_byte_t>(n));
        for (unsigned int i = 0; i < n; ++i){
          j = pind[i];
          if (between<unsigned int>(j, 1U, xn)){
            set_value<r_byte_t>(out, k++, p_x[--j]);
            // If j > n_val then it is a negative 32-bit integer
          } else if (j > na_val){
            SEXP new_indices = SHIELD(exclude_locs(indices, xn));
            SEXP out2 = SHIELD(sset_vec(x, new_indices, false));
            YIELD(3);
            return out2;
          } else if (j != 0U){
            set_value<r_byte_t>(out, k++, r_byte_t{0});
          }
        }
        break;
      }
      case VECSXP: {
        const SEXP *p_x = list_ptr_ro(x);
        out = SHIELD(new_list(n));
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
            SET_VECTOR_ELT(out, k++, r_null);
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

    if (Rf_xlength(x) > r_limits::r_int_max){

      int_fast64_t n = Rf_xlength(indices);

      const double *pind = real_ptr_ro(indices);

      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(r_null);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int *p_x = integer_ptr_ro(x);
        out = SHIELD(internal::new_vec(xtype, n));
        int* RESTRICT p_out = integer_ptr(out);
        OMP_SIMD
        for (int_fast64_t i = 0; i < n; ++i){
          p_out[i] = p_x[static_cast<int_fast64_t>(pind[i] - 1.0)];
        }
        break;
      }
      case REALSXP: {
        const double *p_x = real_ptr_ro(x);
        out = SHIELD(new_vector<double>(n));
        double* RESTRICT p_out = real_ptr(out);
        OMP_SIMD
        for (int_fast64_t i = 0; i < n; ++i){
          p_out[i] = p_x[static_cast<int_fast64_t>(pind[i] - 1.0)];
        }
        break;
      }
      case STRSXP: {
        const r_string_t *p_x = string_ptr_ro(x);
        out = SHIELD(new_vector<r_string_t>(n));
        for (int_fast64_t i = 0; i < n; ++i){
          set_value(out, i, p_x[static_cast<int_fast64_t>(pind[i] - 1.0)]);
        }
        break;
      }
      case CPLXSXP: {
        const r_complex_t *p_x = complex_ptr_ro(x);
        out = SHIELD(new_vector<r_complex_t>(n));
        r_complex_t* RESTRICT p_out = complex_ptr(out);
        for (int_fast64_t i = 0; i < n; ++i){
          set_value(p_out, i, p_x[static_cast<int_fast64_t>(pind[i] - 1.0)]);
        }
        break;
      }
      case RAWSXP: {
        const r_byte_t *p_x = raw_ptr_ro(x);
        out = SHIELD(new_vector<r_byte_t>(n));
        for (int_fast64_t i = 0; i < n; ++i) set_value<r_byte_t>(out, i, p_x[static_cast<int_fast64_t>(pind[i] - 1.0)]);
        break;
      }
      case VECSXP: {
        const SEXP *p_x = list_ptr_ro(x);
        out = SHIELD(new_list(n));
        for (int_fast64_t i = 0; i < n; ++i){
          SET_VECTOR_ELT(out, i, p_x[static_cast<int_fast64_t>(pind[i] - 1.0)]);
        }
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

      const int *pind = integer_ptr_ro(indices);
      switch ( xtype ){

      case NILSXP: {
        out = SHIELD(r_null);
        break;
      }
      case LGLSXP:
      case INTSXP: {
        const int *p_x = integer_ptr_ro(x);
        out = SHIELD(internal::new_vec(xtype, n));
        int* RESTRICT p_out = integer_ptr(out);
        for (int i = 0; i != n; ++i){
          p_out[i] = p_x[pind[i] - 1];
        }
        break;
      }
      case REALSXP: {
        const double *p_x = real_ptr_ro(x);
        out = SHIELD(new_vector<double>(n));
        double* RESTRICT p_out = real_ptr(out);
        for (int i = 0; i != n; ++i){
          p_out[i] = p_x[pind[i] - 1];
        }
        break;
      }
      case STRSXP: {
        const r_string_t *p_x = string_ptr_ro(x);
        out = SHIELD(new_vector<r_string_t>(n));
        for (int i = 0; i != n; ++i){
          set_value(out, i, p_x[pind[i] - 1]);
        }
        break;
      }
      case CPLXSXP: {
        const r_complex_t *p_x = complex_ptr_ro(x);
        out = SHIELD(new_vector<r_complex_t>(n));
        r_complex_t* RESTRICT p_out = complex_ptr(out);
        for (int i = 0; i != n; ++i){
          set_value(p_out, i, p_x[pind[i] - 1]);
        }
        break;
      }
      case RAWSXP: {
        const r_byte_t *p_x = raw_ptr_ro(x);
        out = SHIELD(new_vector<r_byte_t>(n));
        r_byte_t* RESTRICT p_out = raw_ptr(out);
        for (int i = 0; i != n; ++i){
          set_value(p_out, i, p_x[pind[i] - 1]);
        }
        break;
      }
      case VECSXP: {
        const SEXP *p_x = list_ptr_ro(x);
        out = SHIELD(new_list(n));
        for (int i = 0; i != n; ++i){
          SET_VECTOR_ELT(out, i, p_x[pind[i] - 1]);
        }
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
    out = r_null;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
    int* RESTRICT p_out = integer_ptr(out);
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
    double* RESTRICT p_out = real_ptr(out);
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
    const r_string_t *p_out = string_ptr_ro(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      SEXP left = p_out[i];
      set_value(out, i, p_out[k]);
      set_value(out, k, left);
    }
    break;
  }
  case CPLXSXP: {
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
    r_complex_t *p_out = complex_ptr(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      r_complex_t left = p_out[i];
      set_value<r_complex_t>(out, i, p_out[k]);
      set_value<r_complex_t>(out, k, left);
    }
    break;
  }
  case RAWSXP: {
    out = SHIELD(set ? x : cpp_semi_copy(x)); ++NP;
    r_byte_t *p_out = raw_ptr(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      r_byte_t left = p_out[i];
      set_value<r_byte_t>(out, i, p_out[k]);
      set_value<r_byte_t>(out, k, left);
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
  SEXP names = SHIELD(get_old_names(x));
  SHIELD(names = rev(names, set && !ALTREP(x)));
  set_old_names(out, names);
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
    n_rows = df::nrow(x),
    n_locs = Rf_length(locs);

  // Flag to check indices
  bool check = true;

  SEXP names = SHIELD(get_old_names(x)); ++NP;

  SEXP cols;
  int loc_type = TYPEOF(locs);

  switch(loc_type){
  case NILSXP: {
    // If NULL then select all cols
    cols = SHIELD(cpp_seq_len(n_cols)); ++NP;
    n_locs = n_cols;
    check = false;
    break;
  }
  case STRSXP: {
    cols = SHIELD(match(names, locs, na::integer)); ++NP;
    int *match_locs = integer_ptr(cols);
    if (cpp_any_na(cols, false)){
      for (int i = 0; i < Rf_length(cols); ++i){
        if (is_r_na(match_locs[i])){
          const char *bad_loc = utf8_char(get_value<r_string_t>(locs, i));
          YIELD(NP);
          Rf_error("Column %s does not exist", bad_loc);
        }
      }
    }
    check = false;
    break;
  }
  case LGLSXP: {
    // If logical then find locs using `which_()`
    if (Rf_length(locs) != n_cols){
    YIELD(NP);
    Rf_error("`length(j)` must match `ncol(x)` when `j` is a logical vector");
  }
    cols = SHIELD(cpp_which_(locs, false)); ++NP;
    n_locs = Rf_length(cols);
    check = false;
    break;
  }
  default: {
    // Catch-all make sure cols is an int vector
    cols = SHIELD(cast<r_integers_t>(locs, r_null)); ++NP;
    break;
  }
  }

  // Negative subscripting
  if (n_locs > 0 && !is_r_na(integer_ptr(cols)[0]) && integer_ptr(cols)[0] < 0){
    SHIELD(cols = exclude_locs(cols, n_cols)); ++NP;
    n_locs = Rf_length(cols);
    check = false;
  }

  const int *p_cols = integer_ptr_ro(cols);

  SEXP out = SHIELD(new_list(n_locs)); ++NP;
  SEXP out_names = SHIELD(new_vector<r_string_t>(n_locs)); ++NP;

  const SEXP *p_x = list_ptr_ro(x);
  const r_string_t *p_names = string_ptr_ro(names);
  int k = 0;
  int col;

  if (check){
    for (int i = 0; i < n_locs; ++i) {
      col = p_cols[i];
      if (is_r_na(col)){
        YIELD(NP);
        Rf_error("Cannot select `NA` column locations in %s", __func__);
      } else if (col >= 1 && col <= n_cols){
        --col;
        SET_VECTOR_ELT(out, k, p_x[col]);
        set_value(out_names, k, p_names[col]);
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
      set_value(out_names, i, p_names[col - 1]);
    }
  }

  if (check && k != n_locs){
    SHIELD(out = Rf_lengthgets(out, k)); ++NP;
    SHIELD(out_names = Rf_lengthgets(out_names, k)); ++NP;
  }

  // Make a plain data frame
  df::set_row_names(out, n_rows);
  SEXP df_cls = SHIELD(as_vector("data.frame")); ++NP;
  attr::set_old_class(out, df_cls);
  set_old_names(out, out_names);
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
  const SEXP *p_x = list_ptr_ro(x);
  SEXP out = SHIELD(new_list(ncols)); ++NP;

  // Clean indices and get metadata

  int out_size;
  bool check_indices;

  if (check){
    SEXP clean_indices_info = SHIELD(clean_indices(indices, x, true)); ++NP;
    SHIELD(indices = VECTOR_ELT(clean_indices_info, 0)); ++NP;
    out_size = real_ptr(VECTOR_ELT(clean_indices_info, 1))[0];
    check_indices = logical_ptr(VECTOR_ELT(clean_indices_info, 2))[0];
  } else {
    out_size = Rf_length(indices);
    check_indices = false;
  }

  // Subset columns

  for (int j = 0; j < ncols; ++j){
    SET_VECTOR_ELT(out, j, cpp_sset(p_x[j], indices, check_indices));
  }

  SEXP names = SHIELD(get_old_names(x)); ++NP;
  set_old_names(out, names);

  // list to data frame object
  df::set_row_names(out, out_size);
  SEXP df_cls = SHIELD(as_vector("data.frame")); ++NP;
  attr::set_old_class(out, df_cls);
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
SEXP cpp_sset2(SEXP x, SEXP i, SEXP j, bool check, SEXP args){

  int32_t NP = 0;

  SEXP out = r_null;

  if (cheapr_is_simple_vec(x)){

    if (!is_null(j)){
      YIELD(NP);
      Rf_error("`x` does not have cols, please leave `j` as `NULL`");
    }

    if (check){
      SEXP indices_metadata = SHIELD(clean_indices(i, x, false)); ++NP;
      SHIELD(i = VECTOR_ELT(indices_metadata, 0)); ++NP;
      check = logical_ptr(VECTOR_ELT(indices_metadata, 2))[0];
    }

    SEXP names = r_null;

    if (is_compact_seq(i)){
      SEXP seq_data = SHIELD(compact_seq_data(i)); ++NP;
      const double *p_data = real_ptr_ro(seq_data);
      SHIELD(out = cpp_sset_range(x, p_data[0], p_data[1], p_data[2])); ++NP;

      // Subset names
      SHIELD(names = get_old_names(x)); ++NP;
      SHIELD(names = cpp_sset_range(names, p_data[0], p_data[1], p_data[2])); ++NP;
    } else {

      SHIELD(out = sset_vec(x, i, check)); ++NP;

      // Subset names
      SHIELD(names = get_old_names(x)); ++NP;
      SHIELD(names = sset_vec(names, i, check)); ++NP;
    }
    Rf_copyMostAttrib(x, out);
    set_old_names(out, names);
  } else if (is_df(x)){
    if (is_bare_df(x) || is_bare_tbl(x)){
      SHIELD(out = cpp_df_subset(x, i, j, check)); ++NP;
    } else {
      if (Rf_length(args) == 0){
        SHIELD(out = eval_pkg_fun("cheapr_sset", "cheapr", env::base_env, x, i, j)); ++NP;
      } else {
        SEXP usual_args = SHIELD(make_list(
          arg("x") = x,
          arg("i") = i,
          arg("j") = j
        )); ++NP;

        // Combine all args into one list
        SEXP all_args = SHIELD(combine(usual_args, args)); ++NP;
        SEXP cheapr_sset_fn = SHIELD(fn::find_pkg_fun("cheapr_sset", "cheapr", true)); ++NP;
        SHIELD(out = eval_pkg_fun("do.call", "base", env::base_env, cheapr_sset_fn, all_args)); ++NP;
      }
    }
  } else {
    // Fall-back to `cheapr_sset()` S3 methods
    if (Rf_length(args) == 0){
      SHIELD(out = eval_pkg_fun("cheapr_sset", "cheapr", env::base_env, x, i, j)); ++NP;
    } else {
      SEXP usual_args = SHIELD(make_list(
        arg("x") = x,
        arg("i") = i,
        arg("j") = j
      )); ++NP;

      // Combine all args into one list
      SEXP all_args = SHIELD(combine(usual_args, args)); ++NP;
      SEXP cheapr_sset_fn = SHIELD(fn::find_pkg_fun("cheapr_sset", "cheapr", true)); ++NP;
      SHIELD(out = eval_pkg_fun("do.call", "base", env::base_env, cheapr_sset_fn, all_args)); ++NP;
    }
  }
  YIELD(NP);
  return out;
}

// Keep this as it's a handy C function to export
// also keep for legacy reasons as fastplyr directly uses it
[[cpp11::register]]
SEXP cpp_sset(SEXP x, SEXP indices, bool check){
  return cpp_sset2(x, indices, r_null, check, r_null);
}

// scalar subset
SEXP slice_loc(SEXP x, R_xlen_t i){

  if (i < 0){
    Rf_error("`i` must be >= 0");
  }

  if (vec::is_object(x)){
    SEXP loc = SHIELD(as_vector(i));
    SEXP out = SHIELD(cpp_sset(x, loc, true));
    YIELD(2);
    return out;
  }

  if (i == 0){
    return internal::new_vec(TYPEOF(x), 0);
  }

  if (i > Rf_xlength(x)){
    return cpp_na_init(x, 1);
  }

  --i;

  switch ( TYPEOF(x) ){

  case NILSXP: {
    return r_null;
  }
  case LGLSXP: {
    return as_vector(logical_ptr(x)[i]);
  }
  case INTSXP: {
    return as_vector(integer_ptr(x)[i]);
  }
  case REALSXP: {
    return as_vector(real_ptr(x)[i]);
  }
  case STRSXP: {
    return as_vector(STRING_ELT(x, i));
  }
  case CPLXSXP: {
    return as_vector(complex_ptr(x)[i]);
  }
  case RAWSXP: {
    return as_vector(raw_ptr(x)[i]);
  }
  case VECSXP: {
    return VECTOR_ELT(x, i);
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}
