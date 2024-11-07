#include "cheapr.h"

// Internal fast check for valid order values (here they are 0-indexed)
// Where rng = length(x) - 1
void check_order_value(unsigned int x, unsigned int rng, int n_prot){
  if (!(x <= rng)){
    Rf_unprotect(n_prot);
    Rf_error("order must be an integer vector of unique values between 1 and length(x)");
  }
}

[[cpp11::register]]
SEXP cpp_lag(SEXP x, R_xlen_t k, SEXP fill, bool set, bool recursive) {
  R_xlen_t size = Rf_xlength(x);
  R_xlen_t fill_size = Rf_xlength(fill);
  int NP = 0;
  if (fill_size > 1){
    Rf_error("fill size must be NULL or length 1");
  }
  bool set_and_altrep = set && ALTREP(x);
  SEXP out;
  if (ALTREP(x)){
    set = true;
  }
  SEXP xvec = Rf_protect(altrep_materialise(x)); ++NP;
  if (set_and_altrep){
    Rf_warning("Cannot lag an ALTREP by reference, a copy has been made.\n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- lag_(x, set = TRUE)`");
  }
  switch(CHEAPR_TYPEOF(xvec)){
  case NILSXP: {
    out = R_NilValue;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    int fill_value = NA_INTEGER;
    if (fill_size >= 1){
      fill_value = Rf_asInteger(fill);
    }
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec));
    ++NP;
    int *p_out = INTEGER(out);
    int *p_x = INTEGER(xvec);
    if (k >= 0){
      memmove(&p_out[k], &p_x[0], (size - k) * sizeof(int));
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < k; ++i) p_out[i] = fill_value;
    } else {
      memmove(&p_out[0], &p_x[-k], (size + k) * sizeof(int));
      OMP_FOR_SIMD
      for (R_xlen_t i = size - 1; i >= size + k; --i) p_out[i] = fill_value;
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++NP;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++NP;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    long long int fill_value = NA_INTEGER64;
    if (fill_size >= 1){
      SEXP temp_fill = Rf_protect(coerce_vector(fill, CHEAPR_INT64SXP)); ++NP;
      fill_value = INTEGER64_PTR(temp_fill)[0];
    }
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec)); ++NP;
    long long int *p_out = INTEGER64_PTR(out);
    long long int *p_x = INTEGER64_PTR(xvec);

    // sizeof(long long int) == sizeof(double), both are 8 bytes

    if (k >= 0){
      memmove(&p_out[k], &p_x[0], (size - k) * sizeof(long long int));
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < k; ++i) p_out[i] = fill_value;
    } else {
      memmove(&p_out[0], &p_x[-k], (size + k) * sizeof(long long int));
      OMP_FOR_SIMD
      for (R_xlen_t i = size - 1; i >= size + k; --i) p_out[i] = fill_value;
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++NP;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++NP;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case REALSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    double fill_value = NA_REAL;
    if (fill_size >= 1){
      fill_value = Rf_asReal(fill);
    }
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec)); ++NP;
    double *p_out = REAL(out);
    double *p_x = REAL(xvec);
    if (k >= 0){
      memmove(&p_out[k], &p_x[0], (size - k) * sizeof(double));
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < k; ++i) p_out[i] = fill_value;
    } else {
      memmove(&p_out[0], &p_x[-k], (size + k) * sizeof(double));
      OMP_FOR_SIMD
      for (R_xlen_t i = size - 1; i >= size + k; --i) p_out[i] = fill_value;
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++NP;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++NP;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case CPLXSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP fill_sexp = Rf_protect(Rf_allocVector(CPLXSXP, 1)); ++NP;
    Rcomplex *p_fill = COMPLEX(fill_sexp);
    p_fill[0].i = NA_REAL;
    p_fill[0].r = NA_REAL;
    Rcomplex fill_value = fill_size >= 1 ? Rf_asComplex(fill) : COMPLEX(fill_sexp)[0];
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec)); ++NP;
    Rcomplex *p_out = COMPLEX(out);
    Rcomplex *p_x = COMPLEX(xvec);
    if (k >= 0){
      memmove(&p_out[k], &p_x[0], (size - k) * sizeof(Rcomplex));
      for (R_xlen_t i = 0; i < k; ++i) SET_COMPLEX_ELT(out, i, fill_value);
    } else {
      memmove(&p_out[0], &p_x[-k], (size + k) * sizeof(Rcomplex));
      OMP_FOR_SIMD
      for (R_xlen_t i = size - 1; i >= size + k; --i) SET_COMPLEX_ELT(out, i, fill_value);
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++NP;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++NP;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case STRSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP fill_char = Rf_protect(fill_size >= 1 ? Rf_asChar(fill) : NA_STRING);
    ++NP;
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec)); ++NP;
    const SEXP *p_out = STRING_PTR_RO(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(STRSXP, std::abs(k)));
        ++NP;
        SEXP tempv = Rf_protect(Rf_allocVector(STRSXP, 1));
        ++NP;
        const SEXP* __restrict__ p_lag = STRING_PTR_RO(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            SET_STRING_ELT(lag_temp, i, p_out[i]);
            SET_STRING_ELT(out, i, fill_char);
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            SET_STRING_ELT(tempv, 0, p_lag[tempi]);
            SET_STRING_ELT(lag_temp, tempi, p_out[i]);
            SET_STRING_ELT(out, i, STRING_ELT(tempv, 0));
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_STRING_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_STRING_ELT(out, i, fill_char);
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            SET_STRING_ELT(tempv, 0, p_lag[tempi]);
            SET_STRING_ELT(lag_temp, tempi, p_out[i]);
            SET_STRING_ELT(out, i, STRING_ELT(tempv, 0));
          }
        }
      }
    } else {
      const SEXP *p_x = STRING_PTR_RO(xvec);
      if (k >= 0){
        for (R_xlen_t i = 0; i < k; ++i) SET_STRING_ELT(out, i, fill_char);
        for (R_xlen_t i = k; i < size; ++i) SET_STRING_ELT(out, i, p_x[i - k]);
      } else {
        for (R_xlen_t i = size - 1; i >= size + k; --i) SET_STRING_ELT(out, i, fill_char);
        for (R_xlen_t i = size + k - 1; i >= 0; --i) SET_STRING_ELT(out, i, p_x[i - k]);
      }
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++NP;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++NP;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case RAWSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP raw_sexp = Rf_protect(Rf_coerceVector(fill, RAWSXP));
    Rbyte fill_raw = fill_size == 0 ? RAW(Rf_ScalarRaw(0))[0] : RAW(raw_sexp)[0];
    ++NP;
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec));
    ++NP;
    Rbyte *p_out = RAW(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(RAWSXP, std::abs(k)));
        ++NP;
        Rbyte *p_lag = RAW(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            SET_RAW_ELT(lag_temp, i, p_out[i]);
            SET_RAW_ELT(out, i, fill_raw);
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            Rbyte tempv = p_lag[tempi];
            SET_RAW_ELT(lag_temp, tempi, p_out[i]);
            SET_RAW_ELT(out, i, tempv);
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_RAW_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_RAW_ELT(out, i, fill_raw);
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            Rbyte tempv = p_lag[tempi];
            SET_RAW_ELT(lag_temp, tempi, p_out[i]);
            SET_RAW_ELT(out, i, tempv);
          }
        }
      }
    } else {
      const Rbyte *p_x = RAW_RO(xvec);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          SET_RAW_ELT(out, i, i >= k ? p_x[i - k] : fill_raw);
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          SET_RAW_ELT(out, i, (i - size) < k ? p_x[i - k] : fill_raw);
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++NP;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++NP;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    const SEXP *p_x = VECTOR_PTR_RO(xvec);
    out = Rf_protect(Rf_allocVector(VECSXP, size));
    ++NP;
    SHALLOW_DUPLICATE_ATTRIB(out, xvec);
    for (R_xlen_t i = 0; i < size; ++i){
      SET_VECTOR_ELT(out, i, cpp_lag(p_x[i], k, fill, set, true));
    }
  } else {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP fill_value = Rf_protect(Rf_coerceVector(fill_size >= 1 ? fill : R_NilValue, VECSXP));
    ++NP;
    out = Rf_protect(set ? xvec : Rf_allocVector(VECSXP, size));
    ++NP;
    const SEXP *p_out = VECTOR_PTR_RO(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(VECSXP, std::abs(k)));
        ++NP;
        SEXP tempv = Rf_protect(Rf_allocVector(VECSXP, 1));
        ++NP;
        SEXP* __restrict__ p_lag = VECTOR_PTR(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            SET_VECTOR_ELT(lag_temp, i, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(fill_value, 0));
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            SET_VECTOR_ELT(tempv, 0, p_lag[tempi]);
            SET_VECTOR_ELT(lag_temp, tempi, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(tempv, 0));
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_VECTOR_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(fill_value, 0));
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            SET_VECTOR_ELT(tempv, 0, p_lag[tempi]);
            SET_VECTOR_ELT(lag_temp, tempi, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(tempv, 0));
          }
        }
      }
    } else {
      const SEXP *p_x = VECTOR_PTR_RO(xvec);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          SET_VECTOR_ELT(out, i, i >= k ? p_x[i - k] : VECTOR_ELT(fill_value, 0));
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          SET_VECTOR_ELT(out, i, (i - size) < k ? p_x[i - k] : VECTOR_ELT(fill_value, 0));
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++NP;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++NP;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
  }
  break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(xvec)));
  }
  }
  Rf_unprotect(NP);
  return out;
}

[[cpp11::register]]
SEXP cpp_lag2(SEXP x, SEXP lag, SEXP order, SEXP run_lengths, SEXP fill, bool recursive){
  int o_size = Rf_length(order);
  int rl_size = Rf_length(run_lengths);
  int lag_size = Rf_length(lag);
  int fill_size = Rf_length(fill);
  int NP = 0;
  if (fill_size > 1){
    Rf_error("fill size must be NULL or length 1");
  }
  if (lag_size < 1){
    Rf_error("lag must be a non-zero length integer vector");
  }
  bool has_order = !Rf_isNull(order);
  bool has_rl = !Rf_isNull(run_lengths);
  bool recycle_lag = lag_size != 1;

  // When order is NULL we run through x from left to right (as usual)
  // When run_lengths is NULL we run through x without resetting
  // To do this properly we use dummy vectors so that
  // we can statically assign int pointers
  // While keeping order and run_lengths as NULL (if they are NULL)
  // This is mainly done this way because lag2_ is recursive and
  // hence order/run_lengths should remain NULL throughout each recursion
  // if they are NULL

  SEXP dummy_vec1 = Rf_protect(Rf_allocVector(INTSXP, 0));
  ++NP;
  SEXP dummy_vec2 = Rf_protect(Rf_allocVector(INTSXP, 0));
  ++NP;
  Rf_protect(lag = Rf_coerceVector(lag, INTSXP));
  ++NP;
  Rf_protect(order = has_order ? Rf_coerceVector(order, INTSXP) : R_NilValue);
  ++NP;
  Rf_protect(run_lengths = has_rl ? Rf_coerceVector(run_lengths, INTSXP) : R_NilValue);
  ++NP;
  // std::variant<int*, double*> p_o;
  // typedef std::conditional<true, int*, double*>::type p_o;
  int *p_o = INTEGER(has_order ? order : dummy_vec1);
  int *p_rl = INTEGER(has_rl ? run_lengths : dummy_vec2);
  int *p_lag = INTEGER(lag);
  int rl; // Run-length
  int run_start = 0; // Start index of current run
  int run_end = 0; // End index of current run
  int oi; // Indices (specified by order vector) to lag
  int k; // Lag
  int lag1 = p_lag[0];
  // Manually set run rl size to 1 if run_lengths = NULL (checked prior)
  if (!has_rl) rl_size = 1;

  // Macro to lag names(x)
#define CHEAPR_LAG_R_NAMES                                                                                  \
  if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){                                                          \
    SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));                                            \
    ++NP;                                                                                                   \
    SEXP new_names = Rf_protect(cpp_lag2(old_names, lag, order, run_lengths, R_NilValue, recursive));       \
    ++NP;                                                                                                   \
    Rf_setAttrib(out, R_NamesSymbol, new_names);                                                            \
  }                                                                                                         \

// Initialise range of possible order values to fast check
// user-supplied order values
unsigned int o_rng = o_size - 1;
  SEXP out;
  switch(CHEAPR_TYPEOF(x)){
  case NILSXP: {
    out = R_NilValue;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    int *p_x = INTEGER(x);
    int fill_value = NA_INTEGER;
    if (fill_size >= 1){
      fill_value = Rf_asInteger(fill);
    }
    out = Rf_protect(Rf_duplicate(x)); ++NP;
    int *p_out = INTEGER(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      // Manually set run rl if order = NULL
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // int trick to calculate inclusive bounded between(x, lo, hi)
      // unsigned int rng = (run_end - 1) - run_start;
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, NP);
        } else {
          oi = j;
        }
        // Costly to use % if we don't need to
        // k = recycle_lag ? p_lag[j % lag_size] : lag1;
        k = recycle_lag ? p_lag[oi % lag_size] : lag1;
        if (k >= 0){
          p_out[oi] = (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        } else {
          p_out[oi] = (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        }
        // lagi = (k != NA_INTEGER) ? j - k : rng + 1;
        // p_out[oi] = unsigned(lagi - run_start) <= rng ? p_x[has_order ? p_o[lagi] - 1 : lagi] : fill_value;
      }
    }
    if (run_end != size){
      Rf_unprotect(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    CHEAPR_LAG_R_NAMES;
    break;
  }
  case CHEAPR_INT64SXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    long long int *p_x = INTEGER64_PTR(x);
    long long int fill_value = NA_INTEGER64;
    if (fill_size >= 1){
      SEXP temp_fill = Rf_protect(coerce_vector(fill, CHEAPR_INT64SXP)); ++NP;
      fill_value = INTEGER64_PTR(temp_fill)[0];
    }
    out = Rf_protect(Rf_duplicate(x));
    ++NP;
    long long int *p_out = INTEGER64_PTR(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, NP);
        } else {
          oi = j;
        }
        k = recycle_lag ? p_lag[oi % lag_size] : lag1;
        if (k >= 0){
          p_out[oi] = (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        } else {
          p_out[oi] = (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        }
      }
    }
    if (run_end != size){
      Rf_unprotect(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    CHEAPR_LAG_R_NAMES;
    break;
  }
  case REALSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    double *p_x = REAL(x);
    double fill_value = NA_REAL;
    if (fill_size >= 1){
      fill_value = Rf_asReal(fill);
    }
    out = Rf_protect(Rf_duplicate(x)); ++NP;
    double *p_out = REAL(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, NP);
        } else {
          oi = j;
        }
        k = recycle_lag ? p_lag[oi % lag_size] : lag1;
        if (k >= 0){
          p_out[oi] = (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        } else {
          p_out[oi] = (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        }
      }
    }
    if (run_end != size){
      Rf_unprotect(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    CHEAPR_LAG_R_NAMES;
    break;
  }
  case STRSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const SEXP *p_x = STRING_PTR_RO(x);
    SEXP fill_value = Rf_protect(fill_size >= 1 ? Rf_asChar(fill) : NA_STRING);
    ++NP;
    out = Rf_protect(Rf_duplicate(x));
    ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, NP);
        } else {
          oi = j;
        }
        k = recycle_lag ? p_lag[oi % lag_size] : lag1;
        if (k >= 0){
          SET_STRING_ELT(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          SET_STRING_ELT(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      Rf_unprotect(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    CHEAPR_LAG_R_NAMES;
    break;
  }
  case CPLXSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    Rcomplex *p_x = COMPLEX(x);
    SEXP fill_sexp = Rf_protect(Rf_allocVector(CPLXSXP, 1));
    ++NP;
    Rcomplex *p_fill = COMPLEX(fill_sexp);
    p_fill[0].i = NA_REAL;
    p_fill[0].r = NA_REAL;
    Rcomplex fill_value = fill_size >= 1 ? Rf_asComplex(fill) : COMPLEX(fill_sexp)[0];
    out = Rf_protect(Rf_duplicate(x));
    ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, NP);
        } else {
          oi = j;
        }
        k = recycle_lag ? p_lag[oi % lag_size] : lag1;
        if (k >= 0){
          SET_COMPLEX_ELT(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          SET_COMPLEX_ELT(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      Rf_unprotect(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    CHEAPR_LAG_R_NAMES;
    break;
  }
  case RAWSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    Rbyte *p_x = RAW(x);
    SEXP raw_sexp = Rf_protect(Rf_coerceVector(fill, RAWSXP));
    Rbyte fill_value = fill_size == 0 ? RAW(Rf_ScalarRaw(0))[0] : RAW(raw_sexp)[0];
    ++NP;
    out = Rf_protect(Rf_duplicate(x));
    ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, NP);
        } else {
          oi = j;
        }
        k = recycle_lag ? p_lag[oi % lag_size] : lag1;
        if (k >= 0){
          SET_RAW_ELT(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          SET_RAW_ELT(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      Rf_unprotect(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    CHEAPR_LAG_R_NAMES;
    break;
  }
  case VECSXP: {
    if (recursive){
    int size = Rf_length(x);
    const SEXP *p_x = VECTOR_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(VECSXP, size));
    ++NP;
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    for (int i = 0; i < size; ++i){
      SET_VECTOR_ELT(out, i, cpp_lag2(p_x[i], lag, order, run_lengths, fill, true));
    }
  } else {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const SEXP *p_x = VECTOR_PTR_RO(x);
    SEXP fill_value = Rf_protect(VECTOR_ELT(Rf_coerceVector(fill_size >= 1 ? fill : R_NilValue, VECSXP), 0));
    ++NP;
    out = Rf_protect(Rf_allocVector(VECSXP, size));
    ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, NP);
        } else {
          oi = j;
        }
        k = recycle_lag ? p_lag[oi % lag_size] : lag1;
        if (k >= 0){
          SET_VECTOR_ELT(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          SET_VECTOR_ELT(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      Rf_unprotect(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    CHEAPR_LAG_R_NAMES;
  }
  break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(NP);
  return out;
}
