#include "cheapr_cpp.h"

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
  int n_protections = 0;
  if (fill_size > 1){
    Rf_error("fill size must be NULL or length 1");
  }
  bool set_and_altrep = set && ALTREP(x);
  SEXP out;
  SEXP xvec = Rf_protect(altrep_materialise(x));
  ++n_protections;
  if (set_and_altrep){
    Rf_warning("Cannot lag an ALTREP by reference, a copy has been made.\n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- lag_(x, set = TRUE)`");
  }
  switch(TYPEOF(xvec)){
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
    ++n_protections;
    int *p_out = INTEGER(out);
    if (set){
      R_xlen_t tempi;
      int tempv;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(INTSXP, std::abs(k)));
        ++n_protections;
        int* __restrict__ p_lag = INTEGER(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            p_lag[i] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            p_lag[size - i - 1] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
        }

      }
    } else {
      int *p_x = INTEGER(xvec);
      if (k >= 0){
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < k; ++i) p_out[i] = fill_value;
        OMP_FOR_SIMD
        for (R_xlen_t i = k; i < size; ++i) p_out[i] = p_x[i - k];
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = size - 1; i >= size + k; --i) p_out[i] = fill_value;
        OMP_FOR_SIMD
        for (R_xlen_t i = size + k - 1; i >= 0; --i) p_out[i] = p_x[i - k];
      }
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
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
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec));
    ++n_protections;
    double *p_out = REAL(out);
    if (set){
      R_xlen_t tempi;
      double tempv;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(REALSXP, std::abs(k)));
        ++n_protections;
        double* __restrict__ p_lag = REAL(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            p_lag[i] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            p_lag[size - i - 1] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
        }

      }
    } else {
      double *p_x = REAL(xvec);
      if (k >= 0){
        OMP_FOR_SIMD
        for (R_xlen_t i = 0; i < k; ++i) p_out[i] = fill_value;
        OMP_FOR_SIMD
        for (R_xlen_t i = k; i < size; ++i) p_out[i] = p_x[i - k];
      } else {
        OMP_FOR_SIMD
        for (R_xlen_t i = size - 1; i >= size + k; --i) p_out[i] = fill_value;
        OMP_FOR_SIMD
        for (R_xlen_t i = size + k - 1; i >= 0; --i) p_out[i] = p_x[i - k];
      }
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case CPLXSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP fill_sexp = Rf_protect(Rf_allocVector(CPLXSXP, 1));
    ++n_protections;
    Rcomplex *p_fill = COMPLEX(fill_sexp);
    p_fill[0].i = NA_REAL;
    p_fill[0].r = NA_REAL;
    Rcomplex fill_value = fill_size >= 1 ? Rf_asComplex(fill) : COMPLEX(fill_sexp)[0];
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec));
    ++n_protections;
    Rcomplex *p_out = COMPLEX(out);
    if (set){
      R_xlen_t tempi;
      Rcomplex tempv;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(CPLXSXP, std::abs(k)));
        ++n_protections;
        Rcomplex* __restrict__ p_lag = COMPLEX(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            p_lag[i].i = p_out[i].i;
            p_lag[i].r = p_out[i].r;
            SET_COMPLEX_ELT(out, i, fill_value);
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            tempv = p_lag[tempi];
            SET_COMPLEX_ELT(lag_temp, tempi, p_out[i]);
            SET_COMPLEX_ELT(out, i, tempv);
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_COMPLEX_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_COMPLEX_ELT(out, i, fill_value);
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            tempv = p_lag[tempi];
            SET_COMPLEX_ELT(lag_temp, tempi, p_out[i]);
            SET_COMPLEX_ELT(out, i, tempv);
          }
        }

      }
    } else {
      Rcomplex *p_x = COMPLEX(xvec);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          SET_COMPLEX_ELT(out, i, i >= k ? p_x[i - k] : fill_value);
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          SET_COMPLEX_ELT(out, i, (i - size) < k ? p_x[i - k] : fill_value);
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(xvec, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(xvec, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case STRSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP fill_char = Rf_protect(fill_size >= 1 ? Rf_asChar(fill) : NA_STRING);
    ++n_protections;
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec));
    ++n_protections;
    const SEXP *p_out = STRING_PTR_RO(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(STRSXP, std::abs(k)));
        ++n_protections;
        SEXP tempv = Rf_protect(Rf_allocVector(STRSXP, 1));
        ++n_protections;
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
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case RAWSXP: {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP raw_sexp = Rf_protect(Rf_coerceVector(fill, RAWSXP));
    Rbyte fill_raw = fill_size == 0 ? RAW(Rf_ScalarRaw(0))[0] : RAW(raw_sexp)[0];
    ++n_protections;
    out = Rf_protect(set ? xvec : Rf_duplicate(xvec));
    ++n_protections;
    Rbyte *p_out = RAW(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(RAWSXP, std::abs(k)));
        ++n_protections;
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
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    const SEXP *p_x = VECTOR_PTR_RO(xvec);
    out = Rf_protect(Rf_allocVector(VECSXP, size));
    ++n_protections;
    SHALLOW_DUPLICATE_ATTRIB(out, xvec);
    for (R_xlen_t i = 0; i < size; ++i){
      SET_VECTOR_ELT(out, i, cpp_lag(p_x[i], k, fill, set, true));
    }
  } else {
    k = k >= 0 ? std::min(size, k) : std::max(-size, k);
    SEXP fill_value = Rf_protect(Rf_coerceVector(fill_size >= 1 ? fill : R_NilValue, VECSXP));
    ++n_protections;
    out = Rf_protect(set ? xvec : Rf_allocVector(VECSXP, size));
    ++n_protections;
    const SEXP *p_out = VECTOR_PTR_RO(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(VECSXP, std::abs(k)));
        ++n_protections;
        SEXP tempv = Rf_protect(Rf_allocVector(VECSXP, 1));
        ++n_protections;
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
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
  }
  break;
  }
  default: {
    Rf_unprotect(n_protections);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(xvec)));
  }
  }
  Rf_unprotect(n_protections);
  return out;
}

[[cpp11::register]]
SEXP cpp_lag2(SEXP x, SEXP lag, SEXP order, SEXP run_lengths, SEXP fill, bool recursive){
  int o_size = Rf_length(order);
  int rl_size = Rf_length(run_lengths);
  int lag_size = Rf_length(lag);
  int fill_size = Rf_length(fill);
  int n_protections = 0;
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
  ++n_protections;
  SEXP dummy_vec2 = Rf_protect(Rf_allocVector(INTSXP, 0));
  ++n_protections;
  Rf_protect(lag = Rf_coerceVector(lag, INTSXP));
  ++n_protections;
  Rf_protect(order = has_order ? Rf_coerceVector(order, INTSXP) : R_NilValue);
  ++n_protections;
  Rf_protect(run_lengths = has_rl ? Rf_coerceVector(run_lengths, INTSXP) : R_NilValue);
  ++n_protections;
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
#define LAG_R_NAMES                                            \
  if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){             \
    SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));\
    ++n_protections;                                           \
    SEXP new_names = Rf_protect(cpp_lag2(old_names, lag, order, run_lengths, R_NilValue, recursive));\
    ++n_protections;                                           \
    Rf_setAttrib(out, R_NamesSymbol, new_names);               \
  }                                                            \

  // Initialise range of possible order values to fast check
  // user-supplied order values
  unsigned int o_rng = o_size - 1;
  SEXP out;
  switch(TYPEOF(x)){
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
    out = Rf_protect(Rf_duplicate(x));
    ++n_protections;
    int *p_out = INTEGER(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      // Manually set run rl if order = NULL
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(n_protections);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(n_protections);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // int trick to calculate inclusive bounded between(x, lo, hi)
      // unsigned int rng = (run_end - 1) - run_start;
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, n_protections);
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
      Rf_unprotect(n_protections);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    LAG_R_NAMES;
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
    out = Rf_protect(Rf_duplicate(x));
    ++n_protections;
    double *p_out = REAL(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(n_protections);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(n_protections);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, n_protections);
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
      Rf_unprotect(n_protections);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    LAG_R_NAMES;
    break;
  }
  case STRSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const SEXP *p_x = STRING_PTR_RO(x);
    SEXP fill_value = Rf_protect(fill_size >= 1 ? Rf_asChar(fill) : NA_STRING);
    ++n_protections;
    out = Rf_protect(Rf_duplicate(x));
    ++n_protections;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(n_protections);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(n_protections);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, n_protections);
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
      Rf_unprotect(n_protections);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    LAG_R_NAMES;
    break;
  }
  case CPLXSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    Rcomplex *p_x = COMPLEX(x);
    SEXP fill_sexp = Rf_protect(Rf_allocVector(CPLXSXP, 1));
    ++n_protections;
    Rcomplex *p_fill = COMPLEX(fill_sexp);
    p_fill[0].i = NA_REAL;
    p_fill[0].r = NA_REAL;
    Rcomplex fill_value = fill_size >= 1 ? Rf_asComplex(fill) : COMPLEX(fill_sexp)[0];
    out = Rf_protect(Rf_duplicate(x));
    ++n_protections;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(n_protections);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(n_protections);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, n_protections);
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
      Rf_unprotect(n_protections);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    LAG_R_NAMES;
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
    ++n_protections;
    out = Rf_protect(Rf_duplicate(x));
    ++n_protections;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(n_protections);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(n_protections);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, n_protections);
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
      Rf_unprotect(n_protections);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    LAG_R_NAMES;
    break;
  }
  case VECSXP: {
    if (recursive){
    int size = Rf_length(x);
    const SEXP *p_x = VECTOR_PTR_RO(x);
    out = Rf_protect(Rf_allocVector(VECSXP, size));
    ++n_protections;
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
    ++n_protections;
    out = Rf_protect(Rf_allocVector(VECSXP, size));
    ++n_protections;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        Rf_unprotect(n_protections);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        Rf_unprotect(n_protections);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          check_order_value(oi, o_rng, n_protections);
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
      Rf_unprotect(n_protections);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    LAG_R_NAMES;
  }
  break;
  }
  default: {
    Rf_unprotect(n_protections);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(n_protections);
  return out;
}
