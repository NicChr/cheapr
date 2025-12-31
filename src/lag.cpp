#include "cheapr.h"

// Lag vectors and lists recursively, in place, with custom order and run lengths
// Author: Nick Christofides

SEXP lag(SEXP x, R_xlen_t k, SEXP fill, bool set) {
  R_xlen_t size = Rf_xlength(x);
  R_xlen_t fill_size = Rf_xlength(fill);

  int32_t NP = 0;

  if (fill_size > 1){
    Rf_error("fill size must be NULL or length 1");
  }

  if (size == 0){
    return x;
  }

  bool set_and_altrep = set && ALTREP(x);
  if (ALTREP(x)){
    set = true;
  }
  SEXP xvec = SHIELD(altrep_materialise(x)); ++NP;
  if (set_and_altrep){
    Rf_warning("Cannot lag an ALTREP by reference, a copy has been made.\n\tEnsure the result is assigned to an object if used in further calculations\n\te.g. `x <- lag_(x, set = TRUE)`");
  }
  k = k >= 0 ? std::min(size, k) : std::max(-size, k);

  SEXP out = visit_vector(xvec, [&](auto p_x) -> SEXP {

    // Data type of vector elements
    using data_t = std::remove_pointer_t<decltype(p_x)>;

    if constexpr (std::is_same_v<data_t, std::nullptr_t>) {
      if (is_null(xvec)){
        return r_null;
      }
      Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
    } else if constexpr (is_r_ptr_writable_v<data_t>){
      using writable_data_t = std::remove_const_t<data_t>;
      SEXP res = SHIELD(set ? xvec : cpp_semi_copy(xvec)); ++NP;
      auto *p_res = vector_ptr<writable_data_t>(res);
      SHIELD(fill = cast_(get_r_type(xvec), fill, r_null)); ++NP;
      auto fill_value = fill_size > 0 ? get_value<writable_data_t>(fill, 0) : na_value<data_t>();

      if (k >= 0){
        safe_memmove(&p_res[k], &p_x[0], sizeof(data_t) * (size - k));
        r_fill(res, p_res, 0, k, fill_value);
      } else {
        safe_memmove(&p_res[0], &p_x[-k], sizeof(data_t) * (size + k));
        r_fill(res, p_res, size + k, -k, fill_value);
      }
      return res;
    } else {
      SEXP res = SHIELD(set ? xvec : cpp_semi_copy(xvec)); ++NP;
      auto *p_res = vector_ptr<data_t>(res);
      SHIELD(fill = cast_(get_r_type(xvec), fill, r_null)); ++NP;
      auto fill_value = fill_size > 0 ? get_value<data_t>(fill, 0) : na_value<data_t>();

      if (set){
        R_xlen_t tempi;
        // If k = 0 then no lag occurs
        if (std::abs(k) >= 1){
          SEXP lag_temp = SHIELD(new_vector<std::decay_t<data_t>>(std::abs(k))); ++NP;
          std::decay_t<data_t> temp_value;
          auto *p_lag = vector_ptr<data_t>(lag_temp);
          // Positive lags
          if (k >= 0){
            for (R_xlen_t i = 0; i < k; ++i) {
              set_value(lag_temp, i, p_res[i]);
              set_value(res, i, fill_value);
            }
            for (R_xlen_t i = k; i < size; ++i) {
              tempi = ((i - k) % k);
              temp_value = p_lag[tempi];
              set_value(lag_temp, tempi, p_res[i]);
              set_value(res, i, temp_value);
            }
            // Negative lags
          } else {
            for (R_xlen_t i = size - 1; i >= size + k; --i) {
              set_value(lag_temp, size - i - 1, p_res[i]);
              set_value(res, i, fill_value);
            }
            for (R_xlen_t i = size + k - 1; i >= 0; --i) {
              tempi = ( (size - (i - k) - 1) % k);
              temp_value = p_lag[tempi];
              set_value(lag_temp, tempi, p_res[i]);
              set_value(res, i, temp_value);
            }
          }
        }
      } else {
        if (k >= 0){
          for (R_xlen_t i = 0; i < size; ++i) {
            set_value(res, i, i >= k ? p_x[i - k] : fill_value);
          }
        } else {
          for (R_xlen_t i = size - 1; i >= 0; --i) {
            set_value(res, i, (i - size) < k ? p_x[i - k] : fill_value);
          }
        }
      }
      return res;
    }
  });
  YIELD(NP);
  return out;
}

// The reason for having a separate function for doing the recursion is so that
// __restrict__ pointers can be more safely used

[[cpp11::register]]
SEXP cpp_lag(SEXP x, R_xlen_t k, SEXP fill, bool set, bool recursive){
  SEXP out = r_null;
  if (recursive && TYPEOF(x) == VECSXP){
    R_xlen_t size = Rf_xlength(x);
    const SEXP *p_x = list_ptr_ro(x);
    out = SHIELD(new_list(size));
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    for (R_xlen_t i = 0; i < size; ++i){
      SEXP elem = SHIELD(cpp_lag(p_x[i], k, fill, set && !ALTREP(p_x[i]), true));
      SET_VECTOR_ELT(out, i, elem);
      YIELD(1);
    }
    YIELD(1);
  } else {
    out = SHIELD(lag(x, k, fill, set));
    SEXP names = SHIELD(get_old_names(x));
    SEXP lagged_names = SHIELD(lag(names, k, fill, set && !ALTREP(x)));
    set_old_names(out, lagged_names);
    YIELD(3);
  }
  return out;
}

SEXP lag2(SEXP x, SEXP lag, SEXP order, SEXP run_lengths, SEXP fill){
  int o_size = Rf_length(order);
  int rl_size = Rf_length(run_lengths);
  int lag_size = Rf_length(lag);
  int fill_size = Rf_length(fill);
  int32_t NP = 0;
  if (fill_size > 1){
    Rf_error("fill size must be NULL or length 1");
  }
  if (lag_size < 1){
    Rf_error("lag must be a non-zero length integer vector");
  }
  bool has_order = !is_null(order);
  bool has_rl = !is_null(run_lengths);

  // When order is NULL we run through x from left to right (as usual)
  // When run_lengths is NULL we run through x without resetting
  // To do this properly we use dummy vectors so that
  // we can statically assign int pointers
  // While keeping order and run_lengths as NULL (if they are NULL)
  // This is mainly done this way because lag2_ is recursive and
  // hence order/run_lengths should remain NULL throughout each recursion
  // if they are NULL

  SHIELD(lag = vec::coerce_vec(lag, INTSXP)); ++NP;
  SHIELD(order = has_order ? vec::coerce_vec(order, INTSXP) : r_null); ++NP;
  SHIELD(run_lengths = has_rl ? vec::coerce_vec(run_lengths, INTSXP) : r_null); ++NP;
  int foo = 42;
  const int *p_o = &foo;
  const int *p_rl = &foo;
  if (has_order) p_o = integer_ptr(order);
  if (has_rl) p_rl = integer_ptr(run_lengths);
  const int *p_lag = integer_ptr(lag);
  int rl; // Run-length
  int run_start = 0; // Start index of current run
  int run_end = 0; // End index of current run
  int oi; // Indices (specified by order vector) to lag
  int k; // Lag
  // Manually set run rl size to 1 if run_lengths = NULL (checked prior)
  if (!has_rl) rl_size = 1;
  bool recycle_lag = lag_size != 1;
  SEXP out = r_null;
  switch(CHEAPR_TYPEOF(x)){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const int *p_x = integer_ptr(x);
    int fill_value = na::integer;
    if (fill_size >= 1){
      SEXP temp_fill = SHIELD(cast<r_integers_t>(fill, r_null)); ++NP;
      fill_value = integer_ptr(temp_fill)[0];
    }
    out = SHIELD(cpp_semi_copy(x)); ++NP;
    int* RESTRICT p_out = integer_ptr(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      // Manually set run rl if order = NULL
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        YIELD(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        YIELD(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          if (oi < 0 || oi >= size){
            Rf_error("`order` must be an integer vector of unique values between 1 and `length(x)`");
          }
        } else {
          oi = j;
        }
        k = p_lag[recycle_lag ? oi % lag_size : 0];
        if (k >= 0){
          p_out[oi] = (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        } else {
          p_out[oi] = (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        }
      }
    }
    if (run_end != size){
      YIELD(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const int64_t *p_x = integer64_ptr_ro(x);
    int64_t fill_value = na::integer64;
    if (fill_size >= 1){
      SEXP temp_fill = SHIELD(cast<r_integers64_t>(fill, r_null)); ++NP;
      fill_value = integer64_ptr(temp_fill)[0];
    }
    out = SHIELD(cpp_semi_copy(x)); ++NP;
    int64_t* RESTRICT p_out = integer64_ptr(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        YIELD(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        YIELD(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          if (oi < 0 || oi >= size){
            Rf_error("`order` must be an integer vector of unique values between 1 and `length(x)`");
          }
        } else {
          oi = j;
        }
        k = p_lag[recycle_lag ? oi % lag_size : 0];
        if (k >= 0){
          p_out[oi] = (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        } else {
          p_out[oi] = (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        }
      }
    }
    if (run_end != size){
      YIELD(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    break;
  }
  case REALSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const double *p_x = real_ptr(x);
    double fill_value = na::real;
    if (fill_size >= 1){
      SEXP temp_fill = SHIELD(cast<r_doubles_t>(fill, r_null)); ++NP;
      fill_value = real_ptr(temp_fill)[0];
    }
    out = SHIELD(cpp_semi_copy(x)); ++NP;
    double* RESTRICT p_out = real_ptr(out);
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        YIELD(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        YIELD(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          if (oi < 0 || oi >= size){
            Rf_error("`order` must be an integer vector of unique values between 1 and `length(x)`");
          }
        } else {
          oi = j;
        }
        k = p_lag[recycle_lag ? oi % lag_size : 0];
        if (k >= 0){
          p_out[oi] = (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        } else {
          p_out[oi] = (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value;
        }
      }
    }
    if (run_end != size){
      YIELD(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    break;
  }
  case STRSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const r_string_t *p_x = string_ptr_ro(x);
    SEXP temp_fill = SHIELD(cast<r_characters_t>(fill, r_null)); ++NP;
    r_string_t fill_value = fill_size >= 1 ? get_value<r_string_t>(temp_fill, 0) : na::string;
    out = SHIELD(cpp_semi_copy(x)); ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        YIELD(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        YIELD(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          if (oi < 0 || oi >= size){
            Rf_error("`order` must be an integer vector of unique values between 1 and `length(x)`");
          }
        } else {
          oi = j;
        }
        k = p_lag[recycle_lag ? oi % lag_size : 0];
        if (k >= 0){
          set_value<r_string_t>(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          set_value<r_string_t>(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      YIELD(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    break;
  }
  case CPLXSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    r_complex_t *p_x = complex_ptr(x);
    SEXP temp_fill = SHIELD(cast<r_complexes_t>(fill, r_null)); ++NP;
    r_complex_t fill_value = fill_size >= 1 ? get_value<r_complex_t>(temp_fill, 0) : na::complex;
    out = SHIELD(cpp_semi_copy(x)); ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        YIELD(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        YIELD(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          if (oi < 0 || oi >= size){
            Rf_error("`order` must be an integer vector of unique values between 1 and `length(x)`");
          }
        } else {
          oi = j;
        }
        k = p_lag[recycle_lag ? oi % lag_size : 0];
        if (k >= 0){
          set_value<r_complex_t>(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          set_value<r_complex_t>(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      YIELD(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    break;
  }
  case RAWSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    r_byte_t *p_x = raw_ptr(x);
    SEXP raw_sexp = SHIELD(vec::coerce_vec(fill, RAWSXP)); ++NP;
    r_byte_t fill_value = fill_size == 0? r_byte_t{0} : raw_ptr(raw_sexp)[0];
    out = SHIELD(cpp_semi_copy(x)); ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        YIELD(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        YIELD(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          if (oi < 0 || oi >= size){
            Rf_error("`order` must be an integer vector of unique values between 1 and `length(x)`");
          }
        } else {
          oi = j;
        }
        k = p_lag[recycle_lag ? oi % lag_size : 0];
        if (k >= 0){
          set_value<r_byte_t>(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          set_value<r_byte_t>(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      YIELD(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
    }
    break;
  }
  case VECSXP: {
    int size = Rf_length(x);
    if (has_order && (size != o_size)){
      Rf_error("length(order) must equal length(x) (%d)", size);
    }
    const SEXP *p_x = list_ptr_ro(x);
    SEXP fill_value = r_null;
    if (fill_size >= 1){
      SEXP temp_fill = SHIELD(cast<r_list_t>(fill, r_null)); ++NP;
      fill_value = get_value<SEXP>(temp_fill, 0);
    }
    out = SHIELD(new_list(size)); ++NP;
    for (int i = 0; i != rl_size; ++i){
      run_start = run_end; // Start at the end of the previous run
      rl = has_rl ? p_rl[i] : size; // Current run-length
      run_end += rl; // Cumulative run-length

      // If any run-lengths are negative (or NA) we stop
      if (rl < 0){
        YIELD(NP);
        Rf_error("All run lengths must be non-NA and >= 0");
      }

      // If the cumulative run-length exceeds length(x) we stop
      if (run_end > size){
        YIELD(NP);
        Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
      }
      // Loop starting from the end of the previous run-length
      for (int j = run_start; j != run_end; ++j){
        // Check that order value is valid
        if (has_order){
          oi = p_o[j] - 1;
          if (oi < 0 || oi >= size){
            Rf_error("`order` must be an integer vector of unique values between 1 and `length(x)`");
          }
        } else {
          oi = j;
        }
        k = p_lag[recycle_lag ? oi % lag_size : 0];
        if (k >= 0){
          SET_VECTOR_ELT(out, oi, (j - run_start) >= k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        } else {
          SET_VECTOR_ELT(out, oi, (j - run_end) < k ? p_x[has_order ? p_o[j - k] - 1 : j - k] : fill_value);
        }
      }
    }
    if (run_end != size){
      YIELD(NP);
      Rf_error("sum(run_lengths) must be equal to length(x) (%d)", size);
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
SEXP cpp_lag2(SEXP x, SEXP lag, SEXP order, SEXP run_lengths, SEXP fill, bool recursive){
  SEXP out = r_null;
  if (recursive && TYPEOF(x) == VECSXP){
    R_xlen_t size = Rf_xlength(x);
    const SEXP *p_x = list_ptr_ro(x);
    out = SHIELD(new_list(size));
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    for (R_xlen_t i = 0; i < size; ++i){
      SEXP elem = SHIELD(cpp_lag2(p_x[i], lag, order, run_lengths, fill, true));
      SET_VECTOR_ELT(out, i, elem);
      YIELD(1);
    }
    YIELD(1);
  } else {
    SEXP names = SHIELD(get_old_names(x));
    out = SHIELD(lag2(x, lag, order, run_lengths, fill));
    SEXP lagged_names = SHIELD(lag2(names, lag, order, run_lengths, fill));
    set_old_names(out, lagged_names);
    YIELD(3);
  }
  return out;
}
