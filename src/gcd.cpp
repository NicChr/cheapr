#include "cheapr.h"

// Greatest common divisor and least (smallest) common multiple in R
// Very safe, fast, efficient and works for bit64's vectors
// Use `gcd()` for the GCD across a vector of numbers
// and `gcd2()` for the GCD between 2 numbers (or 2 vectors of numbers)

// Author: Nick Christofides

[[cpp11::register]]
SEXP cpp_gcd(SEXP x, double tol, bool na_rm, bool break_early, bool round){

  if (tol < 0 || tol >= 1){
    Rf_error("tol must be >= 0 and < 1");
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out;

  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    const int *p_x = integer_ptr(x);

    out = SHIELD(vec::new_vector<int>((n == 0) ? 0 : 1));

    if (n > 0){
      auto gcd = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = r_gcd(gcd, p_x[i], na_rm);
        if (!na_rm && is_r_na(gcd)){
          break;
        } else if (gcd == 1){
          break;
        }
      }
      set_value(out, 0, r_abs(gcd));
    }
    break;
  }
  case CHEAPR_INT64SXP: {

    const int64_t *p_x = integer64_ptr_ro(x);
    out = SHIELD(new_vector<int64_t>((n == 0) ? 0 : 1));

    if (n > 0){
      auto gcd = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = r_gcd(gcd, p_x[i], na_rm);
        if (!na_rm && is_r_na(gcd)){
          break;
        } else if (gcd == 1){
          break;
        }
      }
      set_value(out, 0, r_abs(gcd));
    }
    break;
  }
  default: {
    const double *p_x = real_ptr(x);
    out = SHIELD(new_vector<double>((n == 0) ? 0 : 1));
    if (n > 0){
      auto gcd = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = r_gcd(gcd, p_x[i], na_rm, tol);
        if (!na_rm && is_r_na(gcd)){
          break;
        }
        if (break_early && gcd > 0.0 && gcd < (tol + tol)){
          gcd = tol * r_sign(gcd);
          break;
        }
      }
      if (round && tol > 0){
        double factor = std::pow(10, std::ceil(std::fabs(std::log10(tol))) + 1);
        gcd = r_round(gcd * factor) / factor;
      }
      set_value(out, 0, r_abs(gcd));
    }
    break;
  }
  }
  YIELD(1);
  return out;
}

// Lowest common multiple using GCD Euclidean algorithm

[[cpp11::register]]
SEXP cpp_lcm(SEXP x, double tol, bool na_rm){
  if (tol < 0 || tol >= 1){
    Rf_error("tol must be >= 0 and < 1");
  }
  R_xlen_t n = Rf_xlength(x);


  bool overflowed = false;

  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {

    if (n == 0){
    return new_vector<int>(0);
  }

    int *p_x = integer_ptr(x);

    // Initialise first value as lcm
    int lcm = p_x[0];

    for (R_xlen_t i = 1; i < n; ++i) {
      if (!na_rm && is_r_na(lcm)){
        break;
      }
      auto res = r_lcm(lcm, p_x[i], na_rm);

      // Only overflowed if result is NA and inputs aren't NA
      overflowed = is_r_na(res) && !is_r_na(lcm) && !is_r_na(p_x[i]);

      if (overflowed){
        i = n; // Terminate the loop
        double lcm_dbl = r_cast<double>(p_x[0]);
        for (R_xlen_t j = 1; j < n; ++j) {
          if (!na_rm && is_r_na(lcm_dbl)){
            break;
          }
          lcm_dbl = r_lcm<double>(lcm_dbl, r_cast<double>(p_x[j]), na_rm, tol);
        }
        return as_vector(lcm_dbl);
      } else {
        lcm = res;
      }
    }
    return as_vector(lcm);
  }
  case CHEAPR_INT64SXP: {

    if (n == 0){
    return new_vector<double>(0);
  }

    int64_t *p_x = integer64_ptr(x);

    // Initialise first value as lcm
    int64_t lcm = p_x[0];

    for (R_xlen_t i = 1; i < n; ++i) {
      if (!na_rm && is_r_na(lcm)){
        break;
      }
      auto res = r_lcm(lcm, p_x[i], na_rm);

      // Only overflowed if result is NA and inputs aren't NA
      overflowed = is_r_na(res) && !is_r_na(lcm) && !is_r_na(p_x[i]);

      if (overflowed){
        i = n; // Terminate the loop
        double lcm_dbl = r_cast<double>(p_x[0]);
        for (R_xlen_t j = 1; j < n; ++j) {
          if (!na_rm && is_r_na(lcm_dbl)){
            break;
          }
          lcm_dbl = r_lcm<double>(lcm_dbl, r_cast<double>(p_x[j]), na_rm, tol);
        }
        return as_vector(lcm_dbl);
      } else {
        lcm = res;
      }
    }
    return as_vector(r_cast<double>(lcm));
  }
  default: {

    if (n == 0){
    return new_vector<double>(0);
  }

    double *p_x = real_ptr(x);

    double lcm = p_x[0];
    for (R_xlen_t i = 1; i < n; ++i) {
      if (!na_rm && is_r_na(lcm)){
        lcm = na::real;
        break;
      }
      lcm = r_lcm(lcm, p_x[i], na_rm, tol);
      if (is_r_inf(lcm)) break;
    }
    return as_vector(lcm);
  }
  }
}

// Vectorised binary gcd

[[cpp11::register]]
SEXP cpp_gcd2_vectorised(SEXP x, SEXP y, double tol, bool na_rm){

  if (tol < 0 || tol >= 1){
    Rf_error("tol must be >= 0 and < 1");
  }

  int32_t NP = 0;
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  R_xlen_t n = std::max(xn, yn);
  if (xn == 0 || yn == 0){
    n = 0;
  }

  if (is_int64(x)){
    SHIELD(x = cpp_int64_to_double(x)); ++NP;
  }
  if (is_int64(y)){
    SHIELD(y = cpp_int64_to_double(y)); ++NP;
  }
  switch(TYPEOF(x)){
  case INTSXP: {
    SHIELD(x = vec::coerce_vec(x, INTSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, INTSXP)); ++NP;
    SEXP out = SHIELD(vec::new_vector<int>(n)); ++NP;
    int* RESTRICT p_out = integer_ptr(out);
    const int *p_x = integer_ptr(x);
    const int *p_y = integer_ptr(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = r_gcd(p_x[xi], p_y[yi], na_rm);
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = vec::coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(new_vector<double>(n)); ++NP;
    double* RESTRICT p_out = real_ptr(out);
    const double *p_x = real_ptr(x);
    const double *p_y = real_ptr(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = r_gcd(p_x[xi], p_y[yi], na_rm, tol);
    }
    YIELD(NP);
    return out;
  }
  }
}

[[cpp11::register]]
SEXP cpp_lcm2_vectorised(SEXP x, SEXP y, double tol, bool na_rm){
  if (tol < 0 || tol >= 1){
    Rf_error("tol must be >= 0 and < 1");
  }
  int32_t NP = 0;
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);
  R_xlen_t n = std::max(xn, yn);
  if (xn == 0 || yn == 0){
    n = 0;
  }

  if (is_int64(x)){
    SHIELD(x = cpp_int64_to_double(x)); ++NP;
  }
  if (is_int64(y)){
    SHIELD(y = cpp_int64_to_double(y)); ++NP;
  }

  switch(TYPEOF(x)){
  case INTSXP: {
    SHIELD(x = vec::coerce_vec(x, INTSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, INTSXP)); ++NP;
    SEXP out = SHIELD(vec::new_vector<int>(n)); ++NP;
    int* RESTRICT p_out = integer_ptr(out);
    const int *p_x = integer_ptr(x);
    const int *p_y = integer_ptr(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = r_cast<int>(r_lcm(p_x[xi], p_y[yi], na_rm));
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = vec::coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(new_vector<double>(n)); ++NP;
    double* RESTRICT p_out = real_ptr(out);
    const double *p_x = real_ptr(x);
    const double *p_y = real_ptr(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = r_cast<double>(r_lcm(p_x[xi], p_y[yi], na_rm, tol));
    }
    YIELD(NP);
    return out;
  }
  }
}
