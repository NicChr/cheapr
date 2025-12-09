#include "cheapr.h"

// Greatest common divisor and least (smallest) common multiple in R
// Very safe, fast, efficient and works for bit64's vectors
// Use `gcd()` for the GCD across a vector of numbers
// and `gcd2()` for the GCD between 2 numbers (or 2 vectors of numbers)

// Author: Nick Christofides

int gcd2_int(int x, int y, bool na_rm){
  bool has_na = is_r_na(x) || is_r_na(y);
  if (!na_rm && has_na){
    return na::integer;
  }
  if (na_rm && has_na){
    if (is_r_na(x)){
      return y;
    } else {
      return x;
    }
  }
  // GCD(0,0)=0
  if (x == 0 && y == 0){
    return 0;
  }
  // GCD(a,0)=a
  if (x == 0){
    return y;
  }
  // GCD(a,0)=a
  if (y == 0){
    return x;
  }
  int r;
  // Taken from number theory lecture notes
  while(y != 0){
    r = x % y;
    x = y;
    y = r;
  }
  return x;
}

int64_t gcd2_int64(int64_t x, int64_t y, bool na_rm){
  int64_t zero = 0;
  bool has_na = is_r_na(x) || is_r_na(y);
  if (!na_rm && has_na){
    return na::integer64;
  }
  if (na_rm && has_na){
    if (is_r_na(x)){
      return y;
    } else {
      return x;
    }
  }
  // GCD(0,0)=0
  if (x == zero && y == zero){
    return zero;
  }
  // GCD(a,0)=a
  if (x == zero){
    return y;
  }
  // GCD(a,0)=a
  if (y == zero){
    return x;
  }
  int64_t r;
  // Taken from number theory lecture notes
  while(y != zero){
    r = x % y;
    x = y;
    y = r;
  }
  return x;
}

double lcm2(double x, double y, double tol, bool na_rm){
  if (na_rm && ( is_r_na(x) || is_r_na(y))){
    return ( is_r_na(x) ? y : x);
  }
  if (x == 0.0 && y == 0.0){
    return 0.0;
  }
  return ( std::fabs(x) / gcd2(x, y, tol, true) ) * std::fabs(y);
}

int64_t lcm2_int64(int64_t x, int64_t y, bool na_rm){
  int num_nas = is_r_na(x) + is_r_na(y);
  if ( num_nas >= 1 ){
    if (na_rm && num_nas == 1){
      return (is_r_na(x) ? y : x);
    } else {
      return na::integer64;
    }
  }
  if (x == 0 && y == 0){
    return 0;
  }
  // 64-bit integer overflow check
  // Make sure not to divide by zero!

  // res can be an int because the gcd ensures the denom
  // divides x by a whole number

  int64_t res = std::llabs(x) / gcd2_int64(x, y, false);
  if (y != 0 && (std::llabs(res) > (limits::r_int64_max / std::llabs(y)))){
    Rf_error("64-bit integer overflow, please use doubles");
  } else {
    return (res * std::llabs(y));
  }
}

double lcm2_int(int x, int y, bool na_rm){
  int num_nas = is_r_na(x) + is_r_na(y);

  if ( num_nas >= 1 ){
    if (na_rm && num_nas == 1){
      return (is_r_na(x) ? y : x);
    } else {
      return na::numeric;
    }
  }
  if (x == 0 && y == 0){
    return 0.0;
  }
  return ( std::fabs(x) / gcd2_int(x, y, false) ) * std::fabs(y);
}

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
    const int *p_x = INTEGER(x);

    out = SHIELD(vec::new_integer((n == 0) ? 0 : 1));

    if (n > 0){
      int gcd = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = gcd2_int(gcd, p_x[i], na_rm);
        if (is_r_na(gcd)){
          if (!na_rm) break;
        } else if (std::abs(gcd) == 1){
          break;
        }
      }
      set_val(out, 0, gcd);
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    const int64_t *p_x = INTEGER64_PTR_RO(x);

    out = SHIELD(internal::new_vec(REALSXP, (n == 0) ? 0 : 1));

    if (n > 0){
      int64_t gcd = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = gcd2_int64(gcd, p_x[i], na_rm);
        if (is_r_na(gcd)){
          if (!na_rm) break;
        } else if (std::abs(gcd) == 1){
          break;
        }
      }
      set_val(out, 0, as_double(gcd));
    }
    break;
  }
  default: {
    const double *p_x = REAL(x);
    out = SHIELD(internal::new_vec(REALSXP, (n == 0) ? 0 : 1));
    if (n > 0){
      double gcd = p_x[0];
      double agcd;
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = gcd2(gcd, p_x[i], tol, na_rm);
        agcd = std::fabs(gcd);
        if (!na_rm && is_r_na(gcd)){
          break;
        }
        if (break_early && agcd > 0.0 && agcd < (tol + tol)){
          gcd = tol * static_cast<double>(sign(gcd));
          break;
        }
      }
      if (round && tol > 0){
        double factor = std::pow(10, std::ceil(std::fabs(std::log10(tol))) + 1);
        gcd = std::round(gcd * factor) / factor;
      }
      set_val(out, 0, gcd);
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

  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);

    if (n > 0){

      // Initialise first value as lcm
      int64_t lcm = as_int64(p_x[0]);

      for (R_xlen_t i = 1; i < n; ++i) {
        if (!na_rm && is_r_na(lcm)){
          break;
        }
        lcm = lcm2_int64(lcm, as_int64(p_x[i]), na_rm);
      }
      return as_vec(lcm);
    } else {
      return vec::new_integer((n == 0) ? 0 : 1);
    }
  }
  case CHEAPR_INT64SXP: {
    int64_t *p_x = INTEGER64_PTR(x);

    if (n > 0){
      // Initialise first value as lcm
      int64_t lcm = p_x[0];

      for (R_xlen_t i = 1; i < n; ++i) {
        if (!na_rm && is_r_na(lcm)){
          break;
        }
        lcm = lcm2_int64(lcm, p_x[i], na_rm);
      }
      return as_vec(lcm);
    } else {
      return vec::new_integer((n == 0) ? 0 : 1);
    }
  }
  default: {
    double *p_x = REAL(x);

    if (n > 0){
      double lcm = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        if (!na_rm && is_r_na(lcm)){
          lcm = na::numeric;
          break;
        }
        lcm = lcm2(lcm, p_x[i], tol, na_rm);
        if (is_r_inf(lcm)) break;
      }
      return as_vec(lcm);
    } else {
      return internal::new_vec(REALSXP, (n == 0) ? 0 : 1);
    }
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
    SEXP out = SHIELD(vec::new_integer(n)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
    const int *p_x = INTEGER(x);
    const int *p_y = INTEGER(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = gcd2_int(p_x[xi], p_y[yi], na_rm);
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = vec::coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(internal::new_vec(REALSXP, n)); ++NP;
    double* RESTRICT p_out = REAL(out);
    const double *p_x = REAL(x);
    const double *p_y = REAL(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
      ++i){
      p_out[i] = gcd2(p_x[xi], p_y[yi], tol, na_rm);
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
    double dbl_lcm;
    int int_lcm;
    double int_max = limits::r_int_max;
    SHIELD(x = vec::coerce_vec(x, INTSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, INTSXP)); ++NP;
    SEXP out = SHIELD(vec::new_integer(n)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
    const int *p_x = INTEGER(x);
    const int *p_y = INTEGER(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      dbl_lcm = lcm2_int(p_x[xi], p_y[yi], na_rm);
      if (is_r_na(dbl_lcm)|| std::fabs(dbl_lcm) > int_max){
        p_out[i] = na::integer;
      } else {
        int_lcm = dbl_lcm;
        p_out[i] = int_lcm;
      }
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = vec::coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(internal::new_vec(REALSXP, n)); ++NP;
    double* RESTRICT p_out = REAL(out);
    const double *p_x = REAL(x);
    const double *p_y = REAL(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = lcm2(p_x[xi], p_y[yi], tol, na_rm);
    }
    YIELD(NP);
    return out;
  }
  }
}
