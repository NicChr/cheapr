#include "cheapr.h"

// Greatest common divisor and least (smallest) common multiple in R
// Very safe, fast, efficient and works for bit64's vectors
// Use `gcd()` for the GCD across a vector of numbers
// and `gcd2()` for the GCD between 2 numbers (or 2 vectors of numbers)

// Author: Nick Christofides

template<typename T> T cpp_sign(T x) {
  return (x > 0) - (x < 0);
}

[[cpp11::register]]
double cpp_gcd2(double x, double y, double tol, bool na_rm){
  double zero = 0.0;
  if (!na_rm && ( !(x == x) || !(y == y) )){
    return NA_REAL;
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
  double r;
  // Taken from number theory lecture notes
  while(std::fabs(y) > tol){
    r = std::fmod(x, y);
    x = y;
    y = r;
  }
  return x;
}

int cpp_gcd2_int(int x, int y, bool na_rm){
  int zero = 0;
  bool has_na = ( x == NA_INTEGER || y == NA_INTEGER );
  if (!na_rm && has_na){
    return NA_INTEGER;
  }
  if (na_rm && has_na){
    if (x == NA_INTEGER){
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
  int r;
  // Taken from number theory lecture notes
  while(y != zero){
    r = x % y;
    x = y;
    y = r;
  }
  return x;
}

int64_t cpp_gcd2_int64(int64_t x, int64_t y, bool na_rm){
  int64_t zero = 0;
  bool has_na = ( x == NA_INTEGER64 || y == NA_INTEGER64 );
  if (!na_rm && has_na){
    return NA_INTEGER64;
  }
  if (na_rm && has_na){
    if (x == NA_INTEGER64){
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

[[cpp11::register]]
double cpp_lcm2(double x, double y, double tol, bool na_rm){
  if (na_rm && ( !(x == x) || !(y == y) )){
    return ( !(x == x) ? y : x);
  }
  if (x == 0.0 && y == 0.0){
    return 0.0;
  }
  return ( std::fabs(x) / cpp_gcd2(x, y, tol, true) ) * std::fabs(y);
}

int64_t cpp_lcm2_int64(int64_t x, int64_t y, bool na_rm){
  int num_nas = (x == NA_INTEGER64) + (y == NA_INTEGER64);
  if ( num_nas >= 1 ){
    if (na_rm && num_nas == 1){
      return (x == NA_INTEGER64 ? y : x);
    } else {
      return NA_INTEGER64;
    }
  }
  if (x == 0 && y == 0){
    return 0;
  }
  // 64-bit integer overflow check
  // Make sure not to divide by zero!

  // res can be an int because the gcd ensures the denom
  // divides x by a whole number

  int64_t res = std::llabs(x) / cpp_gcd2_int64(x, y, false);
  if (y != 0 && (std::llabs(res) > (INTEGER64_MAX / std::llabs(y)))){
    Rf_error("64-bit integer overflow, please use doubles");
  } else {
    return (res * std::llabs(y));
  }
}

double cpp_lcm2_int(int x, int y, bool na_rm){
  int num_nas = (x == NA_INTEGER) + (y == NA_INTEGER);
  if ( num_nas >= 1 ){
    if (na_rm && num_nas == 1){
      return (x == NA_INTEGER ? y : x);
    } else {
      return NA_REAL;
    }
  }
  if (x == 0 && y == 0){
    return 0.0;
  }
  return ( std::fabs(x) / cpp_gcd2_int(x, y, false) ) * std::fabs(y);
}

[[cpp11::register]]
SEXP cpp_gcd(SEXP x, double tol, bool na_rm, bool break_early, bool round){
  if (tol < 0 || tol >= 1){
    Rf_error("tol must be >= 0 and < 1");
  }
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);

  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    const int *p_x = INTEGER(x);
    SEXP out = SHIELD(new_vec(INTSXP, n == 0 ? 0 : 1)); ++NP;
    if (n > 0){
      int gcd = p_x[0];
      int agcd;
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = cpp_gcd2_int(gcd, p_x[i], na_rm);
        if (gcd == NA_INTEGER){
          if (!na_rm) break;
        } else {
          agcd = std::abs(gcd);
          if (agcd > 0 && agcd == 1){
            break;
          }
        }
      }
      INTEGER(out)[0] = gcd;
    }
    YIELD(NP);
    return out;
  }
  case CHEAPR_INT64SXP: {
    const int64_t *p_x = INTEGER64_PTR(x);
    SEXP out = SHIELD(new_vec(REALSXP, n == 0 ? 0 : 1)); ++NP;
    if (n > 0){
      int64_t gcd = p_x[0];
      int64_t agcd;
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = cpp_gcd2_int64(gcd, p_x[i], na_rm);
        if (gcd == NA_INTEGER64){
          if (!na_rm) break;
        } else {
          agcd = std::abs(gcd);
          if (agcd > 0 && agcd == 1){
            break;
          }
        }
      }
      REAL(out)[0] = CHEAPR_INT64_TO_DBL(gcd);
    }
    YIELD(NP);
    return out;
  }
  default: {
    const double *p_x = REAL(x);
    SEXP out = SHIELD(new_vec(REALSXP, n == 0 ? 0 : 1)); ++NP;
    if (n > 0){
      double gcd = p_x[0];
      double agcd;
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = cpp_gcd2(gcd, p_x[i], tol, na_rm);
        agcd = std::fabs(gcd);
        if ((!na_rm && !(gcd == gcd))){
          break;
        }
        if (break_early && agcd > 0.0 && agcd < (tol + tol)){
          gcd = tol * cpp_sign<double>(gcd);
          break;
        }
      }
      if (round && tol > 0){
        double factor = std::pow(10, std::ceil(std::fabs(std::log10(tol))) + 1);
        gcd = std::round(gcd * factor) / factor;
      }
      REAL(out)[0] = gcd;
    }
    YIELD(NP);
    return out;
  }
  }
}

// Lowest common multiple using GCD Euclidean algorithm

[[cpp11::register]]
SEXP cpp_lcm(SEXP x, double tol, bool na_rm){
  if (tol < 0 || tol >= 1){
    Rf_error("tol must be >= 0 and < 1");
  }
  R_xlen_t n = Rf_xlength(x);
  int32_t NP = 0;

  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    int *p_x = INTEGER(x);

    SEXP out;

    if (n > 0){

      // Initialise first value as lcm
      int64_t lcm = CHEAPR_INT_TO_INT64(p_x[0]);

      for (R_xlen_t i = 1; i < n; ++i) {
        if (!na_rm && lcm == NA_INTEGER64){
          break;
        }
        lcm = cpp_lcm2_int64(lcm, CHEAPR_INT_TO_INT64(p_x[i]), na_rm);
      }
      bool is_short = is_na_int64(lcm) || is_integerable<int64_t>(lcm);
      out = SHIELD(new_vec(is_short ? INTSXP : REALSXP, 1)); ++NP;
      if (is_short){
        int temp = CHEAPR_INT64_TO_INT(lcm);
        INTEGER(out)[0] = temp;
      } else {
        double temp = CHEAPR_INT64_TO_DBL(lcm);
        REAL(out)[0] = temp;
      }
    } else {
      out = SHIELD(new_vec(INTSXP, 0)); ++NP;
    }
    YIELD(NP);
    return out;
  }
  case CHEAPR_INT64SXP: {
    int64_t *p_x = INTEGER64_PTR(x);

    SEXP out = SHIELD(new_vec(REALSXP, n == 0 ? 0 : 1)); ++NP;

    if (n > 0){
      // Initialise first value as lcm
      int64_t lcm = p_x[0];

      for (R_xlen_t i = 1; i < n; ++i) {
        if (!na_rm && lcm == NA_INTEGER64){
          break;
        }
        lcm = cpp_lcm2_int64(lcm, p_x[i], na_rm);
      }
      double temp = CHEAPR_INT64_TO_DBL(lcm);
      REAL(out)[0] = temp;
    }
    YIELD(NP);
    return out;
  }
  default: {
    double *p_x = REAL(x);
    SEXP out = SHIELD(new_vec(REALSXP, n == 0 ? 0 : 1)); ++NP;
    if (n > 0){
      double lcm = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        if (!na_rm && !(lcm == lcm)){
          lcm = NA_REAL;
          break;
        }
        lcm = cpp_lcm2(lcm, p_x[i], tol, na_rm);
        if (lcm == R_PosInf || lcm == R_NegInf) break;
      }
      REAL(out)[0] = lcm;
    }
    YIELD(NP);
    return out;
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
    SHIELD(x = coerce_vec(x, INTSXP)); ++NP;
    SHIELD(y = coerce_vec(y, INTSXP)); ++NP;
    SEXP out = SHIELD(new_vec(INTSXP, n)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
    const int *p_x = INTEGER(x);
    const int *p_y = INTEGER(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    xi = (++xi == xn) ? 0 : xi,
      yi = (++yi == yn) ? 0 : yi, ++i){
      p_out[i] = cpp_gcd2_int(p_x[xi], p_y[yi], na_rm);
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(new_vec(REALSXP, n)); ++NP;
    double* RESTRICT p_out = REAL(out);
    const double *p_x = REAL(x);
    const double *p_y = REAL(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    xi = (++xi == xn) ? 0 : xi,
      yi = (++yi == yn) ? 0 : yi, ++i){
      p_out[i] = cpp_gcd2(p_x[xi], p_y[yi], tol, na_rm);
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
    double int_max = INTEGER_MAX;
    SHIELD(x = coerce_vec(x, INTSXP)); ++NP;
    SHIELD(y = coerce_vec(y, INTSXP)); ++NP;
    SEXP out = SHIELD(new_vec(INTSXP, n)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
    const int *p_x = INTEGER(x);
    const int *p_y = INTEGER(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    xi = (++xi == xn) ? 0 : xi,
      yi = (++yi == yn) ? 0 : yi, ++i){
      dbl_lcm = cpp_lcm2_int(p_x[xi], p_y[yi], na_rm);
      if (!(dbl_lcm == dbl_lcm) || std::fabs(dbl_lcm) > int_max){
        p_out[i] = NA_INTEGER;
      } else {
        int_lcm = dbl_lcm;
        p_out[i] = int_lcm;
      }
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(new_vec(REALSXP, n)); ++NP;
    double* RESTRICT p_out = REAL(out);
    const double *p_x = REAL(x);
    const double *p_y = REAL(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    xi = (++xi == xn) ? 0 : xi,
      yi = (++yi == yn) ? 0 : yi, ++i){
      p_out[i] = cpp_lcm2(p_x[xi], p_y[yi], tol, na_rm);
    }
    YIELD(NP);
    return out;
  }
  }
}
