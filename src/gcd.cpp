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
    const int *p_x = INTEGER(x);

    out = SHIELD(vec::new_integer((n == 0) ? 0 : 1));

    if (n > 0){
      int gcd = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = gcd2(gcd, p_x[i], na_rm);
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

    out = SHIELD(new_double((n == 0) ? 0 : 1));

    if (n > 0){
      int64_t gcd = p_x[0];
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = gcd2(gcd, p_x[i], na_rm);
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
    out = SHIELD(new_double((n == 0) ? 0 : 1));
    if (n > 0){
      double gcd = p_x[0];
      double agcd;
      for (R_xlen_t i = 1; i < n; ++i) {
        gcd = gcd2(gcd, p_x[i], na_rm, tol);
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
// SEXP cpp_gcd(SEXP x, double tol, bool na_rm, bool break_early, bool round){
//
//   if (tol < 0 || tol >= 1){
//     Rf_error("tol must be >= 0 and < 1");
//   }
//   R_xlen_t n = Rf_xlength(x);
//
//   SEXP out;
//
//   switch(CHEAPR_TYPEOF(x)){
//   case LGLSXP:
//   case INTSXP: {
//     const int *p_x = INTEGER(x);
//
//     out = SHIELD(vec::new_integer((n == 0) ? 0 : 1));
//
//     if (n > 0){
//       int gcd = p_x[0];
//       if (is_r_na(gcd)){
//         if (na_rm){
//           gcd = 0;
//         } else {
//           set_val(out, 0, gcd);
//           YIELD(1);
//           return out;
//         }
//       }
//       for (R_xlen_t i = 1; i < n; ++i) {
//         if (is_r_na(p_x[i])){
//           if (na_rm){
//             continue;
//           } else {
//             break;
//           }
//         }
//         gcd = gcd2(gcd, p_x[i]);
//
//         if (std::abs(gcd) == 1){
//           break;
//         }
//       }
//       set_val(out, 0, gcd);
//     }
//     break;
//   }
//   case CHEAPR_INT64SXP: {
//     const int64_t *p_x = INTEGER64_PTR_RO(x);
//
//     out = SHIELD(new_double((n == 0) ? 0 : 1));
//
//     if (n > 0){
//       int64_t gcd = p_x[0];
//       if (is_r_na(gcd)){
//         if (na_rm){
//           gcd = 0;
//         } else {
//           set_val(out, 0, as_double(gcd));
//           YIELD(1);
//           return out;
//         }
//       }
//       for (R_xlen_t i = 1; i < n; ++i) {
//         if (is_r_na(p_x[i])){
//           if (na_rm){
//             continue;
//           } else {
//             break;
//           }
//         }
//         gcd = gcd2(gcd, p_x[i]);
//         if (std::abs(gcd) == 1){
//           break;
//         }
//       }
//       set_val(out, 0, as_double(gcd));
//     }
//     break;
//   }
//   default: {
//     const double *p_x = REAL(x);
//     out = SHIELD(new_double((n == 0) ? 0 : 1));
//     if (n > 0){
//       double gcd = p_x[0];
//       if (is_r_na(gcd)){
//         if (na_rm){
//           gcd = 0;
//         } else {
//           set_val(out, 0, gcd);
//           YIELD(1);
//           return out;
//         }
//       }
//       double agcd;
//       for (R_xlen_t i = 1; i < n; ++i) {
//         if (is_r_na(p_x[i])){
//           if (na_rm){
//             continue;
//           } else {
//             break;
//           }
//         }
//         gcd = gcd2(gcd, p_x[i], tol);
//         agcd = std::fabs(gcd);
//         if (break_early && agcd > 0.0 && agcd < (tol + tol)){
//           gcd = tol * static_cast<double>(sign(gcd));
//           break;
//         }
//       }
//       if (round && tol > 0){
//         double factor = std::pow(10, std::ceil(std::fabs(std::log10(tol))) + 1);
//         gcd = std::round(gcd * factor) / factor;
//       }
//       set_val(out, 0, gcd);
//     }
//     break;
//   }
//   }
//   YIELD(1);
//   return out;
// }

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
    return new_integer(0);
  }

    int *p_x = INTEGER(x);

    // Initialise first value as lcm
    int lcm = p_x[0];

    for (R_xlen_t i = 1; i < n; ++i) {
      if (!na_rm && is_r_na(lcm)){
        break;
      }
      auto res = lcm2(lcm, p_x[i], na_rm, 0, overflowed);
      if (overflowed){
        i = n; // Terminate the loop
        double lcm_dbl = r_cast<double>(p_x[0]);
        for (R_xlen_t j = 1; j < n; ++j) {
          if (!na_rm && is_r_na(lcm)){
            break;
          }
          lcm_dbl = lcm2<double>(lcm_dbl, r_cast<double>(p_x[j]), na_rm, 0, overflowed);
        }
        return as_vec(lcm_dbl);
      } else {
        lcm = res;
      }
    }
    return as_vec(lcm);
    // if (is_r_na(lcm) || can_be_int(lcm)){
    //   return as_vec(r_cast<int>(lcm));
    // } else {
    //   return as_vec(r_cast<double>(lcm));
    // }
  }
  case CHEAPR_INT64SXP: {

    if (n == 0){
    return new_double(0);
  }

    int64_t *p_x = INTEGER64_PTR(x);

    // Initialise first value as lcm
    int64_t lcm = p_x[0];

    for (R_xlen_t i = 1; i < n; ++i) {
      if (!na_rm && is_r_na(lcm)){
        break;
      }
      auto res = lcm2(lcm, p_x[i], na_rm, 0, overflowed);
      if (overflowed){
        i = n; // Terminate the loop
        double lcm_dbl = r_cast<double>(p_x[0]);
        for (R_xlen_t j = 1; j < n; ++j) {
          if (!na_rm && is_r_na(lcm)){
            break;
          }
          lcm_dbl = lcm2<double>(lcm_dbl, r_cast<double>(p_x[j]), na_rm, 0, overflowed);
        }
        return as_vec(lcm_dbl);
      } else {
        lcm = res;
      }
    }
    return as_vec(r_cast<double>(lcm));
  }
  default: {

    if (n == 0){
    return new_double(0);
  }

    double *p_x = REAL(x);

    double lcm = p_x[0];
    for (R_xlen_t i = 1; i < n; ++i) {
      if (!na_rm && is_r_na(lcm)){
        lcm = na::numeric;
        break;
      }
      lcm = lcm2(lcm, p_x[i], na_rm, tol);
      if (is_r_inf(lcm)) break;
    }
    return as_vec(lcm);
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
      p_out[i] = gcd2(p_x[xi], p_y[yi], na_rm);
      // if (na_rm){
      //   if (is_r_na(p_x[xi])){
      //     p_out[i] = p_y[yi];
      //   } else if (is_r_na(p_y[yi])){
      //     p_out[i] = p_x[xi];
      //   } else {
      //     p_out[i] = gcd2(p_x[xi], p_y[yi]);
      //   }
      // } else {
      //   if (is_r_na(p_x[xi]) || is_r_na(p_y[yi])){
      //     p_out[i] = na::integer;
      //   } else {
      //     p_out[i] = gcd2(p_x[xi], p_y[yi]);
      //   }
      // }
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = vec::coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(new_double(n)); ++NP;
    double* RESTRICT p_out = REAL(out);
    const double *p_x = REAL(x);
    const double *p_y = REAL(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = gcd2(p_x[xi], p_y[yi], na_rm, tol);
      // if (na_rm){
      //   if (is_r_na(p_x[xi])){
      //     p_out[i] = p_y[yi];
      //   } else if (is_r_na(p_y[yi])){
      //     p_out[i] = p_x[xi];
      //   } else {
      //     p_out[i] = gcd2(p_x[xi], p_y[yi], tol);
      //   }
      // } else {
      //   if (is_r_na(p_x[xi]) || is_r_na(p_y[yi])){
      //     p_out[i] = na::numeric;
      //   } else {
      //     p_out[i] = gcd2(p_x[xi], p_y[yi], tol);
      //   }
      // }
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
    SEXP out = SHIELD(vec::new_integer(n)); ++NP;
    int* RESTRICT p_out = INTEGER(out);
    const int *p_x = INTEGER(x);
    const int *p_y = INTEGER(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = r_cast<int>(lcm2(p_x[xi], p_y[yi], na_rm));
      // dbl_lcm = lcm2(p_x[xi], p_y[yi], 0, na_rm);
      // if (is_r_na(dbl_lcm)|| std::fabs(dbl_lcm) > int_max){
      //   p_out[i] = na::integer;
      // } else {
      //   int_lcm = dbl_lcm;
      //   p_out[i] = int_lcm;
      // }
    }
    YIELD(NP);
    return out;
  }
  default: {
    SHIELD(x = vec::coerce_vec(x, REALSXP)); ++NP;
    SHIELD(y = vec::coerce_vec(y, REALSXP)); ++NP;
    SEXP out = SHIELD(new_double(n)); ++NP;
    double* RESTRICT p_out = REAL(out);
    const double *p_x = REAL(x);
    const double *p_y = REAL(y);
    for (R_xlen_t i = 0, xi = 0, yi = 0; i < n;
    recycle_index(xi, xn),
    recycle_index(yi, yn),
    ++i){
      p_out[i] = r_cast<double>(lcm2(p_x[xi], p_y[yi], na_rm, tol));
    }
    YIELD(NP);
    return out;
  }
  }
}
