#include "cheapr_cpp.h"

int int_div(int x, int y){
  return x / y;
}

[[cpp11::register]]
R_xlen_t cpp_vec_length(SEXP x){
  if (Rf_isFrame(x)){
    return cpp_df_nrow(x);
    // Is x a list?
  } else if (Rf_isVectorList(x)){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return cpp_vec_length(VECTOR_ELT(x, 0));
    } else if (Rf_inherits(x, "POSIXlt")){
      const SEXP *p_x = VECTOR_PTR_RO(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
      // return Rf_xlength(VECTOR_ELT(x, 0));
    } else if (Rf_isObject(x)){
      return Rf_asReal(cpp11::package("base")["length"](x));
    } else {
      return Rf_xlength(x);
    }
    // Catch-all
  } else {
    return Rf_xlength(x);
  }
}

int num_cores(){
  SEXP num_cores = Rf_protect(Rf_GetOption1(Rf_installChar(Rf_mkChar("cheapr.cores"))));
  int out = Rf_asInteger(num_cores);
  Rf_unprotect(1);
  return out >= 1 ? out : 1;
}

SEXP xlen_to_r(R_xlen_t x){
  return x > integer_max_ ? Rf_ScalarReal(x) : Rf_ScalarInteger(x);
}

R_xlen_t cpp_df_nrow(SEXP x){
  return Rf_xlength(Rf_getAttrib(x, R_RowNamesSymbol));
}

// Copy names from source to target
void cpp_copy_names(SEXP source, SEXP target){
  SEXP source_nms = Rf_protect(Rf_getAttrib(source, R_NamesSymbol));
  SEXP target_nms = Rf_protect(Rf_duplicate(source_nms));
  Rf_setAttrib(target, R_NamesSymbol, target_nms);
  Rf_unprotect(2);
}

SEXP r_address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return Rf_mkChar(buf);
}

[[cpp11::register]]
SEXP r_copy(SEXP x){
  return Rf_duplicate(x);
}

// Internal-only function
// Sum of squared-differences
// Used in `cheapr_var()`

[[cpp11::register]]
double var_sum_squared_diff(SEXP x, double mu){
  R_xlen_t n = Rf_xlength(x);

  double out = 0;

  // NA values are always ignored here

  if (!cheapr_is_na_dbl(mu)){
    switch (TYPEOF(x)){

    case INTSXP: {
      int *p_x = INTEGER(x);
      // if (std::abs(mu - std::round(mu)) < std::numeric_limits<double>::epsilon()){
      //   long long temp = 0;
      //   long long temp2;
      //   long long llmu = mu;
      //   long long temp1;
      //   for (R_xlen_t i = 0; i < n; ++i){
      //     if (cheapr_is_na_int(p_x[i])) continue;
      //     temp1 = p_x[i];
      //     temp2 = std::pow(temp1 - llmu, 2);
      //     temp += temp2;
      //   }
      //   out = temp;
      // } else {
        for (R_xlen_t i = 0; i < n; ++i){
          if (cheapr_is_na_int(p_x[i])) continue;
          out += std::pow(p_x[i] - mu, 2);
        }
      break;
    }
    default: {
      double *p_x = REAL(x);
      for (R_xlen_t i = 0; i < n; ++i){
        if (cheapr_is_na_dbl(p_x[i])) continue;
        out += std::pow(p_x[i] - mu, 2);
      }
      break;
    }
    }
  } else {
    out = NA_REAL;
  }
  return out;
}

// Credits to R authors
// Re-purposed .bincode
// The main difference is that codes or breaks can be returned efficiently
// Values outside the (right or left) intervals can be included too

#define CHEAPR_BIN_CODES(_IS_NA_, _NA_VAL_)                                                    \
for (R_xlen_t i = 0; i < n; ++i) {                                                             \
  p_out[i] = _NA_VAL_;                                                                         \
  if (!_IS_NA_(p_x[i])) {                                                                      \
    lo = 0;                                                                                    \
    hi = nb1;                                                                                  \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) || \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                     \
      p_out[i] = (left ? hi : lo) + 1;                                                         \
    }                                                                                          \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                         \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                             \
      while (hi - lo >= 2) {                                                                   \
        cutpoint = (hi + lo)/2;                                                                \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                       \
          lo = cutpoint;                                                                       \
        else                                                                                   \
          hi = cutpoint;                                                                       \
      }                                                                                        \
      p_out[i] = lo + 1 + (right && include_oob);                                              \
    }                                                                                          \
  }                                                                                            \
}                                                                                              \

#define CHEAPR_BIN_NCODES(_IS_NA_, _NA_VAL_)                                                                    \
for (R_xlen_t i = 0; i < n; ++i) {                                                                              \
  p_out[i] = _NA_VAL_;                                                                                          \
  if (!_IS_NA_(p_x[i])) {                                                                                       \
    lo = 0;                                                                                                     \
    hi = nb1;                                                                                                   \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) ||                  \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                                      \
      p_out[i] = p_b[(left ? hi : lo)];                                                                         \
    }                                                                                                           \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                                          \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                                              \
      while (hi - lo >= 2) {                                                                                    \
        cutpoint = (hi + lo)/2;                                                                                 \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                                        \
          lo = cutpoint;                                                                                        \
        else                                                                                                    \
          hi = cutpoint;                                                                                        \
      }                                                                                                         \
      p_out[i] = p_b[lo + (right && include_oob)];                                                                                 \
    }                                                                                                           \
  }                                                                                                             \
}                                                                                                               \

[[cpp11::register]]
SEXP cpp_bin(SEXP x, SEXP breaks, bool codes, bool right,
             bool include_lowest,
             bool include_oob){
  int n = Rf_length(x);
  int lo, hi, cutpoint;
  int nb = Rf_length(breaks);
  int nb1 = nb - 1;
  bool left = !right;
  bool include_border = include_lowest;
  switch(TYPEOF(x)){
  case INTSXP: {
    if (codes){
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
    Rf_protect(breaks = Rf_coerceVector(breaks, REALSXP));
    int *p_x = INTEGER(x);
    double *p_b = REAL(breaks);
    int *p_out = INTEGER(out);
    CHEAPR_BIN_CODES(cheapr_is_na_int, NA_INTEGER);
    Rf_unprotect(2);
    return out;
  } else {
    SEXP out = Rf_protect(Rf_duplicate(x));
    Rf_protect(breaks = Rf_coerceVector(breaks, REALSXP));
    int *p_x = INTEGER(x);
    double *p_b = REAL(breaks);
    int *p_out = INTEGER(out);
    CHEAPR_BIN_NCODES(cheapr_is_na_int, NA_INTEGER);
    Rf_unprotect(2);
    return out;
  }
  }
  default: {
    if (codes){
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
    Rf_protect(breaks = Rf_coerceVector(breaks, REALSXP));
    double *p_x = REAL(x);
    double *p_b = REAL(breaks);
    int *p_out = INTEGER(out);
    CHEAPR_BIN_CODES(cheapr_is_na_dbl, NA_INTEGER);
    Rf_unprotect(2);
    return out;
  } else {
    SEXP out = Rf_protect(Rf_duplicate(x));
    Rf_protect(breaks = Rf_coerceVector(breaks, REALSXP));
    double *p_x = REAL(x);
    double *p_b = REAL(breaks);
    double *p_out = REAL(out);
    CHEAPR_BIN_NCODES(cheapr_is_na_dbl, NA_REAL);
    Rf_unprotect(2);
    return out;
  }
  }
  }
}

// Potentially useful for rolling calculations
// Computes the rolling number of true values in a given
// series of consecutive true values

// SEXP cpp_run_id(SEXP x, bool left_to_right){
//   if (!Rf_isLogical(x)){
//     Rf_error("x must be a logical vector");
//   }
//   R_xlen_t count = 0;
//   R_xlen_t n = Rf_xlength(x);
//   SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
//   int *p_out = INTEGER(out);
//   int *p_x = LOGICAL(x);
//   if (left_to_right){
//     for (R_xlen_t i = 0; i < n; ++i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   } else {
//     for (R_xlen_t i = n - 1; i >= 0; --i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }

// Would use data.table as it is very efficient, but would require extra dependency
// Here x must be an integer vector

// SEXP cpp_between(SEXP x, SEXP lower, SEXP upper){
//   int *p_x = INTEGER(x);
//   int n = Rf_length(x);
//   if (Rf_length(lower) != 1){
//     Rf_error("lower must be of length 1");
//   }
//   if (Rf_length(upper) != 1){
//     Rf_error("upper must be of length 1");
//   }
//   Rf_protect(lower = Rf_coerceVector(lower, INTSXP));
//   Rf_protect(upper = Rf_coerceVector(upper, INTSXP));
//   int lo = Rf_asInteger(lower);
//   int hi = Rf_asInteger(upper);
//   if (hi < lo){
//     int hi2 = hi;
//     hi = lo;
//     lo = hi2;
//   }
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
//   int *p_out = LOGICAL(out);
//   bool do_parallel = n >= 100000;
//   int n_cores = do_parallel ? num_cores() : 1;
//   if (lo == NA_INTEGER || hi == NA_INTEGER){
// #pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
//     for (int i = 0; i < n; ++i){
//       p_out[i] = NA_LOGICAL;
//     }
//   } else {
//     unsigned int rng = hi - lo;
// #pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
//     for (int i = 0; i < n; ++i){
//       int xi = p_x[i];
//       p_out[i] = (xi != NA_INTEGER) ? unsigned(xi - lo) <= rng : NA_LOGICAL;
//     }
//   }
//   Rf_unprotect(3);
//   return out;
// }
