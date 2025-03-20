#include "cheapr.h"

// Miscellaneous functions
// Author: Nick Christofides

int int_div(int x, int y){
  return x / y;
}

SEXP xlen_to_r(R_xlen_t x){
  return x > integer_max_ ? Rf_ScalarReal(x) : Rf_ScalarInteger(x);
}

R_xlen_t vec_length(SEXP x){
  if (is_df(x)){
    return cpp_df_nrow(x);
    // Is x a list?
  } else if (Rf_isVectorList(x)){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return vec_length(VECTOR_ELT(x, 0));
    } else if (Rf_inherits(x, "POSIXlt")){
      const SEXP *p_x = VECTOR_PTR_RO(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
      // return Rf_xlength(VECTOR_ELT(x, 0));
    } else if (Rf_isObject(x)){
      return r_length(x);
    } else {
      return Rf_xlength(x);
    }
    // Catch-all
  } else {
    return Rf_xlength(x);
  }
}

[[cpp11::register]]
SEXP cpp_vector_length(SEXP x){
  return xlen_to_r(vec_length(x));
}

int num_cores(){
  SEXP num_cores = Rf_protect(Rf_GetOption1(Rf_installChar(Rf_mkChar("cheapr.cores"))));
  int out = Rf_asInteger(num_cores);
  Rf_unprotect(1);
  return out >= 1 ? out : 1;
}

R_xlen_t cpp_df_nrow(SEXP x){
  return Rf_xlength(Rf_getAttrib(x, R_RowNamesSymbol));
}

SEXP r_address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return Rf_mkChar(buf);
}

[[cpp11::register]]
SEXP cpp_address(SEXP x){
  return Rf_ScalarString(r_address(x));
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

#define CHEAPR_BIN_NCODES(_IS_NA_, _NA_VAL_)                                                                                       \
for (R_xlen_t i = 0; i < n; ++i) {                                                                                                 \
  p_out[i] = _NA_VAL_;                                                                                                             \
  if (!_IS_NA_(p_x[i])) {                                                                                                          \
    lo = 0;                                                                                                                        \
    hi = nb1;                                                                                                                      \
    if ( (include_oob && !include_border && (left ? p_x[i] == p_b[hi] : p_x[i] == p_b[lo])) ||                                     \
         ((include_oob && (left ? p_x[i] > p_b[hi] : p_x[i] < p_b[lo])))){                                                         \
      p_out[i] = p_b[(left ? hi : lo)];                                                                                            \
    }                                                                                                                              \
    else if (!(p_x[i] < p_b[lo] || p_x[i] > p_b[hi] ||                                                                             \
             (p_x[i] == p_b[left ? hi : lo] && !include_border))){                                                                 \
      while (hi - lo >= 2) {                                                                                                       \
        cutpoint = (hi + lo)/2;                                                                                                    \
        if (p_x[i] > p_b[cutpoint] || (left && p_x[i] == p_b[cutpoint]))                                                           \
          lo = cutpoint;                                                                                                           \
        else                                                                                                                       \
          hi = cutpoint;                                                                                                           \
      }                                                                                                                            \
      p_out[i] = p_b[lo + (right && include_oob)];                                                                                 \
    }                                                                                                                              \
  }                                                                                                                                \
}                                                                                                                                  \

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

[[cpp11::register]]
SEXP cpp_if_else(SEXP condition, SEXP yes, SEXP no, SEXP na){
  int NP = 0; // count num protections
  if (TYPEOF(condition) != LGLSXP){
    Rf_error("condition must be a logical vector");
  }
  if (TYPEOF(yes) != TYPEOF(no)){
    Rf_error("`typeof(yes)` must match `typeof(no)`");
  }
  if (TYPEOF(yes) != TYPEOF(na)){
    Rf_error("`typeof(yes)` must match `typeof(na)`");
  }
  R_xlen_t n = Rf_xlength(condition);
  R_xlen_t yes_size = Rf_xlength(yes);
  R_xlen_t no_size = Rf_xlength(no);
  R_xlen_t na_size = Rf_xlength(na);

  if (yes_size != 1 && yes_size != n){
    Rf_error("`length(yes)` must be 1 or `length(condition)`");
  }
  if (no_size != 1 && no_size != n){
    Rf_error("`length(no)` must be 1 or `length(condition)`");
  }
  if (na_size != 1 && na_size != n){
    Rf_error("`length(na)` must be 1 or `length(condition)`");
  }

  bool yes_scalar = yes_size == 1;
  bool no_scalar = no_size == 1;
  bool na_scalar = na_size == 1;

  int *p_x = LOGICAL(condition);
  SEXP out = Rf_protect(Rf_allocVector(TYPEOF(yes), n)); ++NP;

  switch (TYPEOF(yes)){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_out = INTEGER(out);
    int *p_yes = INTEGER(yes);
    int *p_no = INTEGER(no);
    int *p_na = INTEGER(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      p_out[i] = p_yes[yes_scalar ? 0 : i];
      break;
    }
      case false: {
        p_out[i] = p_no[no_scalar ? 0 : i];
        break;
      }
      default: {
        p_out[i] = p_na[na_scalar ? 0 : i];
        break;
      }
      }
      // p_out[i] = p_x[i] == NA_LOGICAL ? na_val : p_x[i] ?
      // p_yes[yes_scalar ? 0 : i] : p_no[no_scalar ? 0 : i];
    }
    break;
  }
  case REALSXP: {
    double *p_out = REAL(out);
    double *p_yes = REAL(yes);
    double *p_no = REAL(no);
    double *p_na = REAL(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      p_out[i] = p_yes[yes_scalar ? 0 : i];
      break;
    }
      case false: {
        p_out[i] = p_no[no_scalar ? 0 : i];
        break;
      }
      default: {
        p_out[i] = p_na[na_scalar ? 0 : i];
        break;
      }
      }
    }
    break;
  }
  case STRSXP: {
    const SEXP *p_yes = STRING_PTR_RO(yes);
    const SEXP *p_no = STRING_PTR_RO(no);
    const SEXP *p_na = STRING_PTR_RO(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      SET_STRING_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case false: {
        SET_STRING_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_STRING_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_yes = COMPLEX(yes);
    Rcomplex *p_no = COMPLEX(no);
    Rcomplex *p_na = COMPLEX(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      SET_COMPLEX_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case false: {
        SET_COMPLEX_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_COMPLEX_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  case RAWSXP: {
    Rbyte *p_yes = RAW(yes);
    Rbyte *p_no = RAW(no);
    Rbyte *p_na = RAW(na);

    for (R_xlen_t i = 0; i < n; ++i){
      switch(p_x[i]){
      case true: {
      SET_RAW_ELT(out, i, p_yes[yes_scalar ? 0 : i]);
      break;
    }
      case false: {
        SET_RAW_ELT(out, i, p_no[no_scalar ? 0 : i]);
        break;
      }
      default: {
        SET_RAW_ELT(out, i, p_na[na_scalar ? 0 : i]);
        break;
      }
      }
    }
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(yes)));
  }
  }
  Rf_unprotect(NP);
  return out;
}

// Counts number of true, false and NAs in a logical vector in one pass

[[cpp11::register]]
SEXP cpp_lgl_count(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  int *p_x = LOGICAL(x);

  R_xlen_t i;
  R_xlen_t ntrue = 0, nfalse = 0;

  if (n_cores > 1){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:ntrue, nfalse)
    for (i = 0; i < n; ++i){
      ntrue += p_x[i] == TRUE;
      nfalse += p_x[i] == FALSE;
    }
  } else {
    OMP_FOR_SIMD
    for (i = 0; i < n; ++i){
      ntrue += p_x[i] == TRUE;
      nfalse += p_x[i] == FALSE;
    }
  }
  R_xlen_t nna = n - ntrue - nfalse;

  SEXP out = Rf_protect(Rf_allocVector(n > integer_max_ ? REALSXP : INTSXP, 3));
  SEXP names = Rf_protect(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("true"));
  SET_STRING_ELT(names, 1, Rf_mkChar("false"));
  SET_STRING_ELT(names, 2, Rf_mkChar("na"));

  if (n > integer_max_){
    SET_REAL_ELT(out, 0, (double) ntrue);
    SET_REAL_ELT(out, 1, (double) nfalse);
    SET_REAL_ELT(out, 2, (double) nna);
  } else {
    SET_INTEGER_ELT(out, 0, (int) ntrue);
    SET_INTEGER_ELT(out, 1, (int) nfalse);
    SET_INTEGER_ELT(out, 2, (int) nna);
  }

  Rf_setAttrib(out, R_NamesSymbol, names);

  Rf_unprotect(2);
  return out;
}

// Copy atomic elements from source to target

[[cpp11::register]]
void cpp_set_copy_elements(SEXP source, SEXP target){
  if (TYPEOF(source) != TYPEOF(target)){
    Rf_error("`typeof(target)` must match `typeof(source)`");
  }
  R_xlen_t n = Rf_xlength(source);

  if (n != Rf_xlength(target)){
    Rf_error("target and source must have the same length");
  }

  switch (TYPEOF(source)){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    int *p_source = INTEGER(source);
    int *p_target = INTEGER(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(int));
    break;
  }
  case REALSXP: {
    double *p_source = REAL(source);
    double *p_target = REAL(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(double));
    break;
  }
  case STRSXP: {
    const SEXP *p_source = STRING_PTR_RO(source);
    for (R_xlen_t i = 0; i < n; ++i){
      SET_STRING_ELT(target, i, p_source[i]);
    }
    break;
  }
  case CPLXSXP: {
    Rcomplex *p_source = COMPLEX(source);
    Rcomplex *p_target = COMPLEX(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(Rcomplex));
    break;
  }
  case RAWSXP: {
    Rbyte *p_source = RAW(source);
    Rbyte *p_target = RAW(target);
    memmove(&p_target[0], &p_source[0], n * sizeof(Rbyte));
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(source)));
  }
  }
}

// Essentially `x & y` but updates x by reference

// SEXP cpp_set_and(SEXP x, SEXP y){
//   R_xlen_t xn = Rf_xlength(x);
//   R_xlen_t yn = Rf_xlength(y);
//
//   R_xlen_t n = xn == 0 || yn == 0 ? 0 : xn;
//
//   R_xlen_t i, yi;
//
//   int *p_x = LOGICAL(x);
//   int *p_y = LOGICAL(y);
//
//
//   for (i = yi = 0; i < n; yi = (++yi == yn) ? 0 : yi, ++i){
//
//     if (p_x[i] != FALSE){
//       if (p_y[yi] == FALSE){
//         p_x[i] = FALSE;
//       } else if ((p_x[i] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
//         p_x[i] = NA_LOGICAL;
//       } else if (p_x[i] == TRUE && p_y[yi] == TRUE){
//         p_x[i] = TRUE;
//       }
//     }
//   }
//   return x;
// }

// Essentially `x | y` but updates x by reference

[[cpp11::register]]
SEXP cpp_set_or(SEXP x, SEXP y){
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);

  R_xlen_t n = xn == 0 || yn == 0 ? 0 : xn;

  R_xlen_t i, yi;

  int *p_x = LOGICAL(x);
  int *p_y = LOGICAL(y);


  for (i = yi = 0; i < n; yi = (++yi == yn) ? 0 : yi, ++i){

    if (p_x[i] != TRUE){
      if (p_y[yi] == TRUE){
        p_x[i] = TRUE;
      } else if ((p_x[i] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
        p_x[i] = NA_LOGICAL;
      } else if (p_x[i] == TRUE || p_y[yi] == TRUE){
        p_x[i] = TRUE;
      }
    }
  }
  return x;
}

// coerceVector() that accounts for int64

SEXP coerce_vector(SEXP source, SEXPTYPE type){
  if (type == CHEAPR_INT64SXP){
    SEXP temp = Rf_protect(Rf_coerceVector(source, REALSXP));
    SEXP out = Rf_protect(cpp_numeric_to_int64(temp));
    Rf_unprotect(2);
    return out;
  } else if (is_int64(source)){
    SEXP temp = Rf_protect(cpp_int64_to_numeric(source));
    SEXP out = Rf_protect(Rf_coerceVector(temp, type));
    Rf_unprotect(2);
    return out;
  } else {
    return Rf_coerceVector(source, type);
  }
}

// Basic growth rate
// i.e the expected percent change per unit of time
// eqn: (new / old ) ^ 1/n
//
// Which is not the same as: new / old
// of which the latter is mistakenly and commonly referred to
// as a 'growth rate'

// o(1) time and space

double growth_rate(double a, double b, double n){
  if (n < 1.0){
    Rf_error("n must be >= 1");
  }
  if (n == R_PosInf){
    Rf_error("n must be finite positive");
  }
  if (n == 1.0){
    return NA_REAL;
  }
  if (a == 0.0 && b == 0.0){
    return 1.0;
  } else {
    return std::pow(b / a, 1.0 / (n - 1.0));
  }
}

#define R_SCALAR_AS_DOUBLE(x, na_val) ((double) (x == na_val ? NA_REAL : x))

[[cpp11::register]]
SEXP cpp_growth_rate(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  if (n == 0){
    return Rf_allocVector(REALSXP, 0);
  }
  if (n == 1){
    return Rf_ScalarReal(NA_REAL);
  }
  double a, b;
  switch(CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    int x_n = INTEGER(x)[n - 1];
    int x_1 = INTEGER(x)[0];
    a = R_SCALAR_AS_DOUBLE(x_1, NA_INTEGER);
    b = R_SCALAR_AS_DOUBLE(x_n, NA_INTEGER);
    break;
  }
  case CHEAPR_INT64SXP: {
    long long int x_n = INTEGER64_PTR(x)[n - 1];
    long long int x_1 = INTEGER64_PTR(x)[0];
    a = R_SCALAR_AS_DOUBLE(x_1, NA_INTEGER64);
    b = R_SCALAR_AS_DOUBLE(x_n, NA_INTEGER64);
    break;
  }
  case REALSXP: {
    b = REAL(x)[n - 1];
    a = REAL(x)[0];
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  return Rf_ScalarReal(growth_rate(a, b, n));
}

SEXP create_df_row_names(int n){
  if (n > 0){
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, 2));
    INTEGER(out)[0] = NA_INTEGER;
    INTEGER(out)[1] = -n;
    Rf_unprotect(1);
    return out;
  } else {
    return Rf_allocVector(INTSXP, 0);
  }
}
