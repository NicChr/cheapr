#include "cheapr.h"

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

SEXP r_address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return Rf_mkChar(buf);
}

[[cpp11::register]]
SEXP r_copy(SEXP x){
  return Rf_duplicate(x);
}

// Shallow copy here means that a
// new list is created and
// all attributes are deep-copied

// SEXP list_shallow_copy(SEXP x, bool deep_copy_attrs){
//   if (TYPEOF(x) != VECSXP){
//     Rf_error("Can only shallow copy lists");
//   }
//   R_xlen_t n = Rf_xlength(x);
//   SEXP out = Rf_protect(Rf_allocVector(VECSXP, n));
//   cpp_copy_attributes(x, out, deep_copy_attrs);
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//   for (R_xlen_t i = 0; i < n; ++i){
//     SET_VECTOR_ELT(out, i, p_x[i]);
//   }
//   Rf_unprotect(1);
//   return out;
// }

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

// Safe and very efficient reverse in-place (or with copy)
// ~ O(n/2) time

// matrix structure is preserved with cpp_rev

[[cpp11::register]]
SEXP cpp_rev(SEXP x, bool set){
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t half = int_div(n, 2);
  R_xlen_t n2 = n - 1; // Offset n for 0-indexing
  R_xlen_t k;
  int NP = 0;
  bool rev_names = true;
  bool is_altrep = ALTREP(x);
  if (set && is_altrep){
    Rf_warning("Cannot update an ALTREP by reference, a copy has been made.\n\tEnsure the result is assigned to an object if used in further calculations");
  }
  Rf_protect(x = altrep_materialise(x)); ++NP;

  // altrep will have already been materialised so this should be safe
  // and avoids a second copy
  if (is_altrep){
    set = true;
  }
  SEXP out;
  switch (TYPEOF(x)){
  case NILSXP: {
    out = R_NilValue;
    break;
  }
  case LGLSXP:
  case INTSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    int *p_out = INTEGER(out);
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
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    double *p_out = REAL(out);
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
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    const SEXP *p_out = STRING_PTR_RO(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      SEXP left = Rf_protect(p_out[i]);
      SET_STRING_ELT(out, i, p_out[k]);
      SET_STRING_ELT(out, k, left);
      Rf_unprotect(1);
    }
    break;
  }
  case CPLXSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    Rcomplex *p_out = COMPLEX(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      Rcomplex left = p_out[i];
      SET_COMPLEX_ELT(out, i, p_out[k]);
      SET_COMPLEX_ELT(out, k, left);
    }
    break;
  }
  case RAWSXP: {
    out = Rf_protect(set ? x : Rf_duplicate(x)); ++NP;
    Rbyte *p_out = RAW(out);
    for (R_xlen_t i = 0; i < half; ++i) {
      k = n2 - i;
      Rbyte left = p_out[i];
      SET_RAW_ELT(out, i, p_out[k]);
      SET_RAW_ELT(out, k, left);
    }
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  // We can reverse data frames in-place with the below commented-out code
  // It will not work properly though for data.tables

  // case VECSXP: {
  //   if (recursive){
  //   rev_names = false;
  //   out = Rf_protect(Rf_allocVector(VECSXP, n)); ++NP;
  //   SHALLOW_DUPLICATE_ATTRIB(out, x);
  //   for (R_xlen_t i = 0; i < n; ++i){
  //     SET_VECTOR_ELT(out, i, cpp_rev(VECTOR_ELT(x, i), true, set));
  //   }
  //   break;
  // } else if (!recursive && (!Rf_isObject(x) || Rf_isFrame(x))){
  //   out = Rf_protect(set ? x : list_shallow_copy(x, false)); ++NP;
  //   if (!set){
  //     // SHALLOW_DUPLICATE_ATTRIB(out, x);
  //     cpp_copy_names(x, out, true);
  //   }
  //   const SEXP *p_out = VECTOR_PTR_RO(out);
  //   for (R_xlen_t i = 0; i < half; ++i) {
  //     k = n2 - i;
  //     SEXP left = Rf_protect(p_out[i]);
  //     SET_VECTOR_ELT(out, i, p_out[k]);
  //     SET_VECTOR_ELT(out, k, left);
  //     Rf_unprotect(1);
  //   }
  //     break;
  //   }
  // }
  // default: {
  //   // set rev_names to false because catch-all r_rev() will handle that
  //   rev_names = false;
  //   cpp11::function r_rev = cpp11::package("base")["rev"];
  //   if (set){
  //     Rf_unprotect(NP);
  //     Rf_error("Can't reverse in place here");
  //   }
  //   out = Rf_protect(r_rev(x)); ++NP;
  //   break;
  // }
  }
  // // If x has names, reverse them too
  if (rev_names && !Rf_isNull(Rf_getAttrib(out, R_NamesSymbol))){
    SEXP old_names = Rf_protect(Rf_getAttrib(out, R_NamesSymbol)); ++NP;
    // should be okay to rev in-place here because we already copied the names
    Rf_setAttrib(out, R_NamesSymbol, cpp_rev(old_names, true));
  }
  Rf_unprotect(NP);
  return out;
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

// SEXP cpp_c(SEXP x){
//   if (!Rf_isVectorList(x)){
//     Rf_error("x must be a list of vectors");
//   }
//   int NP = 0;
//   R_xlen_t n = Rf_xlength(x);
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//
//   int vector_type = NILSXP;
//   R_xlen_t out_size = 0;
//
//   for (R_xlen_t i = 0; i < n; ++i){
//     vector_type = std::max(vector_type, TYPEOF(p_x[i]));
//     out_size += Rf_xlength(p_x[i]);
//     if (Rf_isFactor(p_x[i])){;
//       return(cpp11::package("cheapr")["combine_factors"](x));
//     }
//   }
//
//   R_xlen_t k = 0;
//
//   switch(vector_type){
//   case NILSXP: {
//     return R_NilValue;
//   }
//   case LGLSXP:
//   case INTSXP: {
//     SEXP out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
//     int *p_out = INTEGER(out);
//
//     for (R_xlen_t i = 0; i < n; ++i){
//       R_xlen_t m = Rf_xlength(p_x[i]);
//       int *p_temp = INTEGER(TYPEOF(p_x[i]) == vector_type ? p_x[i] : Rf_protect(Rf_coerceVector(p_x[i], vector_type)));
//       if (TYPEOF(p_x[i]) != vector_type) ++NP;
//       for (R_xlen_t j = 0; j < m; ++k, ++j){
//         p_out[k] = p_temp[j];
//       }
//     }
//     Rf_unprotect(NP);
//     return out;
//   }
//   case REALSXP: {
//     SEXP out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
//     double *p_out = REAL(out);
//
//     for (R_xlen_t i = 0; i < n; ++i){
//       R_xlen_t m = Rf_xlength(p_x[i]);
//       double *p_temp = REAL(TYPEOF(p_x[i]) == vector_type ? p_x[i] : Rf_protect(Rf_coerceVector(p_x[i], vector_type)));
//       if (TYPEOF(p_x[i]) != vector_type) ++NP;
//       for (R_xlen_t j = 0; j < m; ++k, ++j){
//         p_out[k] = p_temp[j];
//       }
//     }
//     Rf_unprotect(NP);
//     return out;
//   }
//   case STRSXP: {
//     SEXP out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
//
//     for (R_xlen_t i = 0; i < n; ++i){
//       R_xlen_t m = Rf_xlength(p_x[i]);
//       const SEXP *p_temp = STRING_PTR_RO(TYPEOF(p_x[i]) == vector_type ? p_x[i] : Rf_protect(Rf_coerceVector(p_x[i], vector_type)));
//       if (TYPEOF(p_x[i]) != vector_type) ++NP;
//       for (R_xlen_t j = 0; j < m; ++k, ++j){
//         SET_STRING_ELT(out, k, p_temp[j]);
//       }
//     }
//     Rf_unprotect(NP);
//     return out;
//   }
//   case CPLXSXP: {
//     SEXP out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
//
//     for (R_xlen_t i = 0; i < n; ++i){
//       R_xlen_t m = Rf_xlength(p_x[i]);
//       Rcomplex *p_temp = COMPLEX(TYPEOF(p_x[i]) == vector_type ? p_x[i] : Rf_protect(Rf_coerceVector(p_x[i], vector_type)));
//       if (TYPEOF(p_x[i]) != vector_type) ++NP;
//       for (R_xlen_t j = 0; j < m; ++k, ++j){
//         SET_COMPLEX_ELT(out, k, p_temp[j]);
//       }
//     }
//     Rf_unprotect(NP);
//     return out;
//   }
//   case VECSXP: {
//     SEXP out = Rf_protect(Rf_allocVector(vector_type, out_size)); ++NP;
//
//     for (R_xlen_t i = 0; i < n; ++i){
//       R_xlen_t m = Rf_xlength(p_x[i]);
//       const SEXP *p_temp = VECTOR_PTR_RO(TYPEOF(p_x[i]) == vector_type ? p_x[i] : Rf_protect(Rf_coerceVector(p_x[i], vector_type)));
//       if (TYPEOF(p_x[i]) != vector_type) ++NP;
//       for (R_xlen_t j = 0; j < m; ++k, ++j){
//         SET_VECTOR_ELT(out, k, p_temp[j]);
//       }
//     }
//     Rf_unprotect(NP);
//     return out;
//   }
//
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(vector_type));
//   }
//   }
// }

// bool cpp_all_whole_numbers(SEXP x, double tol, bool na_ignore){
//   R_xlen_t n = Rf_xlength(x);
//   double adiff;
//   bool out = false;
//   switch ( TYPEOF(x) ){
//   case LGLSXP:
//   case INTSXP: {
//     out = true;
//     break;
//   }
//   default: {
//     if (is_int64(x)){
//     out = true;
//     break;
//   }
//     out = true;
//     double *p_x = REAL(x);
//     for (R_xlen_t i = 0; i < n; ++i) {
//       adiff = std::fabs(p_x[i] - std::round(p_x[i]));
//       if (cheapr_is_na_dbl(p_x[i])){
//         if (na_ignore){
//           continue;
//         } else {
//           out = false;
//           break;
//         }
//       } else if (adiff > tol){
//         out = false;
//         break;
//       }
//     }
//     break;
//   }
//   }
//   return out;
// }

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

