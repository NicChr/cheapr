#include "cheapr_cpp.h"

// Relational operators
// define CHEAPR_OP_SWITCH
// switch(op){
// case 1: {
//   c_op = equals;
//   break;
// }
// case 2: {
//   c_op = gt;
//   break;
// }
// case 3: {
//   c_op = lt;
//   break;
// }
// case 4: {
//   c_op = gte;
//   break;
// }
// case 5: {
//   c_op = lte;
//   break;
// }
// case 6: {
//   c_op = neq;
//   break;
// }
// default: {
//   Rf_error("Supported relational operations: `==`, `>`, `<`, `>=`, `<=`, `!=`");
// }
// }

// #define equals(a, b) ((int) a == b)
// #define gt(a, b) ((int) a > b);
// #define lt(a, b) ((int) a < b)
// #define gte(a, b) ((int) a >= b)
// #define lte(a, b) ((int) a <= b)
// #define neq(a, b) ((int) a != b)

// template <typename T1, typename T2>
// int equals(T1 a, T2 b) { return a == b; }
// template <typename T1, typename T2>
// int gt(T1 a, T2 b) { return a > b; }
// template <typename T1, typename T2>
// int lt(T1 a, T2 b) { return a < b; }
// template <typename T1, typename T2>
// int gte(T1 a, T2 b) { return a >= b; }
// template <typename T1, typename T2>
// int lte(T1 a, T2 b) { return a <= b; }
// template <typename T1, typename T2>
// int neq(T1 a, T2 b) { return a != b; }

bool is_scalar_na(SEXP x){
  if (Rf_xlength(x) != 1){
    Rf_error("x must be a scalar value");
  }
  switch(TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    return cheapr_is_na_int(Rf_asInteger(x));
  }
  case REALSXP: {
    if (is_int64(x)){
    return cheapr_is_na_int64(INTEGER64_PTR(x)[0]);
  } else {
    return cheapr_is_na_dbl(Rf_asReal(x));
  }
  }
  case STRSXP: {
    return cheapr_is_na_str(Rf_asChar(x));
  }
  case CPLXSXP: {
    return cheapr_is_na_cplx(Rf_asComplex(x));
  }
  case RAWSXP: {
    return false;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}

SEXP coerce_vector(SEXP source, SEXP target){
  if (is_int64(target)){
    if (is_int64(source)){
      return source;
    } else {
      return cpp_numeric_to_int64(Rf_coerceVector(source, REALSXP));
    }
  } else {
    return Rf_coerceVector(source, TYPEOF(target));
  }
}

bool implicit_na_coercion(SEXP x, SEXP target){
  SEXP coerced = Rf_protect(coerce_vector(x, target));
  bool out = !is_scalar_na(x) && is_scalar_na(coerced);
  Rf_unprotect(1);
  return out;
}

R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive){
  if (cpp_vec_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t count = 0;
  int NP = 0;
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  SEXP val_is_na = Rf_protect(cpp_is_na(value)); ++NP;
  if (Rf_length(val_is_na) == 1 && LOGICAL(val_is_na)[0]){
    // Can't count NA > NA for example
    // if (op != 1){
    //   Rf_unprotect(NP);
    //   return 0;
    // } else {
    Rf_unprotect(NP);
    return na_count(x, recursive);
    // }
  }
#define CHEAPR_VAL_COUNT(_val_)                                \
  for (R_xlen_t i = 0; i < n; ++i){                            \
    count += (p_x[i] == _val_);                                \
  }                                                            \
                                                               \

  // Alternative that works for other equality operators
  // _IS_NA_ is a arg that accepts a function like cheapr_is_na_int
  // for (R_xlen_t i = 0; i < n; ++i){
  //   count += (c_op(p_x[i], _val_) && !_IS_NA_(p_x[i]));
  // }


  switch ( TYPEOF(x) ){
  case NILSXP: {
    Rf_unprotect(NP);
    return count;
  }
  case LGLSXP:
  case INTSXP: {
    if (implicit_na_coercion(value, x)) break;
    Rf_protect(value = Rf_coerceVector(value, INTSXP)); ++NP;
    int val = Rf_asInteger(value);
    int *p_x = INTEGER(x);
    // int (*c_op)(int, int);
    // CHEAPR_OP_SWITCH;
    if (n_cores > 1){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      CHEAPR_VAL_COUNT(val)
    } else {
#pragma omp for simd
      CHEAPR_VAL_COUNT(val)
    }
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)) break;
    if (is_int64(x)){
      Rf_protect(value = coerce_vector(value, x)); ++NP;
      long long int val = INTEGER64_PTR(value)[0];
      long long int *p_x = INTEGER64_PTR(x);
      // int (*c_op)(long long int, long long int);
      // CHEAPR_OP_SWITCH;
      if (n_cores > 1){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
        CHEAPR_VAL_COUNT(val)
      } else {
#pragma omp for simd
        CHEAPR_VAL_COUNT(val)
      }
    } else {
      Rf_protect(value = Rf_coerceVector(value, REALSXP)); ++NP;
      double val = Rf_asReal(value);
      double *p_x = REAL(x);
      // int (*c_op)(double, double);
      // CHEAPR_OP_SWITCH;
      if (n_cores > 1){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
        CHEAPR_VAL_COUNT(val)
      } else {
#pragma omp for simd
        CHEAPR_VAL_COUNT(val)
      }
    }
    break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)) break;
    Rf_protect(value = Rf_coerceVector(value, STRSXP)); ++NP;
    SEXP val = Rf_protect(Rf_asChar(value)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);
    // int (*c_op)(SEXP, SEXP);
    // CHEAPR_OP_SWITCH;

    if (n_cores > 1){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:count)
      CHEAPR_VAL_COUNT(val);
    } else {
#pragma omp for simd
      CHEAPR_VAL_COUNT(val);
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    const SEXP *p_x = VECTOR_PTR_RO(x);
    for (R_xlen_t i = 0; i < n; ++i){
      count += scalar_count(p_x[i], value, true);
    }
    break;
  }
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(NP);
  return count;
}

[[cpp11::register]]
SEXP cpp_count_val(SEXP x, SEXP value, bool recursive){
  return xlen_to_r(scalar_count(x, value, recursive));
}

// Quick search to return if x is in y vector

// bool x_in_y(int x, SEXP y){
//   int n = Rf_length(y);
//   bool out;
//   if (n == 1){
//     out = (x == Rf_asInteger(y));
//   } else {
//     int *p_y = INTEGER(y);
//     out = false;
//     for (int i = 0; i < n; ++i){
//       if (x == p_y[i]){
//         out = true;
//         break;
//       }
//     }
//   }
//   return out;
// }


[[cpp11::register]]
SEXP cpp_val_replace(SEXP x, SEXP value, SEXP replace, bool recursive){
  int NP = 0;
  R_xlen_t n = Rf_xlength(x);

  if (Rf_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  if (Rf_length(replace) != 1){
    Rf_error("replace must be a vector of length 1");
  }
  bool val_is_na = cpp_any_na(value, true);
  bool any_eq = false;
  bool eq = false;

  SEXP out = Rf_protect(R_NilValue); ++NP;
  switch ( TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SEXP temp = Rf_protect(Rf_allocVector(INTSXP, 0)); ++NP;
    int *p_out = INTEGER(temp);
    Rf_protect(value = Rf_coerceVector(value, INTSXP)); ++NP;
    Rf_protect(replace = Rf_coerceVector(replace, INTSXP)); ++NP;
    int val = Rf_asInteger(value);
    int repl = Rf_asInteger(replace);
    int *p_x = INTEGER(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = Rf_protect(Rf_duplicate(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER(out);
      }
      if (eq) p_out[i] = repl;
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = Rf_protect(x); ++NP;
    }
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    if (is_int64(x)){
    SEXP temp = Rf_protect(Rf_allocVector(REALSXP, 0)); ++NP;
    long long *p_out = INTEGER64_PTR(temp);

    if (!is_int64(value)){
      Rf_protect(value = Rf_coerceVector(value, REALSXP)); ++NP;
      Rf_protect(value = cpp_numeric_to_int64(value)); ++NP;
    }
    if (!is_int64(replace)){
      Rf_protect(replace = Rf_coerceVector(replace, REALSXP)); ++NP;
      Rf_protect(replace = cpp_numeric_to_int64(replace)); ++NP;
    }
    long long val = INTEGER64_PTR(value)[0];
    long long repl = INTEGER64_PTR(replace)[0];
    long long *p_x = INTEGER64_PTR(x);
    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = Rf_protect(Rf_duplicate(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER64_PTR(out);
      }
      if (eq) p_out[i] = repl;
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = Rf_protect(x); ++NP;
    }
  } else {
    SEXP temp = Rf_protect(Rf_allocVector(REALSXP, 0)); ++NP;
    double *p_out = REAL(temp);
    Rf_protect(value = Rf_coerceVector(value, REALSXP)); ++NP;
    Rf_protect(replace = Rf_coerceVector(replace, REALSXP)); ++NP;
    double val = Rf_asReal(value);
    double repl = Rf_asReal(replace);
    double *p_x = REAL(x);
    if (val_is_na){
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] != p_x[i];
        if (!any_eq && eq){
          any_eq = true;
          out = Rf_protect(Rf_duplicate(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = Rf_protect(x); ++NP;
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!any_eq && eq){
          any_eq = true;
          out = Rf_protect(Rf_duplicate(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = Rf_protect(x); ++NP;
      }
    }
  }
  break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    Rf_protect(value = Rf_coerceVector(value, STRSXP)); ++NP;
    Rf_protect(replace = Rf_coerceVector(replace, STRSXP)); ++NP;
    SEXP val = Rf_protect(Rf_asChar(value)); ++NP;
    SEXP repl = Rf_protect(Rf_asChar(replace)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = Rf_protect(Rf_duplicate(x)); ++NP;
      }
      if (eq) SET_STRING_ELT(out, i, repl);
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = Rf_protect(x); ++NP;
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    out = Rf_protect(Rf_allocVector(VECSXP, n)); ++NP;
    for (R_xlen_t i = 0; i < n; ++i){
      // Initialise each element
      SET_VECTOR_ELT(out, i, VECTOR_ELT(x, i));

      // Once we extract the vector it maybe needs protecting??
      SET_VECTOR_ELT(out, i, cpp_val_replace(VECTOR_ELT(out, i), value, replace, true));
    }
    // Copy attributes from x to y
    SEXP attrs = Rf_protect(Rf_coerceVector(ATTRIB(x), VECSXP)); ++NP;
    cpp_set_add_attributes(out, attrs, true);
    break;
  }
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(NP);
  return out;
}
