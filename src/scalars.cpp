#include "cheapr.h"

// Scalar-optimised functions for R
// Author: Nick Christofides

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

void check_atomic(SEXP x){
  if (!Rf_isVectorAtomic(x)){
    Rf_error("'cheapr' scalar functions can only accept atomic vectors");
  }
}

bool implicit_na_coercion(SEXP x, SEXP target){
  SEXP coerced = SHIELD(coerce_vector(x, CHEAPR_TYPEOF(target)));
  bool out = na_count(x, true) != na_count(coerced, true);
  YIELD(1);
  return out;
}

R_xlen_t scalar_count(SEXP x, SEXP value, bool recursive){
  if (vec_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  R_xlen_t n = Rf_xlength(x);
  R_xlen_t count = 0;
  int32_t NP = 0;
  int n_cores = n >= CHEAPR_OMP_THRESHOLD ? num_cores() : 1;

  SEXP val_is_na = SHIELD(cpp_is_na(value)); ++NP;
  if (Rf_length(val_is_na) == 1 && LOGICAL(val_is_na)[0]){
    // Can't count NA > NA for example
    // if (op != 1){
    //   YIELD(NP);
    //   return 0;
    // } else {
    YIELD(NP);
    return na_count(x, recursive);
    // }
  }
#define CHEAPR_VAL_COUNT(VAL)                                \
  for (R_xlen_t i = 0; i < n; ++i){                            \
    count += (p_x[i] == VAL);                                \
  }                                                            \
                                                               \

// Alternative that works for other equality operators
// _IS_NA_ is a arg that accepts a function like is_na_int
// for (R_xlen_t i = 0; i < n; ++i){
//   count += (c_op(p_x[i], _val_) && !_IS_NA_(p_x[i]));
// }


switch ( CHEAPR_TYPEOF(x) ){
case NILSXP: {
  YIELD(NP);
  return count;
}
case LGLSXP:
case INTSXP: {
  if (implicit_na_coercion(value, x)) break;
  SHIELD(value = coerce_vector(value, INTSXP)); ++NP;
  int val = Rf_asInteger(value);
  const int *p_x = INTEGER(x);
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
  SHIELD(value = coerce_vector(value, REALSXP)); ++NP;
  double val = Rf_asReal(value);
  const double *p_x = REAL(x);
  // int (*c_op)(double, double);
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
case CHEAPR_INT64SXP: {
  if (implicit_na_coercion(value, x)) break;
  SHIELD(value = coerce_vector(value, CHEAPR_INT64SXP)); ++NP;
  int64_t val = INTEGER64_PTR(value)[0];
  const int64_t *p_x = INTEGER64_PTR(x);
  // int (*c_op)(int64_t, int64_t);
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
case STRSXP: {
  if (implicit_na_coercion(value, x)) break;
  SHIELD(value = coerce_vector(value, STRSXP)); ++NP;
  SEXP val = SHIELD(Rf_asChar(value)); ++NP;
  const SEXP *p_x = STRING_PTR_RO(x);
  // int (*c_op)(SEXP, SEXP);
  // CHEAPR_OP_SWITCH;
  CHEAPR_VAL_COUNT(val);
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
  YIELD(NP);
  Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
}
}
  YIELD(NP);
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
  int32_t NP = 0;
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

  SEXP out = SHIELD(R_NilValue); ++NP;
  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SEXP temp = SHIELD(new_vec(INTSXP, 0)); ++NP;
    int *p_out = INTEGER(temp);
    SHIELD(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
    SHIELD(replace = coerce_vector(replace, CHEAPR_TYPEOF(x))); ++NP;
    int val = Rf_asInteger(value);
    int repl = Rf_asInteger(replace);
    int *p_x = INTEGER(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER(out);
      }
      if (eq) p_out[i] = repl;
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = SHIELD(x); ++NP;
    }
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SHIELD(value = coerce_vector(value, REALSXP)); ++NP;
    SHIELD(replace = coerce_vector(replace, REALSXP)); ++NP;
    SEXP temp = SHIELD(new_vec(REALSXP, 0)); ++NP;
    double *p_out = REAL(temp);
    double val = Rf_asReal(value);
    double repl = Rf_asReal(replace);
    double *p_x = REAL(x);
    if (val_is_na){
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] != p_x[i];
        if (!any_eq && eq){
          any_eq = true;
          out = SHIELD(cpp_semi_copy(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = SHIELD(x); ++NP;
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!any_eq && eq){
          any_eq = true;
          out = SHIELD(cpp_semi_copy(x)); ++NP;
          // Change where pointer is pointing to
          p_out = REAL(out);
        }
        if (eq) p_out[i] = repl;
      }
      // Make sure to return x if there were no values to replace
      if (!any_eq){
        out = SHIELD(x); ++NP;
      }
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SHIELD(value = coerce_vector(value, CHEAPR_INT64SXP)); ++NP;
    SHIELD(replace = coerce_vector(replace, CHEAPR_INT64SXP)); ++NP;
    SEXP temp = SHIELD(new_vec(REALSXP, 0)); ++NP;
    int64_t *p_out = INTEGER64_PTR(temp);
    int64_t val = INTEGER64_PTR(value)[0];
    int64_t repl = INTEGER64_PTR(replace)[0];
    int64_t *p_x = INTEGER64_PTR(x);
    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
        // Change where pointer is pointing to
        p_out = INTEGER64_PTR(out);
      }
      if (eq) p_out[i] = repl;
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = SHIELD(x); ++NP;
    }
    break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)){
    out = x;
    break;
  }
    SHIELD(value = coerce_vector(value, STRSXP)); ++NP;
    SHIELD(replace = coerce_vector(replace, STRSXP)); ++NP;
    SEXP val = SHIELD(Rf_asChar(value)); ++NP;
    SEXP repl = SHIELD(Rf_asChar(replace)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);

    for (R_xlen_t i = 0; i < n; ++i){
      eq = p_x[i] == val;
      if (!any_eq && eq){
        any_eq = true;
        out = SHIELD(cpp_semi_copy(x)); ++NP;
      }
      if (eq) SET_STRING_ELT(out, i, repl);
    }
    // Make sure to return x if there were no values to replace
    if (!any_eq){
      out = SHIELD(x); ++NP;
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    out = SHIELD(new_vec(VECSXP, n)); ++NP;
    for (R_xlen_t i = 0; i < n; ++i){
      // Initialise each element
      SET_VECTOR_ELT(out, i, VECTOR_ELT(x, i));

      // Once we extract the vector it maybe needs protecting??
      SET_VECTOR_ELT(out, i, cpp_val_replace(VECTOR_ELT(out, i), value, replace, true));
    }
    // Copy attributes from x to y
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    break;
  }
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
SEXP cpp_val_set_replace(SEXP x, SEXP value, SEXP replace, bool recursive){
  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);

  if (Rf_length(value) != 1){
    Rf_error("value must be a vector of length 1");
  }
  if (Rf_length(replace) != 1){
    Rf_error("replace must be a vector of length 1");
  }
  bool val_is_na = cpp_any_na(value, true);

  if (ALTREP(x)){
    Rf_warning("Cannot update an ALTREP by reference, a copy has been made.\n\tEnsure the result is assigned to an object if used in further calculations");
  }
  SHIELD(x = altrep_materialise(x)); ++NP;

  switch ( CHEAPR_TYPEOF(x) ){
  case NILSXP: {
    break;
  }
  case LGLSXP:
  case INTSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
    SHIELD(replace = coerce_vector(replace, CHEAPR_TYPEOF(x))); ++NP;
    int val = Rf_asInteger(value);
    int repl = Rf_asInteger(replace);
    int *p_x = INTEGER(x);

    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      if (p_x[i] == val) p_x[i] = repl;
    }
    break;
  }
  case REALSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = coerce_vector(value, REALSXP)); ++NP;
    SHIELD(replace = coerce_vector(replace, REALSXP)); ++NP;

    double val = Rf_asReal(value);
    double repl = Rf_asReal(replace);
    double *p_x = REAL(x);
    if (val_is_na){
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] != p_x[i]) p_x[i] = repl;
      }
    } else {
      OMP_FOR_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] == val) p_x[i] = repl;
      }
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = coerce_vector(value, CHEAPR_INT64SXP)); ++NP;
    SHIELD(replace = coerce_vector(replace, CHEAPR_INT64SXP)); ++NP;

    int64_t val = INTEGER64_PTR(value)[0];
    int64_t repl = INTEGER64_PTR(replace)[0];
    int64_t *p_x = INTEGER64_PTR(x);
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i){
      if (p_x[i] == val) p_x[i] = repl;
    }
    break;
  }
  case STRSXP: {
    if (implicit_na_coercion(value, x)) break;
    SHIELD(value = coerce_vector(value, STRSXP)); ++NP;
    SHIELD(replace = coerce_vector(replace, STRSXP)); ++NP;
    SEXP val = SHIELD(Rf_asChar(value)); ++NP;
    SEXP repl = SHIELD(Rf_asChar(replace)); ++NP;
    const SEXP *p_x = STRING_PTR_RO(x);

    for (R_xlen_t i = 0; i < n; ++i){
      for (R_xlen_t i = 0; i < n; ++i){
        if (p_x[i] == val) SET_STRING_ELT(x, i, repl);
      }
    }
    break;
  }
  case VECSXP: {
    if (recursive){
    for (R_xlen_t i = 0; i < n; ++i){
      SET_VECTOR_ELT(x, i, cpp_val_replace(VECTOR_ELT(x, i), value, replace, true));
    }
    break;
  }
  }
  default: {
    YIELD(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  YIELD(NP);
  return x;
}

// At the moment this doesn't coerce what to type of x
// or handle long vectors

[[cpp11::register]]
SEXP cpp_loc_set_replace(SEXP x, SEXP where, SEXP what){
  if (TYPEOF(x) != TYPEOF(what)){
    Rf_error("`typeof(x)` must match `typeof(what)`");
  }
  int *p_where = INTEGER(where);

  if (ALTREP(x)){
    Rf_warning("Cannot update an ALTREP by reference, a copy has been made.\n\tEnsure the result is assigned to an object if used in further calculations");
  }
  SHIELD(x = altrep_materialise(x));

  int64_t xn = Rf_xlength(x);
  int where_size = Rf_length(where);
  int what_size = Rf_length(what);
  if (what_size != 1 && where_size != what_size){
    YIELD(1);
    Rf_error("`what` must be either length 1 or `length(where)`");
  }
  int64_t xi;


#define CHEAPR_REPLACE                                                                                      \
  if (what_size == 1){                                                                                      \
    for (int i = 0; i < where_size; ++i){                                                                   \
      xi = p_where[i];                                                                                      \
      if (xi <= 0 || xi > xn){                                                                   \
        YIELD(1);                                                                                           \
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");                    \
      }                                                                                                     \
      p_x[xi - 1] = p_what[0];                                                                              \
    }                                                                                                       \
  } else {                                                                                                  \
    for (int i = 0; i < where_size; ++i){                                                                   \
      xi = p_where[i];                                                                                      \
      if (xi <= 0 || xi > xn){                                                                   \
        YIELD(1);                                                                                           \
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");                    \
      }                                                                                                     \
      p_x[xi - 1] = p_what[i];                                                                              \
    }                                                                                                       \
  }                                                                                                         \


switch (TYPEOF(x)){
case NILSXP: {
  break;
}
case LGLSXP:
case INTSXP: {
  int *p_x = INTEGER(x);
  int *p_what = INTEGER(what);
  CHEAPR_REPLACE
  break;
}
case REALSXP: {
  double *p_x = REAL(x);
  double *p_what = REAL(what);
  CHEAPR_REPLACE
  break;
}
case STRSXP: {
  const SEXP *p_what = STRING_PTR_RO(what);

  if (what_size == 1){
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_STRING_ELT(x, xi - 1, p_what[0]);
    }
  } else {
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_STRING_ELT(x, xi - 1, p_what[i]);
    }
  }
  break;
}
case CPLXSXP: {
  Rcomplex *p_what = COMPLEX(what);

  if (what_size == 1){
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_COMPLEX_ELT(x, xi - 1, p_what[0]);
    }
  } else {
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_COMPLEX_ELT(x, xi - 1, p_what[i]);
    }
  }
  break;
}
case RAWSXP: {
  Rbyte *p_what = RAW(what);

  if (what_size == 1){
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_RAW_ELT(x, xi - 1, p_what[0]);
    }
  } else {
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_RAW_ELT(x, xi - 1, p_what[i]);
    }
  }
  break;
}
case VECSXP: {
  const SEXP *p_what = VECTOR_PTR_RO(what);

  if (what_size == 1){
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_VECTOR_ELT(x, xi - 1, p_what[0]);
    }
  } else {
    for (int i = 0; i < where_size; ++i){
      xi = p_where[i];
      if (xi <= 0 || xi > xn){
        YIELD(1);
        Rf_error("where must be an integer vector of values between 1 and `length(x)`");
      }
      SET_VECTOR_ELT(x, xi - 1, p_what[i]);
    }
  }
  break;
}
default: {
  YIELD(1);
  Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
}
}
  YIELD(1);
  return x;
}

// Remove elements from a vector very efficiently

[[cpp11::register]]
SEXP cpp_val_remove(SEXP x, SEXP value){
  check_atomic(x);
  int32_t NP = 0;
  R_xlen_t n_vals = scalar_count(x, value, true);
  if (n_vals == 0){
    return x;
  } else if (n_vals == Rf_xlength(x)){
    SEXP out = SHIELD(new_vec(TYPEOF(x), 0)); ++NP;

    cpp_set_add_attributes(out, ATTRIB(x), false);
    YIELD(NP);
    return out;
    return cheapr_sset(x, 0);
  } else {
    R_xlen_t n = Rf_xlength(x);
    R_xlen_t n_keep = n - n_vals;
    bool eq;
    R_xlen_t k = 0;

    SEXP out;
    switch ( CHEAPR_TYPEOF(x) ){
    case NILSXP: {
      out = SHIELD(R_NilValue); ++NP;
      break;
    }
    case LGLSXP:
    case INTSXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
      int val = Rf_asInteger(value);
      int *p_x = INTEGER(x);
      int *p_out = INTEGER(out);

      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!eq){
          p_out[k++] = p_x[i];
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case REALSXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
      double val = Rf_asReal(value);
      double *p_x = REAL(x);
      double *p_out = REAL(out);

      if (cpp_any_na(value, true)){
        for (R_xlen_t i = 0; i < n; ++i){
          eq = is_na_dbl(p_x[i]);
          if (!eq){
            p_out[k++] = p_x[i];
          }
        }
      } else {
        for (R_xlen_t i = 0; i < n; ++i){
          eq = p_x[i] == val;
          if (!eq){
            p_out[k++] = p_x[i];
          }
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case CHEAPR_INT64SXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
      int64_t val = INTEGER64_PTR(value)[0];
      int64_t *p_x = INTEGER64_PTR(x);
      int64_t *p_out = INTEGER64_PTR(out);

      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!eq){
          p_out[k++] = p_x[i];
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    case STRSXP: {
      if (implicit_na_coercion(value, x)){
      out = x;
      break;
    }
      out = SHIELD(new_vec(TYPEOF(x), n_keep)); ++NP;
      SHIELD(value = coerce_vector(value, CHEAPR_TYPEOF(x))); ++NP;
      SEXP val = SHIELD(Rf_asChar(value)); ++NP;
      const SEXP *p_x = STRING_PTR_RO(x);

      for (R_xlen_t i = 0; i < n; ++i){
        eq = p_x[i] == val;
        if (!eq){
          SET_STRING_ELT(out, k++, p_x[i]);
        }
      }
      cpp_set_add_attributes(out, ATTRIB(x), false);
      break;
    }
    default: {
      SEXP sexp_n_vals = SHIELD(Rf_ScalarReal(n_vals)); ++NP;
      SEXP val_locs = SHIELD(cpp_val_find(x, value, true, sexp_n_vals)); ++NP;
      out = SHIELD(cheapr_sset(x, val_locs)); ++NP;
      break;
    }
    }
    YIELD(NP);
    return out;
  }
}
