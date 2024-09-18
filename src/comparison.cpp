#include "cheapr_cpp.h"
#include <string>

// Initialise logical with all TRUE

// SEXP cpp_initialise_lgl(R_xlen_t n){
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
//   int *p_out = LOGICAL(out);
//   OMP_FOR_SIMD
//   for (R_xlen_t i = 0; i < n; ++i){
//     p_out[i] = TRUE;
//   }
//   Rf_unprotect(1);
//   return out;
// }
//
// void cpp_replace_lgl(SEXP lgl, SEXP x, SEXP y){
//   R_xlen_t n1 = Rf_xlength(x);
//   R_xlen_t n2  = Rf_xlength(y);
//   int *p_lgl = LOGICAL(lgl);
//   if (TYPEOF(x) != TYPEOF(y)){
//     Rf_error("x and y must be of the same type");
//   }
//   if (n1 != n2){
//     Rf_error("x and y must be of the same length");
//   }
//   switch ( TYPEOF(x) ){
//   case LGLSXP:
//   case INTSXP: {
//     int a, b;
//     int *p_x = INTEGER(x);
//     int *p_y = INTEGER(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       a = p_x[i];
//       b = p_y[i];
//       p_lgl[i] =
//         (p_lgl[i] == NA_LOGICAL || a == NA_INTEGER || b == NA_INTEGER) ?
//         NA_LOGICAL :
//         p_lgl[i] && a == b;
//     }
//     break;
//   }
//   case REALSXP: {
//     double a, b;
//     double *p_x = REAL(x);
//     double *p_y = REAL(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       a = p_x[i];
//       b = p_y[i];
//       p_lgl[i] =
//         (p_lgl[i] == NA_LOGICAL || a != a || b != b) ?
//         NA_LOGICAL :
//         p_lgl[i] && a == b;
//     }
//     break;
//   }
//   case STRSXP: {
//     const SEXP *p_x = STRING_PTR_RO(x);
//     const SEXP *p_y = STRING_PTR_RO(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       SEXP a = p_x[i];
//       SEXP b = p_y[i];
//       p_lgl[i] =
//         (p_lgl[i] == NA_LOGICAL || a == NA_STRING || b == NA_STRING) ?
//         NA_LOGICAL :
//       p_lgl[i] && a == b;
//     }
//     break;
//   }
//   case CPLXSXP: {
//     Rcomplex a, b;
//     Rcomplex *p_x = COMPLEX(x);
//     Rcomplex *p_y = COMPLEX(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       a = p_x[i];
//       b = p_y[i];
//       p_lgl[i] =
//         (p_lgl[i] == NA_LOGICAL || a.r != a.r|| a.i != a.i || b.r != b.r || b.i != b.i) ?
//         NA_LOGICAL :
//         p_lgl[i] && (a.r == b.r) && (a.i == b.i);
//     }
//     break;
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
// }
//
// SEXP cpp_data_frames_equal(SEXP x, SEXP y){
//   R_xlen_t n1 = Rf_xlength(x);
//   R_xlen_t n2  = Rf_xlength(y);
//   R_xlen_t nrow = cpp_df_nrow(x);
//   if (n1 != n2){
//     Rf_error("x and y must be of the same length");
//   }
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//   const SEXP *p_y = VECTOR_PTR_RO(y);
//
//   SEXP out = Rf_protect(cpp_initialise_lgl(nrow));
//   for (R_xlen_t i = 0; i < n1; ++i){
//     cpp_replace_lgl(out, p_x[i], p_y[i]);
//   }
//   Rf_unprotect(1);
//   return out;
// }
//
// void cpp_replace_lgl2(int* __restrict__ p_lgl, SEXP x, SEXP y){
//   R_xlen_t n1 = Rf_xlength(x);
//   R_xlen_t n2  = Rf_xlength(y);
//   if (TYPEOF(x) != TYPEOF(y)){
//     Rf_error("x and y must be of the same type");
//   }
//   if (n1 != n2){
//     Rf_error("x and y must be of the same length");
//   }
//   switch ( TYPEOF(x) ){
//   case LGLSXP:
//   case INTSXP: {
//     int a, b;
//     int *p_x = INTEGER(x);
//     int *p_y = INTEGER(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       if (p_lgl[i] != FALSE){
//       a = p_x[i];
//       b = p_y[i];
//
//       p_lgl[i] =
//         (a == NA_INTEGER || b == NA_INTEGER || (p_lgl[i] == NA_LOGICAL && a == b)) ?
//         NA_LOGICAL : a == b;
//       }
//     }
//     break;
//   }
//   case REALSXP: {
//     double a, b;
//     double *p_x = REAL(x);
//     double *p_y = REAL(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       if (p_lgl[i] != FALSE){
//       a = p_x[i];
//       b = p_y[i];
//       p_lgl[i] =
//         (a != a || b != b || (p_lgl[i] == NA_LOGICAL && a == b)) ? NA_LOGICAL : a == b;
//       }
//     }
//     break;
//   }
//   case STRSXP: {
//     const SEXP *p_x = STRING_PTR_RO(x);
//     const SEXP *p_y = STRING_PTR_RO(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       if (p_lgl[i] != FALSE){
//       SEXP a = p_x[i];
//       SEXP b = p_y[i];
//       p_lgl[i] =
//         (a == NA_STRING || b == NA_STRING || (p_lgl[i] == NA_LOGICAL && a == b)) ? NA_LOGICAL : a == b;
//       }
//     }
//     break;
//   }
//   case CPLXSXP: {
//     Rcomplex a, b;
//     Rcomplex *p_x = COMPLEX(x);
//     Rcomplex *p_y = COMPLEX(y);
//     // OMP_FOR_SIMD
//     for (R_xlen_t i = 0; i < n1; ++i){
//       if (p_lgl[i] != FALSE){
//       a = p_x[i];
//       b = p_y[i];
//       p_lgl[i] = (
//         ( (a.r != a.r) || (a.i != a.i) || (b.r != b.r) || (b.i != b.i) ) ||
//           (p_lgl[i] == NA_LOGICAL && a.r == b.r && a.i == b.i) ?
//           NA_LOGICAL :
//         (a.r == b.r) && (a.i == b.i)
//       );
//       }
//     }
//     break;
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
// }
//
// SEXP cpp_data_frames_equal2(SEXP x, SEXP y){
//   R_xlen_t n1 = Rf_xlength(x);
//   R_xlen_t n2  = Rf_xlength(y);
//   R_xlen_t nrow = cpp_df_nrow(x);
//   if (n1 != n2){
//     Rf_error("x and y must be of the same length");
//   }
//   const SEXP *p_x = VECTOR_PTR_RO(x);
//   const SEXP *p_y = VECTOR_PTR_RO(y);
//
//   SEXP out = Rf_protect(cpp_initialise_lgl(nrow));
//   int* __restrict__ p_out = LOGICAL(out);
//   for (R_xlen_t i = 0; i < n1; ++i){
//     cpp_replace_lgl2(p_out, p_x[i], p_y[i]);
//   }
//   Rf_unprotect(1);
//   return out;
// }
//
// bool cpp_all_equal(SEXP x, SEXP y){
//   R_xlen_t n1 = Rf_xlength(x);
//   R_xlen_t n2  = Rf_xlength(y);
//   if (TYPEOF(x) != TYPEOF(y)){
//     Rf_error("x and y must be of the same type");
//   }
//   if (n1 != n2){
//     Rf_error("x and y must be of the same length");
//   }
//   // int out = 1;
//   bool out = true;
//   // bool either_na;
//   // R_xlen_t na_count = 0;
//   switch ( TYPEOF(x) ){
//   case LGLSXP:
//   case INTSXP: {
//     int a, b;
//     int *p_x = INTEGER(x);
//     int *p_y = INTEGER(y);
//     for (R_xlen_t i = 0; i < n1; ++i){
//       a = p_x[i];
//       b = p_y[i];
//       // either_na = (a == NA_INTEGER || b == NA_INTEGER);
//       // na_count += either_na;
//       // if (a != b && !either_na){
//       if (a != b){
//         out = false;
//         break;
//       }
//     }
//     break;
//   }
//   case REALSXP: {
//     double a, b;
//     // bool both_na;
//     double *p_x = REAL(x);
//     double *p_y = REAL(y);
//     for (R_xlen_t i = 0; i < n1; ++i){
//       a = p_x[i];
//       b = p_y[i];
//       // either_na = (a != a || b != b);
//       // na_count += either_na;
//       // if (a != b && !either_na){
//       if (a != b && (!(a != a && b != b))){
//         out = false;
//         break;
//       }
//     }
//     break;
//   }
//   case STRSXP: {
//     const SEXP *p_x = STRING_PTR_RO(x);
//     const SEXP *p_y = STRING_PTR_RO(y);
//     for (R_xlen_t i = 0; i < n1; ++i){
//       // either_na = (p_x[i] == NA_STRING || p_y[i] == NA_STRING);
//       // na_count += either_na;
//       // if (p_x[i] != p_y[i] && !either_na){
//       if (p_x[i] != p_y[i]){
//         out = false;
//         break;
//       }
//     }
//     break;
//   }
//   case CPLXSXP: {
//     Rcomplex a, b;
//     Rcomplex *p_x = COMPLEX(x);
//     Rcomplex *p_y = COMPLEX(y);
//     for (R_xlen_t i = 0; i < n1; ++i){
//       a = p_x[i];
//       b = p_y[i];
//       // either_na = cheapr_is_na_cplx(a) || cheapr_is_na_cplx(b);
//       // na_count += either_na;
//       // if (( (a.r != a.r) || (a.i != a.i) ) && !either_na){
//       if (( (a.r != a.r) || (a.i != a.i) )){
//         out = false;
//         break;
//       }
//     }
//     break;
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
//   // if (out == 1 && na_count > 0){
//   //   out = NA_INTEGER;
//   // }
//   // return Rf_ScalarLogical(out);
//   return out;
// }
//


// Relational operators
template <typename T1, typename T2>
int equals(T1 a, T2 b) { return a == b; }
template <typename T1, typename T2>
int gt(T1 a, T2 b) { return a > b; }
template <typename T1, typename T2>
int lt(T1 a, T2 b) { return a < b; }
template <typename T1, typename T2>
int gte(T1 a, T2 b) { return a >= b; }
template <typename T1, typename T2>
int lte(T1 a, T2 b) { return a <= b; }
template <typename T1, typename T2>
int neq(T1 a, T2 b) { return a != b; }


// #define cheapr_char_eq(a, b)((bool)(std::strcmp(a, b) == 0))

#define OP_SWITCH                                                                                 \
switch(op){                                                                                       \
case 1: {                                                                                         \
  c_op = equals;                                                                                  \
  break;                                                                                          \
}                                                                                                 \
case 2: {                                                                                         \
  c_op = gt;                                                                                      \
  break;                                                                                          \
}                                                                                                 \
case 3: {                                                                                         \
  c_op = lt;                                                                                      \
  break;                                                                                          \
}                                                                                                 \
case 4: {                                                                                         \
  c_op = gte;                                                                                     \
  break;                                                                                          \
}                                                                                                 \
case 5: {                                                                                         \
  c_op = lte;                                                                                     \
  break;                                                                                          \
}                                                                                                 \
case 6: {                                                                                         \
  c_op = neq;                                                                                     \
  break;                                                                                          \
}                                                                                                 \
default: {                                                                                        \
  Rf_unprotect(NP);                                                                               \
  Rf_error("Supported relational operations: `==`, `>`, `<`, `>=`, `<=`, `!=`");                  \
}                                                                                                 \
}                                                                                                 \

#define CHEAPR_CHAR_COMPARISON_LOOP(a, b, ISNA1, ISNA2)        \
for (i = xi = yi = 0; i < n; xi = (++xi == xn) ? 0 : xi,       \
     yi = (++yi == yn) ? 0 : yi, ++i){                         \
  p_out[i] = (ISNA1(p_x[xi]) || ISNA2(p_y[yi])) ?              \
  NA_LOGICAL : c_op(a, b);                                     \
}                                                              \

// Below is generally faster but gives wrong results sometimes?

#define CHEAPR_CHAR_COMPARISON_LOOP2(a, b, ISNA1, ISNA2)       \
for (i = xi = yi = 0; i < n; xi = (++xi == xn) ? 0 : xi,       \
     yi = (++yi == yn) ? 0 : yi, ++i){                         \
  p_out[i] = (ISNA1(p_x[xi]) || ISNA2(p_y[yi])) ?              \
  NA_LOGICAL : c_op(std::strcmp(a, b), 0);                     \
}                                                              \

// Fun idea, may explore later
// Adding the number of true values to a logical vector
// as you're creating it can speed up cheapr::which_()
// especially when .n.true == 0

// SEXP cpp_compare(SEXP x, SEXP y, int op){
//   int NP = 0;
//   R_xlen_t xn = Rf_xlength(x);
//   R_xlen_t yn = Rf_xlength(y);
//
//   R_xlen_t n = xn == 0 || yn == 0 ? 0 : std::max(xn, yn);
//
//   R_xlen_t i, xi, yi;
//
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n)); ++NP;
//   int *p_out = LOGICAL(out);
//
//   R_xlen_t n_true = 0;
//   int eq = 0;
//
// switch (TYPEOF(x)){
// case LGLSXP:
// case INTSXP: {
//
//   switch (TYPEOF(y)){
// case LGLSXP:
// case INTSXP: {
//
//   int (*c_op)(int, int);
//   OP_SWITCH;
//
//   int *p_x = INTEGER(x);
//   int *p_y = INTEGER(y);
//
//
//   CHEAPR_COMPARISON_LOOP(p_x[xi], p_y[yi], cheapr_is_na_int, cheapr_is_na_int);
//   break;
// }
// case REALSXP: {
//
//   int (*c_op)(double, double);
//   OP_SWITCH;
//
//   int *p_x = INTEGER(x);
//   double *p_y = REAL(y);
//   CHEAPR_COMPARISON_LOOP(p_x[xi], p_y[yi], cheapr_is_na_int, cheapr_is_na_dbl);
//   break;
// }
// case STRSXP: {
//
//   int (*c_op)(std::string, std::string);
//   OP_SWITCH;
//
//   Rf_protect(x = Rf_coerceVector(x, STRSXP)); ++NP;
//
//   const SEXP *p_x = STRING_PTR_RO(x);
//   const SEXP *p_y = STRING_PTR_RO(y);
//   CHEAPR_COMPARISON_LOOP(Rf_translateCharUTF8(p_x[xi]),
//                          Rf_translateCharUTF8(p_y[yi]),
//                          cheapr_is_na_str,
//                          cheapr_is_na_str);
//   break;
// }
// default: {
//   Rf_unprotect(NP);
//   Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
// }
// }
//   break;
// }
// case REALSXP: {
//
//   switch (TYPEOF(y)){
//   case REALSXP: {
//
//     int (*c_op)(double, double);
//     OP_SWITCH;
//
//     double *p_x = REAL(x);
//     double *p_y = REAL(y);
//     CHEAPR_COMPARISON_LOOP(p_x[xi], p_y[yi], cheapr_is_na_dbl, cheapr_is_na_dbl);
//     break;
//   }
//   case INTSXP: {
//
//     int (*c_op)(double, double);
//     OP_SWITCH;
//
//     double *p_x = REAL(x);
//     int *p_y = INTEGER(y);
//     CHEAPR_COMPARISON_LOOP(p_x[xi], p_y[yi], cheapr_is_na_int, cheapr_is_na_int);
//     break;
//   }
//   case STRSXP: {
//
//     int (*c_op)(std::string, std::string);
//     OP_SWITCH;
//
//     Rf_protect(x = Rf_coerceVector(x, STRSXP)); ++NP;
//
//     const SEXP *p_x = STRING_PTR_RO(x);
//     const SEXP *p_y = STRING_PTR_RO(y);
//     CHEAPR_COMPARISON_LOOP(Rf_translateCharUTF8(p_x[xi]),
//                            Rf_translateCharUTF8(p_y[yi]),
//                            cheapr_is_na_str,
//                            cheapr_is_na_str);
//     break;
//   }
//   default: {
//     Rf_unprotect(1);
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
//   break;
// }
// case STRSXP: {
//
//   switch (TYPEOF(y)){
//
//   // No break in these cases to allow the STRSXP case to be used after coercion
// case LGLSXP:
// case INTSXP:
// case REALSXP: {
//   Rf_protect(y = Rf_coerceVector(y, STRSXP)); ++NP;
// }
// case STRSXP: {
//
//   int (*c_op)(std::string, std::string);
//   OP_SWITCH;
//
//   const SEXP *p_x = STRING_PTR_RO(x);
//   const SEXP *p_y = STRING_PTR_RO(y);
//   CHEAPR_COMPARISON_LOOP(Rf_translateCharUTF8(p_x[xi]),
//                          Rf_translateCharUTF8(p_y[yi]),
//                          cheapr_is_na_str,
//                          cheapr_is_na_str);
//   break;
// }
// default: {
//   Rf_unprotect(NP);
//   Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
// }
// }
//   break;
// }
// default: {
//   Rf_unprotect(NP);
//   Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
// }
// }
//
//   SEXP r_n_true = Rf_protect(Rf_ScalarInteger(n_true)); ++NP;
//
//   SEXP n_true_attrib = Rf_protect(Rf_installChar(Rf_mkChar(".n.true"))); ++NP;
//   Rf_setAttrib(out, n_true_attrib, r_n_true);
//   Rf_unprotect(NP);
//   return out;
// }

[[cpp11::register]]
SEXP cpp_character_compare(SEXP x, SEXP y, int op){
  if (!( Rf_isString(x) || Rf_isString(y) )){
    Rf_error("Either x or y must be a character vector");
  }

  int NP = 0;
  R_xlen_t xn = Rf_xlength(x);
  R_xlen_t yn = Rf_xlength(y);

  R_xlen_t n = xn == 0 || yn == 0 ? 0 : std::max(xn, yn);

  R_xlen_t i, xi, yi;

  SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n)); ++NP;
  int *p_out = LOGICAL(out);

  // int (*c_op)(int, int);
  // int (*c_op)(std::string, std::string);
  // OP_SWITCH;

  double tol = std::sqrt(std::numeric_limits<double>::epsilon());


  // It's faster to coerce to character if x or y is small enough
  // relative to the other

  if (n > 0 && xn <= (n / 2)){
    Rf_protect(x = Rf_coerceVector(x, STRSXP)); ++NP;
  }
  if (n > 0 && yn <= (n / 2)){
    Rf_protect(y = Rf_coerceVector(y, STRSXP)); ++NP;
  }

  switch (TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {

    switch (TYPEOF(y)){
  case LGLSXP:
  case INTSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    int *p_x = INTEGER(x);
    int *p_y = INTEGER(y);


    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(Rf_mkChar(std::to_string(p_x[xi]).c_str())),
      Rf_translateCharUTF8(Rf_mkChar(std::to_string(p_y[yi]).c_str())),
      cheapr_is_na_int,
      cheapr_is_na_int);
    break;
  }
  case REALSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    int *p_x = INTEGER(x);
    double *p_y = REAL(y);


    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(Rf_mkChar(std::to_string(p_x[xi]).c_str())),
      Rf_translateCharUTF8(
        Rf_mkChar(
          std::to_string(
            std::abs(p_y[yi]) < tol ? 0 : p_y[yi]
          ).c_str()
        )
      ),
      cheapr_is_na_int,
      cheapr_is_na_dbl);
    break;
  }
  case STRSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    int *p_x = INTEGER(x);
    const SEXP *p_y = STRING_PTR_RO(y);

    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(Rf_mkChar(std::to_string(p_x[xi]).c_str())),
      Rf_translateCharUTF8(p_y[yi]),
      cheapr_is_na_int,
      cheapr_is_na_str);
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
    break;
  }
  case REALSXP: {

    switch (TYPEOF(y)){
  case REALSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    double *p_x = REAL(x);
    double *p_y = REAL(y);

    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(
        Rf_mkChar(
          std::to_string(
            std::abs(p_x[xi]) < tol ?
    0 : p_x[xi]
          ).c_str()
        )
      ),
      Rf_translateCharUTF8(
        Rf_mkChar(
          std::to_string(
            std::abs(p_y[yi]) < tol ? 0 : p_y[yi]
          ).c_str()
        )
      ),
      cheapr_is_na_dbl,
      cheapr_is_na_dbl);
    break;
  }
  case INTSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    double *p_x = REAL(x);
    int *p_y = INTEGER(y);

    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(
        Rf_mkChar(
          std::to_string(
            std::abs(p_x[xi]) < tol ? 0 : p_x[xi]
          ).c_str()
        )
      ),
      Rf_translateCharUTF8(Rf_mkChar(std::to_string(p_y[yi]).c_str())),
      cheapr_is_na_dbl,
      cheapr_is_na_int);
    break;
  }
  case STRSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    double *p_x = REAL(x);
    const SEXP *p_y = STRING_PTR_RO(y);

    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(
        Rf_mkChar(
          std::to_string(
            std::abs(p_x[xi]) < tol ? 0 : p_x[xi]
          ).c_str()
        )
      ),
      Rf_translateCharUTF8(p_y[yi]),
      cheapr_is_na_dbl,
      cheapr_is_na_str);
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
    break;
  }
  case STRSXP: {

    switch (TYPEOF(y)){

    // No break in these cases to allow the STRSXP case to be used after coercion
  case LGLSXP:
  case INTSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    const SEXP *p_x = STRING_PTR_RO(x);
    int *p_y = INTEGER(y);

    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(p_x[xi]),
      Rf_translateCharUTF8(Rf_mkChar(std::to_string(p_y[yi]).c_str())),
      cheapr_is_na_str,
      cheapr_is_na_int);

    break;
  }
  case REALSXP: {

    int (*c_op)(std::string, std::string);
    OP_SWITCH;

    const SEXP *p_x = STRING_PTR_RO(x);
    double *p_y = REAL(y);

    CHEAPR_CHAR_COMPARISON_LOOP(
      Rf_translateCharUTF8(p_x[xi]),
      Rf_translateCharUTF8(
        Rf_mkChar(
          std::to_string(
            std::abs(p_y[yi]) < std::sqrt(std::numeric_limits<double>::epsilon()) ?
    0 : p_y[yi]
          ).c_str()
        )
      ),
      cheapr_is_na_str,
      cheapr_is_na_dbl);

    break;
  }
  case STRSXP: {

    int (*c_op)(int, int);
    OP_SWITCH;

    const SEXP *p_x = STRING_PTR_RO(x);
    const SEXP *p_y = STRING_PTR_RO(y);

    CHEAPR_CHAR_COMPARISON_LOOP2(
      Rf_translateCharUTF8(p_x[xi]),
      Rf_translateCharUTF8(p_y[yi]),
      cheapr_is_na_str,
      cheapr_is_na_str);
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
    break;
  }
  default: {
    Rf_unprotect(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  Rf_unprotect(NP);
  return out;
}

// #define cheapr_xor(a, b) ((bool) ((a + b) == 1))

// Essentially `x & y`
// SEXP cpp_and(SEXP x, SEXP y){
//   R_xlen_t xn = Rf_xlength(x);
//   R_xlen_t yn = Rf_xlength(y);
//
//   R_xlen_t n = xn == 0 || yn == 0 ? 0 : std::max(xn, yn);
//
//   R_xlen_t i, xi, yi;
//
//   int *p_x = LOGICAL(x);
//   int *p_y = LOGICAL(y);
//
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
//   int *p_out = LOGICAL(out);
//
//   // Copy the values of x over to out
//
//   memmove(p_out, &p_x[0], sizeof(int) * xn);
//
//   for (i = xi = yi = 0; i < n; xi = (++xi == xn) ? 0 : xi,
//        yi = (++yi == yn) ? 0 : yi, ++i){
//
//     // Is x[i] == FALSE then out[i] is always FALSE
//
//     if (p_x[xi] != FALSE){
//       if (p_y[yi] == FALSE){
//         p_out[i] = FALSE;
//       } else if ((p_x[xi] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
//         p_out[i] = NA_LOGICAL;
//       } else {
//         p_out[i] = p_x[xi] == TRUE && p_y[yi] == TRUE;
//       }
//     }
//
//   }
//
//   Rf_unprotect(1);
//   return out;
// }

// Essentially `x | y`

// SEXP cpp_or(SEXP x, SEXP y){
//   R_xlen_t xn = Rf_xlength(x);
//   R_xlen_t yn = Rf_xlength(y);
//
//   R_xlen_t n = xn == 0 || yn == 0 ? 0 : std::max(xn, yn);
//
//   R_xlen_t i, xi, yi;
//
//   int *p_x = LOGICAL(x);
//   int *p_y = LOGICAL(y);
//
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
//   int *p_out = LOGICAL(out);
//
//   // Copy the values of x over to out
//
//   memmove(p_out, &p_x[0], sizeof(int) * xn);
//
//   for (i = xi = yi = 0; i < n; xi = (++xi == xn) ? 0 : xi,
//        yi = (++yi == yn) ? 0 : yi, ++i){
//
//     // Is x[i] == TRUE then out[i] is always TRUE
//
//     if (p_x[xi] != TRUE){
//       if (p_y[yi] == TRUE){
//         p_out[i] = TRUE;
//       } else if ((p_x[xi] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
//         p_out[i] = NA_LOGICAL;
//       } else {
//         p_out[i] = p_x[xi] == TRUE || p_y[yi] == TRUE;
//       }
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }

// Faster (x | y) & !(x & y) found in xor()

// SEXP cpp_xor(SEXP x, SEXP y){
//   R_xlen_t xn = Rf_xlength(x);
//   R_xlen_t yn = Rf_xlength(y);
//
//   R_xlen_t n = xn == 0 || yn == 0 ? 0 : std::max(xn, yn);
//
//   R_xlen_t i, xi, yi;
//
//   int *p_x = LOGICAL(x);
//   int *p_y = LOGICAL(y);
//
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
//   int *p_out = LOGICAL(out);
//
//   for (i = xi = yi = 0; i < n; xi = (++xi == xn) ? 0 : xi,
//        yi = (++yi == yn) ? 0 : yi, ++i){
//
//     if ((p_x[xi] == NA_LOGICAL) || (p_y[yi] == NA_LOGICAL)){
//       p_out[i] = NA_LOGICAL;
//     } else {
//       p_out[i] = p_x[xi] != p_y[yi];
//       // p_out[i] = cheapr_xor(p_x[xi], p_y[yi]);
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }

#undef OP_SWITCH
#undef CHEAPR_CHAR_COMPARISON_LOOP
#undef CHEAPR_CHAR_COMPARISON_LOOP2
// #undef cheapr_xor
