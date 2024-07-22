// #include "cheapr_cpp.h"

// THis is a helper function to assist with doing multiple comparisons
// across binary pairs of vectors between 2 data frames
// The idea is that rows are equal when all values across columns are equal
// And so we create a logical outside this function for the comparisons
// in the first column, and then we simply replace our logical vector
// With comparisons for the rest of the vector pairs, without making more copies.
// For large data frames this should be quite cheap


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
// [[cpp11::register]]
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
//       // either_na = is_na_cplx(a) || is_na_cplx(b);
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
