#include "cheapr_cpp.h"
#include <R_ext/Altrep.h>

// Altrep utils

bool is_altrep(SEXP x){
  return ALTREP(x);
}
SEXP alt_class(SEXP x){
  if (is_altrep(x)){
    SEXP alt_attribs = Rf_protect(Rf_coerceVector(ATTRIB(ALTREP_CLASS(x)), VECSXP));
    SEXP out = Rf_protect(Rf_coerceVector(VECTOR_ELT(alt_attribs, 0), STRSXP));
    Rf_unprotect(2);
    return out;
  } else {
    return R_NilValue;
  }
}
SEXP alt_pkg(SEXP x){
  if (is_altrep(x)){
    SEXP alt_attribs = Rf_protect(Rf_coerceVector(ATTRIB(ALTREP_CLASS(x)), VECSXP));
    SEXP out = Rf_protect(Rf_coerceVector(VECTOR_ELT(alt_attribs, 1), STRSXP));
    Rf_unprotect(2);
    return out;
  } else {
    return R_NilValue;
  }
}

SEXP alt_data1(SEXP x){
  if (is_altrep(x)){
    return R_altrep_data1(x);
  } else {
    return R_NilValue;
  }
}

[[cpp11::register]]
bool is_compact_seq(SEXP x){
  if (!is_altrep(x)) return false;
  SEXP alt_class_nm = Rf_protect(alt_class(x));
  SEXP alt_pkg_nm = Rf_protect(alt_pkg(x));
  SEXP intseq_char = Rf_protect(Rf_mkChar("compact_intseq"));
  SEXP realseq_char = Rf_protect(Rf_mkChar("compact_realseq"));
  SEXP base_char = Rf_protect(Rf_mkChar("base"));
  bool out = (STRING_ELT(alt_class_nm, 0) == intseq_char ||
              STRING_ELT(alt_class_nm, 0) == realseq_char) &&
              STRING_ELT(alt_pkg_nm, 0) == base_char;
  Rf_unprotect(5);
  return out;
}

[[cpp11::register]]
SEXP compact_seq_data(SEXP x){
  if (!is_compact_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = Rf_protect(Rf_coerceVector(alt_data1(x), REALSXP));
  double alt_size = REAL(alt_data)[0];
  double alt_from = REAL(alt_data)[1];
  double alt_by = REAL(alt_data)[2];
  double alt_to = (std::fmax(alt_size - 1.0, 0.0) * alt_by) + alt_from;
  SEXP out = Rf_protect(Rf_allocVector(REALSXP, 4));
  double *p_out = REAL(out);
  p_out[0] = alt_from;
  p_out[1] = alt_to;
  p_out[2] = alt_by;
  p_out[3] = alt_size;
  Rf_unprotect(2);
  return out;
}

SEXP altrep_materialise(SEXP x) {
  if (ALTREP(x)){
    return Rf_duplicate(x);
  } else {
    return x;
  }
  // Even after using DATAPTR, ALTREP(x) == TRUE ?
  // if (ALTREP(x)) DATAPTR(x);
  // return x;
}

// SEXP altrep_materialise(SEXP x) {
//   if (!ALTREP(x)) return x;
//   SEXP out;
//   R_xlen_t n = Rf_xlength(x);
//   switch (TYPEOF(x)) {
//   case INTSXP: {
//     out = Rf_protect(Rf_allocVector(INTSXP, n));
//     INTEGER_GET_REGION(x, 0, n, INTEGER(out));
//     break;
//   }
//   case REALSXP: {
//     out = Rf_protect(Rf_allocVector(REALSXP, n));
//     REAL_GET_REGION(x, 0, n, REAL(out));
//     break;
//   }
//   default: {
//     // This is a pretty hacky way of doing it..
//     cpp11::function base_subset = cpp11::package("base")["["];
//     out = Rf_protect(base_subset(x));
//   }
//   }
//   Rf_unprotect(1);
//   return out;
// }
