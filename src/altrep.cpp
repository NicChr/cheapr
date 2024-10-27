#include "cheapr.h"

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

// bool is_compact_seq(SEXP x){
//   if (!is_altrep(x)) return false;
//   SEXP alt_attribs = Rf_protect(Rf_coerceVector(ATTRIB(ALTREP_CLASS(x)), VECSXP));
//   std::string alt_class_str = CHAR(PRINTNAME(VECTOR_ELT(alt_attribs, 0)));
//   std::string alt_pkg_str = CHAR(PRINTNAME(VECTOR_ELT(alt_attribs, 1)));
//   std::string intseq_str = "compact_intseq";
//   std::string realseq_str = "compact_realseq";
//   std::string base_str = "base";
//   bool out = alt_pkg_str == base_str && (alt_class_str == intseq_str || alt_class_str == realseq_str);
//   Rf_unprotect(1);
//   return out;
// }

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
  // if (ALTREP(x)){
  //   return TYPEOF(x) == VECSXP ? list_shallow_copy(x, false) : Rf_duplicate(x);
  // } else {
  //  return x;
  // }
  return ALTREP(x) ? Rf_duplicate(x) : x;
  // Even after using DATAPTR, ALTREP(x) == TRUE ?
  // if (ALTREP(x)) DATAPTR(x);
  // return x;
}
