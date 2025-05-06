#include "cheapr.h"

// Altrep utils

inline bool is_altrep(SEXP x){
  return ALTREP(x);
}

// Symbols
static SEXP CHEAPR_COMPACT_INTSEQ = NULL;
static SEXP CHEAPR_COMPACT_REALSEQ = NULL;
static SEXP CHEAPR_BASE = NULL;

SEXP alt_class(SEXP x){
  if (is_altrep(x)){
    return CAR(ATTRIB(ALTREP_CLASS(x)));
  } else {
    return R_NilValue;
  }
}
SEXP alt_pkg(SEXP x){
  if (is_altrep(x)){
    return CADR(ATTRIB(ALTREP_CLASS(x)));
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

bool is_compact_seq(SEXP x){
  if (!is_altrep(x)) return false;
  SEXP alt_class_sym = alt_class(x);
  SEXP alt_pkg_sym = alt_pkg(x);

  if (CHEAPR_COMPACT_INTSEQ == NULL){
    CHEAPR_COMPACT_INTSEQ = install_utf8("compact_intseq");
  }
  if (CHEAPR_COMPACT_REALSEQ == NULL){
    CHEAPR_COMPACT_REALSEQ = install_utf8("compact_realseq");
  }
  if (CHEAPR_BASE == NULL){
    CHEAPR_BASE = install_utf8("base");
  }
  return (alt_class_sym == CHEAPR_COMPACT_INTSEQ ||
          alt_class_sym == CHEAPR_COMPACT_REALSEQ) &&
          alt_pkg_sym == CHEAPR_BASE;
}

SEXP compact_seq_data(SEXP x){
  if (!is_compact_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = SHIELD(coerce_vec(alt_data1(x), REALSXP));
  double alt_size = REAL(alt_data)[0];
  double alt_from = REAL(alt_data)[1];
  double alt_by = REAL(alt_data)[2];
  double alt_to = (std::fmax(alt_size - 1.0, 0.0) * alt_by) + alt_from;
  SEXP out = SHIELD(new_vec(REALSXP, 4));
  double *p_out = REAL(out);
  p_out[0] = alt_from;
  p_out[1] = alt_to;
  p_out[2] = alt_by;
  p_out[3] = alt_size;
  YIELD(2);
  return out;
}

SEXP altrep_materialise(SEXP x) {
  // if (ALTREP(x)){
  //   return TYPEOF(x) == VECSXP ? list_shallow_copy(x, false) : Rf_duplicate(x);
  // } else {
  //  return x;
  // }
  return is_altrep(x) ? cpp_semi_copy(x) : x;
  // Even after using DATAPTR, ALTREP(x) == TRUE ?
  // if (ALTREP(x)) DATAPTR(x);
  // return x;
}
