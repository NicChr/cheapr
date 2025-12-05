#include "cheapr.h"

// Altrep utils

// Symbols
static SEXP CHEAPR_COMPACT_INTSEQ = r_null;
static SEXP CHEAPR_COMPACT_REALSEQ = r_null;
static SEXP CHEAPR_BASE = r_null;

SEXP alt_class(SEXP x){
  if (is_altrep(x)){
    return CAR(ATTRIB(ALTREP_CLASS(x)));
  } else {
    return r_null;
  }
}
SEXP alt_pkg(SEXP x){
  if (is_altrep(x)){
    return CADR(ATTRIB(ALTREP_CLASS(x)));
  } else {
    return r_null;
  }
}

SEXP alt_data1(SEXP x){
  if (is_altrep(x)){
    return R_altrep_data1(x);
  } else {
    return r_null;
  }
}

bool is_compact_seq(SEXP x){
  if (!is_altrep(x)) return false;
  SEXP alt_class_sym = alt_class(x);
  SEXP alt_pkg_sym = alt_pkg(x);

  if (is_null(CHEAPR_COMPACT_INTSEQ)){
    CHEAPR_COMPACT_INTSEQ = install_utf8("compact_intseq");
  }
  if (is_null(CHEAPR_COMPACT_REALSEQ)){
    CHEAPR_COMPACT_REALSEQ = install_utf8("compact_realseq");
  }
  if (is_null(CHEAPR_BASE)){
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
  SEXP alt_data = SHIELD(vec::coerce_vec(alt_data1(x), REALSXP));
  double *p_alt_data = REAL(alt_data);
  double size = p_alt_data[0];
  double from = p_alt_data[1];
  double by = p_alt_data[2];
  double to = (std::fmax(size - 1.0, 0.0) * by) + from;
  SEXP out = SHIELD(new_r_vec(from, to, by, size));
  YIELD(2);
  return out;
}

SEXP altrep_materialise(SEXP x) {
  return is_altrep(x) ? cpp_semi_copy(x) : x;
}
