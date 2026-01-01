#include "cheapr.h"
#include <R_ext/Altrep.h>

// Altrep utils

// Symbols
static SEXP CHEAPR_COMPACT_INTSEQ = r_null;
static SEXP CHEAPR_COMPACT_REALSEQ = r_null;
static SEXP CHEAPR_BASE = r_null;

[[cpp11::register]]
SEXP alt_class(SEXP x){
  if (altrep::is_altrep(x) && has_attrs(ALTREP_CLASS(x))){
    return VECTOR_ELT(get_attrs(ALTREP_CLASS(x)), 0);
  } else {
    return r_null;
  }
}

[[cpp11::register]]
SEXP foo2(SEXP x){
  SEXP altrep_cls = SHIELD(ALTREP_CLASS(x));
  bool alt_has_attrs = has_attrs(altrep_cls);

  if (altrep::is_altrep(x) && alt_has_attrs){
    SEXP attrs = SHIELD(get_attrs(altrep_cls));
    SEXP out = SHIELD(VECTOR_ELT(attrs, 0));
    YIELD(3);
    return out;
  } else {
    YIELD(1);
    return r_null;
  }
}


[[cpp11::register]]
SEXP foo3(SEXP x){
  SEXP expr = Rf_protect(Rf_lang2(Rf_install("attributes"), x));
  SEXP out = Rf_protect(Rf_eval(expr, R_BaseEnv));
  Rf_unprotect(2);
  return out;
}


[[cpp11::register]]
bool foo4(SEXP x){
  SEXP altrep_cls = SHIELD(ALTREP_CLASS(x));
  bool alt_has_attrs = has_attrs(altrep_cls);
  YIELD(1);
  return true;
}

[[cpp11::register]]
bool foo5(SEXP x){
  SEXP altrep_cls = SHIELD(ALTREP_CLASS(x));
  bool alt_has_attrs = has_attrs(altrep_cls);
  SEXP attrs = SHIELD(get_attrs(altrep_cls));
  YIELD(2);
  return true;
}
[[cpp11::register]]
SEXP alt_pkg(SEXP x){
  if (altrep::is_altrep(x) && has_attrs(ALTREP_CLASS(x))){
    return VECTOR_ELT(get_attrs(ALTREP_CLASS(x)), 1);
  } else {
    return r_null;
  }
}
[[cpp11::register]]
SEXP alt_data1(SEXP x){
  if (altrep::is_altrep(x)){
    return R_altrep_data1(x);
  } else {
    return r_null;
  }
}

bool is_compact_seq(SEXP x){
  if (!altrep::is_altrep(x)) return false;
  SEXP alt_class_sym = alt_class(x);
  SEXP alt_pkg_sym = alt_pkg(x);

  if (is_null(CHEAPR_COMPACT_INTSEQ)){
    CHEAPR_COMPACT_INTSEQ = r_cast<r_symbol_t>("compact_intseq");
  }
  if (is_null(CHEAPR_COMPACT_REALSEQ)){
    CHEAPR_COMPACT_REALSEQ = r_cast<r_symbol_t>("compact_realseq");
  }
  if (is_null(CHEAPR_BASE)){
    CHEAPR_BASE = r_cast<r_symbol_t>("base");
  }
  return (alt_class_sym == CHEAPR_COMPACT_INTSEQ ||
          alt_class_sym == CHEAPR_COMPACT_REALSEQ) &&
          alt_pkg_sym == CHEAPR_BASE;
}

[[cpp11::register]]
SEXP compact_seq_data(SEXP x){
  if (!is_compact_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = SHIELD(vec::coerce_vec(alt_data1(x), REALSXP));
  double *p_alt_data = real_ptr(alt_data);
  double size = p_alt_data[0];
  double from = p_alt_data[1];
  double by = p_alt_data[2];
  double to = (std::fmax(size - 1.0, 0.0) * by) + from;
  SEXP out = SHIELD(combine(from, to, by, size));
  YIELD(2);
  return out;
}

[[cpp11::register]]
SEXP foobarfoo(SEXP x){
  if (!is_compact_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = SHIELD(vec::coerce_vec(alt_data1(x), REALSXP));
  double *p_alt_data = real_ptr(alt_data);
  double size = p_alt_data[0];
  double from = p_alt_data[1];
  double by = p_alt_data[2];
  double to = (std::fmax(size - 1.0, 0.0) * by) + from;
  YIELD(1);
  return r_null;
}

SEXP altrep_materialise(SEXP x) {
  return altrep::is_altrep(x) ? cpp_semi_copy(x) : x;
}
