#ifndef CHEAPR_R_ATTRS_H
#define CHEAPR_R_ATTRS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_symbols.h>

namespace cheapr {

namespace attr {

inline bool inherits1(SEXP x, const char *r_cls){
  return Rf_inherits(x, r_cls);
}

// Attributes of x as a list
inline SEXP get_attrs(SEXP x){
  if (internal::BASE_ATTRIBUTES == NULL){
    internal::BASE_ATTRIBUTES = Rf_install("attributes");
  }
  SEXP expr = SHIELD(Rf_lang2(internal::BASE_ATTRIBUTES, x));
  SEXP out = SHIELD(Rf_eval(expr, R_BaseEnv));
  YIELD(2);
  return out;
}

inline bool has_attrs(SEXP x){
  return !is_null(get_attrs(x));
}

inline SEXP get_attr(SEXP x, r_symbol_t sym){
  return Rf_getAttrib(x, static_cast<SEXP>(sym));
}

inline void set_attr(SEXP x, r_symbol_t sym, SEXP value){
  Rf_setAttrib(x, static_cast<SEXP>(sym), value);
}

inline void set_old_names(SEXP x, SEXP names){
  if (is_null(names)){
    attr::set_attr(x, symbol::names_sym, r_null); 
  } else {
    Rf_namesgets(x, names);
  }
}
inline SEXP get_old_names(SEXP x){
  return attr::get_attr(x, symbol::names_sym);
}
inline bool has_r_names(SEXP x){
  SEXP names = SHIELD(get_old_names(x));
  bool out = !is_null(names);
  YIELD(1);
  return out;
}

inline SEXP get_old_class(SEXP x){
  return get_attr(x, symbol::class_sym);
}
inline void set_old_class(SEXP x, SEXP cls){
  Rf_classgets(x, cls); 
}

inline bool inherits(SEXP x, SEXP classes){
  R_xlen_t n = Rf_xlength(classes);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (inherits1(x, CHAR(STRING_ELT(classes, i)))){
      return true;
    }
  }
  return false;
}

}

}

#endif
