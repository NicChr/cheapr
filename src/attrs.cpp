#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>

// Adding and removing attributes in-place
// There is a check to ensure that attributes are copied when they are the same
// object as x

[[cpp11::register]]
SEXP cpp_set_rm_attributes(SEXP x){
  SEXP attrs = Rf_protect(ATTRIB(x));
  SEXP names = Rf_protect(Rf_getAttrib(attrs, R_NamesSymbol));
  int n = Rf_length(attrs);
  for (int i = 0; i < n; ++i){
    SEXP attrib_nm = Rf_protect(Rf_installChar(STRING_ELT(names, i)));
    Rf_setAttrib(x, attrib_nm, R_NilValue);
  }
  Rf_unprotect(n + 2);
  return x;
}

// Add attribute onto existing attributes

[[cpp11::register]]
SEXP cpp_set_add_attr(SEXP x, SEXP which, SEXP value) {
  int n_protect;
  Rf_protect(x = x);
  Rf_protect(which = which);
  Rf_protect(value = value);
  SEXP attr_char = Rf_protect(Rf_install(CHAR(STRING_ELT(which, 0))));
  if (cpp_obj_address(x) == cpp_obj_address(value)){
    Rf_protect(value = Rf_duplicate(value));
    n_protect = 5;
  } else {
    n_protect = 4;
  }
  Rf_setAttrib(x, attr_char, value);
  Rf_unprotect(n_protect);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_rm_attr(SEXP x, SEXP which) {
  Rf_protect(x = x);
  Rf_protect(which = which);
  SEXP attr_char = Rf_protect(Rf_installChar(STRING_ELT(which, 0)));
  Rf_setAttrib(x, attr_char, R_NilValue);
  Rf_unprotect(3);
  return x;
}

// Set attributes of x in-place, when add = F, attrs of x are first removed

[[cpp11::register]]
SEXP cpp_set_attributes(SEXP x, SEXP attributes, bool add) {
  int n_protect;
  if (add){
    Rf_protect(x = x);
  } else {
    Rf_protect(x = cpp_set_rm_attributes(x));
  }
  SEXP names = Rf_protect(Rf_getAttrib(attributes, R_NamesSymbol));
  n_protect = 2;
  if (!Rf_isVectorList(attributes) || Rf_isNull(names)){
    Rf_unprotect(n_protect);
    Rf_error("attributes must be a named list");
  }
  const SEXP *p_attributes = VECTOR_PTR_RO(attributes);
  SEXP *p_names = STRING_PTR(names);
  int n = Rf_length(attributes);
  for (int i = 0; i < n; ++i){
    SEXP attr_nm = Rf_protect(Rf_installChar(p_names[i]));
    ++n_protect;
    if (cpp_obj_address(x) == cpp_obj_address(p_attributes[i])){
      SEXP dup_attr = Rf_protect(Rf_duplicate(p_attributes[i]));
      ++n_protect;
      Rf_setAttrib(x, attr_nm, dup_attr);
    } else {
      Rf_setAttrib(x, attr_nm, p_attributes[i]);
    }
  }
  Rf_unprotect(n_protect);
  return x;
}
