#include "cheapr_cpp.h"

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
  SEXP x2 = Rf_protect(x);
  SEXP which2 = Rf_protect(which);
  SEXP value2 = Rf_protect(value);
  SEXP attr_char = Rf_protect(Rf_install(CHAR(STRING_ELT(which2, 0))));
  if (r_address(x2) == r_address(value2)){
    Rf_protect(value2 = Rf_duplicate(value2));
    n_protect = 5;
  } else {
    n_protect = 4;
  }
  Rf_setAttrib(x2, attr_char, value2);
  Rf_unprotect(n_protect);
  return x2;
}

[[cpp11::register]]
SEXP cpp_set_rm_attr(SEXP x, SEXP which) {
  SEXP x2 = Rf_protect(x);
  SEXP which2 = Rf_protect(which);
  SEXP attr_char = Rf_protect(Rf_installChar(STRING_ELT(which2, 0)));
  Rf_setAttrib(x2, attr_char, R_NilValue);
  Rf_unprotect(3);
  return x2;
}

// Set attributes of x in-place, when add = F, attrs of x are first removed

[[cpp11::register]]
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add) {
  int n_protect = 0;
  SEXP x2 = Rf_protect(x);
  ++n_protect;
  if (!add){
    Rf_protect(x2 = cpp_set_rm_attributes(x2));
    ++n_protect;
  }
  SEXP names = Rf_protect(Rf_getAttrib(attributes, R_NamesSymbol));
  ++n_protect;
  if (!Rf_isVectorList(attributes) || Rf_isNull(names)){
    Rf_unprotect(n_protect);
    Rf_error("attributes must be a named list");
  }
  const SEXP *p_attributes = VECTOR_PTR_RO(attributes);
  const SEXP *p_names = STRING_PTR_RO(names);
  int n = Rf_length(attributes);
  for (int i = 0; i < n; ++i){
    SEXP attr_nm = Rf_protect(Rf_installChar(p_names[i]));
    ++n_protect;
    if (r_address(x) == r_address(p_attributes[i])){
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
