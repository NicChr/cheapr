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
  SEXP attr_char = Rf_protect(Rf_installChar(STRING_ELT(which, 0)));
  SEXP value2 = Rf_protect(r_address(x) == r_address(value) ? Rf_duplicate(value) : value);
  Rf_setAttrib(x, attr_char, value2);
  Rf_unprotect(2);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_rm_attr(SEXP x, SEXP which) {
  SEXP attr_char = Rf_protect(Rf_installChar(STRING_ELT(which, 0)));
  Rf_setAttrib(x, attr_char, R_NilValue);
  Rf_unprotect(1);
  return x;
}

// Set attributes of x in-place, when add = F, attrs of x are first removed

[[cpp11::register]]
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add) {
  int NP = 0;
  if (Rf_isNull(attributes) ||
      // is.null or empty list?
      (TYPEOF(attributes) == VECSXP && Rf_length(attributes) == 0)){
    if (add){
      return x;
    } else {
      return cpp_set_rm_attributes(x);
    }
  }
  if (!add){
    cpp_set_rm_attributes(x);
  }
  SEXP names = Rf_protect(Rf_getAttrib(attributes, R_NamesSymbol)); ++NP;
  if (!Rf_isVectorList(attributes) || Rf_isNull(names)){
    Rf_unprotect(NP);
    Rf_error("attributes must be a named list");
  }
  const SEXP *p_attributes = VECTOR_PTR_RO(attributes);
  const SEXP *p_names = STRING_PTR_RO(names);
  int n = Rf_length(attributes);
  for (int i = 0; i < n; ++i){
    SEXP attr_nm = Rf_protect(Rf_installChar(p_names[i])); ++NP;
    if (r_address(x) == r_address(p_attributes[i])){
      SEXP dup_attr = Rf_protect(Rf_duplicate(p_attributes[i])); ++NP;
      Rf_setAttrib(x, attr_nm, dup_attr);
    } else {
      Rf_setAttrib(x, attr_nm, p_attributes[i]);
    }
  }
  Rf_unprotect(NP);
  return x;
}

void cpp_copy_attributes(SEXP source, SEXP target){
  SEXP source_attrs = Rf_protect(ATTRIB(source));
  SEXP target_attrs = Rf_protect(Rf_duplicate(Rf_coerceVector(source_attrs, VECSXP)));
  cpp_set_add_attributes(target, target_attrs, false);
  Rf_unprotect(2);
}
