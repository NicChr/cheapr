#include "cheapr.h"

// Adding and removing attributes in-place
// There is a check to ensure that attributes are copied when they are the same
// object as x

[[cpp11::register]]
SEXP cpp_set_rm_attributes(SEXP x){
  SEXP attrs = SHIELD(ATTRIB(x));
  SEXP names = SHIELD(Rf_getAttrib(attrs, R_NamesSymbol));
  int n = Rf_length(attrs);
  for (int i = 0; i < n; ++i){
    Rf_setAttrib(x, Rf_installChar(STRING_ELT(names, i)), R_NilValue);
  }
  YIELD(2);
  return x;
}

// Add attribute onto existing attributes

[[cpp11::register]]
SEXP cpp_set_add_attr(SEXP x, SEXP which, SEXP value) {
  SEXP attr_char = SHIELD(Rf_install(CHAR(Rf_asChar(which))));
  SEXP value2 = SHIELD(r_address(x) == r_address(value) ? Rf_duplicate(value) : value);
  Rf_setAttrib(x, attr_char, value2);
  YIELD(2);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_rm_attr(SEXP x, SEXP which) {
  Rf_setAttrib(x, Rf_installChar(Rf_asChar(which)), R_NilValue);
  return x;
}

// Set attributes of x in-place, when add = F, attrs of x are first removed

[[cpp11::register]]
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add) {
  int NP = 0;
  SEXP attrs;
  if (Rf_isPairList(attributes)){
    attrs = SHIELD(coerce_vec(attributes, VECSXP)); ++NP;
  } else {
   attrs = attributes;
  }
  int n = Rf_length(attrs);
  bool attrs_are_a_list = TYPEOF(attrs) == VECSXP;
  if (Rf_isNull(attrs) ||
      // is.null or empty list?
      (attrs_are_a_list && n == 0)){
    if (add){
      YIELD(NP);
      return x;
    } else {
      YIELD(NP);
      return cpp_set_rm_attributes(x);
    }
  }
  SEXP names = SHIELD(Rf_getAttrib(attrs, R_NamesSymbol)); ++NP;
  if (!attrs_are_a_list || Rf_isNull(names)){
    YIELD(NP);
    Rf_error("attributes must be a named list");
  }
  if (!add) cpp_set_rm_attributes(x);
  const SEXP *p_attributes = VECTOR_PTR_RO(attrs);
  const SEXP *p_names = STRING_PTR_RO(names);
  for (int i = 0; i < n; ++i){
    SEXP attr_nm = SHIELD(Rf_installChar(p_names[i])); ++NP;
    if (r_address(x) == r_address(p_attributes[i])){
      SEXP dup_attr = SHIELD(Rf_duplicate(p_attributes[i])); ++NP;
      Rf_setAttrib(x, attr_nm, dup_attr);
    } else {
      Rf_setAttrib(x, attr_nm, p_attributes[i]);
    }
  }
  YIELD(NP);
  return x;
}

// Copy attributes from source to target
void cpp_copy_attributes(SEXP source, SEXP target, bool deep_copy){
  SEXP target_attrs = SHIELD(deep_copy ? Rf_duplicate(ATTRIB(source)) : ATTRIB(source));
  cpp_set_add_attributes(target, target_attrs, false);
  YIELD(1);
}

// Copy names from source to target
void cpp_copy_names(SEXP source, SEXP target, bool deep_copy){
  SEXP source_nms = SHIELD(Rf_getAttrib(source, R_NamesSymbol));
  SEXP target_nms = SHIELD(deep_copy ? Rf_duplicate(source_nms) : source_nms);
  if (!Rf_isNull(source_nms)){
    Rf_setAttrib(target, R_NamesSymbol, target_nms);
  }
  YIELD(2);
}

[[cpp11::register]]
void cpp_shallow_duplicate_attrs(SEXP source, SEXP target){
  SHALLOW_DUPLICATE_ATTRIB(target, source);
}
[[cpp11::register]]
void cpp_copy_most_attrs(SEXP source, SEXP target){
  Rf_copyMostAttrib(source, target);
}
