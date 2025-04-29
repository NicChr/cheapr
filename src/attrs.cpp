#include "cheapr.h"

// Adding and removing attributes in-place
// There is a check to ensure that attributes are copied when they are the same
// object as x

[[cpp11::register]]
SEXP cpp_set_rm_attributes(SEXP x){
  SEXP current = ATTRIB(x);
  while (current != R_NilValue){
    Rf_setAttrib(x, TAG(current), R_NilValue);
    current = CDR(current);
  }
  return x;
}

// Add attribute onto existing attributes

[[cpp11::register]]
SEXP cpp_set_add_attr(SEXP x, SEXP which, SEXP value) {
  if (Rf_length(which) != 1){
    Rf_error("`which` must be a character vector of length 1 in %s", __func__);
  }
  SEXP attr_char = SHIELD(Rf_install(utf8_char(STRING_ELT(which, 0))));
  SEXP value2 = SHIELD(address_equal(x, value) ? Rf_duplicate(value) : value);
  Rf_setAttrib(x, attr_char, value2);
  YIELD(2);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_rm_attr(SEXP x, SEXP which){
  Rf_setAttrib(x, Rf_installChar(STRING_ELT(which, 0)), R_NilValue);
  return x;
}

// Set attributes of x in-place, when add = F, attrs of x are first removed

[[cpp11::register]]
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add) {

  if (!add) cpp_set_rm_attributes(x);

  int NP = 0;

  if (is_null(attributes)){
    return x;
  } else if (TYPEOF(attributes) == VECSXP){
    if (Rf_length(attributes) == 0) return x;
    SEXP names = Rf_getAttrib(attributes, R_NamesSymbol);
    if (is_null(names)){
      Rf_error("attributes must be a named list");
    }
    const SEXP *p_attributes = VECTOR_PTR_RO(attributes);
    const SEXP *p_names = STRING_PTR_RO(names);

    SEXP attr_nm = R_NilValue;

    for (int i = 0; i < Rf_length(names); ++i){
      if (p_names[i] != R_BlankString){
        attr_nm = Rf_install(utf8_char(p_names[i]));
        if (address_equal(x, p_attributes[i])){
          SEXP dup_attr = SHIELD(Rf_duplicate(p_attributes[i])); ++NP;
          Rf_setAttrib(x, attr_nm, dup_attr);
        } else {
          Rf_setAttrib(x, attr_nm, p_attributes[i]);
        }
      }
    }
    YIELD(NP);
    return x;
  } else if (TYPEOF(attributes) == LISTSXP){
    SEXP current = attributes;

    while (!is_null(current)){
      if (address_equal(x, CAR(current))){
        SEXP dup_attr = SHIELD(Rf_duplicate(CAR(current))); ++NP;
        Rf_setAttrib(x, TAG(current), dup_attr);
      } else {
        Rf_setAttrib(x, TAG(current), CAR(current));
      }
     current = CDR(current);
    }
    YIELD(NP);
    return x;
  } else {
    Rf_error("`attributes` must be a named list");
  }
}

// Copy attributes from source to target
void cpp_copy_attributes(SEXP source, SEXP target, bool deep_copy){
  SEXP target_attrs = SHIELD(deep_copy ? Rf_duplicate(ATTRIB(source)) : ATTRIB(source));
  cpp_set_add_attributes(target, target_attrs, false);
  YIELD(1);
}

// Copy names from source to target
void cpp_copy_names(SEXP source, SEXP target, bool deep_copy){
  SEXP source_nms = Rf_getAttrib(source, R_NamesSymbol);
  SEXP target_nms = SHIELD(deep_copy ? Rf_duplicate(source_nms) : source_nms);
  if (!is_null(source_nms)){
    Rf_setAttrib(target, R_NamesSymbol, target_nms);
  }
  YIELD(1);
}
