#include "cheapr.h"

// Adding and removing attributes in-place
// There is a check to ensure that attributes are copied when they are the same
// object as x


// Taken from Writing R Extensions
// void CLEAR_ATTRIB(SEXP x){
//   SET_ATTRIB(x, r_null);
//   SET_OBJECT(x, 0);
//   UNSET_S4_OBJECT(x);
// }


// Can't use the above because UNSET_S4_OBJECT is non-API
// and cant use CLEAR_ATTRIB because it's only available in latest R
// Also can't remove S4 bit in-place using any method
// The only way to basically do it is to copy R's header files and
// manipulate SEXP info directly

[[cpp11::register]]
SEXP cpp_set_rm_attributes(SEXP x){
  SEXP current = ATTRIB(x);
  while (!is_null(current)){
    set_attrib(x, TAG(current), r_null);
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
  SEXP attr_char = SHIELD(make_symbol(CHAR(STRING_ELT(which, 0))));
  SEXP value2 = SHIELD(address_equal(x, value) ? vec::deep_copy(value) : value);
  set_attrib(x, attr_char, value2);
  YIELD(2);
  return x;
}

[[cpp11::register]]
SEXP cpp_set_rm_attr(SEXP x, SEXP which){
  if (Rf_length(which) != 1){
    Rf_error("`which` must be a length 1 character vector");
  }
  if (TYPEOF(which) != STRSXP){
    Rf_error("`which` must be a length 1 character vector");
  }
  set_attrib(x, make_symbol(CHAR(STRING_ELT(which, 0))), r_null);
  return x;
}

// Set attributes of x in-place, when add = F, attrs of x are first removed

[[cpp11::register]]
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add) {

  if (!add) attr::clear_attrs(x);

  int32_t NP = 0;

  if (is_null(attributes)){
    return x;
  } else if (TYPEOF(attributes) == VECSXP){
    if (Rf_length(attributes) == 0) return x;
    SEXP names = SHIELD(get_r_names(attributes)); ++NP;
    if (is_null(names)){
      YIELD(NP);
      Rf_error("attributes must be a named list");
    }
    const SEXP *p_attributes = LIST_PTR_RO(attributes);
    const SEXP *p_names = STRING_PTR_RO(names);

    SEXP attr_nm = r_null;

    for (int i = 0; i < Rf_length(names); ++i){
      if (p_names[i] != R_BlankString){
        attr_nm = make_symbol(CHAR(p_names[i]));
        if (address_equal(x, p_attributes[i])){
          SEXP dup_attr = SHIELD(vec::deep_copy(p_attributes[i])); ++NP;
          set_attrib(x, attr_nm, dup_attr);
        } else {
          set_attrib(x, attr_nm, p_attributes[i]);
        }
      }
    }
    YIELD(NP);
    return x;
  } else if (TYPEOF(attributes) == LISTSXP){
    SEXP current = attributes;

    while (!is_null(current)){
      if (address_equal(x, CAR(current))){
        SEXP dup_attr = SHIELD(vec::deep_copy(CAR(current))); ++NP;
        set_attrib(x, TAG(current), dup_attr);
      } else {
        set_attrib(x, TAG(current), CAR(current));
      }
     current = CDR(current);
    }
    YIELD(NP);
    return x;
  } else {
    Rf_error("`attributes` must be a named list");
  }
}
