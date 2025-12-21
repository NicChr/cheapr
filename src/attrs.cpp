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
  cheapr::attr::clear_attrs(x);
  return x;
}

// Add attribute onto existing attributes

[[cpp11::register]]
SEXP cpp_set_add_attr(SEXP x, SEXP which, SEXP value) {
  if (Rf_length(which) != 1){
    Rf_error("`which` must be a character vector of length 1 in %s", __func__);
  }
  r_symbol_t attr_char = r_cast<r_symbol_t>(get_value<r_string_t>(which, 0));
  SEXP value2 = SHIELD(address_equal(x, value) ? vec::deep_copy(value) : value);
  set_attr(x, attr_char, value2);
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
  set_attr(x, r_cast<r_symbol_t>(get_value<r_string_t>(which, 0)), r_null);
  return x;
}

// Set attributes of x in-place, when add = F, attrs of x are first removed

[[cpp11::register]]
SEXP cpp_set_add_attributes(SEXP x, SEXP attributes, bool add) {

  if (!add) attr::clear_attrs(x);
  internal::add_attrs(x, attributes);
  return x;
}
