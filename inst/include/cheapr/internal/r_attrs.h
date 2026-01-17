#ifndef CHEAPR_R_ATTRS_H
#define CHEAPR_R_ATTRS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_symbols.h>
#include <cheapr/internal/r_vector.h>

namespace cheapr {

namespace attr {

inline bool inherits1(r_sexp x, const char *r_cls){
  return Rf_inherits(x, r_cls);
}

// Attributes of x as a list
inline r_vec<r_sexp> get_attrs(r_sexp x){
  if (internal::BASE_ATTRIBUTES == NULL){
    internal::BASE_ATTRIBUTES = Rf_install("attributes");
  }
  r_sexp expr = r_sexp(Rf_lang2(internal::BASE_ATTRIBUTES, x));
  r_sexp out = r_sexp(Rf_eval(expr, R_BaseEnv));
  return r_vec<r_sexp>(out);
}

inline bool has_attrs(r_sexp x){
  return !get_attrs(x).is_null();
}

inline r_sexp get_attr(r_sexp x, r_sym sym){
  return r_sexp(Rf_getAttrib(x, sym));
}

inline void set_attr(r_sexp x, r_sym sym, r_sexp value){
  if (!x.is_null()){
    Rf_setAttrib(x, sym, value);
  }
}

inline void set_old_names(r_sexp x, r_vec<r_str> names){
  if (names.is_null()){
    attr::set_attr(x, symbol::names_sym, r_null);
  } else if (names.length() != x.length()){
    cpp11::stop("`length(names)` must equal `length(x)`");
  } else {
    Rf_namesgets(x, names);
  }
}
inline r_vec<r_str> get_old_names(r_sexp x){
  return r_vec<r_str>(get_attr(x, symbol::names_sym));
}
inline bool has_r_names(r_sexp x){
  auto names = get_old_names(x);
  return !names.is_null();
}

inline r_vec<r_str> get_old_class(r_sexp x){
  return r_vec<r_str>(get_attr(x, symbol::class_sym));
}
inline void set_old_class(r_sexp x, r_vec<r_str> cls){
  r_vec<r_str>(Rf_classgets(x, cls));
}

inline bool inherits(r_sexp x, r_vec<r_str> classes){
  r_size_t n = classes.length();
  for (r_size_t i = 0; i < n; ++i) {
    if (inherits1(x, classes.get(i).c_str())){
      return true;
    }
  }
  return false;
}

inline void clear_attrs(r_sexp x){

  auto attrs = get_attrs(x);

  if (attrs.is_null()){
    return;
  }
  auto names = attr::get_old_names(attrs.sexp);

  int n = attrs.length();

  for (r_size_t i = 0; i < n; ++i){
    r_sym target_sym = internal::as_r<r_sym>(names.get(i));
    set_attr(x, target_sym, r_null);
  }
}

}

namespace internal {

inline void modify_attrs_impl(r_sexp x, cheapr::r_vec<r_sexp> attrs) {

  if (x.is_null()){
    cpp11::stop("Cannot add attributes to `NULL`");
  }

  if (attrs.is_null()){
    return;
  }

  auto names = attr::get_old_names(attrs.sexp);

  if (names.is_null()){
    cpp11::stop("attributes must be a named list");
  }

  r_sym attr_nm;

  int n = names.length();

  for (int i = 0; i < n; ++i){
    if (!(names.get(i) == blank_r_string)){
      attr_nm = internal::as_r<r_sym>(names.get(i));
      // We can't add an object as its own attribute in-place (as this will crash R)
      if (x.address() == attrs.get(i).address()){
        r_sexp dup_attr = r_sexp(Rf_duplicate(attrs.get(i)));
        attr::set_attr(x, attr_nm, dup_attr);
      } else {
        attr::set_attr(x, attr_nm, attrs.get(i));
      }
    }
  }
}

}

namespace attr {

inline void set_attrs(cheapr::r_sexp x, r_vec<r_sexp> attrs){
  if (x.is_null()){
    clear_attrs(x);
    internal::modify_attrs_impl(x, attrs);
  }
}

}

}

#endif

