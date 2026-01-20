#ifndef CHEAPR_R_ATTRS_H
#define CHEAPR_R_ATTRS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_symbols.h>
#include <cheapr/internal/r_vec.h>

namespace cheapr {

namespace attr {

template <RObject T>
inline bool can_have_attributes(T x){
  if constexpr (any<T, r_factors, r_df> || RVector<T>){
    return true; 
  } else if constexpr (is_sexp<T>){
    switch (TYPEOF(x)){
      case LGLSXP:
      case INTSXP: 
      case REALSXP: 
      case STRSXP: 
      case CPLXSXP: 
      case RAWSXP: 
      case VECSXP: 
      case LISTSXP:
      case ENVSXP:
      case CLOSXP:
      case BUILTINSXP: 
      case SPECIALSXP: 
      case LANGSXP:
      case EXPRSXP: 
      case EXTPTRSXP: 
      case OBJSXP: {
        return true;
      }
      default: {
        return false;
      }
    }
  } else {
    return false;
  }
}

template <RObject T>
inline bool can_have_names(T x){
  if constexpr (any<T, r_factors, r_df> || RVector<T>){
    return true; 
  } else if constexpr (is_sexp<T>){
    switch (TYPEOF(x)){
      case LGLSXP:
      case INTSXP: 
      case REALSXP: 
      case STRSXP: 
      case CPLXSXP: 
      case RAWSXP: 
      case VECSXP: 
      case LISTSXP:
      case ENVSXP:
      case LANGSXP:
      case EXPRSXP: {
        return true;
      }
      default: {
        return false;
      }
    }
  } else {
    return false;
  }
}

template <RObject T>
inline bool inherits1(T x, const char *r_cls){
  return Rf_inherits(x, r_cls);
}

// Attributes of x as a list
template <RObject T>
inline r_vec<r_sexp> get_attrs(T x){
  r_sexp attrs = r_sexp(ATTRIB(x)); // Pairlist
  r_size_t n = attrs.length();
  if (n == 0){
    return r_vec<r_sexp>(r_null);
  }
  r_vec<r_sexp> out(n);
  r_vec<r_str> names(n);

  SEXP current = attrs.value;

  for (r_size_t i = 0; i < n; ++i){
    out.set(i, r_sexp(CAR(current), internal::read_only_tag{}));
    names.set(i, internal::as_r<r_str>(symbol::tag(current)));
    current = CDR(current);
  }
  out.set_names(names);
  return out;
}
template <RObject T>
inline bool has_attrs(T x){
  return !get_attrs(x).is_null();
}
template <RObject T>
inline r_sexp get_attr(T x, r_sym sym){
  return r_sexp(Rf_getAttrib(x, sym));
}
template <RObject T, RObject U> 
inline void set_attr(T x, r_sym sym, U value){
  if (can_have_attributes(x)){
    Rf_setAttrib(x, sym, value); 
  }
}
template <RObject T>
inline void set_old_names(T x, r_vec<r_str> names){
  if (x.is_null()){
    return;
  } else if (names.is_null()){
    attr::set_attr(x, symbol::names_sym, r_null);
  } else if (names.length() != x.length()){
    cpp11::stop("`length(names)` must equal `length(x)`");
  } else if (can_have_names(x)){
    Rf_namesgets(x, names);
  } else {
    cpp11::stop("Cannot add names to unsupported type");
  }
}
template <RObject T>
inline r_vec<r_str> get_old_names(T x){
  return r_vec<r_str>(get_attr(x, symbol::names_sym));
}
template <RObject T>
inline bool has_r_names(T x){
  return !get_old_names(x).is_null();
}
template <RObject T>
inline r_vec<r_str> get_old_class(T x){
  return r_vec<r_str>(get_attr(x, symbol::class_sym));
}
template <RObject T>
inline void set_old_class(T x, r_vec<r_str> cls){
  if (can_have_attributes(x)){
    Rf_classgets(x, cls);
  }
}
template <RObject T>
inline bool inherits(T x, r_vec<r_str> classes){
  r_size_t n = classes.length();
  for (r_size_t i = 0; i < n; ++i) {
    if (inherits1(x, classes.get(i).c_str())){
      return true;
    }
  }
  return false;
}
template <RObject T>
inline void clear_attrs(T x){

  auto attrs = get_attrs(x);

  if (attrs.is_null()){
    return;
  }
  auto names = attr::get_old_names(attrs);

  int n = attrs.length();

  for (r_size_t i = 0; i < n; ++i){
    r_sym target_sym = internal::as_r<r_sym>(names.get(i));
    set_attr(x, target_sym, r_null);
  }
}

}

namespace internal {

template <RObject T>
inline void modify_attrs_impl(T x, r_vec<r_sexp> attrs) {

  if (x.is_null()){
    cpp11::stop("Cannot add attributes to `NULL`");
  }

  if (attrs.is_null()){
    return;
  }

  auto names = attr::get_old_names(attrs);

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

template <RObject T>
inline void set_attrs(T x, r_vec<r_sexp> attrs){
  clear_attrs(x);
  internal::modify_attrs_impl(x, attrs);
}

}

}

#endif

