#ifndef CHEAPR_R_MAKE_VEC_H
#define CHEAPR_R_MAKE_VEC_H


#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_coerce.h>
#include <variant>

namespace cheapr {

// Named argument
struct arg {
  const char* name;
  using storage_t = std::variant<
    std::monostate,
    r_lgl,
    r_int,
    r_int64,
    r_dbl,
    r_str,
    r_cplx,
    r_raw,
    r_sexp
  >;

  storage_t storage;

  // Default: no value yet
  explicit arg(const char* n)
    : name(n), storage(std::monostate{}) {}

  template <typename U>
  arg(const char* n, U&& v)
    : name(n), storage(internal::as_r_type(std::forward<U>(v))) {}

  template <typename U>
  arg operator=(U&& v) const {
    return arg{name, std::forward<U>(v)};
  }
};


template <RScalar T>
inline T as(const arg& a) {
  return std::visit(
    [](auto const& x) -> T {
      using X = std::remove_cvref_t<decltype(x)>;
      if constexpr (is<X, std::monostate>) {
        return T();
      } else {
        return as<T>(x);
      }
    },
    a.storage
  );
}


template<RScalar T, typename... Args>
inline r_vec<T> make_vec(Args... args) {

  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return r_vec<T>(0);
  } else {

    auto out = r_vec<T>(n);

    // Are any args named?
    constexpr bool any_named = (is<Args, arg> || ...);

    auto nms = any_named ? r_vec<r_str>(n) : r_vec<r_str>(r_null);

    int i = 0;
    (([&]() {
      if constexpr (is<Args, arg>) {
        out.set(i, as<T>(args));
        nms.set(i, as<r_str>(args.name));
      } else {
        out.set(i, as<T>(args));
      }
      ++i;
    }()), ...);

    attr::set_old_names(out.sexp, nms);
    return out;
  }
}

// template <RScalar T>
// template<typename... Args>
// void r_vec<T>::modify_attrs(Args... args) {

//   if (this->sexp.is_null()){
//     cpp11::stop("Cannot add attributes to `NULL`");
//   }


//   auto attrs = make_vec<r_sexp>(args...);
//   auto names = attrs.get_names();

//   if (attrs.).is_null(){
//     return;
//   }

//   if (names.).is_null(){
//     cpp11::stop("attributes must be a named list");
//   }

//   r_sym attr_nm;

//   int n = names.length();

//   for (int i = 0; i < n; ++i){
//     if (!(names.get(i) == blank_r_string)){
//       attr_nm = internal::as_r<r_sym>(names.get(i));
//       if (this->address() == internal::address(attrs.get(i))){
//         SEXP dup_attr = Rf_duplicate(attrs.get(i));
//         attr::set_attr(sexp, attr_nm, dup_attr);
//       } else {
//         attr::set_attr(sexp, attr_nm, attrs.get(i));
//       }
//     }
//   }
// }

}

#endif
