#ifndef CHEAPR_R_LIST_H
#define CHEAPR_R_LIST_H

#include <cheapr/internal/r_vector.h>
#include <cheapr/internal/r_coerce.h>

namespace cheapr {
// Named argument

struct arg {
  const char* name;
  SEXP value;

  explicit arg(const char* n) : name(n), value(R_NilValue) {}
  explicit arg(const char* n, SEXP v) : name(n), value(v) {}

  template<typename T>
  arg operator=(T v) const {
    return arg(name, as<r_sexp>(v)); 
  }
};

// Variadic list constructor
template<typename... Args>
inline r_vec<r_sexp> make_list(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return r_vec<r_sexp>(0);
  } else {

    auto out = SHIELD(r_vec<r_sexp>(n));

    // Are any args named?
    constexpr bool any_named = (is<Args, arg> || ...);

    auto nms = any_named ? r_vec<r_str>(n) : r_vec<r_str>();
    SHIELD(nms);

    int i = 0;
    (([&]() {
      if constexpr (is<Args, arg>) { 
        out.set(i, args.value);
        nms.set(i, as<r_str>(args.name));
      } else {
        out.set(i, as<r_sexp>(args));
      }
      ++i;
    }()), ...);

    attr::set_old_names(out, nms);
    YIELD(2);
    return out;
  }
}

}

#endif
