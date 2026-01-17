#ifndef CHEAPR_R_EXPRS_H
#define CHEAPR_R_EXPRS_H

#include <cheapr/internal/r_make_vec.h>

namespace cheapr {

inline r_sexp eval(r_sexp expr, r_sexp env){
  return r_sexp(cpp11::safe[Rf_eval](expr, env));
}

template<typename... Args>
inline r_sexp make_pairlist(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return r_sexp(Rf_allocList(0));
  } else {
    r_sexp out = r_sexp(cpp11::safe[Rf_allocList](n));

    SEXP current = out;

    (([&]() {
      if constexpr (is<Args, arg>) {
        SETCAR(current, as<r_sexp>(args.value));
        SET_TAG(current, as<r_sym>(args.name));
      } else {
        SETCAR(current, as<r_sexp>(args));
      }
      current = CDR(current);
    }()), ...);

    return out;
  }
}


template<typename... Args>
inline r_sexp make_expr(Args... args) { 

  r_sexp pairlist = make_pairlist(args...);
  return r_sexp(Rf_lcons(r_null, pairlist));
}

}

#endif
