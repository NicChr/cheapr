#ifndef CHEAPR_R_EXPRS_H
#define CHEAPR_R_EXPRS_H

#include <cheapr/internal/r_vector.h>

namespace cheapr {

inline SEXP eval(SEXP expr, SEXP env){
  return Rf_eval(expr, env);
}

template<typename... Args>
inline SEXP make_pairlist(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return Rf_allocList(0); 
  } else {
    SEXP out = SHIELD(Rf_allocList(n));

    SEXP current = out;

    (([&]() {
      if constexpr (is<Args, arg>) {
        SETCAR(current, args.value);
        SET_TAG(current, as<r_sym>(args.name));
      } else {
        SETCAR(current, as<r_sexp>(args));
      }
      current = CDR(current);
    }()), ...);

    YIELD(1);
    return out;
  }
}


template<typename... Args>
inline SEXP make_expr(Args... args) { 

  SEXP pairlist = SHIELD(make_pairlist(args...));
  SEXP out = SHIELD(Rf_lcons(r_null, pairlist)); 

  YIELD(2);
  return out;
}

}

#endif
