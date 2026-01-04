#ifndef CHEAPR_R_EXPRS_H
#define CHEAPR_R_EXPRS_H

#include <cheapr/internal/r_setup.h>

namespace cheapr {

inline SEXP eval(SEXP expr, SEXP env){
  return Rf_eval(expr, env);
}

template<typename... Args>
inline SEXP make_expr(Args... args) { 
  constexpr int n = sizeof...(args);

  SEXP pairlist = SHIELD(make_pairlist(args...));
  SEXP out = SHIELD(Rf_lcons(r_null, pairlist)); 

  YIELD(2);
  return out;
}

}

#endif
