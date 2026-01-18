#ifndef CHEAPR_R_FNS_H
#define CHEAPR_R_FNS_H

#include <cheapr/internal/r_env.h>
#include <cheapr/internal/r_symbols.h>
#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_exprs.h>

namespace cheapr {

namespace fn {
// Return R function from a specified package
inline r_sexp find_pkg_fun(const char *name, const char *pkg, bool all_fns){

  r_sexp expr = r_null;

  if (all_fns){
    expr = make_expr(symbol::triple_colon_sym, as<r_sym>(pkg), as<r_sym>(name));
  } else {
    expr = make_expr(symbol::double_colon_sym, as<r_sym>(pkg), as<r_sym>(name));
  }
  return eval(expr, env::base_env);
}

template<typename... Args>
inline r_sexp eval_fn(r_sexp r_fn, r_sexp envir, Args... args){
  // Expression
  r_sexp expr = make_expr(r_fn, args...);
  // Evaluate expression
  return eval(expr, envir);
}

}

}

#endif
