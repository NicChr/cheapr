#ifndef CHEAPR_R_FNS_H
#define CHEAPR_R_FNS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_env.h>
#include <cheapr/internal/r_symbols.h>
#include <cheapr/internal/r_exprs.h>
#include <cheapr/internal/r_coerce.h>

namespace cheapr {

namespace fn {
// Return R function from a specified package
inline SEXP find_pkg_fun(const char *name, const char *pkg, bool all_fns){

  SEXP expr = r_null;

  if (all_fns){
    expr = SHIELD(Rf_lang3(symbol::triple_colon_sym, as<r_symbol_t>(pkg), as<r_symbol_t>(name)));
  } else {
    expr = SHIELD(Rf_lang3(symbol::double_colon_sym, as<r_symbol_t>(pkg), as<r_symbol_t>(name)));
  }
  SEXP out = SHIELD(eval(expr, env::base_env));
  YIELD(2);
  return out; 
}

template<typename... Args>
inline SEXP eval_fn(SEXP r_fn, SEXP envir, Args... args){
  // Expression
  SEXP expr = SHIELD(make_expr(r_fns, args...));
  // Evaluate expression
  SEXP out = SHIELD(eval(call, envir));
  YIELD(2);
  return out;
}

}

}

#endif
