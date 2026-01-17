#ifndef CHEAPR_R_ENV_H
#define CHEAPR_R_ENV_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_exprs.h>

namespace cheapr {

namespace env {

inline const r_sexp empty_env = r_sexp(R_EmptyEnv);
inline const r_sexp base_env = r_sexp(R_BaseEnv);

inline r_sexp get(r_sym sym, r_sexp env, bool inherits = true){

  if (TYPEOF(env) != ENVSXP){
    cpp11::stop("second argument to '%s' must be an environment", __func__);
  }

  r_sexp val = r_sexp(inherits ? Rf_findVar(sym, env) : Rf_findVarInFrame(env, sym));

  if (val == static_cast<SEXP>(symbol::missing_arg)){
    cpp11::stop("arg `sym` cannot be missing");
  } else if (val == static_cast<SEXP>(symbol::unbound_value)){
    return r_null;
  } else if (TYPEOF(val) == PROMSXP){
    val = eval(val, env);
  }
  return val;
}
}

}

#endif
