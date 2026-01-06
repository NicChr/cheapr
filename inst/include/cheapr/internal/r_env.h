#ifndef CHEAPR_R_ENV_H
#define CHEAPR_R_ENV_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_exprs.h>

namespace cheapr {

namespace env {

inline const SEXP empty_env = R_EmptyEnv;
inline const SEXP base_env = R_BaseEnv;

inline SEXP get(r_sym sym, SEXP env, bool inherits = true){

  if (TYPEOF(env) != ENVSXP){
    Rf_error("second argument to '%s' must be an environment", __func__);
  }

  SEXP val = inherits ? Rf_findVar(sym, env) : Rf_findVarInFrame(env, sym);

  if (val == static_cast<SEXP>(symbol::missing_arg)){
    Rf_error("arg `sym` cannot be missing");
  } else if (val == static_cast<SEXP>(symbol::unbound_value)){
    return r_null;
  } else if (TYPEOF(val) == PROMSXP){
    SHIELD(val);
    val = eval(val, env);
    YIELD(1);
  }
  return val;
}
}

}

#endif
