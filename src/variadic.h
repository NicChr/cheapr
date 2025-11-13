#ifndef CHEAPR_VARIADIC_H
#define CHEAPR_VARIADIC_H

#include "cheapr_core.h"
#include "decls.h"

template<typename... Args>
inline SEXP r_combine(Args... args){
  SEXP objs = SHIELD(make_r_list(args...));
  SEXP out = SHIELD(cpp_c(objs));
  YIELD(2);
  return out;
}

template<typename... Args>
inline SEXP r_paste(SEXP sep, SEXP collapse, Args... args){
  SEXP objs = SHIELD(make_r_list(args...));
  SEXP out = SHIELD(cpp_paste(objs, sep, collapse));
  YIELD(2);
  return out;
}

#endif
