#ifndef CHEAPR_VARIADIC_H
#define CHEAPR_VARIADIC_H

#include <core.h>
#include "decls.h"

template<typename... Args>
inline SEXP r_combine(Args... args){
  SEXP objs = cheapr::SHIELD(cheapr::new_r_list(args...));
  SEXP out = cheapr::SHIELD(cpp_c(objs));
  cheapr::YIELD(2);
  return out;
}

template<typename... Args>
inline SEXP r_paste(SEXP sep, SEXP collapse, Args... args){
  SEXP objs = cheapr::SHIELD(cheapr::new_r_list(args...));
  SEXP out = cheapr::SHIELD(cpp_paste(objs, sep, collapse));
  cheapr::YIELD(2);
  return out;
}

#endif
