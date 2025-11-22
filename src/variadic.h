#ifndef CHEAPR_VARIADIC_H
#define CHEAPR_VARIADIC_H

#include <core.h>
#include "declarations.h"

template<typename... Args>
inline SEXP r_paste(SEXP sep, SEXP collapse, Args... args){
  SEXP objs = cheapr::SHIELD(cheapr::new_r_list(args...));
  SEXP out = cheapr::SHIELD(cpp_paste(objs, sep, collapse));
  cheapr::YIELD(2);
  return out;
}

// R vector constructor from C++ types and SEXP
template<typename... Args>
SEXP new_r_vec(Args... args) {
  SEXP out = cheapr::SHIELD(cheapr::new_r_list(args...));
  cheapr::SHIELD(out = cpp_c(out));
  cheapr::YIELD(2);
  return out;
}

#endif
