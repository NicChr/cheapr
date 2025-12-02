#ifndef CHEAPR_VARIADIC_H
#define CHEAPR_VARIADIC_H

#include <core.h>
#include "declarations.h"

namespace cheapr {

template<typename... Args>
inline SEXP r_paste(SEXP sep, SEXP collapse, Args... args){
  SEXP objs = SHIELD(cheapr::new_r_list(args...));
  SEXP out = SHIELD(cpp_paste(objs, sep, collapse));
  YIELD(2);
  return out;
}

// Powerful R vector constructor from C++ types and SEXP
template<typename... Args>
SEXP new_r_vec(Args... args) {
  SEXP out = SHIELD(cheapr::new_r_list(args...));
  SHIELD(out = cpp_c(out));
  YIELD(2);
  return out;
}

}

#endif
