#ifndef CHEAPR_VARIADIC_H
#define CHEAPR_VARIADIC_H

#include <cheapr/internal/c_core.h>
#include "declarations.h"

namespace cheapr {

namespace vec {

template<typename... Args>
inline SEXP r_paste(SEXP sep, SEXP collapse, Args... args){
  SEXP objs = SHIELD(cheapr::vec::make_list(args...));
  SEXP out = SHIELD(cpp_paste(objs, sep, collapse));
  YIELD(2);
  return out;
}

// Powerful R vector constructor from C++ types and SEXP
template<typename... Args>
SEXP combine(Args... args) {
  SEXP out = SHIELD(cheapr::vec::make_list(args...));
  SHIELD(out = cpp_c(out));
  YIELD(2);
  return out;
}

}
}

#endif
