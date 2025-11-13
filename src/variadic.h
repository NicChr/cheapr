#ifndef CHEAPR_VARIADIC_H
#define CHEAPR_VARIADIC_H

#include "cheapr_core.h"
#include "decls.h"

template<typename... Args>
inline SEXP make_r_list(Args... args){
  constexpr int n = sizeof...(args);
  SEXP out = SHIELD(new_vec(VECSXP, n));
  int i = 0;
  int dummy[] = {(SET_VECTOR_ELT(out, i++, args), 0)...};
  static_cast<void>(dummy);
  YIELD(1);
  return out;
}

// Make a character vec from const char ptrs
template<typename... Args>
inline SEXP make_r_chars(Args... args){
  constexpr int n = sizeof...(args);
  SEXP out = SHIELD(new_vec(STRSXP, n));
  int i = 0;
  int dummy[] = {(SET_STRING_ELT(out, i++, make_utf8_char(args)), 0)...};
  static_cast<void>(dummy);
  YIELD(1);
  return out;
}

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
