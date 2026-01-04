#ifndef CHEAPR_R_UTF8_H
#define CHEAPR_R_UTF8_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>

namespace cheapr {

namespace internal {
// UTF-8 helpers

inline const char* utf8_char(r_string_t x){
  return Rf_translateCharUTF8(static_cast<SEXP>(x));
}

inline SEXP make_utf8_charsxp(const char *x){
  return Rf_mkCharCE(x, CE_UTF8);
}

inline SEXP make_utf8_strsxp(const char *x){
  return Rf_ScalarString(make_utf8_charsxp(x));
}

inline const char* char_as_utf8(const char *x){
  return CHAR(make_utf8_charsxp(x));
}

// Memory address
inline r_string_t address(SEXP x) {
  char buf[1000];
  std::snprintf(buf, 1000, "%p", static_cast<void*>(x));
  return static_cast<r_string_t>(internal::make_utf8_charsxp(buf));
}

}

}

#endif
