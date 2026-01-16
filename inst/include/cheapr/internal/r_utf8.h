#ifndef CHEAPR_R_UTF8_H
#define CHEAPR_R_UTF8_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>

namespace cheapr {

namespace internal {
// UTF-8 helpers

inline const char* utf8_char(r_str x){
  return Rf_translateCharUTF8(static_cast<r_sexp>(x));
}

// inline r_str make_utf8_charsxp(const char *x){
//   return r_str(Rf_mkCharCE(x, CE_UTF8));
// }


inline const char* char_as_utf8(const char *x){
  return r_str(x).c_str();
}

}

}

#endif
