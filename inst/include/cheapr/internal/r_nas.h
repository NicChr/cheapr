#ifndef CHEAPR_R_NAS
#define CHEAPR_R_NAS

#include <limits>
#include <r_types.h>

namespace cheapr {
    
// NAs

namespace na {
inline constexpr r_bool_t logical = r_na;
inline constexpr r_int_t integer = r_int_t{std::numeric_limits<int>::min()};
inline constexpr r_int64_t integer64 = r_int64_t(std::numeric_limits<int64_t>::min());
inline const r_double_t real = r_double_t{NA_REAL};
inline const r_complex_t complex = r_complex_t{real, real};
inline constexpr r_byte_t raw = r_byte_t{0};
inline const r_string_t string = r_string_t{NA_STRING};
inline const sexp_t nil = r_null;
}

}

#endif
