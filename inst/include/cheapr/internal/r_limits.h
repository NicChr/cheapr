#ifndef CHEAPR_R_LIMITS_H
#define CHEAPR_R_LIMITS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <limits>

namespace cheapr {
    
namespace r_limits { 
inline constexpr r_int_t r_int_min = r_int_t{-std::numeric_limits<int>::max()};
inline constexpr r_int_t r_int_max = r_int_t{std::numeric_limits<int>::max()};
inline constexpr r_int64_t r_int64_min = r_int64_t(-std::numeric_limits<int64_t>::max());
inline constexpr r_int64_t r_int64_max = r_int64_t(std::numeric_limits<int64_t>::max());
inline constexpr r_double_t r_pos_inf = r_double_t(std::numeric_limits<double>::infinity());
inline constexpr r_double_t r_neg_inf = r_double_t(-std::numeric_limits<double>::infinity());
}

}

#endif
