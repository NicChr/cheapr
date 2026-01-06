#ifndef CHEAPR_R_LIMITS_H
#define CHEAPR_R_LIMITS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <limits>

namespace cheapr {
    
namespace r_limits { 
inline constexpr int r_int_min = r_int{-std::numeric_limits<int>::max()};
inline constexpr int r_int_max = r_int{std::numeric_limits<int>::max()};
inline constexpr int64_t r_int64_min = r_int64(-std::numeric_limits<int64_t>::max());
inline constexpr int64_t r_int64_max = r_int64(std::numeric_limits<int64_t>::max());
inline constexpr double r_pos_inf = r_dbl(std::numeric_limits<double>::infinity());
inline constexpr double r_neg_inf = r_dbl(-std::numeric_limits<double>::infinity());
}

}

#endif
