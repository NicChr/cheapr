#ifndef CHEAPR_R_LIMITS_H
#define CHEAPR_R_LIMITS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <limits>

namespace cheapr {
    

template <RMathType T>
struct r_limits;

template <>
struct r_limits<r_lgl>{

    static constexpr r_lgl min() noexcept {
        return r_false;
    }
    static constexpr r_lgl max() noexcept {
        return r_true;
    }
    static constexpr r_lgl epsilon() noexcept {
        return r_false;
    }
    static constexpr r_lgl tolerance() noexcept {
        return r_false;
    }
};
template <>
struct r_limits<r_int>{
    static constexpr r_int min() noexcept {
        return r_int(-std::numeric_limits<int>::max());
    }
    static constexpr r_int max() noexcept {
        return r_int(std::numeric_limits<int>::max());
    }
    static constexpr r_int epsilon() noexcept {
        return r_int(0);
    }
    static constexpr r_int tolerance() noexcept {
        return r_int(0);
    }
};
template <>
struct r_limits<r_int64>{
    static constexpr r_int64 min() noexcept {
        return r_int64(-std::numeric_limits<int64_t>::max());
    }
    static constexpr r_int64 max() noexcept {
        return r_int64(std::numeric_limits<int64_t>::max());
    }
    static constexpr r_int64 epsilon() noexcept {
        return r_int64(0);
    }
    static constexpr r_int64 tolerance() noexcept {
        return r_int64(0);
    }
};
template <>
struct r_limits<r_dbl>{
    static constexpr r_dbl min() noexcept {
        return r_dbl(-std::numeric_limits<double>::infinity());
    }
    static constexpr r_dbl max() noexcept {
        return r_dbl(std::numeric_limits<double>::infinity());
    }
    static constexpr r_dbl epsilon() noexcept {
        return r_dbl(std::numeric_limits<double>::epsilon());
    }
    static constexpr r_dbl tolerance() noexcept {
        return r_dbl(std::sqrt(std::numeric_limits<double>::epsilon()));
    }
};

}

#endif
