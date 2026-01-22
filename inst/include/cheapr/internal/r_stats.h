#ifndef CHEAPR_R_STATS_H
#define CHEAPR_R_STATS_H

#include <cheapr/internal/r_make_vec.h>

namespace cheapr {
    
template <RMathType T>
r_dbl sum(r_vec<T> x, bool na_rm = false){
    r_size_t n = x.length();
    r_dbl out(0);
    if (na_rm){
        for (r_size_t i = 0; i < n; ++i){ 
            if (is_na(x.get(i))){
                continue;
            } else {
                out += x.get(i);
            }
        }
    } else {
        for (r_size_t i = 0; i < n; ++i){
            out += x.get(i);
        }
    }
    return out;
}

// Optimisation for r_dbl
template <>
r_dbl sum(r_vec<r_dbl> x, bool na_rm){
    r_size_t n = x.length();
    r_dbl out(0);
    if (na_rm){
        for (r_size_t i = 0; i < n; ++i){ 
            if (is_na(x.get(i))){
                continue;
            } else {
                out += x.get(i);
            }
        }
    } else {
        double out_ = out;
        const auto* RESTRICT p_x = x.data();
        #pragma omp simd reduction(+:out_)
        for (r_size_t i = 0; i < n; ++i){
            out_ += p_x[i].value;
        }
        out = r_dbl(out_);
    }
    return out;
}

// Integer specific sum (user must accept there may be overflow)
template <RIntegerType T>
r_int64 sum_int(r_vec<T> x, bool na_rm = false){
    r_size_t n = x.length();
    r_int64 out(0);
    if (na_rm){
        for (r_size_t i = 0; i < n; ++i){ 
            if (is_na(x.get(i))){
                continue;
            } else {
                out += x.get(i);
            }
        }
    } else {
        for (r_size_t i = 0; i < n; ++i){
            out += x.get(i);
        }
    }
    return out;
}

template <RMathType T>
r_vec<T> range(r_vec<T> x, bool na_rm = false){
    
    r_size_t n = x.length();

    T lo = r_limits<T>::max();
    T hi = r_limits<T>::min();

    // Use underlying types for SIMD vectorisation
    auto lo_ = lo.value;
    auto hi_ = hi.value;
    const auto* RESTRICT p_x = x.data(); 

    if (na_rm){
    #pragma omp simd reduction(std::min:lo_) reduction(std::max:hi_)
    for (r_size_t i = 0; i < n; ++i){
        if (is_na(p_x[i].value)){
            continue;
        } else {
            lo_ = std::min(lo_, p_x[i].value); 
            hi_ = std::max(hi_, p_x[i].value); 
        }
    }
    lo = T(lo_);
    hi = T(hi_);
    } else {
    // Can't use SIMD, `cheapr::min/max` checks for NAs automatically
    for (r_size_t i = 0; i < n; ++i){
        lo = min(lo, x.get(i)); 
        hi = max(hi, x.get(i));
    }
    }
    return make_vec<T>(lo, hi);
}

// SIMD optimisation for integer types
template <RIntegerType T>
r_vec<T> range(r_vec<T> x, bool na_rm){
    
    r_size_t n = x.length();

    T lo = r_limits<T>::max();
    T hi = r_limits<T>::min();

    auto lo_ = lo.value;
    auto hi_ = hi.value;

    const auto* RESTRICT p_x = x.data(); 

    if (na_rm){
        #pragma omp simd reduction(std::min:lo_) reduction(std::max:hi_)
        for (r_size_t i = 0; i < n; ++i){
            if (is_na(p_x[i].value)){
                continue;
            } else {
                lo_ = std::min(lo_, p_x[i].value); 
                hi_ = std::max(hi_, p_x[i].value); 
            }
        }
        lo = T(lo_);
        hi = T(hi_);
    } else {
        #pragma omp simd reduction(std::min:lo_) reduction(std::max:hi_)
        for (r_size_t i = 0; i < n; ++i){
            lo_ = std::min(lo_, p_x[i].value); 
            hi_ = std::max(hi_, p_x[i].value);
        }
        lo = T(lo_);
        hi = T(hi_);

        // We use the fact that if there were NAs then min(x) would be NA
        // Only works for R's integer types
        bool has_nas = is_na(lo);

        if (has_nas){
            lo = na_value<T>();
            hi = na_value<T>();
        }
    }
    return make_vec<T>(lo, hi);
}

template <RMathType T>
T min(r_vec<T> x, bool na_rm = false){
    return range(x, na_rm).get(0);
}
template <RMathType T>
T max(r_vec<T> x, bool na_rm = false){
    return range(x, na_rm).get(1);
}

} 

#endif
