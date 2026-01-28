#ifndef CHEAPR_R_STATS_H
#define CHEAPR_R_STATS_H

#include <cheapr/internal/r_make_vec.h>

namespace cheapr {
    
template <RMathType T>
r_dbl sum(r_vec<T> x, bool na_rm = false){
    r_size_t n = x.length();
    double out_ = 0;
    const auto* RESTRICT p_x = x.data();

    if (na_rm){
        #pragma omp simd reduction(+:out_)
        for (r_size_t i = 0; i < n; ++i){
            out_ += (is_na(x.get(i))) ? 0 : unwrap(p_x[i]);
        }
    } else {
        #pragma omp simd reduction(+:out_)
        for (r_size_t i = 0; i < n; ++i){
            out_ = (is_na(x.get(i)) || is_na(as_r_val(out_))) ? unwrap(na_value<r_dbl>()) : (out_ + unwrap(p_x[i]));
        }
    }
    return r_dbl(out_);
}

// Optimisation for r_dbl
template <>
r_dbl sum(r_vec<r_dbl> x, bool na_rm){
    r_size_t n = x.length();
    double out_ = 0;
    const auto* RESTRICT p_x = x.data();

    if (na_rm){
        #pragma omp simd reduction(+:out_)
        for (r_size_t i = 0; i < n; ++i){
            out_ += (is_na(x.get(i))) ? 0 : unwrap(p_x[i]);
        }
    } else {
        #pragma omp simd reduction(+:out_)
        for (r_size_t i = 0; i < n; ++i){
            out_ += unwrap(p_x[i]);
        }
        
    }
    return r_dbl(out_);
}

// Integer specific sum (user must accept there may be overflow)
template <RIntegerType T>
auto sum_int(r_vec<T> x, bool na_rm = false){
    r_size_t n = x.length();
    int64_t out_ = 0;
    const auto* RESTRICT p_x = x.data();

    if (na_rm){
        #pragma omp simd reduction(+:out_)
        for (r_size_t i = 0; i < n; ++i){
            out_ += (is_na(x.get(i))) ? int64_t(0) : unwrap(p_x[i]);
        }
    } else {
        r_int64 temp(out_);
        #pragma omp simd reduction(+:out_)
        for (r_size_t i = 0; i < n; ++i){
            if (is_na(x.get(i)) || is_na(temp)){
                temp = na_value<r_int64>();
                // Move underlying value of temp into out_ directly
                out_ = std::move(temp.value);
            } else {
                out_ += unwrap(p_x[i]);
            }
        }
    }
    return out_;
}

template <RMathType T>
r_vec<T> range(r_vec<T> x, bool na_rm = false){
    
    r_size_t n = x.length();

    T lo = r_limits<T>::max();
    T hi = r_limits<T>::min();

    // Can't use SIMD, `cheapr::min/max` checks for NAs automatically
    if (na_rm){
        for (r_size_t i = 0; i < n; ++i){
            const auto v = x.get(i);
            if (is_na(v)){
                continue;
            } else {
                lo = min(lo, v);
                hi = max(hi, v);
            }
        }
    } else {
        for (r_size_t i = 0; i < n; ++i){
            const auto v = x.get(i);
            lo = min(lo, v); 
            hi = max(hi, v);
        }
    }
    return make_vec<T>(lo, hi);
}

// SIMD optimisation for integer types
template <RIntegerType T>
r_vec<T> range(r_vec<T> x, bool na_rm){
    
    r_size_t n = x.length();

    T max_val = r_limits<T>::max();
    T min_val = r_limits<T>::min();

    T lo = max_val;
    T hi = min_val;

    auto lo_ = unwrap(lo);
    auto hi_ = unwrap(hi);

    const auto* RESTRICT p_x = x.data();

    if (na_rm){ 
        #pragma omp simd reduction(std::min:lo_) reduction(std::max:hi_)
        for (r_size_t i = 0; i < n; ++i){
            // Ignore NA for min()
            lo_ = is_na(x.get(i)) ? lo_ : std::min(lo_, unwrap(p_x[i]));
            // No need to ignore NA for max() because NA is defined as lowest representable value
            hi_ = std::max(hi_, unwrap(p_x[i]));
        }        
        lo = T(lo_);
        hi = T(hi_);
    } else {
        #pragma omp simd reduction(std::min:lo_) reduction(std::max:hi_)
        for (r_size_t i = 0; i < n; ++i){
            lo_ = std::min(lo_, unwrap(p_x[i])); 
            hi_ = std::max(hi_, unwrap(p_x[i]));
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

template <RMathType T>
r_vec<T> abs(r_vec<T> x){
    r_size_t n = x.length();
    r_vec<T> out(n);
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
        OMP_PARALLEL_FOR_SIMD(n_threads)
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, abs(x.get(i)));
        }
    } else {
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, abs(x.get(i)));
        }
    }
    return out;
}

} 

#endif
