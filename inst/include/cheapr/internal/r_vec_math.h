#ifndef CHEAPR_R_VEC_MATH_H
#define CHEAPR_R_VEC_MATH_H

#include <cheapr/internal/r_vec.h>

namespace cheapr {

    template<typename T, typename U>
    requires (RVector<T> || RVector<U>)
    inline auto operator+(const T& lhs, const U& rhs) {
    
        if constexpr (RVector<T> && RVector<U>){
            if (lhs.length() == 1){
                return lhs.get(0) + rhs;
            } else if (rhs.length() == 1){
                return lhs + rhs.get(0);
            } else {
                // Slower recycling approach
                r_size_t lhs_size = lhs.length();
                r_size_t rhs_size = rhs.length();
    
                r_size_t n = std::max(lhs_size, rhs_size);
                if (lhs_size == 0 || rhs_size == 0){
                    n = 0;
                }
                using common_t = common_r_math_t<typename T::data_type, typename U::data_type>;
                r_vec<common_t> out(n);
                for (r_size_t i = 0, lhsi = 0, rhsi = 0; i < n; 
                    recycle_index(lhsi, lhs_size),
                    recycle_index(rhsi, rhs_size), 
                    ++i){
                    out.set(i, lhs.get(lhsi) + rhs.get(rhsi));
                }
                return out;
            }
        } else if constexpr (RVector<T>){
            using common_t = common_r_math_t<typename T::data_type, U>;
            r_size_t n = lhs.length();
            r_vec<common_t> out(n);
            OMP_SIMD
            for (r_size_t i = 0; i < n; ++i){
                out.set(i, lhs.get(i) + rhs);
            }
            return out;
        } else {
            using common_t = common_r_math_t<T, typename U::data_type>;
            r_size_t n = rhs.length();
            r_vec<common_t> out(n);
            OMP_SIMD
            for (r_size_t i = 0; i < n; ++i){
                out.set(i, rhs.get(i) + lhs);
            }
            return out;
        }
    }

}


#endif
