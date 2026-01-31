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

template<typename T, typename U>
requires (RVector<T> || RVector<U>)
inline auto operator-(const T& lhs, const U& rhs) {

    if constexpr (RVector<T> && RVector<U>){
        if (lhs.length() == 1){
            return lhs.get(0) - rhs;
        } else if (rhs.length() == 1){
            return lhs - rhs.get(0);
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
                out.set(i, lhs.get(lhsi) - rhs.get(rhsi));
            }
            return out;
        }
    } else if constexpr (RVector<T>){
        using common_t = common_r_math_t<typename T::data_type, U>;
        r_size_t n = lhs.length();
        r_vec<common_t> out(n);
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, lhs.get(i) - rhs);
        }
        return out;
    } else {
        using common_t = common_r_math_t<T, typename U::data_type>;
        r_size_t n = rhs.length();
        r_vec<common_t> out(n);
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, rhs.get(i) - lhs);
        }
        return out;
    }
}

template<typename T, typename U>
requires (RVector<T> || RVector<U>)
inline auto operator*(const T& lhs, const U& rhs) {

    if constexpr (RVector<T> && RVector<U>){
        if (lhs.length() == 1){
            return lhs.get(0) * rhs;
        } else if (rhs.length() == 1){
            return lhs * rhs.get(0);
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
                out.set(i, lhs.get(lhsi) * rhs.get(rhsi));
            }
            return out;
        }
    } else if constexpr (RVector<T>){
        using common_t = common_r_math_t<typename T::data_type, U>;
        r_size_t n = lhs.length();
        r_vec<common_t> out(n);
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, lhs.get(i) * rhs);
        }
        return out;
    } else {
        using common_t = common_r_math_t<T, typename U::data_type>;
        r_size_t n = rhs.length();
        r_vec<common_t> out(n);
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, rhs.get(i) * lhs);
        }
        return out;
    }
}

template<typename T, typename U>
requires (RVector<T> || RVector<U>)
inline auto operator/(const T& lhs, const U& rhs) {

    if constexpr (RVector<T> && RVector<U>){
        if (lhs.length() == 1){
            return lhs.get(0) / rhs;
        } else if (rhs.length() == 1){
            return lhs / rhs.get(0);
        } else {
            // Slower recycling approach
            r_size_t lhs_size = lhs.length();
            r_size_t rhs_size = rhs.length();

            r_size_t n = std::max(lhs_size, rhs_size);
            if (lhs_size == 0 || rhs_size == 0){
                n = 0;
            }
            r_vec<r_dbl> out(n);
            for (r_size_t i = 0, lhsi = 0, rhsi = 0; i < n; 
                recycle_index(lhsi, lhs_size),
                recycle_index(rhsi, rhs_size), 
                ++i){
                out.set(i, lhs.get(lhsi) / rhs.get(rhsi));
            }
            return out;
        }
    } else if constexpr (RVector<T>){
        r_size_t n = lhs.length();
        r_vec<r_dbl> out(n);
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, lhs.get(i) / rhs);
        }
        return out;
    } else {
        r_size_t n = rhs.length();
        r_vec<r_dbl> out(n);
        OMP_SIMD
        for (r_size_t i = 0; i < n; ++i){
            out.set(i, rhs.get(i) / lhs);
        }
        return out;
    }
}


template<RIntegerType T>
T gcd(const r_vec<T> &x, bool na_rm = false, T tol = r_limits<T>::tolerance()){
  if (tol < 0 || tol >= 1){
    abort("`tol` must be >= 0 and < 1");
  }
  r_size_t n = x.length();

  if (n == 0){
    return na_value<T>();
  }

  auto out = x.get(0);
  for (r_size_t i = 1; i < n; ++i) {
      out = gcd(out, x.get(i), na_rm);
      if (!na_rm && is_na(out)){
          break;
      } else if (out == 1){
          break;
      }
  }
  return out;
}
 
template<RMathType T>
T gcd(const r_vec<T> &x, bool na_rm = false, T tol = r_limits<T>::tolerance()){
  if (tol < 0 || tol >= 1){
    abort("`tol` must be >= 0 and < 1");
  }
  r_size_t n = x.length();

  if (n == 0){
    return na_value<T>();
  }

  auto out = x.get(0);
  for (r_size_t i = 1; i < n; ++i) {
      out = gcd(out, x.get(i), na_rm);
      if (!na_rm && is_na(out)){
          break;
      }
      if (out > T(0.0) && out < (tol + tol)){
        out = tol;
        break;
      }
  }
  return out;
}

template<RMathType T>
T lcm(const r_vec<T> &x, bool na_rm = false, T tol = r_limits<T>::tolerance()){
    if (tol < 0 || tol >= 1){
      Rf_error("tol must be >= 0 and < 1");
    }
    r_size_t n = x.length();
    if (n == 0){
        return na_value<T>();
    }

    // Initialise first value as lcm
    T out = x.get(0);

    for (R_xlen_t i = 1; i < n; ++i) {
    if (!na_rm && is_na(out)){
        break;
    } else if (is_r_pos_inf(out)){
        break;
    }
    out = lcm(out, x.get(i), na_rm, tol);
    }
    return out;
  }

}


#endif
