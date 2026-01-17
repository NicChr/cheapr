#ifndef CHEAPR_R_MATH_H
#define CHEAPR_R_MATH_H

#include <cheapr/internal/r_methods.h>
#include <cheapr/internal/r_limits.h>
#include <cheapr/internal/r_coerce.h> 

// R math fns
// Propagate NA values correctly

namespace cheapr {

namespace internal {

inline r_dbl round_to_even(r_dbl x){
  return x - r_dbl{std::remainder(x.value, 1.0)};
}

}

namespace math {
template <typename T>
inline constexpr bool is_r_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_inf<r_dbl>(const r_dbl x){
  return x.value == r_limits::r_pos_inf.value || x.value == r_limits::r_neg_inf.value;
}

template <typename T>
inline constexpr bool is_r_pos_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_pos_inf<r_dbl>(const r_dbl x){
  return x.value == r_limits::r_pos_inf.value;
}

template <typename T>
inline constexpr bool is_r_neg_inf(const T x){
  return false;
}

template <>
inline constexpr bool is_r_neg_inf<r_dbl>(const r_dbl x){
  return x.value == r_limits::r_neg_inf.value;
}

template<RMathType T>
inline constexpr T abs(T x){
  return is_r_na(x) ? x : T{std::abs(x.value)};
}

template<RMathType T>
inline T floor(T x){
  return is_r_na(x) ? x : T{std::floor(x.value)};
}
template<>
inline r_dbl floor(r_dbl x){
  return r_dbl(std::floor(x.value));
}

template<RMathType T>
inline T ceiling(T x){
  return is_r_na(x) ? x : T{std::ceil(x.value)};
}
template<>
inline r_dbl ceiling(r_dbl x){
  return r_dbl(std::ceil(x.value));
}

template<RMathType T>
inline T trunc(T x){
  return is_r_na(x) ? x : T{std::trunc(x.value)};
}

template <>
inline r_dbl trunc(r_dbl x){
  return r_dbl(std::trunc(x.value) + 0.0);
}

template <RMathType T>
inline r_int sign(T x) {
  return is_r_na(x) ? na::integer : (T(0) < x) - (x < T(0));
}

template<RMathType T>
inline T negate(T x){
  return -x;
}

template<RMathType T>
inline r_dbl sqrt(T x){
  return r_dbl(std::sqrt(as<r_dbl>(x).value));
}

template<typename T, typename U>
  requires (AtLeastOneRMathType<T, U>)
inline r_dbl pow(T x, U y){
  if (y == 0) return r_dbl(1.0);
  if (x == 1) return r_dbl(1.0);
  if (y == 2){
    r_dbl left = as<r_dbl>(x);
    return left * left;
  }
  return r_dbl(std::pow(as<r_dbl>(x), as<r_dbl>(y)));
}

template<RMathType T>
inline r_dbl log10(T x){
  return r_dbl(std::log10(as<r_dbl>(x).value));
}

template<RMathType T>
inline r_dbl exp(T x){
  return r_dbl(std::exp(as<r_dbl>(x).value));
}

template<typename T, typename U>
requires (AtLeastOneRMathType<T, U>)
inline r_dbl log(T x, U base){
  return r_dbl(std::log(as<r_dbl>(x)) / std::log(as<r_dbl>(base)));
}
template<RMathType T>
inline r_dbl log(T x){
  return r_dbl(std::log(as<r_dbl>(x).value));
}
inline r_cplx log(r_cplx x){
  if (is_r_na(x)){
    return x;
  }
  r_dbl real = as<r_dbl>(0.5 * (log(pow(x.re(), 2.0) + pow(x.im(), 2.0))));
  r_dbl imag = as<r_dbl>(std::atan2(as<r_dbl>(x.im()), as<r_dbl>(x.re())));
  return r_cplx{real, imag};
}


template<typename T, typename U>
requires (AtLeastOneRMathType<T, U>)
inline r_dbl round(T x, const U digits){
  if (is_r_na(x)){
    return as<r_dbl>(x);
  } else if (is_r_na(digits)){
    return na::real;
  } else if (is_r_inf(x)){
    return x;
  } else if (is_r_neg_inf(digits)){
    return r_dbl(0.0);
  } else if (is_r_pos_inf(digits)){
    return x;
  } else {
    r_dbl scale = std::pow(10.0, as<r_dbl>(digits));
    return internal::round_to_even(as<r_dbl>(x) * scale) / scale;
  }
}

template<RMathType T>
inline r_dbl round(T x){
  if (is_r_na(x)){
    return as<r_dbl>(x);
  } else if (is_r_inf(x)){
    return x;
  } else {
    return internal::round_to_even(as<r_dbl>(x));
  }
}

template<typename T, typename U>
requires (AtLeastOneRMathType<T, U>)
inline r_dbl signif(T x, const U digits){
  U new_digits = std::max(U(1), digits);
  if (is_r_na(x)){
    return as<r_dbl>(x);
  } else if (is_r_na(new_digits)){
    return na::real;
  } else if (is_r_pos_inf(digits)){
    return x;
  } else {
    new_digits -= std::ceil(std::log10(std::abs(x)));
    r_dbl scale = std::pow(10, new_digits);
    return internal::round_to_even(scale * x) / scale;
  }
}

template<MathType T, MathType U>
requires (AtLeastOneRMathType<T, U>)
inline T abs_diff(const T x, const U y){
  return abs(x - y);
}

inline r_lgl is_whole_number(const r_dbl x, const r_dbl tolerance){
  return is_r_na(x) || is_r_na(tolerance) ? na::logical : r_lgl(abs_diff(x, round(x)) <= tolerance);
}


// Greatest common divisor
template<MathType T>
  inline T gcd(T x, T y, bool na_rm = true, T tol = std::sqrt(std::numeric_limits<T>::epsilon())){

    if (is_r_na(x) || is_r_na(y)){
      if (na_rm){ 
        if (is_r_na(x)){
          return abs(y);
        } else {
          return abs(x);
        }
      } else {
        return na_value<decltype(x)>();
      }
    }

    T ax = std::abs(x);
    T ay = std::abs(y);

    if constexpr (IntegerType<T>){

      // Taken from number theory lecture notes

      // GCD(0,0)=0
      if (ax == 0 && ay == 0){
        return 0;
      }
      // GCD(a,0)=a
      if (ax == 0){
        return ay;
      }
      // GCD(a,0)=a
      if (ay == 0){
        return ax;
      }

      T r;
      while(ay != 0){
        r = ax % ay;
        ax = ay;
        ay = r;
      }
      return ax;
    } else {

      // GCD(0,0)=0
      if (ax <= tol && ay <= tol){
        return 0.0;
      }
      // GCD(a,0)=a
      if (ax <= tol){
        return ay;
      }
      // GCD(a,0)=a
      if (ay <= tol){
        return ax;
      }

      T r;
      while(ay > tol){
        r = std::fmod(ax, ay);
        ax = ay;
        ay = r;
      }
      return ax;
    }
  }


// Lowest common multiple
template<MathType T>
  inline T lcm(
      T x, T y, bool na_rm = true, T tol = std::sqrt(std::numeric_limits<T>::epsilon())
  ){
    if (is_r_na(x) || is_r_na(y)){
      if (na_rm){
        if (is_r_na(x)){
          return y;
        } else {
          return x;
        }
      } else {
        return na_value<decltype(x)>();
      }
    }

    if constexpr (IntegerType<T>){
      if (x == 0 && y == 0){
        return 0;
      }
      T res = std::abs(x) / gcd(x, y, na_rm);
      if (y != 0 && (std::abs(res) > (std::numeric_limits<T>::max() / std::abs(y)))){
        return na_value<decltype(x)>();
      }
      return res * std::abs(y);
    } else {
      if (std::fabs(x) <= tol && std::fabs(y) <= tol){
        return 0.0;
      }
      return ( std::fabs(x) / gcd(x, y, na_rm, tol) ) * std::fabs(y);
    }
  }

}

}

#endif
