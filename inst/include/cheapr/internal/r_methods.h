#ifndef CHEAPR_R_METHODS_H
#define CHEAPR_R_METHODS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_nas.h>

namespace cheapr {

// Methods for custom R types

// operators for r_lgl
inline constexpr r_lgl operator!(r_lgl x) {
  return is_na(x) ? r_na : r_lgl(x.value == 0);
}

// r_true = 1, r_false = 0, r_na = INT_MIN

// ---------------------------------------------------------
// OPTIMIZED OR (||) for r_lgl
// If LSB is set (1), return 1
// Otherwise return (a|b).
// ---------------------------------------------------------
inline constexpr r_lgl operator||(r_lgl lhs, r_lgl rhs) {
    int val = lhs.value | rhs.value;
    return (val & 1) ? r_true : r_lgl(val);
}

// ---------------------------------------------------------
// OPTIMIZED AND (&&) for r_lgl
// If either is 0, return 0.
// if either is NA (negative), return NA.
// otherwise return 1.
// ---------------------------------------------------------
inline constexpr r_lgl operator&&(r_lgl lhs, r_lgl rhs) {
    if (lhs.value == 0 || rhs.value == 0) {
        return r_false;
    }
    return (lhs.value | rhs.value) < 0 ? r_na : r_true;
}

template<RVal T, RVal U>
inline constexpr r_lgl operator==(const T &lhs, const U &rhs) {
  return (is_na(lhs) || is_na(rhs)) ? r_na : r_lgl{unwrap(lhs) == unwrap(rhs)};
}

template<RVal T, CppScalar U>
inline constexpr r_lgl operator==(const T &lhs, const U &rhs) {
  return is_na(lhs) ? r_na : r_lgl{unwrap(lhs) == rhs};
}

template<CppScalar T, RVal U>
inline constexpr r_lgl operator==(const T &lhs, const U &rhs) {
  return is_na(rhs) ? r_na : r_lgl{lhs == unwrap(rhs)};
}

// Other comparison operators
template<RVal T, CppScalar U>
inline constexpr r_lgl operator!=(const T &lhs, const U &rhs) {
  r_lgl eq = (lhs == rhs);
  return eq.is_na() ? r_na : r_lgl(eq.is_false());
}
template<CppScalar T, RVal U>
inline constexpr r_lgl operator!=(const T &lhs, const U &rhs) {
  r_lgl eq = (lhs == rhs);
  return eq.is_na() ? r_na : r_lgl(eq.is_false());
}
template<RVal T, RVal U>
inline constexpr r_lgl operator!=(const T &lhs, const U &rhs) {
  r_lgl eq = (lhs == rhs);
  return eq.is_na() ? r_na : r_lgl(eq.is_false());
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator<(T lhs, U rhs) {
  return (is_na(lhs) || is_na(rhs)) ? r_na : r_lgl{lhs.value < unwrap(rhs)};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator<(T lhs, U rhs) {
  return is_na(lhs) ? r_na : r_lgl{lhs.value < rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator<(T lhs, U rhs) {
  return is_na(rhs) ? r_na : r_lgl{lhs < unwrap(rhs)};
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator<=(T lhs, U rhs) {
  return (is_na(lhs) || is_na(rhs)) ? r_na : r_lgl{lhs.value <= unwrap(rhs)};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator<=(T lhs, U rhs) {
  return is_na(lhs) ? r_na : r_lgl{lhs.value <= rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator<=(T lhs, U rhs) {
  return is_na(rhs) ? r_na : r_lgl{lhs <= unwrap(rhs)};
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator>(T lhs, U rhs) {
  return (is_na(lhs) || is_na(rhs)) ? r_na : r_lgl{lhs.value > unwrap(rhs)};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator>(T lhs, U rhs) {
  return is_na(lhs) ? r_na : r_lgl{lhs.value > rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator>(T lhs, U rhs) {
  return is_na(rhs) ? r_na : r_lgl{lhs > unwrap(rhs)};
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator>=(T lhs, U rhs) {
  return (is_na(lhs) || is_na(rhs)) ? r_na : r_lgl{lhs.value >= unwrap(rhs)};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator>=(T lhs, U rhs) {
  return is_na(lhs) ? r_na : r_lgl{lhs.value >= rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator>=(T lhs, U rhs) {
  return is_na(rhs) ? r_na : r_lgl{lhs >= unwrap(rhs)};
}

template<RMathType T, RMathType U>
inline constexpr T& operator+=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value += unwrap(rhs);
  }
  return lhs;
}
template<RMathType T, CppMathType U>
inline constexpr T& operator+=(T &lhs, U rhs) {
  if (is_na(lhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value += rhs;
  }
  return lhs;
}
template<CppMathType T, RMathType U>
inline constexpr T& operator+=(T &lhs, U rhs) {
  if (is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs += unwrap(rhs);
  }
  return lhs;
}

// Fast specialisation for r_dbl
template<>
inline constexpr r_dbl& operator+=(r_dbl &lhs, r_dbl rhs) {
  lhs.value += unwrap(rhs);
  return lhs;
}

template<MathType T, MathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr auto operator+(T lhs, U rhs) {

  using common_math_t = common_r_math_t<T, U>;

  if constexpr (is<common_math_t, r_dbl>){
    return r_dbl(static_cast<double>(unwrap(lhs)) + static_cast<double>(unwrap(rhs)));
  } else {
    return ( is_na(lhs) || is_na(rhs) ) ? 
    na_value<common_math_t>() : 
    common_math_t(static_cast<unwrap_t<common_math_t>>(unwrap(lhs)) + static_cast<unwrap_t<common_math_t>>(unwrap(rhs)));
  }
}

template<MathType T, MathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr auto operator-(T lhs, U rhs) {

  using common_math_t = common_r_math_t<T, U>;

  if constexpr (is<common_math_t, r_dbl>){
    return r_dbl(static_cast<double>(unwrap(lhs)) - static_cast<double>(unwrap(rhs)));
  } else {
    return ( is_na(lhs) || is_na(rhs) ) ? 
    na_value<common_math_t>() : 
    common_math_t(static_cast<unwrap_t<common_math_t>>(unwrap(lhs)) - static_cast<unwrap_t<common_math_t>>(unwrap(rhs)));
  }
}

template<MathType T, MathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr auto operator*(T lhs, U rhs) {

  using common_math_t = common_r_math_t<T, U>;

  if constexpr (is<common_math_t, r_dbl>){
    return r_dbl(static_cast<double>(unwrap(lhs)) * static_cast<double>(unwrap(rhs)));
  } else {
    return ( is_na(lhs) || is_na(rhs) ) ? 
    na_value<common_math_t>() : 
    common_math_t(static_cast<unwrap_t<common_math_t>>(unwrap(lhs)) * static_cast<unwrap_t<common_math_t>>(unwrap(rhs)));
  }
}

template<MathType T, MathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr r_dbl operator/(T lhs, U rhs) {
  return ( is_na(lhs) || is_na(rhs) ) ? na_value<r_dbl>() : r_dbl(static_cast<double>(unwrap(lhs)) / static_cast<double>(unwrap(rhs)));
}

template<RMathType T, RMathType U>
inline constexpr T& operator-=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value -= unwrap(rhs);
  }
  return lhs;
}

template<RMathType T, CppMathType U>
inline constexpr T& operator-=(T &lhs, U rhs) {
  if (is_na(lhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value -= rhs;
  }
  return lhs;
}

template<CppMathType T, RMathType U>
inline constexpr T& operator-=(T &lhs, U rhs) {
  if (is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs -= unwrap(rhs);
  }
  return lhs;
}
template<>
inline constexpr r_dbl& operator-=(r_dbl &lhs, r_dbl rhs) {
  lhs.value -= unwrap(rhs);
  return lhs;
}

template<RMathType T, RMathType U>
inline constexpr T& operator*=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value *= unwrap(rhs);
  }
  return lhs;
}

template<RMathType T, CppMathType U>
inline constexpr T& operator*=(T &lhs, U rhs) {
  if (is_na(lhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value *= rhs;
  }
  return lhs;
}

template<CppMathType T, RMathType U>
inline constexpr T& operator*=(T &lhs, U rhs) {
  if (is_na(rhs)){
    lhs = na_value<T>();
  } else {
    lhs *= unwrap(rhs);
  }
  return lhs;
}

template<>
inline constexpr r_dbl& operator*=(r_dbl &lhs, r_dbl rhs) { 
  lhs.value *= unwrap(rhs);
  return lhs;
}

template<RMathType T, RMathType U>
inline constexpr T& operator/=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value /= unwrap(rhs);
  }
  return lhs;
}

template<RMathType T, CppMathType U>
inline constexpr T& operator/=(T &lhs, U rhs) {
  if (is_na(lhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value /= rhs;
  }
  return lhs;
}

template<CppMathType T, RMathType U>
inline constexpr T& operator/=(T &lhs, U rhs) {
  if (is_na(rhs)){
    lhs = na_value<T>();
  } else {
    lhs /= unwrap(rhs);
  }
  return lhs;
}

template<>
inline constexpr r_dbl& operator/=(r_dbl &lhs, r_dbl rhs) {
  lhs.value /= unwrap(rhs);
  return lhs;
}

template<RMathType T>
inline constexpr T operator-(T x) {
  return is_na(x) ? x : T{-unwrap(x)};
}
template<>
inline constexpr r_dbl operator-(r_dbl x) {
  return r_dbl{-unwrap(x)};
}

template <MathType T, MathType U>
inline constexpr r_lgl between(T x, U lo, U hi){
  return x >= lo && x <= hi;
}

}

#endif
