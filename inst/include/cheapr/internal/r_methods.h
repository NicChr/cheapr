#ifndef CHEAPR_R_METHODS
#define CHEAPR_R_METHODS

#include <r_types.h>
#include <r_concepts.h>
#include <r_nas.h>

namespace cheapr {

// Methods for custom R types

  // ! operator for r_bool_t
inline constexpr r_bool_t operator!(r_bool_t x) {
  if (is_r_na(x)) {
    return r_na;
  }
  return r_bool_t{!x.value};
}

// r_complex_t methods

// Unary minus
inline constexpr r_complex_t operator-(const r_complex_t& x) {
  return r_complex_t{-x.re(), -x.im()};
}

// Binary arithmetic operators
inline constexpr r_complex_t operator+(const r_complex_t& lhs, const r_complex_t& rhs) {
  return r_complex_t{lhs.re() + rhs.re(), lhs.im() + rhs.im()};
}

inline constexpr r_complex_t operator-(const r_complex_t& lhs, const r_complex_t& rhs) {
  return r_complex_t{lhs.re() - rhs.re(), lhs.im() - rhs.im()};
}

inline constexpr r_complex_t operator*(const r_complex_t& lhs, const r_complex_t& rhs) {
  // (a+bi) * (c+di) = (ac-bd) + (ad+bc)i
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  return r_complex_t{a*c - b*d, a*d + b*c};
}

inline constexpr r_complex_t operator/(const r_complex_t& lhs, const r_complex_t& rhs) {
  // (a+bi) / (c+di) = [(ac+bd)/(c²+d²) + (bc-ad)/(c²+d²)i]
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  r_double_t denom = c*c + d*d;
  return r_complex_t{(a*c + b*d) / denom, (b*c - a*d) / denom};
}

// Compound assignment operators
inline constexpr r_complex_t& operator+=(r_complex_t& lhs, const r_complex_t& rhs) {
  lhs.value.r += rhs.value.r;
  lhs.value.i += rhs.value.i;
  return lhs;
}

inline constexpr r_complex_t& operator-=(r_complex_t& lhs, const r_complex_t& rhs) {
  lhs.value.r -= rhs.value.r;
  lhs.value.i -= rhs.value.i;
  return lhs;
}

inline constexpr r_complex_t& operator*=(r_complex_t& lhs, const r_complex_t& rhs) {
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  lhs.re() = a*c - b*d;
  lhs.im() = a*d + b*c;
  return lhs;
}

inline constexpr r_complex_t& operator/=(r_complex_t& lhs, const r_complex_t& rhs) {
  r_double_t a = lhs.re(), b = lhs.im();
  r_double_t c = rhs.re(), d = rhs.im();
  r_double_t denom = c*c + d*d;
  lhs.re() = (a*c + b*d) / denom;
  lhs.im() = (b*c - a*d) / denom;
  return lhs;
}

template<RType T, RType U>
inline constexpr r_bool_t operator==(const T lhs, const U rhs) {

  // Check for NA in either operand
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return na::logical;
  }

  if constexpr (is<T, r_complex_t> && is<U, r_complex_t>){
    // Compare complex types by components
    return r_bool_t{lhs.re() == rhs.re() && lhs.im() == rhs.im()};
  } else {
    return r_bool_t{lhs.value == rhs.value};
  }
}

template<RType T, CppType U>
inline constexpr r_bool_t operator==(const T lhs, const U rhs) {

  // Check for NA in either operand
  if (is_r_na(lhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value == rhs};
}

template<CppType T, RType U>
inline constexpr r_bool_t operator==(const T lhs, const U rhs) {

  // Check for NA in either operand
  if (is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs == rhs.value};
}

// Other comparison operators
template<typename T, typename U>
inline constexpr r_bool_t operator!=(const T lhs, const U rhs) {
  return r_bool_t{!(lhs == rhs)};
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator<(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value < rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_bool_t operator<(const T lhs, const U rhs) {
  if (is_r_na(lhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value < rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_bool_t operator<(const T lhs, const U rhs) {
  if (is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs < rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator<=(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value <= rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_bool_t operator<=(const T lhs, const U rhs) {
  if (is_r_na(lhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value <= rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_bool_t operator<=(const T lhs, const U rhs) {
  if (is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs <= rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator>(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value > rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_bool_t operator>(const T lhs, const U rhs) {
  if (is_r_na(lhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value > rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_bool_t operator>(const T lhs, const U rhs) {
  if (is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs > rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_bool_t operator>=(const T lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value >= rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_bool_t operator>=(const T lhs, const U rhs) {
  if (is_r_na(lhs)) {
    return na::logical;
  }
  return r_bool_t{lhs.value >= rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_bool_t operator>=(const T lhs, const U rhs) {
  if (is_r_na(rhs)) {
    return na::logical;
  }
  return r_bool_t{lhs >= rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr T operator+=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value += rhs.value;
  }
  return lhs;
}
template<RMathType T, CppMathType U>
inline constexpr T operator+=(T &lhs, const U rhs) {
  if (is_r_na(lhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value += rhs;
  }
  return lhs;
}
template<CppMathType T, RMathType U>
inline constexpr T operator+=(T &lhs, const U rhs) {
  if (is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs += rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_double_t operator+=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value += rhs.value;
  return lhs;
}

template<typename T, typename U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator+(const T lhs, const U rhs) {
  auto res = lhs;
  res += rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T operator-=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value -= rhs.value;
  }
  return lhs;
}

template<RMathType T, CppMathType U>
inline constexpr T operator-=(T &lhs, const U rhs) {
  if (is_r_na(lhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value -= rhs;
  }
  return lhs;
}

template<CppMathType T, RMathType U>
inline constexpr T operator-=(T &lhs, const U rhs) {
  if (is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs -= rhs.value;
  }
  return lhs;
}
template<>
inline constexpr r_double_t operator-=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value -= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator-(const T lhs, const U rhs) {
  auto res = lhs;
  res -= rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T operator*=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value *= rhs.value;
  }
  return lhs;
}

template<RMathType T, CppMathType U>
inline constexpr T operator*=(T &lhs, const U rhs) {
  if (is_r_na(lhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value *= rhs;
  }
  return lhs;
}

template<CppMathType T, RMathType U>
inline constexpr T operator*=(T &lhs, const U rhs) {
  if (is_r_na(rhs)){
    lhs = na_value<T>();
  } else {
    lhs *= rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_double_t operator*=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value *= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator*(const T lhs, const U rhs) {
  auto res = lhs;
  res *= rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T operator/=(T &lhs, const U rhs) {
  if (is_r_na(lhs) || is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value /= rhs.value;
  }
  return lhs;
}

template<RMathType T, CppMathType U>
inline constexpr T operator/=(T &lhs, const U rhs) {
  if (is_r_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value /= rhs;
  }
  return lhs;
}

template<CppMathType T, RMathType U>
inline constexpr T operator/=(T &lhs, const U rhs) {
  if (is_r_na(rhs)){
    lhs = na_value<T>();
  } else {
    lhs /= rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_double_t operator/=(r_double_t &lhs, const r_double_t rhs) {
  lhs.value /= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator/(const T lhs, const U rhs) {
  auto res = lhs;
  res /= rhs;
  return res;
}

template<RMathType T>
inline constexpr T operator-(const T x) {
  if (is_r_na(x)){
    return x;
  }
  return T{-x.value};
}
template<>
inline constexpr r_double_t operator-(const r_double_t x) {
  return r_double_t{-x.value};
}

}

#endif
