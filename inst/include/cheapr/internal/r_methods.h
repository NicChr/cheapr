#ifndef CHEAPR_R_METHODS_H
#define CHEAPR_R_METHODS_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_nas.h>

namespace cheapr {

// Methods for custom R types

// ! operator for r_lgl
inline constexpr r_lgl operator!(r_lgl x) {
  if (is_na(x)) {
    return r_na;
  }
  return r_lgl{x.value == 0};
}

// r_cplx methods

// Unary minus
inline constexpr r_cplx operator-(const r_cplx& x) {
  return r_cplx{r_dbl(-x.re().value), r_dbl(-x.im().value)};
}

// Binary arithmetic operators
inline constexpr r_cplx operator+(const r_cplx& lhs, const r_cplx& rhs) {
  return r_cplx{r_dbl(lhs.re() + rhs.re()), r_dbl(lhs.im() + rhs.im())};
}

inline constexpr r_cplx operator-(const r_cplx& lhs, const r_cplx& rhs) {
  return r_cplx{r_dbl(lhs.re() - rhs.re()), r_dbl(lhs.im() - rhs.im())};
}

inline constexpr r_cplx operator*(const r_cplx& lhs, const r_cplx& rhs) {
  // (a+bi) * (c+di) = (ac-bd) + (ad+bc)i
  double a = lhs.re(), b = lhs.im();
  double c = rhs.re(), d = rhs.im();
  return r_cplx{r_dbl(a*c - b*d), r_dbl(a*d + b*c)};
}

inline constexpr r_cplx operator/(const r_cplx& lhs, const r_cplx& rhs) {
  // (a+bi) / (c+di) = [(ac+bd)/(c²+d²) + (bc-ad)/(c²+d²)i]
  double a = lhs.re(), b = lhs.im();
  double c = rhs.re(), d = rhs.im();
  double denom = c*c + d*d; 
  return r_cplx{r_dbl((a*c + b*d) / denom), r_dbl((b*c - a*d) / denom)};
}

// Compound assignment operators
inline constexpr r_cplx& operator+=(r_cplx& lhs, const r_cplx& rhs) {
  lhs.value.r += rhs.value.r;
  lhs.value.i += rhs.value.i;
  return lhs;
}

inline constexpr r_cplx& operator-=(r_cplx& lhs, const r_cplx& rhs) {
  lhs.value.r -= rhs.value.r;
  lhs.value.i -= rhs.value.i;
  return lhs;
}

inline constexpr r_cplx& operator*=(r_cplx& lhs, const r_cplx& rhs) {
  double a = lhs.re(), b = lhs.im();
  double c = rhs.re(), d = rhs.im();
  lhs.re() = r_dbl(a*c - b*d);
  lhs.im() = r_dbl(a*d + b*c);
  return lhs;
}

inline constexpr r_cplx& operator/=(r_cplx& lhs, const r_cplx& rhs) {
  double a = lhs.re(), b = lhs.im();
  double c = rhs.re(), d = rhs.im();
  double denom = c*c + d*d;
  lhs.re() = r_dbl((a*c + b*d) / denom);
  lhs.im() = r_dbl((b*c - a*d) / denom); 
  return lhs;
}

template<RVal T, RVal U>
inline constexpr r_lgl operator==(const T &lhs, const U &rhs) {

  // Check for NA in either operand
  if (is_na(lhs) || is_na(rhs)) {
    return na::logical;
  }

  if constexpr (is<T, r_cplx> && is<U, r_cplx>){
    // Compare complex types by components
    return r_lgl{lhs.re() == rhs.re() && lhs.im() == rhs.im()};
  } else {
    return r_lgl{lhs.value == rhs.value};
  }
}

template<RVal T, CppScalar U>
inline constexpr r_lgl operator==(const T &lhs, const U &rhs) {

  // Check for NA in either operand
  if (is_na(lhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value == rhs};
}

template<CppScalar T, RVal U>
inline constexpr r_lgl operator==(const T &lhs, const U &rhs) {

  // Check for NA in either operand
  if (is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs == rhs.value};
}

// Other comparison operators
template<RVal T, CppScalar U>
inline constexpr r_lgl operator!=(const T &lhs, const U &rhs) {
  r_lgl eq = (lhs == rhs);
  if (eq.is_na()) return r_na;
  return r_lgl(eq.is_false());
}
template<CppScalar T, RVal U>
inline constexpr r_lgl operator!=(const T &lhs, const U &rhs) {
  r_lgl eq = (lhs == rhs);
  if (eq.is_na()) return r_na;
  return r_lgl(eq.is_false());
}
template<RVal T, RVal U>
inline constexpr r_lgl operator!=(const T &lhs, const U &rhs) {
  r_lgl eq = (lhs == rhs);
  if (eq.is_na()) return r_na;
  return r_lgl(eq.is_false());
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator<(T lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value < rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator<(T lhs, U rhs) {
  if (is_na(lhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value < rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator<(T lhs, U rhs) {
  if (is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs < rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator<=(T lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value <= rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator<=(T lhs, U rhs) {
  if (is_na(lhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value <= rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator<=(T lhs, U rhs) {
  if (is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs <= rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator>(T lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value > rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator>(T lhs, U rhs) {
  if (is_na(lhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value > rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator>(T lhs, U rhs) {
  if (is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs > rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr r_lgl operator>=(T lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value >= rhs.value};
}
template<RMathType T, CppMathType U>
inline constexpr r_lgl operator>=(T lhs, U rhs) {
  if (is_na(lhs)) {
    return na::logical;
  }
  return r_lgl{lhs.value >= rhs};
}
template<CppMathType T, RMathType U>
inline constexpr r_lgl operator>=(T lhs, U rhs) {
  if (is_na(rhs)) {
    return na::logical;
  }
  return r_lgl{lhs >= rhs.value};
}

template<RMathType T, RMathType U>
inline constexpr T& operator+=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value += rhs.value;
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
    lhs += rhs.value;
  }
  return lhs;
}

// Fast specialisation for r_dbl
template<>
inline constexpr r_dbl& operator+=(r_dbl &lhs, r_dbl rhs) {
  lhs.value += rhs.value;
  return lhs;
}

template<typename T, typename U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator+(T lhs, U rhs) {
  auto res = lhs;
  res += rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T& operator-=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value -= rhs.value;
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
    lhs -= rhs.value;
  }
  return lhs;
}
template<>
inline constexpr r_dbl& operator-=(r_dbl &lhs, r_dbl rhs) {
  lhs.value -= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator-(T lhs,  U rhs) {
  auto res = lhs;
  res -= rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T& operator*=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value *= rhs.value;
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
    lhs *= rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_dbl& operator*=(r_dbl &lhs, r_dbl rhs) { 
  lhs.value *= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator*(T lhs, U rhs) {
  auto res = lhs;
  res *= rhs;
  return res;
}

template<RMathType T, RMathType U>
inline constexpr T& operator/=(T &lhs, U rhs) {
  if (is_na(lhs) || is_na(rhs)) {
    lhs = na_value<T>();
  } else {
    lhs.value /= rhs.value;
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
    lhs /= rhs.value;
  }
  return lhs;
}

template<>
inline constexpr r_dbl& operator/=(r_dbl &lhs, r_dbl rhs) {
  lhs.value /= rhs.value;
  return lhs;
}

template<RMathType T, RMathType U>
  requires (AtLeastOneRMathType<T, U>)
inline constexpr T operator/(T lhs, U rhs) {
  auto res = lhs;
  res /= rhs;
  return res;
}

template<RMathType T>
inline constexpr T operator-(T x) {
  if (is_na(x)){
    return x;
  }
  return T{-x.value};
}
template<>
inline constexpr r_dbl operator-(r_dbl x) {
  return r_dbl{-x.value};
}

// template <typename T, typename U>
// inline constexpr r_lgl between(T x, U lo, U hi){
//   return x >= lo && x <= hi;
// }

}

#endif
