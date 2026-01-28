#ifndef CHEAPR_R_COERCE_H
#define CHEAPR_R_COERCE_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_dates.h>
#include <cheapr/internal/r_posixcts.h>
#include <cheapr/internal/r_factor.h>

namespace cheapr {

// Powerful and flexible coercion function that can handle many types and convert to R-spcific C++ types and R vectors
template<typename T, typename U>
  requires (RVal<T> || RVector<T> || any<T, r_sexp, SEXP, r_factors>)
inline T as(U x) {
  if constexpr (is<U, T>){
    return x;
  } else if constexpr (RVector<U> && is<T, r_sexp>){
    return x.sexp;
  } else if constexpr (RVector<T> && any<U, SEXP, r_sexp>){
    return internal::visit_vector(x, [&](auto xvec) -> T {
      // This will trigger the branch that checks that both are RVector
      return as<T>(xvec);
    });
  } else if constexpr (RVal<T> && any<U, SEXP, r_sexp>){
    return internal::visit_vector(x, [&](auto xvec) -> T {
      // Use branch below current branch
      return as<T>(xvec);
    });
  } else if constexpr (RVector<U> && RVal<T>){
    if (x.length() != 1){
      abort("Vector must be length-1 to be coerced to requested scalar type");
    }
    return internal::as_r<T>(x.get(0));
  } else if constexpr (RVector<T> && RVal<U>){
    using data_t = typename T::data_type;
    return r_vec<data_t>(1, internal::as_r<data_t>(x));
  } else if constexpr (RVector<U> && RVector<T>){
    r_size_t n = x.length();
    auto out = T(n);
    using data_t = typename T::data_type;
    if constexpr (RPtrWritableType<data_t>){
      OMP_SIMD
      for (r_size_t i = 0; i < n; ++i){
      out.set(i, internal::as_r<data_t>(x.get(i)));
      }
    } else {
      for (r_size_t i = 0; i < n; ++i){
      out.set(i, internal::as_r<data_t>(x.get(i)));
      }
    }
    return out;
  } else if constexpr (is<T, r_factors>){
    auto str_vec = as<r_vec<r_str>>(x);
    auto out = r_factors(str_vec);
    return out; 
  } else if constexpr (RVal<T> && !RVector<U>) {
    return internal::as_r<T>(x);
    // If input is not an R type or an R vector type
  } else if constexpr (!RVal<U> && !RVector<U>){
    return as<T>(as_r_val(x));
  } else {
    static_assert(always_false<T>, "Unsupported type for `as`");
  }
}

// Convert any C obj to an r_vec<>
template<typename T>
inline auto as_vector(T x){
  if constexpr (RVector<T>){
    return x;
  } else if constexpr (is_sexp<T>){
    static_assert(always_false<T>, "Can't convert `SEXP/r_sexp` to `r_vec<>`, please use `as<>` to convert");
    return T();
  } else {
    auto rt_val = as_r_val(x);
    return r_vec<decltype(rt_val)>(1, rt_val);
  }
}

}

#endif
