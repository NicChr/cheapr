#ifndef CHEAPR_R_COERCE_H
#define CHEAPR_R_COERCE_H

#include <cheapr/internal/r_vector.h>
#include <cheapr/internal/r_factor.h>

namespace cheapr {

// Powerful and flexible coercion function that can handle many types and convert to R-spcific C++ types and R vectors
template<typename T, typename U>
  requires (RType<T> || RVectorType<T> || any<T, SEXP, r_factors>)
inline T as(U x) {
  if constexpr (is<U, T>){
    return x;
  } else if constexpr (RType<U> && is<T, SEXP>){
    return static_cast<SEXP>(internal::as_r<r_sexp>(x));
  } else if constexpr (RVectorType<U> && is<T, SEXP>){
    return x.value;
  } else if constexpr (RVectorType<U> && is<T, r_sexp>){
    return r_sexp(x.value);
  } else if constexpr (RVectorType<T> && any<U, SEXP, r_sexp>){
    return internal::visit_vector(x, [&](auto xvec) -> T {
      // This will trigger the branch that checks that both are RVectorType
      return as<T>(xvec);
    });
  } else if constexpr (RType<T> && any<U, SEXP, r_sexp>){
    return internal::visit_vector(x, [&](auto xvec) -> T {
      // Use branch below current branch
      return as<T>(xvec);
    });
  } else if constexpr (RVectorType<U> && RType<T>){
    if (x.length() != 1){
      Rf_error("Vector must be length-1 to be coerced to requested type");
    }
    return internal::as_r<T>(x.get(0));
  } else if constexpr (RVectorType<T> && RType<U>){  
    using data_t = typename T::data_type;
    return r_vec<data_t>(1, internal::as_r<data_t>(x));
  } else if constexpr (RVectorType<U> && RVectorType<T>){
    r_size_t n = x.length();
    auto out = SHIELD(T(n));
    using data_t = typename T::data_type;
    OMP_SIMD
    for (r_size_t i = 0; i < n; ++i){
    out.set(i, internal::as_r<data_t>(x.get(i)));
    }
    YIELD(1);
    return out;
  } else if constexpr (is<T, r_factors>){
    auto str_vec = SHIELD(as<r_vec<r_str>>(x));
    auto out = SHIELD(r_factors(str_vec));
    YIELD(2);
    return out; 
  } else if constexpr (RType<T> && !RVectorType<U>) {
    return internal::as_r<T>(x);
  } else {
    static_assert(always_false<T>, "Unsupported type for `as`");
  }
}

// Convert any C obj to an r_vec<>
template<typename T>
inline auto as_vector(T x){
  if constexpr (RVectorType<T>){
    return x;
  } else if constexpr (any<T, SEXP, r_sexp>){
    static_assert(always_false<T>, "Can't convert `SEXP/r_sexp` to `r_vec<>`, please use `as<>` to convert");
  } else {
    auto rt_val = internal::as_r_type(x);
    return r_vec<decltype(rt_val)>(1, rt_val);
  }
}

}

#endif
