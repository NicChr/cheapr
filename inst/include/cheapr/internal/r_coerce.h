#ifndef CHEAPR_R_COERCE_H
#define CHEAPR_R_COERCE_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_dates.h>
#include <cheapr/internal/r_posixcts.h>
#include <cheapr/internal/r_factor.h>

namespace cheapr {

template<typename T>
concept RVectorType = internal::is_r_vector_v<T> || is<T, r_dates> || is<T, r_posixcts>;

// Coerce to an R type based on the C type (useful for RType templates)
namespace internal {

template<typename T>
inline auto as_r_type(T x) {
  if constexpr (RType<T>){
    return x;
  } else if constexpr (is<T, bool>){
    return r_lgl(x);
  } else if constexpr (MathType<T>){
    if constexpr (internal::can_be_int<T>){
      return r_int(static_cast<int>(x));
    } else {
      return r_dbl(static_cast<double>(x));
    }
  } else if constexpr (is<T, const char*>){
    return internal::as_r<r_str>(x);
  } else if constexpr (is<T, SEXP>){
    return r_sexp(static_cast<SEXP>(x));
  } else if constexpr (RVectorType<T>){
    return x.sexp;
  } else {    
    static_assert(
      always_false<T>,
      "Unsupported type for `as_r_type`"
    );
  } 
}

}

// Powerful and flexible coercion function that can handle many types and convert to R-spcific C++ types and R vectors
template<typename T, typename U>
  requires (RType<T> || RVectorType<T> || any<T, SEXP, r_factors>)
inline T as(U x) {
  if constexpr (is<U, T>){
    return x;
  } else if constexpr (RVectorType<U> && is<T, r_sexp>){
    return x.sexp;
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
      cpp11::stop("Vector must be length-1 to be coerced to requested type");
    }
    return internal::as_r<T>(x.get(0));
  } else if constexpr (RVectorType<T> && RType<U>){
    using data_t = typename T::data_type;
    return r_vec<data_t>(1, internal::as_r<data_t>(x));
  } else if constexpr (RVectorType<U> && RVectorType<T>){
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
  } else if constexpr (RType<T> && !RVectorType<U>) {
    return internal::as_r<T>(x);
    // If input is not an R type or an R vector type
  } else if constexpr (!RType<U> && !RVectorType<U>){
    return as<T>(internal::as_r_type(x));
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
    return T();
  } else {
    auto rt_val = internal::as_r_type(x);
    return r_vec<decltype(rt_val)>(1, rt_val);
  }
}

}

#endif
