#ifndef CHEAPR_R_RTYPE_COERCE_H
#define CHEAPR_R_RTYPE_COERCE_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_limits.h>
#include <cheapr/internal/r_nas.h>
#include <cheapr/internal/r_vector_utils.h>
#include <charconv>

namespace cheapr {

namespace internal {

// Assumes no NAs at all
template<typename T>
inline constexpr bool can_be_int(T x){
  using xt = std::remove_cvref_t<T>;
  if constexpr (CppIntegerType<xt> && sizeof(xt) <= sizeof(int)){
    // Check if unsigned type's max exceeds signed int range
    if constexpr (
        std::is_unsigned_v<xt> && std::numeric_limits<xt>::max() > static_cast<xt>(std::numeric_limits<int>::max())){
      return x <= static_cast<xt>(std::numeric_limits<int>::max());
    } else {
      return true;  // Small types can safely cast to int
    }
    // Larger types that can safely cast to T
  } else if constexpr (CppMathType<xt>){
    return between<xt>(x, r_limits::r_int_min, r_limits::r_int_max);
  } else if constexpr (RMathType<xt>){
    return between(x.value, static_cast<decltype(x.value)>(r_limits::r_int_min), static_cast<decltype(x.value)>(r_limits::r_int_max));
  } else {
    return false;
  }
}
template<typename T>
inline constexpr bool can_be_int64(T x){
  using xt = std::remove_cvref_t<T>;
  if constexpr (CppIntegerType<xt> && sizeof(xt) <= sizeof(int64_t)){
    // Check if unsigned type's max exceeds signed int64 range
    if constexpr (
        std::is_unsigned_v<xt> &&
          std::numeric_limits<xt>::max() > static_cast<xt>(std::numeric_limits<int64_t>::max())){
      return x <= static_cast<xt>(std::numeric_limits<int64_t>::max());
    } else {
      return true;  // Small types can safely cast to int64
    }
  } else if constexpr (CppMathType<xt>){
    return between<xt>(x, r_limits::r_int64_min, r_limits::r_int64_max);
  } else if constexpr (RMathType<xt>){
    return between(x.value, static_cast<decltype(x.value)>(r_limits::r_int64_min), static_cast<decltype(x.value)>(r_limits::r_int64_max));
  } else {
    return false;
  }
}

// Coerce functions that account for NA
template<typename T>
inline r_lgl as_bool(T x){
  if constexpr (is<T, int> || is<T, r_lgl>){
    return static_cast<r_lgl>(x);
  } else if constexpr (MathType<T>){
    return is_r_na(x) ? na::logical : static_cast<r_lgl>(static_cast<bool>(x));
  } else {
    return na::logical;
  }
}
template<typename T>
inline r_int as_int(T x){
  if constexpr (is<T, int> || is<T, r_int>){
    return static_cast<r_int>(x);
  } else if constexpr (MathType<T>){
    return is_r_na(x) || !internal::can_be_int(x) ? na::integer : static_cast<r_int>(x);
  } else {
    return na::integer;
  }
}
template<typename T>
inline r_int64 as_int64(T x){
  if constexpr (is<T, r_int64>){
    return x;
  } else if constexpr (MathType<T>){
    return is_r_na(x) || !internal::can_be_int64(x) ? na::integer64 : static_cast<r_int64>(x);
  } else {
    return na::integer64;
  }
}
template<typename T>
inline r_dbl as_double(T x){
  if constexpr (is<T, double> || is<T, r_dbl>){
    return static_cast<r_dbl>(x);
  } else if constexpr (MathType<T>){
    return is_r_na(x) ? na::real : static_cast<r_dbl>(x);
  } else {
    return na::real;
  }
}
template<typename T>
inline r_cplx as_complex(T x){
  if constexpr (is<T, r_cplx>){
    return x;
  } else if constexpr (MathType<T>){
    return r_cplx{as_double(x), r_dbl(0.0)};
  } else {
    return na::complex;
  }
}
template<typename T>
inline r_raw as_raw(T x){
  if constexpr (is<T, r_raw>){
    return x;
  } else if constexpr (IntegerType<T> && sizeof(T) <= sizeof(int8_t)){
    return is_r_na(x) || x < 0 ? na::raw : static_cast<r_raw>(x);
  } else if constexpr (IntegerType<T>){
    return is_r_na(x) || !between(x, static_cast<T>(0), static_cast<T>(255)) ? na::raw : static_cast<r_raw>(x);
  } else {
    return na::raw;
  }
}
// As CHARSXP
template<typename T>
inline r_str as_r_string(T x){
  if constexpr (is<T, r_str>){
    return x;
  } else if constexpr (is<T, const char *>){
    return static_cast<r_str>(r_str(x));
  } else if constexpr (is<T, r_sym>){
    return static_cast<r_str>(PRINTNAME(static_cast<SEXP>(x)));
  } else if constexpr (is<T, r_lgl>){
    if (is_r_na(x)){
      return na::string;
    } else if (x == r_true){
      return as_r_string("TRUE");
    } else {
      return as_r_string("FALSE");
    }
  } else if constexpr (RMathType<T>){
    if (is_r_na(x)){
      return na::string;
    }
    char buffer[48];
    auto result = std::to_chars(buffer, buffer + sizeof(buffer), x.value);
    if (result.ec != std::errc{}) {
      Rf_error("Internal error, increase buffer size for string conversion");
    }
    *result.ptr = '\0';  // Null-terminate
    return as_r_string(static_cast<const char *>(buffer));
  } else if constexpr (CppMathType<T>){
    char buffer[48];
    auto result = std::to_chars(buffer, buffer + sizeof(buffer), x);
    if (result.ec != std::errc{}) {
      Rf_error("Internal error, increase buffer size for string conversion");
    }
    *result.ptr = '\0';  // Null-terminate
    return as_r_string(static_cast<const char *>(buffer));
  } else if constexpr (is<T, r_cplx>){
    if (is_r_na(x)){
      return na::string;
    }
    double re = x.re();
    double im = x.im();

    char buffer[96];
    if (im >= 0){
      snprintf(buffer, sizeof(buffer), "%g+%gi", re, im);
    } else {
      snprintf(buffer, sizeof(buffer), "%g%gi", re, im);
    }
    return as_r_string(static_cast<const char *>(buffer));
  } else if constexpr (is<T, r_raw>){
    char buffer[8];
    snprintf(buffer, sizeof(buffer), "%02x", x.value);
    return as_r_string(static_cast<const char *>(buffer));
  } else if constexpr (is<T, SEXP> || is<T, r_sexp>){
    if (Rf_length(x) != 1){
      Rf_error("`x` is a non-scalar vector and cannot be converted to an `r_str` in %s", __func__);
    }
    r_sexp str = r_sexp(Rf_coerceVector(x, STRSXP));
    r_str out = r_str(STRING_ELT(str, 0));
    return out;
  } else {
    static_assert(always_false<T>, "Unsupported type for `as_r_string`");
  }
}

// As SYMSXP
template<typename T>
inline r_sym as_r_sym(T x){
  if constexpr (is<T, r_sym>){
    return x;
  } else if constexpr (is<T, const char *>){
    SEXP str = r_str(x);
    r_sym out = r_sym(Rf_installChar(str));
    return out;
  } else if constexpr (is<T, r_str>){
    return r_sym(Rf_installChar(x));
  } else {
    r_str str = as_r_string(x);
    r_sym out = as_r_sym(str);
    return out;
  }
}

// CHARSXP is always converted to STRSXP here, see `r_types.h` for info
template<typename T>
inline r_sexp as_sexp(T x){
  if constexpr (is<T, r_sexp>){
    return x;
  } else if constexpr (std::is_convertible_v<T, SEXP>){
    return r_sexp(static_cast<SEXP>(x));
  } else if constexpr (RType<T>){
    return r_sexp(new_scalar_vector(x));
  } else if constexpr (MathType<T>){
    if constexpr (internal::can_be_int<T>){
      return r_sexp(new_scalar_vector(r_int(static_cast<int>(x))));
    } else {
      return r_sexp(new_scalar_vector(r_dbl(static_cast<double>(x))));
    }
  } else {
    static_assert(
      always_false<T>,
      "Unimplemented `as_sexp` specialisation"
    );
    return r_null;
  }
}

template<>
inline r_sexp as_sexp<bool>(bool x){
  return r_sexp(new_scalar_vector(r_lgl(static_cast<int>(x))));
}
template<>
inline r_sexp as_sexp<const char *>(const char *x){
  return r_sexp(Rf_ScalarString(r_str(x)));
}
template<>
inline r_sexp as_sexp<r_sym>(r_sym x){
  return r_sexp(static_cast<SEXP>(x));
}

template<>
inline r_sexp as_sexp<SEXP>(SEXP x){ 
  return r_sexp(x);
}

// R version of static_cast
template<typename T, typename U>
struct as_impl {
  static T cast(U x) {
    static_assert(
      always_false<T>,
      "Can't `as` this type, use `static_cast`"
    );
    return T{};
  }
};

// Specializations for each target type

template<typename U>
struct as_impl<r_lgl, U> {
  static constexpr r_lgl cast(U x) {
    return as_bool(x);
  }
};

template<typename U>
struct as_impl<r_int, U> {
  static constexpr r_int cast(U x) {
    return as_int(x);
  }
};

template<typename U>
struct as_impl<r_int64, U> {
  static constexpr r_int64 cast(U x) {
    return as_int64(x);
  }
};

template<typename U>
struct as_impl<r_dbl, U> {
  static constexpr r_dbl cast(U x) {
    return as_double(x);
  }
};

template<typename U>
struct as_impl<r_cplx, U> {
  static constexpr r_cplx cast(U x) {
    return as_complex(x);
  }
};

template<typename U>
struct as_impl<r_raw, U> {
  static constexpr r_raw cast(U x) {
    return as_raw(x);
  }
};

template<typename U>
struct as_impl<r_str, U> {
  static r_str cast(U x) {
    return as_r_string(x);
  }
};

template<typename U>
struct as_impl<r_sym, U> {
  static r_sym cast(U x) {
    return as_r_sym(x);
  }
};

template<typename U>
struct as_impl<r_sexp, U> {
  static r_sexp cast(U x) {
    return as_sexp(x);
  }
};

template<RType T, typename U>
inline T as_r(U x) {
  if constexpr (is<U, T>){ 
    return x;
  } else {
    using r_t = std::remove_cvref_t<T>;
    return internal::as_impl<r_t, U>::cast(x);
  } 
}

}

}

#endif
