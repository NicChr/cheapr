#ifndef CHEAPR_R_COERCE_H
#define CHEAPR_R_COERCE_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_vector.h>

namespace cheapr {

namespace internal {

// Assumes no NAs at all
template<typename T>
inline constexpr bool can_be_int(T x){
  using xt = std::remove_cvref_t<T>;
  if constexpr (std::is_integral_v<xt> && sizeof(xt) <= sizeof(int)){
    // Check if unsigned type's max exceeds signed int range
    if constexpr (
        std::is_unsigned_v<xt> && std::numeric_limits<xt>::max() > static_cast<xt>(std::numeric_limits<int>::max())){
      return between<xt>(x, 0, std::numeric_limits<int>::max());
    } else {
      return true;  // Small types can safely cast to int
    }
    // Larger types that can safely cast to T
  } else if constexpr (MathType<xt>){
    return between<xt>(x, r_limits::r_int_min, r_limits::r_int_max);
  } else {
    return false;
  }
}
template<typename T>
inline constexpr bool can_be_int64(T x){
  using xt = std::remove_cvref_t<T>;
  if constexpr (std::is_integral_v<xt> && sizeof(xt) <= sizeof(int64_t)){
    // Check if unsigned type's max exceeds signed int64 range
    if constexpr (
        std::is_unsigned_v<xt> &&
          std::numeric_limits<xt>::max() > static_cast<xt>(std::numeric_limits<int64_t>::max())){
      return x <= std::numeric_limits<int64_t>::max();
    } else {
      return true;  // Small types can safely cast to int64
    }
  } else if constexpr (MathType<xt>){
    return between<xt>(x, r_limits::r_int64_min, r_limits::r_int64_max);
  } else {
    return false;
  }
}

}

namespace internal {
template<typename T>
inline SEXP as_r_obj(const T x){
  if constexpr (is_r_vector_v<T>){
    return static_cast<SEXP>(x);
  } else if constexpr (std::is_convertible_v<T, SEXP>){
    return static_cast<SEXP>(x);
  } else {
    return static_cast<SEXP>(vec::as_vector(x));
  }
}

template<>
inline SEXP as_r_obj<r_string_t>(const r_string_t x){
  return vec::as_vector<r_string_t>(x);
}
template<>
inline SEXP as_r_obj<r_symbol_t>(const r_symbol_t x){
  return static_cast<SEXP>(x);
}


// Coerce functions that account for NA
template<typename T>
inline constexpr r_bool_t as_bool(T x){
  if constexpr (is<T, int> || is<T, r_bool_t>){
    return static_cast<r_bool_t>(x);
  } else if constexpr (MathType<T>){
    return is_r_na(x) ? na::logical : static_cast<r_bool_t>(static_cast<bool>(x));
  } else {
    return na::logical;
  }
}
template<typename T>
inline constexpr r_int_t as_int(T x){
  if constexpr (is<T, int> || is<T, r_int_t>){
    return static_cast<r_int_t>(x);
  } else if constexpr (MathType<T>){
    return is_r_na(x) || !internal::can_be_int(x) ? na::integer : static_cast<r_int_t>(x);
  } else {
    return na::integer;
  }
}
template<typename T>
inline constexpr r_int64_t as_int64(T x){
  if constexpr (is<T, r_int64_t>){
    return x;
  } else if constexpr (MathType<T>){
    return is_r_na(x) || !internal::can_be_int64(x) ? na::integer64 : static_cast<r_int64_t>(x);
  } else {
    return na::integer64;
  }
}
template<typename T>
inline constexpr r_double_t as_double(T x){
  if constexpr (is<T, double> || is<T, r_double_t>){
    return static_cast<r_double_t>(x);
  } else if constexpr (MathType<T>){
    return is_r_na(x) ? na::real : static_cast<r_double_t>(x);
  } else {
    return na::real;
  }
}
template<typename T>
inline constexpr r_complex_t as_complex(T x){
  if constexpr (is<T, r_complex_t>){
    return x;
  } else if constexpr (MathType<T>){
    return r_complex_t{as_double(x), 0.0};
  } else {
    return na::complex;
  }
}
template<typename T>
inline constexpr r_byte_t as_raw(T x){
  if constexpr (is<T, r_byte_t>){
    return x;
  } else if constexpr (is_r_or_cpp_integral_v<T> && sizeof(T) <= sizeof(int8_t)){
    return is_r_na(x) || x < 0 ? na::raw : static_cast<r_byte_t>(x);
  } else if constexpr (std::is_convertible_v<T, r_byte_t>){
    return is_r_na(x) || !between(x, static_cast<T>(0), static_cast<T>(255)) ? na::raw : static_cast<r_byte_t>(x);
  } else {
    return na::raw;
  }
}
// As CHARSXP
template<typename T>
inline r_string_t as_r_string(T x){
  if constexpr (is<T, r_string_t>){
    return x;
  } else if constexpr (is<T, const char *>){
    return static_cast<r_string_t>(internal::make_utf8_charsxp(x));
  } else if constexpr (is<T, r_symbol_t>){
    return static_cast<r_string_t>(PRINTNAME(static_cast<SEXP>(x)));
  } else {
    SEXP scalar = SHIELD(as_r_obj(x));
    if (Rf_length(scalar) != 1){
      YIELD(1);
      Rf_error("`x` is a non-scalar vector and cannot be convered to an `r_string_t` in %s", __func__);
    }
    if (is_r_na(x)){
      YIELD(1);
      return na::string;
    }
    SEXP str = SHIELD(internal::coerce_vec(scalar, STRSXP));
    r_string_t out = internal::get_value<r_string_t>(str, 0);
    YIELD(2);
    return out;
  }
}

// As SYMSXP
template<typename T>
inline r_symbol_t as_r_sym(T x){
  if constexpr (is<T, r_symbol_t>){
    return x;
  } else if constexpr (is<T, const char *>){
    return static_cast<r_symbol_t>(Rf_install(x));
  } else if constexpr (is<T, r_string_t>){
    return as_r_sym(CHAR(static_cast<SEXP>(x)));
  } else {
    SEXP scalar = SHIELD(as_r_obj(x));
    if (Rf_length(scalar) != 1){
      YIELD(1);
      Rf_error("`x` is a non-scalar vector and cannot be convered to an `r_symbol_t` in %s", __func__);
    }
    SEXP str = SHIELD(internal::coerce_vec(scalar, STRSXP));
    r_symbol_t out = as_r_sym(internal::get_value<r_string_t>(str, 0));
    YIELD(2);
    return out;
  }
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
struct as_impl<r_bool_t, U> {
  static constexpr r_bool_t cast(U x) {
    return as_bool(x);
  }
};

template<typename U>
struct as_impl<r_int_t, U> {
  static constexpr r_int_t cast(U x) {
    return as_int(x);
  }
};

template<typename U>
struct as_impl<r_int64_t, U> {
  static constexpr r_int64_t cast(U x) {
    return as_int64(x);
  }
};

template<typename U>
struct as_impl<r_double_t, U> {
  static constexpr r_double_t cast(U x) {
    return as_double(x);
  }
};

template<typename U>
struct as_impl<r_complex_t, U> {
  static constexpr r_complex_t cast(U x) {
    return as_complex(x);
  }
};

template<typename U>
struct as_impl<r_byte_t, U> {
  static constexpr r_byte_t cast(U x) {
    return as_raw(x);
  }
};

template<typename U>
struct as_impl<r_string_t, U> {
  static r_string_t cast(U x) {
    return as_r_string(x);
  }
};

template<typename U>
struct as_impl<r_symbol_t, U> {
  static r_symbol_t cast(U x) {
    return as_r_sym(x);
  }
};

template<typename U>
struct as_impl<SEXP, U> {
  static SEXP cast(U x) {
    return as_r_obj(x);
  }
};
}

template<typename T, typename U>
inline constexpr T as(U x) {
  if constexpr (is<U, T>){ 
    return x;
  } else {
  using r_t = std::remove_cvref_t<T>;
  // if constexpr (is_r_vector_v<r_t>){    
  // using data_t = std::remove_pointer_t<std::remove_cvref_t<decltype(x.data())>>;
  // return as_vector(as<data_t>(x));
  // } else {
  return internal::as_impl<r_t, U>::cast(x);
  } 
}


// Get the R type from the C type (useful for mixed type operators)

namespace internal {
template<typename T>
inline auto as_r_type(T x) {
  if constexpr (is<T, bool>) {
    return r_bool_t{x};
  } else if constexpr (is<T, int>) {
    return r_int_t{x};
  } else if constexpr (is<T, double>) {
    return r_double_t{x};
  } else if constexpr (is<T, r_string_t>) {
    return r_string_t{x};
  } else if constexpr (is<T, Rcomplex>) {
    return r_complex_t{x};
  }  else if constexpr (is<T, Rbyte>) {
    return r_byte_t{x};
  } else {
  static_assert(
    always_false<T>,
    "Unsupported type for `as_r_type`"
  );
  return T{};
  } 
}

}

}

#endif
