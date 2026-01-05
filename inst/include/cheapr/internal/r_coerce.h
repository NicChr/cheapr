#ifndef CHEAPR_R_COERCE_H
#define CHEAPR_R_COERCE_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_limits.h>
#include <cheapr/internal/r_nas.h>
#include <cheapr/internal/r_vector_utils.h>

namespace cheapr {

namespace internal {

// Coerce to an R type from the C type (useful for RType templates)
// template<typename T>
// inline auto as_r_type(T x) {
//   if constexpr (RType<T>){
//     return x;
//   } else if (CppMathType>){
//     if constexpr (internal::can_be_int<T>){
//       return r_int_t(static_cast<int>(x));
//     } else {
//       r_double_t(static_cast<double>(x));
//     }
//   } else {
//   static_assert(
//     always_false<T>,
//     "Unsupported type for `as_r_type`"
//   );
//   } 
// }

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
  } else if constexpr (IntegerType<T> && sizeof(T) <= sizeof(int8_t)){
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
  } else if constexpr (is<T, r_bool_t>){
    if (x == r_true){
      return as_r_string("TRUE");
    } else if (x == r_false){
      return as_r_string("FALSE");
    } else {
      return na::string;
    }
  } else if constexpr (RMathType<T>){
    if (is_r_na(x)){
      return na::string;
    }
    std::string str = std::to_string(x.value);
    return as_r_string(str.c_str());
  } else if constexpr (MathType<T>){
    std::string str = std::to_string(x);
    return as_r_string(str.c_str());
  } else if constexpr (is<T, r_complex_t>){
    if (is_r_na(x)){
      return na::string;
    }
    double re = x.re();
    double im = x.im();

    std::string str;
    if (im >= 0){
      str = std::to_string(re) + "+" + std::to_string(im) + "i";
    } else {
      str = std::to_string(re) + std::to_string(im) + "i";
    }
    return as_r_string(str.c_str());
  } else if constexpr (is<T, SEXP>){
    if (Rf_length(x) != 1){
      Rf_error("`x` is a non-scalar vector and cannot be converted to an `r_string_t` in %s", __func__);
    }
    SEXP str = SHIELD(Rf_coerceVector(x, STRSXP));
    r_string_t out = r_string_t(STRING_ELT(str, 0));
    YIELD(1);
    return out;
  } else {
    static_assert(always_false<T>, "Unsupported type for `as_r_string`");
  }
}

// As SYMSXP
template<typename T>
inline r_symbol_t as_r_sym(T x){
  if constexpr (is<T, r_symbol_t>){
    return x;
  } else if constexpr (is<T, const char *>){
    SEXP str = SHIELD(internal::make_utf8_charsxp(x));
    r_symbol_t out = r_symbol_t(Rf_installChar(str));
    YIELD(1);
    return out;
  } else if constexpr (is<T, r_string_t>){
    return r_symbol_t(Rf_installChar(x));
  } else {
    r_string_t str = SHIELD(as_r_string(x));
    r_symbol_t out = as_r_sym(str);
    YIELD(1);
    return out;
  }
}

// CHARSXP is always converted to STRSXP here, see `r_types.h` for info
template<typename T>
inline sexp_t as_sexp(T x){
  if constexpr (is<T, sexp_t>){
    return x;
  } else if constexpr (RType<T>){
    return sexp_t(new_scalar_vector(x));
  } else if constexpr (IntegerType<T>){
    if constexpr (internal::can_be_int<T>){
      return sexp_t(new_scalar_vector(r_int_t(static_cast<int>(x))));
    } else {
      return sexp_t(new_scalar_vector(r_double_t(static_cast<double>(x))));
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
inline sexp_t as_sexp<const char *>(const char *x){
  return sexp_t(internal::make_utf8_strsxp(x));
}
template<>
inline sexp_t as_sexp<r_symbol_t>(r_symbol_t x){
  return sexp_t(static_cast<SEXP>(x));
}

template<>
inline sexp_t as_sexp<SEXP>(SEXP x){ 
  return sexp_t(x);
}

// template<typename T>
// inline sexp_t as_sexp(const T x){
//   if constexpr (is<T, SEXP> || std::is_convertible_v<T, SEXP>){
//     return sexp_t(static_cast<SEXP>(x));
//   } else {
//     return sexp_t(static_cast<SEXP>(vec::as_vector(x)));
//   }
// }

// template<>
// inline sexp_t as_sexp<r_string_t>(const r_string_t x){
//   return sexp_t(Rf_ScalarString(x));
// }
// template<>
// inline sexp_t as_sexp<r_symbol_t>(const r_symbol_t x){
//   return sexp_t(static_cast<SEXP>(x));
// }

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
struct as_impl<sexp_t, U> {
  static sexp_t cast(U x) {
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
