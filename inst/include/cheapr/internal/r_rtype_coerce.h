#ifndef CHEAPR_R_RTYPE_COERCE_H
#define CHEAPR_R_RTYPE_COERCE_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_limits.h>
#include <cheapr/internal/r_nas.h>
#include <cheapr/internal/r_vec_utils.h>
#include <charconv>

namespace cheapr {

namespace internal {

// Assumes no NAs at all
template<typename T>
inline constexpr bool can_be_int(T x){
  constexpr int max_int = std::numeric_limits<int>::max();
  constexpr int min_int = -max_int; // Doesn't include lowest int (reserved for NA)

  if constexpr (can_definitely_be_int<T>()){
    return true;
  } else if constexpr (MathType<T>){
    // This should be a 'practical' way to get the wider type of the 2
    using common_t = std::common_type_t<unwrap_t<T>, int>;
    return internal::between_impl<common_t>(unwrap(x), min_int, max_int);
  } else {
    return false;
  }
}
template<typename T>
inline constexpr bool can_be_int64(T x){
  constexpr int64_t max_int64 = std::numeric_limits<int64_t>::max();
  constexpr int64_t min_int64 = -max_int64; // Doesn't include lowest int (reserved for NA)

  if constexpr (can_definitely_be_int64<T>()){
    return true;
  } else if constexpr (MathType<T>){
    using common_t = std::common_type_t<unwrap_t<T>, int64_t>;
    return internal::between_impl<common_t>(unwrap(x), min_int64, max_int64);
  } else {
    return false;
  }
}

// Coerce functions that account for NA
template<typename T>
inline r_lgl as_bool(T x){
  if constexpr (is<T, int> || is<T, r_lgl>){
    return r_lgl(unwrap(x));
  } else if constexpr (MathType<T>){
    return is_na(x) ? na_value<r_lgl>() : r_lgl(static_cast<bool>(unwrap(x)));
  } else {
    return na_value<r_lgl>();
  }
}
template<typename T>
inline r_int as_int(T x){
  if constexpr (is<T, int> || is<T, r_int>){
    return r_int(unwrap(x));
  } else if constexpr (MathType<T>){
    return is_na(x) || !internal::can_be_int(x) ? na_value<r_int>() : r_int(static_cast<int>(unwrap(x)));
  } else {
    return na_value<r_int>();
  }
}
template<typename T>
inline r_int64 as_int64(T x){
  if constexpr (is<T, int64_t> || is<T, r_int64>){
    return r_int64(unwrap(x));
  } else if constexpr (MathType<T>){
    return is_na(x) || !internal::can_be_int64(x) ? na_value<r_int64>() : r_int64(static_cast<int64_t>(unwrap(x)));
  } else {
    return na_value<r_int64>();
  }
}
template<typename T>
inline r_dbl as_double(T x){
  if constexpr (is<T, double> || is<T, r_dbl>){
    return r_dbl(unwrap(x));
  } else if constexpr (MathType<T>){
    return is_na(x) ? na_value<r_dbl>() : r_dbl(static_cast<double>(unwrap(x)));
  } else {
    return na_value<r_dbl>();
  }
}
template<typename T>
inline r_cplx as_complex(T x){
  if constexpr (is<T, std::complex<double>> || is<T, r_cplx>){
    return r_cplx(unwrap(x));
  } else if constexpr (MathType<T>){
    return r_cplx{as_double(x), r_dbl(0.0)};
  } else {
    return na_value<r_cplx>();
  }
}
template<typename T>
inline r_raw as_raw(T x){
  if constexpr (is<T, Rbyte> || is<T, r_raw>){
    return r_raw(unwrap(x));
  } else if constexpr (IntegerType<T> && sizeof(T) <= sizeof(int8_t)){
    return is_na(x) || x < 0 ? na_value<r_raw>() : r_raw(static_cast<Rbyte>(unwrap(x)));
  } else if constexpr (MathType<T>){
    using r_t = unwrap_t<T>;
    return is_na(x) || !internal::between_impl(unwrap(x), r_t(0), r_t(255)) ? na_value<r_raw>() : r_raw(static_cast<Rbyte>(unwrap(x)));
  } else {
    return na_value<r_raw>();
  }
}
// As CHARSXP
template<typename T>
inline r_str as_r_string(T x){
  if constexpr (is<T, r_str>){
    return x;
  } else if constexpr (is<T, const char *>){
    return r_str(x);
  } else if constexpr (is<T, std::string>){
    return r_str(x.c_str());
  } else if constexpr (is<T, r_sym>){
    return r_str(PRINTNAME(static_cast<SEXP>(x)));
  } else if constexpr (is<T, r_lgl>){
    if (is_na(x)){
      return na_value<r_str>();
    } else if (x == r_true){
      return as_r_string("TRUE");
    } else {
      return as_r_string("FALSE");
    }
  } else if constexpr (RMathType<T>){
    if (is_na(x)){
      return na_value<r_str>();
    }
    char buffer[48];
    auto result = std::to_chars(buffer, buffer + sizeof(buffer), x.value);
    if (result.ec != std::errc{}) {
      abort("Internal error, increase buffer size for string conversion");
    }
    *result.ptr = '\0';  // Null-terminate
    return as_r_string(static_cast<const char *>(buffer));
  } else if constexpr (CppMathType<T>){
    char buffer[48];
    auto result = std::to_chars(buffer, buffer + sizeof(buffer), x);
    if (result.ec != std::errc{}) {
      abort("Internal error, increase buffer size for string conversion");
    }
    *result.ptr = '\0';  // Null-terminate
    return as_r_string(static_cast<const char *>(buffer));
  } else if constexpr (is<T, r_cplx>){
    if (is_na(x)){
      return na_value<r_str>();
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
      abort("`x` is a non-scalar vector and cannot be converted to an `r_str` in %s", __func__);
    }
    r_sexp str = r_sexp(Rf_coerceVector(x, STRSXP));
    r_str out = r_str(STRING_ELT(str, 0));
    return out;
  } else {
    static_assert(always_false<T>, "Unsupported type for `as_r_string`");
  }
}

// Unprotected variant of above (only used to assign into char vec)
// template<typename T>
// inline r_str as_maybe_unprotected_r_str(T x){
//   if constexpr (is<T, r_str>){
//     return x;
//   } else if constexpr (is<T, const char *>){
//     return r_str(r_sexp(Rf_mkCharCE(x, CE_UTF8), internal::read_only_tag{}));
//   } else if constexpr (is<T, r_sym>){
//     return static_cast<r_str>(PRINTNAME(static_cast<SEXP>(x)));
//   } else if constexpr (is<T, r_lgl>){
//     if (is_na(x)){
//       return na_value<r_str>();
//     } else if (x == r_true){
//       return as_r_string("TRUE");
//     } else {
//       return as_r_string("FALSE");
//     }
//   } else if constexpr (RMathType<T>){
//     if (is_na(x)){
//       return na_value<r_str>();
//     }
//     char buffer[48];
//     auto result = std::to_chars(buffer, buffer + sizeof(buffer), x.value);
//     if (result.ec != std::errc{}) {
//       abort("Internal error, increase buffer size for string conversion");
//     }
//     *result.ptr = '\0';  // Null-terminate
//     return as_r_string(static_cast<const char *>(buffer));
//   } else if constexpr (CppMathType<T>){
//     char buffer[48];
//     auto result = std::to_chars(buffer, buffer + sizeof(buffer), x);
//     if (result.ec != std::errc{}) {
//       abort("Internal error, increase buffer size for string conversion");
//     }
//     *result.ptr = '\0';  // Null-terminate
//     return as_r_string(static_cast<const char *>(buffer));
//   } else if constexpr (is<T, r_cplx>){
//     if (is_na(x)){
//       return na_value<r_str>();
//     }
//     double re = x.re();
//     double im = x.im();

//     char buffer[96];
//     if (im >= 0){
//       snprintf(buffer, sizeof(buffer), "%g+%gi", re, im);
//     } else {
//       snprintf(buffer, sizeof(buffer), "%g%gi", re, im);
//     }
//     return as_r_string(static_cast<const char *>(buffer));
//   } else if constexpr (is<T, r_raw>){
//     char buffer[8];
//     snprintf(buffer, sizeof(buffer), "%02x", x.value);
//     return as_r_string(static_cast<const char *>(buffer));
//   } else if constexpr (is<T, SEXP> || is<T, r_sexp>){
//     if (Rf_length(x) != 1){
//       abort("`x` is a non-scalar vector and cannot be converted to an `r_str` in %s", __func__);
//     }
//     r_sexp str = r_sexp(Rf_coerceVector(x, STRSXP));
//     r_str out = as_r_string(CHAR(STRING_ELT(str, 0)));
//     return out;
//   } else {
//     static_assert(always_false<T>, "Unsupported type for `as_r_string`");
//   }
// }

// As SYMSXP
template<typename T>
inline r_sym as_r_sym(T x){
  if constexpr (is<T, r_sym>){
    return x;
  } else if constexpr (is<T, const char *>){
    return r_sym(Rf_install(x));
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
  } else if constexpr (RVector<T>){
    return x.sexp;
  } else if constexpr (std::is_convertible_v<T, SEXP>){
    return r_sexp(static_cast<SEXP>(x));
  } else if constexpr (RVal<T>){
    return r_sexp(new_scalar_vec(x));
  } else {
    return new_scalar_vec(as_r_val(x)); 
  }
}

template<>
inline r_sexp as_sexp<bool>(bool x){
  return r_sexp(new_scalar_vec(r_lgl(static_cast<int>(x))));
}
template<>
inline r_sexp as_sexp<const char *>(const char *x){
  return new_scalar_vec(r_str(x));
}
template<>
inline r_sexp as_sexp<r_sym>(r_sym x){
  return x.value;
}
template<>
inline r_sexp as_sexp<r_str>(r_str x){
  return x.value;
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

template<RVal T, typename U>
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
