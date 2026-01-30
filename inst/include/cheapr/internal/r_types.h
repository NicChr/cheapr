#ifndef CHEAPR_R_TYPES_H
#define CHEAPR_R_TYPES_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_concepts.h>
#include <limits>
#include <string>

// R types

namespace cheapr {

// General SEXP, reserved for everything except CHARSXP and SYMSXP
// Wrapper around cpp11::sexp to benefit from automatic protection (cpp11-managed linked list)
// All credits go to cpp11 authors/maintainers for `cpp11::sexp`
struct r_sexp {
  cpp11::sexp value;
  using value_type = cpp11::sexp;
  
  // Default constructor
  r_sexp() : value(R_NilValue) {}

  // Constructor from SEXP
  explicit r_sexp(SEXP s) : value(s) {}

  // Implicit conversion to SEXP
  operator SEXP() const { return value.data(); }

  r_size_t length() const noexcept {
    return Rf_xlength(value.data());
  }

  r_size_t size() const noexcept {
    return length();
  }

  bool is_null() const { return value.data() == R_NilValue; }
  
  r_str address() const;
};

// bool type, similar to Rboolean
// Implicit coercion to bool (not int) provided no NA
struct r_lgl {
  int value;
  using value_type = int;
  r_lgl() : value{0} {}
  explicit constexpr r_lgl(int x) : value{x} {}
  explicit constexpr r_lgl(bool x) : value{x} {}  
  explicit constexpr operator int() const { return value; }

  explicit operator bool() const;
  constexpr bool is_true() const;
  constexpr bool is_false() const;
  constexpr bool is_na() const;
};

// The 3 possible values of r_lgl
inline constexpr r_lgl r_true{1};
inline constexpr r_lgl r_false{0};
inline constexpr r_lgl r_na{std::numeric_limits<int>::min()};

  inline constexpr bool r_lgl::is_true() const {
    return value == 1;
  }

  inline constexpr bool r_lgl::is_false() const {
    return value == 0;
  }

  inline constexpr bool r_lgl::is_na() const {
    return value == r_na.value;
  }

  inline r_lgl::operator bool() const {
    if (is_na()){
    abort("Cannot implicitly convert NA to bool, please check");
    }
    return static_cast<bool>(value);
  }

// is_na is defined later as a template

// R integer
struct r_int {
  int value;
  using value_type = int;
  r_int() : value{0} {}
  template <CppMathType T>
  requires (internal::can_definitely_be_int<T>())
  explicit constexpr r_int(T x) : value{static_cast<int>(x)} {}
  constexpr operator int() const { return value; }
};
// R double
struct r_dbl {
  double value;
  using value_type = double;
  r_dbl() : value{0} {}
  template <CppMathType T>
  explicit constexpr r_dbl(T x) : value{static_cast<double>(x)} {}
  constexpr operator double() const { return value; }
};
// R integer64 (closely mimicking how bit64 defines it)
struct r_int64 {
  int64_t value;
  using value_type = int64_t;
  r_int64() : value{0} {}
  template <CppMathType T>
  requires (internal::can_definitely_be_int64<T>())
  explicit constexpr r_int64(T x) : value{static_cast<int64_t>(x)} {}
  constexpr operator int64_t() const { return value; }
};

// Alias type for CHARSXP
// r_str must never be converted to `SEXP`/`r_sexp`
// all templates assume that `SEXP`/`r_sexp` is reserved for objects that can safely fit into an R list vector
// Furthermore CHARSXP is a special case because it is essentially the only SEXP that already fits into a non-list vector: a character vector
struct r_str {
  r_sexp value;
  using value_type = r_sexp;
  r_str() : value{R_BlankString} {}
  // Explicit SEXP/const char* -> r_str
  explicit r_str(SEXP x) : value{x} {}
  explicit r_str(r_sexp x) : value(std::move(x)) {}
  explicit r_str(const char *x) : value(Rf_mkCharCE(x, CE_UTF8)) {}
  // Implicit r_str -> SEXP 
  operator SEXP() const { return value; }

  const char *c_str() const {
    return CHAR(value);
  }

  std::string cpp_str() const {
    return static_cast<std::string>(c_str());
  }
};

// Alias type for SYMSXP
struct r_sym {
  r_sexp value;
  using value_type = r_sexp;
  r_sym() : value{R_MissingArg} {}
  explicit r_sym(SEXP x) : value{x} {}
  explicit r_sym(r_sexp x) : value(std::move(x)) {} 
  explicit r_sym(const char *x) : value(Rf_installChar(r_str(x))) {}
  operator SEXP() const { return value; }
};


// Alias type for Rcomplex
struct r_cplx {
  Rcomplex value;
  using value_type = Rcomplex;

  // Constructors
  constexpr r_cplx() : value{0.0, 0.0} {}
  constexpr r_cplx(r_dbl r, r_dbl i) : value{r, i} {}

  // Conversion handling
  explicit constexpr r_cplx(Rcomplex x) : value{x} {}
  constexpr operator Rcomplex() const { return value; }

  // Get real and imaginary parts
  constexpr r_dbl re() const { return r_dbl{value.r}; }
  constexpr r_dbl im() const { return r_dbl{value.i}; }
};

// Alias type for r_raw
struct r_raw {
  Rbyte value;
  using value_type = Rbyte;

  // Constructors
  constexpr r_raw() : value{static_cast<Rbyte>(0)} {}

  // Conversion handling
  explicit constexpr r_raw(Rbyte x) : value{x} {}
  constexpr operator Rbyte() const { return value; }
};

inline r_str r_sexp::address() const {
  char buf[1000];
  std::snprintf(buf, 1000, "%p", static_cast<void*>(value));
  return r_str(buf);
}

namespace internal {

template <typename T>
struct unwrapped_type {
    using type = T;
};

template <>
struct unwrapped_type<cpp11::sexp> {
    using type = SEXP;
};

template <RVal T>
struct unwrapped_type<T> {
    // Recursively call unwrapped_type on the inner type
    using type = typename unwrapped_type<typename T::value_type>::type;
};

}

template <typename T>
using unwrap_t = typename internal::unwrapped_type<T>::type;


// Important (recursive) helper to extract the underlying NON-RVal value
// Recursively unwrap until we hit a primitive type
template <typename T>
inline constexpr auto unwrap(const T& x){
if constexpr (RVal<T>){
    return unwrap(x.value);
  } else if constexpr (is<T, cpp11::sexp>){
    return x.data();
  } else {
    return x;
  }
}

// Constants

// R C NULL constant
inline const r_sexp r_null = r_sexp();
// Blank string ''
inline const r_str blank_r_string = r_str();
  
// Coerce to an R type based on the C type (useful for RVal templates)
template<typename T>
inline constexpr auto as_r_val(T x) { 
  if constexpr (RVal<T>){
    return x;
  } else if constexpr (CastableToRVal<T>){
    return static_cast<as_r_val_t<T>>(x);
  } else {
    static_assert(
      always_false<T>,
      "Unsupported type for `as_r_val`"
    );
    return r_null;
  } 
}

template<typename T>
inline constexpr auto as_r_scalar(T x) {
  if constexpr (RVector<T>){
    if (x.length() != 1){
      abort("Vector must be length-1 to be coerced to a scalar");
    }
    auto out = x.get(0);
    
    // Only happens if x is a list
    if (!RScalar<decltype(out)>){
      abort("`x` cannot be coereced to a scalar, first list-element is not a scalar");
    }
    return out;
  }
  else {
    return as_r_val(x);
  } 
}

}

#endif
