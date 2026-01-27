#ifndef CHEAPR_R_TYPES_H
#define CHEAPR_R_TYPES_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_concepts.h>
#include <limits>
#include <string>

// R types

namespace cheapr {

// Internal struct to provide a means to convert SEXP -> r_sexp directly without extra protection
namespace internal {
struct read_only_tag {};
}

// General SEXP, reserved for everything except CHARSXP and SYMSXP
// Wrapper around cpp11::sexp to benefit from automatic protection (cpp11-managed linked list)
// All credits go to cpp11 authors/maintainers for `cpp11::sexp`
struct r_sexp {
  
// value must be declared first (before protector_) as it is the primary data
public:
  SEXP value;

private: 
  // cpp11::sexp will automatically protect underlying SEXP
  //  We never use this directly
  cpp11::sexp protector_;

public:
  // Constructor from SEXP
  explicit r_sexp(SEXP s) : value(s), protector_(s) {}

  // Default constructor
  r_sexp() : value(R_NilValue) {}

  // Optimized constructor
  // convert SEXP -> r_sexp directly without extra protection
  r_sexp(SEXP s, internal::read_only_tag) : value(s) {}

  // Copy Constructor
  r_sexp(const r_sexp& other) : value(other.value), protector_(other.value) {}

  // Move Constructor
  r_sexp(r_sexp&& other) noexcept 
  : value(other.value), protector_(std::move(other.protector_)) {
    other.value = R_NilValue;
  }

  // Copy Assignment
  r_sexp& operator=(const r_sexp& other) {
      if (this != &other) {
          value = other.value;
          // protector_ = other.value; // Previously was this line (investigate)
          protector_ = other.protector_;
      }
      return *this;
  }

  // Move Assignment
  r_sexp& operator=(r_sexp&& other) noexcept {
      if (this != &other) {
          value = other.value;
          protector_ = std::move(other.protector_);
          other.value = R_NilValue;
      }
      return *this;
  }

  // Implicit conversion to SEXP
  constexpr operator SEXP() const { return value; }

  r_size_t length() const noexcept {
    return Rf_xlength(value);
  }

  r_size_t size() const noexcept {
    return length();
  }

  bool is_null() const { return value == R_NilValue; }
  
  r_str address() const;
};

// bool type, similar to Rboolean
// Implicit coercion to bool (not int) provided no NA
struct r_lgl {
  int value;
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
    cpp11::stop("Cannot implicitly convert NA to bool, please check");
    }
    return static_cast<bool>(value);
  }

// is_na is defined later as a template

// R integer
struct r_int {
  int value;
  r_int() : value{0} {}
  explicit constexpr r_int(int x) : value{x} {}
  constexpr operator int() const { return value; }
};
// R double
struct r_dbl {
  double value;
  r_dbl() : value{0} {}
  explicit constexpr r_dbl(double x) : value{x} {}
  constexpr operator double() const { return value; }
};
// R integer64 (closely mimicking how bit64 defines it)
struct r_int64 {
  int64_t value;
  r_int64() : value{0} {}
  explicit constexpr r_int64(int64_t x) : value{x} {}
  constexpr operator int64_t() const { return value; }
};

// Alias type for CHARSXP
// r_str must never be converted to `SEXP`/`r_sexp`
// all templates assume that `SEXP`/`r_sexp` is reserved for objects that can safely fit into an R list vector
// Furthermore CHARSXP is a special case because it is essentially the only SEXP that already fits into a non-list vector: a character vector
struct r_str {
  r_sexp value;
  r_str() : value{R_BlankString} {}
  // Explicit SEXP/const char* -> r_str
  explicit r_str(SEXP x) : value{x} {}
  explicit r_str(r_sexp x) : value(std::move(x)) {}
  explicit r_str(const char *x) : value(Rf_mkCharCE(x, CE_UTF8)) {}
  // Implicit r_str -> SEXP 
  constexpr operator SEXP() const { return value.value; }

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
  r_sym() : value{R_MissingArg} {}
  explicit r_sym(r_sexp x) : value(std::move(x)) {} 
  explicit r_sym(SEXP x) : value{std::move(r_sexp(x, internal::read_only_tag{}))} {} // Assume symbols are already protected
  constexpr operator SEXP() const { return value.value; }
};


// Alias type for Rcomplex
struct r_cplx {
  Rcomplex value;

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


// TO-DO: The same thing as below but returning a typename

// This is the simplest way to write unwrap() but we're relying on compiler to properly optimise it
// template <typename T>
// inline constexpr auto unwrap(const T& x){
//   if constexpr (RVal<T>){
//       return unwrap(x.value);
//   } else {
//     return x;
//   }
// }

// Important (recursive) helper to extract the underlying NON-RVal value
// Recursively unwrap until we hit a primitive type
template <typename T>
inline constexpr auto unwrap(const T& x){
  if constexpr (RVal<T>){
    if constexpr (!RVal<decltype(x.value)>){
      return x.value;
    } else if constexpr (!RVal<decltype(x.value.value)>){
      return x.value.value;
    } else {
      return unwrap(x.value.value);
    }
  } else {
    return x;
  }
}

// Constants

// R C NULL constant
inline const r_sexp r_null = r_sexp();
// Blank string ''
inline const r_str blank_r_string = r_str();

namespace internal {

template<typename T>
inline constexpr bool can_definitely_be_int(){

  constexpr int max_int = std::numeric_limits<int>::max();

  using xt = std::remove_cvref_t<T>;
  if constexpr (CppIntegerType<xt> && sizeof(xt) <= sizeof(int)){
    // Check if unsigned type's max exceeds signed int range
    if constexpr (std::is_unsigned_v<xt> && std::numeric_limits<xt>::max() <= static_cast<xt>(max_int)){
      return true; // Small types can safely cast to int
    } else {
      return false;
    }
  } else {
    return false;
  }
}

template<typename T>
inline constexpr bool can_definitely_be_int64(){

  constexpr int64_t max_int64 = std::numeric_limits<int64_t>::max();

  using xt = std::remove_cvref_t<T>;
  if constexpr (CppIntegerType<xt> && sizeof(xt) <= sizeof(int64_t)){
    // Check if unsigned type's max exceeds signed int64 range
    if constexpr (std::is_unsigned_v<xt> && std::numeric_limits<xt>::max() <= static_cast<xt>(max_int64)){
      return true; // Small types can safely cast to int64
    } else {
      return false;
    }
  } else {
    return false;
  }
}


// Assumes no NAs at all
template<typename T>
inline constexpr bool can_be_int(T x){
  constexpr int max_int = std::numeric_limits<int>::max();
  constexpr int min_int = -max_int; // Doesn't include lowest int (reserved for NA)

  if constexpr (can_definitely_be_int<T>()){
    return true;
 } else if constexpr (CppMathType<T>){
    using data_t = decltype(x);
    return internal::between_impl<data_t>(x, min_int, max_int);
  } else if constexpr (RMathType<T>){
    using data_t = decltype(x.value);
    return internal::between_impl<data_t>(x.value, min_int, max_int);
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
 } else if constexpr (CppMathType<T>){
    using data_t = decltype(x);
    return internal::between_impl<data_t>(x, min_int64, max_int64);
  } else if constexpr (RMathType<T>){
    using data_t = decltype(x.value);
    return internal::between_impl<data_t>(x.value, min_int64, max_int64);
  } else {
    return false;
  }
}

}
  
// Coerce to an R type based on the C type (useful for RVal templates)
template<typename T>
inline constexpr auto as_r_val(T x) { 
  if constexpr (RVal<T>){
    return x;
  } else if constexpr (ConstructibleToRVal<T>){
    return to_r_val_t<T>(x);
  } else if constexpr (MathType<T>){
    if constexpr (internal::can_definitely_be_int<T>()){
      return r_int(static_cast<int>(x));
    } else {
      return r_dbl(static_cast<double>(x));
    }
  } else if constexpr (RVector<T>){
    return x.sexp;
  } else if constexpr (is<T, SEXP>){
    return r_sexp(x);
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
      cpp11::stop("Vector must be length-1 to be coerced to a scalar");
    }
    auto out = x.get(0);
    
    // Only happens if x is a list
    if (!RScalar<decltype(out)>){
      cpp11::stop("`x` cannot be coereced to a scalar, first list-element is not a scalar");
    }
    return out;
  }
  else {
    return as_r_val(x);
  } 
}

}

#endif
