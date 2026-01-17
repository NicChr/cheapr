#ifndef CHEAPR_R_TYPES_H
#define CHEAPR_R_TYPES_H

#include <cheapr/internal/r_setup.h>

// R types

namespace cheapr {

struct r_str; // Forward declaration

// General SEXP, reserved for everything except CHARSXP and SYMSXP
// Wrapper around cpp11::sexp to benefit from automatic protection (cpp11-managed linked list)
// All credits go to cpp11 authors/maintainers for `cpp11::sexp`
#include <cpp11.hpp>

struct r_sexp {

  SEXP value;

private:
    // 2. The manager that keeps 'value' protected
    // We don't use this directly
    cpp11::sexp protector_;

public:
    // Constructor from SEXP
    explicit r_sexp(SEXP s) : value(s), protector_(s) {}

    // Default constructor
    r_sexp() : value(R_NilValue), protector_(R_NilValue) {}

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
            protector_ = other.value;
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

    // Implicit conversion to SEXP (Zero overhead, returns raw pointer)
    constexpr operator SEXP() const { return value; }

    bool is_null() const { return value == R_NilValue; }
    
    r_str address() const;
};

// struct r_sexp {
//   cpp11::sexp value;

//   // Constructor from SEXP
//   explicit r_sexp(SEXP s) : value(s) {}

//   // Default constructor
//   r_sexp() : value(R_NilValue) {}

//   // Implicit conversion to SEXP
//   operator SEXP() const { return value.data(); }

//   bool is_null() const {
//     return value.data() == R_NilValue;
//   }

//   r_str address() const;

// };

// R C NULL constant
inline const r_sexp r_null = r_sexp();

// bool type, similar to Rboolean
// Implicit coercion to bool (not int) provided no NA
struct r_lgl {
  int value;
  r_lgl() : value{0} {}
  explicit constexpr r_lgl(int x) : value{x} {}
  explicit constexpr operator int() const { return value; }

  operator bool() const;
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
    Rf_error("Cannot convert NA to bool, please check");
  }
  return static_cast<bool>(value);
}

// is_r_na is defined later as a template

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
  explicit r_str(SEXP x) : value{r_sexp(x)} {}
  explicit r_str(const r_sexp &x) : value{x} {}
  explicit r_str(const char *x) : value(r_sexp(Rf_mkCharCE(x, CE_UTF8))) {}
  // Implicit r_str -> SEXP
  operator SEXP() const { return value; }

  const char *c_str() const {
    return CHAR(value);
  }

};

inline const r_str blank_r_string = r_str();

// Alias type for SYMSXP
struct r_sym {
  SEXP value;
  r_sym() : value{R_MissingArg} {}
  explicit constexpr r_sym(SEXP x) : value{x} {}
  constexpr operator SEXP() const { return value; }
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

r_str r_sexp::address() const {
  char buf[1000];
  std::snprintf(buf, 1000, "%p", static_cast<void*>(value));
  return r_str(buf);
}

}

#endif
