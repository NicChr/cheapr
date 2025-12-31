#include <cpp11.hpp>

// R types


// General SEXP, reserved for everything except R vectors, CHARSXP, and SYMSXP
struct sexp_t {
  SEXP value;

  sexp_t() : value(R_NilValue) {}
  // Implicit coercion both ways
  explicit constexpr sexp_t(SEXP s) : value(s) {}
  constexpr operator SEXP() const { return value; }
};

// R C NULL constant
inline const sexp_t r_null = sexp_t();

// bool type, similar to Rboolean
// Contains bember functions to extract bool
struct r_bool_t {
  int value;
  r_bool_t() : value{0} {}
  explicit constexpr r_bool_t(int x) : value{x} {}
  constexpr operator int() const { return value; }


  constexpr bool as_bool() const {
    return static_cast<bool>(value);
  }

  constexpr bool is_true() const {
    return value == 1;
  }

  constexpr bool is_false() const {
    return value == 0;
  }
};

// The 3 possible values of r_bool_t
inline constexpr r_bool_t r_true{1};
inline constexpr r_bool_t r_false{0};
inline constexpr r_bool_t r_na{std::numeric_limits<int>::min()};

inline bool is_r_true(r_bool_t x){
  return x.is_true();
}
inline bool is_r_false(r_bool_t x){
  return x.is_false();
}
// is_r_na is defined later as a template

// R integer 
struct r_int_t {
  int value;
  r_int_t() : value{0} {}
  explicit constexpr r_int_t(int x) : value{x} {}
  constexpr operator int() const { return value; }
};
// R double
struct r_double_t {
  double value;
  r_double_t() : value{0} {}
  explicit constexpr r_double_t(double x) : value{x} {}
  constexpr operator double() const { return value; }
};
// R integer64 (closely mimicking how bit64 defines it)
struct r_int64_t {
  int64_t value;
  r_int64_t() : value{0} {}
  explicit constexpr r_int64_t(int64_t x) : value{x} {}
  constexpr operator int64_t() const { return value; }
};

// Alias type for CHARSXP
struct r_string_t {
  SEXP value;
  r_string_t() : value{R_BlankString} {}
  // Explicit SEXP -> r_string_t
  explicit constexpr r_string_t(SEXP x) : value{x} {}
  // Implicit r_string_t -> SEXP
  constexpr operator SEXP() const { return value; }
};

inline const r_string_t blank_r_string = r_string_t();

// Alias type for SYMSXP
struct r_symbol_t {
  SEXP value;
  r_symbol_t() : value{R_MissingArg} {}
  explicit constexpr r_symbol_t(SEXP x) : value{x} {}
  constexpr operator SEXP() const { return value; }
};


// Alias type for Rcomplex
struct r_complex_t {
  Rcomplex value;

  // Constructors
  constexpr r_complex_t() : value{0.0, 0.0} {}
  constexpr r_complex_t(r_double_t r, r_double_t i) : value{r, i} {}

  // Conversion handling
  explicit constexpr r_complex_t(Rcomplex x) : value{x} {}
  constexpr operator Rcomplex() const { return value; }

  // Get real and imaginary parts
  constexpr r_double_t re() const { return r_double_t{value.r}; }
  constexpr r_double_t im() const { return r_double_t{value.i}; }
};

// Alias type for r_byte_t
struct r_byte_t {
  Rbyte value;

  // Constructors
  constexpr r_byte_t() : value{static_cast<Rbyte>(0)} {}

  // Conversion handling
  explicit constexpr r_byte_t(Rbyte x) : value{x} {}
  constexpr operator Rbyte() const { return value; }
};
