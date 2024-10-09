#include "cheapr_cpp.h"

bool is_int64(SEXP x){
  return Rf_isReal(x) && Rf_inherits(x, "integer64");
}

// Convert 64-bit integer vec to 32-bit integer vec

[[cpp11::register]]
SEXP cpp_int64_to_int(SEXP x){
  if (!is_int64(x)){
    Rf_error("x must be an integer64");
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
  int *p_out = INTEGER(out);

  long long *p_x = INTEGER64_PTR(x);

  int repl;
  long long int int_max = integer_max_;
  for (R_xlen_t i = 0; i < n; ++i){
    repl = cheapr_is_na_int64(p_x[i]) || std::llabs(p_x[i]) > int_max ? NA_INTEGER : p_x[i];
    p_out[i] = repl;
  }
  Rf_unprotect(1);
  return out;
}

// Convert 64-bit integer vec to double vec

[[cpp11::register]]
SEXP cpp_int64_to_double(SEXP x){
  if (!is_int64(x)){
    Rf_error("x must be an integer64");
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
  double *p_out = REAL(out);

  long long *p_x = INTEGER64_PTR(x);

  double repl;
  for (R_xlen_t i = 0; i < n; ++i){
    repl = cheapr_is_na_int64(p_x[i]) ? NA_REAL : p_x[i];
    p_out[i] = repl;
  }
  Rf_unprotect(1);
  return out;
}

// Can all numbers be safely converted to 32-bit int?

bool cpp_all_integerable(SEXP x, int shift = 0){
  R_xlen_t n = Rf_xlength(x);
  bool out = true;

  switch (TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    break;
  }
  case REALSXP: {
    if (is_int64(x)){
    long long int *p_x = INTEGER64_PTR(x);
    long long int int_max = integer_max_;
    long long int shift_ = shift;
    for (R_xlen_t i = 0; i < n; ++i){
      if (!cheapr_is_na_int64(p_x[i]) && ( (std::llabs(p_x[i]) + shift_) > int_max )){
        out = false;
        break;
      }
    }
  } else {
    double *p_x = REAL(x);
    double int_max = integer_max_;
    double shift_ = shift;
    for (R_xlen_t i = 0; i < n; ++i){
      if (!cheapr_is_na_dbl(p_x[i]) && ( (std::fabs(p_x[i]) + shift_) > int_max )){
        out = false;
        break;
      }
    }
  }
  break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  return out;
}

// Convert 64-bit integer to 32-bit int if possible, otherwise double

[[cpp11::register]]
SEXP cpp_int64_to_numeric(SEXP x){
  if (!is_int64(x)){
    Rf_error("x must be an integer64");
  }
  return cpp_all_integerable(x, 0) ? cpp_int64_to_int(x) : cpp_int64_to_double(x);
}

// The reverse operation
// Convert any numeric vector into 64-integer vec
// Very much an internal function as cheapr never explicitly outputs
// a 64-bit integer

[[cpp11::register]]
SEXP cpp_numeric_to_int64(SEXP x){

  if (is_int64(x)){
    return x;
  }

  R_xlen_t n = Rf_xlength(x);

  switch (TYPEOF(x)){
  case INTSXP: {
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
    long long *p_out = INTEGER64_PTR(out);
    int *p_x = INTEGER(x);
    long long int repl;
    for (R_xlen_t i = 0; i < n; ++i){
      repl = cheapr_is_na_int(p_x[i]) ? NA_INTEGER64 : p_x[i];
      p_out[i] = repl;
    }
    Rf_classgets(out, Rf_mkString("integer64"));
    Rf_unprotect(1);
    return out;
  }
  default: {
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
    long long *p_out = INTEGER64_PTR(out);
    double *p_x = REAL(x);

    long long int repl;
    double temp;
    for (R_xlen_t i = 0; i < n; ++i){
      temp = p_x[i];
      if (cheapr_is_na_dbl(temp) || (temp == R_PosInf) || temp == R_NegInf){
        repl = NA_INTEGER64;
      } else {
        repl = temp;
      }
      p_out[i] = repl;
    }
    Rf_classgets(out, Rf_mkString("integer64"));
    Rf_unprotect(1);
    return out;
  }
  }
}

// Found here stackoverflow.com/questions/347949
template<typename ... Args>
std::string string_format( const std::string& format, Args ... args){
  int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
  if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
  auto size = static_cast<size_t>( size_s );
  std::unique_ptr<char[]> buf( new char[ size ] );
  std::snprintf( buf.get(), size, format.c_str(), args ... );
  return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

// Simple format large integers
// same as `format(x, scientific = FALSE, trim = TRUE)`
// Except NA values become NA_character_

[[cpp11::register]]
SEXP cpp_format_numeric_as_int64(SEXP x){
  R_xlen_t n = Rf_xlength(x);

  switch (TYPEOF(x)){
  case INTSXP: {
    SEXP out = Rf_protect(Rf_allocVector(STRSXP, n));
    int *p_x = INTEGER(x);

    for (R_xlen_t i = 0; i < n; ++i){
      if (cheapr_is_na_int(p_x[i])){
        SET_STRING_ELT(out, i, NA_STRING);
      } else {
        long long temp = p_x[i];
        std::string s = string_format("%lld", temp);
        SET_STRING_ELT(out, i, Rf_mkChar(s.c_str()));
      }
    }
    Rf_unprotect(1);
    return out;
  }
  case REALSXP: {
    SEXP out = Rf_protect(Rf_allocVector(STRSXP, n));
    if (is_int64(x)){
      long long *p_x = INTEGER64_PTR(x);

      for (R_xlen_t i = 0; i < n; ++i){
        if (cheapr_is_na_int64(p_x[i])){
          SET_STRING_ELT(out, i, NA_STRING);
        } else {
          long long temp = p_x[i];
          std::string s = string_format("%lld", temp);
          SET_STRING_ELT(out, i, Rf_mkChar(s.c_str()));
        }
      }
    } else {
      double *p_x = REAL(x);

      for (R_xlen_t i = 0; i < n; ++i){
        if (cheapr_is_na_dbl(p_x[i])){
          SET_STRING_ELT(out, i, NA_STRING);
        } else {
          long long temp = p_x[i];
          std::string s = string_format("%lld", temp);
          SET_STRING_ELT(out, i, Rf_mkChar(s.c_str()));
        }
      }
    }
    Rf_unprotect(1);
    return out;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}
