#include "cheapr.h"

// Convert 64-bit integer vec to 32-bit integer vec

[[cpp11::register]]
SEXP cpp_int64_to_int(SEXP x){
  if (!is_int64(x)){
    Rf_error("x must be an integer64");
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = SHIELD(new_vec(INTSXP, n));
  int* RESTRICT p_out = INTEGER(out);

  const int64_t *p_x = INTEGER64_PTR(x);

  int64_t int_max = integer_max_;

  for (R_xlen_t i = 0; i < n; ++i){
    p_out[i] = is_na_int64(p_x[i]) || std::llabs(p_x[i]) > int_max ? NA_INTEGER : p_x[i];
  }
  YIELD(1);
  return out;
}

// Convert 64-bit integer vec to double vec

[[cpp11::register]]
SEXP cpp_int64_to_double(SEXP x){
  if (!is_int64(x)){
    Rf_error("x must be an integer64");
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = SHIELD(new_vec(REALSXP, n));
  double* RESTRICT p_out = REAL(out);

  const int64_t *p_x = INTEGER64_PTR(x);

  double repl;
  for (R_xlen_t i = 0; i < n; ++i){
    repl = is_na_int64(p_x[i]) ? NA_REAL : p_x[i];
    p_out[i] = repl;
  }
  YIELD(1);
  return out;
}

// Can all numbers be safely converted to 32-bit int?

bool cpp_all_integerable(SEXP x, int shift = 0){
  R_xlen_t n = Rf_xlength(x);
  bool out = true;

  switch (CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    break;
  }
  case CHEAPR_INT64SXP: {
    const int64_t *p_x = INTEGER64_PTR(x);
    int64_t int_max = integer_max_;
    int64_t shift_ = shift;
    for (R_xlen_t i = 0; i < n; ++i){
      if (!is_na_int64(p_x[i]) && ( (std::llabs(p_x[i]) + shift_) > int_max )){
        out = false;
        break;
      }
    }
    break;
  }
  case REALSXP: {
    const double *p_x = REAL(x);
    double int_max = integer_max_;
    double shift_ = shift;
    for (R_xlen_t i = 0; i < n; ++i){
      if (!is_na_dbl(p_x[i]) && ( (std::fabs(p_x[i]) + shift_) > int_max )){
        out = false;
        break;
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
  R_xlen_t n = Rf_xlength(x);

  SEXP out;
  int64_t repl;

  switch (CHEAPR_TYPEOF(x)){
  case INTSXP: {
    int *p_x = INTEGER(x);
    out = SHIELD(new_vec(REALSXP, n));
    int64_t *p_out = INTEGER64_PTR(out);
    for (R_xlen_t i = 0; i < n; ++i){
      repl = is_na_int(p_x[i]) ? NA_INTEGER64 : p_x[i];
      p_out[i] = repl;
    }
    Rf_classgets(out, make_utf8_str("integer64"));
    break;
  }
  case CHEAPR_INT64SXP: {
    out = SHIELD(x);
    break;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    out = SHIELD(new_vec(REALSXP, n));
    int64_t *p_out = INTEGER64_PTR(out);
    double temp;
    for (R_xlen_t i = 0; i < n; ++i){
      temp = p_x[i];
      if (is_na_dbl(temp) || temp == R_PosInf || temp == R_NegInf){
        repl = NA_INTEGER64;
      } else {
        repl = temp;
      }
      p_out[i] = repl;
    }
    Rf_classgets(out, make_utf8_str("integer64"));
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  YIELD(1);
  return out;
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

  SEXP out;
  std::string s;

  switch (CHEAPR_TYPEOF(x)){
  case INTSXP: {
    out = SHIELD(new_vec(STRSXP, n));
    int *p_x = INTEGER(x);

    for (R_xlen_t i = 0; i < n; ++i){
      if (is_na_int(p_x[i])){
        SET_STRING_ELT(out, i, NA_STRING);
      } else {
        int64_t temp = p_x[i];
        s = string_format("%lld", temp);
        SET_STRING_ELT(out, i, make_utf8_char(s.c_str()));
      }
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    out = SHIELD(new_vec(STRSXP, n));
    int64_t *p_x = INTEGER64_PTR(x);

    for (R_xlen_t i = 0; i < n; ++i){
      if (is_na_int64(p_x[i])){
        SET_STRING_ELT(out, i, NA_STRING);
      } else {
        int64_t temp = p_x[i];
        s = string_format("%lld", temp);
        SET_STRING_ELT(out, i, make_utf8_char(s.c_str()));
      }
    }
    break;
  }
  case REALSXP: {
    out = SHIELD(new_vec(STRSXP, n));
    double *p_x = REAL(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (is_na_dbl(p_x[i])){
        SET_STRING_ELT(out, i, NA_STRING);
      } else {
        int64_t temp = p_x[i];
        s = string_format("%lld", temp);
        SET_STRING_ELT(out, i, make_utf8_char(s.c_str()));
      }
    }
    break;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  YIELD(1);
  return out;
}
