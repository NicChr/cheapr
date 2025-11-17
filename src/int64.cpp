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

  const int64_t *p_x = INTEGER64_PTR_RO(x);

  int64_t int_max = INTEGER_MAX;

  for (R_xlen_t i = 0; i < n; ++i){
    p_out[i] = is_r_na(p_x[i]) || std::llabs(p_x[i]) > int_max ? NA_INTEGER : p_x[i];
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

  const int64_t *p_x = INTEGER64_PTR_RO(x);

  double repl;
  for (R_xlen_t i = 0; i < n; ++i){
    repl = is_r_na(p_x[i]) ? NA_REAL : p_x[i];
    p_out[i] = repl;
  }
  YIELD(1);
  return out;
}

// Can all numbers be safely converted to 32-bit int?

bool cpp_all_integerable(SEXP x){
  R_xlen_t n = Rf_xlength(x);
  bool out = true;

  switch (CHEAPR_TYPEOF(x)){
  case LGLSXP:
  case INTSXP: {
    break;
  }
  case CHEAPR_INT64SXP: {
    const int64_t *p_x = INTEGER64_PTR_RO(x);

    int64_t int_min = INTEGER_MIN;
    int64_t int_max = INTEGER_MAX;

    for (R_xlen_t i = 0; i < n; ++i){
      if (!(is_r_na(p_x[i]) || between(p_x[i], int_min, int_max))){
        out = false;
        break;
      }
    }
    break;
  }
  case REALSXP: {
    const double *p_x = REAL_RO(x);

    double int_min = INTEGER_MIN;
    double int_max = INTEGER_MAX;

    for (R_xlen_t i = 0; i < n; ++i){
      if (!(is_r_na(p_x[i]) || between(p_x[i], int_min, int_max))){
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
  return cpp_all_integerable(x) ? cpp_int64_to_int(x) : cpp_int64_to_double(x);
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
      p_out[i] = as_int64(p_x[i]);
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
      if (is_r_na(temp) || temp == R_PosInf || temp == R_NegInf){
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
      if (is_r_na(p_x[i])){
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
      if (is_r_na(p_x[i])){
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
      if (is_r_na(p_x[i])){
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


[[cpp11::register]]
SEXP cpp_sset_int64(SEXP x, SEXP locs){

  int32_t NP = 0;

  const int64_t* p_x = INTEGER64_PTR_RO(x);

  SEXP clean_locs = SHIELD(clean_indices(locs, x, false)); ++NP;
  SHIELD(locs = VECTOR_ELT(clean_locs, 0)); ++NP;

  SEXP out = SHIELD(new_vec(REALSXP, Rf_xlength(locs))); ++NP;
  int64_t* RESTRICT p_out = INTEGER64_PTR(out);

  SEXP names = SHIELD(get_names(x)); ++NP;
  SEXP out_names = SHIELD(sset_vec(names, locs, true)); ++NP;

  if (Rf_xlength(x) > INTEGER_MAX){

    int_fast64_t xn = Rf_xlength(x);

    int_fast64_t
    n = Rf_xlength(locs), k = 0, j;

    const double* pind = REAL_RO(locs);

    for (int_fast64_t i = 0; i < n; ++i){
      j = pind[i];
      if (j < 0){
        SEXP new_i = SHIELD(exclude_locs(locs, xn)); ++NP;
        SEXP out2 = SHIELD(cpp_sset_int64(x, new_i)); ++NP;
        YIELD(NP);
        return out2;
      } else if (j != 0){
        p_out[k++] = (is_r_na(pind[i]) || j > xn) ? NA_INTEGER64 : p_x[j - 1];
      }
    }

    // Resize if necessary (only when locs contains zeroes)

    if (k != n){
      SHIELD(out = Rf_xlengthgets(out, k)); ++NP;
    }

    set_names(out, out_names);

    YIELD(NP);
    return out;

  } else {

    unsigned int
    xn = Rf_length(x),
      n = Rf_xlength(locs),
      k = 0,
      na_val = NA_INTEGER,
      j;

    const int *pind = INTEGER_RO(locs);

    for (unsigned int i = 0; i < n; ++i){
      j = pind[i];
      if (between(j, 1U, xn)){
        p_out[k++] = p_x[--j];
        // If j > n_val then it is a negative 32-bit integer
      } else if (j > na_val){
        SEXP new_i = SHIELD(exclude_locs(locs, xn)); ++NP;
        SEXP out2 = SHIELD(cpp_sset_int64(x, new_i)); ++NP;
        YIELD(NP);
        return out2;
      } else if (j != 0U){
        p_out[k++] = NA_INTEGER64;
      }
    }

    if (k != n){
      SHIELD(out = Rf_lengthgets(out, k)); ++NP;
    }

    set_names(out, out_names);

    YIELD(NP);
    return out;
  }
}
