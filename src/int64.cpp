#include "cheapr.h"

// Convert 64-bit integer vec to 32-bit integer vec

[[cpp11::register]]
SEXP cpp_int64_to_int(SEXP x){
  if (!is_int64(x)){
    Rf_error("x must be an integer64");
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = SHIELD(new_vector<int>(n));
  int* RESTRICT p_out = integer_ptr(out);

  const int64_t *p_x = integer64_ptr_ro(x);

  int64_t int_max = r_limits::r_int_max;

  for (R_xlen_t i = 0; i < n; ++i){
    p_out[i] = is_r_na(p_x[i]) || std::llabs(p_x[i]) > int_max ? na::integer : p_x[i];
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

  SEXP out = SHIELD(new_vector<double>(n));
  double* RESTRICT p_out = real_ptr(out);

  const int64_t *p_x = integer64_ptr_ro(x);

  for (R_xlen_t i = 0; i < n; ++i){
    p_out[i] = r_cast<double>(p_x[i]);
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
    const int64_t *p_x = integer64_ptr_ro(x);

    int64_t int_min = r_limits::r_int_min;
    int64_t int_max = r_limits::r_int_max;

    for (R_xlen_t i = 0; i < n; ++i){
      if (!(is_r_na(p_x[i]) || between(p_x[i], int_min, int_max))){
        out = false;
        break;
      }
    }
    break;
  }
  case REALSXP: {
    const double *p_x = real_ptr_ro(x);

    double int_min = r_limits::r_int_min;
    double int_max = r_limits::r_int_max;

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

  int32_t NP = 0;
  R_xlen_t n = Rf_xlength(x);

  SEXP out = r_null;

  switch (CHEAPR_TYPEOF(x)){
  case INTSXP: {
    int *p_x = integer_ptr(x);
    out = SHIELD(new_vector<int64_t>(n)); ++NP;
    int64_t *p_out = integer64_ptr(out);
    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = r_cast<int64_t>(p_x[i]);
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    out = x;
    break;
  }
  case REALSXP: {
    double *p_x = real_ptr(x);
    out = SHIELD(new_vector<int64_t>(n)); ++NP;
    int64_t *p_out = integer64_ptr(out);
    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = r_cast<int64_t>(p_x[i]);
    }
    break;
  }
  default: {
    YIELD(NP);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
  YIELD(NP);
  return out;
}

// Found here stackoverflow.com/questions/347949
template<typename ... Args>
std::string string_format(const char *format, Args ... args){
  int size_s = std::snprintf( nullptr, 0, format, args ... ) + 1; // Extra space for '\0'
  if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
  auto size = static_cast<size_t>( size_s );
  std::unique_ptr<char[]> buf( new char[ size ] );
  std::snprintf( buf.get(), size, format, args ... );
  return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

// Simple format large integers
// same as `format(x, scientific = FALSE, trim = TRUE)`
// Except NA values become NA_character_

[[cpp11::register]]
SEXP cpp_format_numeric_as_int64(SEXP x){

  R_xlen_t n = Rf_xlength(x);

  SEXP out;

  switch (CHEAPR_TYPEOF(x)){
  case INTSXP: {
    std::string s;
    out = SHIELD(new_vector<r_string_t>(n));
    int *p_x = integer_ptr(x);

    for (R_xlen_t i = 0; i < n; ++i){
      if (is_r_na(p_x[i])){
        set_value<r_string_t>(out, i, na::string);
      } else {
        int64_t temp = p_x[i];
        s = string_format("%lld", temp);
        set_value<r_string_t>(out, i, r_cast<r_string_t>(s.c_str()));
      }
    }
    break;
  }
  case CHEAPR_INT64SXP: {
    std::string s;
    out = SHIELD(new_vector<r_string_t>(n));
    int64_t *p_x = integer64_ptr(x);

    for (R_xlen_t i = 0; i < n; ++i){
      if (is_r_na(p_x[i])){
        set_value<r_string_t>(out, i, na::string);
      } else {
        int64_t temp = p_x[i];
        s = string_format("%lld", temp);
        set_value<r_string_t>(out, i, r_cast<r_string_t>(s.c_str()));
      }
    }
    break;
  }
  case REALSXP: {
    std::string s;
    out = SHIELD(new_vector<r_string_t>(n));
    double *p_x = real_ptr(x);
    for (R_xlen_t i = 0; i < n; ++i){
      if (is_r_na(p_x[i])){
        set_value<r_string_t>(out, i, na::string);
      } else {
        int64_t temp = p_x[i];
        s = string_format("%lld", temp);
        set_value<r_string_t>(out, i, r_cast<r_string_t>(s.c_str()));
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

  const int64_t* p_x = integer64_ptr_ro(x);

  SEXP clean_locs = SHIELD(clean_indices(locs, x, false)); ++NP;
  SHIELD(locs = VECTOR_ELT(clean_locs, 0)); ++NP;

  SEXP out = SHIELD(new_vector<double>(Rf_xlength(locs))); ++NP;
  int64_t* RESTRICT p_out = integer64_ptr(out);

  SEXP names = SHIELD(get_old_names(x)); ++NP;
  SEXP out_names = SHIELD(sset_vec(names, locs, true)); ++NP;

  if (Rf_xlength(x) > r_limits::r_int_max){

    int_fast64_t xn = Rf_xlength(x);

    int_fast64_t
    n = Rf_xlength(locs), k = 0, j;

    const double* pind = real_ptr_ro(locs);

    for (int_fast64_t i = 0; i < n; ++i){
      j = pind[i];
      if (j < 0){
        SEXP new_i = SHIELD(exclude_locs(locs, xn)); ++NP;
        SEXP out2 = SHIELD(cpp_sset_int64(x, new_i)); ++NP;
        YIELD(NP);
        return out2;
      } else if (j != 0){
        p_out[k++] = (is_r_na(pind[i]) || j > xn) ? na::integer64 : p_x[j - 1];
      }
    }

    // Resize if necessary (only when locs contains zeroes)

    if (k != n){
      SHIELD(out = Rf_xlengthgets(out, k)); ++NP;
    }

    set_old_names(out, out_names);

    YIELD(NP);
    return out;

  } else {

    unsigned int
    xn = Rf_length(x),
      n = Rf_xlength(locs),
      k = 0,
      na_val = na::integer,
      j;

    const int *pind = integer_ptr_ro(locs);

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
        p_out[k++] = na::integer64;
      }
    }

    if (k != n){
      SHIELD(out = Rf_lengthgets(out, k)); ++NP;
    }

    set_old_names(out, out_names);

    YIELD(NP);
    return out;
  }
}
