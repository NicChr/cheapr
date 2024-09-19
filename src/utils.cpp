#include "cheapr_cpp.h"

int int_div(int x, int y){
  return x / y;
}

[[cpp11::register]]
R_xlen_t cpp_vec_length(SEXP x){
  if (Rf_isFrame(x)){
    return cpp_df_nrow(x);
    // Is x a list?
  } else if (Rf_isVectorList(x)){
    if (Rf_inherits(x, "vctrs_rcrd")){
      return cpp_vec_length(VECTOR_ELT(x, 0));
    } else if (Rf_inherits(x, "POSIXlt")){
      const SEXP *p_x = VECTOR_PTR_RO(x);
      R_xlen_t out = 0;
      for (int i = 0; i != 10; ++i){
        out = std::max(out, Rf_xlength(p_x[i]));
      }
      return out;
      // return Rf_xlength(VECTOR_ELT(x, 0));
    } else if (Rf_isObject(x)){
      return Rf_asReal(cpp11::package("base")["length"](x));
    } else {
      return Rf_xlength(x);
    }
    // Catch-all
  } else {
    return Rf_xlength(x);
  }
}

int num_cores(){
  SEXP num_cores = Rf_protect(Rf_GetOption1(Rf_installChar(Rf_mkChar("cheapr.cores"))));
  int out = Rf_asInteger(num_cores);
  Rf_unprotect(1);
  return out >= 1 ? out : 1;
}

SEXP xlen_to_r(R_xlen_t x){
  return x > integer_max_ ? Rf_ScalarReal(x) : Rf_ScalarInteger(x);
}

R_xlen_t cpp_df_nrow(SEXP x){
  return Rf_xlength(Rf_getAttrib(x, R_RowNamesSymbol));
}

// Copy names from source to target
void cpp_copy_names(SEXP source, SEXP target){
  SEXP source_nms = Rf_protect(Rf_getAttrib(source, R_NamesSymbol));
  SEXP target_nms = Rf_protect(Rf_duplicate(source_nms));
  Rf_setAttrib(target, R_NamesSymbol, target_nms);
  Rf_unprotect(2);
}

SEXP r_address(SEXP x) {
  static char buf[1000];
  snprintf(buf, 1000, "%p", (void*) x);
  return Rf_mkChar(buf);
}

[[cpp11::register]]
SEXP r_copy(SEXP x){
  return Rf_duplicate(x);
}

bool is_int64(SEXP x){
  return Rf_isReal(x) && Rf_inherits(x, "integer64");
}

// We almost never want to convert back to 32-bit int

[[cpp11::register]]
SEXP cpp_int64_to_double(SEXP x){
  if (!is_int64(x)){
    Rf_error("x must be an integer64");
  }
  R_xlen_t n = Rf_xlength(x);

  SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
  double *p_out = REAL(out);

  long long *p_x = INTEGER64_PTR(x);

  for (R_xlen_t i = 0; i < n; ++i){
    p_out[i] = cheapr_is_na_int64(p_x[i]) ? NA_REAL : (double) p_x[i];
  }
  Rf_unprotect(1);
  return out;
}

// The reverse operation

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

    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = cheapr_is_na_int(p_x[i]) ? NA_INTEGER64 : (long long) p_x[i];
    }
    Rf_classgets(out, Rf_mkString("integer64"));
    Rf_unprotect(1);
    return out;
  }
  default: {
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
    long long *p_out = INTEGER64_PTR(out);
    double *p_x = REAL(x);

    for (R_xlen_t i = 0; i < n; ++i){
      p_out[i] = cheapr_is_na_dbl(p_x[i]) ? NA_INTEGER64 : (long long) p_x[i];
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

[[cpp11::register]]
SEXP cpp_format_double_as_int64(SEXP x){
  R_xlen_t n = Rf_xlength(x);

  SEXP out = Rf_protect(Rf_allocVector(STRSXP, n));
  double *p_x = REAL(x);
  SEXP na_char = Rf_protect(Rf_mkChar("NA"));
  for (R_xlen_t i = 0; i < n; ++i){
    if (cheapr_is_na_dbl(p_x[i])){
      SET_STRING_ELT(out, i, na_char);
    } else {
      long long temp = p_x[i];
      std::string s = string_format("%lld", temp);
      SET_STRING_ELT(out, i, Rf_mkChar(s.c_str()));
    }
  }
  Rf_unprotect(2);
  return out;
}

// Potentially useful for rolling calculations
// Computes the rolling number of true values in a given
// series of consecutive true values

// SEXP cpp_run_id(SEXP x, bool left_to_right){
//   if (!Rf_isLogical(x)){
//     Rf_error("x must be a logical vector");
//   }
//   R_xlen_t count = 0;
//   R_xlen_t n = Rf_xlength(x);
//   SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
//   int *p_out = INTEGER(out);
//   int *p_x = LOGICAL(x);
//   if (left_to_right){
//     for (R_xlen_t i = 0; i < n; ++i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   } else {
//     for (R_xlen_t i = n - 1; i >= 0; --i){
//       count = (count + p_x[i]) * p_x[i];
//       p_out[i] = count;
//     }
//   }
//   Rf_unprotect(1);
//   return out;
// }

// Would use data.table as it is very efficient, but would require extra dependency
// Here x must be an integer vector

// SEXP cpp_between(SEXP x, SEXP lower, SEXP upper){
//   int *p_x = INTEGER(x);
//   int n = Rf_length(x);
//   if (Rf_length(lower) != 1){
//     Rf_error("lower must be of length 1");
//   }
//   if (Rf_length(upper) != 1){
//     Rf_error("upper must be of length 1");
//   }
//   Rf_protect(lower = Rf_coerceVector(lower, INTSXP));
//   Rf_protect(upper = Rf_coerceVector(upper, INTSXP));
//   int lo = Rf_asInteger(lower);
//   int hi = Rf_asInteger(upper);
//   if (hi < lo){
//     int hi2 = hi;
//     hi = lo;
//     lo = hi2;
//   }
//   SEXP out = Rf_protect(Rf_allocVector(LGLSXP, n));
//   int *p_out = LOGICAL(out);
//   bool do_parallel = n >= 100000;
//   int n_cores = do_parallel ? num_cores() : 1;
//   if (lo == NA_INTEGER || hi == NA_INTEGER){
// #pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
//     for (int i = 0; i < n; ++i){
//       p_out[i] = NA_LOGICAL;
//     }
//   } else {
//     unsigned int rng = hi - lo;
// #pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
//     for (int i = 0; i < n; ++i){
//       int xi = p_x[i];
//       p_out[i] = (xi != NA_INTEGER) ? unsigned(xi - lo) <= rng : NA_LOGICAL;
//     }
//   }
//   Rf_unprotect(3);
//   return out;
// }
