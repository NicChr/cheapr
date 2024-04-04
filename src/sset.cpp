#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>
// #include <vector>
// using namespace cpp11;

// SEXP cpp_sset(SEXP x, SEXP indices){
//   int *pi = INTEGER(indices);
//   int xn = Rf_xlength(x);
//   int n = Rf_xlength(indices);
//   int n_protections = 0;
//   int zero_count = 0;
//   int pos_count = 0;
//   int oob_count = 0;
//   int na_count = 0;
//   int out_size;
//   bool do_parallel = n >= 10000;
//   int n_cores = do_parallel ? num_cores() : 1;
//   do_parallel = do_parallel && n_cores > 1;
//
//   // Counting the number of:
//   // Zeroes
//   // Out-of-bounds indices
//   // Positive indices
//   // From this we can also work out the number of negatives
//
//   if (do_parallel){
// #pragma omp parallel for simd num_threads(n_cores) reduction(+:zero_count,pos_count,oob_count,na_count)
//     for (int j = 0; j < n; ++j){
//       zero_count += pi[j] == 0;
//       pos_count += pi[j] > 0;
//       na_count += pi[j] == NA_INTEGER;
//       oob_count += std::abs(pi[j]) > xn;
//     }
//   } else {
// #pragma omp for simd
//     for (int j = 0; j < n; ++j){
//       zero_count += (pi[j] == 0);
//       pos_count += (pi[j] > 0);
//       na_count += pi[j] == NA_INTEGER;
//       oob_count += (std::abs(pi[j]) > xn);
//     }
//   }
//   bool neg_count = n - pos_count - zero_count - na_count;
//   if ( (pos_count + zero_count) > 0 && neg_count > 0){
//     Rf_error("Cannot mix positive and negative indices");
//   }
//   bool simple_sset = zero_count == 0 && oob_count == 0 && na_count == 0 && pos_count == n;
//
//   // Convert negative index vector to positive
//
//   if (neg_count > 0){
//     SEXP indices2 = Rf_protect(cpp11::package("cheapr")["neg_indices_to_pos"](indices, xn));
//     ++n_protections;
//     int *pi2 = INTEGER(indices2);
//     pi = pi2;
//     out_size = Rf_xlength(indices2);
//     n = out_size;
//     simple_sset = true;
//   } else {
//     out_size = n - zero_count;
//   }
//   switch ( TYPEOF(x) ){
//   int i;
//   case NILSXP: {
//     return R_NilValue;
//   }
//   case LGLSXP: {
//     int *p_x = LOGICAL(x);
//     SEXP out = Rf_protect(Rf_allocVector(LGLSXP, out_size));
//     ++n_protections;
//     zero_count = 0;
//     int *p_out = LOGICAL(out);
//     if (simple_sset){
//       if (do_parallel){
// #pragma omp parallel for simd num_threads(n_cores) private(i)
//         for (i = 0; i < n; ++i){
//           p_out[i] = p_x[pi[i] - 1];
//         }
//       } else {
// #pragma omp for simd
//         for (i = 0; i < n; ++i){
//           p_out[i] = p_x[pi[i] - 1];
//         }
//       }
//     } else {
//       for (i = 0; i < n; ++i){
//         if (pi[i] == 0){
//           ++zero_count;
//         } else {
//           p_out[i - zero_count] = (pi[i] <= xn && pi[i] != NA_INTEGER) ? p_x[pi[i] - 1] : NA_LOGICAL;
//         }
//         // p_out[i - zero_count] = (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_LOGICAL;
//         // p_out[i - ( pi[i] == 0 ? zero_count++ : zero_count)] = (pi[i] > 0 && pi[i] <= xn) ? p_x[pi[i] - 1] : NA_INTEGER;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case INTSXP: {
//     int *p_x = INTEGER(x);
//     SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
//     ++n_protections;
//     zero_count = 0;
//     int *p_out = INTEGER(out);
//     if (simple_sset){
//       if (do_parallel){
// #pragma omp parallel for simd num_threads(n_cores) private(i)
//         for (i = 0; i < n; ++i){
//           p_out[i] = p_x[pi[i] - 1];
//         }
//       } else {
// #pragma omp for simd
//         for (i = 0; i < n; ++i){
//           p_out[i] = p_x[pi[i] - 1];
//         }
//       }
//     } else {
//       for (i = 0; i < n; ++i){
//         if (pi[i] == 0){
//           ++zero_count;
//         } else {
//           p_out[i - zero_count] = (pi[i] <= xn && pi[i] != NA_INTEGER) ? p_x[pi[i] - 1] : NA_INTEGER;
//         }
//         // p_out[i - zero_count] = (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_INTEGER;
//         // p_out[i - ( pi[i] == 0 ? zero_count++ : zero_count)] = (pi[i] > 0 && pi[i] <= xn) ? p_x[pi[i] - 1] : NA_INTEGER;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case REALSXP: {
//     double *p_x = REAL(x);
//     SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
//     ++n_protections;
//     zero_count = 0;
//     double *p_out = REAL(out);
//     if (simple_sset){
//       if (do_parallel){
// #pragma omp parallel for simd num_threads(n_cores) private(i)
//         for (i = 0; i < n; ++i){
//           p_out[i] = p_x[pi[i] - 1];
//         }
//       } else {
// #pragma omp for simd
//         for (i = 0; i < n; ++i){
//           p_out[i] = p_x[pi[i] - 1];
//         }
//       }
//     } else {
//       for (i = 0; i < n; ++i){
//         if (pi[i] == 0){
//           ++zero_count;
//         } else {
//           p_out[i - zero_count] = (pi[i] <= xn && pi[i] != NA_INTEGER) ? p_x[pi[i] - 1] : NA_REAL;
//         }
//         // p_out[i - ( pi[i] == 0 ? zero_count++ : zero_count)] = (pi[i] > 0 && pi[i] <= xn) ? p_x[pi[i] - 1] : NA_REAL;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case STRSXP: {
//     SEXP *p_x = STRING_PTR(x);
//     SEXP out = Rf_protect(Rf_allocVector(STRSXP, out_size));
//     ++n_protections;
//     zero_count = 0;
//     if (simple_sset){
//       for (i = 0; i < n; ++i){
//         SET_STRING_ELT(out, i, p_x[pi[i] - 1]);
//       }
//     } else {
//       for (i = 0; i < n; ++i){
//         if (pi[i] == 0){
//           ++zero_count;
//         } else {
//           SET_STRING_ELT(out, i - zero_count,
//                          (pi[i] <= xn && pi[i] != NA_INTEGER) ? p_x[pi[i] - 1] : NA_STRING);
//         }
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case RAWSXP: {
//     Rbyte *p_x = RAW(x);
//     SEXP out = Rf_protect(Rf_allocVector(RAWSXP, out_size));
//     ++n_protections;
//     zero_count = 0;
//     if (simple_sset){
//       for (i = 0; i < n; ++i){
//         SET_RAW_ELT(out, i, p_x[pi[i] - 1]);
//       }
//     } else {
//       for (i = 0; i < n; ++i){
//         if (pi[i] == 0){
//           ++zero_count;
//         } else {
//           SET_RAW_ELT(out, i - zero_count,
//                       (pi[i] <= xn && pi[i] != NA_INTEGER) ? p_x[pi[i] - 1] : 0);
//         }
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   case VECSXP: {
//     const SEXP *p_x = VECTOR_PTR_RO(x);
//     SEXP out = Rf_protect(Rf_allocVector(VECSXP, out_size));
//     ++n_protections;
//     zero_count = 0;
//     if (simple_sset){
//       for (i = 0; i < n; ++i){
//         SET_VECTOR_ELT(out, i, p_x[pi[i] - 1]);
//       }
//     } else {
//       for (i = 0; i < n; ++i){
//         if (pi[i] == 0){
//           ++zero_count;
//         } else {
//           SET_VECTOR_ELT(out, i - zero_count,
//                          (pi[i] <= xn && pi[i] != NA_INTEGER) ? p_x[pi[i] - 1] : R_NilValue);
//         }
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
// }

// Altrep utils

bool is_altrep(SEXP x){
  return ALTREP(x);
}
SEXP alt_class(SEXP x){
  if (is_altrep(x)){
    SEXP alt_attribs = Rf_protect(Rf_coerceVector(ATTRIB(ALTREP_CLASS(x)), VECSXP));
    SEXP out = Rf_protect(Rf_coerceVector(VECTOR_ELT(alt_attribs, 0), STRSXP));
    Rf_unprotect(2);
    return out;
  } else {
    return R_NilValue;
  }
}
SEXP alt_pkg(SEXP x){
  if (is_altrep(x)){
    SEXP alt_attribs = Rf_protect(Rf_coerceVector(ATTRIB(ALTREP_CLASS(x)), VECSXP));
    SEXP out = Rf_protect(Rf_coerceVector(VECTOR_ELT(alt_attribs, 1), STRSXP));
    Rf_unprotect(2);
    return out;
  } else {
    return R_NilValue;
  }
}
[[cpp11::register]]
SEXP alt_data1(SEXP x){
  if (is_altrep(x)){
    return R_altrep_data1(x);
  } else {
    return R_NilValue;
  }
}
[[cpp11::register]]
bool is_alt_int_seq(SEXP x){
  if (!is_altrep(x)) return false;
  SEXP alt_class_nm = Rf_protect(alt_class(x));
  SEXP alt_pkg_nm = Rf_protect(alt_pkg(x));
  SEXP intseq_char = Rf_protect(Rf_mkChar("compact_intseq"));
  SEXP base_char = Rf_protect(Rf_mkChar("base"));
  bool out = STRING_ELT(alt_class_nm, 0) == intseq_char &&
    STRING_ELT(alt_pkg_nm, 0) == base_char;
  Rf_unprotect(4);
  return out;
}
double alt_int_seq_start(SEXP x){
  if (!is_alt_int_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = Rf_protect(Rf_coerceVector(alt_data1(x), REALSXP));
  double out = REAL(alt_data)[1];
  Rf_unprotect(1);
  return out;
}
double alt_int_seq_increment(SEXP x){
  if (!is_alt_int_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = Rf_protect(Rf_coerceVector(alt_data1(x), REALSXP));
  double out = REAL(alt_data)[2];
  Rf_unprotect(1);
  return out;
}

double alt_int_seq_end(SEXP x){
  if (!is_alt_int_seq(x)){
    Rf_error("x must be an altrep compact_intseq");
  }
  SEXP alt_data = Rf_protect(Rf_coerceVector(alt_data1(x), REALSXP));
  double alt_size = REAL(alt_data)[0];
  double alt_from = REAL(alt_data)[1];
  double alt_increment = REAL(alt_data)[2];
  double out = (std::fmax(alt_size - 1.0, 0.0) * alt_increment) + alt_from;
  Rf_unprotect(1);
  return out;
}

// Subset with no checks, indices vector must be pre-curated
[[cpp11::register]]
SEXP cpp_sset_simple(SEXP x, SEXP indices){
  int *pi = INTEGER(indices);
  int n = Rf_xlength(indices);
  int n_protections = 0;
  int out_size = Rf_length(indices);
  bool do_parallel = n >= 100000;
  int n_cores = do_parallel ? num_cores() : 1;

#define SIMPLE_SSET                                            \
  for (i = 0; i < n; ++i){                                     \
    p_out[i] = p_x[pi[i] - 1];                                 \
  }

switch ( TYPEOF(x) ){

int i;
case NILSXP: {
  return R_NilValue;
}
case LGLSXP: {
  int *p_x = LOGICAL(x);
  SEXP out = Rf_protect(Rf_allocVector(LGLSXP, out_size));
  ++n_protections;
  int *p_out = LOGICAL(out);
  if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores)
    SIMPLE_SSET;
  } else {
#pragma omp for simd
    SIMPLE_SSET;
  }
  Rf_unprotect(n_protections);
  return out;
}
case INTSXP: {
  int *p_x = INTEGER(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
  ++n_protections;
  int *p_out = INTEGER(out);
  if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores)
    SIMPLE_SSET;
  } else {
#pragma omp for simd
SIMPLE_SSET;
}
  for (i = 0; i < n; ++i){
    p_out[i] = p_x[pi[i] - 1];
  }
  Rf_unprotect(n_protections);
  return out;
}
case REALSXP: {
  double *p_x = REAL(x);
  SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
  ++n_protections;
  double *p_out = REAL(out);
  if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores)
    SIMPLE_SSET;
  } else {
#pragma omp for simd
    SIMPLE_SSET;
  }
  Rf_unprotect(n_protections);
  return out;
}
case STRSXP: {
  SEXP *p_x = STRING_PTR(x);
  SEXP out = Rf_protect(Rf_allocVector(STRSXP, out_size));
  ++n_protections;
#pragma omp for simd
  for (i = 0; i < n; ++i){
    SET_STRING_ELT(out, i, p_x[pi[i] - 1]);
  }
  Rf_unprotect(n_protections);
  return out;
}
case CPLXSXP: {
  Rcomplex *p_x = COMPLEX(x);
  SEXP out = Rf_protect(Rf_allocVector(CPLXSXP, out_size));
  ++n_protections;
#pragma omp for simd
  for (i = 0; i < n; ++i){
    SET_COMPLEX_ELT(out, i, p_x[pi[i] - 1]);
  }
  Rf_unprotect(n_protections);
  return out;
}
case RAWSXP: {
  Rbyte *p_x = RAW(x);
  SEXP out = Rf_protect(Rf_allocVector(RAWSXP, out_size));
  ++n_protections;
#pragma omp for simd
  for (i = 0; i < n; ++i){
    SET_RAW_ELT(out, i, p_x[pi[i] - 1]);
  }
  Rf_unprotect(n_protections);
  return out;
}
case VECSXP: {
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, out_size));
  ++n_protections;
#pragma omp for simd
  for (i = 0; i < n; ++i){
    SET_VECTOR_ELT(out, i, p_x[pi[i] - 1]);
  }
  Rf_unprotect(n_protections);
  return out;
}
default: {
  Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
}
}
}

// A range-based subset method
// Can be readily used when indices are an altrep compact integer sequence
// Also works with negative-indexing

[[cpp11::register]]
SEXP cpp_sset_range(SEXP x, int from, int to, int by){
  int n = Rf_length(x);
  if (by != 1 && by != -1){
    Rf_error("by increment must be 1 or -1");
  }
  if ( (from > 0 && to < 0) ||
       (from < 0 && to > 0) ){
    Rf_error("Cannot mix positive and negative indices");
  }
  if ( (from > to && by > 0) || (from < to && by < 0)){
    Rf_error("Wrong increment sign in by arg");
  }
  int istart = from;
  int iend = to;
  int out_size, istart1, istart2, iend1, iend2;
  bool double_loop = false;
  // Negative indexing is complicated
  // Because there are 7 scenarios...
  // Scenario 1 and 2 are out-of-bound scenarios
  // Scenario 3 (prob jsut falls back to scenario 5): n:n - We just exclude x[n]
  // Scenario 4: -1:length(x) - We return an empty vector
  // Scenario 5: -1:nE{x:x< length(x)} - We subset from n+1:length(x)
  // Scenario 6: m:nE{n = length(x)} - We subset from m:length(x)
  // Scenario 7: m:kE{m > 1 & k < length(x)} - This requires 2 loops
  if (istart == 0 && iend == 0){
    out_size = 0;
  } else if (istart < 0 || iend < 0){
    if (istart == 0){
      istart = -1;
    }
    if (iend == 0){
      iend = -1;
    }
    // We first switch them
    if (istart < iend){
      int iend_temp = iend;
      iend = istart;
      istart = iend_temp;
    }
    // Out-of-bounds adjustments
    if (std::abs(istart) <= n && std::abs(iend) > n){
      iend = std::abs(istart) - 1;
      istart = 1;
      by = 1;
      out_size = (iend - istart) + 1;
    } else if (std::abs(istart) > n && std::abs(iend) > n){
      istart = 1;
      iend  = n;
      by = 1;
      out_size = n;
    } else if (istart == -1 && iend == -n){
      istart = n;
      iend = n;
      by = 1;
      out_size = 0;
      // Scenario 3
    } else if (istart == -1 && std::abs(iend) < n){
      istart = std::abs(iend) + 1;
      iend = n;
      by = 1;
      out_size = (iend - istart) + 1;
      // Scenario 4
    } else if (std::abs(istart) < n && std::abs(iend) == n){
      iend = std::abs(istart) - 1;
      istart = 1;
      by = 1;
      out_size = (iend - istart) + 1;
    } else {
      // Scenario 7
      double_loop = true;
      istart1 = 1;
      iend1 = std::abs(istart) - 1;
      istart2 = std::abs(iend) + 1;
      iend2 = n;
      by = 1;
      out_size = (iend1 - istart1) + (iend2 - istart2) + 2;
    }
  } else {
    if (istart == 0){
      istart = 1;
    }
    if (iend == 0){
      iend = 1;
    }
    out_size = ((iend - istart) * by) + 1;
  }
  bool do_parallel = out_size >= 100000;
  int n_cores = do_parallel ? num_cores() : 1;
  unsigned int k = 0;

switch ( TYPEOF(x) ){
case NILSXP: {
  return R_NilValue;
}
case LGLSXP: {
  int *p_x = LOGICAL(x);
  SEXP out = Rf_protect(Rf_allocVector(LGLSXP, out_size));
  int *p_out = LOGICAL(out);
  if (double_loop){
#pragma omp parallel num_threads(n_cores) if (do_parallel)
#pragma omp for simd
    for (int i = istart1 - 1; i < iend1; ++i){
      p_out[i - istart1 + 1] = i < n ? p_x[i] : NA_LOGICAL;
    }
#pragma omp for simd
    for (int i = istart2 - 1; i < iend2; ++i){
      p_out[i - istart2 + 1] = i < n ? p_x[i] : NA_LOGICAL;
    }
  } else {
    if (by > 0){
#pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
      for (int i = istart - 1; i < iend; ++i){
        p_out[i - istart + 1] = i < n ? p_x[i] : NA_LOGICAL;
      }
    } else {
#pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
      for (int i = istart - 1; i >= iend - 1; --i){
        p_out[istart - i - 1] = i < n ? p_x[i] : NA_LOGICAL;
      }
    }
  }
  Rf_unprotect(1);
  return out;
}
case INTSXP: {
  int *p_x = INTEGER(x);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
  int *p_out = INTEGER(out);
  if (double_loop){
#pragma omp parallel num_threads(n_cores) if (do_parallel)
#pragma omp for simd
    for (int i = istart1 - 1; i < iend1; ++i){
      p_out[i - istart1 + 1] = i < n ? p_x[i] : NA_INTEGER;
    }
#pragma omp for simd
    for (int i = istart2 - 1; i < iend2; ++i){
      p_out[i - istart2 + 1] = i < n ? p_x[i] : NA_INTEGER;
    }
  } else {
    if (by > 0){
#pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
      for (int i = istart - 1; i < iend; ++i){
        p_out[i - istart + 1] = i < n ? p_x[i] : NA_INTEGER;
      }
    } else {
#pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
      for (int i = istart - 1; i >= iend - 1; --i){
        p_out[istart - i - 1] = i < n ? p_x[i] : NA_INTEGER;
      }
    }
  }
  Rf_unprotect(1);
  return out;
}
case REALSXP: {
  double *p_x = REAL(x);
  SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
  double *p_out = REAL(out);
  if (double_loop){
#pragma omp parallel num_threads(n_cores) if (do_parallel)
#pragma omp for simd
    for (int i = istart1 - 1; i < iend1; ++i){
      p_out[i - istart1 + 1] = i < n ? p_x[i] : NA_REAL;
    }
#pragma omp for simd
    for (int i = istart2 - 1; i < iend2; ++i){
      p_out[i - istart2 + 1] = i < n ? p_x[i] : NA_REAL;
    }
  } else {
    if (by > 0){
#pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
      for (int i = istart - 1; i < iend; ++i){
        p_out[i - istart + 1] = i < n ? p_x[i] : NA_REAL;
      }
    } else {
#pragma omp parallel for simd num_threads(n_cores) if (do_parallel)
      for (int i = istart - 1; i >= iend - 1; --i){
        p_out[istart - i - 1] = i < n ? p_x[i] : NA_REAL;
      }
    }
  }
  Rf_unprotect(1);
  return out;
}
case STRSXP: {
  SEXP *p_x = STRING_PTR(x);
  SEXP out = Rf_protect(Rf_allocVector(STRSXP, out_size));
  if (double_loop){
    for (int i = istart1 - 1; i < iend1; ++i){
      SET_STRING_ELT(out, k++, i < n ? p_x[i] : NA_STRING);
    }
    for (int j = istart2 - 1; j < iend2; ++j){
      SET_STRING_ELT(out, k++, j < n ? p_x[j] : NA_STRING);
    }
  } else {
    if (by > 0){
      for (int i = istart - 1; i < iend; ++i){
        SET_STRING_ELT(out, k++, i < n ? p_x[i] : NA_STRING);
      }
    } else {
      for (int i = istart - 1; i >= iend - 1; --i){
        SET_STRING_ELT(out, k++, i < n ? p_x[i] : NA_STRING);
      }
    }
  }
  Rf_unprotect(1);
  return out;
}
case CPLXSXP: {
  Rcomplex *p_x = COMPLEX(x);
  SEXP out = Rf_protect(Rf_allocVector(CPLXSXP, out_size));
  SEXP na_complex_sexp = Rf_protect(Rf_allocVector(CPLXSXP, 1));
  Rcomplex *p_na_complex = COMPLEX(na_complex_sexp);
  p_na_complex[0].i = NA_REAL;
  p_na_complex[0].r = NA_REAL;
  Rcomplex na_complex = Rf_asComplex(na_complex_sexp);
  if (double_loop){
    for (int i = istart1 - 1; i < iend1; ++i){
      SET_COMPLEX_ELT(out, k++, i < n ? p_x[i] : na_complex);
    }
    for (int j = istart2 - 1; j < iend2; ++j){
      SET_COMPLEX_ELT(out, k++, j < n ? p_x[j] : na_complex);
    }
  } else {
    if (by > 0){
      for (int i = istart - 1; i < iend; ++i){
        SET_COMPLEX_ELT(out, k++, i < n ? p_x[i] : na_complex);
      }
    } else {
      for (int i = istart - 1; i >= iend - 1; --i){
        SET_COMPLEX_ELT(out, k++, i < n ? p_x[i] : na_complex);
      }
    }
  }
  Rf_unprotect(2);
  return out;
}
case RAWSXP: {
  Rbyte *p_x = RAW(x);
  SEXP out = Rf_protect(Rf_allocVector(RAWSXP, out_size));
  if (double_loop){
    for (int i = istart1 - 1; i < iend1; ++i){
      SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
    }
    for (int j = istart2 - 1; j < iend2; ++j){
      SET_RAW_ELT(out, k++, j < n ? p_x[j] : 0);
    }
  } else {
    if (by > 0){
      for (int i = istart - 1; i < iend; ++i){
        SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
      }
    } else {
      for (int i = istart - 1; i >= iend - 1; --i){
        SET_RAW_ELT(out, k++, i < n ? p_x[i] : 0);
      }
    }
  }
  Rf_unprotect(1);
  return out;
}
case VECSXP: {
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, out_size));
  if (double_loop){
    for (int i = istart1 - 1; i < iend1; ++i){
      if (i < n) SET_VECTOR_ELT(out, k++, p_x[i]);
    }
    for (int j = istart2 - 1; j < iend2; ++j){
      if (j < n) SET_VECTOR_ELT(out, k++, p_x[j]);
    }
  } else {
    if (by > 0){
      for (int i = istart - 1; i < iend; ++i){
        if (i < n) SET_VECTOR_ELT(out, k++, p_x[i]);
      }
    } else {
      for (int i = istart - 1; i >= iend - 1; --i){
        if (i < n) SET_VECTOR_ELT(out, k++, p_x[i]);
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

[[cpp11::register]]
SEXP cpp_sset_df(SEXP x, SEXP indices){
  int *pi = INTEGER(indices);
  int xn = cpp_df_nrow(x);
  int ncols = Rf_length(x);
  int n = Rf_xlength(indices);
  int n_protections = 0;
  int zero_count = 0;
  int pos_count = 0;
  int oob_count = 0;
  int na_count = 0;
  bool do_parallel = n >= 10000;
  int n_cores = do_parallel ? num_cores() : 1;
  // Counting the number of:
  // Zeroes
  // Out-of-bounds indices
  // Positive indices
  // From this we can also work out the number of negatives

  // If indices is a special type of ALTREP compact int sequence, we can infer
  // this information
  if (is_alt_int_seq(indices) &&
      std::fabs(alt_int_seq_increment(indices)) == 1.0 &&
      alt_int_seq_start(indices) > 0 &&
      alt_int_seq_start(indices) <= xn &&
      alt_int_seq_end(indices) > 0 &&
      alt_int_seq_end(indices) <= xn){
    zero_count = 0;
    na_count = 0;
    oob_count = 0;
    pos_count = n;
  } else {
    if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:zero_count,pos_count,oob_count,na_count)
      for (int j = 0; j < n; ++j){
        zero_count += (pi[j] == 0);
        pos_count += (pi[j] > 0);
        oob_count += (std::fabs(pi[j]) > xn);
        na_count += (pi[j] == NA_INTEGER);
      }
    } else {
#pragma omp for simd
      for (int j = 0; j < n; ++j){
        zero_count += (pi[j] == 0);
        pos_count += (pi[j] > 0);
        oob_count += (std::fabs(pi[j]) > xn);
        na_count += (pi[j] == NA_INTEGER);
      }
    }
  }
  int neg_count = n - pos_count - zero_count - na_count;
  if ( (pos_count > 0 && neg_count > 0) ||
       (neg_count > 0 && na_count > 0)){
    Rf_error("Cannot mix positive and negative indices");
  }
  bool simple_sset = zero_count == 0 && oob_count == 0 && na_count == 0 && pos_count == n;
  cpp11::function cheapr_sset = cpp11::package("cheapr")["sset"];

  // Convert negative index vector to positive
  const SEXP *p_x = VECTOR_PTR_RO(x);
  SEXP out = Rf_protect(Rf_allocVector(VECSXP, ncols));
  ++n_protections;
  // Index vector is clean, we can use fast subset
  if (simple_sset){
    for (int j = 0; j < ncols; ++j){
      if (Rf_isObject(p_x[j])){
        SET_VECTOR_ELT(out, j, cheapr_sset(p_x[j], indices));
      } else {
        SET_VECTOR_ELT(out, j, cpp_sset_simple(p_x[j], indices));
      }
    }
    // Negative indexing
  } else if (neg_count > 0){
    SEXP indices2 = Rf_protect(cpp11::package("cheapr")["neg_indices_to_pos"](indices, xn));
    ++n_protections;
    for (int j = 0; j < ncols; ++j){
      if (Rf_isObject(p_x[j])){
        SET_VECTOR_ELT(out, j, cheapr_sset(p_x[j], indices2));
      } else {
        SET_VECTOR_ELT(out, j, cpp_sset_simple(p_x[j], indices2));
      }
    }
    // If index vector is clean except for existence of zeroes
  } else if (zero_count > 0 && oob_count == 0 && na_count == 0){
    SEXP r_zero = Rf_protect(Rf_ScalarInteger(0));
    ++n_protections;
    SEXP indices2 = Rf_protect(cpp11::package("cheapr")["val_rm"](indices, r_zero));
    ++n_protections;
    for (int j = 0; j < ncols; ++j){
      if (Rf_isObject(p_x[j])){
        SET_VECTOR_ELT(out, j, cheapr_sset(p_x[j], indices2));
      } else {
        SET_VECTOR_ELT(out, j, cpp_sset_simple(p_x[j], indices2));
      }
    }
  } else {
    for (int j = 0; j < ncols; ++j){
      SET_VECTOR_ELT(out, j, cheapr_sset(p_x[j], indices));
    }
  }
  SEXP names = Rf_protect(Rf_duplicate(Rf_getAttrib(x, R_NamesSymbol)));
  ++n_protections;
  Rf_setAttrib(out, R_NamesSymbol, names);
  Rf_protect(out = cpp_list_as_df(out));
  ++n_protections;
  Rf_unprotect(n_protections);
  return out;
}

// SEXP cpp_sset(SEXP x, SEXP indices){
//   if (!Rf_isObject(x) && Rf_isNull(Rf_getAttrib(x, R_NamesSymbol)) && is_alt_int_seq(indices)){
//     SEXP int_seq_data = Rf_protect(Rf_coerceVector(alt_data1(indices), INTSXP));
//     int size = INTEGER(int_seq_data)[0];
//     int from = INTEGER(int_seq_data)[1];
//     int by = INTEGER(int_seq_data)[2];
//     int to = from + (std::max(size - 1, 0) * by);
//     Rf_unprotect(1);
//     return cpp_sset_range(x, from, to, by);
//   }
//   int *pi = INTEGER(indices);
//   int xn = Rf_length(x);
//   int n = Rf_length(indices);
//   int n_protections = 0;
//   int zero_count = 0;
//   int pos_count = 0;
//   int oob_count = 0;
//   int na_count = 0;
//   int k = 0;
//   int out_size;
//   // bool do_parallel = n >= 10000;
//   // int n_cores = do_parallel ? num_cores() : 1;
//
//   // Counting the number of:
//   // Zeroes
//   // Out-of-bounds indices
//   // Positive indices
//   // From this we can also work out the number of negatives
//
// //   if (do_parallel){
// // #pragma omp parallel for simd num_threads(n_cores) reduction(+:zero_count,pos_count,oob_count,na_count)
//   //   for (int j = 0; j < n; ++j){
//   //     zero_count += (pi[j] == 0);
//   //     pos_count += (pi[j] > 0);
//   //     na_count += (pi[j] == NA_INTEGER);
//   //     oob_count += (std::fabs(pi[j]) > xn);
//   //   }
//   // } else {
// // #pragma omp for simd
//     for (int j = 0; j < n; ++j){
//       zero_count += (pi[j] == 0);
//       pos_count += (pi[j] > 0);
//       na_count += pi[j] == NA_INTEGER;
//       oob_count += (std::abs(pi[j]) > xn);
//     }
//   // }
//   int neg_count = n - pos_count - zero_count - na_count;
//   if ( pos_count > 0 && neg_count > 0){
//     Rf_error("Cannot mix positive and negative indices");
//   }
//   bool simple_sset = zero_count == 0 && oob_count == 0 && na_count == 0 && pos_count == n;
//
//   // Convert negative index vector to positive
//
//   if (neg_count > 0){
//     SEXP indices2 = Rf_protect(cpp11::package("cheapr")["neg_indices_to_pos"](indices, xn));
//     // ++n_protections;
//     // int *pi2 = INTEGER(indices2);
//     // pi = pi2;
//     // out_size = Rf_xlength(indices2);
//     // n = out_size;
//     // ++n_protections;
//     Rf_unprotect(n_protections + 1);
//     return cpp_sset_simple(x, indices2);
//   } else {
//     out_size = n - zero_count;
//   }
//   if (simple_sset){
//     Rf_unprotect(n_protections);
//     return cpp_sset_simple(x, indices);
//   }
//   switch ( TYPEOF(x) ){
//   case NILSXP: {
//     return R_NilValue;
//   }
//   case INTSXP: {
//     int *p_x = INTEGER(x);
//     SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
//     ++n_protections;
//     int *p_out = INTEGER(out);
//     for (R_xlen_t i = 0; i < n; ++i){
//       int xi = pi[i];
//       if (xi != 0){
//         p_out[k++] = (xi <= xn && xi != NA_INTEGER) ? p_x[xi - 1] : NA_INTEGER;
//       }
//     }
//     Rf_unprotect(n_protections);
//     return out;
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
// }


// A subset method using c++ vectors

// list cpp_sset(SEXP x, integers i){
//   int xn = Rf_xlength(x);
//   int n = i.size();
//   switch ( TYPEOF(x) ){
//   case INTSXP: {
//     std::vector<int> out;
//     int *p_x = INTEGER(x);
//     out.reserve(n);
//     for (int j = 0; j < n; ++j){
//       if (i[j] > 0 && i[j] <= xn){
//         int val = p_x[i[j] - 1];
//         out.push_back(val);
//       } else {
//         out.push_back(NA_INTEGER);
//       }
//     }
//       return writable::list({
//         "out"_nm = out
//       });
//   }
//   default: {
//     Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
//   }
//   }
// }
