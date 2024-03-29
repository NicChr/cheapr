#include "cheapr_cpp.h"
#include <cpp11.hpp>
#include <Rinternals.h>
// #include <vector>
// using namespace cpp11;

[[cpp11::register]]
SEXP cpp_sset(SEXP x, SEXP indices){
  int *pi = INTEGER(indices);
  int xn = Rf_xlength(x);
  int n = Rf_xlength(indices);
  int n_protections = 0;
  int zero_count = 0;
  int pos_count = 0;
  int oob_count = 0;
  int out_size;
  bool do_parallel = n >= 10000;
  int n_cores = do_parallel ? num_cores() : 1;
  do_parallel = do_parallel && n_cores > 1;

  // Counting the number of:
  // Zeroes
  // Out-of-bounds indices
  // Positive indices
  // From this we can also work out the number of negatives

  if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) reduction(+:zero_count,pos_count,oob_count)
    for (int j = 0; j < n; ++j){
      zero_count += pi[j] == 0;
      pos_count += pi[j] > 0;
      oob_count += std::abs(pi[j]) > xn;
    }
  } else {
#pragma omp for simd
    for (int j = 0; j < n; ++j){
      zero_count += (pi[j] == 0);
      pos_count += (pi[j] > 0);
      oob_count += (std::abs(pi[j]) > xn);
    }
  }
  bool neg_count = n - pos_count - zero_count;
  if ( (pos_count + zero_count) > 0 && neg_count > 0){
    Rf_error("Cannot mix positive and negative indices");
  }
  bool simple_sset = zero_count == 0 && oob_count == 0 && pos_count == n;

  // Convert negative index vector to positive

  if (neg_count > 0){
    SEXP indices2 = Rf_protect(cpp11::package("cheapr")["neg_indices_to_pos"](xn, indices));
    ++n_protections;
    int *pi2 = INTEGER(indices2);
    pi = pi2;
    out_size = Rf_xlength(indices2);
    n = out_size;
    simple_sset = true;
  } else {
    out_size = n - zero_count;
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
    zero_count = 0;
    int *p_out = LOGICAL(out);
    if (simple_sset){
      if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) private(i)
        for (i = 0; i < n; ++i){
          p_out[i] = p_x[pi[i] - 1];
        }
      } else {
#pragma omp for simd
        for (i = 0; i < n; ++i){
          p_out[i] = p_x[pi[i] - 1];
        }
      }
    } else {
      for (i = 0; i < n; ++i){
        if (pi[i] == 0){
          ++zero_count;
        } else {
          p_out[i - zero_count] = (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_LOGICAL;
        }
        // p_out[i - zero_count] = (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_LOGICAL;
        // p_out[i - ( pi[i] == 0 ? zero_count++ : zero_count)] = (pi[i] > 0 && pi[i] <= xn) ? p_x[pi[i] - 1] : NA_INTEGER;
      }
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case INTSXP: {
    int *p_x = INTEGER(x);
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
    ++n_protections;
    zero_count = 0;
    int *p_out = INTEGER(out);
    if (simple_sset){
      if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) private(i)
        for (i = 0; i < n; ++i){
          p_out[i] = p_x[pi[i] - 1];
        }
      } else {
#pragma omp for simd
        for (i = 0; i < n; ++i){
          p_out[i] = p_x[pi[i] - 1];
        }
      }
    } else {
      for (i = 0; i < n; ++i){
        if (pi[i] == 0){
          ++zero_count;
        } else {
          p_out[i - zero_count] = (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_INTEGER;
        }
        // p_out[i - zero_count] = (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_INTEGER;
        // p_out[i - ( pi[i] == 0 ? zero_count++ : zero_count)] = (pi[i] > 0 && pi[i] <= xn) ? p_x[pi[i] - 1] : NA_INTEGER;
      }
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case REALSXP: {
    double *p_x = REAL(x);
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
    ++n_protections;
    zero_count = 0;
    double *p_out = REAL(out);
    if (simple_sset){
      if (do_parallel){
#pragma omp parallel for simd num_threads(n_cores) private(i)
        for (i = 0; i < n; ++i){
          p_out[i] = p_x[pi[i] - 1];
        }
      } else {
#pragma omp for simd
        for (i = 0; i < n; ++i){
          p_out[i] = p_x[pi[i] - 1];
        }
      }
    } else {
      for (i = 0; i < n; ++i){
        if (pi[i] == 0){
          ++zero_count;
        } else {
          p_out[i - zero_count] = (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_REAL;
        }
        // p_out[i - ( pi[i] == 0 ? zero_count++ : zero_count)] = (pi[i] > 0 && pi[i] <= xn) ? p_x[pi[i] - 1] : NA_REAL;
      }
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case STRSXP: {
    SEXP *p_x = STRING_PTR(x);
    SEXP out = Rf_protect(Rf_allocVector(STRSXP, out_size));
    ++n_protections;
    zero_count = 0;
    if (simple_sset){
      for (i = 0; i < n; ++i){
        SET_STRING_ELT(out, i, p_x[pi[i] - 1]);
      }
    } else {
      for (i = 0; i < n; ++i){
        if (pi[i] == 0){
          ++zero_count;
        } else {
          SET_STRING_ELT(out, i - zero_count,
                         (pi[i] <= xn) ? p_x[pi[i] - 1] : NA_STRING);
        }
      }
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case RAWSXP: {
    Rbyte *p_x = RAW(x);
    SEXP out = Rf_protect(Rf_allocVector(RAWSXP, out_size));
    ++n_protections;
    zero_count = 0;
    if (simple_sset){
      for (i = 0; i < n; ++i){
        SET_RAW_ELT(out, i, p_x[pi[i] - 1]);
      }
    } else {
      for (i = 0; i < n; ++i){
        if (pi[i] == 0){
          ++zero_count;
        } else {
          SET_RAW_ELT(out, i - zero_count,
                      (pi[i] <= xn) ? p_x[pi[i] - 1] : 0);
        }
      }
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case VECSXP: {
    const SEXP *p_x = VECTOR_PTR_RO(x);
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, out_size));
    ++n_protections;
    zero_count = 0;
    if (simple_sset){
      for (i = 0; i < n; ++i){
        SET_VECTOR_ELT(out, i, p_x[pi[i] - 1]);
      }
    } else {
      for (i = 0; i < n; ++i){
        if (pi[i] == 0){
          ++zero_count;
        } else {
          SET_VECTOR_ELT(out, i - zero_count,
                         (pi[i] <= xn) ? p_x[pi[i] - 1] : R_NilValue);
        }
      }
    }
    Rf_unprotect(n_protections);
    return out;
  }
  default: {
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}


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
