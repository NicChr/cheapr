#include "cheapr_cpp.h"

[[cpp11::register]]
SEXP cpp_lag(SEXP x, int k, SEXP fill, bool set, bool recursive) {
  R_xlen_t size = Rf_xlength(x);
  int fill_size = Rf_length(fill);
  int n_protections = 0;
  if (fill_size > 1){
    Rf_error("fill size must be NULL or length 1");
  }
  R_xlen_t klag = k;
  switch(TYPEOF(x)){
  case NILSXP: {
    return R_NilValue;
  }
  case LGLSXP:
  case INTSXP: {
    k = k >= 0 ? std::min(size, klag) : std::max(-size, klag);
    int fill_value = NA_INTEGER;
    if (fill_size >= 1){
      fill_value = Rf_asInteger(fill);
    }
    SEXP out = Rf_protect(set ? x : Rf_duplicate(x));
    ++n_protections;
    int *p_out = INTEGER(out);
    if (set){
      R_xlen_t tempi;
      int tempv;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(INTSXP, std::abs(k)));
        ++n_protections;
        int* __restrict__ p_lag = INTEGER(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            p_lag[i] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            p_lag[size - i - 1] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
        }

      }
    } else {
      int *p_x = INTEGER(x);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          p_out[i] = i >= k ? p_x[i - k] : fill_value;
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          p_out[i] = (i - size) < k ? p_x[i - k] : fill_value;
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case REALSXP: {
    k = k >= 0 ? std::min(size, klag) : std::max(-size, klag);
    double fill_value = NA_REAL;
    if (fill_size >= 1){
      fill_value = Rf_asReal(fill);
    }
    SEXP out = Rf_protect(set ? x : Rf_duplicate(x));
    ++n_protections;
    double *p_out = REAL(out);
    if (set){
      R_xlen_t tempi;
      double tempv;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(REALSXP, std::abs(k)));
        ++n_protections;
        double* __restrict__ p_lag = REAL(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            p_lag[i] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            p_lag[size - i - 1] = p_out[i];
            p_out[i] = fill_value;
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            tempv = p_lag[tempi];
            p_lag[tempi] = p_out[i];
            p_out[i] = tempv;
          }
        }

      }
    } else {
      double *p_x = REAL(x);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          p_out[i] = i >= k ? p_x[i - k] : fill_value;
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          p_out[i] = (i - size) < k ? p_x[i - k] : fill_value;
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case CPLXSXP: {
    k = k >= 0 ? std::min(size, klag) : std::max(-size, klag);
    SEXP fill_sexp = Rf_protect(Rf_allocVector(CPLXSXP, 1));
    ++n_protections;
    Rcomplex *p_fill = COMPLEX(fill_sexp);
    p_fill[0].i = NA_REAL;
    p_fill[0].r = NA_REAL;
    Rcomplex fill_value = fill_size >= 1 ? Rf_asComplex(fill) : COMPLEX(fill_sexp)[0];
    SEXP out = Rf_protect(set ? x : Rf_duplicate(x));
    ++n_protections;
    Rcomplex *p_out = COMPLEX(out);
    if (set){
      R_xlen_t tempi;
      Rcomplex tempv;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(CPLXSXP, std::abs(k)));
        ++n_protections;
        Rcomplex* __restrict__ p_lag = COMPLEX(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            p_lag[i].i = p_out[i].i;
            p_lag[i].r = p_out[i].r;
            SET_COMPLEX_ELT(out, i, fill_value);
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            tempv = p_lag[tempi];
            SET_COMPLEX_ELT(lag_temp, tempi, p_out[i]);
            SET_COMPLEX_ELT(out, i, tempv);
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_COMPLEX_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_COMPLEX_ELT(out, i, fill_value);
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            tempv = p_lag[tempi];
            SET_COMPLEX_ELT(lag_temp, tempi, p_out[i]);
            SET_COMPLEX_ELT(out, i, tempv);
          }
        }

      }
    } else {
      Rcomplex *p_x = COMPLEX(x);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          SET_COMPLEX_ELT(out, i, i >= k ? p_x[i - k] : fill_value);
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          SET_COMPLEX_ELT(out, i, (i - size) < k ? p_x[i - k] : fill_value);
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case STRSXP: {
    k = k >= 0 ? std::min(size, klag) : std::max(-size, klag);
    SEXP fill_char = Rf_protect(fill_size >= 1 ? Rf_asChar(fill) : NA_STRING);
    ++n_protections;
    SEXP out = Rf_protect(set ? x : Rf_duplicate(x));
    ++n_protections;
    SEXP *p_out = STRING_PTR(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(STRSXP, std::abs(k)));
        ++n_protections;
        SEXP tempv = Rf_protect(Rf_allocVector(STRSXP, 1));
        ++n_protections;
        SEXP* __restrict__ p_lag = STRING_PTR(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            SET_STRING_ELT(lag_temp, i, p_out[i]);
            SET_STRING_ELT(out, i, fill_char);
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            SET_STRING_ELT(tempv, 0, p_lag[tempi]);
            SET_STRING_ELT(lag_temp, tempi, p_out[i]);
            SET_STRING_ELT(out, i, STRING_ELT(tempv, 0));
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_STRING_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_STRING_ELT(out, i, fill_char);
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            SET_STRING_ELT(tempv, 0, p_lag[tempi]);
            SET_STRING_ELT(lag_temp, tempi, p_out[i]);
            SET_STRING_ELT(out, i, STRING_ELT(tempv, 0));
          }
        }
      }
    } else {
      SEXP *p_x = STRING_PTR(x);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          SET_STRING_ELT(out, i, i >= k ? p_x[i - k] : fill_char);
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          SET_STRING_ELT(out, i, (i - size) < k ? p_x[i - k] : fill_char);
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case RAWSXP: {
    k = k >= 0 ? std::min(size, klag) : std::max(-size, klag);
    SEXP raw_sexp = Rf_protect(Rf_coerceVector(fill, RAWSXP));
    Rbyte fill_raw = fill_size == 0 ? RAW(Rf_ScalarRaw(0))[0] : RAW(raw_sexp)[0];
    ++n_protections;
    SEXP out = Rf_protect(set ? x : Rf_duplicate(x));
    ++n_protections;
    Rbyte *p_out = RAW(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(RAWSXP, std::abs(k)));
        ++n_protections;
        Rbyte *p_lag = RAW(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            SET_RAW_ELT(lag_temp, i, p_out[i]);
            SET_RAW_ELT(out, i, fill_raw);
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            Rbyte tempv = p_lag[tempi];
            SET_RAW_ELT(lag_temp, tempi, p_out[i]);
            SET_RAW_ELT(out, i, tempv);
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_RAW_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_RAW_ELT(out, i, fill_raw);
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            Rbyte tempv = p_lag[tempi];
            SET_RAW_ELT(lag_temp, tempi, p_out[i]);
            SET_RAW_ELT(out, i, tempv);
          }
        }
      }
    } else {
      const Rbyte *p_x = RAW_RO(x);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          SET_RAW_ELT(out, i, i >= k ? p_x[i - k] : fill_raw);
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          SET_RAW_ELT(out, i, (i - size) < k ? p_x[i - k] : fill_raw);
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  case VECSXP: {
    if (recursive){
    const SEXP *p_x = VECTOR_PTR_RO(x);
    SEXP out = Rf_protect(Rf_allocVector(VECSXP, size));
    ++n_protections;
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    for (R_xlen_t i = 0; i < size; ++i){
      SET_VECTOR_ELT(out, i, cpp_lag(p_x[i], k, fill, set, true));
    }
    Rf_unprotect(n_protections);
    return out;
  } else {
    k = k >= 0 ? std::min(size, klag) : std::max(-size, klag);
    SEXP fill_value = Rf_protect(Rf_coerceVector(fill_size >= 1 ? fill : R_NilValue, VECSXP));
    ++n_protections;
    SEXP out = Rf_protect(set ? x : Rf_allocVector(VECSXP, size));
    ++n_protections;
    const SEXP *p_out = VECTOR_PTR_RO(out);
    if (set){
      R_xlen_t tempi;
      // If k = 0 then no lag occurs
      if (std::abs(k) >= 1){
        SEXP lag_temp = Rf_protect(Rf_allocVector(VECSXP, std::abs(k)));
        ++n_protections;
        SEXP tempv = Rf_protect(Rf_allocVector(VECSXP, 1));
        ++n_protections;
        SEXP* __restrict__ p_lag = VECTOR_PTR(lag_temp);
        // Positive lags
        if (k >= 0){
          for (R_xlen_t i = 0; i < k; ++i) {
            SET_VECTOR_ELT(lag_temp, i, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(fill_value, 0));
          }
          for (R_xlen_t i = k; i < size; ++i) {
            tempi = ((i - k) % k);
            SET_VECTOR_ELT(tempv, 0, p_lag[tempi]);
            SET_VECTOR_ELT(lag_temp, tempi, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(tempv, 0));
          }
          // Negative lags
        } else {
          for (R_xlen_t i = size - 1; i >= size + k; --i) {
            SET_VECTOR_ELT(lag_temp, size - i - 1, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(fill_value, 0));
          }
          for (R_xlen_t i = size + k - 1; i >= 0; --i) {
            tempi = ( (size - (i - k) - 1) % k);
            SET_VECTOR_ELT(tempv, 0, p_lag[tempi]);
            SET_VECTOR_ELT(lag_temp, tempi, p_out[i]);
            SET_VECTOR_ELT(out, i, VECTOR_ELT(tempv, 0));
          }
        }
      }
    } else {
      const SEXP *p_x = VECTOR_PTR_RO(x);
      if (k >= 0){
        for (R_xlen_t i = 0; i < size; ++i) {
          SET_VECTOR_ELT(out, i, i >= k ? p_x[i - k] : VECTOR_ELT(fill_value, 0));
        }
      } else {
        for (R_xlen_t i = size - 1; i >= 0; --i) {
          SET_VECTOR_ELT(out, i, (i - size) < k ? p_x[i - k] : VECTOR_ELT(fill_value, 0));
        }
      }
    }
    if (!Rf_isNull(Rf_getAttrib(x, R_NamesSymbol))){
      SEXP old_names = Rf_protect(Rf_getAttrib(x, R_NamesSymbol));
      ++n_protections;
      SEXP new_names = Rf_protect(cpp_lag(old_names, k, R_NilValue, set, recursive));
      ++n_protections;
      Rf_setAttrib(out, R_NamesSymbol, new_names);
    }
    Rf_unprotect(n_protections);
    return out;
  }
  }
  default: {
    Rf_unprotect(n_protections);
    Rf_error("%s cannot handle an object of type %s", __func__, Rf_type2char(TYPEOF(x)));
  }
  }
}
