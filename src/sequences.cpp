#include "cheapr.h"

double r_sum(SEXP x, bool na_rm = false){
  cpp11::function base_sum = cpp11::package("base")["sum"];
  return Rf_asReal(base_sum(x, cpp11::named_arg("na.rm") = na_rm));
}

double r_min(SEXP x){
  cpp11::function base_min = cpp11::package("base")["min"];
  double out = R_PosInf;
  if (Rf_xlength(x) > 0){
    out = Rf_asReal(base_min(x));
  }
  return out;
}

// My version of base::sequence()

[[cpp11::register]]
SEXP cpp_int_sequence(SEXP size, SEXP from, SEXP by) {
  int size_n = Rf_length(size);
  int from_n = Rf_length(from);
  int by_n = Rf_length(by);
  if (size_n > 0 && (from_n <= 0 || by_n <= 0)){
    Rf_error("from and by must both have length > 0");
  }
  double out_size = r_sum(size, false);
  double min_size = r_min(size);
  if (!(out_size == out_size)){
    Rf_error("size must not contain NA values");
  }
  if (min_size < 0){
    Rf_error("size must be a vector of non-negative integers");
  }
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, out_size));
  int *p_out = INTEGER(out);
  R_xlen_t index = 0, interrupt_counter = 0;
  int fj;
  int bj;
  int start;
  int increment;
  int seq_size;
  double seq_end;
  if (size_n > 0){
    int *p_size = INTEGER(size);
    int *p_from = INTEGER(from);
    int *p_by = INTEGER(by);
    for (int j = 0; j < size_n; ++j){
      seq_size = p_size[j];
      fj = j % from_n;
      bj = j % by_n;
      start = p_from[fj];
      increment = p_by[bj];
      // Throw error if integer overflow
      seq_end = ( (std::fmax(seq_size - 1, 0.0)) * increment ) + start;
      if (std::fabs(seq_end) > integer_max_){
        Rf_unprotect(1);
        Rf_error("Integer overflow value of %g in sequence %d", seq_end, j + 1);
      }
      if (start == NA_INTEGER){
        Rf_unprotect(1);
        Rf_error("from contains NA values");
      }
      if (increment == NA_INTEGER){
        Rf_unprotect(1);
        Rf_error("by contains NA values");
      }
      for (int i = 0; i < seq_size; ++i){
        if (interrupt_counter == 100000000){
          R_CheckUserInterrupt();
          interrupt_counter = 0;
        }
        p_out[index] = start;
        start += increment;
        ++index;
        ++interrupt_counter;
      }
    }
  }
  Rf_unprotect(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_dbl_sequence(SEXP size, SEXP from, SEXP by) {
  int size_n = Rf_length(size);
  int from_n = Rf_length(from);
  int by_n = Rf_length(by);
  if (size_n > 0 && (from_n <= 0 || by_n <= 0)){
    Rf_error("from and by must both have length > 0");
  }
  // To recycle we would need to do sum * remainder of the sum over n
  double out_size = r_sum(size, false);
  double min_size = r_min(size);
  if (!(out_size == out_size)){
    Rf_error("size must not contain NA values");
  }
  if (min_size < 0){
    Rf_error("size must be a vector of non-negative integers");
  }
  SEXP out = Rf_protect(Rf_allocVector(REALSXP, out_size));
  double *p_out = REAL(out);
  R_xlen_t index = 0, interrupt_counter = 0;
  int fj;
  int bj;
  int seq_size;
  double start;
  double increment;
  if (size_n > 0){
    int *p_size = INTEGER(size);
    double *p_from = REAL(from);
    double *p_by = REAL(by);
    for (int j = 0; j < size_n; ++j){

      seq_size = p_size[j];
      fj = j % from_n;
      bj = j % by_n;
      start = p_from[fj];
      increment = p_by[bj];
      if (!(start == start)){
        Rf_unprotect(1);
        Rf_error("from contains NA values");
      }
      if (!(increment == increment)){
        Rf_unprotect(1);
        Rf_error("by contains NA values");
      }
      for (int i = 0; i < seq_size; ++i){
        if (interrupt_counter == 100000000){
          R_CheckUserInterrupt();
          interrupt_counter = 0;
        }
        p_out[index] = ( start + (i * increment) );
        ++index;
        ++interrupt_counter;
      }
    }
  }
  Rf_unprotect(1);
  return out;
}

[[cpp11::register]]
SEXP cpp_sequence(SEXP size, SEXP from, SEXP by) {
  int size_n = Rf_length(size);
  int from_n = Rf_length(from);
  int by_n = Rf_length(by);
  switch (TYPEOF(from)){
  case INTSXP: {
    switch (TYPEOF(by)){
  case INTSXP: {
    int n = std::max(std::max(size_n, from_n), by_n);
    double seq_end;
    bool out_is_integer = true;
    int *p_size = INTEGER(size);
    int *p_from = INTEGER(from);
    int *p_by = INTEGER(by);

    // Checking that the sequence values are integers
    // Only do the loop if vectors are not zero-length
    if (size_n > 0 && from_n > 0 && by_n > 0){
      for (int i = 0; i < n; ++i){
        seq_end = (std::fmax(p_size[i % size_n] - 1.0, 0.0) * p_by[i % by_n]) + p_from[i % from_n];
        if (seq_end > integer_max_){
          out_is_integer = false;
          break;
        }
      }
    }
    // If all sequence values are < 2^31 then we can safely use cpp_int_sequence
    if (out_is_integer){
      return cpp_int_sequence(size, from, by);
    } else {
      Rf_protect(from = Rf_coerceVector(from, REALSXP));
      Rf_protect(by = Rf_coerceVector(by, REALSXP));
      SEXP out = Rf_protect(cpp_dbl_sequence(size, from, by));
      Rf_unprotect(3);
      return out;
    }
  }
  case REALSXP: {
    Rf_protect(from = Rf_coerceVector(from, REALSXP));
    SEXP out = Rf_protect(cpp_dbl_sequence(size, from, by));
    Rf_unprotect(2);
    return out;
  }
  default: {
    Rf_error("by must have type integer or double in %s", __func__);
  }
  }
    break;
  }
  case REALSXP: {
    switch (TYPEOF(by)){
  case INTSXP: {
    Rf_protect(by = Rf_coerceVector(by, REALSXP));
    SEXP out = Rf_protect(cpp_dbl_sequence(size, from, by));
    Rf_unprotect(2);
    return out;
  }
  case REALSXP: {
    return cpp_dbl_sequence(size, from, by);
  }
  default: {
    Rf_error("by must have type integer or double in %s", __func__);
  }
  }
  }
  default: {
    Rf_error("from must have type integer or double in %s", __func__);
  }
  }
}

[[cpp11::register]]
SEXP cpp_window_sequence(SEXP size,
                         double k,
                         bool partial = true,
                         bool ascending = true) {
  int size_n = Rf_length(size);
  SEXP size_sexp = Rf_protect(Rf_coerceVector(size, INTSXP));
  if (r_min(size_sexp) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  k = std::fmax(k, 0);
  R_xlen_t N = r_sum(size_sexp, false);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, N));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size_sexp);
  R_xlen_t index = 0;
  if (ascending){
    // right aligned window sequences
    if (partial){
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          if (i < k){
            p_out[index] = i + 1;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    } else {
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          if (i < (k - 1)){
            p_out[index] = NA_INTEGER;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    }
  } else {
    // left aligned window sequences
    int idiff;
    if (partial){
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          idiff = p_size[j] - i - 1;
          if (idiff < k){
            p_out[index] = idiff + 1;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    } else {
      for (int j = 0; j < size_n; ++j){
        for (int i = 0; i < p_size[j]; ++i){
          idiff = p_size[j] - i - 1;
          if (idiff < (k - 1)){
            p_out[index] = NA_INTEGER;
          } else {
            p_out[index] = k;
          }
          ++index;
        }
      }
    }
  }
  Rf_unprotect(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_lag_sequence(SEXP size, double k, bool partial = false) {
  Rf_protect(size = Rf_coerceVector(size, INTSXP));
  if (r_min(size) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  int size_n = Rf_length(size);
  k = std::fmax(k, 0);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, r_sum(size, false)));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size);
  R_xlen_t index = 0;
  if (partial){
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        if (i < k){
          p_out[index] = i;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  } else {
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        if (i < k){
          p_out[index] = NA_INTEGER;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  }
  Rf_unprotect(2);
  return out;
}
[[cpp11::register]]
SEXP cpp_lead_sequence(SEXP size, double k, bool partial = false) {
  Rf_protect(size = Rf_coerceVector(size, INTSXP));
  if (r_min(size) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  int size_n = Rf_length(size);
  k = std::fmax(k, 0);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, r_sum(size, false)));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size);
  R_xlen_t index = 0;
  int idiff;
  if (partial){
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        idiff = p_size[j] - i - 1;
        if (idiff < k){
          p_out[index] = idiff;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  } else {
    for (int j = 0; j < size_n; ++j){
      for (int i = 0; i < p_size[j]; ++i){
        idiff = p_size[j] - i - 1;
        if (idiff < k){
          p_out[index] = NA_INTEGER;
        } else {
          p_out[index] = k;
        }
        ++index;
      }
    }
  }
  Rf_unprotect(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_sequence_id(SEXP size){
  int size_n = Rf_length(size);
  SEXP size_sexp = Rf_protect(Rf_coerceVector(size, INTSXP));
  if (r_min(size_sexp) < 0){
    Rf_unprotect(1);
    Rf_error("size must be a vector of non-negative integers");
  }
  R_xlen_t N = r_sum(size_sexp, false);
  SEXP out = Rf_protect(Rf_allocVector(INTSXP, N));
  int *p_out = INTEGER(out);
  int *p_size = INTEGER(size_sexp);
  R_xlen_t k = 0;
  int seq_size;
  for (int i = 0; i < size_n; ++i){
    seq_size = p_size[i];
    for (int j = 0; j < seq_size; ++j){
      p_out[k++] = i + 1;
    }
  }
  Rf_unprotect(2);
  return out;
}

SEXP cpp_seq_len(R_xlen_t n){
  if (n > integer_max_){
    SEXP out = Rf_protect(Rf_allocVector(REALSXP, n));
    double *p_out = REAL(out);
    OMP_FOR_SIMD
    for (R_xlen_t i = 0; i < n; ++i) p_out[i] = 1.0 + i;
    Rf_unprotect(1);
    return out;
  } else {
    SEXP out = Rf_protect(Rf_allocVector(INTSXP, n));
    int *p_out = INTEGER(out);
    OMP_FOR_SIMD
    for (int i = 0; i < n; ++i) p_out[i] = i + 1;
    Rf_unprotect(1);
    return out;
  }
}

bool is_whole_number(double x, double tolerance){
  return (std::fabs(x - std::round(x)) < tolerance);
}

double log10_scale(double x){
  return std::floor(std::log10(std::fabs(x == 0.0 ? 1.0 : x)));
}
double nearest_floor(double x, double n){
  return std::floor(x / n) * n;
}
double nearest_ceiling(double x, double n){
  return std::ceil(x / n) * n;
}

double pretty_ceiling(double x){
  return nearest_ceiling(x, (std::pow(10.0, log10_scale(x))));
}

// Not sure yet if approximate division is needed

// double aprx_floor(double x, double tol){
//   return is_whole_number(x, tol) ? std::round(x) : std::floor(x);
// }
// double aprx_ceil(double x, double tol){
//   return is_whole_number(x, tol) ? std::round(x) : std::ceil(x);
// }
// double aprx_trunc(double x, double tol){
//   return is_whole_number(x, tol) ? std::round(x) : std::trunc(x);
// }

// truncated division with a tolerance
// aka integer division towards zero
// double int_div_towards_zero(double x, double y, double tol){
//   return aprx_trunc(x / y, tol);
// }

// approximate ceiling division
// double ceil_div(double x, double y, double tol){
//   return aprx_ceil(x / y, tol);
// }

// integer division away from zero
// double int_div_away_from_zero(double x, double y, double tolerance){
//   double out = x / y;
//   if (is_whole_number(out, tolerance)){
//     return std::round(out);
//   } else {
//     double sign = (out > 0) - (out < 0);
//     return std::ceil(std::fabs(out)) * sign;
//   }
// }

bool can_be_int(double x, double tolerance){
  return is_whole_number(x, tolerance) && std::fabs(x) <= integer_max_;
}

double seq_end(double size, double from, double by){
  return from + (std::fmax(size - 1.0, 0.0) * by);
}
double seq_width(double size, double from, double to){
  double del = (to - from);
  double out = del / std::fmax(size - 1.0, 0.0);
  return del == 0.0 ? 0.0 : out;
}

bool is_infinite(double x){
  return (x == R_NegInf) || (x == R_PosInf);
}

[[cpp11::register]]
SEXP cpp_fixed_width_breaks(double start, double end, double n,
                            bool pretty, bool expand_min, bool expand_max){
  if (cheapr_is_na_dbl(n)){
    Rf_error("n must not be `NA`");
  }
  if (n < 1){
    Rf_error("n must be >= 1");
  }
  if (n >= R_PosInf){
    Rf_error("n must be finite");
  }
  if (cheapr_is_na_dbl(start) || cheapr_is_na_dbl(end) ||
      is_infinite(start) || is_infinite(end)){
    return Rf_ScalarReal(NA_REAL);
  }
  // Switch them if needed
  if (start > end){
    double temp = end;
    end = start;
    start = temp;
  }
  double rng_width = end - start;
  bool spans_zero = (start < 0.0 && end > 0.0) || (start > 0.0 && end < 0.0);
  bool zero_range = rng_width == 0.0;

  // e^2/3 which is basically saying 1/3 of the representable digits
  // aren't compared
  // e^1/2 says 1/2 of them aren't to be compared so e^2/3 is a bit more strict
  // as the tolerance is smaller
  double tol = std::pow(std::numeric_limits<double>::epsilon(), 2.0/3.0);

  double n_breaks,
  bin_width, adj_width,
  adj_start, adj_end, adj_rng_width,
  scale, scale_diff,
  zero_adjustment, n_rm, n_add;
  double scale_adj = 1;
  bool lands_on_zero;

  bool scale_up = false;

  if (zero_range){
    if (start == 0.0){
      rng_width = 1.0;
    } else {
      rng_width = std::fabs(start);
    }
    adj_start = start - (rng_width / 1000.0);
    adj_end = end + (rng_width / 1000.0);
    SEXP size = Rf_protect(Rf_ScalarInteger(n + 1.0));
    SEXP from = Rf_protect(Rf_ScalarReal(adj_start));
    SEXP by = Rf_protect(Rf_ScalarReal(seq_width(n + 1.0, adj_start, adj_end)));
    SEXP out = Rf_protect(cpp_dbl_sequence(size, from, by));
    Rf_unprotect(4);
    return out;
  }

  if (!pretty){
    bin_width = rng_width / n;
    int out_size = n + 1.0;
    if (expand_min){
      start -= bin_width;
      ++out_size;
    }
    if (expand_max){
      ++out_size;
    }
    SEXP size_sexp = Rf_protect(Rf_ScalarInteger(out_size));
    SEXP start_sexp = Rf_protect(Rf_ScalarReal(start));
    SEXP width_sexp = Rf_protect(Rf_ScalarReal(bin_width));
    SEXP out = Rf_protect(cpp_dbl_sequence(size_sexp, start_sexp, width_sexp));
    Rf_unprotect(4);
    return out;

  } else {

    // Making breaks prettier

    // This is the number of orders of magnitude the data spans
    scale_diff = log10_scale(rng_width);

    // If large range & relatively small starting value
    // floor start to the nearest difference in orders of magnitude

    // If range spans across zero and start val is small
    if (scale_diff >= 1.0 && spans_zero && std::fabs(start) < 1.0){
      adj_start = nearest_floor(start, std::pow(10.0, log10_scale(end)));
    } else {
      adj_start = nearest_floor(start, std::pow(10.0, scale_diff));
    }

    // Calculate bin-width (guaranteed to span end-points inclusively)
    adj_rng_width = end - adj_start;
    bin_width = adj_rng_width / n;

    // Make width look nicer
    if (bin_width > 2.0 && bin_width < 5.0){
      adj_width = 5.0;
    } else if (bin_width > 5.0 && bin_width < 10.0){
      adj_width = 10.0;
    } else if (bin_width < 1.0){
      adj_width = nearest_ceiling(bin_width, std::pow(10.0, log10_scale(bin_width)) * 5.0);
    } else {
      adj_width = pretty_ceiling(bin_width);
    }

    n_breaks = n;

    // We scale up to work with whole numbers
    scale_up = adj_width < 1.0;
    if (scale_up){

      // The below only works because
      // a) width won't contain decimals for widths >= 1
      // b) It is a decimal of the form 0.0..n for >= 0 zeros
      // so -log10(width) tells us the number of decimals here

      scale = std::fabs(log10_scale(adj_width));
      scale_adj = std::pow(10.0, scale);
      adj_width = round_nearest_even(adj_width * scale_adj);

      // Not sure if adj_start needs rounding here
      adj_start *= scale_adj;

      if (is_whole_number(adj_start, tol)){
        adj_start = round_nearest_even(adj_start);
      }
      start *= scale_adj;
      end *= scale_adj;
    }

    // If breaks span zero make sure they actually land on zero
    if (spans_zero){
      zero_adjustment = adj_start + (adj_width * std::ceil(std::fabs(adj_start) / adj_width));
      lands_on_zero = std::fabs(zero_adjustment) < tol;
      if (!lands_on_zero){
        adj_start -= zero_adjustment;
        n_breaks += std::trunc(zero_adjustment / adj_width);
        // n_breaks += int_div_towards_zero(zero_adjustment, adj_width, tol);
      }
    }

    // Final break?
    adj_end = seq_end(n_breaks, adj_start, adj_width);

    // If too many breaks, reduce
    if (adj_end > end){
      n_rm = std::ceil(adj_end - end / adj_width);
      n_breaks -= n_rm;
      adj_end -= (adj_width * n_rm);
    }

    // adj_end might also have too few breaks
    if (adj_end < end){
      // n_add = int_div_towards_zero(end - adj_end, adj_width, tol);
      n_add = std::trunc((end - adj_end) / adj_width);
      n_breaks += n_add;
      adj_end += (adj_width * n_add);
    }

    if (adj_start < start){
      // n_rm = int_div_towards_zero(start - adj_start, adj_width, tol);
      n_rm = std::trunc((start - adj_start) / adj_width);
      n_breaks -= n_rm;
      adj_start += (adj_width * n_rm);
    }

    // At this point, adj_end will be in range (end - adj_width, end]
    // If we want last break to extend beyond data, add adj_width to it

    if (expand_max && adj_end <= end){
      n_breaks += 1.0;
      adj_end += adj_width;
    }

    if (expand_min && adj_start >= start){
      n_breaks += 1.0;
      adj_start -= adj_width;
    }

    // Work with integers where possible

    // If we have to scale up it means width was a decimal and
    // hence result is not an integer

    SEXP out, seq_size, seq_from, seq_width;
    out = R_NilValue, seq_size = R_NilValue;
    seq_from = R_NilValue, seq_width = R_NilValue;

    if (!scale_up && can_be_int(adj_width, tol) && can_be_int(adj_start, tol) &&
        std::fabs(adj_end) <= integer_max_){
      adj_width = round_nearest_even(adj_width);
      adj_start = round_nearest_even(adj_start);

      seq_size = Rf_protect(Rf_ScalarInteger(n_breaks));
      seq_from = Rf_protect(Rf_ScalarInteger(adj_start));
      seq_width = Rf_protect(Rf_ScalarInteger(adj_width));

      out = Rf_protect(cpp_int_sequence(seq_size, seq_from, seq_width));
    }

    if (scale_up){
      seq_size = Rf_protect(Rf_ScalarInteger(n_breaks));
      seq_from = Rf_protect(Rf_ScalarReal(adj_start));
      seq_width = Rf_protect(Rf_ScalarReal(adj_width));

      out = Rf_protect(cpp_dbl_sequence(seq_size, seq_from, seq_width));
      int seq_n = n_breaks;
      double *p_out = REAL(out);
      OMP_FOR_SIMD
      for (int i = 0; i < seq_n; ++i) p_out[i] = p_out[i] / scale_adj;
    }

    Rf_unprotect(4);
    return out;
  }
}
