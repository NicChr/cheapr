#ifndef CHEAPR_R_DATES_H
#define CHEAPR_R_DATES_H

#include <cheapr/internal/r_vector.h>
#include <cheapr/internal/r_attrs.h>

namespace cheapr {

// R Date vector
struct r_dates : public r_vec<r_int> {

  // Constructors
  r_dates() : r_vec<r_int>() {}

  explicit r_dates(SEXP x) : r_vec<r_int>(x) {
    if (!is_null() && !attr::inherits1(sexp, "Date")){
      Rf_error("`SEXP` must be a Date");
    }
  }

  explicit r_dates(r_size_t n) : r_vec<r_int>(n) { 
    attr::set_old_class(this->sexp, r_vec<r_str>(1, r_str("Date")));
  }
};

}

#endif
