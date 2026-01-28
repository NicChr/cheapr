#ifndef CHEAPR_R_DATES_H
#define CHEAPR_R_DATES_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_attrs.h>

namespace cheapr {

// R Date vector
struct r_dates : public r_vec<r_int> {

private:
  void init_date_attrs() {
      attr::set_old_class(*this, r_vec<r_str>(1, r_str("Date")));
  }

public: 
  // Constructors
  r_dates() : r_vec<r_int>() {
    init_date_attrs();
  }

  explicit r_dates(SEXP x) : r_vec<r_int>(x) {
    if (!is_null() && !attr::inherits1(this->sexp, "Date")){
      abort("`SEXP` must be a Date");
    }
  }

  explicit r_dates(r_size_t n) : r_vec<r_int>(n) { 
    init_date_attrs();
  }
  template<typename U>
  explicit r_dates(r_size_t n, U default_value) : r_vec<r_int>(n, default_value) { 
    init_date_attrs();
  }
};

}

#endif
