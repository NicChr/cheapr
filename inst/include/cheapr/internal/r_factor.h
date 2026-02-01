#ifndef CHEAPR_R_FACTOR_H
#define CHEAPR_R_FACTOR_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_match.h>
#include <cheapr/internal/r_unique.h>
#include <cheapr/internal/r_attrs.h>
#include <vector>

namespace cheapr {

struct r_factors : public r_vec<r_int> { 

  public: 

  r_vec<r_str> levels() const {
    return r_vec<r_str>(attr::get_attr(*this, r_sym("levels")));
  }

  void set_levels(const r_vec<r_str>& levels) {
    attr::set_attr(*this, r_sym("levels"), levels);
  }

  private:
  void init_factor_attrs(const r_vec<r_str>& levels) {
      // Set class
      attr::set_attr(*this, r_sym("class"), r_vec<r_str>(1, "factor"));
      // Set levels
      set_levels(levels);
  }

  public: 

  // Constructors
  r_factors() : r_vec<r_int>() {
    init_factor_attrs(r_vec<r_str>());
  }

  explicit r_factors(SEXP x) : r_vec<r_int>(x) {
    if (!is_null() && !attr::inherits1(*this, "factor")){
      abort("`SEXP` must be a factor");
    }
  }

  template <RScalar T>
  explicit r_factors(const r_vec<T>& x, const r_vec<T>& levels){
    auto fct = match(x, levels);
    r_vec<r_int>::operator=(std::move(fct));
    r_vec<r_str> str_levels;
    if constexpr (is<T, r_str>){
      str_levels = levels;
    } else {
      r_size_t n = levels.length();
      str_levels = r_vec<r_str>(n);
      for (r_size_t i = 0; i < n; ++i){
        str_levels.set(i, internal::as_r<r_str>(levels.get(i)));
      }
    }
    init_factor_attrs(str_levels);
  }

  template <RScalar T>
  explicit r_factors(const r_vec<T>& x) : r_factors(x, unique(x)) {}

};

}

#endif
