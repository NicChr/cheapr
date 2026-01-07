#ifndef CHEAPR_R_FACTOR_H
#define CHEAPR_R_FACTOR_H

#include <cheapr/internal/r_vector.h>
#include <cheapr/internal/r_coerce.h>
#include <ankerl/unordered_dense.h> // For unique strings and matching


namespace cheapr {

namespace internal {

r_vec<r_str> unique_strings(r_vec<r_str> x) {
  R_xlen_t n = x.length();
  
  // Hash set for O(n) de-duplication
  ankerl::unordered_dense::set<SEXP> seen;
  seen.reserve(n);
  
  auto unique_vec = SHIELD(r_vec<r_str>(n));
  R_xlen_t n_unique = 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    r_str str = x.get(i);
    if (seen.insert(str).second) {
        unique_vec.set(n_unique++, str);
    }
  }

  if (n_unique == n){
    YIELD(1);
    return unique_vec;
  }

  auto out = SHIELD(r_vec<r_str>(n_unique));

  for (R_xlen_t i = 0; i < n_unique; ++i){
    out.set(i, unique_vec.get(i));
  }

  YIELD(2);
  return out;
}


// Match needle strings to first occurrence in haystack
r_vec<r_int> string_match(r_vec<r_str> needles, r_vec<r_str> haystack) { 
  
  R_xlen_t n_needles = needles.length();
  R_xlen_t n_haystack = haystack.length();

  if (n_haystack > r_limits::r_int_max){
    Rf_error("Cannot match to a long vector, please use a short character vector");
  }
  
  // Build hash table
  ankerl::unordered_dense::map<SEXP, int> lookup;
  lookup.reserve(n_haystack);
  
  for (R_xlen_t i = 0; i < n_haystack; ++i) {
    r_str str = haystack.get(i);
    lookup.try_emplace(str, static_cast<int>(i) + 1);
  }
  
  // Match needles
  auto out = SHIELD(r_vec<r_int>(n_needles));
  for (R_xlen_t i = 0; i < n_needles; ++i) {
    r_str needle = needles.get(i);
    auto it = lookup.find(needle);
    out.set(i, it != lookup.end() ? r_int(it->second) : na::integer);
  }
  
  YIELD(1);
  return out;
}

}


struct r_factors : public r_vec<r_int> {
  
  // Constructors
  r_factors() : r_vec<r_int>() {}
  
  explicit r_factors(SEXP x) : r_vec<r_int>(x) {
    if (!(is_null() || attr::inherits1(x, "factor"))){ 
      Rf_error("`SEXP` must be a factor");
    }
  }
  
  explicit r_factors(r_vec<r_str> x, r_vec<r_str> levels){

    auto fct = SHIELD(internal::string_match(x, levels));
    auto cls = SHIELD(internal::new_scalar_vector(internal::as_r<r_str>("factor")));
    // Set class
    attr::set_old_class(fct, cls);
    // Set levels
    attr::set_attr(fct, internal::as_r<r_sym>("levels"), levels);
    this->value = fct.value;
    YIELD(2);
  }

  explicit r_factors(r_vec<r_str> x){
    auto levels = SHIELD(internal::unique_strings(x));
    auto fct = SHIELD(r_factors(x, levels));
    this->value = fct.value;
    YIELD(2);
  }

  r_vec<r_str> levels_get() const {
    return r_vec<r_str>(attr::get_attr(this->value, internal::as_r<r_sym>("levels")));
  }

  void levels_set(r_vec<r_str> levels) {
    attr::set_attr(this->value, internal::as_r<r_sym>("levels"), levels);
  }

};

}

#endif
