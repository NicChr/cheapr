#ifndef CHEAPR_R_FACTOR_H
#define CHEAPR_R_FACTOR_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_attrs.h>
#include <ankerl/unordered_dense.h> // For unique strings and matching
#include <vector>


namespace cheapr {

namespace internal {

r_vec<r_str> unique_strings(r_vec<r_str> x) {
  r_size_t n = x.length();
  
  // Hash set for O(n) de-duplication
  ankerl::unordered_dense::set<SEXP> seen;
  seen.reserve(n);
  
  // Store uniques in insertion order using raw SEXP pointers
  std::vector<SEXP> uniques;
  uniques.reserve(n);
  
  auto *p_x = x.data();
  
  for (r_size_t i = 0; i < n; ++i) {
      SEXP str = p_x[i];
      if (seen.insert(str).second) {
          uniques.push_back(str);
      }
  }
  
  r_size_t n_unique = uniques.size();
  
  if (n_unique == n) {
      return x;
  }
  
  auto out = r_vec<r_str>(n_unique);
  for (r_size_t i = 0; i < n_unique; ++i) {
      SET_STRING_ELT(out, i, uniques[i]);
  }
  
  return out;
}


// Match needle strings to first occurrence in haystack
r_vec<r_int> string_match(r_vec<r_str> needles, r_vec<r_str> haystack) {

  r_size_t n_needles = needles.length();
  r_size_t n_haystack = haystack.length();

  if (n_haystack > r_limits<r_int>::max()){
    abort("Cannot match to a long vector, please use a short character vector");
  }

  // Build hash table
  ankerl::unordered_dense::map<SEXP, int> lookup;
  lookup.reserve(n_haystack);

  for (r_size_t i = 0; i < n_haystack; ++i) {
    r_str str = haystack.get(i);
    lookup.try_emplace(str, static_cast<int>(i) + 1);
  }

  // Match needles
  auto out = r_vec<r_int>(n_needles);
  for (r_size_t i = 0; i < n_needles; ++i) {
    r_str needle = needles.get(i);
    auto it = lookup.find(needle);
    out.set(i, it != lookup.end() ? r_int(it->second) : na::integer);
  }

  return out;
}

}


struct r_factors : public r_vec<r_int> { 

  // Constructors
  r_factors() : r_vec<r_int>() {}

  explicit r_factors(SEXP x) : r_vec<r_int>(x) {
    if (!is_null() && !attr::inherits1(*this, "factor")){
      abort("`SEXP` must be a factor");
    }
  }

  explicit r_factors(r_vec<r_str> x, r_vec<r_str> levels){

    auto fct = internal::string_match(x, levels);
    auto cls = r_vec<r_str>(1, internal::as_r<r_str>("factor"));
    // Set class
    attr::set_old_class(fct, cls);
    // Set levels
    attr::set_attr(fct, internal::as_r<r_sym>("levels"), levels.sexp);
    this->sexp = fct.sexp;
  }

  explicit r_factors(r_vec<r_str> x){
    auto levels = internal::unique_strings(x);
    auto fct = r_factors(x, levels);
    this->sexp = fct.sexp;
  }

  r_vec<r_str> levels() const {
    return r_vec<r_str>(attr::get_attr(this->sexp, internal::as_r<r_sym>("levels")));
  }

  void set_levels(r_vec<r_str> levels) {
    attr::set_attr(this->sexp, internal::as_r<r_sym>("levels"), levels.sexp);
  }

};

}

#endif
