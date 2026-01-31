#ifndef CHEAPR_R_FACTOR_H
#define CHEAPR_R_FACTOR_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_attrs.h>
#include <ankerl/unordered_dense.h> // For unique strings and matching
#include <vector>

namespace cheapr {

template <RScalar T>
r_vec<r_int> match(r_vec<T> needles, r_vec<T> haystack) {

  r_size_t n_needles = needles.length();
  r_size_t n_haystack = haystack.length();

  if (n_haystack > r_limits<r_int>::max()){
    abort("Cannot match to a long vector, please use a short vector");
  }

  // Build hash table
  ankerl::unordered_dense::map<unwrap_t<T>, int> lookup;
  lookup.reserve(n_haystack);

  for (r_size_t i = 0; i < n_haystack; ++i) {
    auto elem = unwrap(haystack.get(i));
    lookup.try_emplace(elem, static_cast<int>(i) + 1);
  }

  // Match needles
  auto out = r_vec<r_int>(n_needles);
  for (r_size_t i = 0; i < n_needles; ++i) {
    auto needle = unwrap(needles.get(i));
    auto it = lookup.find(needle);
    out.set(i, it != lookup.end() ? r_int(it->second) : na_value<r_int>());
  }

  return out;
}

template <RScalar T>
r_vec<T> unique(r_vec<T> x) {
  r_size_t n = x.length();
  
  // Hash set for O(n) de-duplication
  ankerl::unordered_dense::set<unwrap_t<T>> seen;
  seen.reserve(n);
  
  // Store uniques in insertion order using raw SEXP pointers
  std::vector<unwrap_t<T>> uniques;
  uniques.reserve(n);
  
  auto *p_x = x.data();
  
  for (r_size_t i = 0; i < n; ++i) {
      auto elem = p_x[i];
      if (seen.insert(elem).second) {
          uniques.push_back(elem);
      }
  }
  
  r_size_t n_unique = uniques.size();
  
  if (n_unique == n) {
      return x;
  }
  
  auto out = r_vec<T>(n_unique);
  for (r_size_t i = 0; i < n_unique; ++i) {
      out.set(i, T(uniques[i]));
  }
  
  return out;
}

// Specialisation for strings
r_vec<r_str> unique(r_vec<r_str> x) {
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
  requires (!is<T, r_cplx>)
  explicit r_factors(r_vec<T> x) : r_factors(x, unique(x)) {}
};

}

#endif
