#ifndef CHEAPR_R_UNIQUE_H
#define CHEAPR_R_UNIQUE_H

#include <cheapr/internal/r_vec.h>
#include <ankerl/unordered_dense.h> // For unique strings and matching
#include <vector>

namespace cheapr {

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

}

#endif
