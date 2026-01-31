#ifndef CHEAPR_R_MATCH_H
#define CHEAPR_R_MATCH_H

#include <cheapr/internal/r_vec.h>
#include <ankerl/unordered_dense.h> // For unique strings and matching

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

}

#endif
