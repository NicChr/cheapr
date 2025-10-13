#include "cheapr.h"

[[cpp11::register]]
SEXP cpp_group_starts(SEXP group_id, int n_groups){

  int n = Rf_length(group_id);

  SEXP out = SHIELD(new_vec(INTSXP, n_groups));
  const int* p_group_id = INTEGER_RO(group_id);
  int* RESTRICT p_out = INTEGER(out);

  int fill_value = std::numeric_limits<int>::max();

  if (n < fill_value){
    // Initialise start locations
    std::fill(p_out, p_out + n_groups, fill_value);
    for (int i = 0; i < n; ++i){
      p_out[p_group_id[i] - 1] = std::min(p_out[p_group_id[i] - 1], i + 1);
    }

    for (int i = 0; i < n_groups; ++i){
      if (p_out[i] == fill_value){
        p_out[i] = 0;
      }
    }
  } else {

    // Slightly slower method than above
    // this can handle the edge-case where a group's start location
    // happens to be at .Machine$integer.max

    // Initialise start locations
    std::fill(p_out, p_out + n_groups, 0);
    for (int i = 0; i < n; ++i){
      p_out[p_group_id[i] - 1] = static_cast<int>(
        std::min(
          static_cast<unsigned int>(p_out[p_group_id[i] - 1]) - 1U,
          static_cast<unsigned int>(i)
        ) + 1U
      );
    }
  }
  YIELD(1);
  return out;
}

