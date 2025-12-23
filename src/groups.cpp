#include "cheapr.h"

[[cpp11::register]]
SEXP cpp_group_starts(SEXP group_id, int n_groups){

  int n = Rf_length(group_id);

  SEXP out = SHIELD(vec::new_vector<int>(n_groups));
  const int* p_group_id = integer_ptr_ro(group_id);
  int* RESTRICT p_out = integer_ptr(out);

  int fill_value = r_limits::r_int_max;

  bool sorted = true;

  if (n < fill_value){
    // Initialise start locations
    r_fill(out, p_out, 0, n_groups, fill_value);

    for (int i = 0; i < n; ++i){
      int curr_group = p_group_id[i] - 1;
      int curr_group_start = p_out[p_group_id[i] - 1];
      p_out[curr_group] = std::min(curr_group_start, i + 1);
      sorted = sorted && p_out[curr_group] >= p_out[std::max(curr_group, 1) - 1];
    }

    // This will set groups with no start locations to 0
    // (e.g. undropped factor levels)
    r_replace(out, p_out, 0, n_groups, fill_value, 0);
  } else {

    // Slightly slower method than above
    // this can handle the edge-case where a group's start location
    // happens to be at .Machine$integer.max

    // Initialise start locations
    r_fill(out, p_out, 0, n_groups, 0);

    for (int i = 0; i < n; ++i){
      int curr_group = p_group_id[i] - 1;
      unsigned int curr_group_start = p_out[curr_group];


      p_out[curr_group] = static_cast<int>(
        std::min(curr_group_start - 1U, static_cast<unsigned int>(i)) + 1U
      );
      sorted = sorted && p_out[curr_group] >= p_out[std::max(curr_group, 1) - 1];
    }
  }

  SEXP r_sorted = SHIELD(as_vector(sorted));
  set_attr(out, r_cast<r_symbol_t>("sorted"), r_sorted);
  YIELD(2);
  return out;
}

SEXP cpp_group_counts(SEXP group_id, int n_groups){

  int n = Rf_length(group_id);

  // Counts initialised to zero
  SEXP out = SHIELD(vec::new_vector<int>(n_groups, 0));
  const int* RESTRICT p_group_id = integer_ptr_ro(group_id);
  int* RESTRICT p_out = integer_ptr(out);

  // Count groups
  for (int i = 0; i < n; ++i) p_out[p_group_id[i] - 1]++;

  YIELD(1);
  return out;
}
