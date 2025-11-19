#include "cheapr.h"

[[cpp11::register]]
SEXP cpp_group_starts(SEXP group_id, int n_groups){

  int n = Rf_length(group_id);

  SEXP out = SHIELD(new_vec(INTSXP, n_groups));
  const int* p_group_id = INTEGER_RO(group_id);
  int* RESTRICT p_out = INTEGER(out);

  int fill_value = std::numeric_limits<int>::max();

  bool sorted = true;

  if (n < fill_value){
    // Initialise start locations
    std::fill(p_out, p_out + n_groups, fill_value);

    for (int i = 0; i < n; ++i){
      int curr_group = p_group_id[i] - 1;
      int curr_group_start = p_out[p_group_id[i] - 1];
      p_out[curr_group] = std::min(curr_group_start, i + 1);
      sorted = sorted && p_out[curr_group] >= p_out[std::max(curr_group, 1) - 1];
    }

    // This will set groups with no start locations to 0
    // (e.g. undropped factor levels)
    std::replace(p_out, p_out + n_groups, fill_value, 0);
  } else {

    // Slightly slower method than above
    // this can handle the edge-case where a group's start location
    // happens to be at .Machine$integer.max

    // Initialise start locations
    std::fill(p_out, p_out + n_groups, 0);

    for (int i = 0; i < n; ++i){
      int curr_group = p_group_id[i] - 1;
      unsigned int curr_group_start = p_out[curr_group];


      p_out[curr_group] = static_cast<int>(
        std::min(curr_group_start - 1U, static_cast<unsigned int>(i)) + 1U
      );
      sorted = sorted && p_out[curr_group] >= p_out[std::max(curr_group, 1) - 1];
    }
  }

  SEXP r_sorted = SHIELD(as_r_scalar(sorted));
  Rf_setAttrib(out, Rf_installChar(make_utf8_char("sorted")), r_sorted);
  YIELD(2);
  return out;
}

[[cpp11::register]]
SEXP cpp_group_counts(SEXP group_id, int n_groups){

  int n = Rf_length(group_id);

  SEXP out = SHIELD(new_vec(INTSXP, n_groups));
  const int* RESTRICT p_group_id = INTEGER_RO(group_id);
  int* RESTRICT p_out = INTEGER(out);

  // Initialise counts
  std::fill(p_out, p_out + n_groups, 0);
  // Count groups
  for (int i = 0; i < n; ++i) p_out[p_group_id[i] - 1]++;

  YIELD(1);
  return out;
}
