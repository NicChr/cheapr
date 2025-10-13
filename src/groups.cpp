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

    // Manually do first iteration
    // We have to look at lagged value to check if group IDs are sorted
    int curr_group = p_group_id[0] - 1;
    int curr_group_start = p_out[curr_group];
    p_out[curr_group] = std::min(curr_group_start, 1);

    for (int i = 1; i < n; ++i){
      curr_group = p_group_id[i] - 1;
      curr_group_start = p_out[p_group_id[i] - 1];
      p_out[curr_group] = std::min(curr_group_start, i + 1);
      sorted = sorted && curr_group > (p_group_id[i - 1] - 1);
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

    int curr_group = p_group_id[0] - 1;
    unsigned int curr_group_start = p_out[curr_group];

    p_out[curr_group] = static_cast<int>(
      std::min(curr_group_start - 1U, static_cast<unsigned int>(0)) + 1U
    );

    for (int i = 1; i < n; ++i){
      curr_group = p_group_id[i] - 1;
      curr_group_start = p_out[curr_group];


      p_out[curr_group] = static_cast<int>(
        std::min(curr_group_start - 1U, static_cast<unsigned int>(i)) + 1U
      );
      sorted = sorted && curr_group > (p_group_id[i - 1] - 1);
    }
  }

  SEXP r_sorted = SHIELD(Rf_ScalarLogical(sorted));
  Rf_setAttrib(out, Rf_installChar(make_utf8_char("sorted")), r_sorted);
  YIELD(2);
  return out;
}
