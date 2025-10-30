#include "cheapr.h"

// R Function imports

cpp11::function cheapr_sset = cpp11::package("cheapr")["cheapr_sset"];
cpp11::function cheapr_is_na = cpp11::package("cheapr")["is_na"];
cpp11::function cheapr_factor = cpp11::package("cheapr")["factor_"];
cpp11::function base_rep = cpp11::package("base")["rep"];
cpp11::function base_do_call = cpp11::package("base")["do.call"];
cpp11::function base_as_character = cpp11::package("base")["as.character"];
cpp11::function base_paste0 = cpp11::package("base")["paste0"];
cpp11::function cheapr_fast_match = cpp11::package("cheapr")["fast_match"];
cpp11::function cheapr_fast_unique = cpp11::package("cheapr")["fast_unique"];
cpp11::function cheapr_rebuild = cpp11::package("cheapr")["rebuild"];
cpp11::function cheapr_as_df = cpp11::package("cheapr")["as_df"];
cpp11::function base_cast = cpp11::package("cheapr")["base_cast"];
cpp11::function base_assign = cpp11::package("cheapr")["base_assign_at"];
