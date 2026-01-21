#include <cheapr/internal/c_core.h>

namespace cheapr {

// Defining custom R types
// for casting, initialising, combining and assigning
// Here R type must be an object that directly interfaces with R
// Things like `NULL`, vectors, factors, data frames, etc
// Which is different from the lower-level RVal concept for C/C++ we have built

namespace internal {
using r_type = uint8_t;
}

using internal::r_type;

struct r_unknown_t {};

// r type constants
enum : r_type {
  R_null = 0,
    R_lgl = 1,
    R_int = 2,
    R_int64 = 3,
    R_dbl = 4,
    R_cplx = 5,
    R_raw = 6,
    R_date = 7,
    R_pxt = 8,
    R_chr = 9,
    R_fct = 10,
    R_list = 11,
    R_df = 12,
    R_unk = 13
};

inline constexpr int n_types = 14;

// R type chars
inline constexpr const char* r_type_names[14] = {
  "NULL",        // 0
  "logical",     // 1
  "integer",     // 2
  "integer64",   // 3
  "numeric",     // 4
  "complex",     // 5
  "raw",         // 6
  "Date",        // 7
  "POSIXct",     // 8
  "character",  // 9
  "factor",     // 10
  "list",       // 11
  "data.frame", // 12
  "unknown"     // 13
};

// An n x n matrix of r types and their common cast type

inline constexpr r_type r_type_pairs[14][14] = {
  /*                0-NULL  1-LGL   2-INT   3-I64   4-DBL    5-CPLX  6-RAW   7-DATE  8-PXCT  9-CHR   10-FCT  11-LIST 12-DF 13-UNK */
  /* 0 - NULL */  { R_null, R_lgl,  R_int,  R_int64, R_dbl,  R_cplx, R_raw,  R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 1 - LGL  */  { R_lgl,  R_lgl,  R_int,  R_int64, R_dbl,  R_cplx, R_raw,  R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 2 - INT  */  { R_int,  R_int,  R_int,  R_int64, R_dbl,  R_cplx, R_raw,  R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 3 - I64  */  { R_int64,R_int64, R_int64, R_int64, R_dbl, R_cplx, R_raw, R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 4 - DBL  */  { R_dbl,  R_dbl,  R_dbl,  R_dbl,   R_dbl,  R_cplx, R_raw, R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 5 - CPLX */  { R_cplx, R_cplx, R_cplx, R_cplx,  R_cplx, R_cplx, R_raw, R_date,  R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 6 - RAW  */  { R_raw,  R_raw,  R_raw,  R_raw, R_raw,  R_raw,  R_raw, R_unk, R_unk, R_chr,  R_fct, R_list, R_df, R_unk },
  /* 7 - DATE */  { R_date, R_date, R_date, R_date,  R_date, R_date,  R_unk, R_date, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 8 - PXCT */  { R_pxt, R_pxt, R_pxt, R_pxt,  R_pxt, R_pxt, R_unk, R_pxt, R_pxt, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 9 - CHR  */  { R_chr,  R_chr,  R_chr,  R_chr, R_chr,  R_chr,  R_chr,  R_chr,  R_chr, R_chr,  R_fct,  R_list, R_df, R_unk },
  /* 10 - FCT  */  { R_fct,  R_fct,  R_fct,  R_fct, R_fct,  R_fct,  R_fct, R_fct,  R_fct, R_fct,  R_fct,  R_list, R_df, R_unk },
  /* 11 - LIST */  { R_list, R_list, R_list, R_list,  R_list, R_list, R_list, R_list, R_list, R_list, R_list, R_list, R_df, R_unk },
  /* 12 - DF */    { R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_df, R_unk },
  /* 13 - Unknown */ { R_unk,  R_unk,  R_unk,  R_unk, R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk,  R_unk }
};

inline r_type common_type(const r_type a, const r_type b) {
  return r_type_pairs[a][b];
}

template <typename T>
inline r_type get_r_type(T x) {
  static_assert(always_false<T>, "Unsupported type for `get_r_type`");
  return R_unk;
}

template <>
inline r_type get_r_type<r_vec<r_lgl>>(r_vec<r_lgl> x) {
    return R_lgl;
}
template <>
inline r_type get_r_type<r_vec<r_int>>(r_vec<r_int> x) {
    return R_int;
}
template <>
inline r_type get_r_type<r_vec<r_int64>>(r_vec<r_int64> x) {
    return R_int64;
}

}
