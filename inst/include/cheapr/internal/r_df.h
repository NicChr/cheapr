#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_attrs.h>
#include <cheapr/internal/r_coerce.h>
#include <cheapr/internal/r_make_vec.h>

namespace cheapr {

namespace internal {
r_vec<r_sexp> new_df_impl(r_vec<r_sexp> cols){

    int n = cols.length();
    int nrows;

    if (n == 0){
        nrows = 0;
    } else {
        nrows = cols.get(0).length();
    }

    r_vec<r_sexp> out(n);

    r_vec<r_str> names = cols.names();
    if (names.is_null()){
        names = r_vec<r_str>(n, blank_r_string);
    }

    // Copy cols to new list
    for (int i = 0; i < n; ++i){
        out.set(i, cols.get(i));
    }

    r_vec<r_int> row_names;

    if (nrows != 0){ 
        row_names = make_vec<r_int>(na::integer, -nrows);
    }
    out.set_names(names); 
    attr::set_attr(out, symbol::class_sym, r_vec<r_str>(1, r_str("data.frame")));
    attr::set_attr(out, symbol::row_names_sym, row_names);
    return out;
}
r_vec<r_sexp> new_df_impl(int nrows){

    r_vec<r_sexp> out{};
    r_vec<r_int> row_names;

    if (nrows != 0){
        row_names = make_vec<r_int>(na::integer, -nrows);
    }
    out.set_names(r_vec<r_str>(1, blank_r_string));
    attr::set_attr(out, symbol::class_sym, r_vec<r_str>(1, r_str("data.frame")));
    attr::set_attr(out, symbol::row_names_sym, row_names);
    return out;
}
}

// Forward declare to allow r_df to contain r_df cols
struct r_df; 

struct r_df {

    r_sexp sexp;
    
    // Default constructor (empty data frame)
    r_df() : sexp(std::move(internal::new_df_impl(0).sexp)) {}

    // Constructor from existing SEXP
    explicit r_df(SEXP s) : sexp(s) {
        if (!sexp.is_null() && !attr::inherits1(sexp, "data.frame")) {
            abort("Object is not a data.frame");
        }
    }
    explicit r_df(r_sexp s) : sexp(std::move(s)) {
        if (!sexp.is_null() && !attr::inherits1(sexp, "data.frame")) {
            abort("Object is not a data.frame");
        }
    }

    // Constructor creating a new DF with N rows and M columns
//     r_df(r_size_t n_rows, r_size_t n_cols) 
//         : sexp(internal::new_df_impl(n_rows, n_cols)) {}

    // Implicit conversion to SEXP
    operator SEXP() const { return sexp; }

    bool is_null() const { return sexp.is_null(); }

    r_vec<r_str> rownames() const {
        return as<r_vec<r_str>>(attr::get_attr(sexp, symbol::row_names_sym));
    }
    r_vec<r_str> colnames() const {
        return attr::get_old_names(sexp);
    }
    
    r_size_t nrow() const noexcept {
        return rownames().length();
    }
    r_size_t ncol() const noexcept {
        return sexp.length();
    }
};

}
