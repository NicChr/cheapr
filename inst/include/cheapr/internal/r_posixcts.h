#ifndef CHEAPR_R_POSIXCTS_H
#define CHEAPR_R_POSIXCTS_H

#include <cheapr/internal/r_vec.h>
#include <cheapr/internal/r_attrs.h>

namespace cheapr {

// R POSIXct vector
struct r_posixcts : public r_vec<r_dbl> {

    private:
    void init_posixct_attrs() {
        auto cls = r_vec<r_str>(2);
        auto tz = r_vec<r_str>(1);
        cls.set(0, "POSIXct"); cls.set(1, "POSIXt");
        // Set class
        attr::set_old_class(*this, cls);
        // Set timezone
        attr::set_attr(*this, internal::as_r<r_sym>("tzone"), tz);
    }

    public: 

    // Constructors
    r_posixcts() : r_vec<r_dbl>() {
        init_posixct_attrs();
    }

    explicit r_posixcts(SEXP x) : r_vec<r_dbl>(x) {
    if (!is_null() && !(attr::inherits1(this->sexp, "POSIXct") && attr::inherits1(this->sexp, "POSIXt"))){
        abort("`SEXP` must be a POSIXct");
      }
    }

    explicit r_posixcts(r_size_t n) : r_vec<r_dbl>(n) {
        init_posixct_attrs();
    }
    template <typename U>
    explicit r_posixcts(r_size_t n, U default_value) : r_vec<r_dbl>(n, default_value) {
        init_posixct_attrs();
    }

    r_str tzone() const {
    auto tz = r_vec<r_str>(attr::get_attr(*this, internal::as_r<r_sym>("tzone")));
    if (tz.is_null() || tz.length() == 0){
        return blank_r_string;
    } else {
        return tz.get(0);
    }
    }

    void set_tzone(r_str tzone) {
    auto tz = r_vec<r_str>(attr::get_attr(*this, internal::as_r<r_sym>("tzone")));
    if (tz.length() != 0){
        tz.set(0, tzone);
    } else {
        auto new_tz = r_vec<r_str>(1, tzone);
        attr::set_attr(*this, internal::as_r<r_sym>("tzone"), new_tz);
    }
    }

};

}

#endif
