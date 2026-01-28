#include <cheapr/internal/core.h>
#include <cheapr/internal/test.h>

using namespace cheapr;
using namespace vec;

[[cpp11::register]]
SEXP foo(SEXP x) {
  // auto ok = internal::as_r_string(r_true);
  // return as_vector(as<r_string_t>(r_true));
  // return r_vector_t<r_lgl>(x);
// return as<r_vector_t<r_string_t>>(r_vector_t<r_lgl>(x));

  SEXP out = internal::visit_vector(x, [&](auto xvec) -> SEXP {
    return as<r_vec<r_str>>(xvec);
  });
  return out;

  // return as<r_vector_t<r_string_t>>(r_vector_t<r_lgl>(x));

  // // return make_list(as<r_string_t>(1.5), as<r_string_t>(na::real));
  // return as_vector(as<r_string_t>(vec::new_vector<r_complex_t>(1, 123)));
}

[[cpp11::register]]
SEXP bar(SEXP x) {
  SEXP out = internal::visit_vector(x, [&](auto xvec) -> SEXP {
    return as<r_vec<r_int>>(xvec);
  });
  return out;
}


[[cpp11::register]]
SEXP foobar(SEXP x) {
    auto out = internal::visit_vector(x, [&](auto xvec) -> r_vec<r_lgl> {
      return xvec.is_na();
  });
  return out;
}



[[cpp11::register]]
SEXP foofoo() {
  auto out = r_posixcts(10);
  out.fill(0, out.length(), 0);
  out.set_tzone(as<r_str>("UTC"));

  return out;
}



[[cpp11::register]]
SEXP ok() {
  r_lgl out = r_na;
  return make_vec<r_sexp>(out, out.is_false(), out.is_true(), out.is_na());
}

[[cpp11::register]]
SEXP yeah(SEXP x) {
  auto out = internal::visit_vector(x, [&](auto xvec) -> r_factors {
    return as<r_factors>(xvec);
  });

  return out;
}

[[cpp11::register]]
SEXP yeah2(SEXP x) {
  return internal::unique_strings(r_vec<r_str>(x));
}


[[cpp11::register]]
SEXP foo2() {
  return as<r_vec<r_int>>(r_int(1));
}


[[cpp11::register]]
SEXP foo3(SEXP x) {
  auto x2 = internal::visit_vector(x, [&](auto xvec) -> r_vec<r_str> {
    return as<r_vec<r_str>>(xvec);
  });
  auto out = internal::unique_strings(x2);

  return out;
}


[[cpp11::register]]
SEXP foo4(SEXP x) {
  auto out = internal::visit_vector(x, [&](auto xvec) -> r_vec<r_str> {
    return as<r_vec<r_str>>(xvec);
  });

  return out;
}

[[cpp11::register]]
SEXP foo5(SEXP x) {

  auto xvec = r_vec<r_str>(x);
  auto out = as<r_vec<r_dbl>>(xvec);

  return out;
}


[[cpp11::register]]
SEXP foo6(SEXP x){
    auto out = internal::visit_vector(x, [&](auto xvec) -> r_posixcts {
    return as<r_posixcts>(xvec);
  });

  return out;
}


[[cpp11::register]]
SEXP foo7(SEXP x){
    auto out = internal::visit_vector(x, [&](auto xvec) -> r_vec<r_str> {
    return xvec.names();
  });

  return out;
}


[[cpp11::register]]
SEXP foo8(){
  auto x = r_vec<r_sexp>(1);
  auto y = r_vec<r_str>(3);
  x.set(0, y);

  return x;
}


[[cpp11::register]]
SEXP foo9(){
  auto x = make_vec<r_int>(1, 2, 3);
  auto y = make_vec<r_dbl>(arg("x") = 2.5);
  auto out = make_vec<r_sexp>(arg("first") = x, arg("second") = y);

  return out;
}

[[cpp11::register]]
SEXP foo10(SEXP x){
  auto y = r_vec<r_int>(x);
  return y.resize(10);
}


[[cpp11::register]]
SEXP foo11(SEXP x){
  auto y = r_vec<r_int>(x);
  return attr::get_attrs(y.sexp);
}


[[cpp11::register]]
SEXP foo12(SEXP x, SEXP attrs){
  auto y = r_vec<r_int>(x);
  attr::set_attrs(y.sexp, as<r_vec<r_sexp>>(attrs));
  return y;
}

[[cpp11::register]]
SEXP foo13(SEXP x){
  auto y = r_vec<r_int>(x);
  attr::modify_attrs(y.sexp, arg("ok") = make_vec<r_int>(1, 2, 3));
  return y;
}

// Go from an r_sexp to an r_vec<>

[[cpp11::register]]
SEXP foo14(){
  r_sexp x = r_sexp(r_vec<r_int>(10));
  auto out = as<r_vec<r_str>>(x);

  return out;
}


[[cpp11::register]]
SEXP foo15(){
  return r_vec<r_str>(10, "Hello");
}

[[cpp11::register]]
SEXP foo16(SEXP x){
  return as_vector(as<r_int>(x));
}

[[cpp11::register]]
SEXP foo17(SEXP x){
  return r_vec<r_sexp>(1, r_sexp(x));
}

[[cpp11::register]]
SEXP foo18(SEXP x){
  return x;
  // auto y = r_vec<r_int>(10, 3);
  // auto out = r_list(1, y);
  // auto out2 = out.get(0);
  //
  // return out2;
}


// Compilation error (expected)
// SEXP foo19(SEXP x){
//   return as_vector(x);
// }

[[cpp11::register]]
double foo20(){
  int x = std::abs(na::integer);
  return x == NA_INTEGER ? na::real : as<r_dbl>(x);
}

[[cpp11::register]]
SEXP foo21(int n){
  return r_vec<r_int>(n, 0);
}
[[cpp11::register]]
SEXP foo22(SEXP cols){
  auto new_cols = r_vec<r_sexp>(cols);
  return new_cols.get(1);
  // return internal::new_df_impl(r_vec<r_sexp>(cols));
  // r_vec<r_int> row_names;
  // return row_names;
}
[[cpp11::register]]
SEXP foo23(){
  return as<r_vec<r_str>>("data.frame");
}

[[cpp11::register]]
SEXP foo24(SEXP cols){
  auto new_cols = r_vec<r_sexp>(cols);
  return new_cols.get(1);
}

[[cpp11::register]]
SEXP foo25(SEXP cols, SEXP x){
  auto new_cols = r_vec<r_sexp>(cols);
  new_cols.set(1, r_sexp(x));
  return new_cols;
}

[[cpp11::register]]
SEXP foo26(){
  return make_vec<r_sexp>(sizeof(r_sexp), sizeof(SEXP), sizeof(cpp11::sexp), sizeof(r_str), 
arg("dbl") = sizeof(r_dbl), arg("int") = sizeof(r_int), arg("lgl") = sizeof(r_lgl), arg("cplx") = sizeof(r_cplx));
}


[[cpp11::register]]
SEXP foo27(SEXP cols){
  return internal::new_df_impl(r_vec<r_sexp>(cols));
}


[[cpp11::register]]
SEXP foo28(SEXP x, SEXP y){
  r_vec<r_sexp> x_ = r_vec<r_sexp>(x);
  r_sexp y_ = r_sexp(y);
  int n = x_.length();
  for (int i = 0; i < n; ++i){
    x_.set(i, r_sexp(y_));
  }
  return x_;
}

[[cpp11::register]]
SEXP foo29(SEXP x){
  return internal::visit_vector(x, [&](auto xvec) -> SEXP {
    return xvec.is_na();
  });
}

[[cpp11::register]]
int foo30(){
  r_lgl x = r_true;
  // return !x;
  if (x != r_false){
    return 1;
  } else {
    return 2;
  }
}

[[cpp11::register]]
SEXP foo31(){
  return make_vec<r_dbl>(arg("x") = 2.5, arg("y") = 3.5);
}
[[cpp11::register]]
SEXP foo32(){
  return make_vec<r_int>(1, 2, 3);
}

[[cpp11::register]]
SEXP foo33() {
  return r_dates();
}
[[cpp11::register]]
SEXP foo34() {
  return r_dates(10);
}
[[cpp11::register]]
SEXP foo35() {
  return r_dates(10, 3);
}

[[cpp11::register]]
SEXP foo36(int n) {
  return r_dates(n, 0);
}

[[cpp11::register]]
SEXP foo37(SEXP x, int n) {
  return internal::visit_vector(x, [&](auto xvec) -> SEXP {
    return xvec.rep_len(n);
  });
}


[[cpp11::register]]
SEXP foo38(SEXP x, bool na_rm){
  r_vec<r_dbl> y = r_vec<r_dbl>(x);
  return as_vector(sum(y, na_rm));
}


[[cpp11::register]]
SEXP foo39(SEXP x, bool na_rm){
  return internal::visit_vector(x, [&](auto xvec) -> SEXP {
    using t = decltype(xvec);
    if constexpr (is<t, r_vec<r_int>> || is<t, r_vec<r_dbl>> || is<t, r_vec<r_lgl>>){
      return range(xvec, na_rm);
    } else {
      return r_null;
    }
  });
}

[[cpp11::register]]
SEXP foo40(SEXP x, bool na_rm){
  return internal::visit_vector(x, [&](auto xvec) -> SEXP {
    using t = decltype(xvec);
    if constexpr (is<t, r_vec<r_int>> || is<t, r_vec<r_dbl>> || is<t, r_vec<r_lgl>>){
      return as_vector(min(xvec, na_rm));
    } else {
      Rf_error("error");
      return r_null;
    }
  });
} 

[[cpp11::register]]
SEXP foo41(SEXP x){
  r_vec<r_int> y = r_vec<r_int>(x);
  return as<r_vec<r_dbl>>(sum_int(y));
}


[[cpp11::register]]
SEXP foo42(SEXP x) {
  auto xvec = r_vec<r_dbl>(x);

    r_size_t n = xvec.length();

    r_dbl sum(0);

    for (r_size_t i = 0; i < n; ++i){
      sum += xvec.get(i);
    }
    return as_vector(sum);
}


[[cpp11::register]]
SEXP foo43(SEXP x) {
  auto xvec = r_vec<r_dbl>(x);

    r_size_t n = xvec.length();
    const auto *p_x = xvec.data();

    r_dbl sum(0);

    for (r_size_t i = 0; i < n; ++i){
      sum += p_x[i];
    }
    return as_vector(sum);
}


[[cpp11::register]]
SEXP foo44(SEXP x) {
  auto xvec = r_vec<r_dbl>(x);

    r_size_t n = xvec.length();
    const auto* RESTRICT p_x = xvec.data();

    r_dbl sum(0);

    for (r_size_t i = 0; i < n; ++i){
      sum += p_x[i];
    }
    return as_vector(sum);
}


[[cpp11::register]]
SEXP foo45(SEXP x) {
  auto xvec = r_vec<r_dbl>(x);

    r_size_t n = xvec.length();
    const auto* RESTRICT p_x = xvec.data();

    double sum(0);

    for (r_size_t i = 0; i < n; ++i){
      sum += p_x[i].value;
    }
    return as_vector(r_dbl(sum));
}

[[cpp11::register]]
SEXP foo46() {
  r_str x("hi");
  auto y = unwrap(x);
  return y;
}

[[cpp11::register]]
SEXP foo47(int n) {
  return r_vec<r_int>(n, 0);
}



[[cpp11::register]]
SEXP foo48(SEXP x) {
  auto xvec = r_vec<r_int>(x);

    r_size_t n = xvec.length();

    r_int min_ = r_limits<r_int>::max();

    for (r_size_t i = 0; i < n; ++i){
      min_ = min(min_, xvec.get(i)); 
    }
    return as_vector(min_);
}

[[cpp11::register]]
SEXP foo49(SEXP x, bool na_rm){
  return internal::visit_vector(x, [&](auto xvec) -> SEXP {
    using t = decltype(xvec);
    if constexpr (is<t, r_vec<r_int>> || is<t, r_vec<r_dbl>> || is<t, r_vec<r_lgl>>){
      return as_vector(sum(xvec, na_rm));
    } else {
      Rf_error("error");
      return r_null;
    }
  });
} 

[[cpp11::register]]
SEXP foo50(SEXP x, bool na_rm){
  r_vec<r_int> x_ = r_vec<r_int>(x);

  return as<r_vec<r_dbl>>(sum_int(x_, na_rm));
} 


[[cpp11::register]]
SEXP foo51(SEXP x){
  r_vec<r_int> x_ = r_vec<r_int>(x);
  return abs(x_); 
} 


[[cpp11::register]]
bool foo52(){
  if constexpr (ConstructibleToRVal<decltype("yes")>){
    return true;
  } else {
   return false; 
  }
}

[[cpp11::register]]
bool foo53(){
  if constexpr (ConstructibleToRVal<int16_t>){
    return true;
  } else {
   return false; 
  }
}


[[cpp11::register]]
SEXP foo54(SEXP x){
  auto y = r_vec<r_int>(x);
  return as_r_val(y);
}


[[cpp11::register]]
SEXP foo55(SEXP x, SEXP y){
  auto x_ = r_vec<r_int>(x);
  auto y_ = r_vec<r_int>(y);  
  r_size_t n = x_.length();
  auto z = r_vec<r_int>(n);

  OMP_SIMD
  for (r_size_t i = 0; i < n; ++i){
    z.set(i, x_.get(i) + y_.get(i));
  }
  return z;
}


[[cpp11::register]]
SEXP foo56(SEXP x, SEXP y){
  auto x_ = r_vec<r_dbl>(x);
  auto y_ = r_vec<r_dbl>(y);  
  r_size_t n = x_.length();
  auto z = r_vec<r_dbl>(n);

  OMP_SIMD
  for (r_size_t i = 0; i < n; ++i){
    z.set(i, x_.get(i) + y_.get(i));
  }
  return z;
}

[[cpp11::register]]
SEXP foo57(SEXP x, SEXP y){
  auto x_ = r_vec<r_dbl>(x);
  auto y_ = r_vec<r_dbl>(y);  
  r_size_t n = x_.length();
  auto z = r_vec<r_dbl>(n);

  auto* RESTRICT p_x = x_.data();
  auto* RESTRICT p_y = y_.data();
  auto* RESTRICT p_z = z.data();

  OMP_SIMD
  for (r_size_t i = 0; i < n; ++i){
    p_z[i].value = p_x[i].value + p_y[i].value;
  }
  return z;
}
