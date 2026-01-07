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
  auto out = SHIELD(r_posixcts(10));
  out.fill(0, out.length(), 0);
  out.tzone_set(as<r_str>("UTC"));
  YIELD(1);
  return out;
}



[[cpp11::register]]
SEXP ok() {
  r_lgl out = r_na;
  return make_list(out, out.is_false(), out.is_true(), out.is_na());
}

[[cpp11::register]]
SEXP yeah(SEXP x) {
  auto out = internal::visit_vector(x, [&](auto xvec) -> r_factors { 
    return SHIELD(as<r_factors>(xvec));
  });
  YIELD(1);
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
    return SHIELD(as<r_vec<r_str>>(xvec));
  });
  auto out = SHIELD(internal::unique_strings(x2));
  YIELD(2);
  return out;
}


[[cpp11::register]]
SEXP foo4(SEXP x) {
  auto out = internal::visit_vector(x, [&](auto xvec) -> r_vec<r_str> { 
    return SHIELD(as<r_vec<r_str>>(xvec));
  });
  YIELD(1);
  return out;
}

[[cpp11::register]]
SEXP foo5(SEXP x) {

  auto xvec = r_vec<r_str>(x);
  auto out = SHIELD(as<r_vec<r_dbl>>(xvec));
  YIELD(1);
  return out;
}


[[cpp11::register]]
SEXP foo6(SEXP x){
    auto out = internal::visit_vector(x, [&](auto xvec) -> r_posixcts { 
    return SHIELD(as<r_posixcts>(xvec));
  });
  YIELD(1);
  return out;
}
