#ifndef CHEAPR_R_LIST_H
#define CHEAPR_R_LIST_H

#include <cheapr/internal/r_vector.h>
#include <cheapr/internal/r_coerce.h>

namespace cheapr {

// struct r_list {
//   SEXP value;
  
//   // Constructors
//   r_list() : value(R_NilValue) {}
  
//   explicit r_list(SEXP x) : value(x) {
//     if (!(is_null() || (TYPEOF(x) == VECSXP && vec::is_bare(x)))){
//       Rf_error("`SEXP` must be a plain list");
//     }
//   }
  
//   explicit r_list(r_size_t n) {
//     value = r_vec<r_sexp>(n);
//   }
  
//   // Length
//   r_size_t length() const {
//     return Rf_xlength(value);
//   }
  
//   bool is_null() const {
//     return value == R_NilValue;
//   }
  
//   // Get element as SEXP
//   SEXP get(r_size_t i) const {
//     return VECTOR_ELT(value, i);
//   }
  
//   // Set element from SEXP
//   void set(r_size_t i, SEXP val) {
//     SET_VECTOR_ELT(value, i, val);
//   }
  
//   // Set element from r_vec<> or any wrapper with .value member
//   template<typename T>
//   void set(r_size_t i, const T& vec) {
//     SET_VECTOR_ELT(value, i, vec.value);
//   }
  
//   // Get element as specific type
//   template<typename T>
//   T get_as(r_size_t i) const {
//     return T(VECTOR_ELT(value, i));
//   }
  
//   // Conversion operator
//   operator SEXP() const {
//     return value;
//   }
// };


// Variadic list constructor
// template<typename... Args>
// inline r_vec<r_sexp> make_list(Args... args) {
//   constexpr int n = sizeof...(args);

//   if constexpr (n == 0){
//     return r_vec<r_sexp>(0);
//   } else {

//     auto out = SHIELD(r_vec<r_sexp>(n));

//     // Are any args named?
//     constexpr bool any_named = (is<Args, arg> || ...);

//     auto nms = any_named ? r_vec<r_str>(n) : r_vec<r_str>();
//     SHIELD(nms);

//     int i = 0;
//     (([&]() {
//       if constexpr (is<Args, arg>) { 
//         out.set(i, args.value);
//         nms.set(i, as<r_str>(args.name));
//       } else {
//         out.set(i, as<r_sexp>(args));
//       }
//       ++i;
//     }()), ...);

//     attr::set_old_names(out, nms);
//     YIELD(2);
//     return out;
//   }
// }

}

#endif
