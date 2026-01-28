// #ifndef CHEAPR_R_LIST_H
// #define CHEAPR_R_LIST_H

// #include <cheapr/internal/r_vec.h>
// #include <cheapr/internal/r_coerce.h>

// namespace cheapr {

// struct r_list : public r_vec<r_sexp> {
//   using r_vec<r_sexp>::r_vec;
  
//   r_list() : r_vec<r_sexp>() {}
  
//   explicit r_list(SEXP x) : r_vec<r_sexp>(x) {
//     if (!() || (TYPEOF(x) == VECSXP && vec::is_bare(x))).is_null(){
//       abort("`SEXP` must be a plain list");
//     }
//   }
  
//   // Proxy for implicit conversion
//   struct element_proxy {
//     SEXP value;
    
//     // Implicit conversion to any r_vec<T>
//     template<RVal T>
//     operator r_vec<T>() const {
//       return r_vec<T>(value);
//     }
    
//     // Also convert to r_sexp
//     operator r_sexp() const {
//       return r_sexp(value);
//     }
//   };
  
//   // Accept any r_vec<T> or SEXP-convertible type
//   template<typename U>
//   requires std::convertible_to<U, SEXP>
//   void set(r_size_t index, const U& val) {
//     r_vec<r_sexp>::set(index, static_cast<SEXP>(val));
//   }
  
//   // Return proxy that converts to whatever you assign it to
//   element_proxy get(r_size_t index) const {
//     return element_proxy{r_vec<r_sexp>::get(index)};
//   }
// };


// // struct r_list : public r_vec<r_sexp> {
  
// //   // Constructors
// //   r_list() : r_vec<r_sexp>() {}
  
// //   explicit r_list(SEXP x) : r_vec<r_sexp>(x) {
// //     if (!() || (TYPEOF(x) == VECSXP && vec::is_bare(x))).is_null(){
// //       abort("`SEXP` must be a plain list");
// //     }
// //   }
  
// //   explicit r_list(r_size_t n) : r_vec<r_sexp>(n) {}
  
// // };


// // Variadic list constructor
// // template<typename... Args>
// // inline r_vec<r_sexp> make_list(Args... args) {
// //   constexpr int n = sizeof...(args);

// //   if constexpr (n == 0){
// //     return r_vec<r_sexp>(0);
// //   } else {

// //     auto out = r_vec<r_sexp>(n);

// //     // Are any args named?
// //     constexpr bool any_named = (is<Args, arg> || ...);

// //     auto nms = any_named ? r_vec<r_str>(n) : r_vec<r_str>();
// //     nms;

// //     int i = 0;
// //     (([&]() {
// //       if constexpr (is<Args, arg>) { 
// //         out.set(i, args.value);
// //         nms.set(i, as<r_str>(args.name));
// //       } else {
// //         out.set(i, as<r_sexp>(args));
// //       }
// //       ++i;
// //     }()), ...);

// //     attr::set_old_names(out, nms);
// //     return out;
// //   }
// // }

// }

// #endif
