#ifndef CHEAPR_R_VECTOR_H
#define CHEAPR_R_VECTOR_H

#include <cheapr/internal/r_setup.h>
#include <cheapr/internal/r_utf8.h>
#include <cheapr/internal/r_types.h>
#include <cheapr/internal/r_concepts.h>
#include <cheapr/internal/r_attrs.h>
#include <cheapr/internal/r_vector_utils.h>
#include <cheapr/internal/r_coerce.h>
#include <unordered_map> // For string vector matching
#include <unordered_set> // For unique strings
#include <vector> // Also for unique strings

namespace cheapr {

template<RType T>
struct r_vec {
  SEXP value = R_NilValue;
  const T* const_ptr = nullptr;  // Always created
  T* ptr = nullptr;              // Only initialized if writable
  using data_type = T;              // Type of data vec contains

  // Constructor that wraps new_vector<T>
  explicit r_vec(R_xlen_t size)
    : value(internal::new_vector_impl<std::remove_cvref_t<T>>(size))
    , const_ptr(internal::vector_ptr<const T>(value))
    , ptr(nullptr)
  {
    if constexpr (RPtrWritableType<T>) {
      ptr = internal::vector_ptr<T>(value);
    }
  }

  // Constructor from existing SEXP
  explicit r_vec(SEXP s = R_NilValue)
    : value(s)
    , const_ptr(s == R_NilValue ? nullptr : internal::vector_ptr<const T>(s))
    , ptr(nullptr)
  {
    if constexpr (RPtrWritableType<T>){
      if (const_ptr != nullptr){
        ptr = internal::vector_ptr<T>(s);
      }
    }
  }

  bool is_null() const {
    return value == R_NilValue;
  }

    // Direct pointer access
  T* data() requires RPtrWritableType<T>{
    return ptr;
  }

  const T* data() const {
    return const_ptr;
  }

  // Implicit conversion to SEXP
  operator SEXP() const {
    return value;
  }

  // Get size
  R_xlen_t size() const {
    return Rf_xlength(value);
  }
  R_xlen_t length() const {
    return Rf_xlength(value);
  }

  const r_str address() const {
    return internal::address(value);
  }

r_vec<r_lgl> is_na() const {
  R_xlen_t n = length();
  auto out = SHIELD(r_vec<r_lgl>(n));
  
  if constexpr (RPtrWritableType<T>){
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1){
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i){
        out.set(i, is_r_na(get(i)));
      }
    } else {
      OMP_SIMD
      for (R_xlen_t i = 0; i < n; ++i){
        out.set(i, is_r_na(get(i)));
      }
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i){
      out.set(i, is_r_na(get(i)));
    }
  }
  
  YIELD(1);
  return out;
}

    // get uses const pointer
  T get(R_xlen_t index) const {
    return const_ptr[index];
  }

  // Set only available if writable
  // We use flexible template to be able to coerce it to an RType
  template <typename U>
   void set(R_xlen_t index, U val) {
    T val2 = internal::as_r<T>(val);
    if constexpr (RPtrWritableType<T>){
      ptr[index] = val2;
    } else {
      internal::set_value<T>(value, index, val2);
    }
  }

  template <typename U>
  void fill(R_xlen_t start, R_xlen_t n, const U val){
  auto val2 = internal::as_r<T>(val);
  if constexpr (RPtrWritableType<T>){
    int n_threads = internal::calc_threads(n);
    auto *p_target = data();
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        p_target[start + i] = val2;
      }
    } else {
      std::fill_n(p_target + start, n, val2);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      set(start + i, val2);
    }
  }
}

template <typename U1, typename U2> 
void replace(R_xlen_t start, R_xlen_t n, const U1 old_val, const U2 new_val){
  auto old_val2 = internal::as_r<T>(old_val);
  auto new_val2 = internal::as_r<T>(new_val);
  bool implicit_na_coercion = !is_r_na(old_val) && is_r_na(old_val2);
  if (!implicit_na_coercion){
    if constexpr (RPtrWritableType<T>){
      int n_threads = internal::calc_threads(n);
      auto *p_target = data();
      if (n_threads > 1) {
        OMP_PARALLEL_FOR_SIMD(n_threads)
        for (R_xlen_t i = 0; i < n; ++i) {
          if (p_target[start + i] == old_val2){
            p_target[start + i] = new_val2;
          }
        }
      } else {
        std::replace(data() + start, data() + start + n, old_val2, new_val2);
      }
    } else {
      for (R_xlen_t i = 0; i < n; ++i) {
        R_xlen_t idx = start + i;
        if (get(idx) == old_val2){
          set(idx, new_val2);
        }
      }
    }
}
}

};


// R Date vector
struct r_dates : public r_vec<r_int> {
  
  // Constructors
  r_dates() : r_vec<r_int>() {}
  
  explicit r_dates(SEXP x) : r_vec<r_int>(x) {
    if (!is_null() && !attr::inherits1(x, "Date")){ 
      Rf_error("`SEXP` must be a Date");
    }
  }
  
  explicit r_dates(R_xlen_t n) : r_vec<r_int>(n) {
    attr::set_old_class(this->value, internal::make_utf8_strsxp("Date"));
  }
};

// R POSIXct vector
struct r_posixcts : public r_vec<r_dbl> {
  
  // Constructors
  r_posixcts() : r_vec<r_dbl>() {}
  
  explicit r_posixcts(SEXP x) : r_vec<r_dbl>(x) {
    if (!(is_null() || (attr::inherits1(x, "POSIXct") && attr::inherits1(x, "POSIXt")))){ 
      Rf_error("`SEXP` must be a POSIXct");
    }
  }
  
  explicit r_posixcts(R_xlen_t n) : r_vec<r_dbl>(n) {
    auto cls = SHIELD(r_vec<r_str>(2));
    auto tz = SHIELD(internal::new_scalar_vector(blank_r_string));
    cls.set(0, "POSIXct"); cls.set(1, "POSIXt");
    // Set class
    attr::set_old_class(this->value, cls);
    // Set timezone
    attr::set_attr(this->value, internal::as_r<r_sym>("tzone"), tz);
    YIELD(2);
  }

  r_str tzone_get() const {
    auto tz = r_vec<r_str>(attr::get_attr(this->value, internal::as_r<r_sym>("tzone")));
    if (tz.is_null() || tz.length() == 0){
      return blank_r_string;
    } else {
      return tz.get(0);
    }
  }

    void tzone_set(r_str tzone) {
    auto tz = SHIELD(r_vec<r_str>(attr::get_attr(this->value, internal::as_r<r_sym>("tzone"))));
    if (tz.length() != 0){
      tz.set(0, tzone);
      YIELD(1);
    } else {
      auto new_tz = SHIELD(internal::new_scalar_vector(tzone));
      attr::set_attr(this->value, internal::as_r<r_sym>("tzone"), new_tz);
      YIELD(2);
    }
  }

};

r_vec<r_str> unique_strings(r_vec<r_str> x) {
  R_xlen_t n = x.length();
  
  // Hash set for O(n) de-duplication
  std::unordered_set<SEXP> seen;
  seen.reserve(n);
  
  std::vector<r_str> unique_vec;
  unique_vec.reserve(n / 2);
  
  for (R_xlen_t i = 0; i < n; ++i) {
    r_str str = x.get(i);
    if (seen.insert(str).second) {
      unique_vec.push_back(str);
    }
  }
  
  // Result vector
  R_xlen_t n_unique = unique_vec.size();
  auto out = SHIELD(r_vec<r_str>(n_unique));
  
  for (R_xlen_t i = 0; i < n_unique; ++i) {
    out.set(i, unique_vec[i]);
  }
  
  YIELD(1);
  return out;
}


// Match needle strings to first occurrence in haystack
r_vec<r_int> string_match(r_vec<r_str> needles, r_vec<r_str> haystack) {
  
  R_xlen_t n_needles = needles.length();
  R_xlen_t n_haystack = haystack.length();
  
  // Build hash table: O(m) where m = haystack size
  std::unordered_map<SEXP, R_xlen_t> lookup;
  lookup.reserve(n_haystack);
  
  for (R_xlen_t i = 0; i < n_haystack; ++i) {
    r_str str = haystack.get(i);
    // Only store first occurrence
    if (lookup.find(str) == lookup.end()) {
      lookup[str] = i + 1;
    }
  }
  
  // Match needles: O(n) where n = needles size
  auto out = SHIELD(r_vec<r_int>(n_needles));
  for (R_xlen_t i = 0; i < n_needles; ++i) {
    r_str needle = needles.get(i);
    auto it = lookup.find(needle);
    out.set(i, it != lookup.end() ? r_int(it->second) : na::integer);
  }
  
  YIELD(1);
  return out;
}


struct r_factors : public r_vec<r_int> {
  
  // Constructors
  r_factors() : r_vec<r_int>() {}
  
  explicit r_factors(SEXP x) : r_vec<r_int>(x) {
    if (!(is_null() || attr::inherits1(x, "factor"))){ 
      Rf_error("`SEXP` must be a factor");
    }
  }
  
  explicit r_factors(r_vec<r_str> x, r_vec<r_str> levels){

    auto fct = SHIELD(string_match(x, levels));
    auto cls = SHIELD(internal::new_scalar_vector(internal::as_r<r_str>("factor")));
    // Set class
    attr::set_old_class(fct, cls);
    // Set levels
    attr::set_attr(fct, internal::as_r<r_sym>("levels"), levels);
    this->value = fct.value;
    YIELD(2);
  }

  explicit r_factors(r_vec<r_str> x){
    auto levels = SHIELD(unique_strings(x));
    auto fct = SHIELD(r_factors(x, levels));
    this->value = fct.value;
    YIELD(2);
  }

  r_vec<r_str> levels_get() const {
    return r_vec<r_str>(attr::get_attr(this->value, internal::as_r<r_sym>("levels")));
  }

  void levels_set(r_vec<r_str> levels) {
    attr::set_attr(this->value, internal::as_r<r_sym>("levels"), levels);
  }

};

// R data frame
// struct r_df : public r_vec<r_sexp> {
  
//   // Constructors
//   r_df() : r_vec<r_sexp>() {}
  
//   explicit r_df(SEXP x) : r_vec<r_sexp>(x) {
//     if (!is_null() && !attr::inherits1(x, "Date")){ 
//       Rf_error("`SEXP` must be a Date");
//     }
//   }
  
//   explicit r_df(R_xlen_t n) : r_vec<r_sexp>(n) {
//     attr::set_old_class(this->value, internal::make_utf8_strsxp("Date"));
//   }
// };

namespace internal {

template<typename T>
struct is_r_vector : std::false_type {};

template<typename T>
struct is_r_vector<r_vec<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_r_vector_v = is_r_vector<std::remove_cvref_t<T>>::value;

}

template<typename T>
concept RVectorType = internal::is_r_vector_v<T> || is<T, r_dates> || is<T, r_posixcts>;

template <RType T>
inline void r_copy_n(r_vec<T> &target, const T *p_source, R_xlen_t target_offset, R_xlen_t n){
  if constexpr (RPtrWritableType<T>){
    int n_threads = internal::calc_threads(n);
    if (n_threads > 1) {
      OMP_PARALLEL_FOR_SIMD(n_threads)
      for (R_xlen_t i = 0; i < n; ++i) {
        target.set(target_offset + i, p_source[i]);
      }
    } else {
      std::copy_n(p_source, n, target.data() + target_offset);
    }
  } else {
    for (R_xlen_t i = 0; i < n; ++i) {
      target.set(target_offset + i, p_source[i]);
    }
  }
}

namespace vec {

// Templates for creating new vectors (can also be done via r_vec)
template <RType T>
inline r_vec<T> new_vector(R_xlen_t n){
  return r_vec<T>(n);
}

template <RType T, typename U>
inline r_vec<T> new_vector(R_xlen_t n, const U default_value){
  auto out = SHIELD(r_vec<T>(n));
  out.fill(0, n, default_value); 
  YIELD(1);
  return out;
}

inline bool is_object(SEXP x){
  return Rf_isObject(x);
}

inline bool is_atomic(SEXP x){
  return Rf_isVectorAtomic(x);
}

inline bool is_vec(SEXP x){
  return Rf_isVector(x);
}

inline bool is_bare(SEXP x){
  return !is_object(x);
}

inline bool is_logical(SEXP x){
  return TYPEOF(x) == LGLSXP;
}
inline bool is_integer(SEXP x){
  return TYPEOF(x) == INTSXP;
}
inline bool is_double(SEXP x){
  return TYPEOF(x) == REALSXP;
}
inline bool is_integer64(SEXP x){
  return is_double(x) && attr::inherits1(x, "integer64");
}
inline bool is_character(SEXP x){
  return TYPEOF(x) == STRSXP;
}
inline bool is_list(SEXP x){
  return TYPEOF(x) == VECSXP;
}
inline bool is_complex(SEXP x){
  return TYPEOF(x) == CPLXSXP;
}
inline bool is_raw(SEXP x){
  return TYPEOF(x) == RAWSXP;
}
inline bool is_date(SEXP x){
  return attr::inherits1(x, "Date");
}
inline bool is_datetime(SEXP x){
  return attr::inherits1(x, "POSIXt");
}
inline bool is_factor(SEXP x){
  return is_integer(x) && attr::inherits1(x, "factor");
}


template<typename T>
inline auto as_vector(T x){
  if constexpr (RVectorType<T>){
    return x;
  } else if constexpr (RType<T>){
    return new_vector<T>(1, x);
  } else if constexpr (IntegerType<T>){
    if constexpr (internal::can_be_int<T>){
      return new_vector<r_int>(1, r_int(static_cast<int>(x)));
    } else {
      return new_vector<r_dbl>(1, r_dbl(static_cast<double>(x)));
    }
  } else {
    static_assert(
      always_false<T>,
      "Unimplemented `as_vector` specialisation"
    );
  }
}
template<>
inline auto as_vector<bool>(bool x){
  return as_vector(r_lgl(x));
}
template<>
inline auto as_vector<int>(int x){
  return as_vector(r_int(x));
}
template<>
inline auto as_vector<double>(double x){
  return as_vector(r_dbl(x));
}
template<>
inline auto as_vector<const char *>(const char *x){
  return as_vector(r_str(internal::make_utf8_charsxp(x))); 
}

template<>
inline auto as_vector<SEXP>(SEXP x){ 
  switch (TYPEOF(x)){
  case LGLSXP:
  case INTSXP:
  case REALSXP:
  case STRSXP:
  case VECSXP:
  case CPLXSXP:
  case RAWSXP: {
    return x;
  }
  default: {
    return static_cast<SEXP>(new_vector<r_sexp>(1, r_sexp(x)));
  }
  }
}

}

// Powerful and flexible coercion function that can handle many types and convert to R-spcific C++ types and R vectors
template<typename T, typename U>
requires (RType<T> || RVectorType<T> || is<T, SEXP>)
inline T as(U x) {
if constexpr (is<U, T>){
  return x;
} else if constexpr (RType<U> && is<T, SEXP>){
  return static_cast<SEXP>(internal::as_r<r_sexp>(x));
} else if constexpr (RVectorType<U> && is<T, SEXP>){
  return x.value;
} else if constexpr (RVectorType<U> && is<T, r_sexp>){
  return r_sexp(x.value);
} else if constexpr (RVectorType<U> && RType<T>){
  if (x.length() != 1){
    Rf_error("Vector must be length-1 to be coerced to requested type");
  }
  return internal::as_r<T>(x.get(0));
} else if constexpr (RVectorType<T> && RType<U>){  
  using data_t = typename T::data_type;
  return vec::new_vector<data_t>(1, internal::as_r<data_t>(x));
} else if constexpr (RVectorType<U> && RVectorType<T>){
  R_xlen_t n = x.length();
  auto out = SHIELD(T(n));
  using data_t = typename T::data_type;
  OMP_SIMD
  for (R_xlen_t i = 0; i < n; ++i){
  out.set(i, internal::as_r<data_t>(x.get(i)));
  }
  YIELD(1);
  return out;
} else if constexpr (RType<T> && !RVectorType<U>) {
  return internal::as_r<T>(x);
} else {
  static_assert(always_false<T>, "Unsupported type for `as`");
}
}

// Named argument

struct arg {
  const char* name;
  SEXP value;

  explicit arg(const char* n) : name(n), value(R_NilValue) {}
  explicit arg(const char* n, SEXP v) : name(n), value(v) {}

  template<typename T>
  arg operator=(T v) const {
    return arg(name, as<r_sexp>(v)); 
  }
};

// Variadic list constructor
template<typename... Args>
inline r_vec<r_sexp> make_list(Args... args) {
  constexpr int n = sizeof...(args);

  if constexpr (n == 0){
    return r_vec<r_sexp>(0);
  } else {

    auto out = SHIELD(r_vec<r_sexp>(n));

    // Are any args named?
    constexpr bool any_named = (is<Args, arg> || ...);

    auto nms = any_named ? r_vec<r_str>(n) : r_vec<r_str>();
    SHIELD(nms);

    int i = 0;
    (([&]() {
      if constexpr (is<Args, arg>) { 
        out.set(i, args.value);
        nms.set(i, as<r_str>(args.name));
      } else {
        out.set(i, as<r_sexp>(args));
      }
      ++i;
    }()), ...);

    attr::set_old_names(out, nms);
    YIELD(2);
    return out;
  }
}

}

#endif
