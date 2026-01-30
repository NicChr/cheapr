#ifndef CHEAPR_DISPATCH_HPP
#define CHEAPR_DISPATCH_HPP

#include <Rinternals.h>
#include <cheapr/internal/r_coerce.h>
#include <tuple>
#include <utility>

namespace cheapr {

namespace detail {

template <typename T>
T as_cpp(SEXP x) {
    return as<T>(x);
}

// TRAITS: Deduce argument types from a function pointer
template <typename> struct fn_traits;

template <typename Ret, typename... Args>
struct fn_traits<Ret(*)(Args...)> {
    using return_type = Ret;
    using args_tuple = std::tuple<Args...>;
    static constexpr size_t arity = sizeof...(Args);
};

// CALLER: Unpacks SEXPs, casts them, calls Fn, wraps result
template <auto Fn, typename... Args, size_t... Is>
SEXP invoke_impl(SEXP* sexp_args, std::index_sequence<Is...>) {
    // Fold expression to cast each SEXP at index I to the target type Arg
    // Note: We use the raw SEXP array from the generated wrapper
    if constexpr (std::is_void_v<typename fn_traits<decltype(Fn)>::return_type>) {
        Fn(as_cpp<std::decay_t<Args>>(sexp_args[Is])...);
        return R_NilValue;
    } else {
        return unwrap(as_sexp(
            Fn(as_cpp<std::decay_t<Args>>(sexp_args[Is])...)
        ));
    }
}

}

// THE PUBLIC INTERFACE
// Generated code calls this single template
template <auto Fn, typename... SexpArgs>
SEXP dispatch(SexpArgs... args) {
    using Traits = detail::fn_traits<decltype(Fn)>;
    using ArgsTuple = typename Traits::args_tuple;
    
    // Pack SEXPs into an array for indexed access
    SEXP arg_array[] = {args...};
    
    return std::apply([&](auto... target_type_dummies) {
        return detail::invoke_impl<Fn, decltype(target_type_dummies)...>(
            arg_array, 
            std::make_index_sequence<sizeof...(SexpArgs)>{}
        );
    }, ArgsTuple{});
}
}

#endif
