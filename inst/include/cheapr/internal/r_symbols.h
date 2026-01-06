#ifndef CHEAPR_R_SYMBOLS_H
#define CHEAPR_R_SYMBOLS_H

#include <cheapr/internal/r_types.h>

namespace cheapr {

namespace internal {
inline SEXP BASE_ATTRIBUTES = NULL;
inline SEXP BASE_LENGTH = NULL;
inline SEXP r_as_lgl = NULL;
inline SEXP r_as_int = NULL;
inline SEXP r_as_dbl = NULL;
inline SEXP r_as_char = NULL;
inline SEXP r_as_cplx = NULL;
inline SEXP r_as_raw = NULL;
inline SEXP r_as_date = NULL;
inline SEXP r_as_posixct = NULL;
inline SEXP r_as_list = NULL;
}

namespace symbol {

inline r_sym class_sym = r_sym(R_ClassSymbol);
inline r_sym names_sym = r_sym(R_NamesSymbol);
inline r_sym dim_sym = r_sym(R_DimSymbol);
inline r_sym dim_names_sym = r_sym(R_DimNamesSymbol);
inline r_sym row_names_sym = r_sym(R_RowNamesSymbol);
inline r_sym levels_sym = r_sym(R_LevelsSymbol);
inline r_sym double_colon_sym = r_sym(R_DoubleColonSymbol);
inline r_sym triple_colon_sym = r_sym(R_TripleColonSymbol);
inline r_sym dollar_sym = r_sym(R_DollarSymbol);
inline r_sym bracket_sym = r_sym(R_BracketSymbol);
inline r_sym double_brackets_sym = r_sym(R_Bracket2Symbol);
inline r_sym brace_sym = r_sym(R_BraceSymbol);
inline r_sym dots_sym = r_sym(R_DotsSymbol);
inline r_sym tsp_sym = r_sym(R_TspSymbol);
inline r_sym name_sym = r_sym(R_NameSymbol);
inline r_sym base_sym = r_sym(R_BaseSymbol);
inline r_sym quote_sym = r_sym(R_QuoteSymbol);
inline r_sym function_sym = r_sym(R_FunctionSymbol);
inline r_sym namespace_env_sym = r_sym(R_NamespaceEnvSymbol);
inline r_sym package_sym = r_sym(R_PackageSymbol);
inline r_sym seeds_sym = r_sym(R_SeedsSymbol);
inline r_sym na_rm_sym = r_sym(R_NaRmSymbol);
inline r_sym source_sym = r_sym(R_SourceSymbol);
inline r_sym mode_sym = r_sym(R_ModeSymbol);
inline r_sym device_sym = r_sym(R_DeviceSymbol);
inline r_sym last_value_sym = r_sym(R_LastvalueSymbol);
inline r_sym spec_sym = r_sym(R_SpecSymbol);
inline r_sym previous_sym = r_sym(R_PreviousSymbol);
inline r_sym sort_list_sym = r_sym(R_SortListSymbol);
inline r_sym eval_sym = r_sym(R_EvalSymbol);
inline r_sym drop_sym = r_sym(R_DropSymbol);
inline r_sym missing_arg = r_sym(R_MissingArg);
inline r_sym unbound_value = r_sym(R_UnboundValue);

inline r_sym tag(SEXP x){
  return static_cast<r_sym>(TAG(x));
}

}

}

#endif
