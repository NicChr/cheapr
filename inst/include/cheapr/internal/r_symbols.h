#ifndef CHEAPR_R_SYMBOLS_H
#define CHEAPR_R_SYMBOLS_H

#include <r_types.h>

namespace cheapr {

namespace symbol {

inline r_symbol_t class_sym = r_symbol_t(R_ClassSymbol);
inline r_symbol_t names_sym = r_symbol_t(R_NamesSymbol);
inline r_symbol_t dim_sym = r_symbol_t(R_DimSymbol);
inline r_symbol_t dim_names_sym = r_symbol_t(R_DimNamesSymbol);
inline r_symbol_t row_names_sym = r_symbol_t(R_RowNamesSymbol);
inline r_symbol_t levels_sym = r_symbol_t(R_LevelsSymbol);
inline r_symbol_t double_colon_sym = r_symbol_t(R_DoubleColonSymbol);
inline r_symbol_t triple_colon_sym = r_symbol_t(R_TripleColonSymbol);
inline r_symbol_t dollar_sym = r_symbol_t(R_DollarSymbol);
inline r_symbol_t bracket_sym = r_symbol_t(R_BracketSymbol);
inline r_symbol_t double_brackets_sym = r_symbol_t(R_Bracket2Symbol);
inline r_symbol_t brace_sym = r_symbol_t(R_BraceSymbol);
inline r_symbol_t dots_sym = r_symbol_t(R_DotsSymbol);
inline r_symbol_t tsp_sym = r_symbol_t(R_TspSymbol);
inline r_symbol_t name_sym = r_symbol_t(R_NameSymbol);
inline r_symbol_t base_sym = r_symbol_t(R_BaseSymbol);
inline r_symbol_t quote_sym = r_symbol_t(R_QuoteSymbol);
inline r_symbol_t function_sym = r_symbol_t(R_FunctionSymbol);
inline r_symbol_t namespace_env_sym = r_symbol_t(R_NamespaceEnvSymbol);
inline r_symbol_t package_sym = r_symbol_t(R_PackageSymbol);
inline r_symbol_t seeds_sym = r_symbol_t(R_SeedsSymbol);
inline r_symbol_t na_rm_sym = r_symbol_t(R_NaRmSymbol);
inline r_symbol_t source_sym = r_symbol_t(R_SourceSymbol);
inline r_symbol_t mode_sym = r_symbol_t(R_ModeSymbol);
inline r_symbol_t device_sym = r_symbol_t(R_DeviceSymbol);
inline r_symbol_t last_value_sym = r_symbol_t(R_LastvalueSymbol);
inline r_symbol_t spec_sym = r_symbol_t(R_SpecSymbol);
inline r_symbol_t previous_sym = r_symbol_t(R_PreviousSymbol);
inline r_symbol_t sort_list_sym = r_symbol_t(R_SortListSymbol);
inline r_symbol_t eval_sym = r_symbol_t(R_EvalSymbol);
inline r_symbol_t drop_sym = r_symbol_t(R_DropSymbol);
inline r_symbol_t missing_arg = r_symbol_t(R_MissingArg);
inline r_symbol_t unbound_value = r_symbol_t(R_UnboundValue);

}

}

#endif
