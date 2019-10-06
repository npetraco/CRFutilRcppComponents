// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ff_C
arma::Mat<int> ff_C(int x);
RcppExport SEXP _CRFutilRcppComponents_ff_C(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(ff_C(x));
    return rcpp_result_gen;
END_RCPP
}
// fix_node_and_edge_par
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par);
RcppExport SEXP _CRFutilRcppComponents_fix_node_and_edge_par(SEXP node_parSEXP, SEXP edge_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Cube<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    rcpp_result_gen = Rcpp::wrap(fix_node_and_edge_par(node_par, edge_par));
    return rcpp_result_gen;
END_RCPP
}
// row_match
arma::uvec row_match(arma::Mat<int> x, arma::Mat<int> table);
RcppExport SEXP _CRFutilRcppComponents_row_match(SEXP xSEXP, SEXP tableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type table(tableSEXP);
    rcpp_result_gen = Rcpp::wrap(row_match(x, table));
    return rcpp_result_gen;
END_RCPP
}
// phi_features_C
arma::Mat<int> phi_features_C(arma::Mat<int> config, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default);
RcppExport SEXP _CRFutilRcppComponents_phi_features_C(SEXP configSEXP, SEXP edge_matSEXP, SEXP node_parSEXP, SEXP edge_parSEXP, SEXP num_params_defaultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type config(configSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    Rcpp::traits::input_parameter< int >::type num_params_default(num_params_defaultSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_features_C(config, edge_mat, node_par, edge_par, num_params_default));
    return rcpp_result_gen;
END_RCPP
}
// compute_model_matrix
arma::Mat<int> compute_model_matrix(arma::Mat<int> configs, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default);
RcppExport SEXP _CRFutilRcppComponents_compute_model_matrix(SEXP configsSEXP, SEXP edge_matSEXP, SEXP node_parSEXP, SEXP edge_parSEXP, SEXP num_params_defaultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type configs(configsSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    Rcpp::traits::input_parameter< int >::type num_params_default(num_params_defaultSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_model_matrix(configs, edge_mat, node_par, edge_par, num_params_default));
    return rcpp_result_gen;
END_RCPP
}
// get_par_off
int get_par_off(arma::Mat<int> config, Rcpp::Nullable<int> i_in, Rcpp::Nullable<int> j_in, Rcpp::Nullable<IntegerMatrix> node_par_in, Rcpp::Nullable<List> edge_par_in, Rcpp::Nullable<IntegerMatrix> edge_mat_in, bool printQ);
RcppExport SEXP _CRFutilRcppComponents_get_par_off(SEXP configSEXP, SEXP i_inSEXP, SEXP j_inSEXP, SEXP node_par_inSEXP, SEXP edge_par_inSEXP, SEXP edge_mat_inSEXP, SEXP printQSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type i_in(i_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type j_in(j_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<IntegerMatrix> >::type node_par_in(node_par_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<List> >::type edge_par_in(edge_par_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<IntegerMatrix> >::type edge_mat_in(edge_mat_inSEXP);
    Rcpp::traits::input_parameter< bool >::type printQ(printQSEXP);
    rcpp_result_gen = Rcpp::wrap(get_par_off(config, i_in, j_in, node_par_in, edge_par_in, edge_mat_in, printQ));
    return rcpp_result_gen;
END_RCPP
}
// phi_component
int phi_component(arma::Mat<int> config, Rcpp::Nullable<int> i_in, Rcpp::Nullable<int> j_in, Rcpp::Nullable<IntegerMatrix> node_par_in, Rcpp::Nullable<List> edge_par_in, Rcpp::Nullable<IntegerMatrix> edge_mat_in);
RcppExport SEXP _CRFutilRcppComponents_phi_component(SEXP configSEXP, SEXP i_inSEXP, SEXP j_inSEXP, SEXP node_par_inSEXP, SEXP edge_par_inSEXP, SEXP edge_mat_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type i_in(i_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type j_in(j_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<IntegerMatrix> >::type node_par_in(node_par_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<List> >::type edge_par_in(edge_par_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<IntegerMatrix> >::type edge_mat_in(edge_mat_inSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_component(config, i_in, j_in, node_par_in, edge_par_in, edge_mat_in));
    return rcpp_result_gen;
END_RCPP
}
// symbolic_conditional_energy
arma::Mat<int> symbolic_conditional_energy(arma::Mat<int> config, int condition_element_number, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default);
RcppExport SEXP _CRFutilRcppComponents_symbolic_conditional_energy(SEXP configSEXP, SEXP condition_element_numberSEXP, SEXP edge_matSEXP, SEXP node_parSEXP, SEXP edge_parSEXP, SEXP num_params_defaultSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::Mat<int> >::type config(configSEXP);
    Rcpp::traits::input_parameter< int >::type condition_element_number(condition_element_numberSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type edge_mat(edge_matSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int> >::type node_par(node_parSEXP);
    Rcpp::traits::input_parameter< List >::type edge_par(edge_parSEXP);
    Rcpp::traits::input_parameter< int >::type num_params_default(num_params_defaultSEXP);
    rcpp_result_gen = Rcpp::wrap(symbolic_conditional_energy(config, condition_element_number, edge_mat, node_par, edge_par, num_params_default));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CRFutilRcppComponents_ff_C", (DL_FUNC) &_CRFutilRcppComponents_ff_C, 1},
    {"_CRFutilRcppComponents_fix_node_and_edge_par", (DL_FUNC) &_CRFutilRcppComponents_fix_node_and_edge_par, 2},
    {"_CRFutilRcppComponents_row_match", (DL_FUNC) &_CRFutilRcppComponents_row_match, 2},
    {"_CRFutilRcppComponents_phi_features_C", (DL_FUNC) &_CRFutilRcppComponents_phi_features_C, 5},
    {"_CRFutilRcppComponents_compute_model_matrix", (DL_FUNC) &_CRFutilRcppComponents_compute_model_matrix, 5},
    {"_CRFutilRcppComponents_get_par_off", (DL_FUNC) &_CRFutilRcppComponents_get_par_off, 7},
    {"_CRFutilRcppComponents_phi_component", (DL_FUNC) &_CRFutilRcppComponents_phi_component, 6},
    {"_CRFutilRcppComponents_symbolic_conditional_energy", (DL_FUNC) &_CRFutilRcppComponents_symbolic_conditional_energy, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_CRFutilRcppComponents(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
