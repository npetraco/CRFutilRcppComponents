#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace std;


//===============================================
// Tests and prototypeing
//===============================================
//
//
// [[Rcpp::export]]
void testemptymat(arma::Mat<int> emptymat) {
  
  Rcout << emptymat(0,0) << endl;
  
}