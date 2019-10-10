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


// [[Rcpp::export]]
void testslicedetect(arma::Cube<int> amat) {
  
  // doesnt work
  
  //string s = typeid(amat).name();
  //Rcout << s << endl;
  Rcout << amat.n_slices << endl;
  
}