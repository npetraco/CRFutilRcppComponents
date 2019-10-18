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
  
  // DOESN'T WORK. USE SEXP INSTEAD????
  
  //string s = typeid(amat).name();
  //Rcout << s << endl;
  Rcout << amat.n_slices << endl;
  
}


// [[Rcpp::export]]
void testslicedetect2(SEXP amat) {
  
  // WORKS BUT TYPENAME OUTPUT IS IDENTICAL FOR MATS AND ARRAYS
  string s = typeid(amat).name();
  Rcout << s << endl;
  //Rcout << amat.n_slices << endl;
  
}


// [[Rcpp::export]]
NumericMatrix testslicedetect3(SEXP amat) {
  
  // DOESN'T WORK FOR ARRAYS. 
  
  string s = typeid(amat).name();
  Rcout << s << endl;
  
  return(NumericMatrix(amat));
  
}

// [[Rcpp::export]]
void testslicedetect4(SEXP amat) {
  
  // DOESN'T WORK
  // Rcout << Dimension(amat) << endl;
  
  //string s = typeid(amat).name();
  //Rcout << amat.size() << endl;
  
  // DOESN'T WORK for matrices
  arma::Cube<int> aamat = as<arma::Cube<int>>(amat);
  aamat.print();
  
  
}


// [[Rcpp::export]]
void testslicedetect5(NumericVector amat) {
  
  // DOESNT WORK........

  //string s = typeid(amat).name();
  Rcout << amat.size() << endl; // Works for both mats and arrays but just spits out num elems
  
  //Rcout << amat.dims() << endl; // Doesn't work with and without ()
  Rcpp::Dimension d(amat);
  Rcout << d[1] << endl; // Non informative
  
}

