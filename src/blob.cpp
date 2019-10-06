#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace std;

//===============================================
// Utility functions:
//===============================================
//
//
//============================================================
// Two state spinor function for internal use.
// We Don't want to pass in R functions so use this C version
// x: 1 or 2 ONLY
//============================================================
// [[Rcpp::export]]
arma::Mat<int> ff_C(int x){
  arma::Mat<int> aspinor(2,1);

  aspinor[0] = (x == 1);
  aspinor[1] = (x == 2);

  return(aspinor);

}

//===================================================================================
// CRF default node_par and edge_par are 3D "Cubes" with 1 slice (3 index R arrays).
// This function makes them armadillo integer matrices (2D) to head of annoying type 
// conflicts and code bloat from re-use. For now, we export to R as well for testing.
// Ultimately intended for internal use.
// Armadillo version
//===================================================================================
// [[Rcpp::export]]
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par){
  
  arma::Mat<int> node_par_new(node_par.n_rows,2);
  arma::Mat<int> aedge_par_mat(2,2);

  node_par_new = node_par.slice(0);

  // Remove third index from each of the edge_par
  IntegerVector tmpe(4);
  List edge_par_new(edge_par.size());
  for(int i=0; i<edge_par.size(); ++i) {

    tmpe = as<IntegerVector>(edge_par(i));
    aedge_par_mat(0,0) = tmpe(0);
    aedge_par_mat(0,1) = tmpe(1);
    aedge_par_mat(1,0) = tmpe(2);
    aedge_par_mat(1,1) = tmpe(3);

    edge_par_new(i) = aedge_par_mat;

  }

  List theta_pars;
  theta_pars["node_par"] = node_par_new;
  theta_pars["edge_par"] = edge_par_new;

  return theta_pars;

}

//===============================================
// row.match in C. Much faster than in R and we
// don't want to pass in R functions
//===============================================
// [[Rcpp::export]]
arma::uvec row_match(arma::Mat<int> x, arma::Mat<int> table){
  
  arma::uvec matching_row_idxs;

  // Note: Below assume x is a effectively a vector with one row or one column
  // We pass it in as a matrix though. Noticed Armadillo routenes often ran faster that way......??

  // Get common first set of possible row match indices   
  matching_row_idxs = intersect(find(table.col(0) == x(0)), find(table.col(1) == x(1)));

  // Look over the rest and update (pairdown) the possible matching rows 
  for(int i=2; i<x.size(); ++i){
    matching_row_idxs = intersect(matching_row_idxs, find(table.col(i) == x(i)));
  }

  // Note: a size 0 return means no matching row was found
  return matching_row_idxs;
}


//=====================================================
// features_util.R function ports:
//=====================================================
//
//
//===============================================
// phi.features port to C
//===============================================
//[[Rcpp::export]]
arma::Mat<int> phi_features_C(arma::Mat<int> config, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default=0) {
  
  int num_nodes = config.size();
  int num_edges = edge_mat.n_rows;
  
  // The original R function determines this everytime it is called. This is stupid and will
  // slow things down. To keep parity with the R function signature, by default the number of
  // parameters will be calculated. However the user can also prespcify it and change the default
  // value from 0 so the loop below isn't excecuted over an over when calling the function
  // multiple times.
  int num_params;
  if(num_params_default == 0) {
    
    num_params = node_par.max();
    int amax       = 0;
    for(int i=0; i<edge_par.size(); ++i) {
      amax = as<arma::Mat<int>>( edge_par(i) ).max(); // a copy performed with this as<>() ????
      if(amax > num_params){
        num_params = amax;
      } else {
        //num_params = num_params_default;
        stop("num_pars specification broken...");
      }
    }
    
  } 
  
  // Initialize a phi (row) vector. Choose row vector. Thats what the R code in compute.model.matrix assumes.
  // We do this here to keep parity with the R code.
  arma::Mat<int> phi_vec(1,num_params);
  phi_vec.zeros();
  
  // Nodes: \phi_i({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \tau}_i, 0}
  int phi_off = -1;                                        // -1 bec we want offsets NOT indices
  for(int i=0; i<num_nodes; ++i) {
    
    // Usually reversed but these are rowvec . colvec. Doing it this way avoids a transpose 
    phi_off = dot( node_par.row(i), ff_C(config(i)) ) - 1; // -1 bec we want offsets NOT indices
    
    // Kronecker delta part:
    if(phi_off != -1) {
      phi_vec(phi_off) = 1; 
    }
    
  }
  
  // Edges: \phi_{k_{[ij]}}({\bf X}) = 1-\delta_{{\bf f}^{\dagger}(X_i) {\boldsymbol \omega}_{ij} {\bf f}(X_j) , 0}
  arma::Mat<int> aepm;
  for(int i=0; i<num_edges; ++i) {
    
    // Check and see if we reached the end of phi. No point in doing the rest of the edges if we did:
    // NOTE: assumes parameters have unique consecutive indices, so no wierd parameterizations.
    // Sticking to "standard" or "flexible" should be safe.
    if(phi_off == num_params) {
      break;
    }
    
    aepm = as<arma::Mat<int>>(edge_par(i));
    int left_off  = edge_mat(i,0) - 1;                                    // offset NOT index, do -1
    int right_off = edge_mat(i,1) - 1;                                    // offset NOT index, do -1
    
    arma::Mat<int> tmp = ff_C(config(left_off)).t() * aepm * ff_C(config(right_off));
    
    phi_off = tmp(0,0) - 1;                                               // offset NOT index, do -1
    
    // Kronecker delta part:
    if(phi_off != -1) {
      phi_vec(phi_off) = 1; 
    }
    
  }
  
  return phi_vec;
  
}

//===============================================
// compute.model.matrix port
//===============================================
// [[Rcpp::export]]
arma::Mat<int> compute_model_matrix(arma::Mat<int> configs, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default = 0) {
  
  int num_configs = configs.n_rows;
  int num_params;
  
  if(num_params_default == 0) {
    
    num_params = node_par.max();
    int amax       = 0;
    for(int i=0; i<edge_par.size(); ++i) {      // ********* RcppParallel HERE???????
      amax = as<arma::Mat<int>>( edge_par(i) ).max(); // a copy performed with this as<>() ????
      if(amax > num_params){
        num_params = amax;
      } else {
        //num_params = num_params_default;
        stop("num_pars specification broken...");
      }
    }
    
  }
  
  arma::Mat<int> model_mat(num_configs, num_params);
  
  for(int i=0; i<num_configs; ++i){
    model_mat.row(i) = phi_features_C(configs.row(i), edge_mat, node_par, edge_par);
  }
  
  return(model_mat);
  
  
}

//===============================================
// get.par.idx port Now an offset however  ******* R-Nullables may be causing PROBLEMS?
//===============================================
// [[Rcpp::export]]
int get_par_off(arma::Mat<int>                config, 
                Rcpp::Nullable<int>           i_in        = R_NilValue, 
                Rcpp::Nullable<int>           j_in        = R_NilValue, 
                Rcpp::Nullable<IntegerMatrix> node_par_in = R_NilValue,
                Rcpp::Nullable<List>          edge_par_in = R_NilValue,
                Rcpp::Nullable<IntegerMatrix> edge_mat_in = R_NilValue,
                bool                          printQ      = false) {
  
  int i, j;
  List edge_par;
  arma::Mat<int> edge_mat;
  arma::Mat<int> node_par;
  int par_off = -1;          // parameter offset NOT index
  
  if(i_in.isNotNull()) {
    i = as<int>(i_in);
    if(j_in.isNotNull()) {
      
      // An edge was input
      j = as<int>(j_in);
      
      // Need edge_par
      if(edge_par_in.isNotNull()) {
        edge_par = edge_par_in;
        // for(int i=0; i<edge_par.size(); ++i){
        //   Rcout << as<arma::Mat<int>>(edge_par(i)) << endl;
        // }
      } else {
        stop("Edge param queried but no edge par input.");
      }
      
      // Need edge mat
      if(edge_mat_in.isNotNull()) {
        edge_mat = as<arma::Mat<int>>(edge_mat_in);
        //Rcout << edge_mat << endl;
      } else {
        stop("Edge param queried but no edge mat input.");
      }
      
      // Check and see if the edge indices are together in the edge matrix
      // Note: actually get the edge offset, not index
      arma::Mat<int> avec(1,2);
      avec(0,0) = i;
      avec(0,1) = j;
      
      // edge offset in the edge matrix:
      arma::uvec edge_off = row_match(avec, edge_mat);
      
      if(edge_off.size() == 0) {
        Rcout << "Input edge indices: i=" << i  << " j=" << j << endl;
        stop("Input edge indices not found in edge mat");
      }
      if(edge_off.size() > 1) {
        Rcout << "Input edge indices: i=" << i << " j=" << j << endl;
        stop("Something is wierd. Mutiple instances of this edge found in edge mat.");
      }
      
      // If all looks ok, compute parameter offset (not index!) associated with edge
      arma::Mat<int> aepm;
      aepm = as<arma::Mat<int>>(edge_par(edge_off(0)));
      int left_off  = edge_mat(edge_off(0),0) - 1;                                    // offset NOT index, do -1
      int right_off = edge_mat(edge_off(0),1) - 1;                                    // offset NOT index, do -1
      
      arma::Mat<int> tmp;                                                             // to hold product
      tmp = ff_C(config(left_off)).t() * aepm * ff_C(config(right_off));
      
      par_off = tmp(0,0) - 1;                                                         // offset NOT index, do -1
      
    } else {
      
      //A node was input
      if(node_par_in.isNotNull()) {
        node_par = as<arma::Mat<int>>(node_par_in);
        //Rcout << node_par << endl;
      } else {
        stop("Node param queried but no node par input.");
      }
      
      // Need to get node row (offset) out of node par
      // Note: assumes each node only has one parameter. One parameter can be shared
      // between many nodes however.
      
      // If all looks ok, compute parameter offset (not index!) associated with node
      int node_off = i-1;                                                  // Just for readability for when I forget how I did this...
      par_off = dot( node_par.row(node_off), ff_C(config(node_off)) ) - 1; // -1 bec we want offsets NOT indices
    }
    
  } else {
    stop("No node i index entered!");
  }
  
  // Note: If offset is -1 that means parameter is not associated with this node/edge and spin set
  return par_off;
}


//===============================================
// phi.component port  **** R nullables
//===============================================
// [[Rcpp::export]]
int phi_component(arma::Mat<int>                config,
                  Rcpp::Nullable<int>           i_in        = R_NilValue, 
                  Rcpp::Nullable<int>           j_in        = R_NilValue, 
                  Rcpp::Nullable<IntegerMatrix> node_par_in = R_NilValue,
                  Rcpp::Nullable<List>          edge_par_in = R_NilValue,
                  Rcpp::Nullable<IntegerMatrix> edge_mat_in = R_NilValue) {
  
  int par_off = get_par_off(config, i_in, j_in, node_par_in, edge_par_in, edge_mat_in);
  
  int swtch;
  if(par_off == -1) {
    swtch = 1;
  } else {
    swtch = 0;
  }
  
  int comp = 1 - swtch;
  
  return comp;
  
}


//=====================================================
// energy_util.R function ports:
//=====================================================
//
//
//===============================================
// symbolic.conditional.energy port to C
//symbolic.conditional.energy(config, condition.element.number, crf, ff, format="tex", printQ=FALSE)
//===============================================
//[[Rcpp::export]]
arma::Mat<int> symbolic_conditional_energy(arma::Mat<int> config, int condition_element_number, arma::Mat<int> edge_mat, arma::Mat<int> node_par, List edge_par, int num_params_default=0) {
  
  // Same as for phi_features_C:
  // The original R function determines this everytime it is called. This is stupid and will
  // slow things down. To keep parity with the R function signature, by default the number of
  // parameters will be calculated. However the user can also prespcify it and change the default
  // value from 0 so the loop below isn't excecuted over an over when calling the function
  // multiple times.
  int num_params;
  if(num_params_default == 0) {
    
    num_params = node_par.max();
    int amax       = 0;
    for(int i=0; i<edge_par.size(); ++i) {
      amax = as<arma::Mat<int>>( edge_par(i) ).max(); // a copy performed with this as<>() ????
      if(amax > num_params){
        num_params = amax;
      } else {
        //num_params = num_params_default;
        stop("num_pars specification broken...");
      }
    }
    
  } 
  
  // Initialize a out_eq (row) vector. Choose row vector. 
  arma::Mat<int> out_eq(1,num_params);
  out_eq.zeros();
  

  return out_eq;
}


