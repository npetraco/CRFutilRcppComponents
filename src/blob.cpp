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
arma::Mat<int> ff_C(int x) {
  
  arma::Mat<int> aspinor(2,1);

  aspinor[0] = (x == 1);
  aspinor[1] = (x == 2);

  return(aspinor);

}

//===================================================================================
// CRF default node_par and edge_par are 3D "Cubes" with 1 slice (3 index R arrays).
// This function makes them armadillo integer matrices (2D) to head of annoying type 
// conflicts and code bloat from re-use. For now, we export to R as well for testing.
// Ultimately intended for internal use only.
//===================================================================================
// [[Rcpp::export]]
List fix_node_and_edge_par(arma::Cube<int> node_par, List edge_par) {
  
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

  // Node pars and Edge pars with third dimension removed:
  List theta_pars;
  theta_pars["node_par"] = node_par_new;
  theta_pars["edge_par"] = edge_par_new;

  return theta_pars;

}


//===============================================
// row.match in C. Much faster than in R and we
// don't want to pass in R functions anyway
//===============================================
// [[Rcpp::export]]
arma::uvec row_match(arma::Mat<int> x, arma::Mat<int> table) {
  
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


//===============================================
// Function to extract number of parameters if not
// input. Intended for internal use only, but 
// exported for testing
//===============================================
//[[Rcpp::export]]
int get_num_params(arma::Mat<int> node_par, List edge_par) {
  
  int num_params;
  
  num_params = node_par.max();
  int amax   = 0;
  for(int i=0; i<edge_par.size(); ++i) {
    amax = as<arma::Mat<int>>( edge_par(i) ).max(); // a copy performed with this as<>() ????
    if(amax > num_params){
      num_params = amax;
    } else {
      stop("num_pars specification broken...");
    }
  }
  
  return num_params;
}


//=====================================================
// features_util.R function ports:
//=====================================================
//
//
//===============================================
// phi.features port to C  ******* CLEAN THIS UP!!!!!!!
//===============================================
//[[Rcpp::export]]
arma::Mat<int> phi_features_C(arma::Mat<int> config, 
                              arma::Mat<int> edge_mat, 
                              arma::Mat<int> node_par, 
                              List           edge_par, 
                              int            num_params_in=0) 
{
  int num_nodes = config.size();
  int num_edges = edge_mat.n_rows;
  
  // The original R function determines the number of parameters everytime it is called. 
  // This is stupid and will slow things down. However, to keep parity with the R function signature, 
  // by default the number of parameters will be calculated. The user can also pre-specify the
  // number of parameters so that the extra loop the calculation requires is not executed. 
  // The default value of 0 for num_params_in indicates the number of parameters will be computed.
  int num_params;
  if(num_params_in == 0) {
    num_params = get_num_params(node_par, edge_par);
  } else {
    num_params = num_params_in;
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
    // Sticking to "standard" or "flexible" should be safe. ****** REMOVE AT SOME POINT!!!!!!!!!!!
    if(phi_off == num_params) {
      break;
    }
    
    aepm = as<arma::Mat<int>>(edge_par(i));
    int left_off  = edge_mat(i,0) - 1; // offset NOT index, do -1
    int right_off = edge_mat(i,1) - 1; // offset NOT index, do -1
    
    arma::Mat<int> tmp = ff_C(config(left_off)).t() * aepm * ff_C(config(right_off));
    
    phi_off = tmp(0,0) - 1;            // offset NOT index, do -1
    
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
arma::Mat<int> compute_model_matrix(arma::Mat<int> configs, 
                                    arma::Mat<int> edge_mat, 
                                    arma::Mat<int> node_par, 
                                    List           edge_par, 
                                    int            num_params_in = 0) 
{
  int num_configs = configs.n_rows;
  int num_params;
  
  // Check and see if number of parameters was input. If not compute them:
  if(num_params_in == 0) {
    num_params = get_num_params(node_par, edge_par);
  } else {
    num_params = num_params_in;
  }
  
  // Compute model matrix:  
  arma::Mat<int> model_mat(num_configs, num_params);
  
  for(int i=0; i<num_configs; ++i){  // ********* RcppParallel HERE???????
    model_mat.row(i) = phi_features_C(configs.row(i), edge_mat, node_par, edge_par);
  }
  
  return(model_mat);
}


//===============================================
// get.par.idx port Now an offset AND ALMOST NO 
// MORE R-Nullables and almost NO DEFAULT ARGUEMENTS
// * We can't do default arma arguements since they
// are not supported in Rcpp
// * So instead pass in an array object that can be
// cast to an arma Mat and has a -1 as the (0,0) 
// element.
// * i_in and j_in are indicated as empty by passing
// in -1
//===============================================
// [[Rcpp::export]]
int get_par_off(arma::Mat<int>       config, 
                int                  i_in,  
                int                  j_in, 
                arma::Mat<int>       node_par_in,
                Rcpp::Nullable<List> edge_par_in,
                arma::Mat<int>       edge_mat_in,
                bool                 printQ = false) 
{
  int i, j;
  List edge_par;
  arma::Mat<int> edge_mat;
  arma::Mat<int> node_par;
  int par_off = -1;          // parameter offset NOT index
  
  if(i_in != -1) {
    i = i_in;
    if(j_in != -1) {
      
      // A non -1 i_in and j_in passed in indicate and edge was input
      j = j_in;
      
      // Need edge_par for an edge
      if(edge_par_in.isNotNull()) {
        edge_par = edge_par_in;
      } else {
        stop("Edge param queried but no edge par input.");
      }
      
      // Need edge mat for an edge
      if(edge_mat_in(0,0) != -1) { // -1 at (0,0) means empty arguement
        edge_mat = edge_mat_in;
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
        Rcout << "Input edge indices: i=" << i << " j=" << j << endl;
        stop("Input edge indices not found in edge mat");
      }
      if(edge_off.size() > 1) {
        Rcout << "Input edge indices: i=" << i << " j=" << j << endl;
        stop("Something is wierd. Mutiple instances of this edge found in edge mat.");
      }
      
      // If all looks ok, compute parameter offset (not index!) associated with edge
      arma::Mat<int> aepm;
      aepm = as<arma::Mat<int>>(edge_par(edge_off(0))); // parameter index matrix for the edge
      int left_off  = edge_mat(edge_off(0),0) - 1;      // offset NOT index, so do -1
      int right_off = edge_mat(edge_off(0),1) - 1;      // offset NOT index, so do -1
      
      arma::Mat<int> tmp;                               // to hold product
      tmp = ff_C(config(left_off)).t() * aepm * ff_C(config(right_off));
      
      par_off = tmp(0,0) - 1;                           // offset NOT index, so do -1
      
    } else {
      
      //A node was input
      if(node_par_in(0,0) != -1) { //Indicate empty node_par_in by passing in -1 0,0 element 
        node_par = node_par_in;
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
// phi.component port  **** NOW WITH ALMOST NO R nullables and NO DEFAULT ARGUEMENTS
// See get_par_off2 for what "empty" arguements should be passed in as
//===============================================
// [[Rcpp::export]]
int phi_component(arma::Mat<int>        config,
                   int                  i_in, 
                   int                  j_in, 
                   arma::Mat<int>       node_par_in,
                   Rcpp::Nullable<List> edge_par_in,
                   arma::Mat<int>       edge_mat_in) 
{
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
// Port of guts (ONLY) of symbolic.conditional.energy
// We are really just interested in computing the alpha vector in C.
// Recall an alpha vector is (among other things) the number 
// of times each parameter appears in the expression for a 
// conditional energy E(X_i|X/X_i). 
// The alpha vectors are needed to form the Delta-alpha matrix
//===============================================
//[[Rcpp::export]]
arma::Mat<int> alpha_vector(arma::Mat<int> config, 
                            int            condition_element_number, 
                            arma::Mat<int> edge_mat, 
                            arma::Mat<int> node_par, 
                            List           edge_par, 
                            List           adj_nodes, 
                            int            num_params_in=0) 
{
  
  // Check and see if number of parameters was input. If not compute them: 
  int num_params;
  if(num_params_in == 0) {
    num_params = get_num_params(node_par, edge_par);
  } else {
    num_params = num_params_in;
  }
  
  // To indicate an empty arma matrix arguement: ***** IS THERE A BETTER WAY????  
  arma::Mat<int> emptym(1,1);
  emptym(0,0) = -1;
  
  // Initialize an alpha vector (param.num.vec). Choose row vector.
  // This implementation is a little different than the in the R code. There param.num.vec acumulated by appending
  // Here we will allocate the full length of the vector instead and increment the elements as needed.
  // Note: an alpha vector contains the number of times each theta_k appears in E(X_i|X/X_i)(amon other things)
  arma::Mat<int> param_num_vec(1,num_params); // This will be the alpha vector
  param_num_vec.zeros();
  
  // Parameter (if any) associated with conditioned node
  int l = get_par_off(config, condition_element_number, -1, node_par, R_NilValue, emptym, false);
  if(l >= 0) {
    param_num_vec(l) ++;
  }
  
  // MAKE adj_nodes NULLABLE LATER in case we input a model with unconnected nodes *********
  // adj_nodes is pased in as a List. Didn't bother to make nullable because we probably wouldn't use this function 
  // if a node in the model has no attached neighbors.
  IntegerVector adj_nodes_loc = (IntegerVector)adj_nodes(condition_element_number-1); // -1 for offset conversion
  IntegerVector edge_nods(2);
  
  for(int i=0; i<adj_nodes_loc.size(); ++i) {

    edge_nods(0) = condition_element_number;
    edge_nods(1) = adj_nodes_loc(i);
    std::sort(edge_nods.begin(), edge_nods.end());
    
    int k = get_par_off(config, edge_nods(0), edge_nods(1), emptym, edge_par, edge_mat, false);
    if(k >= 0) {
      param_num_vec(k) ++;
    }

  }

  return param_num_vec;
}


//=====================================================
// logistic_util.R function ports:
//=====================================================
//
//
//===============================================
// Port of delta.alpha
//===============================================
//[[Rcpp::export]]
arma::Mat<int> delta_alpha(arma::Mat<int> samples, 
                           arma::Mat<int> node_par, 
                           List           edge_par, 
                           arma::Mat<int> edge_mat,
                           List           adj_nodes,
                           int            num_params_in=0) 
{
  // Check and see if number of parameters was input. If not compute them:  
  int num_params;
  if(num_params_in == 0) {
    num_params = get_num_params(node_par, edge_par);
  } else {
    num_params = num_params_in;
  }
  
  int num_nodes = samples.n_cols; // For readability
  
  // Initialize some memory to be used and re-used
  arma::Mat<int> Da_mat(samples.n_rows*num_nodes, num_params, arma::fill::zeros);
  arma::Mat<int> X_cfg(1, samples.n_cols,arma::fill::zeros);
  arma::Mat<int> Xc_cfg(1, samples.n_cols,arma::fill::zeros);
  
  // Compute Delta-alpha matrix:
  int count = 0;
  for(int i=0; i<num_nodes; ++i) {
    for(int n=0; n<samples.n_rows; ++n) {
      
      X_cfg = samples.row(n);
      Xc_cfg = samples.row(n);
      X_cfg(0,i) = 1;
      Xc_cfg(0,i) = 2;

      Da_mat.row(count) =  alpha_vector(X_cfg,  i+1, edge_mat, node_par, edge_par, adj_nodes, num_params_in)
                         - alpha_vector(Xc_cfg, i+1, edge_mat, node_par, edge_par, adj_nodes, num_params_in);

      count++;
      
    }
  } 
  
  return Da_mat;
  
}