#include <RcppArmadillo.h>
#include <math.h>  
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp ;

// [[Rcpp::export]]
arma::sp_mat get_transition_multiplex_rcpp(int Number_Proteins, int Number_Layers, double lambda, arma::sp_mat SupraAdjacencyMatrix, arma::sp_mat SupraBipartiteMatrix)  {
  
  int nrow = Number_Proteins*Number_Layers, ncol =Number_Proteins*Number_Layers;
  arma::sp_mat Transition_Multiplex_Network(nrow, ncol);
  
  arma::sp_mat Col_Sum_Multiplex = sum(SupraAdjacencyMatrix,0) ;
  arma::sp_mat Row_Sum_Bipartite = sum(SupraBipartiteMatrix,1) ;

//  for(int i = 0; i < 10; i++){
    for (int j = 0; j < ncol; j++){
//      printf ("j: %d ", j);
      if (Row_Sum_Bipartite(j) != 0) { 
        Transition_Multiplex_Network.col(j) = ((1-lambda)*SupraAdjacencyMatrix.col(j))/Col_Sum_Multiplex(j);
      } else {
        Transition_Multiplex_Network.col(j) = SupraAdjacencyMatrix.col(j) /Col_Sum_Multiplex(j);
      }
    }
//  }
  
  return Transition_Multiplex_Network ;
}

