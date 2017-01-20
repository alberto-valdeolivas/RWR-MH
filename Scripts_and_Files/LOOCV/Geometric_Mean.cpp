#include <Rcpp.h>
#include <iostream>     
#include <algorithm>  
#include <vector>  

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector Geometric_Mean(NumericVector Scores, int L, int N) {
  
  Rcpp::NumericVector FinalScores(N,1.0);
  
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < L; j++) {
      FinalScores[i] = FinalScores[i]*Scores[i+ N*j] ;
      
    }
    FinalScores[i] = pow(FinalScores[i], 1.0/L) ;
  }
  return(FinalScores);
}
