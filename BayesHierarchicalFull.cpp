#include <RcppArmadillo.h>
#include <cmath>
#include <time.h>
#include <updateThetaCpp.h>
#include <mysetdiff.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

mat rwishart(int const& nu, mat const& V){

  int m = V.n_rows;
  mat T = zeros(m,m);
  mat R = zeros(m,m);
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]); 
  }
   
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = rnorm(1)[0]; 
    }}
  
  mat C = trans(T)*chol(V);
  
    return R = trans(C)*C;
}

mat rinvwishart(int const& nu, mat const& V){

  int m = V.n_rows;
  mat T = zeros(m,m);
  mat IR = zeros(m,m);
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(rchisq(1,nu-i)[0]);
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {  
      T(i,j) = rnorm(1)[0]; 
    }}
  
  mat C = trans(T)*chol(V);
  mat CI = solve(trimatu(C),eye(m,m)); 

    return IR = CI*trans(CI);
}

// [[Rcpp::export]] 
List BayesUnmixing(mat const& M, mat const& Y, double const& nu, mat const& Omega, double const& gamma, uword const& nIter){
  Rcpp::Rcout << " ****************************************************** " << std::endl;
  Rcpp::Rcout << " Bayesian hierarchical model for unmixing compositions  " << std::endl;
  Rcpp::Rcout << " ****************************************************** " << std::endl;
  Rcpp::Rcout << "                                                        " << std::endl;
  uword r             = M.n_rows;
  uword p             = M.n_cols;
  uword nAnim         = Y.n_cols;
  Rcpp::Rcout << " Number of plant-wax markers: " << r << std::endl;
  Rcpp::Rcout << " Number of forages/plants " << p <<std::endl;
  Rcpp::Rcout << " Number of observations " << nAnim << std::endl;
  Rcpp::Rcout << " Starting Gibbs sampler... " << std::endl;
  uvec plants;
  uvec tmp;
  uvec lastj;
  uvec inTheta;
  uvec estim;
  uvec last;
  vec Theta_p;
  mat Mstar;
  mat mp;
  mat Psi;
  vec tmp1          = zeros(1);
  mat One                 = ones <mat> (p - 1, 1);
  plants                  = sequence(1, p);
  double ttime      = 0;
  //double pct        = 0;
  clock_t start;
  clock_t end;
  vec init(p);
  mat Theta         = ones(p, nAnim)/p;
  /*for (uword t = 0; t < nAnim; t++){
   init          = runif(p);
   Theta.col(t) = init/accu(init);
  }*/
  mat store_Sigma   = zeros(nIter, r);
  mat store_Xi      = zeros(nIter, r);
  cube store_Theta  = zeros(nIter, p, nAnim);
  List Out;
  
  /* For updating Theta_i */
  vec ThetaStar;
  
  /* For updating Sigma */
  mat Sigma = rinvwishart(nAnim + nu, Omega);   
  mat Xi    = rwishart(gamma, Omega);
  mat R     = zeros <mat> (r, r);
  
  start = clock();
  
  for(uword iter = 0; iter < nIter; iter++){
    tmp                     = shuffle(plants);
    lastj                   = tmp(0);
    inTheta                 = mysetdiff(tmp, lastj);
    estim                   = inTheta -  1;
    last                    = lastj - 1;
    Mstar                   = M.cols(estim);
    mp                      = M.cols(last); 
    Psi                     = trans(Mstar - mp*trans(One))*solve(Sigma, Mstar - mp*trans(One));
    for(uword i = 0; i < nAnim; i++){
      ThetaStar                                   = updateTheta(Theta.col(i), M, Sigma, Y.col(i), Psi, estim, last, Mstar, mp, One);
      store_Theta(span(iter), span::all, span(i)) = ThetaStar;
      Theta.col(i)                                = ThetaStar;
      R                                           = R  + (Y.col(i) - M*ThetaStar)*trans(Y.col(i) - M*ThetaStar);
      }

      Sigma                 = rinvwishart(nAnim + nu, solve(R + Xi, eye(r, r)));
      Xi                    = rwishart(nu + gamma, Omega + Sigma);
      store_Sigma.row(iter) = trans(Sigma.diag()); 
      store_Xi.row(iter)    = trans(Xi.diag()); 
      R                     = zeros(r, r);
      if(iter % 1000 == 0){
        //pct = iter/nIter;
        Rcpp::Rcout << " Iteration " << iter << "/" << nIter << std::endl;//" (" << (pct*100.00) << "%) "<< std::endl;
      }  
    }
    
    end = clock();
    ttime = ((double) (end - start)) / CLOCKS_PER_SEC;
    Rcpp::Rcout << "                                               " << std::endl;
    Rcpp::Rcout << " Wrapping up! " << std::endl;
    Rcpp::Rcout << "                                               " << std::endl;
    Rcpp::Rcout << nIter << " iterations in " << ttime << " seconds" << std::endl;
  
    Out["Theta"]        = store_Theta;
    Out["Sigma"]        = store_Sigma;
    Out["Xi"]           = store_Xi;
    Out["Elapsed time"] = ttime;
    return(wrap(Out));    
}