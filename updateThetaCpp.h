#include <RcppArmadillo.h>
#include <cmath>
#include <rmtnormsimplex.h>
#include <trnorm.h>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]] 
arma::vec updateTheta(arma::vec ThetaStar, arma::mat M, arma::mat Sigma, arma::vec y, arma::mat Psi, arma::uvec estim, arma::uvec last, 
                      arma::mat Mstar, arma::vec mp, arma::vec One){
    uword p                 = M.n_cols;
    vec theta_old(1);
    vec tmp1 = zeros(1);
    vec Mean                = zeros(p - 1);
    vec Mu                  = zeros(p - 1);
    vec Theta_p;
    if (p == 2){
        theta_old               = ThetaStar.elem(estim);
        Mean                    = (1/Psi[0,0])*trans(Mstar - mp*trans(One))*solve(Sigma, y - mp);
        //Mu                      = rtnorm(0.0, 1.0, Mean[0], sqrt(1/Psi[0, 0]));
        Mu                      = dtrandn(theta_old[0], Mean[0], sqrt(1/Psi[0, 0]), 0.0, 1.0);
        ThetaStar.elem(estim)   = Mu;
        ThetaStar.elem(last)    = 1 - Mu;
    }
    else {
        Mean                    = solve(Psi, eye(p - 1, p - 1))*trans(Mstar - mp*trans(One))*solve(Sigma, y - mp);
        Mu                      = rmtnorm_simplex(ThetaStar.elem(estim), Mean, solve(Psi, eye(p - 1, p - 1)));
        Theta_p << 1 - accu(Mu);
        ThetaStar.elem(estim)   = Mu;
        ThetaStar.elem(last)    = max(Theta_p, tmp1);
    }

    return (ThetaStar);
}