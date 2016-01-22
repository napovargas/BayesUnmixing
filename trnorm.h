#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double dtrandn(double x_old, double mu, double sigma, double mum, double mup)
{
  double sigma2;
  int accept;
  int compt;
  double x;
  double z;
  double d0;
  sigma2 = sigma*sigma;
  mup = (mup - mu)/sigma;
  mum = (mum - mu)/sigma;

  accept = 0;
  compt = 0;
  x = x_old;
  while ((accept == 0) && (compt < 200)) {
    compt++;
    z = runif(1)[0]*(mup - mum) + mum;
    if (0.0 < mum) {
      d0 = exp((mum * mum - z * z) / 2.0);
    } else if (mup < 0.0) {
      d0 = exp((mup * mup - z * z) / 2.0);
    } else {
      d0 = exp(-(z * z) / 2.0);
    }

    if (runif(1)[0] < d0) {
      x = z;
      accept = 1;
    }
  }

  return (x * sqrt(sigma2) + mu) * (double)accept + x_old * (1.0 - (double)accept);
}
