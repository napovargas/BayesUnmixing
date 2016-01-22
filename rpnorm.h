#include <limits>
#include <cmath>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

double rpnorm(double mean, double dev){
const double  A  = 1.136717791056118;
const double pi = M_PI;
double var = dev*dev;
double v = std::numeric_limits<double>::quiet_NaN();
double mA = (1 - A*A)/A*dev;
double mC = dev*std::sqrt(pi/2.0);
double a = 0;
double z = 0;
double rho = 0;
double r = 0;
double u = 0;
double g = 0;

    while(std::isnan(v)){
        if (mean < mA){
            a = (-mean + sqrt(mean*mean + 4.0*var))/2/var;
            z = -log(1.0 - runif(1)[0])/a;
            rho = exp(-(z - mean)*(z - mean)/2.0/var - a*(mean - z + a*var/2.0));
        }
        
        else if(mean <= 0){
            z = std::abs(rnorm(1)[0])*dev + mean;
			rho = (z >= 0.0)?1.0:0.0;
        }
        
        else if(mean < mC){
            r = (runif(1)[0] < mean/(mean + std::sqrt(pi/2.0)*dev))?1.0:0.0;
			u = runif(1)[0]*mean;
			g = std::abs(rnorm(1)[0]*dev) + mean;
			z = r*u + (1.0 - r)*g;
			rho = r*exp(-(z - mean)*(z - mean)/2.0/var) + (1.0 - r);
        }
        
        else {
            z = rnorm(1)[0]*dev + mean;
			rho = (z >= 0)?1.0:0.0;
        }
        
        if(runif(1)[0] < rho){
            v = z;
        }
    }
    return (v);
}

