# BayesUnmixing
Bayesian hierarchical model for estimating (unmixing) diet composition using a Gibbs sampler. Uses a methodology similar to the one 
describedby Yu (2015). Its intended use is for delineation (estimation) of diet composition in ruminants through plant-wax markers, 
alkanes or long chain alohols. 
The posterior means for the forages proportions are sampled from a multivariate normal distribution truncated on a simplex using the
methodology described in Dobigeon et al, (2007). Posterior means for the covariance matrix and its scale matrix are sampled from an 
inverse Wishart and a Wishart distribution respectively.
Written in C++ through Rcpp and RcppArmadillo.
