# dplfit.R - Density, distribution function, quantile function, random
# generation, and maximum-likelihood estimation for the self-similar
# discrete power-law distribution
#
# (c) Mitchell Newberry 2019 - licensed under GPL version 3.0, see LICENSE
#
library(ggplot2)

require(VGAM)
# Fix a bug in VGAM by redefining pzeta if VGAM version is 1.0.5 or less.
# Thomas Yee says it will be fixed in the next release.
# https://www.stat.auckland.ac.nz/~yee/VGAM/prerelease/
pzeta = ifelse( 
  compareVersion(as.character(packageVersion("VGAM")), "1.0.5") > 0,
  VGAM::pzeta,
  function(d,shape,...) {
  return(ifelse(d < 11, VGAM::pzeta(d,shape,...),
    ifelse(d == 11, VGAM::pzeta(10,shape,...) + VGAM::dzeta(11,shape),
      VGAM::pzeta(d-1, shape,...)))) })
# Test that pzeta is working, modulo floating-point error:
stopifnot(dzeta(9:13,2) - (pzeta(9:13,2)-pzeta(8:12,2)) < 10^-15)

# We define our own continuous (Pareto) distribution with P(x) \propto x^-a, so
# don't use the one from VGAM.

# There are two conventions for the exponent in the Pareto distribution,
# according to whether corresponds to the log-log slope of the probability
# density function or the cumulative distribution function. The two conventions
# differ by 1.

# Clauset et al. style exponent definition
dpareto = function(x, a=2, xm=1) (a - 1)*xm^(a - 1)*x^-a
ppareto = function(x, a=2, xm=1) (x > xm)*(1-(x/xm)^(1 - a))
qpareto = function(u, a=2, xm=1) xm*(1-u)^(1/(1 - a))
rpareto = function(n, a=2, xm=1) qpareto(runif(n),a,xm) 

# Pareto exponent definition
dpareto=function(x, a=0.5, b=1) a*b^a/x^(a+1) 
ppareto=function(x, a=0.5, b=1) (x > b)*(1-(b/x)^a)
qpareto=function(u, a=0.5, b=1) b/(1-u)^(1/a)
rpareto=function(n, a=0.5, b=1) qpareto(runif(n),a,b)

# Self-similar discrete power law
# distribution, density, quantile and sampler functions
ddpl = function(x, a=2, l=sqrt(2), xm=1) {
  if(any(abs((log(x/xm)/log(l)) - round(log(x/xm)/log(l))) > 5e-15)) {
    warning(sprintf("x = %f is not in the support of dpl", x)) }
  return(ifelse(
    abs((log(x/xm)/log(l)) - round(log(x/xm)/log(l))) > 5e-15,
    0, ((1 - l^(-a))*(xm/x)^a))) }
pdpl = function(x, a=2, l=sqrt(2), xm=1) 
  (1 - l^(-a*(floor(log(x/xm)/log(l)) + 1)))
qdpl = function (u, a=2, l=sqrt(2), xm=1) {
  lma = l^(-a)
  return(xm * l^floor(log(1 - u/xm^a)/(-a*log(l)))) } 
rdpl = function (n, a=2, l=sqrt(2), xm=1) qdpl(runif(n),a,l,xm)

# Maximum-likelihood fits to the self-similar discrete power law

# --- UNSAFE --- internal functions
# functions tagged _ mean unsafe, because they do not do bin the input.
dplfit_ = function(xs,lambda,xm) {
  xs = xs[xs >= xm]
  logl = function(x) { return(logb(x, lambda)) }
  return(logl(1 + length(xs)/sum(logl(xs/xm)))) }

dplfit_se_ = function(xs,lambda,xm) {
  xs = xs[xs >= xm]
  n = length(xs)
  logl = function(x) { return(logb(x, lambda)) }
  alpha = logl(1 + length(xs)/sum(logl(xs/xm)))
  var = (lambda^alpha - 1)^2 / (n * lambda^alpha * log(lambda)^2)
  return(data.frame(alpha = alpha, se = sqrt(var))) }
# --- END UNSAFE ---

#
# Pareto maximum-likelihood estimator
#
paretofit = function(xs,xm=min(xs)) {
  xs = xs[xs >= xm]
  return (length(xs)/sum(log(xs/xm))) }

# Pareto maximum-likelihood estimator with standard error, as a data frame
paretofit_se = function(xs,xm) {
  xs = xs[xs >= xm]
  alpha = length(xs)/sum(log(xs/xm))
  n = length(xs)
  var = alpha^2 / n
  return (data.frame(alpha = alpha, se = sqrt(var))) }

# logarithmically bin data
log_bin = function(xs,lambda,xm=min(xs)) {
  xs = xs[xs >= xm]
  return(xm * lambda^floor(log(xs/xm)/log(lambda))) }

#
# Self-similar discrete power law maximum-likelihood estimator
#
# Fit the discrete power law after binning the input
#
dplfit_bin = function(xs,lambda,xm=min(xs)) {
  # "bin" by rounding down to the nearest xm\lambda^k
  return(dplfit_(log_bin(xs,lambda,xm),lambda,xm)) }

dplfit_bin_se = function(xs,lambda,xm) {
  return(dplfit_se_(log_bin(xs,lambda,xm),lambda,xm)) }
