####################################################################################
# Jim Gatheral, 2021
#
# The following functions use the Lewis representation of the option price
# in terms of the characteristic function, equation (5.6) of 
# The Volatility Surface.
####################################################################################

source("BlackFormula.R")

option.OTM.raw <- function (phi, k, t) {
  
  integrand <- function(u) {
    Re(exp(-1i * u * k) * phi(u - 1i/2, t)/(u^2 + 1/4))
  }
  k.minus <- (k < 0) * k
  
  res <- exp(k.minus) - exp(k/2)/pi * integrate(integrand, lower = 0, upper = Inf,rel.tol = 1e-8)$value #
  return(res)
}

option.OTM <- Vectorize(option.OTM.raw,vectorize.args="k")

impvol.phi <- function (phi) 
{
  function(k, t) {
    return( ivOTM(1, exp(k), t, option.OTM(phi,k, t)))
  }
}