#####################################################################################
# The following code is a translation to R of Matlab code
# written by Fabio Baschetti, Giacomo Bormetti, Pietro Rossi, 
# and Silvia Romagnoli at the University of Bologna (2020).
#####################################################################################

#####################################################################################
# h(.) using the Adams scheme
hA.GB <- function(H,rho,nu,a,n,bigT){
  
  DELTA <- bigT/n
  alpha <- H + 1/2
  
  hA <- numeric(n+1)
  
  F <- function(u,h){ -0.5*u*(u+1i)+1i*rho*nu*u*h+0.5*nu^2*h^2}
  
  a01    <- DELTA^alpha / gamma(alpha+2) * alpha
  a11    <- DELTA^alpha / gamma(alpha+2)
  b01    <- DELTA^alpha / gamma(alpha+1)
  # hA[1] is zero, boundary condition
  hP     <- b01*F(a,hA[1])
  hA[2]  <- a01*F(a,hA[1]) + a11*F(a,hP)
  # ------
  
  k      <- c(1:(n-1))
  a0     <- DELTA^alpha / gamma(alpha+2) * (k^(alpha+1) - (k-alpha)*(k+1)^alpha)
  q      <- c(0:(n-1))
  aj     <- DELTA^alpha / gamma(alpha+2) * ((q+2)^(alpha+1) + (q)^(alpha+1) - 2*(q+1)^(alpha+1))
  r      <- c(0:n)
  bj     <- DELTA^alpha / gamma(alpha+1) * ((r+1)^alpha - (r)^alpha)
  
  akp1   <- DELTA^alpha / gamma(alpha+2)
  
  for (k in 1:(n-1)){
    f    <- F(a,hA[1:(k+1)])
    aA   <- c(a0[k], rev(aj[1:k]))
    hP   <- t(f) %*% rev(bj[1:(k+1)])
    hA[k+2] <- t(f) %*% aA + akp1*F(a,hP)
  }
  return(hA)
}

dhA.GB.raw <- function(H, rho, nu, n, a, bigT){
  DELTA <- bigT/n
  alpha <- H + 1/2
  
  hA <- numeric(n+1)
  dhA <- numeric(n+1)  
  
  F <- function(u,h){ -0.5*u*(u+1i)+1i*rho*nu*u*h+0.5*nu^2*h^2}
  
  a01    <- DELTA^alpha / gamma(alpha+2) * alpha
  a11    <- DELTA^alpha / gamma(alpha+2)
  b01    <- DELTA^alpha / gamma(alpha+1)
  dhA[1] <- F(a,hA[1])
  # hA[1] is zero, boundary condition
  hP     <- b01*F(a,hA[1])
  hA[2]  <- a01*F(a,hA[1]) + a11*F(a,hP)
  dhA[2] <- F(a,hA[2])
  # ------
  
  k      <- c(1:(n-1))
  a0     <- DELTA^alpha / gamma(alpha+2) * (k^(alpha+1) - (k-alpha)*(k+1)^alpha)
  q      <- c(0:(n-1))
  aj     <- DELTA^alpha / gamma(alpha+2) * ((q+2)^(alpha+1) + (q)^(alpha+1) - 2*(q+1)^(alpha+1))
  r      <- c(0:n)
  bj     <- DELTA^alpha / gamma(alpha+1) * ((r+1)^alpha - (r)^alpha)
  
  akp1   <- DELTA^alpha / gamma(alpha+2)
  
  for (k in 1:(n-1)){
    f    <- F(a,hA[1:(k+1)])
    aA   <- c(a0[k], rev(aj[1:k]))
    hP   <- t(f) %*% rev(bj[1:(k+1)])
    hA[k+2] <- t(f) %*% aA + akp1*F(a,hP)
    dhA[k+2] <- F(a,hA[k+2])
  }
  
  
  return(dhA)
}
# Vectorize the above function
dhA <- Vectorize(dhA.GB.raw,vectorize.args = "a")
#####################################################################################


#####################################################################################
# The following code takes the above and computes the rough Heston
# characteristic function.
#####################################################################################
phiRoughHeston.raw <- function (params, xiCurve, n) 
  function(u, t) {
    rho <- params$rho
    al <- params$al
    H <- al-1/2
    nu <- params$nu
    yy <-sqrt(1-rho)
    ti <- (0:n)/n * t
    xi <- xiCurve(ti)
    dah <- dhA(H=H, rho=rho, nu=nu, n=n, a=u, bigT=t)
    return(exp(t(dah) %*% rev(xi) * t/n))
  }
# Vectorize the above function
phiRoughHeston <- function (params, xiCurve, nSteps){
  Vectorize(phiRoughHeston.raw(params, xiCurve, nSteps),vectorize.args = "u")
}


#####################################################################################
# The following code is a translation to R of Python code 
# written by Omar El Euch of École Polytechnique Paris (2017).

# All errors in the translation are mine.

# The computation is explained in Appendix B of Roughening Heston.
#####################################################################################

#####################################################################################
# Compute upper bound for Amax to truncate the Fourier inversion in 
# Lewis option price formula with ivol error eps*sqrt(maturity). 
#####################################################################################

ik.alpha <- function(H,tau,xiCurve,k){
  alpha <- H+1/2
  integ <- function(s){xiCurve(s)*(tau-s)^(k*alpha)}
  res <- integrate(integ,lower=0,upper=tau)$value
  return(res)
}

# y.alpha is for amax computation
y.alpha <- function(H,tau,xiCurve){
  alpha <- 1/2+H
  integ <- function(s){xiCurve(s)/(tau-s)^alpha}
  res <- integrate(integ,lower=0,upper=tau)$value/gamma(1-alpha)
  return(res)
}

findzero<- function(c){
  f <- function (x) {exp(-x) - c*(x^2+1/4)}
  xmin <- uniroot(f, c(0, 1/sqrt(c)), tol = 0.0001)$root
  return(xmin)
}

amax.cf <- function(params, xiCurve, k, tau, eps.sig){ # eps.sig is roughly error in ivol
  rho <- params$rho
  al <- params$al
  H <- al-1/2
  nu <- params$nu
  yy <-sqrt(1-rho^2)/nu * y.alpha(H,tau,xiCurve) 
  c <- sqrt(pi/2)*exp(-k/2)*eps.sig#*sqrt(tau) # Note we divided by sqrt(tau)
  f <- function (a){exp(-yy*a) - c*yy*(a^2+1/4)}
  amax <- uniroot(f, c(0, 1/sqrt(c*yy)), tol = 0.0001)$root
  return(amax)
}

#####################################################################################
# Now we compute rough Heston option prices
#####################################################################################
otmRoughHeston.raw <- function (params, xiCurve, nSteps) 
  function(k, t) {
    k.minus <- (k < 0) * k
    phi <- phiRoughHeston(params, xiCurve, n = nSteps)
    a.max <- amax.cf(params, xiCurve, k, t, eps.sig = 1e-04)
    integrand <- function(u) {
      Re(exp(-(0 + (0+1i)) * u * k) * phi(u - (0 + (0+1i))/2, 
                                          t)/(u^2 + 1/4))
    }
    res <- exp(k.minus) - exp(k/2)/pi * integrate(integrand, lower = 0, upper = a.max, rel.tol = 1e-08)$value
    return(res)
  }
# Vectorize the above function
otmRoughHeston <- function (params, xiCurve, nSteps){
  Vectorize(otmRoughHeston.raw(params, xiCurve, nSteps),vectorize.args = "k")
}
#####################################################################################

#####################################################################################
# Finally compute implied volatilities
#####################################################################################
impliedVolRoughHeston <- function (params, xiCurve, nSteps) 
  function(k, t) {
    ivOTM(1, exp(k), t, otmRoughHeston(params, xiCurve, nSteps)(k, t))
  }
#####################################################################################