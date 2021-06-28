BlackCall <- function(S0, K, T, sigma)
{
k <- log(K/S0)
sig <- sigma*sqrt(T)
d1 <- -k/sig+sig/2
d2 <- d1 - sig
return( S0*pnorm(d1) - K*pnorm(d2))
}

BlackPut <- function (S0, K, T, sigma) 
{
  k <- log(K/S0)
  sig <- sigma * sqrt(T)
  d1 <- -k/sig + sig/2
  d2 <- d1 - sig
  return( K * pnorm(-d2) - S0 * pnorm(-d1))
}

# This function works with vectors of strikes and option values
ivCall <- function(S0, K, T, C)
{
nK <- length(K)
sigmaL <- rep(1e-10,nK)
CL <- BlackCall(S0, K, T, sigmaL)
sigmaH <- rep(10,nK)
CH <- BlackCall(S0, K, T, sigmaH)
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2
    CM <- BlackCall(S0, K, T, sigma)
    CL <- CL + (CM < C)*(CM-CL)
    sigmaL <- sigmaL + (CM < C)*(sigma-sigmaL)
    CH <- CH + (CM >= C)*(CM-CH)
    sigmaH <- sigmaH + (CM >= C)*(sigma-sigmaH)
  }
  return(sigma)
}

# This function also works with vectors of strikes and option values  
ivPut <- function (S0, K, T, P) 
{
  nK <- length(K)
  sigmaL <- 1e-10
  PL <- BlackPut(S0, K, T, sigmaL)
  sigmaH <- 10
  PH <- BlackPut(S0, K, T, sigmaH)
  while (mean(sigmaH - sigmaL) > 1e-10) {
    sigma <- (sigmaL + sigmaH)/2
    PM <- BlackPut(S0, K, T, sigma)
    PL <- PL + (PM < P) * (PM - PL)
    sigmaL <- sigmaL + (PM < P) * (sigma - sigmaL)
    PH <- PH + (PM >= P) * (PM - PH)
    sigmaH <- sigmaH + (PM >= P) * (sigma - sigmaH)
  }
  return(sigma)
}

ivOTM <- function(S0, K, T, V){
  res <- ifelse(K>S0,ivCall(S0, K, T, V), ivPut(S0, K, T, V))
  return(res)
}

# Function to compute option prices and implied vols given vector of final values of underlying
ivS <- function (Sf, T, AK) 
{
  nK <- length(AK)
  N <- length(Sf)
  Sfbar <- mean(Sf)
  V <- numeric(nK)
  ivBlack <- numeric(nK)
  for (j in 1:nK) {
    payoff <- (Sf - AK[j]) * (Sf > AK[j])
    V <- mean(payoff)
    V.OTM <- ifelse(Sfbar<AK[j], V, V + AK[j] - Sfbar)
    ivBlack[j] <- ivOTM(Sfbar, AK[j], T, V.OTM)
  }
  return(ivBlack)
}
