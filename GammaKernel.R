#######################################################
# Gamma kernel functions
#######################################################

# Gamma kernel
gGamma <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  return(sqrt(2*H)*tau^{al-1}*exp(-lam*tau))
}

# Gamma variance G00
G00 <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  H2 <- 2*H
  lam <- params$lam
  
  prefactor <- H2/((2*lam)^H2)
  bkt <- gamma(H2)- gamma_inc(H2,2*lam*tau)
  res2 <- tau^(2*H)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
}

# G11
G11 <- function(params) function(tau){
  al <- params$al
  H <- al - 1/2
  H2 <- 2*H
  lam <- params$lam
  
  prefactor <- H2/((2*lam)^H2)
  bkt <- gamma_inc(H2,2*lam*tau)- gamma_inc(H2,4*lam*tau)
  res2 <- tau^(2*H)*(2^H2-1)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
}

# G0
G0 <- function(params) function(dt){
  
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  
  prefactor <- sqrt(2*H)/(lam^al)
  bkt <- gamma(al)- gamma_inc(al,lam*dt)
  res2 <- sqrt(2*H)/al*dt^(al)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
  
}

G1 <- function(params) function(dt){
  
  al <- params$al
  H <- al - 1/2
  lam <- params$lam
  
  prefactor <- sqrt(2*H)/(lam^al)
  bkt <- gamma_inc(al,lam*dt)- gamma_inc(al,2*lam*dt)
  res2 <- sqrt(2*H)/al*dt^(al)*(2^al-1)
  res <- ifelse(lam>0,prefactor * bkt,res2)
  return(res)
  
}

# Gamma covariance
G0k <- function(params,k)function(t){
  
  gp <- gGamma(params)
  eps <- 0
  integr <- function(s){gp(s)*gp(s+k*t)}
  res <- integrate(integr, lower=0,upper=t)$value
  return(res)
  
}

# Gamma first order covariance
G01 <- function(params)function(t){
  
  gp <- gGamma(params)
  eps <- 0
  integr <- function(s){gp(s)*gp(s+t)}
  res <- integrate(integr, lower=0,upper=t)$value
  return(res)
  
}

