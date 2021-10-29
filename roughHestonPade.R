#####################################################################################
# Jim Gatheral, 2021
# Various rational approximations to the rough Heston solution
#####################################################################################

########################################################################
# Pade approximations to h(a,x)
########################################################################

h.Pade33 <- function(params)function(a,x){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+.5
  
  aa <- sqrt(a * (a + (0+1i)) - rho^2 * a^2)
  rm <- -(0+1i) * rho * a - aa
  rp <- -(0+1i) * rho * a + aa
  
  b1 <- -a*(a+1i)/(2 * gamma(1+al))
  b2 <- (1-a*1i) * a^2 * rho/(2* gamma(1+2*al))               
  b3 <- gamma(1+2*al)/gamma(1+3*al) * 
    (a^2*(1i+a)^2/(8*gamma(1+al)^2)+(a+1i)*a^3*rho^2/(2*gamma(1+2*al)))
  
  g0 <- rm
  g1 <- -rm/(aa*gamma(1-al))
  g2 <- rm/aa^2/gamma(1-2*al) * (1 + rm/(2*aa)*gamma(1-2*al)/gamma(1-al)^2)
  
  den <- g0^3 +2*b1*g0*g1-b2*g1^2+b1^2*g2+b2*g0*g2
  
  p1 <- b1
  p2 <- (b1^2*g0^2 + b2*g0^3 + b1^3*g1 + b1*b2*g0*g1 - b2^2*g1^2 +b1*b3*g1^2 +b2^2*g0*g2 - b1*b3*g0*g2)/den
  q1 <- (b1*g0^2 + b1^2*g1 - b2*g0*g1 + b3*g1^2 - b1*b2*g2 -b3*g0*g2)/den
  q2 <- (b1^2*g0 + b2*g0^2 - b1*b2*g1 - b3*g0*g1 + b2^2*g2 - b1*b3*g2)/den
  q3 <- (b1^3 + 2*b1*b2*g0 + b3*g0^2 -b2^2*g1 +b1*b3*g1 )/den
  p3 <- g0*q3
  
  y <- x^al
  
  h.pade <- (p1*y + p2*y^2 + p3*y^3)/(1 + q1*y + q2*y^2 + q3*y^3)
  
  #res <- 1/2*(h.pade-rm)*(h.pade-rp)
  
  return(h.pade)
}

########################################################################

h.Pade44 <- function (params)function(a, x) {
  H <- params$H
  rho <- params$rho
  al <- H + 0.5
  aa <- sqrt(a * (a + (0 + (0+1i))) - rho^2 * a^2)
  rm <- -(0 + (0+1i)) * rho * a - aa
  rp <- -(0 + (0+1i)) * rho * a + aa
  b1 <- -a * (a + (0+1i))/(2 * gamma(1 + al))
  b2 <- (1 - a * (0+1i)) * a^2 * rho/(2 * gamma(1 + 2 * al))
  b3 <- gamma(1 + 2 * al)/gamma(1 + 3 * al) * (a^2 * (0+1i + 
                                                        a)^2/(8 * gamma(1 + al)^2) + (a + (0+1i)) * a^3 * rho^2/(2 * 
                                                                                                                   gamma(1 + 2 * al)))
  b4 <- ((a^2 * (0+1i + a)^2)/(8 * gamma(1 + al)^2) + ((0+1i) * 
                                                         rho^2 * (1 - (0+1i) * a) * a^3)/(2 * gamma(1 + 2 * al))) * 
    gamma(1 + 2 * al)/gamma(1 + 3 * al)
  g0 <- rm
  g1 <- -rm/(aa * gamma(1 - al))
  g2 <- rm/aa^2/gamma(1 - 2 * al) * (1 + rm/(2 * aa) * gamma(1 - 
                                                               2 * al)/gamma(1 - al)^2)
  g3 <- (rm * (-1 - (rm * gamma(1 - 2 * al))/(2 * aa * gamma(1 - 
                                                               al)^2) - (rm * gamma(1 - 3 * al) * (1 + (rm * gamma(1 - 
                                                                                                                     2 * al))/(2 * aa * gamma(1 - al)^2)))/(aa * gamma(1 - 
                                                                                                                                                                         2 * al) * gamma(1 - al))))/(aa^3 * gamma(1 - 3 * al))
  den <- (g0^4 + 3 * b1 * g0^2 * g1 + b1^2 * g1^2 - 2 * b2 * 
            g0 * g1^2 + b3 * g1^3 + 2 * b1^2 * g0 * g2 + 2 * b2 * 
            g0^2 * g2 - 2 * b1 * b2 * g1 * g2 - 2 * b3 * g0 * g1 * 
            g2 + b2^2 * g2^2 - b1 * b3 * g2^2 + b1^3 * g3 + 2 * b1 * 
            b2 * g0 * g3 + b3 * g0^2 * g3 - b2^2 * g1 * g3 + b1 * 
            b3 * g1 * g3)
  p1 <- b1
  p2 <- (b1^2 * g0^3 + b2 * g0^4 + 2 * b1^3 * g0 * g1 + 2 * 
           b1 * b2 * g0^2 * g1 - b1^2 * b2 * g1^2 - 2 * b2^2 * g0 * 
           g1^2 + b1 * b3 * g0 * g1^2 + b2 * b3 * g1^3 - b1 * b4 * 
           g1^3 + b1^4 * g2 + 2 * b1^2 * b2 * g0 * g2 + 2 * b2^2 * 
           g0^2 * g2 - b1 * b3 * g0^2 * g2 - b1 * b2^2 * g1 * g2 + 
           b1^2 * b3 * g1 * g2 - 2 * b2 * b3 * g0 * g1 * g2 + 2 * 
           b1 * b4 * g0 * g1 * g2 + b2^3 * g2^2 - 2 * b1 * b2 * 
           b3 * g2^2 + b1^2 * b4 * g2^2 + b1 * b2^2 * g0 * g3 - 
           b1^2 * b3 * g0 * g3 + b2 * b3 * g0^2 * g3 - b1 * b4 * 
           g0^2 * g3 - b2^3 * g1 * g3 + 2 * b1 * b2 * b3 * g1 * 
           g3 - b1^2 * b4 * g1 * g3)/den
  p3 <- (b1^3 * g0^2 + 2 * b1 * b2 * g0^3 + b3 * g0^4 + b1^4 * 
           g1 + 2 * b1^2 * b2 * g0 * g1 - b2^2 * g0^2 * g1 + 2 * 
           b1 * b3 * g0^2 * g1 - 2 * b1 * b2^2 * g1^2 + 2 * b1^2 * 
           b3 * g1^2 - b2 * b3 * g0 * g1^2 + b1 * b4 * g0 * g1^2 + 
           b3^2 * g1^3 - b2 * b4 * g1^3 + b1 * b2^2 * g0 * g2 - 
           b1^2 * b3 * g0 * g2 + b2 * b3 * g0^2 * g2 - b1 * b4 * 
           g0^2 * g2 + b2^3 * g1 * g2 - 2 * b1 * b2 * b3 * g1 * 
           g2 + b1^2 * b4 * g1 * g2 - 2 * b3^2 * g0 * g1 * g2 + 
           2 * b2 * b4 * g0 * g1 * g2 - b2^3 * g0 * g3 + 2 * b1 * 
           b2 * b3 * g0 * g3 - b1^2 * b4 * g0 * g3 + b3^2 * g0^2 * 
           g3 - b2 * b4 * g0^2 * g3)/den
  q1 <- (b1 * g0^3 + 2 * b1^2 * g0 * g1 - b2 * g0^2 * g1 - 
           2 * b1 * b2 * g1^2 + b3 * g0 * g1^2 - b4 * g1^3 + b1^3 * 
           g2 - b3 * g0^2 * g2 + b2^2 * g1 * g2 + b1 * b3 * g1 * 
           g2 + 2 * b4 * g0 * g1 * g2 - b2 * b3 * g2^2 + b1 * b4 * 
           g2^2 - b1^2 * b2 * g3 - b2^2 * g0 * g3 - b1 * b3 * g0 * 
           g3 - b4 * g0^2 * g3 + b2 * b3 * g1 * g3 - b1 * b4 * g1 * 
           g3)/den
  q2 <- (b1^2 * g0^2 + b2 * g0^3 + b1^3 * g1 - b3 * g0^2 * 
           g1 + b1 * b3 * g1^2 + b4 * g0 * g1^2 - b1^2 * b2 * g2 + 
           b2^2 * g0 * g2 - 3 * b1 * b3 * g0 * g2 - b4 * g0^2 * 
           g2 - b2 * b3 * g1 * g2 + b1 * b4 * g1 * g2 + b3^2 * g2^2 - 
           b2 * b4 * g2^2 + b1 * b2^2 * g3 - b1^2 * b3 * g3 + b2 * 
           b3 * g0 * g3 - b1 * b4 * g0 * g3 - b3^2 * g1 * g3 + b2 * 
           b4 * g1 * g3)/den
  q3 <- (b1^3 * g0 + 2 * b1 * b2 * g0^2 + b3 * g0^3 - b1^2 * 
           b2 * g1 - 2 * b2^2 * g0 * g1 - b4 * g0^2 * g1 + b2 * 
           b3 * g1^2 - b1 * b4 * g1^2 + b1 * b2^2 * g2 - b1^2 * 
           b3 * g2 + b2 * b3 * g0 * g2 - b1 * b4 * g0 * g2 - b3^2 * 
           g1 * g2 + b2 * b4 * g1 * g2 - b2^3 * g3 + 2 * b1 * b2 * 
           b3 * g3 - b1^2 * b4 * g3 + b3^2 * g0 * g3 - b2 * b4 * 
           g0 * g3)/den
  q4 <- (b1^4 + 3 * b1^2 * b2 * g0 + b2^2 * g0^2 + 2 * b1 * 
           b3 * g0^2 + b4 * g0^3 - 2 * b1 * b2^2 * g1 + 2 * b1^2 * 
           b3 * g1 - 2 * b2 * b3 * g0 * g1 + 2 * b1 * b4 * g0 * 
           g1 + b3^2 * g1^2 - b2 * b4 * g1^2 + b2^3 * g2 - 2 * b1 * 
           b2 * b3 * g2 + b1^2 * b4 * g2 - b3^2 * g0 * g2 + b2 * 
           b4 * g0 * g2)/den
  p4 <- g0 * q4
  y <- x^al
  h.pade <- (p1 * y + p2 * y^2 + p3 * y^3 + p4 * y^4)/(1 + 
                                                         q1 * y + q2 * y^2 + q3 * y^3 + q4 * y^4)
  #    res <- 1/2 * (h.pade - rm) * (h.pade - rp)
  return(h.pade)
}

########################################################################

h.Pade22 <- function(params)function(a,x){
  
  H <- params$H
  rho <- params$rho
  al <- H+.5
  
  aa <- sqrt(a * (a + (0+1i)) - rho^2 * a^2)
  rm <- -(0+1i) * rho * a - aa
  rp <- -(0+1i) * rho * a + aa
  
  b1 <- -a*(a+1i)/(2 * gamma(1+al))
  b2 <- (1-a*1i) * a^2 * rho/(2* gamma(1+2*al))               
  
  g0 <- rm
  g1 <- -rm/(aa*gamma(1-al))
  
  den <- g0^2 +b1*g1
  
  p1 <- b1
  p2 <- (b1^2*g0+b2*g0^2)/den
  q1 <- (b1*g0-b2*g1)/den
  q2 <- (b1^2+b2*g0)/den
  
  y <- x^al
  
  h.pade <- (p1*y + p2*y^2)/(1 + q1*y + q2*y^2)
  
  #res <- 1/2*(h.pade-rm)*(h.pade-rp) # F[h] = D^alpha h
  
  return(h.pade)
}

########################################################################
# Pade approximations to D^\alpha h(a,x)
########################################################################

d.h.Pade22 <- function(params)function(a,x){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+.5
  
  aa <- sqrt(a * (a + (0+1i)) - rho^2 * a^2)
  rm <- -(0+1i) * rho * a - aa
  rp <- -(0+1i) * rho * a + aa
  
  b1 <- -a*(a+1i)/(2 * gamma(1+al))
  b2 <- (1-a*1i) * a^2 * rho/(2* gamma(1+2*al))               
  
  g0 <- rm
  g1 <- -rm/(aa*gamma(1-al))
  
  den <- g0^2 +b1*g1
  
  p1 <- b1
  p2 <- (b1^2*g0+b2*g0^2)/den
  q1 <- (b1*g0-b2*g1)/den
  q2 <- (b1^2+b2*g0)/den
  
  y <- x^al
  
  h.pade <- (p1*y + p2*y^2)/(1 + q1*y + q2*y^2)
  
  res <- 1/2*(h.pade-rm)*(h.pade-rp) # F[h] = D^alpha h
  
  return(res)
}


########################################################################
d.h.Pade33 <- function(params)function(a,x){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+.5
  
  aa <- sqrt(a * (a + (0+1i)) - rho^2 * a^2)
  rm <- -(0+1i) * rho * a - aa
  rp <- -(0+1i) * rho * a + aa
  
  b1 <- -a*(a+1i)/(2 * gamma(1+al))
  b2 <- (1-a*1i) * a^2 * rho/(2* gamma(1+2*al))               
  b3 <- gamma(1+2*al)/gamma(1+3*al) * 
    (a^2*(1i+a)^2/(8*gamma(1+al)^2)+(a+1i)*a^3*rho^2/(2*gamma(1+2*al)))
  
  g0 <- rm
  g1 <- -rm/(aa*gamma(1-al))
  g2 <- rm/aa^2/gamma(1-2*al) * (1 + rm/(2*aa)*gamma(1-2*al)/gamma(1-al)^2)
  
  den <- g0^3 +2*b1*g0*g1-b2*g1^2+b1^2*g2+b2*g0*g2
  
  p1 <- b1
  p2 <- (b1^2*g0^2 + b2*g0^3 + b1^3*g1 + b1*b2*g0*g1 - b2^2*g1^2 +b1*b3*g1^2 +b2^2*g0*g2 - b1*b3*g0*g2)/den
  q1 <- (b1*g0^2 + b1^2*g1 - b2*g0*g1 + b3*g1^2 - b1*b2*g2 -b3*g0*g2)/den
  q2 <- (b1^2*g0 + b2*g0^2 - b1*b2*g1 - b3*g0*g1 + b2^2*g2 - b1*b3*g2)/den
  q3 <- (b1^3 + 2*b1*b2*g0 + b3*g0^2 -b2^2*g1 +b1*b3*g1 )/den
  p3 <- g0*q3
  
  y <- x^al
  
  h.pade <- (p1*y + p2*y^2 + p3*y^3)/(1 + q1*y + q2*y^2 + q3*y^3)
  
  res <- 1/2*(h.pade-rm)*(h.pade-rp)
  
  return(res)
}

########################################################################

d.h.Pade44 <- function(params)function(a,x){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+.5
  
  aa <- sqrt(a * (a + (0+1i)) - rho^2 * a^2)
  rm <- -(0+1i) * rho * a - aa
  rp <- -(0+1i) * rho * a + aa
  
  b1 <- -a*(a+1i)/(2 * gamma(1+al))
  b2 <- (1-a*1i) * a^2 * rho/(2* gamma(1+2*al))               
  b3 <- gamma(1+2*al)/gamma(1+3*al) * 
    (a^2*(1i+a)^2/(8*gamma(1+al)^2)+(a+1i)*a^3*rho^2/(2*gamma(1+2*al)))
  
  b4 <- ((a^2*(1i+a)^2)/(8*gamma(1+al)^2) + (1i*rho^2*(1-1i*a)*a^3)/(2*gamma(1+2*al))) * gamma(1+2*al)/gamma(1+3*al)
  
  g0 <- rm
  g1 <- -rm/(aa*gamma(1-al))
  g2 <- rm/aa^2/gamma(1-2*al) * (1 + rm/(2*aa)*gamma(1-2*al)/gamma(1-al)^2)
  
  g3 <- (rm*(-1 - (rm*gamma(1 - 2*al))/(2.*aa*gamma(1 - al)^2) - 
               (rm*gamma(1 - 3*al)*(1 + (rm*gamma(1 - 2*al))/(2.*aa*gamma(1 - al)^2)))/(aa*gamma(1 - 2*al)*gamma(1 - al))))/(aa^3*gamma(1 - 3*al))    
  
  #     g4 <-    (rm*(1 + (rm*((Gamma(1 - 4*al)*(1 + (rm*Gamma(1 - 2*al))/(2.*aa*Gamma(1 - al)^2))^2)/Gamma(1 - 2*al)^2 - 
  #                 (2*Gamma(1 - 4*al)*(-1 - (rm*Gamma(1 - 2*al))/(2.*aa*Gamma(1 - al)^2) - 
  #                      (rm*Gamma(1 - 3*al)*(1 + (rm*Gamma(1 - 2*al))/(2.*aa*Gamma(1 - al)^2)))/(aa*Gamma(1 - 2*al)*Gamma(1 - al))))/
  #                  (Gamma(1 - 3*al)*Gamma(1 - al))))/(2.*aa) + (rm*Gamma(1 - 2*al))/(2.*aa*Gamma(1 - al)^2) + 
  #            (rm*Gamma(1 - 3*al)*(1 + (rm*Gamma(1 - 2*al))/(2.*aa*Gamma(1 - al)^2)))/(aa*Gamma(1 - 2*al)*Gamma(1 - al))))/(aa^4*Gamma(1 - 4*al))
  
  den <- (g0^4 + 3*b1*g0^2*g1 + b1^2*g1^2 - 2*b2*g0*g1^2 + b3*g1^3 + 
            2*b1^2*g0*g2 + 2*b2*g0^2*g2 - 2*b1*b2*g1*g2 - 2*b3*g0*g1*g2 + 
            b2^2*g2^2 - b1*b3*g2^2 + b1^3*g3 + 2*b1*b2*g0*g3 + b3*g0^2*g3 - 
            b2^2*g1*g3 + b1*b3*g1*g3)
  
  p1 <- b1
  p2 <- (b1^2*g0^3 + b2*g0^4 + 2*b1^3*g0*g1 + 2*b1*b2*g0^2*g1 - 
           b1^2*b2*g1^2 - 2*b2^2*g0*g1^2 + b1*b3*g0*g1^2 + b2*b3*g1^3 - 
           b1*b4*g1^3 + b1^4*g2 + 2*b1^2*b2*g0*g2 + 2*b2^2*g0^2*g2 - 
           b1*b3*g0^2*g2 - b1*b2^2*g1*g2 + b1^2*b3*g1*g2 - 2*b2*b3*g0*g1*g2 + 
           2*b1*b4*g0*g1*g2 + b2^3*g2^2 - 2*b1*b2*b3*g2^2 + b1^2*b4*g2^2 + 
           b1*b2^2*g0*g3 - b1^2*b3*g0*g3 + b2*b3*g0^2*g3 - b1*b4*g0^2*g3 - 
           b2^3*g1*g3 + 2*b1*b2*b3*g1*g3 - b1^2*b4*g1*g3)/den
  
  p3 <- (b1^3*g0^2 + 2*b1*b2*g0^3 + b3*g0^4 + b1^4*g1 + 2*b1^2*b2*g0*g1 - b2^2*g0^2*g1 + 2*b1*b3*g0^2*g1 - 
           2*b1*b2^2*g1^2 + 2*b1^2*b3*g1^2 - b2*b3*g0*g1^2 + b1*b4*g0*g1^2 + b3^2*g1^3 - b2*b4*g1^3 + 
           b1*b2^2*g0*g2 - b1^2*b3*g0*g2 + b2*b3*g0^2*g2 - b1*b4*g0^2*g2 + b2^3*g1*g2 - 2*b1*b2*b3*g1*g2 + 
           b1^2*b4*g1*g2 - 2*b3^2*g0*g1*g2 + 2*b2*b4*g0*g1*g2 - b2^3*g0*g3 + 2*b1*b2*b3*g0*g3 - b1^2*b4*g0*g3 + 
           b3^2*g0^2*g3 - b2*b4*g0^2*g3)/den
  
  
  q1 <- (b1*g0^3 + 2*b1^2*g0*g1 - b2*g0^2*g1 - 2*b1*b2*g1^2 + b3*g0*g1^2 - b4*g1^3 + b1^3*g2 - 
           b3*g0^2*g2 + b2^2*g1*g2 + b1*b3*g1*g2 + 2*b4*g0*g1*g2 - b2*b3*g2^2 + b1*b4*g2^2 - b1^2*b2*g3 - b2^2*g0*g3 - 
           b1*b3*g0*g3 - b4*g0^2*g3 + b2*b3*g1*g3 - b1*b4*g1*g3)/den
  
  q2 <- (b1^2*g0^2 + b2*g0^3 + b1^3*g1 - b3*g0^2*g1 + b1*b3*g1^2 + b4*g0*g1^2 - b1^2*b2*g2 + b2^2*g0*g2 - 
           3*b1*b3*g0*g2 - b4*g0^2*g2 - b2*b3*g1*g2 + b1*b4*g1*g2 + b3^2*g2^2 - b2*b4*g2^2 + b1*b2^2*g3 - b1^2*b3*g3 + 
           b2*b3*g0*g3 - b1*b4*g0*g3 - b3^2*g1*g3 + b2*b4*g1*g3)/den
  
  q3 <- (b1^3*g0 + 2*b1*b2*g0^2 + b3*g0^3 - b1^2*b2*g1 - 2*b2^2*g0*g1 - b4*g0^2*g1 + b2*b3*g1^2 - b1*b4*g1^2 + 
           b1*b2^2*g2 - b1^2*b3*g2 + b2*b3*g0*g2 - b1*b4*g0*g2 - b3^2*g1*g2 + b2*b4*g1*g2 - b2^3*g3 + 2*b1*b2*b3*g3 - 
           b1^2*b4*g3 + b3^2*g0*g3 - b2*b4*g0*g3)/den
  
  q4 <- (b1^4 + 3*b1^2*b2*g0 + b2^2*g0^2 + 2*b1*b3*g0^2 + b4*g0^3 - 2*b1*b2^2*g1 + 2*b1^2*b3*g1 - 
           2*b2*b3*g0*g1 + 2*b1*b4*g0*g1 + b3^2*g1^2 - b2*b4*g1^2 + b2^3*g2 - 2*b1*b2*b3*g2 + b1^2*b4*g2 - b3^2*g0*g2 + 
           b2*b4*g0*g2)/den
  
  p4 <- g0*q4
  
  y <- x^al
  
  h.pade <- (p1*y + p2*y^2+ p3*y^3 +p4*y^4)/(1 + q1*y + q2*y^2 + q3*y^3 + q4*y^4)
  
  res <- 1/2*(h.pade-rm)*(h.pade-rp)
  
  return(res)
}



########################################################################
# Pade approximations to D^\alpha h(a,x)
########################################################################

dh.Pade12 <- function(params)function(a,x){
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H+.5
  
  aa <- sqrt(a * (a + (0+1i)) - rho^2 * a^2)
  rm = -(0+1i) * rho * a - aa
  rp = -(0+1i) * rho * a + aa
  
  b0 <- -a*(a+1i)/2
  b1 <- (1-a*1i) * a^2 * rho/(2* gamma(1+al))               
  
  g1 <- rm/gamma(1-al)
  g2 <- -rm/(aa*gamma(1-2*al))
  
  p0 <- b0
  p1 <- (b0^2*g1+b1*g1^2)/(g1^2+b0*g2)
  q1 <- (b0*g1-b1*g2)/(g1^2+b0*g2)
  q2 <- (b0^2+b1*g1)/(g1^2+b0*g2)
  
  y <- x^al
  
  res <- (p0 + p1*y)/(1 + q1*y + q2*y^2)

  return(res)
}

########################################################################

dh.Pade23 <- function (params) function(a, x) {
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  
  al <- H + 0.5
  aa <- sqrt(a * (a + (0 + (0+1i))) - rho^2 * a^2)
  rm = -(0 + (0+1i)) * rho * a - aa
  
  b0 <- -a * (a + (0+1i))/2
  b1 <- (1 - a * (0+1i)) * a^2 * rho/(2 * gamma(1 + al))
  b2 <-  (a^2 * (1i + a)^2)/(8 * gamma(1 + al)^2) - (((0+1i)/2) * rho^2 * (1 - (0+1i) * a) * a^3)/gamma(1 + 2 * al)
  
  g1 <- rm/gamma(1 - al)
  g2 <- -rm/(aa * gamma(1 - 2*al))
  gam2 <- (1 + rm/(2*aa)*gamma(1-2*al)/gamma(1-al)^2) # Fix this formula!!!
  #g3 <- gam2*rm/(aa^2 * gamma(1 - 3*al))
  
  g3 <- (rm*(1 + (rm*gamma(1 - 2*al))/(2.*aa*gamma(1 - al)^2)))/(aa^2*gamma(1 - 3*al))
  
  den <- g1^3 + 2*b0*g1*g2 - b1*g2^2 + b0^2*g3 + b1*g1*g3
  
  p0 <- b0
  p1 <- (b0^2*g1^2 + b1*g1^3 + b0^3*g2 + b0*b1*g1*g2 - b1^2*g2^2 + b0*b2*g2^2 + b1^2*g1*g3 - b0*b2*g1*g3)/den
  q1 <- (b0*g1^2 + b0^2*g2 - b1*g1*g2 + b2*g2^2 - b0*b1*g3 - b2*g1*g3)/den
  q2 <- (b0^2*g1 + b1*g1^2 - b0*b1*g2 - b2*g1*g2 + b1^2*g3 - b0*b2*g3)/den
  q3 <- (b0^3 + 2*b0*b1*g1 + b2*g1^2 - b1^2*g2 + b0*b2*g2)/den
  p2 <- g1 * q3
  
  y <- x^al
  
  res <- (p0 + p1*y + p2*y^2 )/(1 + q1*y + q2*y^2 + q3*y^3)
  
  return(res)
}

########################################################################
#
# Characteristic function using Padé approximation
#
########################################################################

phiRoughHestonDhApprox.raw <- function(params, xiCurve, dh.approx, n=100) function(u, t) {
  
  H <- params$H
  rho <- params$rho
  nu <- params$nu
  al <- H + 1/2
  
  ti <- (0:n)/n * t
  x <- nu^(1/al)*ti
  xi <- xiCurve(ti)
  dah <- dh.approx(params)(u, x) 
  
  return(exp(t(dah) %*% rev(xi) * t/n))
}

phiRoughHestonDhApprox <- function(params, xiCurve, dh.approx, n=100) function(u, t){
  phi1 <- function(u){ifelse(u==0,1,phiRoughHestonDhApprox.raw(params, xiCurve, dh.approx, n)(u,t))}
  return(sapply(u,phi1))
}