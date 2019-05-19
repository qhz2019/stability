
library(Matrix) 

## simulation 1
law <- function(alpha, N=1, omni.i=NA, omni.j=NA, omega=NULL){
  a<-NULL
  b<-NULL
  out.a <- NULL
  out.b <-NULL
  out.c <- NULL
  out.d <- NULL
  out.e <- NULL
  out.f <- NULL
  out.g <- NULL
  out.h <-NULL
  out.i <- NULL
  MM <- NULL
  MMb <- NULL
  muuu<- NULL
  big<- NULL
  MI<-NULL
  MIS<-NULL
  
  S <- nrow(alpha)
  k<-length(muu)
  
  if(is.na(omni.i)) {
    out <- matrix(NA, nrow=N, ncol=2)
    
    colnames(out) <- c("DomEig", "Im")
    
    for(n in 1:N) out[n,] <- {
      muuuu<-function(muu){
        muuu<-muu
        muuu[1:6]<-1
        return(muuu)
      }
      mu<-muuuu(muu)
      
      runn<- function (alpha) {
        al <- kk[n,1] * alpha 
        diag(al[1:6,1:6]) <- -1
        diag(al[7:9,7:9]) <- -0.1
        return(as.matrix(al))
      }
      alp<-runn(alpha)
      
      
      linearequations<-function(alp,mu){
        a <- rbind(c(alp[1,1], alp[1,2], alp[1,3], alp[1,4], alp[1,5], alp[1,6], alp[1,7], alp[1,8], alp[1,9]), 
                   c(alp[2,1], alp[2,2], alp[2,3], alp[2,4], alp[2,5], alp[2,6], alp[2,7], alp[2,8], alp[2,9]),
                   c(alp[3,1], alp[3,2], alp[3,3], alp[3,4], alp[3,5], alp[3,6], alp[3,7], alp[3,8], alp[3,9]),
                   c(alp[4,1], alp[4,2], alp[4,3], alp[4,4], alp[4,5], alp[4,6], alp[4,7], alp[4,8], alp[4,9]),
                   c(alp[5,1], alp[5,2], alp[5,3], alp[5,4], alp[5,5], alp[5,6], alp[5,7], alp[5,8], alp[5,9]),
                   c(alp[6,1], alp[6,2], alp[6,3], alp[6,4], alp[6,5], alp[6,6], alp[6,7], alp[6,8], alp[6,9]),
                   c(alp[7,1], alp[7,2], alp[7,3], alp[7,4], alp[7,5], alp[7,6], alp[7,7], alp[7,8], alp[7,9]),
                   c(alp[8,1], alp[8,2], alp[8,3], alp[8,4], alp[8,5], alp[8,6], alp[8,7], alp[8,8], alp[8,9]),
                   c(alp[9,1], alp[9,2], alp[9,3], alp[9,4], alp[9,5], alp[9,6], alp[9,7], alp[9,8], alp[9,9]))
        b <- c(-mu[1], -mu[2], -mu[3], -mu[4], -mu[5], -mu[6], -mu[7], -mu[8], -mu[9])
        dd<-solve(a, b) 
        return(dd)}
      d<-linearequations(alp,mu)
      
      big<-d[which.min(d)]
      bi<-d[which.max(d)]
      
      jaco<- function (mu,alp,d) {
        
        out.a1 <- mu[1]+alp[1,1]*d[1]*2+alp[1,2]*d[2]+alp[1,3]*d[3]+alp[1,4]*d[4]+alp[1,5]*d[5]+alp[1,6]*d[6]+alp[1,7]*d[7]+alp[1,8]*d[8]+alp[1,9]*d[9]
        out.a2 <- alp[1,2]*d[1]
        out.a3 <- alp[1,3]*d[1]
        out.a4 <- alp[1,4]*d[1]
        out.a5 <- alp[1,5]*d[1]
        out.a6 <- alp[1,6]*d[1]
        out.a7 <- alp[1,7]*d[1]
        out.a8 <- alp[1,8]*d[1]
        out.a9 <- alp[1,9]*d[1]
        
        
        out.b1 <- alp[2,1]*d[2]
        out.b2 <- mu[2]+alp[2,1]*d[1]+alp[2,2]*d[2]*2+alp[2,3]*d[3]+alp[2,4]*d[4]+alp[2,5]*d[5]+alp[2,6]*d[6]+alp[2,7]*d[7]+alp[2,8]*d[8]+alp[2,9]*d[9]
        out.b3 <- alp[2,3]*d[2]
        out.b4 <- alp[2,4]*d[2]
        out.b5 <- alp[2,5]*d[2]
        out.b6 <- alp[2,6]*d[2]
        out.b7 <- alp[2,7]*d[2]
        out.b8 <- alp[2,8]*d[2]
        out.b9 <- alp[2,9]*d[2]
        
        
        out.c1 <- alp[3,1]*d[3]
        out.c2 <- alp[3,2]*d[3]
        out.c3 <- mu[3]+alp[3,1]*d[1]+alp[3,2]*d[2]+alp[3,3]*d[3]*2+alp[3,4]*d[4]+alp[3,5]*d[5]+alp[3,6]*d[6]+alp[3,7]*d[7]+alp[3,8]*d[8]+alp[3,9]*d[9]
        out.c4 <- alp[3,4]*d[3]
        out.c5 <- alp[3,5]*d[3]
        out.c6 <- alp[3,6]*d[3]
        out.c7 <- alp[3,7]*d[3]
        out.c8 <- alp[3,8]*d[3]
        out.c9 <- alp[3,9]*d[3]
        
        out.d1 <- alp[4,1]*d[4]
        out.d2 <- alp[4,2]*d[4]
        out.d3 <- alp[4,3]*d[4]
        out.d4 <- mu[4]+alp[4,1]*d[1]+alp[4,2]*d[2]+alp[4,3]*d[3]+alp[4,4]*d[4]*2+alp[4,5]*d[5]+alp[4,6]*d[6]+alp[4,7]*d[7]+alp[4,8]*d[8]+alp[4,9]*d[9]
        out.d5 <- alp[4,5]*d[4]
        out.d6 <- alp[4,6]*d[4]
        out.d7 <- alp[4,7]*d[4]
        out.d8 <- alp[4,8]*d[4]
        out.d9 <- alp[4,9]*d[4]
        
        
        out.e1 <- alp[5,1]*d[4]
        out.e2 <- alp[5,2]*d[5]
        out.e3 <- alp[5,3]*d[5]
        out.e4 <- alp[5,4]*d[5]
        out.e5 <- mu[5]+alp[5,1]*d[1]+alp[5,2]*d[2]+alp[5,3]*d[3]+alp[5,4]*d[4]+alp[5,5]*d[5]*2+alp[5,6]*d[6]+alp[5,7]*d[7]+alp[5,8]*d[8]+alp[5,9]*d[9]
        out.e6 <- alp[5,6]*d[5]
        out.e7 <- alp[5,7]*d[5]
        out.e8 <- alp[5,8]*d[5]
        out.e9 <- alp[5,9]*d[5]
        
        out.f1 <- alp[6,1]*d[4]
        out.f2 <- alp[6,2]*d[6]
        out.f3 <- alp[6,3]*d[6]
        out.f4 <- alp[6,4]*d[6]
        out.f5 <- alp[6,5]*d[6]
        out.f6 <- mu[6]+alp[6,1]*d[1]+alp[6,2]*d[2]+alp[6,3]*d[3]+alp[6,4]*d[4]+alp[6,5]*d[5]+alp[6,6]*d[6]*2+alp[6,7]*d[7]+alp[6,8]*d[8]+alp[6,9]*d[9]
        out.f7 <- alp[6,7]*d[6]
        out.f8 <- alp[6,8]*d[6]
        out.f9 <- alp[6,9]*d[6]
        
        
        out.g1 <- alp[7,1]*d[7]
        out.g2 <- alp[7,2]*d[7]
        out.g3 <- alp[7,3]*d[7]
        out.g4 <- alp[7,4]*d[7]
        out.g5 <- alp[7,5]*d[7]
        out.g6 <- alp[7,6]*d[7]
        out.g7 <- mu[7]+alp[7,1]*d[1]+alp[7,2]*d[2]+alp[7,3]*d[3]+alp[7,4]*d[4]+alp[7,5]*d[5]+alp[7,6]*d[6]+alp[7,7]*d[7]*2+alp[7,8]*d[8]+alp[7,9]*d[9]
        out.g8 <- alp[7,8]*d[7]
        out.g9 <- alp[7,9]*d[7]
        
        
        out.h1 <- alp[8,1]*d[8]
        out.h2 <- alp[8,2]*d[8]
        out.h3 <- alp[8,3]*d[8]
        out.h4 <- alp[8,4]*d[8]
        out.h5 <- alp[8,5]*d[8]
        out.h6 <- alp[8,6]*d[8]
        out.h7 <- alp[8,7]*d[8]
        out.h8 <- mu[8]+alp[8,1]*d[1]+alp[8,2]*d[2]+alp[8,3]*d[3]+alp[8,4]*d[4]+alp[8,5]*d[5]+alp[8,6]*d[6]+alp[8,7]*d[7]+alp[8,8]*d[8]*2+alp[8,9]*d[9]
        out.h9 <- alp[8,9]*d[8]
        
        
        
        out.i1 <- alp[9,1]*d[9]
        out.i2 <- alp[9,2]*d[9]
        out.i3 <- alp[9,3]*d[9]
        out.i4 <- alp[9,4]*d[9]
        out.i5 <- alp[9,5]*d[9]
        out.i6 <- alp[9,6]*d[9]
        out.i7 <- alp[9,7]*d[9]
        out.i8 <- alp[9,8]*d[9]
        out.i9 <- mu[9]+alp[9,1]*d[1]+alp[9,2]*d[2]+alp[9,3]*d[3]+alp[9,4]*d[4]+alp[9,5]*d[5]+alp[9,6]*d[6]+alp[9,7]*d[7]+alp[9,8]*d[8]+alp[9,9]*d[9]*2
        
        MM<-rbind(c(out.a1,out.a2, out.a3, out.a4, out.a5, out.a6, out.a7, out.a8, out.a9),
                  c(out.b1,out.b2, out.b3, out.b4, out.b5, out.b6, out.b7, out.b8, out.b9),
                  c(out.c1,out.c2, out.c3, out.c4, out.c5, out.c6, out.c7, out.c8, out.c9),
                  c(out.d1,out.d2, out.d3, out.d4, out.d5, out.d6, out.d7, out.d8, out.d9),
                  c(out.e1,out.e2, out.e3, out.e4, out.e5, out.e6, out.e7, out.e8, out.e9),
                  c(out.f1,out.f2, out.f3, out.f4, out.f5, out.f6, out.f7, out.f8, out.f9),
                  c(out.g1,out.g2, out.g3, out.g4, out.g5, out.g6, out.g7, out.g8, out.g9),
                  c(out.h1,out.h2, out.h3, out.h4, out.h5, out.h6, out.h7, out.h8, out.h9),
                  c(out.i1,out.i2, out.i3, out.i4, out.i5, out.i6, out.i7, out.i8, out.i9))
        return(as.matrix(MM))
      }
      
      MMb<-jaco(mu,alp,d)
      eigs <- eigen(MMb)[["values"]]
      mx <- which.max(Re(eigs))
      deM <- Re(eigs)[mx]
      deM<-ifelse(big<0 |bi ==Inf, NA,  deM)    
      deMi <- Im(eigs)[mx]
      c(deM, deMi)
    } } else { print("exist omnivores")}
  return(as.data.frame(out))
}

out.E <- law(alpha, N = 10000)
summary(out.E)
RT.E <- -1/(out.E[["DomEig"]]) 
summary(RT.E)

## simulation2 
law <- function(alpha, N=1, omni.i=NA, omni.j=NA, omega=NULL){
  
  a<-NULL
  b<-NULL
  out.a <- NULL
  out.b <-NULL
  out.c <- NULL
  out.d <- NULL
  out.e <- NULL
  out.f <- NULL
  out.g <- NULL
  out.h <-NULL
  out.i <- NULL
  MM <- NULL
  MMb <- NULL
  muuu<- NULL
  
  S <- nrow(alpha)
  k<-length(muu)
  if(is.na(omni.i)) {
    out <- matrix(NA, nrow=N, ncol=2)
    colnames(out) <- c("DomEig", "Im")
    
    
     for(n in 1:N) out[n,] <- {
    
    muuuu<-function(muu){
      muuu<-muu
      muuu[1:6]<-1
      return(muuu)
    }
    mu<-muuuu(muu)
    
    runn<- function (alpha) {
      al <- kk[n,1] * alpha  
      diag(al[1:6,1:6]) <- -1
      diag(al[7:10,7:10]) <- -0.1
      
      return(as.matrix(al))
    }
    alp<-runn(alpha)
    
    linearequations<-function(alp,mu){
      a <- rbind(c(alp[1,1], alp[1,2], alp[1,3], alp[1,4], alp[1,5], alp[1,6], alp[1,7], alp[1,8], alp[1,9], alp[1,10]), 
                 c(alp[2,1], alp[2,2], alp[2,3], alp[2,4], alp[2,5], alp[2,6], alp[2,7], alp[2,8], alp[2,9], alp[2,10]),
                 c(alp[3,1], alp[3,2], alp[3,3], alp[3,4], alp[3,5], alp[3,6], alp[3,7], alp[3,8], alp[3,9], alp[3,10]),
                 c(alp[4,1], alp[4,2], alp[4,3], alp[4,4], alp[4,5], alp[4,6], alp[4,7], alp[4,8], alp[4,9], alp[4,10]),
                 c(alp[5,1], alp[5,2], alp[5,3], alp[5,4], alp[5,5], alp[5,6], alp[5,7], alp[5,8], alp[5,9], alp[5,10]),
                 c(alp[6,1], alp[6,2], alp[6,3], alp[6,4], alp[6,5], alp[6,6], alp[6,7], alp[6,8], alp[6,9], alp[6,10]),
                 c(alp[7,1], alp[7,2], alp[7,3], alp[7,4], alp[7,5], alp[7,6], alp[7,7], alp[7,8], alp[7,9], alp[7,10]),
                 c(alp[8,1], alp[8,2], alp[8,3], alp[8,4], alp[8,5], alp[8,6], alp[8,7], alp[8,8], alp[8,9], alp[8,10]),
                 c(alp[9,1], alp[9,2], alp[9,3], alp[9,4], alp[9,5], alp[9,6], alp[9,7], alp[9,8], alp[9,9], alp[9,10]),
                 c(alp[10,1], alp[10,2], alp[10,3], alp[10,4], alp[10,5], alp[10,6], alp[10,7], alp[10,8], alp[10,9], alp[10,10]))
      b <- c(-mu[1], -mu[2], -mu[3], -mu[4], -mu[5], -mu[6], -mu[7], -mu[8], -mu[9], -mu[10])
      dd<-solve(a, b) 
      return(dd)}
    d<-linearequations(alp,mu)
  
    big<-d[which.min(d)]
    bi<-d[which.max(d)]
    
    
    jaco<- function (mu,alp,d) {
      
      out.a1 <- mu[1]+alp[1,1]*d[1]*2+alp[1,2]*d[2]+alp[1,3]*d[3]+alp[1,4]*d[4]+alp[1,5]*d[5]+alp[1,6]*d[6]+alp[1,7]*d[7]+alp[1,8]*d[8]+alp[1,9]*d[9]+alp[1,10]*d[10]
      out.a2 <- alp[1,2]*d[1]
      out.a3 <- alp[1,3]*d[1]
      out.a4 <- alp[1,4]*d[1]
      out.a5 <- alp[1,5]*d[1]
      out.a6 <- alp[1,6]*d[1]
      out.a7 <- alp[1,7]*d[1]
      out.a8 <- alp[1,8]*d[1]
      out.a9 <- alp[1,9]*d[1]
      out.a10 <- alp[1,10]*d[1]
      
      
      out.b1 <- alp[2,1]*d[2]
      out.b2 <- mu[2]+alp[2,1]*d[1]+alp[2,2]*d[2]*2+alp[2,3]*d[3]+alp[2,4]*d[4]+alp[2,5]*d[5]+alp[2,6]*d[6]+alp[2,7]*d[7]+alp[2,8]*d[8]+alp[2,9]*d[9]+alp[2,10]*d[10]
      out.b3 <- alp[2,3]*d[2]
      out.b4 <- alp[2,4]*d[2]
      out.b5 <- alp[2,5]*d[2]
      out.b6 <- alp[2,6]*d[2]
      out.b7 <- alp[2,7]*d[2]
      out.b8 <- alp[2,8]*d[2]
      out.b9 <- alp[2,9]*d[2]
      out.b10 <- alp[2,10]*d[2]
      
      
      out.c1 <- alp[3,1]*d[3]
      out.c2 <- alp[3,2]*d[3]
      out.c3 <- mu[3]+alp[3,1]*d[1]+alp[3,2]*d[2]+alp[3,3]*d[3]*2+alp[3,4]*d[4]+alp[3,5]*d[5]+alp[3,6]*d[6]+alp[3,7]*d[7]+alp[3,8]*d[8]+alp[3,9]*d[9]+alp[3,10]*d[10]
      out.c4 <- alp[3,4]*d[3]
      out.c5 <- alp[3,5]*d[3]
      out.c6 <- alp[3,6]*d[3]
      out.c7 <- alp[3,7]*d[3]
      out.c8 <- alp[3,8]*d[3]
      out.c9 <- alp[3,9]*d[3]
      out.c10 <- alp[3,10]*d[3]
      
      
      out.d1 <- alp[4,1]*d[4]
      out.d2 <- alp[4,2]*d[4]
      out.d3 <- alp[4,3]*d[4]
      out.d4 <- mu[4]+alp[4,1]*d[1]+alp[4,2]*d[2]+alp[4,3]*d[3]+alp[4,4]*d[4]*2+alp[4,5]*d[5]+alp[4,6]*d[6]+alp[4,7]*d[7]+alp[4,8]*d[8]+alp[4,9]*d[9]+alp[4,10]*d[10]
      out.d5 <- alp[4,5]*d[4]
      out.d6 <- alp[4,6]*d[4]
      out.d7 <- alp[4,7]*d[4]
      out.d8 <- alp[4,8]*d[4]
      out.d9 <- alp[4,9]*d[4]
      out.d10 <- alp[4,10]*d[4]
      
      
      out.e1 <- alp[5,1]*d[4]
      out.e2 <- alp[5,2]*d[5]
      out.e3 <- alp[5,3]*d[5]
      out.e4 <- alp[5,4]*d[5]
      out.e5 <- mu[5]+alp[5,1]*d[1]+alp[5,2]*d[2]+alp[5,3]*d[3]+alp[5,4]*d[4]+alp[5,5]*d[5]*2+alp[5,6]*d[6]+alp[5,7]*d[7]+alp[5,8]*d[8]+alp[5,9]*d[9]+alp[5,10]*d[10]
      out.e6 <- alp[5,6]*d[5]
      out.e7 <- alp[5,7]*d[5]
      out.e8 <- alp[5,8]*d[5]
      out.e9 <- alp[5,9]*d[5]
      out.e10 <- alp[5,10]*d[5]
      
      out.f1 <- alp[6,1]*d[4]
      out.f2 <- alp[6,2]*d[6]
      out.f3 <- alp[6,3]*d[6]
      out.f4 <- alp[6,4]*d[6]
      out.f5 <- alp[6,5]*d[6]
      out.f6 <- mu[6]+alp[6,1]*d[1]+alp[6,2]*d[2]+alp[6,3]*d[3]+alp[6,4]*d[4]+alp[6,5]*d[5]+alp[6,6]*d[6]*2+alp[6,7]*d[7]+alp[6,8]*d[8]+alp[6,9]*d[9]+alp[6,10]*d[10]
      out.f7 <- alp[6,7]*d[6]
      out.f8 <- alp[6,8]*d[6]
      out.f9 <- alp[6,9]*d[6]
      out.f10 <- alp[6,10]*d[6]
      
      
      out.g1 <- alp[7,1]*d[7]
      out.g2 <- alp[7,2]*d[7]
      out.g3 <- alp[7,3]*d[7]
      out.g4 <- alp[7,4]*d[7]
      out.g5 <- alp[7,5]*d[7]
      out.g6 <- alp[7,6]*d[7]
      out.g7 <- mu[7]+alp[7,1]*d[1]+alp[7,2]*d[2]+alp[7,3]*d[3]+alp[7,4]*d[4]+alp[7,5]*d[5]+alp[7,6]*d[6]+alp[7,7]*d[7]*2+alp[7,8]*d[8]+alp[7,9]*d[9]+alp[7,10]*d[10]
      out.g8 <- alp[7,8]*d[7]
      out.g9 <- alp[7,9]*d[7]
      out.g10 <- alp[7,10]*d[7]
      
      
      out.h1 <- alp[8,1]*d[8]
      out.h2 <- alp[8,2]*d[8]
      out.h3 <- alp[8,3]*d[8]
      out.h4 <- alp[8,4]*d[8]
      out.h5 <- alp[8,5]*d[8]
      out.h6 <- alp[8,6]*d[8]
      out.h7 <- alp[8,7]*d[8]
      out.h8 <- mu[8]+alp[8,1]*d[1]+alp[8,2]*d[2]+alp[8,3]*d[3]+alp[8,4]*d[4]+alp[8,5]*d[5]+alp[8,6]*d[6]+alp[8,7]*d[7]+alp[8,8]*d[8]*2+alp[8,9]*d[9]+alp[8,10]*d[10]
      out.h9 <- alp[8,9]*d[8]
      out.h10 <- alp[8,10]*d[8]
      
      
      out.i1 <- alp[9,1]*d[9]
      out.i2 <- alp[9,2]*d[9]
      out.i3 <- alp[9,3]*d[9]
      out.i4 <- alp[9,4]*d[9]
      out.i5 <- alp[9,5]*d[9]
      out.i6 <- alp[9,6]*d[9]
      out.i7 <- alp[9,7]*d[9]
      out.i8 <- alp[9,8]*d[9]
      out.i9 <- mu[9]+alp[9,1]*d[1]+alp[9,2]*d[2]+alp[9,3]*d[3]+alp[9,4]*d[4]+alp[9,5]*d[5]+alp[9,6]*d[6]+alp[9,7]*d[7]+alp[9,8]*d[8]+alp[9,9]*d[9]*2+alp[9,10]*d[10]
      out.i10 <- alp[9,10]*d[9]
      
      out.j1 <- alp[10,1]*d[10]
      out.j2 <- alp[10,2]*d[10]
      out.j3 <- alp[10,3]*d[10]
      out.j4 <- alp[10,4]*d[10]
      out.j5 <- alp[10,5]*d[10]
      out.j6 <- alp[10,6]*d[10]
      out.j7 <- alp[10,7]*d[10]
      out.j8 <- alp[10,8]*d[10]
      out.j9 <- alp[10,9]*d[10]
      out.j10 <- mu[10]+alp[10,1]*d[1]+alp[10,2]*d[2]+alp[10,3]*d[3]+alp[10,4]*d[4]+alp[10,5]*d[5]+alp[10,6]*d[6]+alp[10,7]*d[7]+alp[10,8]*d[8]+alp[10,9]*d[9]+alp[10,10]*d[10]*2
      
      MM<-rbind(c(out.a1,out.a2, out.a3, out.a4, out.a5, out.a6, out.a7, out.a8, out.a9, out.a10),
                c(out.b1,out.b2, out.b3, out.b4, out.b5, out.b6, out.b7, out.b8, out.b9, out.b10),
                c(out.c1,out.c2, out.c3, out.c4, out.c5, out.c6, out.c7, out.c8, out.c9, out.c10),
                c(out.d1,out.d2, out.d3, out.d4, out.d5, out.d6, out.d7, out.d8, out.d9, out.d10),
                c(out.e1,out.e2, out.e3, out.e4, out.e5, out.e6, out.e7, out.e8, out.e9, out.e10),
                c(out.f1,out.f2, out.f3, out.f4, out.f5, out.f6, out.f7, out.f8, out.f9, out.f10),
                c(out.g1,out.g2, out.g3, out.g4, out.g5, out.g6, out.g7, out.g8, out.g9, out.g10),
                c(out.h1,out.h2, out.h3, out.h4, out.h5, out.h6, out.h7, out.h8, out.h9, out.h10),
                c(out.i1,out.i2, out.i3, out.i4, out.i5, out.i6, out.i7, out.i8, out.i9, out.i10),
                c(out.j1,out.j2, out.j3, out.j4, out.j5, out.j6, out.j7, out.j8, out.j9, out.j10))
      return(as.matrix(MM))
    }
    
    MMb<-jaco(mu,alp,d)
    eigs <- eigen(MMb)[["values"]]
    mx <- which.max(Re(eigs))
    deM <- Re(eigs)[mx]
    
    deM<-ifelse(big<0 |bi ==Inf, NA,  deM) 
    
    deMi <- Im(eigs)[mx]
    c(deM, deMi)
  } } else { print("exist omnivores")}
return(as.data.frame(out))
}

out.E <- law(alpha, N = 10000)
summary(out.E)
RT.E<- -1/(out.E[["DomEig"]]) 
summary(RT.E)


