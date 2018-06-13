library(rTensor)
library(tidyverse)

# Zhou et al. Acceleratin Online CP Decomposition

# Algorithm 1

getKhatriRaoList <- function(A){
  
  N <- length(A)
  
  left <- list(A[[N-1]])
  right <- list(A[[1]])
  
  if (N > 3){
    
    for (n in 2:(N-2)){
      left[[n]] <- khatri_rao(left[[n-1]], A[[N-n]])
      right[[n]] <- khatri_rao(A[[n]], right[[n-1]])
    }
    
  }
  
  K <- list(left[[N-2]])
  K[[N-1]] <- right[[N-2]]
  
  if (N > 3){
    
    for (n in 2:(N-2)){
      K[[n]] <- khatri_rao(left[[N-n-1]], right[[n-1]])
    }
    
  }
  
  return(K)
  
}


# Algorithm 2

initializeOnlineCP <- function(initX, R, A=NULL){
  
  if(is.null(A)){
    
    estInitX <- cp(initX, R, tol = 1e-8)
    
    A <- estInitX$U
    
    A[[length(A)]] <- A[[length(A)]] %*% diag(estInitX$lambdas)
    
  }
  
  N <- initX@num_modes
  
  H <- hadamard_list(map(A, function(x){ return(t(x) %*% x) }))
  
  K <- getKhatriRaoList(A)
  
  P <- list()
  Q <- list()
  
  for (n in 1:(N-1)){
    
    X <- k_unfold(initX, n)@data
    P[[n]] <- X %*% khatri_rao(A[[N]], K[[n]])
    Q[[n]] <- H / (t(A[[n]]) %*% A[[n]])
    
  }
  
  return(list(A=A,P=P,Q=Q))
  
}


initX <- as.tensor(array(rnorm(60),dim=c(4,5,3)))

initPQ <- initializeOnlineCP(initX, 2)

A <- initPQ$A
P <- initPQ$P
Q <- initPQ$Q

newX <- as.tensor(array(rnorm(40),dim=c(4,5,2)))

onlineCP <- function(newX, A, P, Q){
  
  N <- length(A)
  R <- ncol(A[[1]])
  
  dims <- newX@modes
  
  batchsize <- dims %>% tail(1)
  
  K <- getKhatriRaoList(A)
  
  H <- hadamard_list(map(A[-N], function(x){ return(t(x) %*% x) }))
  
  # Update Mode-N
  
  KN <- khatri_rao(K[[1]],A[[1]])
  AN <- (k_unfold(newX,N)@data %*% KN) / H
  
  A[[N]] <- rbind(A[[N]], AN)
  
  for (n in 1:(N-1)){
    newXn <- k_unfold(newX,n)@data
    P[[n]] <- P[[n]] + newXn %*% khatri_rao(AN, K[[n]])
    Hn = H / (t(A[[n]]) %*% A[[n]])
    Q[[n]] <- Q[[n]] + ((t(AN)%*%AN) %*% Hn)
    A[[n]] <- P[[n]] %*% solve(Q[[n]])
  }
  
  return(list(A=A,P=P,Q=Q))
  
}

J <- onlineCP(newX, A, P, Q)
