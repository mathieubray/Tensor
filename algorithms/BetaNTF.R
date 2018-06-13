library(dplyr)
library(purrr)
library(rTensor)

# Beta NTF. Algorithm 2. Flatz 2013

# Input: Non-negative N-way tensor Y
# Output: Component matrices A, objective values at each iteration
betaNTF <- function(Y, A, beta, tol=1e-5, max_k=10000){
  
  N <- Y@num_modes
  
  normalize.component <- function(Q){
    
    R <- t(rep(1,nrow(Q))) %*% Q %>%
      as.vector %>%
      diag %>%
      solve
    
    return(Q %*% R)
  }
  
  k <- 0
  diff <- 1
  prev.diff <- 1
  objs <- rep(0,max_k)
  
  while (k < max_k & diff > tol & prev.diff > tol){
    
    k <- k + 1
    
    Y.hat <- A
    
    for(n in 1:N){
      
      Y.hat.n <- Y.hat[[n]] %*% t(khatri_rao_list(Y.hat[-n],reverse=TRUE))
      kr.A <- khatri_rao_list(A[-n],reverse=TRUE)
      
      A[[n]] <- A[[n]] * ((k_unfold(Y,n)@data / (Y.hat.n)^(beta-1)) %*% kr.A) / ((Y.hat.n)^(beta) %*% kr.A)
      
      if (n != N){
        A[[n]] <- A[[n]] %>% normalize.component
      }
      
    }
    
    diff <- norm(k_unfold(Y,1)@data - A[[1]] %*% t(khatri_rao_list(A[-1],reverse=TRUE)),"F")
    
    objs[k] <- diff
    
    if (k >= 2){
      prev.diff <- abs(objs[k] - objs[k-1])
    }
    
  }
  
  return(list(A=A, objs=objs[1:k]))
  
}


