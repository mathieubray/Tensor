library(dplyr)
library(rTensor)

# Xu, Yin 2013. Algorithm 2 (p. 17)

# Input: Non-negative N-way tensor Y
# Output: Component matrices A, objective values at each iteration
BCD <- function(Y, A, delta, tol=1e-5, max_k=10000){
  
  N <- Y@num_modes
  
  A.prev <- A
  t.prev <- 1
  L.prev <- rep(1,N)
  L <- rep(1,N)
  
  k <- 0
  diff <- 1
  prev.diff <- 1
  objs <- rep(0,max_k)
  
  #obj.prev <- 0.5 * norm(k_unfold(Y,1)@data - A[[1]]%*%t(khatri_rao_list(A[-1],reverse=TRUE)), "2")^2
  obj.prev <- norm(k_unfold(Y,1)@data - A[[1]]%*%t(khatri_rao_list(A[-1],reverse=TRUE)), "F")
  
  while (k < max_k & diff > tol & prev.diff > tol){
    
    k <- k + 1
    
    t <- 0.5 * (1 + sqrt(1 + 4*t.prev^2))
    w.hat <- (t.prev-1)/t
    
    for (n in 1:N){
      
      B <- khatri_rao_list(A[-n],reverse=TRUE)
      
      L[n] <- norm(t(B)%*%B, "2") # Spectral Norm
      
      w <- ifelse(k==1,
                  min(w.hat, delta),
                  min(w.hat, delta * sqrt(L.prev[n]/L[n])))
      
      A.hat <- A[[n]] + w*(A[[n]] - A.prev[[n]])
      
      G <- (A.hat%*%t(B) - k_unfold(Y,n)@data)%*%B
      
      A[[n]] <- pmax(A.hat - G/L[n], 0)
      
    }
    
    #obj <- 0.5 * norm(k_unfold(Y,1)@data - A[[1]]%*%t(khatri_rao_list(A[-1],reverse=TRUE)), "2")^2
    obj <- norm(k_unfold(Y,1)@data - A[[1]]%*%t(khatri_rao_list(A[-1],reverse=TRUE)), "F")
    
    if (obj >= obj.prev){
      
      for (n in 1:N){
        
        B <- khatri_rao_list(A[-n],reverse=TRUE)
        
        G <- (A.prev[[n]]%*%t(B) - k_unfold(Y,n)@data)%*%B
        
        A[[n]] <- pmax(A.prev[[n]] - G/L[n], 0)
        
      }
    }
    
    #diff <- 0.5 * norm(k_unfold(Y,1)@data - A[[1]]%*%t(khatri_rao_list(A[-1],reverse=TRUE)), "2")^2
    diff <- norm(k_unfold(Y,1)@data - A[[1]]%*%t(khatri_rao_list(A[-1],reverse=TRUE)), "F")
    
    objs[k] <- diff
    
    if (k >= 2){
      prev.diff <- abs(objs[k] - objs[k-1])
    }
    
    A.prev <- A
    t.prev <- t
    obj.prev <- diff
   
  }
  
  return(list(A=A, objs=objs[1:k]))
  
}





