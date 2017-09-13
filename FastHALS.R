library(dplyr)
library(rTensor)

# FAST HALS NTF. Algorithm 4. Flatz 2013

# Input: Non-negative N-way tensor Y
# Output: Component matrices A, objective values at each iteration
fastHALS <- function(Y, A, tol=1e-5, max_k=10000){
  
  N <- Y@num_modes
  
  J <- ncol(A[[1]])
  
  I <- array(rep(0,J^3),dim=c(J,J,J)) %>% as.tensor
  for (j in 1:J){
    I[j,j,j] <- 1
  }
  
  normalize.matrix.columns <- function(Q){
    return(apply(Q, 2, function(x){ x/sqrt(sum(x^2)) }))
  }
  
  A.l <- A %>% map(normalize.matrix.columns)
  A[1:(N-1)] <- A.l[1:(N-1)]
  
  k <- 0
  diff <- 1
  prev.diff <- 1
  objs <- rep(0,max_k)
  
  T.1 <-  map(A, function(X) { return(t(X) %*% X) }) %>% hadamard_list
  
  while (k <= max_k & diff > tol & prev.diff > tol){
    
    k <- k+1
    
    gamma <- t(A[[N]]) %*% A[[N]] %>% diag
    
    for(n in 1:N){
      if (n == N){
        gamma <- rep(1,J)
      }
      
      T.2 <- k_unfold(Y,n)@data %*% khatri_rao_list(A[-n],reverse=TRUE)
      
      T.3 <- T.1 / (t(A[[n]]) %*% A[[n]])
      
      for(j in 1:J){
        
        b <- A[[n]][,j] * gamma[j] + T.2[,j] - A[[n]] %*% T.3[,j]
        
        b <- ifelse(b>0,b,0)
        
        if (n != N){
          b <- b %>% normalize.matrix.columns
        }
        
        A[[n]][,j] <- b
      }
      
      T.1 <- T.3 *(t(A[[n]]) %*% A[[n]])
      
    }
    
    Y.hat <- I %>% ttl(A,ms=c(1,2,3))
    
    diff <- fnorm(Y - Y.hat)
    
    objs[k] <- diff
    
    if (k >= 2){
      prev.diff <- abs(objs[k] - objs[k-1])
    }
  }
  
  return(list(A=A, objs=objs[1:k]))
}





