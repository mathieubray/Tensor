
set.seed(7900)


#Y <- runif(64,0,10) %>%
# array(dim=c(4,4,4)) %>%
#  as.tensor()



# Input: Non-negative N-way tensor M and rank r
Y.components <- list(matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2))

Y <- I %>% ttl(Y.components,ms=c(1,2,3))



# Random Initialization of A
A <- list(matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2)) # Last component not normalized for some reason

J <- 2






results <- alphaNTF(Y,A,J)

B <- results$A

X.1 <- as.vector(B[[1]][,1]) %o% as.vector(B[[2]][,1]) %o% as.vector(B[[3]][,1]) %>% as.tensor
X.2 <- as.vector(B[[1]][,2]) %o% as.vector(B[[2]][,2]) %o% as.vector(B[[3]][,2]) %>% as.tensor

X <- X.1 + X.2

X@data

X - Y -> l


l@data


library(ggplot2)

norm.frame <- data.frame(Time=1:length(results$objs),Norm=results$objs) %>%
  filter(Time > 100)

ggplot(data=norm.frame, aes(x=Time,y=Norm)) +
  geom_path() +
  theme_bw()




results <- fastHALS(Y,A,J)


B <- results$A

X.1 <- as.vector(B[[1]][,1]) %o% as.vector(B[[2]][,1]) %o% as.vector(B[[3]][,1]) %>% as.tensor
X.2 <- as.vector(B[[1]][,2]) %o% as.vector(B[[2]][,2]) %o% as.vector(B[[3]][,2]) %>% as.tensor

X <- X.1 + X.2

X@data



norm.frame <- data.frame(Time=1:length(results$objs),Norm=results$objs) %>%
  filter(Time > 5)

ggplot(data=norm.frame,aes(x=Time,y=Norm)) +
  geom_path() +
  theme_bw()


set.seed(7932)

I <- array(c(1,0,0,0,0,0,0,1),dim=c(2,2,2)) %>%
  as.tensor

# Input: Non-negative N-way tensor M and rank r
Y.components <- list(matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2))

Y <- I %>% ttl(Y.components,ms=c(1,2,3)) + rand_tensor(modes=c(4,4,4))


A <- list(matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2)) # Last component not normalized for some reason



J <- 2







set.seed(7900)


#Y <- runif(64,0,10) %>%
# array(dim=c(4,4,4)) %>%
#  as.tensor()


I <- array(c(1,0,0,0,0,0,0,1),dim=c(2,2,2)) %>%
  as.tensor

# Input: Non-negative N-way tensor M and rank r
Y.components <- list(matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2))

Y <- I %>% ttl(Y.components,ms=c(1,2,3))



# Random Initialization of A
A <- list(matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2)) # Last component not normalized for some reason

J <- 2





# Compare to General Tensor Decomposition Case
# Kalman Filter, Online Algorithm for new matrix slices
# 






















outer_product <- function(A,r){
  
  len <- length(A)
  
  d <- purrr::map_int(A,nrow)
  
  prod.tensor <- rep(0,prod(d)) %>%
    array(dim=d) %>%
    as.tensor()
  
  for (i in 1:r){
    
  }
  
}

list.to.tensor <- function(Q){
  
  i <- nrow(Q[[1]])
  j <- ncol(Q[[1]])
  k <- length(Q)
  
  R <- Q %>%
    unlist %>%
    array(dim=c(i,j,k)) %>%
    as.tensor
  
  return(R)
}

X <- array(1:24,dim=c(3,4,2)) %>% as.tensor

U <- matrix(1:6,nrow=2)

Y <- X %>% ttm(U,m=1)

k_unfold(Y,1) == U%*%k_unfold(X,1)@data



A <- matrix(1:8, nrow=4)

B <- matrix(1:6, nrow=3)

C <- matrix(1:4, nrow=2)

X.1 <- as.vector(A[,1]) %o% as.vector(B[,1]) %o% as.vector(C[,1]) %>% as.tensor
X.2 <- as.vector(A[,2]) %o% as.vector(B[,2]) %o% as.vector(C[,2]) %>% as.tensor

X <- X.1 + X.2

k_unfold(X,1)

A %*% t(khatri_rao(C,B))

k_unfold(X,2)

B %*% t(khatri_rao(C,A))

k_unfold(X,3)

C %*% t(khatri_rao(B,A))

D <- matrix(8:1, nrow=4)

lambda <- c(2,0.5)
X.1 <- lambda[1] * as.vector(A[,1]) %o% as.vector(B[,1]) %o% as.vector(C[,1]) %o% as.vector(D[,1]) %>% as.tensor
X.2 <- lambda[2] * as.vector(A[,2]) %o% as.vector(B[,2]) %o% as.vector(C[,2]) %o% as.vector(D[,2]) %>% as.tensor

X = X.1 + X.2 

comp.list <- list(A,B,C,D)



diag(lambda)

k_unfold(X,1) == comp.list[[1]] %*% diag(lambda) %*% t(khatri_rao_list(comp.list[-1],reverse=TRUE))



node.frame <- rbind(data.frame(U[[1]],nodetype=as.factor(nodes),u="Donor"),data.frame(U[[2]],nodetype=as.factor(nodes),u="Candidate"))


node.frame <- cbind(node.frame,node = c(1:100,1:100))

ggplot(data=node.frame, aes(x=X1,y=X2,color=nodetype)) +
  facet_wrap(~u) +
  geom_point(size=4,alpha=0.5) +
  theme_bw()

K <- data.frame(Time=1:25,Component=U[[3]]) %>%
  gather(key=Component,value=Value,Component.1,Component.2)


ggplot(data=K,aes(x=Time,y=Value,color=Component)) + geom_line() + theme_bw()


####


ggplot(data=node.frame, aes(x=X1,y=X2,color=nodetype)) +
  facet_wrap(~u) +
  geom_point(size=4,alpha=0.5) +
  theme_bw()

K <- data.frame(Time=1:25,Component=B[[3]]) %>%
  gather(key=Component,value=Value,Component.1,Component.2)


ggplot(data=K,aes(x=Time,y=Value,color=Component)) + geom_line() + theme_bw()



B <- decomp.hals$A
B <- decomp.alpha$A









































library(igraph)
library(visNetwork)
library(dplyr)
library(RColorBrewer)

set.seed(900707)

block.matrix <- matrix(nrow=3,ncol=3,0.05)
diag(block.matrix) <- c(0.8,0.8,0.8)
block.sizes <- c(10,10,10)

true.label <- c(rep(1,10),rep(2,10),rep(3,10))

graph <- sample_sbm(sum(block.sizes),block.matrix,block.sizes,directed=F,loops=F)


plot(graph)


clp  <- cluster_label_prop(graph)

clp


edges <- igraph::as_data_frame(graph, what="edges")
vertices <- data.frame(id=1:sum(block.sizes),true=true.label) %>%
  mutate(label = id)

colors <- brewer.pal(3,"Set1")
true.colors <- brewer.pal(3,"Set2")


vertices$color.background <- true.colors[vertices$true]

vertices$color.border <- colors[clp$membership]

visNetwork(vertices,edges,main="Generated Network") %>%
  visNodes(shadow = T,borderWidth = 5)


set.seed(900707)

G <- as_adjacency_matrix(graph,sparse=F)

X <- matrix(runif(90),nrow=30,ncol=3)

numIts <- 1
tol <- 1
while(tol > 0.00001 & numIts < 1000){
  
  GX <- G%*%X
  XXTX <- X%*%t(X)%*%X
  
  Xnew <- X * (0.5 + GX/(2*XXTX))
  
  tol <- sqrt(sum(t(Xnew - X)%*%(Xnew - X)))
  
  X <- Xnew
  
  numIts <- numIts + 1
  
}

X <- X/rowSums(X)

vertices$color.border <- colors[max.col(X)]

visNetwork(vertices,edges,main="Generated Network") %>%
  visNodes(shadow = T,borderWidth = 5)





###


set.seed(900707)

block.matrix <- matrix(nrow=3,ncol=3,0.05)
diag(block.matrix) <- c(0.8,0.8,0.8)
block.sizes <- c(10,10,10)

true.label <- c(rep(1,10),rep(2,10),rep(3,10))

graph <- sample_sbm(sum(block.sizes),block.matrix,block.sizes,directed=T,loops=F)

plot(graph)


clp  <- cluster_label_prop(graph)
clp

edges <- igraph::as_data_frame(graph, what="edges")
vertices <- data.frame(id=1:sum(block.sizes),true=true.label) %>%
  mutate(label = id)

colors <- brewer.pal(3,"Set1")
true.colors <- brewer.pal(3,"Set2")


vertices$color.background <- true.colors[vertices$true]

vertices$color.border <- colors

edges$arrows="to"

visNetwork(vertices,edges,main="Generated Network") %>%
  visNodes(shadow = T,borderWidth = 5)


set.seed(900707)

A <- as_adjacency_matrix(graph,sparse=F)

X <- matrix(runif(90),nrow=30,ncol=3)
S <- matrix(runif(9),nrow=3,ncol=3)

numIts <- 1
tol <- 1
while(tol > 0.00001 & numIts < 1000){
  
  ATXS <- t(A)%*%X%*%S
  AXST <- A%*%X%*%t(S)
  XSXTXST <- X%*%S%*%t(X)%*%X%*%t(S)
  XSTXTXS <- X%*%t(S)%*%t(X)%*%X%*%S
  XTAX <- t(X)%*%A%*%X
  XTXSXTX <- t(X)%*%X%*%S%*%t(X)%*%X
  
  
  Xnew <- X * ((ATXS + AXST)/(XSXTXST + XSTXTXS))^(0.25)
  Snew <- S * (XTAX/XTXSXTX)^(0.25)
  
  tol <- sqrt(sum(t(Xnew - X)%*%(Xnew - X))) +(sum(t(Snew - S)%*%(Snew - S)))
  
  X <- Xnew
  S <- Snew
  
  numIts <- numIts + 1
  
}

X <- X/rowSums(X)

vertices$color.border <- colors[max.col(X)]

visNetwork(vertices,edges,main="Generated Network") %>%
  visNodes(shadow = T,borderWidth = 5)







###################

# Tensor Networks 
# 





















set.seed(900)

#M <- runif(64,0,10) %>%
# array(dim=c(4,4,4)) %>%
#as.tensor()


I <- array(c(1,0,0,0,0,0,0,1),dim=c(2,2,2)) %>%
  as.tensor

# Input: Non-negative N-way tensor M and rank r
Y.components <- list(matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2),
                     matrix(runif(16,0,3),nrow=4,ncol=2))

M <- I %>% ttl(Y.components,ms=c(1,2,3)) + rand_tensor(modes=c(4,4,4))


# Initialization: postive # delta < 1, random A non-negative
A <- list(matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2),
          matrix(runif(16,0,10),nrow=4,ncol=2))


r <- 2

results <- BCD(M,A,r)


B <- results$A

X.1 <- as.vector(B[[1]][,1]) %o% as.vector(B[[2]][,1]) %o% as.vector(B[[3]][,1]) %>% as.tensor
X.2 <- as.vector(B[[1]][,2]) %o% as.vector(B[[2]][,2]) %o% as.vector(B[[3]][,2]) %>% as.tensor

X <- X.1 + X.2

X@data



norm.frame <- data.frame(Time=1:length(results$objs),Norm=results$objs) %>%
  filter(Time > 100)

ggplot(data=norm.frame,aes(x=Time,y=Norm)) +
  geom_path() +
  theme_bw()







































per_mr <- 4 # Number of Nodes Arriving between each Match Run
num_mrs <- 25 # Number of Match Runs
num_nodes <- per_mr * num_mrs # Number of Total Nodes

prob_type <- c(0.4,0.1,0.1,0.4)
prob_connect <- matrix(c(0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05,0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05),nrow=4)

node_types <- sample(1:length(prob_type),num_nodes,replace=T,prob=prob_type) # Randomly Assign Type to Each Node

crossmatch_prob <- c(0.9,0.75,0.75,0.5) # Probability of Connection Remaining After Evaluation


alphaNTF <- function(Y, A, J, alpha=2, tol=1e-5, max_k=10000){
  
  N <- Y@num_modes
  
  normalize.component <- function(Q){
    
    R <- t(rep(1,nrow(Q))) %*% Q %>%
      as.vector %>%
      diag %>%
      solve
    
    return(Q %*% R)
  }
  
  # Normalize to unit length
  A.l <- A %>% map(normalize.component)
  A[1:(N-1)] <- A.l[1:(N-1)]
  
  k <- 0
  diff <- 1
  prev.diff <- 1
  objs <- rep(0,max_k)
  
  while (k < max_k & diff > tol & prev.diff > tol){
    
    k <- k + 1
    
    Y.hat <- A
    
    for(n in 1:N){
      
      Y.hat.n <- Y.hat[[n]] %*% t(khatri_rao_list(Y.hat[-n],reverse=TRUE))
      
      A[[n]] <- A[[n]] * ((k_unfold(Y,n)@data / Y.hat.n)^alpha %*% khatri_rao_list(A.l[-n],reverse=TRUE))^(1/alpha)
      
      A.l[[n]] <- A[[n]] %>% normalize.component
      
      if (n != N){
        A[[n]] <- A.l[[n]]
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

fastHALS <- function(Y, A, J, tol=1e-5, max_k=10000){
  
  N <- Y@num_modes
  
  I <- array(rep(0,J^3),dim=c(J,J,J)) %>% as.tensor
  for (j in 1:J){
    I[j,j,j] <- 1
  }
  
  normalize.matrix.columns <- function(Q){
    return(apply(Q, 2, function(x){ x/sqrt(sum(x^2)) }))
  }
  
  #Normalize
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



