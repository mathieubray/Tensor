library(rTensor)

generate.kpd <- function(per_mr, num_mrs, prob_type, prob_connect, xm_prob){
  
  num_nodes <- per_mr * num_mrs # Number of Total Nodes
  node_types <- sample(1:length(prob_type),num_nodes,replace=T,prob=prob_type) # Randomly Assign Type to Each Node
  
  matrix.list <- array(0,dim=c(num_nodes,num_nodes,num_mrs)) # Initialize a Matrix List
  
  old.matrix <- matrix(0,nrow=num_nodes,ncol=num_nodes) # Initialize First Matrix
  
  for (it in 1:num_mrs){
    
    new.matrix <- old.matrix
    
    # Renege
    
    if (it > 1){
      for (i in 1:(per_mr*(it-1))){
        for (j in 1:(per_mr*(it-1))){
          if (new.matrix[i,j] == 1){
            
            q <- rbinom(1,1,xm_prob[node_types[j]])
            
            new.matrix[i,j] <- q
          }
        }
      }
    }
    
    # New Connections
    
    for (i in 1:(per_mr*it)){
      for (j in 1:(per_mr*it)){
        if (i != j & !(i <= per_mr*(it-1) & j <= per_mr*(it-1))){
          q <- rbinom(1,1,prob_connect[node_types[i],node_types[j]])
          
          new.matrix[i,j] <- q 
        }
      }
    }
    
    matrix.list[,,it] <- new.matrix
    
    old.matrix <- new.matrix
    
  }
  
  kpd <- as.tensor(matrix.list)
  
  return(list(kpd=kpd,node_types=node_types))
}


initialize.matrix <- function(per_mr, num_mrs, num_components){
  
  num_nodes <- per_mr * num_mrs
  
  init.matrix <- list(matrix(runif(num_components * num_nodes,0,1),nrow=num_nodes,ncol=num_components),
                      matrix(runif(num_components * num_nodes,0,1),nrow=num_nodes,ncol=num_components),
                      matrix(runif(num_components * num_mrs,0,1),nrow=num_mrs,ncol=num_components))
  
  return(init.matrix)
}