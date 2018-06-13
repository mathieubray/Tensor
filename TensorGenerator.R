library(rTensor)
library(tidyverse)
library(igraph)

library(ompr) # Mixed integer programming
library(ompr.roi)
library(ROI) # R optimization interface
library(ROI.plugin.glpk) # Will use solver 'glpk'



initialize.nodes <- function(nodes.per.matchrun, number.of.matchruns, true.cluster.prob){
  
  number.of.nodes <- nodes.per.matchrun * number.of.matchruns # Number of Total Nodes
  
  node.types <- sample(1:length(true.cluster.prob), number.of.nodes, replace=T, prob=true.cluster.prob) # Randomly Assign Type to Each Node
  
  node.joined <- rep(1:number.of.matchruns, each=nodes.per.matchrun)
  
  node.list <- tibble(Node = 1:number.of.nodes,
                      Type = node.types,
                      MatchRun = node.joined,
                      Status = "Joined")
  
  return(node.list)
}

virtual.crossmatch <- function(match.run, node.info, cluster.connect.prob){
  
  
  simulate.match <- function(donor.type, candidate.type, prob.connect){
    
    match <- rbinom(1, 1, prob.connect[donor.type, candidate.type])
    
    return(match)
  } 
  
  new.nodes <- node.info %>%
    filter(MatchRun == match.run, Status == "Joined") %>%
    .$Node
  
  available.nodes <- node.info %>%
    filter(MatchRun <= match.run) %>%
    group_by(Node,Type) %>%
    summarize(Status = last(Status)) %>%
    ungroup %>%
    filter(Status != "Transplanted") %>%
    select(-Status)
  
  new.matches <- tidyr::crossing(Donor = available.nodes$Node, Candidate = available.nodes$Node) %>%
    left_join(available.nodes, by=c("Donor"="Node")) %>%
    rename(DonorType = Type) %>%
    left_join(available.nodes, by=c("Candidate"="Node")) %>%
    rename(CandidateType = Type) %>%
    filter(Donor != Candidate,
           Donor %in% new.nodes | Candidate %in% new.nodes) %>%
    rowwise %>%
    mutate(Match = simulate.match(DonorType, CandidateType, cluster.connect.prob)) %>%
    ungroup %>%
    filter(Match == 1) %>%
    select(Donor, Candidate) %>%
    mutate(MatchRun = match.run,
           Status = "Virtual Crossmatch")
  
  return(new.matches)
  
}

get.adjacency.matrix <- function(match.info, number.of.nodes){
  
  matches <- match.info %>%
    group_by(Donor, Candidate) %>%
    summarize(Status = last(Status)) %>%
    ungroup %>%
    filter(Status %in% c("Virtual Crossmatch","Successful Lab Crossmatch")) 
  
  adj.matrix <- matrix(0, nrow=number.of.nodes, ncol=number.of.nodes)
  
  for (m in 1:nrow(matches)){
    adj.matrix[matches$Donor[m], matches$Candidate[m]] <- 1
  }
  
  return(adj.matrix)
}

collect.cycles <- function(incidence.matrix, nodes, max.cycle.size){
  
  cycle_list <- list() # Collect found cycles
  cycles <- 0 # Count cycles
  
  # Retrieve first unvisited child of 'current' node (after 'from' node)
  # Returns -1 if no child found
  get_child <- function(from, current, visited){
    
    if (from != nodes){
      for (j in (from + 1):nodes){
        if (incidence.matrix[current,j] & !visited[j]){
          return(j)
        }
      }
    }
    
    return(-1)
  }
  
  # Use stack to store and evaluate current set of nodes
  node_stack <- numeric()
  
  # Keep note of the first node in the stack (head node)
  current_head_node <- 1
  
  # Keep track of whether each node has been explored already
  visited <- rep(FALSE, nodes)
  
  # Iterate over all nodes
  while (current_head_node <= nodes){
    
    # Place new head node in the stack and mark as visited
    node_stack <- c(node_stack, current_head_node)
    visited[current_head_node] <- TRUE 
    
    # Get the first child of this node
    v <- get_child(current_head_node, current_head_node, visited)
    
    while (length(node_stack) > 0){
      
      # If there are no children to explore
      if (v == -1){
        
        # Pop node from stack
        backtrack_node <- node_stack %>% tail(1)
        node_stack <- node_stack %>% head(-1)
        
        if (backtrack_node == current_head_node){
          # We've popped the head node from the stack, which is now empty
          # Break from the loop and start a new stack with a new head node
          visited[backtrack_node] <- FALSE
          break
        }
        
        # Mark the popped node as not visited
        visited[backtrack_node] <- FALSE
        
        # Get new node from top of stack and get its next unvisited child
        new_top_node <- node_stack %>% tail(1)
        v <- get_child(backtrack_node, new_top_node, visited)
        
        
      } else {
        
        # Place new node on top of stack and mark as visited
        node_stack <- c(node_stack, v)
        visited[v] <- TRUE
        
        # If this node connects back to the head node, we've found a cycle!
        if (incidence.matrix[v,current_head_node]){
          
          # Check size constraints (may be unnecessary...)
          if (length(node_stack) <= max.cycle.size){
            
            # Add cycle to list
            cycles <- cycles + 1
            cycle_list[[cycles]] <- node_stack
          }
          
        }
        
        if (length(node_stack) >= max.cycle.size){
          # If stack has reached the maximum cycle size, signal that there are no
          # more children to explore
          v <- -1
        } else {
          # Else, grow stack with child of current node
          v <- get_child(current_head_node, v, visited)
        }
      }
    }
    
    # Advance head node 
    current_head_node <- current_head_node + 1
    
  }
  
  return(cycle_list)
}

select.cycles <- function(cycles){
  
  if (length(cycles) > 1){
    
    cycle_contains_node <- function(cycle,node){
      map2_lgl(cycle,node, ~ .y %in% cycles[[.x]]) # A little confusing, but .x = 'cycle', .y = 'node' here
    }
    
    utility <- map_int(cycles,length)
    number_of_cycles <- cycles %>% length
    number_of_nodes <- cycles %>% unlist %>% unique %>% max
    
    model <- MIPModel() %>%
      add_variable(y[i], i = 1:number_of_cycles, type = "binary") %>%
      add_variable(x[i,j], i = 1:number_of_cycles, j = 1:number_of_nodes, type = "binary") %>%
      set_objective(sum_expr(utility[i] * y[i], i = 1:number_of_cycles), "max") %>%
      add_constraint(sum_expr(x[i,j], i = 1:number_of_cycles) <= 1, j = 1:number_of_nodes) %>%
      add_constraint(x[i,j] == y[i], i = 1:number_of_cycles, j = 1:number_of_nodes, cycle_contains_node(i,j)) %>%
      suppressWarnings() %>%
      solve_model(with_ROI(solver = "glpk", verbose=F))
    
    response <- model %>%
      get_solution(y[i]) %>%
      filter(value==1) %>%
      .$i
    
    selected.cycles <- cycles[response]
    
    return(selected.cycles)
    
  } else {
    return(cycles)
  }
  
}

perform.lab.crossmatches <- function(match.run, selected.cycles, node.info, match.info, cluster.fail.prob){
  
  get.cycle.edges <- function(cycle.node.list){
    
    edge.list <- data.frame(Donor=numeric(), Candidate=numeric())
    
    cycle.length <- length(cycle.node.list)
    
    for(node in 1:(cycle.length - 1)){
      edge.list <- bind_rows(edge.list, data.frame(Donor = cycle.node.list[node], Candidate = cycle.node.list[node+1]))
    }
    
    edge.list <- bind_rows(edge.list, data.frame(Donor = cycle.node.list[cycle.length], Candidate = cycle.node.list[1]))
    
    return(edge.list)
  }
  
  confirm.match <- function(candidate.type, match.status, cluster.fail.prob){
    
    match <- if_else(match.status == "Successful Lab Crossmatch", 1.0, as.double(rbinom(1, 1, 1-cluster.fail.prob[candidate.type])))
    
    return(match)
    
  } 
  
  node.info.reduced <- node.info %>%
    group_by(Node,Type) %>%
    summarize(Status = last(Status))
  
  match.info.reduced <- match.info %>%
    group_by(Donor, Candidate) %>%
    summarize(Status = last(Status))
  
  test.matches <- map_df(selected.cycles, get.cycle.edges, .id="Cycle") %>%
    left_join(node.info.reduced, by=c("Candidate"="Node")) %>%
    dplyr::select(-Status)  %>%
    rename(CandidateType = Type) %>%
    left_join(match.info.reduced, by=c("Donor"="Donor","Candidate"="Candidate")) %>%
    rename(MatchStatus = Status) %>%
    rowwise %>%
    mutate(Confirmed = confirm.match(CandidateType, MatchStatus, cluster.fail.prob)) %>%
    ungroup
    
  successful.cycles <- test.matches %>%
    group_by(Cycle) %>%
    summarize(N = n(), Confirmations = sum(Confirmed)) %>%
    ungroup %>%
    mutate(Success = N == Confirmations) %>%
    filter(Success) %>%
    .$Cycle
    
  final.test.matches <- test.matches %>%
    mutate(MatchRun = match.run,
           Status = case_when(Cycle %in% successful.cycles ~ "Transplanted",
                                   Confirmed == 1 ~ "Successful Lab Crossmatch",
                                   TRUE ~ "Failed Lab Crossmatch")) %>%
    select(Donor,Candidate,MatchRun,Status)
  
  return(final.test.matches)
}

set.transplant.info <- function(match.run, lab.crossmatches, node.info, match.info){
  
  node.info.reduced <- node.info %>%
    group_by(Node,Type) %>%
    summarize(Status = last(Status))
  
  transplant.info <- node.info.reduced %>%
    select(Node, Type) %>%
    right_join(lab.crossmatches, by=c("Node" = "Candidate")) %>%
    select(-Donor) %>%
    mutate(MatchRun = match.run, 
           Status = if_else(Status == "Transplanted","Transplanted","Selected; Not Transplanted"))
  
  
  transplanted.nodes <- transplant.info %>%
    filter(Status == "Transplanted") %>%
    .$Node
  
  transplanted.matches <- lab.crossmatches %>%
    filter(Status == "Transplanted") %>%
    select(Donor, Candidate)
  
  match.info.reduced <- match.info %>%
    group_by(Donor, Candidate) %>%
    summarize(Status = last(Status))
  
  corrected.matches <- match.info.reduced %>%
    filter(Donor %in% transplanted.nodes | Candidate %in% transplanted.nodes) %>%
    anti_join(transplanted.matches, by=c("Donor"="Donor","Candidate"="Candidate")) %>%
    mutate(MatchRun = match.run,
           Status = "Donor or Candidate was Transplanted")
  
  return(list(Nodes = transplant.info, Matches = corrected.matches))
  
}

check.renege <- function(match.run, node.info, match.info){ # Check for Matches that Renege
  
  renege.match <- function(candidate.type, cluster.fail.prob){
    
    match <- rbinom(1, 1, cluster.fail.prob[candidate.type])
    
    return(match)
    
  } 
  
  node.info.reduced <- node.info %>%
    group_by(Node,Type) %>%
    summarize(Status = last(Status))
  
  match.info.reduced <- match.info %>%
    group_by(Donor, Candidate) %>%
    summarize(Status = last(Status))
  
  reneged.matches <- match.info.reduced %>% 
    filter(Status %in% c("Virtual Crossmatch","Successful Lab Crossmatch")) %>%
    select(-Status) %>%
    left_join(node.info.reduced, by=c("Candidate"="Node")) %>%
    select(-Status) %>%
    rename(CandidateType = Type) %>%
    rowwise %>%
    mutate(Renege = renege.match(CandidateType, cluster.fail.prob)) %>%
    ungroup %>%
    filter(Renege == 1) %>%
    select(-Renege,-CandidateType) %>%
    mutate(MatchRun = match.run,
           Status = "Reneged")
  
  return(reneged.matches)
  
}
