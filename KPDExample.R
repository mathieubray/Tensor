
source("TensorGenerator.R")

set.seed(90707)

### Initialize Nodes

nodes.per.matchrun <- 10 # Number of Nodes Arriving before each Match Run
number.of.matchruns <- 20 # Number of Match Runs
number.of.nodes <- number.of.matchruns * nodes.per.matchrun
cluster.connect.prob <- matrix(c(0.25,0.15,0.15,0.05),nrow=2)
cluster.fail.prob <- c(0.1,0.25) # Probability of Connection Remaining After Evaluation
max.cycle.size <- 3

node.info <- initialize.nodes(nodes.per.matchrun, number.of.matchruns, c(0.5,0.5))

matrix.list <- list() # Initialize a Matrix List



### Get Initial Matches

match.info <- virtual.crossmatch(1, node.info, cluster.connect.prob)

current.matrix <- get.adjacency.matrix(match.info, number.of.nodes)

matrix.list[[1]] <- current.matrix



### Perform Match Runs

for (match.run in 1:number.of.matchruns){
  
  print(paste0("Match Run ", match.run, " Ready"))
  
  nodes <- node.info %>%
    filter(MatchRun <= match.run) %>%
    .$Node %>%
    max
  
  # Collect Cycles
  cycles <- collect.cycles(current.matrix, nodes, max.cycle.size)
  
  # Select Transplants
  selected.cycles <- select.cycles(cycles)
  
  if (length(selected.cycles) > 0){
    
    # Verify Selected Matches and Perform Transplants
    lab.crossmatches <- perform.lab.crossmatches(match.run, selected.cycles, node.info, match.info, cluster.fail.prob)
    transplants <- set.transplant.info(match.run, lab.crossmatches, node.info, match.info)
    
    match.info <- bind_rows(match.info, lab.crossmatches, transplants$Matches)
    node.info <- bind_rows(node.info, transplants$Nodes)
  
  }
  
  if (match.run != number.of.matchruns){
    
    # Verify Matches that Renege
    renege.nodes <- check.renege(match.run + 1, node.info, match.info)
    
    match.info <- bind_rows(match.info, renege.nodes)
    
    # Prepare Next Crossmatch
    new.matches <- virtual.crossmatch(match.run + 1, node.info, cluster.connect.prob)
    
    match.info <- bind_rows(match.info, new.matches)
    
    current.matrix <- get.adjacency.matrix(match.info, number.of.nodes)
    
    matrix.list[[match.run + 1]] <- current.matrix
    
  }
  
}

kpd <- as.tensor(matrix.list)

















##################################

#kpd <- kpd.info$kpd
#node.list.new <- kpd.info$node.list
#match.list.new <- kpd.info$match.list


node.list.joined <- node.info %>%
  filter(Status == "Joined") %>%
  rename(TimeJoined = MatchRun) %>%
  select(-Status)


node.list.transplanted <- node.info %>%
  filter(Status == "Transplanted") %>%
  rename(TimeTransplanted = MatchRun) %>%
  select(-Status)


node.list.new <- node.list.joined %>%
  left_join(node.list.transplanted, by=c("Node"="Node","Type"="Type")) %>%
  mutate(Censored = is.na(TimeTransplanted),
         TimeTransplanted = if_else(Censored,20.0,as.double(TimeTransplanted)))

node.list.selected <- node.info %>%
  filter(Status == "Selected; Not Transplanted") %>%
  select(-Status)


ggplot(data = node.list.new) +
  geom_point(aes(x=Node,y=TimeJoined,color=as.factor(Type))) +
  geom_point(aes(x=Node,y=TimeTransplanted,color=as.factor(Type),shape=Censored)) +
  geom_point(data=node.list.selected,aes(x=Node,y=MatchRun,color=as.factor(Type)), shape=4) +
  geom_segment(aes(x=Node,xend=Node,y=TimeJoined,yend=TimeTransplanted,color=as.factor(Type))) +
  ylim(0,20) +
  coord_flip() +
  theme_bw()





