
source("AlphaNTF.R")
source("BetaNTF.R")
source("FastHALS.R")
source("BCD.R")

source("TensorGenerators.R")

set.seed(90707)


# Generate Adjacency Tensor

per_mr <- 4 # Number of Nodes Arriving between each Match Run
num_mrs <- 25 # Number of Match Runs

prob_type <- c(0.4,0.1,0.1,0.4)
prob_connect <- matrix(c(0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05,0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05),nrow=4)

xm_prob <- c(0.9,0.75,0.75,0.5) # Probability of Connection Remaining After Evaluation


kpd.info <- generate.kpd(per_mr, num_mrs, prob_type, prob_connect, xm_prob)

kpd <- kpd.info$kpd
node_types <- kpd.info$node_types



# Generate Initial Matrix and Perform Decompositions

num_components <- 2

init <- initialize.matrix(per_mr, num_mrs, num_components)

maxits <- 1000
tol <- 1e-5

decomp.cp <- cp(kpd, num_components, max_iter = maxits, tol = tol)
decomp.hals <- fastHALS(kpd, init, max_k = maxits, tol = tol)
decomp.alpha <- alphaNTF(kpd, init, alpha = 1.5, max_k = maxits, tol = tol)
decomp.beta <- betaNTF(kpd, init, beta = 0.5, max_k = maxits, tol=tol)
decomp.bcd <- BCD(kpd, init, delta = 0.5, max_k = maxits, tol=tol)


# CP
B <- decomp.cp$U
norms <- decomp.cp$all_resids
value <- decomp.cp$fnorm_resid

# HALS
B <- decomp.hals$A
norms <- decomp.hals$objs
value <- decomp.hals$objs %>% last

# Alpha
B <- decomp.alpha$A
norms <- decomp.alpha$objs
value <- decomp.alpha$objs %>% last

# Beta
B <- decomp.beta$A
norms <- decomp.beta$objs
value <- decomp.beta$objs %>% last

# BCD
B <- decomp.bcd$A
norms <- decomp.bcd$objs
value <- decomp.bcd$objs %>% last


# Component Plot

library(ggplot2)

node.frame <- rbind(data.frame(B[[1]],Node_Type=as.factor(node_types),u="Donor"),
                    data.frame(B[[2]],Node_Type=as.factor(node_types),u="Candidate"))


ggplot(data=node.frame, aes(x=X1,y=X2,color=Node_Type)) +
  facet_wrap(~u) +
  geom_point(size=4,alpha=0.5) +
  theme_bw() +
  xlab("Component 1") +
  ylab("Component 2") +
  ggtitle("Component Plot")



## Convergence Plot

norm.table <- data.frame(Iteration=1:length(norms),Objective=norms) #%>%
  #filter (Iteration > 200)

ggplot(data=norm.table,aes(x=Iteration,y=Objective)) +
  geom_path() +
  theme_bw() +
  ggtitle("Convergence Plot") +
  geom_hline(yintercept=value,color="red",linetype="dashed") +
  geom_label(aes(label=paste0(round(value,3)),x=10,y=value, color="red"),show.legend=FALSE)





