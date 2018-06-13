
source("AlphaNTF.R")
source("BetaNTF.R")
source("FastHALS.R")
source("BCD.R")

source("TensorGenerators.R")

set.seed(90707)


# Generate Adjacency Tensor

#per_mr <- 4 # Number of Nodes Arriving between each Match Run
#num_mrs <- 25 # Number of Match Runs

#prob_type <- c(0.4,0.1,0.1,0.4)
#prob_connect <- matrix(c(0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05,0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05),nrow=4)

#xm_prob <- c(0.9,0.75,0.75,0.5) # Probability of Connection Remaining After Evaluation


#kpd.info <- generate.kpd(per_mr, num_mrs, prob_type, prob_connect, xm_prob)

#kpd <- kpd.info$kpd
#node_types <- kpd.info$node_types


# Generate Adjacency Tensor

per_mr<- 10 # Number of Nodes Arriving between each Match Run
num_mrs <- 20 # Number of Match Runs

prob_type <- c(0.5,0.5)



prob_connect <- matrix(c(0.25,0.15,0.15,0.05),nrow=2)

xm_prob <- c(0.9,0.75) # Probability of Connection Remaining After Evaluation


kpd.info <- generate.kpd(per_mr, num_mrs, prob_type, prob_connect, xm_prob)

kpd <- kpd.info$kpd
node_types <- kpd.info$node_types







# Generate Initial Matrix and Perform Decompositions

num_components <- 2

init <- initialize.matrix(per_mr, num_mrs, num_components)

maxits <- 1000
tol <- 1e-5
eps <- 1e-5

decomp.cp <- cp(kpd, num_components, max_iter = maxits, tol = tol)
decomp.hals <- fastHALS(kpd, init, max_k = maxits, tol = tol)
decomp.alpha <- alphaNTF(kpd, init, alpha = 1.5, max_k = maxits, tol = tol)
decomp.beta <- betaNTF(kpd, init, beta = 0.5, max_k = maxits, tol=tol)
decomp.bcd <- BCD(kpd, init, delta = 0.5, max_k = maxits, tol=tol)


# CP
#B <- decomp.cp$U
#norms <- decomp.cp$all_resids
#value <- decomp.cp$fnorm_resid

# HALS
B.hals <- decomp.hals$A
norms.hals <- decomp.hals$objs
value.hals <- decomp.hals$objs %>% last

# Alpha
B.alpha <- decomp.alpha$A
norms.alpha <- decomp.alpha$objs
value.alpha <- decomp.alpha$objs %>% last

# Beta
B.beta <- decomp.beta$A
norms.beta <- decomp.beta$objs
value.beta <- decomp.beta$objs %>% last

# BCD
#B <- decomp.bcd$A
#norms <- decomp.bcd$objs
#value <- decomp.bcd$objs %>% last


# Component Plot

library(ggplot2)

node.frame <- rbind(data.frame(B.hals[[1]],Node_Type=as.factor(node_types),u="Donor"),
                    data.frame(B.hals[[2]],Node_Type=as.factor(node_types),u="Candidate"))


ggplot(data=node.frame, aes(x=X1,y=X2,color=Node_Type)) +
  facet_wrap(~u) +
  geom_point(size=4,alpha=0.5) +
  theme_bw() +
  xlab("Component 1") +
  ylab("Component 2")

node.frame <- rbind(data.frame(B.hals[[1]],Node_Type=as.factor(node_types),u="Donor"),
                    data.frame(B.hals[[2]],Node_Type=as.factor(node_types),u="Candidate"))

# 5 10 15 20 sequential match runs
# Side by side with KPD
# PEnalized group optimization
# New penalty function...
# Coefficient in front of utility for the group


ggplot(data=node.frame, aes(x=X1,y=X2,color=Node_Type)) +
  facet_wrap(~u) +
  geom_point(size=4,alpha=0.5) +
  theme_bw() +
  xlab("Component 1") +
  ylab("Component 2")



## Convergence Plot

norm.table <- data.frame(Iteration=1:length(norms.hals),Objective=norms.hals) #%>%
  #filter (Iteration > 200)

ggplot(data=norm.table,aes(x=Iteration,y=Objective)) +
  geom_path() +
  theme_bw() +
  ggtitle("Convergence Plot") +
  geom_hline(yintercept=value.hals,color="red",linetype="dashed") +
  geom_label(aes(label=paste0(round(value.hals,3)),x=10,y=value.hals, color="red"),show.legend=FALSE)



#################################################


num_components <- 3

init <- initialize.matrix(per_mr, num_mrs, num_components)

maxits <- 1000
tol <- 1e-5
eps <- 1e-5

#decomp.cp <- cp(kpd, num_components, max_iter = maxits, tol = tol)
decomp.hals <- fastHALS(kpd, init, max_k = maxits, tol = tol)
#decomp.alpha <- alphaNTF(kpd, init, alpha = 1.5, max_k = maxits, tol = tol)
#decomp.beta <- betaNTF(kpd, init, beta = 0.5, max_k = maxits, tol=tol)
#decomp.bcd <- BCD(kpd, init, delta = 0.5, max_k = maxits, tol=tol)


# CP
#B <- decomp.cp$U
#norms <- decomp.cp$all_resids
#value <- decomp.cp$fnorm_resid

# HALS
B.hals <- decomp.hals$A
norms.hals <- decomp.hals$objs
value.hals <- decomp.hals$objs %>% last

# Alpha
#B.alpha <- decomp.alpha$A
#norms.alpha <- decomp.alpha$objs
#value.alpha <- decomp.alpha$objs %>% last

# Beta
#B.beta <- decomp.beta$A
#norms.beta <- decomp.beta$objs
#value.beta <- decomp.beta$objs %>% last

# BCD
#B <- decomp.bcd$A
#norms <- decomp.bcd$objs
#value <- decomp.bcd$objs %>% last


# Component Plot

library(ggplot2)
library(scatterplot3d)

node.frame <- rbind(data.frame(B.hals[[1]],Node_Type=as.factor(node_types),u="Donor"),
                    data.frame(B.hals[[2]],Node_Type=as.factor(node_types),u="Candidate"))

node.frame.donor <- node.frame %>% filter(u=="Donor")

colors <- c("red","#EE69F00")
colors <- colors[as.numeric(node.frame.donor$Node_Type)]

scatterplot3d(node.frame.donor,
              pch=16,
              type="h",
              color=colors)

node.frame.candidates <- node.frame %>% filter(u=="Candidate")

colors <- c("red","#EE69F00")
colors <- colors[as.numeric(node.frame.donor$Node_Type)]

scatterplot3d(node.frame.candidates,
              pch=16,
              type="h",
              color=colors)



