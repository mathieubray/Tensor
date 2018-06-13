library(tidyverse)
library(ggraph)
library(ggforce)
library(gganimate)
library(ggplot2)
library(igraph)

set.seed(900709)

number.of.nodes <- 30
floor.val <- 12
runoff <- 5

vertex.arrival <- rexp(number.of.nodes,0.5) %>% cumsum %>% floor %>% if_else(. <= floor.val, floor.val, .) - floor.val

vertex.list <- tibble(name = as.character(1:number.of.nodes),
                      arrival = vertex.arrival)

edge.list <- map_df(1:max(vertex.arrival + runoff), function(i) {
  
  max.val <- vertex.arrival[vertex.arrival <= i] %>% length
  
  number.of.edges <- rpois(1,1.5)
  
  if (number.of.edges > 0){
    edges <- tibble(from = as.character(sample.int(n = max.val, size = number.of.edges, replace=T)),
                        to = as.character(sample.int(n = max.val, size = number.of.edges, replace=T)),
                        time = i) %>% unique()
    
    return(edges)
  } else {
    return(NULL)
  }
})

# Layout the vertices
subGr <- graph_from_data_frame(edge.list,directed=T,vertices=vertex.list)
lay <- create_layout(subGr, 'auto')


edgesAnim <- map_df(1:max(edge.list$time), function(i) {edge.list$timebin <- i; edge.list$appear <- if_else(edge.list$time <= i, "Appear", "Do Not Appear"); edge.list}) %>% 
  mutate(label = timebin)
verticesAnim <- map_df(1:max(edge.list$time), function(i) {lay$timebin <- i; lay}) %>%
  mutate(appear = if_else(arrival <= timebin, "Appear", "Do Not Appear"),
         label = paste0(name,": ",timebin))

edgesGraph <- graph_from_data_frame(edgesAnim, directed = T, vertices=vertex.list)


# Then we reassign the full graph with edge trails
attr(lay, 'graph') <- edgesGraph


# Now we create the graph with timebins as frame
p <- ggraph(graph = lay) +
  geom_edge_link(aes(frame = timebin, alpha = appear, size = appear),
                 edge_colour = "black",
                 start_cap=circle(5,"mm"),
                 end_cap=circle(5,"mm"),
                 arrow=arrow(angle=20, 
                             length = unit(0.1,"inches"),
                             type="closed"), 
                 data = get_edges()) +
  geom_node_point(data=verticesAnim,
                  mapping=aes(x=x,y=y,color=appear,frame=timebin),
                  size=12,alpha=0.6,show.legend=FALSE) +
  #geom_node_text(data=verticesAnim,
  #               mapping=aes(x=x, y=y,label=label,frame=timebin),
  #               size=6,show.legend=FALSE) +
  scale_alpha_manual(values=c(1,0)) +
  scale_color_manual(values=c("red","white")) +
  scale_edge_alpha_manual(values = c(1,0)) + 
  scale_edge_size_manual(values = c(1,0)) +
  theme_graph() + 
  theme(legend.position="none")

#p

# And then we animate
animation::ani.options(interval=0.5)
gg_animate(p, 'animation.gif', title_frame = FALSE)
gg_animate(p, 'animation.mp4', title_frame = FALSE)
