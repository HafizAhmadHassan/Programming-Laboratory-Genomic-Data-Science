## Exercise 1
# - Create an empty directed graph with 10 nodes. Set color of all nodes to green and shape to square
# - Add 15 random edges to the graph
# - Add a vertex to the graph, color it to red and add edges: 3->6, 6->5. Set vertex shape to sphere.
# - Replace the first edge with the edge  6->3.
# - Name vertices with letters A-K. List all vertices and edges.
# - Plot the graph. The size of a vertex should depend on the number of incoming edges.
# - Plot the distribution of vertices degrees.
# - Create a heatmap from the adjacency matrix.

library(igraph)
library(ggplot2)

eg <- make_empty_graph(10, directed = T)
V(eg)$color <- "green"
V(eg)$shape <- "square"
plot(eg)

eg2 <- add.edges(eg, sample(1:10,30, replace = TRUE))

eg3 <- add.vertices(eg2, 1, color="red", shape="sphere") %>%
  add_edges(c(3,6, 6,5))
plot(eg3)

eg3[3,6] <- 0
eg3[6,3] <- 1
V(eg3)$label <- c(LETTERS[1:11])
list(V(eg3), E(eg3))
plot(eg3)

deg <-degree(eg3, mode="all")
V(eg3)$size <- deg*4
plot(eg3)

deg.dist <- degree_distribution(eg3, cumulative = T, mode="all")
plot(x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="green", 
     xlab="Degree", ylab="Cumulative Frequency")

netm <- get.adjacency(eg3, sparse=F)

palf <- colorRampPalette(c("gold", "dark orange")) 
heatmap(netm, Rowv = NA, Colv = NA, col = palf(20), 
        scale="none", margins=c(10,10))

## Exercise 2
# - Create an Erdos-Renyi random graph of 100 nodes and 120 edges and plot it.
# - Add weights to edges by sampling them from the range [1,20]
# - Compute the diameter and plot the network assigning two colors to the nodes 
# (one for the nodes belonging to the diameter and one for all the others), set also 
# a width to edges in function of the weight and set red color for edges having weight > 10,
# blue otherwise 
# - Remove edges that have weigth <= 10 and plot again
# - Remove disconnected nodes (you can do this in two ways: figure out the first one, for the second one try using the function "components")
# - Find communities based on greedy optimization of modularity

wei=sample(1:20,120,replace=TRUE)

er <- sample_gnm(n=100, m=120) %>% set_edge_attr("weight", value = wei)

plot(er)

#edge_attr(er)

diameter(er, directed=F)
diam=get_diameter(er,directed=F)
V(er)$color <- "yellow"
V(er)$color[diam] <- "green"
E(er)$width=E(er)$weight
E(er)$color='blue'
E(er)$color[E(er)$weight>10]='red'
plot(er, edge.arrow.mode=0)

er2 <- delete.edges(er, E(er)[E(er)$weight<=10])
plot(er2, edge.arrow.mode=0)

deg <- degree(er2, mode="all")
er3 <- delete.vertices(er2, deg[]==0)
plot(er3, edge.arrow.mode=0)

comps <- components(er3)$membership
colbar <- rainbow(max(comps)+1)
V(er3)$color <- colbar[comps+1]
plot(er3, layout=layout_with_fr, vertex.size=5, vertex.label=NA)

cfg <- cluster_fast_greedy(as.undirected(er3))
plot(cfg, as.undirected(er3))

## Exercise 3
# - Create an undirected star graph with 5 nodes and plot the graph.
# - Add edges between the following node pairs: 8 and 5, 6 and 3, 6 and 10, 5 and 2, 4 and 8,
# 2 and 7, 6 and 5, 7 and 9

# - Re-plot the graph, this time with a circular layout, ordered by vertex number.

# - Verify that the graph is connected.

# - Is the graph bipartite?

# - Calculate the graph's diameter.

# - Find the size of each clique in the graph. Hint: Use the sapply function.
# - Make a new plot of the graph, with the node size being relative to the nodes closeness, multiplied by 500.
# - Color the nodes of the graph: even nodes blue, odd nodes red. Hint: Using an if-else statement will make this more concise.

st <- make_star(10, mode="undirected")  #with 5 nodes is impossible to add the edge 8-5
st
plot(st, vertex.size=10, vertex.label=NA)

st=add.edges(st,edges=c(8,5,6,3,6,10,5,2,4,8,2,7,6,5,7,9))
st
plot(st, vertex.size=10, layout=layout_in_circle, layout=sort)

#plot(st, vertex.size=10, layout=layout_in_circle, layout=gorder(st))

is.connected(st)
is.bipartite(st)
diam <- get_diameter(st, directed=F)
diam
sapply(cliques(st), length)

#V(st)$color=rep(c('blue','red'), length(V(st)))

media <- ifelse(!as.logical(V(st)%%2), "blue", "red")

plot(st, vertex.size=closeness(st, mode="all", weights=NA)*500, layout=layout_in_circle, layout=sort, vertex.color=media)

## Exercise 4
# - Load the data from the csv "sociogram-employees-un.csv"
# and create an undirected graph. Name nodes as letters A to Y. Set node color to orange and shape to circle. 
# Set edge's color to blue and arrow size to 0.2. Plot the graph.
# - Find the largest cliques in the group.
# - How many maximal cliques are there?
# - Calculate the network cohesion.
# - Find the clusters based on betweenness.
# - Find the loop edges.

m = as.matrix(read.csv("D:\\uni/magistrale/programmazione/network/exercises/sociogram-employees-un.csv", header = FALSE))

g1<-graph.adjacency(adjmatrix=m,mode="undirected",weighted=TRUE,diag=FALSE)
V(g1)$label <- LETTERS[1:ncol(m)]
V(g1)$color='orange'
V(g1)$shape='circle'
plot(g1)

E(g1)$color='blue'
E(g1)$arrow.size=0.2
plot(g1)

largest_cliques(g1)
length(largest_cliques(g1))
cohesion(g1)
cluster_edge_betweenness(g1)
has.multiple(g1)#nope
which_multiple(g1, eids = E(g1))

## Exercise 5
# - Create an undirected graph from "net_nodes.csv" and "net_edges.csv" files
# - Compute the clustering coefficient of the network
# - Create a ring graph of 10 nodes and rename vertices with letters from k-t
# - Rewire the vertices of the previous graph connecting them in a neighborhood of 3
# - Create another graph by disjoint union of the previous two and add 3 random edges to connect them
# (use sample function)
# - Compute again the clustering coefficient of the new network

nodes <- read.csv("D:\\uni/magistrale/programmazione/network/exercises/net_nodes.csv", header=T, as.is=T)
links <- read.csv("D:\\uni/magistrale/programmazione/network/exercises/net_edges.csv", header=T, as.is=T)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
transitivity(net, type="global")

rn <- make_ring(10)
V(rn)$name <- letters[11:20]
plot(rn)

rn = connect.neighborhood(rn, 3)
plot(rn)
du = graph.disjoint.union(rn,net)
#x %du% y
du <- du + edges(sample(V(du), 3*2, replace = T))
transitivity(du, type = "global")
plot(du)

