#new prototype analysis

require(synapseClient)
synapseLogin()

a <- synGet('syn4553095')
load(a@filePath)

#sparse graph
require(igraph)
graph2 <- graph_from_adjacency_matrix(sparseNetwork,mode='undirected')
res <- cluster_fast_greedy(graph2)
spar_mods <- membership(res)
keep<-which(table(spar_mods)>20)
spar_mods <- spar_mods[which(spar_mods%in%keep)]

table(spar_mods)



b <- synGet('syn4550165')
key <- read.csv(b@filePath)






net <- as.matrix(sparseNetwork)
ne <- rowSums(net)


ne2 <- which.max(ne)
getNeighb <- function(index,x){
  return(which(x[index,]))
}

neighb <- getNeighb(ne2,net)
neighb2 <- unique(unlist(lapply(neighb,getNeighb,net)))
neighb3 <- unique(c(ne2,neighb2,neighb))
neighb4 <- unique(unlist(lapply(neighb3,getNeighb,net)))
neighb5 <- unique(c(neighb4,neighb3))
#neighb <- which(net[9171,])

b <- synGet('syn4550165')
key <- read.csv(b@filePath)



cat(as.matrix(key)[neighb5,2],file='~/Desktop/test.csv',sep='\n')
#b <- kmeans(net[1:1000,1:1000],centers = 20)
