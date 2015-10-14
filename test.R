#
require(synapseClient)
synapseLogin()
a <- synGet('syn4810725')
load(a@filePath)
consensus <- consensus/149

b <- synGet('syn4550165')
key <- read.csv(b@filePath,stringsAsFactors = FALSE)

require(igraph)
network2 <- (consensus>90/149)
graph2 <- graph_from_adjacency_matrix(network2,mode='undirected')
res <- cluster_fast_greedy(graph2)
spar_mods <- membership(res)
spar_mods2 <- membership(res)
keep<-which(table(spar_mods)>20)
spar_mods <- spar_mods[which(spar_mods%in%keep)]

table(spar_mods)
names(spar_mods2) <- key$geneName


ppisyn <- synGet('syn4552965')
load(ppisyn@filePath)



a <- synGet('syn4553095')
load(a@filePath)


hubness <- rowSums(sparseNetwork)

