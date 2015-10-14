##schizophrenia network analysis


getAllNetworks <- function(x){
  y <- synGet(x)
  load(y@filePath)
  return(sparseNetwork)
}

require(synapseClient)
synapseLogin()

synResult <- synQuery('select name,id from file where projectId==\'syn3455058\' and networkStorageType==\'sparse\' and disease==\'Schizophrenia\' and normalization==\'None\'')

allSczNets <- lapply(synResult$file.id,getAllNetworks)
consensus <- matrix(0,16423,16423)
for (i in 1:length(allSczNets)){
  cat('i: ',i,'\n')
  consensus <- consensus+as.matrix(allSczNets[[i]])
  gc()
}

load('~/Desktop/consensusScz.rda')
require(synapseClient)
synapseLogin()

b <- synGet('syn4550165')
key <- read.csv(b@filePath,stringsAsFactors = FALSE)
consensus <- consensus/149

require(igraph)
graph2 <- graph_from_adjacency_matrix(consensus,mode='undirected',weighted=TRUE)
res <- cluster_fast_greedy(graph2)
spar_mods <- membership(res)
spar_mods2 <- membership(res)
keep<-which(table(spar_mods)>20)
spar_mods <- spar_mods[which(spar_mods%in%keep)]

table(spar_mods)


