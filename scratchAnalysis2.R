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
spar_mods2 <- membership(res)
keep<-which(table(spar_mods)>20)
spar_mods <- spar_mods[which(spar_mods%in%keep)]

table(spar_mods)


b <- synGet('syn4550165')
key <- read.csv(b@filePath,stringsAsFactors = FALSE)
newName <- (key$geneName)[key$ensemblId %in% names(spar_mods)]
names(spar_mods) <- newName
spar_mods <- spar_mods[-which(names(spar_mods)=='.')]

getClusters <- lapply(unique(spar_mods),function(x,y) return(names(y[which(y==x)])),spar_mods)
enrichRes<- sapply(getClusters,function(list2,list1,list3) return(enrichment(list1,list2,list3)), schizophreniaHits, key$geneName)

