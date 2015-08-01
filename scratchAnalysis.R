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

library(metanetwork)
e3 <- rankedEdgeList(net,symmetric=TRUE)
gc()
e3$var1 <- key$geneName[e3$var1]
e3$var2 <- key$geneName[e3$var2]

noName<-which(e3=='.',T)
rmIndex <- unique(noName[,1])

e3 <- e3[-rmIndex,]
write.table(e3,col.names = F,row.names = F,sep='\t',file='~/Desktop/sparrowNet',quote=F)
modules <- readLines('~/Desktop/cluster')
modules <- lapply(modules,function(x) strsplit(x,'\t')[[1]])
moduleIndex <- lapply(modules,function(x,y) which(y%in%x), key$geneName)

netModularized <- makeModuleNetwork(moduleIndex,net)
netModularized <- (netModularized + t(netModularized)) != 0
ee1 <- rankedEdgeList(netModularized,symmetric = TRUE)
ee1$var1 <- colnames(netModularized)[ee1$var1]
ee1$var2 <- colnames(netModularized)[ee1$var2]
write.table(ee1,col.names=F,row.names=F,sep='\t',file='~/Desktop/sparrowNetM',quote=F)

mods <- readLines('~/Desktop/clusterM')
names(modules) <- paste0('m',1:length(modules))
mods <- lapply(mods,function(x) strsplit(x,'\t')[[1]])
modIndex <- lapply(mods,function(x,y) which(y%in%x), colnames(netModularized))

netModularized2 <- makeModuleNetwork(modIndex,netModularized)
netModularized2 <- (netModularized2 + t(netModularized2))
heatmap(netModularized2)


metaModules <- lapply(mods,function(x,y) return(unlist(unique(y[x]))),modules)
metaModuleIndex <- lapply(metaModules, function(x,y) which(y%in%x), key$geneName)

writeNewEdges <- function(net,name){
  e4 <- rankedEdgeList(net,symmetric = TRUE)
  gc()
  e4$var1 <- key$geneName[e4$var1]
  e4$var2 <- key$geneName[e4$var2]

  noName<- which(e4=='.',T)
  rmIndex <- unique(noName[,1])

  e4 <- e4[-rmIndex,]
  write.table(e4,col.names= F, row.names = F, sep='\t',file=paste0('~/Desktop/',name),quote=F)
}
modules2 <- readLines('~/Desktop/cluster2')
modules2 <- lapply(modules2,function(x) strsplit(x,'\t')[[1]])
moduleIndex2 <- lapply(modules2,function(x,y) which(y%in%x), key$geneName)
netModularized2 <- makeModuleNetwork(moduleIndex2,netModularized)

writeNewEdges(netModularized2,'sparrowNet3')

modules3 <- readLines('~/Desktop/cluster3')
modules3 <- lapply(modules3,function(x) strsplit(x,'\t')[[1]])
moduleIndex3 <- lapply(modules3,function(x,y) which(y%in%x), key$geneName)

#which(key$geneName=='KCNN2')
#which(key$geneName=='FBXL13')

##########make heatmap of drivers
heatDriver <- function(expr,network,driverIndex){
  require(gplots)
  neighbors <- which(network[driverIndex,])
  covarianceNeighbors <- cor(expr[,neighbors],expr)
  expandNeighbors <- unique(c(neighbors,which(covarianceNeighbors>sqrt(.5),T)[,2]))
  require(RColorBrewer)
  mypalette <- colorRampPalette(c('blue','yellow'))
  heatmap.2(cor(expr[,expandNeighbors]),trace='none',scale='none',col=mypalette)
  return(hclust(dist(t(expr[,expandNeighbors]))))
}


ptpn11<-heatDriver(as.matrix(sczData),net,4436)
ptpn11mod2<-cutree(ptpn11,3)
cat(as.character(key$geneName)[key$ensemblId%in%names(which(ptpn11mod2==1))],sep='\n',file='~/Desktop/ptpn11mod1')
