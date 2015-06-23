makeMultiSparseNetwork <- function(sparsitySyn, networkSyn, geneSyn,uploadFolder,networkName,nNodes,fileName){
  #move to metanetworkSynapse
  sparsitySyn <- 'syn4546859'
  networkSyn <- 'syn4545862'
  geneSyn <- 'syn4550155'
  
  library(Matrix)
  sparObj <- synGet(sparsitySyn)
  networkSyn <- synGet(networkSyn)
  geneSyn <- synGet(geneSyn)
  
  spar <- read.csv(sparObj@filePath,stringsAsFactors=F,row.names=1)
  network <- read.csv(networkSyn@filePath,stringsAsFactors=F,row.names=1)
  gene <- read.csv(geneSyn@filePath,stringsAsFactors=F,row.names=1)
  
  edgeListToMatrix <- function(spar,edgeList,geneName,nNodes){
    network2 <- matrix(0,nNodes,nNodes)
    print(dim(network2))
    internal <- function(x,y,n){
      return((y-1)*n+x)
    }
    
    #for (i in 1:spar){
    #  network[edgeList[i,1],edgeList[i,2]]<-edgeList[i,3]
    #}
    avec <- internal(edgeList[1:spar,1],edgeList[1:spar,2],nNodes)
    print(avec[1:5])
    print(length(avec))
    print(max(avec))
    print(dim(network2))
    network2[avec] <- as.numeric(edgeList[1:spar,3])
    network <- network+t(network)
    print(dim(network))
    colnames(network2) <- geneName
    rownames(network2) <- geneName
    network <- Matrix(network2,sparse=TRUE)
    return(network)
  }
  
  nNodes <- 16423
  allNetworks<-sapply(spar$V2[spar$V2<nrow(network)],edgeListToMatrix,network,rownames(gene),nNodes)
  save(allNetworks,file=fileName)
  
  synObj <- File(fileName,parentId=uploadFolder)
  anno<- synGetAnnotations(networkSyn)
  
  ###UPDATE ANNOTATIONS!!
  #anno$... = ...
  
  ###Add provenance
  #act <- Activity...
  #act <- storeEntity...
  #generatedBy(synObj) <- act
  
  ###Store
  #synObj <- synStore(synObj)
  
}