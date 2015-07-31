require(synapseClient)
synapseLogin()

res <- synQuery('select name,id,method from file where projectId==\'syn3455058\' and matrixType==\'sparse\'')

loadData <- function(syn){
  synMeta<-synGet(syn)
  load(synMeta@filePath)
  return(network)
}

SVAGeneExpressionSummary <- synGet('syn2757151')
cmcSVAGeneSummary <- read.delim(SVAGeneExpressionSummary@filePath,stringsAsFactors = FALSE)

SVATFExpressionSummary <- synGet('syn2757153')
cmcSVATFSummary <- read.delim(SVATFExpressionSummary@filePath, stringsAsFactors = FALSE)



require(metanetwork)
networks <- sapply(res$file.id,loadData)
names(networks) <- res$file.method
#netStat
netStats <- simplify2array(sapply(networks,networkStatistics))
netStats <- cbind(rownames(netStats),c(cmcSVATFSummary$MAPPED_genes, cmcSVAGeneSummary$MAPPED_genes),netStats)
colnames(netStats) <- c('EnsemblId','GeneId',res$file.method)
netStats <- data.frame(netStats,stringsAsFactors = FALSE)
netStats[,3:9]<-apply(netStats[,3:9],2,as.numeric)
netStats[,3:9] <- apply(netStats[,3:9],2,scale)

tcresult <- as.tableColumns(netStats)
cols <- tcresult$tableColumns
cols[[2]]$maximumSize <- as.integer(80)
fileHandleId <- tcresult$fileHandleId

projectId <- "syn3455058"

schema<-TableSchema(name="SVA Schizophrenia Network Statistics", parent=projectId, columns=cols)
table<-Table(schema, fileHandleId)
table<-synStore(table, retrieveData=TRUE)


geneRanks <- apply(netStats[,3:9],2,order,decreasing=T)
geneRanksName <- apply(geneRanks,2,function(x,g) return(g[x]),netStats$GeneId)

enrichRes <- apply(geneRanksName,2,enrichmentPath,schizophreniaHits)
enrichResPval <- sapply(enrichRes,function(x) return(-log10(as.numeric(x$pval))))


matplot(1:5000,enrichResPval[1:5000,],'l',lwd=3,xlab='Ranked Hubs',ylab='-log10 p-value',main='Enrichment for PGC2 GWAS LOCI based on hubness')
legend('topleft',colnames(enrichResPval),col=c(1:6,1),lty=1:7,lwd=2,cex=.6)


#get differentially expressed genes

deSyn <- synGet('syn3493933')
deSynFile <- read.delim(deSyn@filePath,stringsAsFactors=F)

targetListDE <- as.character(deSynFile$MAPPED_genes[1:488])
targetListDEens <- as.character(deSynFile$genes[1:488])
enrichResDe <- apply(geneRanksName,2,enrichmentPath,targetListDE)
enrichResPvalDE <- sapply(enrichResDe,function(x) return(-log10(as.numeric(x$pval))))

matplot(1:5000,enrichResPvalDE[1:5000,],'l',lwd=3,xlab='Ranked Hubs',ylab='-log10 pvalue',main='Enrichment for DE genes based on hubness')
legend('topright',colnames(enrichResPvalDE),col=c(1:6,1),lty=1:7,lwd=2,cex=.6)

pairs(netStats[,3:9],pch=16,col='grey')
def <- synTableQuery('SELECT * FROM syn3546459')
def@values$isGWAS <- netStats$GeneId%in%schizophreniaHits
def@values$isDE <- netStats$EnsemblId%in%targetListDEens
def <- synStore(def)


def <- synTableQuery('SELECT * FROM syn3546459')
colnames(def@values)
defMethod <- def@values[,c('sparrow','wgcna','aracne','genie3')]
a3 <- rowSums(defMethod)
plot(a3)
def@values <- cbind(def@values,a3)


geneIdForNetwork <- netStats$GeneId
geneIdForNetwork[geneIdForNetwork=='.'] <- netStats$EnsemblId[geneIdForNetwork=='.']

makeSparseNetworkFile(network = networks[[1]]%>%as.matrix(), geneId = geneIdForNetwork, file='~/Desktop/sparrowNet.csv')


a3 <- read.csv('~/Desktop/sparrowNetworkStatistics.csv',row.names=9,stringsAsFactors=FALSE)
a4 <- names(sort(a3$Stress,decreasing=T))
rankedList <- rownames(a3)[order(a3$BetweennessCentrality,decreasing=T)]
rankedList <- c(rankedList,setdiff(netStats$GeneId,rankedList))
blah <- enrichmentPath(rankedList = rankedList,targetList = schizophreniaHits)
plot(-log10(as.numeric(blah$pval))[1:2000])

#####run some quick analyses on what Thanneer generated
res2 <- synTableQuery('select * from syn3546459')


networkStat <- synQuery('select name,id from file where projectId==\'syn3455058\' and networkStatistic==1')

newTable <- res2@values
newTable2 <- vector('list',7)
i <- 1
for (i in 1:nrow(networkStat)){
  synObj <- synGet(networkStat$file.id[i])
  newTable2[[i]] <- read.delim(synObj@filePath,row.names=1)
}
names(newTable2) <- c('aracne','genie3','lasso','ridge','sparrow','tigress','wgcna')


getGeneList<-function(x,g){
  return(g[order(x,decreasing=T)])
}

metaGetGeneList <- function(X,g){
  return(apply(X,2,getGeneList,g))
}

masterList <- lapply(newTable2,metaGetGeneList,newTable$GeneId)

pgc2 <- newTable$GeneId[ newTable$isGWAS==TRUE]

enrichRes2 <- lapply(masterList,function(x,y) apply(x,2,enrichmentPath,pgc2), pgc2)
