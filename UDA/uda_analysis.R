library(synapseClient)
synapseLogin()
foo <- synGet('syn7222231')
load(foo@filePath)

##########make a plot of the bic path
png(file='~/Desktop/bicCMCplot.png',pointsize=44,height=1280,width=1600)
plot(2*(8e4:1e5),log10(bicNetworks$bicPath[8e4:1e5]),'l',lwd=3,col='red',xlab='Number of Edges in the Network',ylab='log10 BIC fit criterion',main='Model selection for GGM Network')
dev.off()





###########download RIMBAnet file
rimbaObj <- synGet('syn7113487')
rimba <- read.delim(rimbaObj@filePath,stringsAsFactors=F,row.names=NULL,header=F)

ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                         dataset = 'hsapiens_gene_ensembl',
                         host='www.ensembl.org')

genes<-biomaRt::getBM(attributes = c('ensembl_gene_id','external_gene_name'),
             filters='ensembl_gene_id',
             values=colnames(bicNetworks$network),
             mart=ensembl)
longForm <- which(bicNetworks$network + t(bicNetworks$network) !=0,T)

longForm2 <- longForm
longForm2[,1] <- colnames(bicNetworks$network)[longForm[,1]]
longForm2[,2] <- colnames(bicNetworks$network)[longForm[,2]]

hgnc <- genes$external_gene_name
names(hgnc) <- genes$ensembl_gene_id
longForm <- longForm2
longForm[,1] <- hgnc[longForm2[,1]]
longForm[,2] <- hgnc[longForm2[,2]]

longForm <- longForm[which(!is.na(longForm[,1])),]
longForm <- longForm[which(!is.na(longForm[,2])),]
dim(rimba)

rimbaEdges <- apply(rimba,1,paste0,collapse='_')
mNedges <- apply(longForm,1,paste0,collapse='_')

length(rimbaEdges)
length(mNedges)
rimbaEdges2 <- rimbaEdges[rimbaEdges%in%mNedges]
rimbaEdges2 <- t(sapply(rimbaEdges2,function(x){return(strsplit(x,'_')[[1]])}))
write.csv(rimbaEdges2,file='~/Desktop/intersectNet.csv',quote=F,row.names=F)
#173810/2
tab1 <- cbind(c(269534113,7009),c(169972,3835))

fisher.test(tab1)

graph1 <- utilityFunctions::convertSparseMatrixToGraph(bicNetworks$network)
library(VennDiagram)
venn.plot <- venn.diagram(list(
  rimbaEdges = rimbaEdges,
  metanetworkEdges = mNedges),
  filename= NULL
)


barplot((c(173810,10844,3835)))


fastGreedy <- igraph::fastgreedy.community(graph1)

#geneSet <- sample(V(graphNew)%>%names,100)

geneSetsObj <- synGet('syn7217924')

#munge gene sets

geneSets <- read.delim(geneSetsObj@filePath,stringsAsFactors=F,row.names=NULL,header=F)
colnames(geneSets) <- c('external_gene_name','geneSet')

geneMappingTable <- utilityFunctions::convertEnsemblToHgnc(names(V(graph1)))

geneMappingTable2 <- merge(geneMappingTable,geneSets,by='external_gene_name')

pgc2 <- dplyr::filter(geneMappingTable2,geneSet=='PGC2') %>% dplyr::select(ensembl_gene_id)


pgc2Dist <- metanetwork::computeDriverDistance(pgc2$ensembl_gene_id,graph1)
pgc2Dist2 <- metanetwork::computeDriverDistancePvalue(pgc2$ensembl_gene_id,graph1,nsamp=100)

pgc2DistTab <- data.frame(pgc2Dist = pgc2Dist,ensembl_gene_id=names(pgc2Dist),stringsAsFactors=F)
pgc2DistTab <- merge(pgc2DistTab,geneMappingTable,by='ensembl_gene_id')
pgc2DistTab <- dplyr::arrange(pgc2DistTab,pgc2Dist)


deg <- dplyr::filter(geneMappingTable2,geneSet=='CMC-DE-Sz') %>% dplyr::select(ensembl_gene_id)

degDist <- metanetwork::computeDriverDistance(deg$ensembl_gene_id,graph1)

degDistTab <- data.frame(degDist = degDist,ensembl_gene_id=names(degDist),stringsAsFactors=F)
degDistTab <- merge(degDistTab,geneMappingTable,by='ensembl_gene_id')
degDistTab <- dplyr::arrange(degDistTab,degDist)


lof <- dplyr::filter(geneMappingTable2,geneSet=='SCZ.LoF') %>% dplyr::select(ensembl_gene_id)

degDist <- metanetwork::computeDriverDistance(deg$ensembl_gene_id,graph1)
lofDist <- metanetwork::computeDriverDistance(lof$ensembl_gene_id,graph1)

lofDistTab <- data.frame(lofDist = lofDist,ensembl_gene_id=names(lofDist),stringsAsFactors=F)
lofDistTab <- merge(lofDistTab,geneMappingTable,by='ensembl_gene_id')
lofDistTab <- dplyr::arrange(lofDistTab,lofDist)



ns <- dplyr::filter(geneMappingTable2,geneSet=='SCZ.NS') %>% dplyr::select(ensembl_gene_id)

nsDist <- metanetwork::computeDriverDistance(ns$ensembl_gene_id,graph1)

nsDistTab <- data.frame(nsDist = nsDist,ensembl_gene_id=names(nsDist),stringsAsFactors=F)
nsDistTab <- merge(nsDistTab,geneMappingTable,by='ensembl_gene_id')
nsDistTab <- dplyr::arrange(nsDistTab,nsDist)


aggr <- dplyr::filter(geneMappingTable2,geneSet=='SCZ.NS' | geneSet=='SCZ.LoF' | geneSet == 'PGC2' | geneSet == 'CMC-DE-Sz') %>% dplyr::select(ensembl_gene_id)

aggrDist <- metanetwork::computeDriverDistance(unique(aggr$ensembl_gene_id),graph1)

aggrDistTab <- data.frame(aggrDist = aggrDist,ensembl_gene_id=names(aggrDist),stringsAsFactors=F)
aggrDistTab <- merge(aggrDistTab,geneMappingTable,by='ensembl_gene_id')
aggrDistTab <- dplyr::arrange(aggrDistTab,aggrDist)

rando <- sample(names(V(graph1)),638)
randoDist <- metanetwork::computeDriverDistance(rando,graph1)

randoDistTab <- data.frame(randoDist = randoDist,ensembl_gene_id=names(randoDist),stringsAsFactors=F)
randoDistTab <- merge(randoDistTab,geneMappingTable,by='ensembl_gene_id')
randoDistTab <- dplyr::arrange(randoDistTab,randoDist)



mergedTable <- merge(degDistTab,pgc2DistTab,by='ensembl_gene_id')
mergedTable <- dplyr::mutate(mergedTable,aggScore = degDist+pgc2Dist) %>% dplyr::arrange(aggScore)

View(mergedTable)


View(degDistTab)
plot(foobar)



