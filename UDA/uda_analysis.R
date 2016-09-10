library(synapseClient)
synapseLogin()
foo <- synGet('syn7222231')
load(foo@filePath)

graph1 <- utilityFunctions::convertSparseMatrixToGraph(bicNetworks$network)


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

pgc2DistTab <- data.frame(pgc2Dist = pgc2Dist,ensembl_gene_id=names(pgc2Dist),stringsAsFactors=F)
pgc2DistTab <- merge(pgc2DistTab,geneMappingTable,by='ensembl_gene_id')
pgc2DistTab <- dplyr::arrange(pgc2DistTab,pgc2Dist)


deg <- dplyr::filter(geneMappingTable2,geneSet=='CMC-DE-Sz') %>% dplyr::select(ensembl_gene_id)

degDist <- metanetwork::computeDriverDistance(deg$ensembl_gene_id,graph1)

degDistTab <- data.frame(degDist = degDist,ensembl_gene_id=names(degDist),stringsAsFactors=F)
degDistTab <- merge(degDistTab,geneMappingTable,by='ensembl_gene_id')
degDistTab <- dplyr::arrange(degDistTab,degDist)


lof <- dplyr::filter(geneMappingTable2,geneSet=='SCZ.LoF') %>% dplyr::select(ensembl_gene_id)

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



