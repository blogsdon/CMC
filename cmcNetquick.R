#get enrichment results
grabNetworkAnalysisResults <- function(method,sparsityMethod,disease,projectId){
  require(synapseClient)
  synapseLogin()
  queryStatement <- paste0('select name,id from file where projectId==\'',projectId,'\' and disease==\'',disease,'\' and sparsityMethod==\'',sparsityMethod,'\' and method ==\'',method,'\'')
  cat(queryStatement,'\n')
  foo <- synQuery(queryStatement)
  bar <- lapply(foo$file.id,synGet)
  load(bar[[1]]@filePath)
  require(dplyr)
  require(data.table)
  sparrow2Sparse <- list(network=sparseNetwork,
                         modules=fread(bar[[2]]@filePath)%>%data.frame,
                         enrichments=fread(bar[[3]]@filePath)%>%data.frame)
  sparrow2Sparse$enrichments <- arrange(sparrow2Sparse$enrichments,pval)
  return(sparrow2Sparse)
}


makeEnrichmentSummary <- function(x){
  require(dplyr)
  require(data.table)
  nMods <- sum(table(x$modules$moduleNumber)>20)
  cat(nMods)
  nameKeep <- filter(x$modules,moduleNumber%in%1:nMods)[,3] %>%
    table %>%
    sort(decreasing=T) %>%
    names
  cat(nameKeep,'\n')


  internal <- function(y,x){
    green <- filter(x$enrichments,ComparisonName==y)
    greenSplit <- lapply(names(table(green$category)),function(x,y) return(filter(y,category==x)),green)
    names(greenSplit) <- names(table(green$category))
    summarygreen <- do.call(rbind,lapply(greenSplit,function(x) return(x[1:5,c('GeneSetName','category','pval','noverlap','OR','fdr','ComparisonName')])))
    try(summarygreen <- arrange(summarygreen,pval),silent=TRUE)
    return(summarygreen)
  }

  res <- lapply(nameKeep,internal,x)
  names(res) <- nameKeep
  return(res)
}


controlResult<- grabNetworkAnalysisResults('rankconsensus','sparrow2Bonferroni','Control','syn3455058')
sczResult<- grabNetworkAnalysisResults('rankconsensus','sparrow2Bonferroni','Schizophrenia','syn3455058')

controlWgcnaModsObj <- synGet('syn3348750')
controlMods <- fread(controlWgcnaModsObj@filePath,data.table=F)

sczWgcnaModsObj <- synGet('syn3348769')
sczMods <- fread(sczWgcnaModsObj@filePath,data.table=F)

sczCombMod <- merge(sczResult$modules,sczMods,by.x = 'GeneIDs', by.y = 'Gene')[,c(1,3,4,5)]
controlCombMod <- merge(controlResult$modules,controlMods,by.x = 'GeneIDs', by.y = 'Gene')[,c(1,3,4,5)]




w1 <- which(table(sczCombMod$modulelabels)>20)
w2 <- which(table(sczCombMod$Module)>20)
name1 <- names(table(sczCombMod$modulelabels))[w1]
name2 <- names(table(controlCombMod$Module))[w2]

outerSapply <- function(FUN,X,Y,...){
  require(dplyr)
  internal <- function(X,Y,FUN,...){
    return(Y%>% sapply(FUN,X,...))
  }
  return(X %>% sapply(internal,Y,FUN,...))
}
require(metanetwork)

internalEnr <- function(x,y,z){
  return(enrichment(x,y,z)$enr)
}

internalPval <- function(x,y,z){
  return(enrichment(x,y,z)$pval)
}

sczFiltered <- dplyr::filter(sczCombMod,modulelabels%in%name1 & Module%in%name2)

mnLt <- lapply(unique(sczFiltered$modulelabels),function(i,x,y){return(x[y==i])},sczFiltered$GeneIDs,sczFiltered$modulelabels)
names(mnLt) <- unique(sczFiltered$modulelabels)

wgLt <- lapply(unique(sczFiltered$Module),function(i,x,y){return(x[y==i])},sczFiltered$GeneIDs,sczFiltered$Module)
names(wgLt) <- unique(sczFiltered$Module)

enrTabScz <- outerSapply(internalEnr,mnLt,wgLt,sczFiltered$GeneIDs)
library(gplots)
library(RColorBrewer)

heatmap.2(sqrt(t(enrTabScz)),trace='none',scale = 'none',margins = c(6,7))
#tab3 <- table(sczFiltered$modulelabels,sczFiltered$Module)
#heatmap(tab3)


controlSummary <- makeEnrichmentSummary(controlResult)
sczSummary <- makeEnrichmentSummary(sczResult)

catSum <- do.call(rbind,lapply(names(sczSummary),function(y,x){return(filter(x[[y]],category=='Cross_Species_Phenotype'))},sczSummary))
catSum <- arrange(catSum,fdr)
catSum[1:10,]



makeFilesForCytoscape <- function(x,networkFile,moduleFile){
  require(dplyr)
  nMods <- sum(table(x$modules$moduleNumber)>20)
  mods <- filter(x$modules,moduleNumber%in%1:nMods)
  allAdj <- x$network[mods$GeneIDs,mods$GeneIDs] %>% as.matrix
  library(metanetwork)
  alledgeList <- rankedEdgeList(allAdj,symmetric=T)
  exprFile <- synGet('syn4259377')
  require(data.table);
  expr <- fread(exprFile@filePath,header=T)
  expr <- data.frame(expr)

  hgnc <- expr$hgnc_symbol
  names(hgnc) <- expr$ensembl_gene_id

  alledgeList$var1 <- hgnc[mods$GeneIDs][alledgeList$var1]
  alledgeList$var2 <- hgnc[mods$GeneIDs][alledgeList$var2]

  alledgeList <- filter(alledgeList,!is.na(var1))
  alledgeList <- filter(alledgeList,!is.na(var2))
  alledgeList <- filter(alledgeList,var1!='')
  alledgeList <- filter(alledgeList,var2!='')
  write.csv(alledgeList,file=networkFile,quote=F,row.names=F)
  mod2 <- dplyr::select(mods,GeneIDs,modulelabels)
  mod2$GeneIDs <- hgnc[mod2$GeneIDs]
  write.csv(mod2,file=moduleFile,quote=F,row.names=F)
}
redSczResult <- dplyr::filter(sczResult$enrichments,ComparisonName=='red')

#yellowSczResult <- dplyr::filter(sczResult$enrichments,ModuleName=='yellow')
require(data.table)

accRes <- fread('ACC.ensembl.DxSCZ.DE.KNOWN_AND_NO_HIDDEN_CONFOUND.ADJUSTED.VOOM_NORMALIZED.LIMMA_FIT.tsv',data.table=F,stringsAsFactors=FALSE)

pfcRes <- fread('DLPFC.ensembl.DxSCZ.DE.KNOWN.ADJUSTED.VOOM_NORMALIZED.LIMMA_FIT.tsv',data.table=F,stringsAsFactors=FALSE)

accAllGenes <- accRes$genes[accRes$adj.P.Val<0.05]
accPosGenes <- accRes$genes[accRes$adj.P.Val<0.05 & sign(accRes$t)==1]
accNegGenes <- accRes$genes[accRes$adj.P.Val<0.05 & sign(accRes$t)==-1]
pfcAllGenes <- pfcRes$genes[pfcRes$adj.P.Val<0.05]
pfcPosGenes <- pfcRes$genes[pfcRes$adj.P.Val<0.05 & sign(pfcRes$t)==1]
pfcNegGenes <- pfcRes$genes[pfcRes$adj.P.Val<0.05 & sign(pfcRes$t)==-1]
accAllUnique <- setdiff(accAllGenes,pfcAllGenes)
pfcAllUnique <- setdiff(pfcAllGenes,accAllGenes)
allShared <- intersect(accAllGenes,pfcAllGenes)

accPosUnique <- setdiff(accPosGenes,pfcPosGenes)
pfcPosUnique <- setdiff(pfcPosGenes,accPosGenes)
posShared <- intersect(accPosGenes,pfcPosGenes)

accNegUnique <- setdiff(accNegGenes,pfcNegGenes)
pfcNegUnique <- setdiff(pfcNegGenes,accNegGenes)
negShared <- intersect(accNegGenes,pfcNegGenes)

require(metanetwork)

sczNeuron <- filter(sczResult$modules,modulelabels=='blue') %>% select(GeneIDs)
sczNeuron <- sczNeuron[,1]

conNeuron <- filter(controlResult$modules,modulelabels=='turquoise') %>% select(GeneIDs)
conNeuron <- conNeuron[,1]
allGenes <- controlResult$modules$GeneIDs


fisherWrapper <- function(moduleGenes,annotationGenes,allGenes){
  a00 <- sum(!(allGenes%in%moduleGenes) & !(allGenes%in%annotationGenes))
  a10 <- sum((allGenes%in%moduleGenes) & !(allGenes%in%annotationGenes))
  a01 <- sum(!(allGenes%in%moduleGenes) & (allGenes%in%annotationGenes))
  a11 <- sum((allGenes%in%moduleGenes) & (allGenes%in%annotationGenes))
  bar <- matrix(c(a00,a10,a01,a11),2,2)
  #print(bar)
  foo <- fisher.test(bar,alternative='greater')
  return(list(p.value=foo$p.value,
              estimate=foo$estimate,
              countMatrix=bar))
}


internal <- function(i,x){
  geneList <- x$GeneIDs[which(x$modulelabels==i)]
}
sczTable <- table(sczResult$modules$modulelabels)
sczTable <- sczTable[which(sczTable>20)]
sczMods <- names(sczTable)

sczModulesList <- lapply(sczMods,internal,sczResult$modules)
names(sczModulesList) <- sczMods


controlTable <- table(controlResult$modules$modulelabels)
controlTable <- controlTable[which(controlTable>20)]
controlMods <- names(controlTable)

controlModulesList <- lapply(controlMods,internal,controlResult$modules)
names(controlModulesList) <- controlMods


pfcAllSczResults <- lapply(sczModulesList,fisherWrapper,negShared,allGenes)
ORs <-sapply(pfcAllSczResults,function(x){return(x$estimate)})
p.values <- sapply(pfcAllSczResults,function(x){return(x$p.value)})
print(ORs[p.values< 0.05/length(p.values)])
print(p.values[p.values< 0.05/length(p.values)])
barplot(ORs[order(p.values)])

