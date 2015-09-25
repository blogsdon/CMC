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
    green <- filter(x$enrichments,ModuleName==y)
    greenSplit <- lapply(names(table(green$category)),function(x,y) return(filter(y,category==x)),green)
    names(greenSplit) <- names(table(green$category))
    summarygreen <- do.call(rbind,lapply(greenSplit,function(x) return(x[1:5,c('GeneSetName','category','pval','noverlap')])))
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

sczSummary <- makeEnrichmentSummary(sczResult)
#yellowSczResult <- dplyr::filter(sczResult$enrichments,ModuleName=='yellow')

