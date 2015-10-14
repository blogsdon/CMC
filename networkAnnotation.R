makeModule <- function(x){
  g = igraph::graph.adjacency(x, mode = 'undirected', weighted = T, diag = F)
  
  # Get modules using fast.greedy method (http://arxiv.org/abs/cond-mat/0408187)
  mod = igraph::fastgreedy.community(g)
  collectGarbage()
  
  # Get individual clusters from the igraph community object
  clust.numLabels = igraph::membership(mod)
  collectGarbage()
  
  # Change cluster number to color labels
  labels = WGCNA::labels2colors(clust.numLabels)
  
  # Get results
  geneModules = data.frame(GeneIDs = igraph::V(g)$name,
                           moduleNumber = as.numeric(clust.numLabels), 
                           modulelabels = labels)
  
  return(geneModules)  
}

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

populateADenrichments <- function(x){
  require(dplyr)
  require(data.table)
  nMods <- sum(table(x$modules$moduleNumber)>20)
  cat(nMods)
  nameKeep <- filter(x$modules,moduleNumber%in%1:nMods)[,3] %>%  
    table %>% 
    sort(decreasing=T) %>% 
    names
  cat(nameKeep,'\n')
  internal <- function(x,y){
    filter(y,ComparisonName==x & category =='ADRelated')
  }
  res <- lapply(nameKeep,internal,x$enrichments)
  names(res) <- nameKeep
  return(res)
}

runEnrichmentAnalysis <- function(MOD,geneset=NULL){
  library(synapseClient)
  library(dplyr)
  library(WGCNA)
  library(tools)
  library(stringr)
  library(igraph)
  library(data.table)
  library(biomaRt)
  filterGeneSets <- function(GeneLists, # List of lists
                             genesInBackground, # background set of genes
                             minSize = 10,
                             maxSize = 1000){
    GeneLists = lapply(GeneLists, 
                       function(x, genesInBackground){
                         x = lapply(x, 
                                    function(x, genesInBackground){
                                      return(intersect(x, genesInBackground))
                                    },
                                    genesInBackground)
                         return(x)
                       }, 
                       genesInBackground)
    
    GeneLists = lapply(GeneLists, 
                       function(x, minSize, maxSize){
                         len = sapply(x, length)
                         x = x[len>minSize & len<maxSize]
                         return(x)
                       },
                       minSize,
                       maxSize)
    len = sapply(GeneLists, length)
    GeneLists = GeneLists[len != 0]
    
    return(GeneLists)
  }
  
  # Function to perform Fishers enrichment analysis
  fisherEnrichment <- function(genesInSignificantSet, # A character vector of differentially expressed or some significant genes to test
                               genesInGeneSet, # A character vector of genes in gene set like GO annotations, pathways etc...
                               genesInBackground # Background genes that are 
  ){
    genesInSignificantSet = intersect(genesInSignificantSet, genesInBackground) # back ground filtering
    genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
    genesInGeneSet = intersect(genesInGeneSet, genesInBackground) # back ground filtering
    genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
    
    pval = fisher.test(
      matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
               length(intersect(genesInGeneSet, genesInNonSignificantSet)),
               length(intersect(genesOutGeneSet, genesInSignificantSet)),
               length(intersect(genesOutGeneSet, genesInNonSignificantSet))), 
             nrow=2, ncol=2),
      alternative="greater")
    OR = (length(intersect(genesInGeneSet, genesInSignificantSet)) * length(intersect(genesOutGeneSet, genesInNonSignificantSet))) / (length(intersect(genesInGeneSet, genesInNonSignificantSet)) * length(intersect(genesOutGeneSet, genesInSignificantSet)))
    return(data.frame(pval = pval$p.value,
                      ngenes = length(genesInGeneSet),
                      noverlap = length(intersect(genesInGeneSet, genesInSignificantSet)),
                      OR = OR 
    )
    )
  }
  
  # Function to convert rownames to first column of a df
  rownameToFirstColumn <- function(DF,colname){
    DF <- as.data.frame(DF)
    DF[,colname] <- row.names(DF)
    DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
    return(DF)
  }
  
  # Download AD related gene sets from synapse
  GL_OBJ = synGet('syn4893059');
  ALL_USED_IDs = c(GL_OBJ$properties$id)
  load(GL_OBJ@filePath)
  GeneSets$PGC2 = geneset
  GeneSets.AD = list(ADRelated = GeneSets)

  
  
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = MOD$GeneIDs, mart = ensembl)
  backGroundGenes = unique(ensg2hgnc$hgnc_symbol)
  
  MOD = merge(MOD, ensg2hgnc, by.x = 'GeneIDs', by.y = 'ensembl_gene_id', all.x=T)
  ############################################################################################################
  
  ############################################################################################################
  #### Filter gene list ####
  GeneSets.AD = filterGeneSets(GeneSets.AD, backGroundGenes, minSize = 1, maxSize = 10000)
  if(!is.null(geneset)){
    GeneSets.Custom = filterGeneSets(geneset,backGroundGenes,minSize=1,maxSize=1e4)
  }
  GeneSets = c(GeneSets.AD,GeneSets.Custom)
  ############################################################################################################
  
  ############################################################################################################
  #### Perform enrichment analysis ####
  # Perform enrichment analysis (for modules greater than 20 genes only)
  enrichResults = list()
  for (name in unique(MOD$modulelabels)){
    genesInModule = unique(MOD$hgnc_symbol[MOD$modulelabels == name])  
    if (length(genesInModule) > 20){
      tmp = lapply(GeneSets,
                   function(x, genesInModule, genesInBackground){
                     tmp = as.data.frame(t(sapply(x, fisherEnrichment, genesInModule, genesInBackground)))
                     tmp = rownameToFirstColumn(tmp,'GeneSetName')
                     return(tmp)
                   },
                   unique(genesInModule), unique(backGroundGenes))
      
      for (name1 in names(tmp))
        tmp[[name1]]$category = name1
      
      enrichResults[[name]] = as.data.frame(rbindlist(tmp))
      enrichResults[[name]]$fdr = p.adjust(enrichResults[[name]]$pval, 'fdr')
    } else {
      enrichResults[[name]] = data.frame(GeneSetName = NA, pval = NA, ngenes = NA, noverlap = NA, OR = NA, category = NA, fdr = NA)
    }
    writeLines(paste0('Completed ',name))  
  }
  
  # Write results to file
  for(name in names(enrichResults))
    enrichResults[[name]]$ComparisonName = name
  enrichmentResults = as.data.frame(rbindlist(enrichResults))
  enrichmentResults$ngenes = unlist(enrichmentResults$ngenes)
  enrichmentResults$noverlap = unlist(enrichmentResults$noverlap)
  enrichmentResults$fdr = unlist(enrichmentResults$fdr)
  enrichmentResults$OR = unlist(enrichmentResults$OR)
  enrichmentResults$pval = unlist(enrichmentResults$pval)
  return(enrichmentResults)
}


runEnrichmentAnalysisSingle <- function(MOD,geneset=NULL){
  library(synapseClient)
  library(dplyr)
  library(WGCNA)
  library(tools)
  library(stringr)
  library(igraph)
  library(data.table)
  library(biomaRt)
  filterGeneSets <- function(GeneLists, # List of lists
                             genesInBackground, # background set of genes
                             minSize = 10,
                             maxSize = 1000){
    GeneLists = lapply(GeneLists, 
                       function(x, genesInBackground){
                         x = lapply(x, 
                                    function(x, genesInBackground){
                                      return(intersect(x, genesInBackground))
                                    },
                                    genesInBackground)
                         return(x)
                       }, 
                       genesInBackground)
    
    GeneLists = lapply(GeneLists, 
                       function(x, minSize, maxSize){
                         len = sapply(x, length)
                         x = x[len>minSize & len<maxSize]
                         return(x)
                       },
                       minSize,
                       maxSize)
    len = sapply(GeneLists, length)
    GeneLists = GeneLists[len != 0]
    
    return(GeneLists)
  }
  
  # Function to perform Fishers enrichment analysis
  fisherEnrichment <- function(genesInSignificantSet, # A character vector of differentially expressed or some significant genes to test
                               genesInGeneSet, # A character vector of genes in gene set like GO annotations, pathways etc...
                               genesInBackground # Background genes that are 
  ){
    genesInSignificantSet = intersect(genesInSignificantSet, genesInBackground) # back ground filtering
    genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
    genesInGeneSet = intersect(genesInGeneSet, genesInBackground) # back ground filtering
    genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
    
    pval = fisher.test(
      matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
               length(intersect(genesInGeneSet, genesInNonSignificantSet)),
               length(intersect(genesOutGeneSet, genesInSignificantSet)),
               length(intersect(genesOutGeneSet, genesInNonSignificantSet))), 
             nrow=2, ncol=2),
      alternative="greater")
    OR = (length(intersect(genesInGeneSet, genesInSignificantSet)) * length(intersect(genesOutGeneSet, genesInNonSignificantSet))) / (length(intersect(genesInGeneSet, genesInNonSignificantSet)) * length(intersect(genesOutGeneSet, genesInSignificantSet)))
    return(data.frame(pval = pval$p.value,
                      ngenes = length(genesInGeneSet),
                      noverlap = length(intersect(genesInGeneSet, genesInSignificantSet)),
                      OR = OR 
    )
    )
  }
  
  # Function to convert rownames to first column of a df
  rownameToFirstColumn <- function(DF,colname){
    DF <- as.data.frame(DF)
    DF[,colname] <- row.names(DF)
    DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
    return(DF)
  }
  
  # Download AD related gene sets from synapse
  #GL_OBJ = synGet('syn4893059');
  #ALL_USED_IDs = c(GL_OBJ$properties$id)
  #load(GL_OBJ@filePath)
  #rm(GeneSets)
  #gc()
  #GeneSets <- list()
  #GeneSets = geneset
  GeneSets.AD = list(SCZRelated = geneset)
  
  
  
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = MOD$GeneIDs, mart = ensembl)
  backGroundGenes = unique(ensg2hgnc$hgnc_symbol)
  
  MOD = merge(MOD, ensg2hgnc, by.x = 'GeneIDs', by.y = 'ensembl_gene_id', all.x=T)
  ############################################################################################################
  
  ############################################################################################################
  #### Filter gene list ####
  GeneSets.AD = filterGeneSets(GeneSets.AD, backGroundGenes, minSize = 1, maxSize = 10000)
  #if(!is.null(geneset)){
  #  GeneSets.Custom = filterGeneSets(geneset,backGroundGenes,minSize=1,maxSize=1e4)
  #}
  GeneSets = c(GeneSets.AD)
  ############################################################################################################
  
  ############################################################################################################
  #### Perform enrichment analysis ####
  # Perform enrichment analysis (for modules greater than 20 genes only)
  enrichResults = list()
  for (name in unique(MOD$modulelabels)){
    genesInModule = unique(MOD$hgnc_symbol[MOD$modulelabels == name])  
    if (length(genesInModule) > 20){
      tmp = lapply(GeneSets,
                   function(x, genesInModule, genesInBackground){
                     tmp = as.data.frame(t(sapply(x, fisherEnrichment, genesInModule, genesInBackground)))
                     tmp = rownameToFirstColumn(tmp,'GeneSetName')
                     return(tmp)
                   },
                   unique(genesInModule), unique(backGroundGenes))
      
      for (name1 in names(tmp))
        tmp[[name1]]$category = name1
      
      enrichResults[[name]] = as.data.frame(rbindlist(tmp))
      enrichResults[[name]]$fdr = p.adjust(enrichResults[[name]]$pval, 'fdr')
    } else {
      enrichResults[[name]] = data.frame(GeneSetName = NA, pval = NA, ngenes = NA, noverlap = NA, OR = NA, category = NA, fdr = NA)
    }
    writeLines(paste0('Completed ',name))  
  }
  
  # Write results to file
  for(name in names(enrichResults))
    enrichResults[[name]]$ComparisonName = name
  enrichmentResults = as.data.frame(rbindlist(enrichResults))
  enrichmentResults$ngenes = unlist(enrichmentResults$ngenes)
  enrichmentResults$noverlap = unlist(enrichmentResults$noverlap)
  enrichmentResults$fdr = unlist(enrichmentResults$fdr)
  enrichmentResults$OR = unlist(enrichmentResults$OR)
  enrichmentResults$pval = unlist(enrichmentResults$pval)
  return(enrichmentResults)
}
#truevec <- trueMat[which(upper.tri(trueMat))]


runEnrichmentAnalysisDouble <- function(MOD,geneset=NULL){
  library(synapseClient)
  library(dplyr)
  library(WGCNA)
  library(tools)
  library(stringr)
  library(igraph)
  library(data.table)
  library(biomaRt)
  filterGeneSets <- function(GeneLists, # List of lists
                             genesInBackground, # background set of genes
                             minSize = 10,
                             maxSize = 1000){
    GeneLists = lapply(GeneLists, 
                       function(x, genesInBackground){
                         x = lapply(x, 
                                    function(x, genesInBackground){
                                      return(intersect(x, genesInBackground))
                                    },
                                    genesInBackground)
                         return(x)
                       }, 
                       genesInBackground)
    
    GeneLists = lapply(GeneLists, 
                       function(x, minSize, maxSize){
                         len = sapply(x, length)
                         x = x[len>minSize & len<maxSize]
                         return(x)
                       },
                       minSize,
                       maxSize)
    len = sapply(GeneLists, length)
    GeneLists = GeneLists[len != 0]
    
    return(GeneLists)
  }
  
  # Function to perform Fishers enrichment analysis
  fisherEnrichment <- function(genesInSignificantSet, # A character vector of differentially expressed or some significant genes to test
                               genesInGeneSet, # A character vector of genes in gene set like GO annotations, pathways etc...
                               genesInBackground # Background genes that are 
  ){
    genesInSignificantSet = intersect(genesInSignificantSet, genesInBackground) # back ground filtering
    genesInNonSignificantSet = base::setdiff(genesInBackground, genesInSignificantSet)
    genesInGeneSet = intersect(genesInGeneSet, genesInBackground) # back ground filtering
    genesOutGeneSet = base::setdiff(genesInBackground,genesInGeneSet)
    
    pval = fisher.test(
      matrix(c(length(intersect(genesInGeneSet, genesInSignificantSet)),             
               length(intersect(genesInGeneSet, genesInNonSignificantSet)),
               length(intersect(genesOutGeneSet, genesInSignificantSet)),
               length(intersect(genesOutGeneSet, genesInNonSignificantSet))), 
             nrow=2, ncol=2),
      alternative="greater")
    OR = (length(intersect(genesInGeneSet, genesInSignificantSet)) * length(intersect(genesOutGeneSet, genesInNonSignificantSet))) / (length(intersect(genesInGeneSet, genesInNonSignificantSet)) * length(intersect(genesOutGeneSet, genesInSignificantSet)))
    return(data.frame(pval = pval$p.value,
                      ngenes = length(genesInGeneSet),
                      noverlap = length(intersect(genesInGeneSet, genesInSignificantSet)),
                      OR = OR 
    )
    )
  }
  
  # Function to convert rownames to first column of a df
  rownameToFirstColumn <- function(DF,colname){
    DF <- as.data.frame(DF)
    DF[,colname] <- row.names(DF)
    DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
    return(DF)
  }
  
  # Download AD related gene sets from synapse
  #GL_OBJ = synGet('syn4893059');
  #ALL_USED_IDs = c(GL_OBJ$properties$id)
  #load(GL_OBJ@filePath)
  #rm(GeneSets)
  #gc()
  #GeneSets <- list()
  #GeneSets = geneset
  GeneSets.AD = geneset
  
  
  
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  ensg2hgnc = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = MOD$GeneIDs, mart = ensembl)
  backGroundGenes = unique(ensg2hgnc$hgnc_symbol)
  
  MOD = merge(MOD, ensg2hgnc, by.x = 'GeneIDs', by.y = 'ensembl_gene_id', all.x=T)
  ############################################################################################################
  
  ############################################################################################################
  #### Filter gene list ####
  GeneSets.AD = filterGeneSets(GeneSets.AD, backGroundGenes, minSize = 1, maxSize = 10000)
  #if(!is.null(geneset)){
  #  GeneSets.Custom = filterGeneSets(geneset,backGroundGenes,minSize=1,maxSize=1e4)
  #}
  GeneSets = c(GeneSets.AD)
  ############################################################################################################
  
  ############################################################################################################
  #### Perform enrichment analysis ####
  # Perform enrichment analysis (for modules greater than 20 genes only)
  enrichResults = list()
  for (name in unique(MOD$modulelabels)){
    genesInModule = unique(MOD$hgnc_symbol[MOD$modulelabels == name])  
    if (length(genesInModule) > 20){
      tmp = lapply(GeneSets,
                   function(x, genesInModule, genesInBackground){
                     tmp = as.data.frame(t(sapply(x, fisherEnrichment, genesInModule, genesInBackground)))
                     tmp = rownameToFirstColumn(tmp,'GeneSetName')
                     return(tmp)
                   },
                   unique(genesInModule), unique(backGroundGenes))
      
      for (name1 in names(tmp))
        tmp[[name1]]$category = name1
      
      enrichResults[[name]] = as.data.frame(rbindlist(tmp))
      enrichResults[[name]]$fdr = p.adjust(enrichResults[[name]]$pval, 'fdr')
    } else {
      enrichResults[[name]] = data.frame(GeneSetName = NA, pval = NA, ngenes = NA, noverlap = NA, OR = NA, category = NA, fdr = NA)
    }
    writeLines(paste0('Completed ',name))  
  }
  
  # Write results to file
  for(name in names(enrichResults))
    enrichResults[[name]]$ComparisonName = name
  enrichmentResults = as.data.frame(rbindlist(enrichResults))
  enrichmentResults$ngenes = unlist(enrichmentResults$ngenes)
  enrichmentResults$noverlap = unlist(enrichmentResults$noverlap)
  enrichmentResults$fdr = unlist(enrichmentResults$fdr)
  enrichmentResults$OR = unlist(enrichmentResults$OR)
  enrichmentResults$pval = unlist(enrichmentResults$pval)
  return(enrichmentResults)
}
#truevec <- trueMat[which(upper.tri(trueMat))]
