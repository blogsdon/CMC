###
require(synapseClient)
synapseLogin()

foo <- synGet('syn4933783')
require(dplyr)
require(data.table)
load(foo@filePath)
####PPI

ppiNets <- synGet('syn4552965')
ppiNets@filePath
load(ppiNets@filePath)
ppiNew <- Matrix(trueMat + t(trueMat) !=0,sparse=TRUE)
diag(ppiNew) <- FALSE
ppiNCIdiff <- sparseNetwork - ppiNew
sum(ppiNCIdiff==1)/2

####cell type markers
source('networkAnnotation.R')
require(WGCNA)
sczModules <- makeModule(sparseNetwork)
sczSpecificEnrichments <- runEnrichmentAnalysis(sczModules)
sczSpecificEnrichments <- arrange(sczSpecificEnrichments,fdr)
bar <- synGet('syn4933752')
load(bar@filePath)
controlModules <- makeModule(sparseNetwork)
controlSpecificEnrichments <- runEnrichmentAnalysis(controlModules)
controlSpecificEnrichments <- arrange(controlSpecificEnrichments,fdr)



microglia <- filter(sczOverall$enrichments,GeneSetName=='Zhang:Microglia')
endothelial <- filter(sczOverall$enrichments,GeneSetName=='Zhang:Endothelial')
astrocyte <- filter(sczOverall$enrichments,GeneSetName=='Zhang:Astrocyte')
myelinOligos <- filter(sczOverall$enrichments,GeneSetName=='Zhang:MyelinOligos')
neuron <- filter(sczOverall$enrichments,GeneSetName=='Zhang:Neuron')
opc <- filter(sczOverall$enrichments,GeneSetName=='Zhang:OPC')

microgliaControl <- filter(controlSpecificEnrichments,GeneSetName=='Zhang:Microglia')
endothelialControl <- filter(controlSpecificEnrichments,GeneSetName=='Zhang:Endothelial')
astrocyteControl <- filter(controlSpecificEnrichments,GeneSetName=='Zhang:Astrocyte')
myelinOligosControl <- filter(controlSpecificEnrichments,GeneSetName=='Zhang:MyelinOligos')
neuronControl <- filter(controlSpecificEnrichments,GeneSetName=='Zhang:Neuron')
opcControl <- filter(controlSpecificEnrichments,GeneSetName=='Zhang:OPC')



#sczOverall <- list(network=sparseNetwork,
#                   modules=sczModules,
#                   enrichments=sczSpecificEnrichments)

controlOverall <- list(network=sparseNetwork,
                   modules=controlModules,
                   enrichments=controlSpecificEnrichments)

makeFilesForCytoscape(sczOverall,'sczNetwork.csv','sczModules.csv')
makeFilesForCytoscape(controlOverall,'controlNetwork.csv','controlModules.csv')

####quick PGC2
namedModule <- fread('sczMinusControlModules.csv') %>% data.frame
library(metanetwork)

enrichRes <- lapply(names(table(namedModule$modulelabels)),function(x) enrichment(schizophreniaHits,namedModule$GeneIDs[namedModule$modulelabels==x],namedModule$GeneIDs))

pvals <- sapply(1:length(enrichRes),function(i,x) return(x[[i]]$pval),enrichRes)
library(gap)
qqunif(pvals)
darkmagenta <- filter(sczOverall$enrichments,ComparisonName=='darkmagenta')
darkmagentaGenes <- filter(namedModule,modulelabels=='darkmagenta')$GeneIDs
cat(darkmagentaGenes,sep='\n',file='~/Desktop/darkmagenta.csv')

####differentialNetwork Analysis

diffNet <- sczOverall$network - controlOverall$network

sczSpecificMods <- makeModule(Matrix(diffNet==1,sparse=TRUE))
#nciSpecificMods <- makeModule(Matrix(diffNet==-1,sparse=TRUE))


sczSpecificEnrichments <- runEnrichmentAnalysis(sczSpecificMods,geneset=schizophreniaHits)
#nciSpecificEnrichments <- runEnrichmentAnalysis(nciSpecificMods)

sczMinusControl <- list(network=Matrix(diffNet==1,sparse=TRUE),
                   modules=sczSpecificMods,
                   enrichments=sczSpecificEnrichments)

#nciMinusAd <- list(network=Matrix(diffNet==-1,sparse=TRUE),
#                   modules=nciSpecificMods,
#                   enrichments=nciSpecificEnrichments)

sczMinusControl$enrichments <- arrange(sczMinusControl$enrichments,fdr)
#nciMinusAd$enrichments <- arrange(nciMinusAd$enrichments,fdr)
makeFilesForCytoscape(sczMinusControl,'sczMinusControlNetwork.csv','sczMinusControlModules.csv')
#makeFilesForCytoscape(nciMinusAd,'~/Projects/AMP-AD-Network-Analysis/nciMinusAdNetwork.csv','~/Projects/AMP-AD-Network-Analysis/nciMinusAdModules.csv')

filter(sczMinusControl$enrichments,GeneSetName=='Zhang:MyelinOligos')[1:5,]


####
sczE <- runEnrichmentAnalysis(sczOverall$modules,geneset=schizophreniaHits)
contE <- runEnrichmentAnalysis(controlOverall$modules,geneset=schizophreniaHits)
sczMiCon <- runEnrichmentAnalysis(sczMinusControl$modules,geneset = schizophreniaHits)

sczE<-arrange(sczE,pval)
contE <- arrange(contE,pval)
sczMiCon <- arrange(sczMiCon,pval)

filter(sczE,GeneSetName=='PGC2')[1:5,]
filter(contE,GeneSetName=='PGC2')[1:5,]
filter(sczMiCon,GeneSetName=='PGC2')[1:5,]
filter(sczE,ComparisonName=='red')[1:5,]


NoSVAGeneExpressionSummary <- synGet('syn2757151')
cmcNoSVAGeneSummary <- fread(NoSVAGeneExpressionSummary@filePath) %>% data.frame

NoSVATFExpressionSummary <- synGet('syn2757153')
cmcNoSVATFSummary <- fread(NoSVATFExpressionSummary@filePath) %>% data.frame

foo <- read.delim(NoSVATFExpressionSummary@filePath)

combinedSummary <- rbind(cmcNoSVATFSummary,cmcNoSVAGeneSummary)
colnames(combinedSummary) <- c('',colnames(foo))
diffGenes <- combinedSummary$MAPPED_genes[combinedSummary$adj.P.Val<0.05]

sczE <- runEnrichmentAnalysis(sczOverall$modules,geneset=diffGenes)
contE <- runEnrichmentAnalysis(controlOverall$modules,geneset=diffGenes)
sczMiCon <- runEnrichmentAnalysis(sczMinusControl$modules,geneset = diffGenes)

sczE<-arrange(sczE,fdr)
contE <- arrange(contE,fdr)
sczMiCon <- arrange(sczMiCon,fdr)

filter(sczE,GeneSetName=='PGC2')[1:5,]
filter(contE,GeneSetName=='PGC2')[1:5,]
filter(sczMiCon,GeneSetName=='PGC2')[1:5,]

####compare to WGCNA
foo <- synGet('syn3348753')
controlNetworkEnrichments <- fread(foo@filePath) %>% data.frame
filter(controlNetworkEnrichments,annot2=='Zhang:MyelinOligos')


bar <- synGet('syn3348772')
sczNetworkEnrichments <- fread(bar@filePath) %>% data.frame
filter(sczNetworkEnrichments,annot2=='Zhang:Neuron')

foobar <- synGet('syn3368830')
abbb <- fread(foobar@filePath) %>% data.frame


####genetic genesets
foo <- synGet('syn3368830')
bar <- fread(foo@filePath,data.table=F) 
#pgc2<-strsplit(bar$symAfterOverlap[250],'\\|')[[1]]
sczLists <- lapply(bar$symAfterOverlap,function(x) return(strsplit(x,'\\|')[[1]]))
names(sczLists) <- bar$Name


sczMods <- names(sort(table(sczOverall$modules$modulelabels),decreasing=T)[1:36])
controlMods <- names(sort(table(controlOverall$modules$modulelabels),decreasing=T)[1:31])
sczControlMods <- names(sort(table(sczMinusControl$modules$modulelabels),decreasing=T)[1:42])


sczE <- runEnrichmentAnalysisSingle(filter(sczOverall$modules,modulelabels%in%sczMods),geneset=sczLists)
contE <- runEnrichmentAnalysisSingle(filter(controlOverall$modules,modulelabels%in%controlMods),geneset=sczLists)
sczMiCon <- runEnrichmentAnalysisSingle(filter(sczMinusControl$modules,modulelabels%in%sczControlMods),geneset = sczLists)

sczE<-arrange(sczE,fdr)
contE <- arrange(contE,fdr)
sczMiCon <- arrange(sczMiCon,fdr)

#filter(sczE,GeneSetName=='PGC2')[1:5,]
bar$Name
filter(sczE,ComparisonName=='purple')[1:10,]
filter(sczE,GeneSetName=='Zhang:Astrocyte')[1:10,]

###Enrichr

GL_OBJ = synGet('syn4867851')
ALL_USED_IDs = GL_OBJ$properties$id
load(GL_OBJ@filePath)

#gsets = c("Achilles_fitness_decrease", "Achilles_fitness_increase", "Allen_Brain_Atlas_down", "Allen_Brain_Atlas_up",
#          "BioCarta", "CMAP_down", "CMAP_up", "Cancer_Cell_Line_Encyclopedia", "ChEA", "Cross_Species_Phenotype",
#          "Disease_Signatures_from_GEO_down", "Disease_Signatures_from_GEO_up", "Drug_Perturbations_from_GEO",
#          "ENCODE_Histone_Modifications_2013", "ESCAPE", "GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function",
#          "GeneSigDB", "Genome_Browser_PWMs.1", "HMDB_Metabolites", "HomoloGene", "Human_Gene_Atlas", "KEGG_2015",
#          "MGI_Mammalian_Phenotype", "MGI_Mammalian_Phenotype_Level_3", "MGI_Mammalian_Phenotype_Level_4", "MSigDB_Computational",
#          "MSigDB_Oncogenic_Signatures", "Mouse_Gene_Atlas", "NCI-60_Cancer_Cell_Lines", "NCI-Nature", 
#          "NURSA_Human_Endogenous_Complexome", "OMIM_Disease", "OMIM_Expanded", "PPI_Hub_Proteins", "Pfam_InterPro_Domains",
#          "Phosphatase_Substrates_from_DEPOD", "Reactome", "SILAC_Phosphoproteomics", "TF-LOF_Expression_from_GEO", 
#          "TargetScan_microRNA", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB",
#          "Transcription_Factor_PPIs", "Virus_Perturbations_from_GEO_down", "Virus_Perturbations_from_GEO_up", "WikiPathways_2015")

gsets = c("BioCarta", "CMAP_down", "CMAP_up", "ChEA",
          "Disease_Signatures_from_GEO_down", "Disease_Signatures_from_GEO_up", "Drug_Perturbations_from_GEO",
          "ENCODE_Histone_Modifications_2013", "GO_Biological_Process", "GO_Cellular_Component", "GO_Molecular_Function",
          "Genome_Browser_PWMs.1", "HMDB_Metabolites","KEGG_2015",
          "Mouse_Gene_Atlas",
          "OMIM_Disease", "OMIM_Expanded", "PPI_Hub_Proteins", "Pfam_InterPro_Domains",
          "Reactome", "TF-LOF_Expression_from_GEO", 
          "TargetScan_microRNA", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB",
          "Transcription_Factor_PPIs", "WikiPathways_2015")

GeneSets.Enrichr = GeneSets[gsets]


sczEnrichr <- runEnrichmentAnalysisDouble(filter(sczOverall$modules,modulelabels%in%sczMods),geneset=GeneSets.Enrichr)
contEnrichr <- runEnrichmentAnalysisDouble(filter(controlOverall$modules,modulelabels%in%controlMods),geneset=GeneSets.Enrichr)
sczMiConEnrichr <- runEnrichmentAnalysisDouble(filter(sczMinusControl$modules,modulelabels%in%sczControlMods),geneset = GeneSets.Enrichr)

sczEnrichr<-arrange(sczEnrichr,fdr)
contEnrichr <- arrange(contEnrichr,fdr)
sczMiConEnrichr <- arrange(sczMiConEnrichr,fdr)


makeModuleNetwork <- function(x,modKeep){
  network <- matrix(0,length(modKeep),length(modKeep))
  #restrict to necessary modules
  mods <- dplyr::filter(x$modules,modulelabels%in%modKeep)
  
  for (i in 1:(length(modKeep)-1)){
    for (j in (i+1):(length(modKeep))){
      ind1 <-  x$modules$GeneIDs[x$modules$modulelabels %in% modKeep[i]]
      ind2 <-  x$modules$GeneIDs[x$modules$modulelabels %in% modKeep[j]]			
      network[i,j] <- mean(x$network[ind1,ind2])
    }
  }
  network <- network + t(network)
  colnames(network) <- modKeep
  rownames(network) <- modKeep
  return(network)
  
}

sczModNetwork <- makeModuleNetwork(sczOverall,sczMods)
sczs<-scale(t(scale(sczModNetwork)))
sczs <- sczs+t(sczs)
el<-metanetwork::rankedEdgeList(sczs,symmetric=TRUE)
el$var1 <- colnames(sczModNetwork)[el$var1]
el$var2 <- colnames(sczModNetwork)[el$var2]
foo<-(filter(sczOverall$modules,modulelabels%in%sczMods))
foo$modulelabels <- as.character(foo$modulelabels)
write.csv(el,file='sczModNet.csv',quote=F)
bar <- (table(foo[,3]))
bar <- cbind(names(bar),bar)
colnames(bar) <- c('module','moduleSize')
write.csv(bar,file='sczModMod.csv',quote=F)
