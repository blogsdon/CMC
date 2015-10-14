#script to grab SVA CMC data

makeDataMatrices <- function(geneSyn,tfSyn,sczFile,controlFile,bpFile){
  require(synapseClient)
  require(data.table)
  require(dplyr)
  synapseLogin()
  
  SVAGeneExpression <- synGet(geneSyn)
  cmcSVAGeneExpression <- read.delim(SVAGeneExpression@filePath)
  #cmcSVAGeneExpression <- fread(SVAGeneExpression@filePath) %>% data.frame
  
  SVATFExpression <- synGet(tfSyn)
  cmcSVATFExpression <- read.delim(SVATFExpression@filePath)
  #cmcSVATFExpression <- fread(SVATFExpression@filePath) %>% data.frame
  
  #combine gene expression data
  cmcExpression <- rbind(cmcSVATFExpression,cmcSVAGeneExpression)
  
  #download metadata
  syncmcMetaData <- synGet('syn2299154')
  metadata <- read.delim(syncmcMetaData@filePath,sep=',',stringsAsFactors=FALSE,check.names=F)
  
  metadata_freeze = metadata[-which(metadata[['DLPFC_RNA_report: Exclude?']] == 1), ]
  
  cmcExpressionFreeze <- cmcExpression[,which(colnames(cmcExpression)%in%as.character(metadata_freeze[['DLPFC_RNA_isolation: Sample RNA ID']]))]
  
  scz_samples <- metadata_freeze[['DLPFC_RNA_isolation: Sample RNA ID']][which(metadata_freeze$Dx == "SCZ")]
  
  control_samples <- metadata_freeze[['DLPFC_RNA_isolation: Sample RNA ID']][which(metadata_freeze$Dx == "Control")]
  
  bp_samples <- metadata_freeze[['DLPFC_RNA_isolation: Sample RNA ID']][which(metadata_freeze$Dx == "BP")] 
  
  #extract schizophrenia cases only
  cmcExpressionFreezeSCZ <- scale(t(data.matrix(cmcExpressionFreeze[,colnames(cmcExpressionFreeze)%in%scz_samples])))
  
  cmcExpressionFreezeControl <- scale(t(data.matrix(cmcExpressionFreeze[,colnames(cmcExpressionFreeze)%in%control_samples])))
  
  cmcExpressionFreezeBP <- scale(t(data.matrix(cmcExpressionFreeze[,colnames(cmcExpressionFreeze)%in%bp_samples])))
  
  write.csv(cmcExpressionFreezeSCZ,file=sczFile,quote=FALSE)
  write.csv(cmcExpressionFreezeControl,file=controlFile,quote=FALSE)
  write.csv(cmcExpressionFreezeBP,file=bpFile,quote=FALSE)
}

geneSyn <- as.character(commandArgs(TRUE)[[1]])
tfSyn <- as.character(commandArgs(TRUE)[[2]])
sczFile <- as.character(commandArgs(TRUE)[[3]])
controlFile <- as.character(commandArgs(TRUE)[[4]])
bpFile <- as.character(commandArgs(TRUE)[[5]])

makeDataMatrices(geneSyn,tfSyn,sczFile,controlFile,bpFile)
