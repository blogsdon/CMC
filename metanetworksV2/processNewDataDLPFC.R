#login to synapse client
require(synapseClient)
library(dplyr)
library(utilityFunctions)
library(bit64)
synapseLogin()

covariatesObj <- synGet('syn5570053')
residualExpressionObj <- synGet('syn5570051')

library(data.table)

cmcCovariates <- fread(covariatesObj@filePath,data.table=F,stringsAsFactors=F)
cmcExpression <- fread(residualExpressionObj@filePath,data.table=F,stringsAsFactors=F)

cmcExpression2 <- t(cmcExpression[,-c(1)])
colnames(cmcExpression2) <- cmcExpression$GeneFeature
cmcCovariates <- cmcCovariates[,-86]

#Bipolar
BipolarCov <- dplyr::filter(cmcCovariates,Dx=='AFF')

cmcExpression2[dplyr::filter(cmcCovariates,Dx=='AFF')$`DLPFC_RNA_isolation: Sample RNA ID`,] %>%
  apply(2,utilityFunctions::winsorize) %>%
  scale %>%
  write.csv(file='cmcDLPFCBipolarRNAseq.csv',quote=F)

#MCI
cmcExpression2[dplyr::filter(cmcCovariates,Dx=='SCZ')$`DLPFC_RNA_isolation: Sample RNA ID`,] %>%
  apply(2,utilityFunctions::winsorize) %>%
  scale %>%
  write.csv(file='cmcDLPFCSchizophreniaRNAseq.csv',quote=F)

#AD
cmcExpression2[dplyr::filter(cmcCovariates,Dx=='Control')$`DLPFC_RNA_isolation: Sample RNA ID`,] %>%
  apply(2,utilityFunctions::winsorize) %>%
  scale %>%
  write.csv(file='cmcDLPFCControlRNAseq.csv',quote=F)

cmcExpression2 <- apply(cmcExpression2,2,utilityFunctions::winsorize)
cmcExpression2 <- scale(cmcExpression2)

write.csv(cmcExpression2,file='cmcDLPFCRNAseq.csv',quote=F)
