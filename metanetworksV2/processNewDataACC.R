#login to synapse client
require(synapseClient)
library(dplyr)
library(utilityFunctions)
library(bit64)
synapseLogin()

covariatesObj <- synGet('syn5570059')
residualExpressionObj <- synGet('syn5570057')

library(data.table)

cmcCovariates <- fread(covariatesObj@filePath,data.table=F,stringsAsFactors=F)
cmcExpression <- fread(residualExpressionObj@filePath,data.table=F,stringsAsFactors=F)

cmcExpression2 <- t(cmcExpression[,-c(1)])
colnames(cmcExpression2) <- cmcExpression$GeneFeature
cmcCovariates <- cmcCovariates[,-136]

#Bipolar
#BipolarCov <- dplyr::filter(cmcCovariates,Dx=='AFF')

cmcExpression2[dplyr::filter(cmcCovariates,Dx=='AFF')$`ACC_RNA_isolation: Sample RNA ID`,] %>%
  apply(2,utilityFunctions::winsorize) %>%
  scale %>%
  write.csv(file='cmcACCBipolarRNAseq.csv',quote=F)

#MCI
cmcExpression2[dplyr::filter(cmcCovariates,Dx=='SCZ')$`ACC_RNA_isolation: Sample RNA ID`,] %>%
  apply(2,utilityFunctions::winsorize) %>%
  scale %>%
  write.csv(file='cmcACCSchizophreniaRNAseq.csv',quote=F)

#AD
cmcExpression2[dplyr::filter(cmcCovariates,Dx=='Control')$`ACC_RNA_isolation: Sample RNA ID`,] %>%
  apply(2,utilityFunctions::winsorize) %>%
  scale %>%
  write.csv(file='cmcACCControlRNAseq.csv',quote=F)

cmcExpression2 <- apply(cmcExpression2,2,utilityFunctions::winsorize)
cmcExpression2 <- scale(cmcExpression2)

write.csv(cmcExpression2,file='cmcACCRNAseq.csv',quote=F)
