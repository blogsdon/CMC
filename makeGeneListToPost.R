makeGeneListToPost <- function(gene_de,tf_de,syn1,name1,disease,testNet){
  require(synapseClient)
  synapseLogin()
  gene_de_obj <- synGet(gene_de);
  tf_de_obj <- synGet(tf_de);
  
  gene_de_obj2 <- read.delim(gene_de_obj@filePath,stringsAsFactors=F)
  tf_de_obj2 <- read.delim(tf_de_obj@filePath,stringsAsFactors=F)
  
  n2 <- c(rownames(tf_de_obj2),rownames(gene_de_obj2))
  
  load(testNet)
  n1 <- rownames(network)
  if(sum(n1!=n2)!=0){
    stop('error will robinson!!')
  }
  
  n3 <- c(tf_de_obj2$MAPPED_genes,gene_de_obj2$MAPPED_genes)
  n4 <- as.character(1:length(n3))
  nameKey <- cbind(n2,n3,n4)
  colnames(nameKey) <- c('ensemblId','geneName','variableNumber')
  #name1 <- 'geneNameKeyNoSVA.csv'
  write.csv(nameKey,quote=F,row.names=F,file=name1)
  #syn1 <- 'syn3526290'
  synObj <- File(name1,parentId=syn1)
  anno <- list(fileType='csv',
               dataType='metaData',
               disease=disease,
               organism='HomoSapiens')
  synSetAnnotations(synObj) <- anno;
  act <- Activity(name='Gene Key for Networks',used=as.list(c(gene_de,tf_de)),executed = as.list(c('https://github.com/blogsdon/CMC/blob/master/makeGeneListToPost.R')))
  act <- storeEntity(act)
  generatedBy(synObj) <- act
  
  synObj <- synStore(synObj)
}

#control
#makeGeneListToPost('syn2757142','syn2757144','syn3526290','geneNameKeyNoSVA.csv','Control','/shared/CMC/controlNetworks/nosva/result_lassoBIC.rda')


#control sva
#makeGeneListToPost('syn2757151','syn2757153','syn3526286','geneNameKeySVA.csv','Control','/shared/CMC/controlNetworks/sva/result_lassoBIC.rda')

#scz
#makeGeneListToPost('syn2757142','syn2757144','syn3526289','geneNameKeyNoSVA.csv','Schizophrenia','/shared/CMC/sczNetworks/nosva/result_lassoBIC.rda')


#scz sva
makeGeneListToPost('syn2757151','syn2757153','syn4549880','geneNameKeySVA.csv','Schizophrenia','/shared/CMC/sczNetworks/sva/result_lassoBIC.rda')
