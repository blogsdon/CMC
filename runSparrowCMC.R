require(synapseClient)
synapseLogin()

NoSVAGeneExpressionSummary <- synGet('syn2757142')
cmcNoSVAGeneSummary <- read.delim(NoSVAGeneExpressionSummary@filePath)

NoSVATFExpressionSummary <- synGet('syn2757144')
cmcNoSVATFSummary <- read.delim(NoSVATFExpressionSummary@filePath)

NoSVAGeneExpression <- synGet('syn2757138')
cmcNoSVAGeneExpression <- read.delim(NoSVAGeneExpression@filePath)

NoSVATFExpression <- synGet('syn2757140')
cmcNoSVATFExpression <- read.delim(NoSVATFExpression@filePath)


SVAGeneExpressionSummary <- synGet('syn2757151')
cmcSVAGeneSummary <- read.delim(SVAGeneExpressionSummary@filePath)

SVATFExpressionSummary <- synGet('syn2757153')
cmcSVATFSummary <- read.delim(SVATFExpressionSummary@filePath)

SVAGeneExpression <- synGet('syn2757147')
cmcSVAGeneExpression <- read.delim(SVAGeneExpression@filePath)

SVATFExpression <- synGet('syn2757149')
cmcSVATFExpression <- read.delim(SVATFExpression@filePath)


sparrow <- function(data,driverIndex,...){
  extractEdges <- function(index,data,driverIndex,...){
    X <- data;
    if(index%%100==0){
      cat(index,'!','\n')
    }
    if (index%in%driverIndex){
      y <- X[,index];
      G <- X[,driverIndex];
      wi <- which(driverIndex%in%index);
      G <- G[,-wi];
      res <- vbsr(y=y,X=G,...);
      pval <- rep(1,length(res$pval)+1);
      pval[-wi] <- res$pval;
      beta <- rep(0,length(res$beta)+1);
      beta[-wi] <- res$beta;
      
    }else{
      y <- X[,index];
      G <- X[,driverIndex];
      res <- vbsr(y=y,X=G,...);
      pval <- res$pval;
      beta <- res$beta;
    }
    return(list(pval=pval,beta=beta));
  }
  ind <- 1:ncol(data);
  edgeMat <- lapply(ind,extractEdges,data=data,driverIndex=driverIndex,...)
  #print(edgeMat)
  edgeMatPval <- sapply(edgeMat,function(x){return(x$pval)});
  edgeMatBeta <- sapply(edgeMat,function(x){return(x$beta)});
  driverMat <- edgeMatPval< 0.05/(ncol(data)*length(driverIndex));
  colnames(driverMat) <- colnames(data);
  rownames(driverMat) <- names(driverIndex)
  colnames(edgeMatBeta) <- colnames(driverMat)
  rownames(edgeMatBeta) <- rownames(driverMat)
  return(list(driverMat=driverMat,edgeMatBeta=edgeMatBeta))
}


sparrow2 <- function(data,driverIndex,...){
  extractEdges <- function(index,data,driverIndex,...){
    X <- data;
    if(index%%100==0){
      cat(index,'!','\n')
    }
    if (index%in%driverIndex){
      y <- X[,index];
      G <- X[,driverIndex];
      wi <- which(driverIndex%in%index);
      G <- G[,-wi];
      res <- vbsr(y=y,X=G,...);
      pval <- rep(1,length(res$pval)+1);
      pval[-wi] <- res$pval;
      beta <- rep(0,length(res$beta)+1);
      beta[-wi] <- res$z;
      
    }else{
      y <- X[,index];
      G <- X[,driverIndex];
      res <- vbsr(y=y,X=G,...);
      pval <- res$pval;
      beta <- res$z;
    }
    return(list(pval=pval,z=beta));
  }
  ind <- 1:ncol(data);
  edgeMat <- lapply(ind,extractEdges,data=data,driverIndex=driverIndex,...)
  #print(edgeMat)
  edgeMatPval <- sapply(edgeMat,function(x){return(x$pval)});
  edgeMatBeta <- sapply(edgeMat,function(x){return(x$z)});
  driverMat <- edgeMatPval< 0.05/(ncol(data)*length(driverIndex));
  colnames(driverMat) <- colnames(data);
  rownames(driverMat) <- names(driverIndex)
  colnames(edgeMatBeta) <- colnames(driverMat)
  rownames(edgeMatBeta) <- rownames(driverMat)
  return(list(driverMat=driverMat,edgeMatBeta=edgeMatBeta))
}

runDriverPathwayEnrichment <- function(){
  #load pathways
  #load network
  #run enrichment
  #run network clustering
}
  
runDriverGWASEnrichment <- function(){
  
}
  
schizophreniaHits <- c('DPYD', 'MIR137', 'ARL3', 'AS3MT', 'C10orf32', 'CNNM2', 'CYP17A1', 'INA', 
                       'NT5C2','PCGF6','PDCD11','SFXN3','TAF5','TRIM8','USMG5','WBP1L','CACNA1C',
                       'TSNARE1','SLC39A8','MAD1L1','ZSWIM6','ABCB9','ARL6IP4','C12orf65','CDK2AP1',
                       'MPHOSPH9','OGFOD2','PITPNM2','RILPL2','SBNO1','SETD8','AC073043.2','C2orf47',
                       'C2orf69','TYW5','FES','FURIN','MAN2A2','TRANK1','AL049840.1','APOPT1','BAG5',
                       'CKB','KLC1','PPP1R13B','TRMT61A','XRCC3','ZFYVE21','AC027228.1','AGPHD1','CHRNA3',
                       'CHRNA5','CHRNB4','IREB2','PSMA4','IMMP2L','SNX19','ZNF804A','CNKSR2','CACNB2',
                       'LRP1','MYO1A','NAB2','NDUFA4L2','NXPH4','R3HDM2','SHMT2','STAC3','STAT6',
                       'TAC3','TMEM194A','LRRIQ3','C2orf82','EFHD1','GIGYF2','KCNJ13','NGEF','ESAM',
                       'MSANTD2','NRGN','VSIG2','TCF4','AMBRA1','ARHGAP1','ATG13','CHRM4','CKAP5','CREB3L1',
                       'DGKZ','F2','HARBI1','MDK','ZNF408','CCDC39','DNAJC19','FXR1','ACTR5','PPP1R16B',
                       'SLC32A1','FANCL','VRK2','ADAMTSL3','GOLGA6L4','ZSCAN2','TCF4','ANKRD44','BOLL','COQ10B',
                       'HSPD1','HSPE1','HSPE1','MARS2','PLCL1','RFTN2','SF3B1','CHADL','EP300','L3MBTL2','RANGAP1',
                       'KCNV1','CNTN4','DRD2','IGSF9B','GLT8D1','GNL3','ITIH1','ITIH3','ITIH4','MUSTN1','NEK4',
                       'NISCH','NT5DC2','PBRM1','SMIM4','SPCS1','STAB1','TMEM110','TMEM110-MUSTN1','ALDOA','ASPHD1',
                       'C16orf92','DOC2A','FAM57B','GDPD3','HIRIP3','INO80E','KCTD13','MAPK3','PPP4C','SEZ6L2','TAOK2',
                       'TBX6','TMEM219','YPEL3','CACNA1I','MSL2','NCK1','PCCB','PPP2R3A','SLC35G2','STAG1','GRIA1','PJA1',
                       'SGSM2','SMG6','SRR','TSR1','GRM3','VPS14C','KDM4A','PTPRF','CILP2','GATAD2A','HAPLN4','MAU2',
                       'NCAN','NDUFA13','PBX4','SUGP1','TM6SF2','TSSK6','ANP32E','APH1A','C1orf51','C1orf54','CA14','OTUD7B',
                       'PLEKHO1','VPS45','SNAP91','PLCH2','ERCC','MLL5','PUS7','SRPK2','RERE','SLC45A1','ATP2A2',
                       'C4orf27','CLCN3','NEK1','FUT9','CENPM','CYP2D6','FAM109B','NAGA','NDUFA6','SEPT3','SHISA8','SMDT1',
                       'SREBF2','TCF20','TNFRSF13C','WBP2NL','BTBD18','C11orf31','CLP1','CTNND1','MED19','SERPING1','TMX2',
                       'YPEL4','ZDHHC5','LUZP2','DGKI','PTN','TLE1','AKT3','SDCCAG8','ANKRD63','PAK6','PLCB2','ZNF536',
                       'MEF2C','TBC1D5','CDC25C','CTNNA1','EGR1','ETF1','FAM53C','GFRA3','HSPA9','KDM3B','REEP2','BCL11B',
                       'AC005477.1','RGS6','HCN1','CA8','CYP26B1','GRAMD1B','SATB2','PCGEM1','GPM6A','CSMD1','CUL3','MMP16','GRIN2A',
                       'PRKD1','ATXN7','C3orf49','PSMD6','THOC7','ACD','C16orf86','CENPT','TRL','DDX28','DPEP2','DPEP3','DUS2L',
                       'EDC4','ENKD1','ESRP2','GFOD2','LCAT','NFATC3','NRN1L','NUTF2','PARD6A','PLA2G15','PSKH','PSMB10','RANBP10',
                       'SLC12A4','SLC7A6','SLC7A6OS','THAP11','TSNAXIP1','EPC2','ATPAF2','DRG2','GID4','LRRC48','MYO15A',
                       'RAI1','SREBF1','TOM1L2','TLE3','CNOT1','SLC38A7','CLU','EPHX2','NLGN4X','RIMS1','DFNA5','MPP6','OSBPL3',
                       'MAN2A1','MIR548AJ2','GALN10','C11orf87','IMMP2L','TMTC1','PODXL','FAM5B','C1orf132','CD46','CR1L',
                       'KCNB1','PTGIS','C12orf79','DPP4','SLC4A10','NOSIP','PRR12','PRRG2','RCN3','RRAS','SCAF1','C12orf42',
                       'AC005609.1','CD14','DND1','HARS','HARS2','IK','NDUFA2','PCDHA1','PCDHA10','PCDHA2','PCDHA3',
                       'PCDHA4','PCDHA5','PCDHA6','PCDHA7','PCDHA8','PCDHA9','TMCO6','WDR55','ZMAT2')





svaMat <- rbind(cmcSVAGeneExpression,cmcSVATFExpression)
noSvaMat <- rbind(cmcNoSVAGeneExpression,cmcNoSVATFExpression)

save(svaMat,noSvaMat,sparrow,file='cmcSparrowDataAndFunction.rda')
library(vbsr)
sdSva <- apply(svaMat,1,sd)
set.seed(7452327)
sparrowSVA <- sparrow(t(svaMat),13838:16423)


set.seed(7452327)
sparrowSVA2 <- sparrow2(t(svaMat),13838:16423)
set.seed(7452327)
sparrowNoSVA2 <- sparrow2(t(noSvaMat),13838:16423)


set.seed(7452327)
sparrowNoSVA <- sparrow(t(noSvaMat),13838:16423)
rownames(sparrowSVA) <- rownames(svaMat)[13838:16423]
rownames(sparrowNoSVA) <- rownames(noSvaMat)[13838:16423]
#rownames(sparrowSVA) <- names(sort(sdSva,decreasing=T)[1:3000])
svaDriver <- rowSums(sparrowSVA)
cmcSVASummary <- rbind(cmcSVAGeneSummary,cmcSVATFSummary)
schizoList <- scan('schizophreniaGeneCards.txt',what='character')






computeFisherEnrichment <- function(list1,list2,list3){
  n2 <- length(list3);
  n1 <- length(list1);
  k <- length(list2);
  a1 <- intersect(list1,list3);
  m <- length(a1);
  n <- n2-m;
  #cat(n2,n1,a1,m,n,'\n')
  
  #s1 <- unlist(lapply(1:n2,getSum,list2,list1));
  s1<-sum(list2%in%list1);
  
  #cat(n2,n1,a1,m,n,s1)
  enr <- (s1/(k))/(m/n2)
  return(enr);
  
}
computeFisherPvalue <- function(list1,list2,list3){
  n2 <- length(list3);
  n1 <- length(list1);
  k <- length(list2);
  a1 <- intersect(list1,list3);
  m <- length(a1);
  n <- n2-m;
  #cat(n2,n1,a1,m,n,'\n')
  
  #s1 <- unlist(lapply(1:n2,getSum,list2,list1));
  s1<-sum(list2%in%list1);
  
  #cat(n2,n1,a1,m,n,s1)
  #model <- list();
  pval <- phyper(q=s1-1,m=m,n=n,k=k,lower.tail=F);
  #model$enr <- (s1/(k))/(m/n2)
  return(pval);
  
}


GOfun <- function(driverMatrix){
  #load msigdb GO annotations
  goGene <- vector('list',1454);
  goName <- rep(0,1454);
  for(i in 1:1454){	
    d1 <- scan('c2.cp.v4.0.symbols.gmt',nlines=1,skip=i-1,what='character');
    goGene[[i]] <- d1[-c(1,2)];
    goName[i] <- d1[1];
  }
  
  #run GO enrichment analysis
  model <- list();
  model$GOpvalMat <- matrix(NA,nrow(driverMatrix),1454);
  model$GOenrMat <- matrix(NA,nrow(driverMatrix),1454);
  
  for (i in 1:nrow(driverMatrix)){
    for (j in 1:1454){
      chrmGene <- colnames(driverMatrix)[driverMatrix[i,]!=0];
      if(length(chrmGene)>0){
        model$GOpvalMat[i,j] <- computeFisherPvalue(goGene[[j]],chrmGene,colnames(driverMatrix))
        model$GOenrMat[i,j] <- computeFisherEnrichment(goGene[[j]],chrmGene,colnames(driverMatrix))
        
      }
      
    }
    print(i)
  }	
  model <- lapply(model,function(x,s){rownames(x)<-s;return(x);},rownames(driverMatrix));
  model <- lapply(model,function(x,s){colnames(x)<-s;return(x);},goName);
  return(model);
}

rownames(sparrowSVA$driverMat) <- cmcSVATFSummary$MAPPED_genes
colnames(sparrowSVA$driverMat) <- c(as.character(cmcSVAGeneSummary$MAPPED_genes),as.character(cmcSVATFSummary$MAPPED_genes))
gores <- GOfun(sparrowSVA$driverMat)

#run sparrow with top 3000 most variable genes




fdrThres <- function (pval, fdr = 0.05) {
  n <- length(pval)
  comp <- sort(pval) < ((fdr/n) * (1:n))
  if (min(pval) < (fdr/n)) {
    w1 <- which(!comp)[1]
    return((fdr/n) * w1)
  }
  else {
    return(fdr/n)
  }
}


sparrowDriver<-as.character(cmcSVATFSummary[as.character(names(sort(svaDriver,decreasing=T))),1]


require(ROCR)


vec1 <- rep(0,length(sparrowDriver))
names(vec1) <- as.character(cmcSVATFSummary[as.character(names(svaDriver)),1])
vec1[names(vec1)%in%schizophreniaHits]<- 1

pred1 <- rowSums(sparrowSVA)
pred2 <- rowSums(sparrowNoSVA)
pred1b <- prediction(pred1,vec1)
pred2b <- prediction(pred2,vec1)
perf1 <- performance(pred1b,'tpr','fpr')
perf2 <- performance(pred2b,'tpr','fpr')
tiff(file='sparrowDriver1.tiff',compression='lzw',height=3.42,width=3.42,pointsize=11,units='in',res = 300)
plot(perf1)
par(new=T);plot(perf2,col='green',main='GWAS Enrichment\n Among SPARROW drivers')
lines(c(0,1),c(0,1),col='red',lwd=2)
legend('bottomright',c('SVA Corrected\nAUC=0.62','No SVA Correction\nAUC=0.52'),lty=1,col=c(1,3),cex=.7)
dev.off()


rownames(sparrowSVA$driverMat) <- rownames(svaMat)[13838:16423]
rownames(sparrowNoSVA$driverMat) <- rownames(svaMat)[13838:16423]


tableFun <- function(x){
  tableFun1 <- which(x,T)
  tableFun2 <- tableFun1
  tableFun2[,1] <- rownames(x)[tableFun1[,1]];
  tableFun2[,2] <- colnames(x)[tableFun1[,2]];
  return(tableFun2)
}

tableFun2 <- function(x,z){
  tableFun1 <- which(x,T)
  tableFun2 <- tableFun1
  tableFun2[,1] <- rownames(x)[tableFun1[,1]];
  tableFun2[,2] <- colnames(x)[tableFun1[,2]];
  y <- rep(0,nrow(tableFun2));
  cat(nrow(tableFun2),'\n')
  
  for(i in 1:nrow(tableFun2)){
    
    y[i] <- z[tableFun1[i,1],tableFun1[i,2]];
  }
  tableFun2 <- cbind(tableFun2,y)
  return(tableFun2)
}


sparrowSVAtable2 <- tableFun(t(driverMat[1:200,]))
rml <- which(sparrowSVAtable2=='.',T)
sparrowSVAtable2 <- sparrowSVAtable2[-rml[,1],]
sparrowSVAt <- tableFun2(sparrowSVA$driverMat,sparrowSVA$edgeMatBeta)[,c(1,3,2)]
sparrowNoSVAt <- tableFun2(sparrowNoSVA$driverMat,sparrowNoSVA$edgeMatBeta)[,c(1,3,2)]


sparrowSVAall <- sparrowSVA$driverMat
sparrowSVAall[1:length(sparrowSVAall)] <- TRUE
sparrowSVAt2 <- tableFun2(sparrowSVAall,sparrowSVA$edgeMatBeta)
sparrowNoSVAall <- sparrowNoSVA$driverMat
sparrowNoSVAall[1:length(sparrowNoSVAall)] <- TRUE
sparrowNoSVAt2 <- tableFun2(sparrowNoSVAall,sparrowNoSVA$edgeMatBeta)

rownames(sparrowSVA2$edgeMatBeta) <- rownames(sparrowSVA$edgeMatBeta)
rownames(sparrowNoSVA2$edgeMatBeta) <- rownames(sparrowNoSVA$edgeMatBeta)
colnames(sparrowSVAt) <- c('Regulator','Weight','Target')
colnames(sparrowNoSVAt) <- c('Regulator','Weight','Target')
colnames(sparrowSVAt2) <- c('Regulator','Weight','Target')
colnames(sparrowNoSVAt2) <- c('Regulator','Weight','Target')
#write.table(sparrowSVAt,file='sva/VBSR_NoRand_Weight.tsv',sep='\t',row.names=F,quote=F)
#write.table(sparrowNoSVAt,file='nosva/VBSR_NoRand_Weight.tsv',sep='\t',row.names=F,quote=F)
write.table(sparrowSVA$edgeMatBeta,file='sva/VBSR_Rand_Weight.tsv',sep='\t',row.names=T,quote=F)

write.table(sparrowSVA$edgeMatBeta,file='nosva/VBSR_NoRand_Weight.tsv',sep='\t',row.names=T,quote=F)

write.table(sparrowSVA2$edgeMatBeta,file='sva/VBSR_Z_Statistic.tsv',sep='\t',row.names=T,quote=F)
write.table(sparrowNoSVA2$edgeMatBeta,file='nosva/VBSR_Z_Statistic.tsv',sep='\t',row.names=T,quote=F)

require(synapseClient)
synapseLogin()
svafolder <- synGet('syn2791415')
svaFile <- File(path=paste(getwd(),'/sva/VBSR_Z_Statistic.tsv',sep=''),parentId='syn2791415')
nosvaFile <- File(path=paste(getwd(),'/nosva/VBSR_Z_Statistic.tsv',sep=''),parentId='syn2791383')
svaFile <- synStore(svaFile)
nosvaFile <- synStore(nosvaFile)


library(igraph)

d1 <- sparrowSVA$driverMat;
d1 <- rbind(d1,matrix(0,16423-2586,16423))
d1[1:16423,1:2586] <- t(sparrowSVA$driverMat)
g1 <- graph.adjacency(d1)

clus <- clusters(g1)
g1 <-delete.vertices(g1,which(clus$membership!=1))
#res <- spinglass.community(g1)
g2 <- as.undirected(g1)
res <- fastgreedy.community(g2)

gggm <- matrix(0,592,max(res$membership))
for (i in 1:max(res$membership)){
  gggm[,i] <- colMeans(svaMat2[res$membership==i,]);

}


s2 <- colnames(sparrowSVA$driverMat)[-which(clus$membership!=1)]

#test

makeDigraph <- function(G,name,graphfile){
  cat('graph',name,'{\n',file=graphfile)
  x <- which(G,T)
  
}

