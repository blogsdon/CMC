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

microglia <- filter(sczSpecificEnrichments,GeneSetName=='Zhang:Microglia')
endothelial <- filter(sczSpecificEnrichments,GeneSetName=='Zhang:Endothelial')
astrocyte <- filter(sczSpecificEnrichments,GeneSetName=='Zhang:Astrocyte')
myelinOligos <- filter(sczSpecificEnrichments,GeneSetName=='Zhang:MyelinOligos')
neuron <- filter(sczSpecificEnrichments,GeneSetName=='Zhang:Neuron')
opc <- filter(sczSpecificEnrichments,GeneSetName=='Zhang:OPC')
black <- filter(sczSpecificEnrichments,ComparisonName=='black')


sczOverall <- list(network=sparseNetwork,
                   modules=sczModules,
                   enrichments=sczSpecificEnrichments)
makeFilesForCytoscape(sczOverall,'sczNetwork.csv','sczModules.csv')


####quick PGC2
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

namedModule <- fread('sczModules.csv') %>% data.frame
library(metanetwork)

enrichRes <- lapply(names(table(namedModule$modulelabels)),function(x) enrichment(schizophreniaHits,namedModule$GeneIDs[namedModule$modulelabels==x],namedModule$GeneIDs))

pvals <- sapply(1:length(enrichRes),function(i,x) return(x[[i]]$pval),enrichRes)
library(gap)
qqunif(pvals)
