#!/bin/bash

rm codeScz.txt
#rm syn.txt
#rm synSVA.txt

#./grabCMCdata.sh
#echo "syn2757147" >> synSVA.txt
#echo "syn2757149" >> synSVA.txt
#echo "syn2299154" >> synSVA.txt

#echo "syn2299154" >> syn.txt
#echo "syn2757138" >> syn.txt
#echo "syn2757140" >> syn.txt

echo "https://github.com/blogsdon/CMC/blob/master/runCMCcontrols.sh" >> codeScz.txt
echo "https://github.com/blogsdon/CMC/blob/master/grabCMCdata.sh" >> codeScz.txt

cd ../metanetworkSynapse/
qsub -v dataFile="../CMC/sczData.csv",pathv="/shared/metanetworkSynapse/",sparrowZ=1,sparrow2Z=1,lassoCV1se=1,lassoCVmin=1,lassoAIC=1,lassoBIC=1,ridgeCV1se=1,ridgeCVmin=1,ridgeAIC=1,ridgeBIC=1,genie3=1,tigress=1,numberCore=159,outputpath="/shared/CMC/sczNetworks/nosva/" -pe orte 159 buildNet.sh
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/metanetworkSynapse/",sparrowZ=1,sparrow2Z=1,lassoCV1se=1,lassoCVmin=1,lassoAIC=1,lassoBIC=1,ridgeCV1se=1,ridgeCVmin=1,ridgeAIC=1,ridgeBIC=1,genie3=1,tigress=1,numberCore=159,outputpath="/shared/CMC/sczNetworks/sva/" -pe orte 159 buildNet.sh
qsub -v dataFile="../CMC/sczData.csv",pathv="/shared/CMC/sczNetworks/nosva/",aracne=1,aracneFull=1,correlation=1,correlationBonferroni=1,correlationFDR=1,wgcna=1,numberCore=1,outputpath="/shared/CMC/sczNetworks/nosva/" -pe orte 1 buildNet.sh
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/CMC/sczNetworks/sva/",aracne=1,aracneFull=1,correlation=1,correlationBonferroni=1,correlationFDR=1,wgcna=1,numberCore=1,outputpath="/shared/CMC/sczNetworks/sva/" -pe orte 1 buildNet.sh
qsub -v dataFile="../CMC/sczData.csv",pathv="/shared/CMC/sczNetworks/nosva/",aracne=1,aracneFull=1,numberCore=1,outputpath="/shared/CMC/sczNetworks/nosva/" -pe orte 1 buildNet.sh
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/CMC/sczNetworks/sva/",aracne=1,aracneFull=1,numberCore=1,outputpath="/shared/CMC/sczNetworks/sva/" -pe orte 1 buildNet.sh
qsub -v dataFile="../CMC/sczData.csv",pathv="/shared/metanetworkSynapse/",sparrowZ=1,sparrow2Z=1,lassoCVmin=1,ridgeCVmin=1,numberCore=159,outputpath="/shared/CMC/sczNetworks/nosva/" -pe orte 159 buildNet.sh


#new submissions
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/metanetworkSynapse/",tigress=1,numberCore=479,outputpath="/shared/CMC/sczNetworks/sva/" -pe orte 479 buildNet.sh
qsub -v s3="s3://metanetworks/CMC/sczNetworks/nosva/",dataFile="/shared/ROSMAP/sczData.csv",pathv="/shared/metanetworkSynapse/",tigressRootN=1,elasticNetAIC=1,elasticNetBIC=1,elasticNetCVmin=1,elasticNetCV1se=1,numberCore=320,outputpath="/shared/CMC/sczNetworks/nosva/" -pe orte 320 -S /bin/bash -V -cwd -N cmcscz -e /shared/CMC/sczNetworks/nosva/error.txt -o /shared/CMC/sczNetworks/nosva/out.txt buildNet.sh
qsub -v s3="s3://metanetworks/CMC/sczNetworks/sva/",dataFile="/shared/ROSMAP/sczDataSVA.csv",pathv="/shared/metanetworkSynapse/",tigressRootN=1,elasticNetAIC=1,elasticNetBIC=1,elasticNetCVmin=1,elasticNetCV1se=1,numberCore=320,outputpath="/shared/CMC/sczNetworks/sva/" -pe orte 320 -S /bin/bash -V -cwd -N cmcsczSVA -e /shared/CMC/sczNetworks/sva/error.txt -o /shared/CMC/sczNetworks/sva/out.txt buildNet.sh


#qsub -v dataFile="../CMC/controlData.csv",pathv="/shared/metanetworkSynapse/ARACNE/",aracne=1,correlation=1,numberCore=8 -pe orte 8 buildNet.sh
#./pushNet.sh -a "syn3526285" -b "../CMC/code.txt" -c "../CMC/syn.txt" -defghijklmnopq -r "SVA" -s "HomoSapiens" -t "Schizophrenia" -u "DorsolateralPrefrontalCortex"
#./cleanNet.sh

