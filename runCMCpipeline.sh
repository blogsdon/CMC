#!/bin/bash



#./grabCMCdata.sh
echo "syn2757147" >> synSVA.txt
echo "syn2757149" >> synSVA.txt
echo "syn2299154" >> synSVA.txt

echo "syn2299154" >> syn.txt
echo "syn2757138" >> syn.txt
echo "syn2757140" >> syn.txt

echo "https://github.com/blogsdon/CMC/blob/master/runCMCpipeline.sh" >> code.txt
echo "https://github.com/blogsdon/CMC/blob/master/grabCMCdata.sh" >> code.txt

cd ../metanetworkSynapse/
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/metanetworkSynapse/",sparrow=1,aracne=0,wgcna=0,lasso=1,ridge=1,genie3=1,tigress=1,numberCore=160 -pe orte 160 buildNet.sh
qsub -v dataFile="../CMC/sczDataSVA.csv",pathv="/shared/metanetworkSynapse/ARACNE/",sparrow=0,aracne=1,wgcna=1,lasso=0,ridge=0,genie3=0,tigress=0,numberCore=8 -pe orte 8 buildNet.sh
./pushNet.sh -p "syn3526285" -c "../CMC/code.txt" -b "../CMC/synSVA.txt" -swlrgta -v "SVA" -o "HomoSapiens" -d "Schizophrenia" -u "DorsolateralPrefrontalCortex"
#./cleanNet.sh

