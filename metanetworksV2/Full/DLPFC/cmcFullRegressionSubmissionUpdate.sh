#!/bin/bash

#number of cores to reserve for job
nthreads=319

#full s3 path where networks will go
s3="s3://metanetworks/CMC/metanetworksV2/Full/DLPFC/"

#location of data file
dataFile="/shared/CMC/metanetworksV2/cmcDLPFCRNAseq.csv"

#location of metanetwork synapse scripts
pathv="/shared/metanetworkSynapse/"

#output path for temporary result file prior to pushing to s3/synapse
outputpath="/shared/CMC/metanetworksV2/Full/DLPFC/"

#path within s3
s3b="CMC/metanetworksV2/Full/DLPFC"

#id of folder on Synapse that files will go to
parentId="syn5816035"

#path to csv file with annotations to add to file on Synapse
annotationFile="/shared/CMC/metanetworksV2/Full/DLPFC/annoFile.txt"

#path to csv file with provenance to add to file on synapse
provenanceFile="/shared/CMC/metanetworksV2/Full/DLPFC/provenanceFileRegressionUpdate.txt"

#path to error output
errorOutput="/shared/CMC/metanetworksV2/Full/DLPFC/Regressionerror.txt"

#path to out output
outOutput="/shared/CMC/metanetworksV2/Full/DLPFC/Regressionout.txt"

#job script name
jobname="CMCDLPFCRegressionFull"


echo "qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,c3net=0,mrnet=0,wgcnaTOM=0,sparrowZ=1,lassoCV1se=1,ridgeCV1se=1,genie3=1,tigress=1,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe orte $nthreads -S /bin/bash -V -cwd -N $jobname -e $errorOutput -o $outOutput $pathv/buildNet.sh"

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,c3net=0,mrnet=0,wgcnaTOM=0,sparrowZ=0,sparrow2Z=1,lassoCVmin=1,lassoAIC=1,lassoBIC=1,lassoCV1se=0,ridgeCV1se=0,ridgeCVmin=1,ridgeAIC=1,ridgeBIC=1,genie3=0,tigress=0,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe orte $nthreads -S /bin/bash -V -cwd -N $jobname -e $errorOutput -o $outOutput $pathv/buildNet.sh
