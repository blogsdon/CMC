#!/bin/bash

#number of cores to reserve for job
nthreads=319

#full s3 path where networks will go
s3="s3://metanetworks/CMC/metanetworksV2/Full/ACC/"

#location of data file
dataFile="/shared/CMC/metanetworksV2/cmcACCRNAseq.csv"

#location of metanetwork synapse scripts
pathv="/shared/metanetworkSynapse/"

#output path for temporary result file prior to pushing to s3/synapse
outputpath="/shared/CMC/metanetworksV2/Full/ACC/"

#path within s3
s3b="CMC/metanetworksV2/Full/ACC"

#id of folder on Synapse that files will go to
parentId="syn5816032"

#path to csv file with annotations to add to file on Synapse
annotationFile="/shared/CMC/metanetworksV2/Full/ACC/annoFile.txt"

#path to csv file with provenance to add to file on synapse
provenanceFile="/shared/CMC/metanetworksV2/Full/ACC/provenanceFileRegression.txt"

#path to error output
errorOutput="/shared/CMC/metanetworksV2/Full/ACC/Regressionerror.txt"

#path to out output
outOutput="/shared/CMC/metanetworksV2/Full/ACC/Regressionout.txt"

#job script name
jobname="CMCACCRegressionFull"


echo "qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,c3net=0,mrnet=0,wgcnaTOM=0,sparrowZ=1,lassoCV1se=1,ridgeCV1se=1,genie3=1,tigress=1,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe orte $nthreads -S /bin/bash -V -cwd -N $jobname -e $errorOutput -o $outOutput $pathv/buildNet.sh"

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,c3net=0,mrnet=0,wgcnaTOM=0,sparrowZ=1,lassoCV1se=1,ridgeCV1se=1,genie3=1,tigress=1,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe orte $nthreads -S /bin/bash -V -cwd -N $jobname -e $errorOutput -o $outOutput $pathv/buildNet.sh
