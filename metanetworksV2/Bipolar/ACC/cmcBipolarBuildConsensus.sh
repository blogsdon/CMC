#!/bin/sh
#number of cores to reserve for job
nthreads=8

#full s3 path where networks will go
s3="s3://metanetworks/CMC/metanetworksV2/Bipolar/ACC/"

#location of data file
dataFile="/shared/CMC/metanetworksV2/cmcACCBipolarRNAseq.csv"

#location of metanetwork synapse scripts
pathv="/shared/metanetworkSynapse/"

#output path for temporary result file prior to pushing to s3/synapse
outputpath="/local/CMC/metanetworksV2/Bipolar/ACC/"

#path within s3
s3b="CMC/metanetworksV2/Bipolar/ACC"

#id of folder with networks to combine
networkFolderId="syn5816026"

#id of folder on Synapse that file will go to
parentId="syn5816026"

#path to csv file with annotations to add to file on Synapse
annotationFile="/shared/CMC/metanetworksV2/Bipolar/ACC/annoFile.txt"

provenanceFile="/shared/CMC/metanetworksV2/Bipolar/ACC/provenanceFile.txt"

#path to error output
errorOutput="/shared/CMC/metanetworksV2/Bipolar/ACC/Aggregationerror.txt"

#path to out output
outOutput="/shared/CMC/metanetworksV2/Bipolar/ACC/Aggregationout.txt"

#job script name
jobname="CMCbipolarAccaggregation"

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile,networkFolderId=$networkFolderId -pe orte $nthreads -S /bin/bash -V -cwd -N $jobname -e $errorOutput -o $outOutput $pathv/buildConsensus.sh
