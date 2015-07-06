#!/bin/sh

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N scz_cmc_networks
#$ -e error_scz.txt
#$ -o out_scz.txt


#cases
#/shared/metanetworkSynapse/pushNet.sh -a "syn3526289" -b "/shared/CMC/codeScz.txt" -c "/shared/CMC/syn.txt" -defghijklmnopqv -r "None" -s "HomoSapiens" -t "Schizophrenia" -u "DLPFC" -x "/shared/metanetworkSynapse/pushNetworkSynapse.R"
/shared/metanetworkSynapse/pushNet.sh -a "syn3526289" -b "/shared/CMC/codeScz.txt" -c "/shared/CMC/syn.txt" -defghijklmnopq -r "None" -s "HomoSapiens" -t "Schizophrenia" -u "DLPFC" -x "/shared/metanetworkSynapse/pushSparseNetworkSynapse.R"