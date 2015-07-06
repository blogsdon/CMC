#!/bin/sh

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N scz_sva_cmc_networks
#$ -e error_scz_sva.txt
#$ -o out_scz_sva.txt


#cases
#/shared/metanetworkSynapse/pushNet.sh -a "syn4549880" -b "/shared/CMC/codeScz.txt" -c "/shared/CMC/synSVA.txt" -defghijklmnopqv -r "SVA" -s "HomoSapiens" -t "Schizophrenia" -u "DLPFC" -x "/shared/metanetworkSynapse/pushNetworkSynapse.R"
/shared/metanetworkSynapse/pushNet.sh -a "syn4549880" -b "/shared/CMC/codeScz.txt" -c "/shared/CMC/synSVA.txt" -defghijklmnopq -r "SVA" -s "HomoSapiens" -t "Schizophrenia" -u "DLPFC" -x "/shared/metanetworkSynapse/pushSparseNetworkSynapse.R"