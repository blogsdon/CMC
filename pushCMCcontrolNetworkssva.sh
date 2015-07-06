#!/bin/sh

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N control_sva_cmc_networks
#$ -e error_control_sva.txt
#$ -o out_control_sva.txt


#controls
#/shared/metanetworkSynapse/pushNet.sh -a "syn3526286" -b "/shared/CMC/codeControl.txt" -c "/shared/CMC/synSVA.txt" -defghijklmnopqv -r "SVA" -s "HomoSapiens" -t "Control" -u "DLPFC" -x "/shared/metanetworkSynapse/pushNetworkSynapse.R"
/shared/metanetworkSynapse/pushNet.sh -a "syn3526286" -b "/shared/CMC/codeControl.txt" -c "/shared/CMC/synSVA.txt" -defghijklmnopq -r "SVA" -s "HomoSapiens" -t "Control" -u "DLPFC" -x "/shared/metanetworkSynapse/pushSparseNetworkSynapse.R"

