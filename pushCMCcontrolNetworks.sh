#!/bin/sh

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N control_cmc_networks
#$ -e error_control.txt
#$ -o out_control.txt
#$ -pe orte 8


#controls
#/shared/metanetworkSynapse/pushNet.sh -a "syn3526290" -b "/shared/CMC/codeControl.txt" -c "/shared/CMC/syn.txt" -defghijklmnopqv -r "None" -s "HomoSapiens" -t "Control" -u "DLPFC" -x "/shared/metanetworkSynapse/pushNetworkSynapse.R"
#/shared/metanetworkSynapse/pushNet.sh -a "syn3526290" -b "/shared/CMC/codeControl.txt" -c "/shared/CMC/syn.txt" -defghijklmnopq -r "None" -s "HomoSapiens" -t "Control" -u "DLPFC" -x "/shared/metanetworkSynapse/pushSparseNetworkSynapse.R"
/shared/metanetworkSynapse/pushNet.sh -a "syn3526290" -b "/shared/CMC/codeControl.txt" -c "/shared/CMC/syn.txt" -z -r "None" -s "HomoSapiens" -t "Control" -u "DLPFC" -x "/shared/metanetworkSynapse/pushSparseNetworkSynapse.R"
