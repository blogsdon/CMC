#!/bin/sh

#scz no sva
echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546233" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N aracne -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546344" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N correlation -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546964" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N genie3 -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546464" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N lassoAIC -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546491" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N lassoBIC -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546372" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N lassoCV1se -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546422" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N lassoCVmin -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546863" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N ridgeAIC -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546865" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N ridgeBIC -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546597" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N ridgeCV1se -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546777" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N ridgeCVmin -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546065" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N sparrow1 -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546164" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N sparrow2 -e error.txt -o out.txt test.sh

echo "Rscript /shared/metanetworkSynapse/makeMultiSparseNetwork.R "syn4546999" "syn4546994" "syn4550165" "syn3526289" "executed.txt"" > test.sh
qsub -S /bin/bash -V -cwd -N tigress -e error.txt -o out.txt test.sh

#scz sva

#control no sva

#control sva