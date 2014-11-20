#!/bin/bash

prefix="test1"
Nrep=10000
n_sample=20
#echo ${prefix} > prefix
hybrid-Lambda -spcu "(A:1,B:1);" -seg -fst -S ${n_sample} ${n_sample} -num ${Nrep} -o ${prefix} -seed 1

echo "rm(list = ls())
source(\"hybridLambda_fst.r\")
prefix = \"${prefix}\"
hybridLambdaFst( ${Nrep}, ${n_sample} ) 
fst = read.table(\"${prefix}_fst\")\$V1 
mean(fst)
sd(fst)
plot (ecdf(fst), add = TRUE, col = \"blue\") 
" > run.r

R CMD BATCH run.r
tail -15 run.r.Rout 

#prefix="test2"
#Nrep=100
#n_sample=20
##echo ${prefix} > prefix
#hybrid-Lambda -spcu "(A:1,B:1);" -seg -fst -S ${n_sample} ${n_sample} -num ${Nrep} -o ${prefix} -seed 1

#echo "rm(list = ls())
#source(\"hybridLambda_fst.r\")
#prefix = \"${prefix}\"
#hybridLambdaFst( ${Nrep}, ${n_sample} ) 
#" > run.r

#R CMD BATCH run.r
#tail -15 run.r.Rout 
