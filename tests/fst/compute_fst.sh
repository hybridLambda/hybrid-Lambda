#!/bin/bash

prefix="test1"
Nrep=10
n_sample=2
#echo ${prefix} > prefix
hybrid-Lambda -spcu "(A:1,B:1);" -seg -fst -S ${n_sample} ${n_sample} -mu 0.01 -num ${Nrep} -o ${prefix}

echo "rm(list = ls())
source(\"hybridLambda_fst.r\")
prefix = \"${prefix}\"
hybridLambdaFst( ${Nrep}, ${n_sample} ) 
" > run.r

R CMD BATCH run.r
