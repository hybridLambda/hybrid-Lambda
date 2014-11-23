#!/bin/bash

fst_test () {
	#echo ${prefix} > prefix
	hybrid-Lambda -spcu "(A:1,B:1);" -seg -fst -S ${n_sample} ${n_sample} -num ${Nrep}  -seed ${seed}
	
	echo "rm(list = ls())
	source(\"hybridLambda_fst.r\")
	prefix = \"OUT\"
	bjarki_fst = hybridLambdaFst( ${Nrep}, ${n_sample} ) 
	fst = read.table(\"OUT_fst\")\$V1 
	fst = fst[!is.nan(fst)]
	hybridlambda_fst = c( mean(fst), var(fst))
	#bjarki_fst == c( mean(fst), var(fst))
	bjarki_fst
	hybridlambda_fst
	result_bool = (bjarki_fst - hybridlambda_fst) < 1e-6
	
	keyword = \"Not \"
	if ( result_bool[1] && result_bool[2] ) 
	keyword = \"\"
	
	print ( paste ( \"Status \", keyword, \"Ok\", sep = \"\"))
	
	#(bjarki_fst - hybridlambda_fst)/ hybridlambda_fst < 1e-7
	#plot (ecdf(fst), add = TRUE, col = \"blue\") 
	" > run.r
	
	R CMD BATCH run.r 
	grep "Status Ok" run.r.Rout 	
}

for i in $(seq 1 9); do 
	Nrep=${i}00
	n_sample=2${i}
	seed=${i}
	fst_test || exit 1
done

for i in $(seq 1 9); do 
	Nrep=${i}00
	n_sample=1${i}
	seed=${i}
	fst_test || exit 1
done

for i in $(seq 1 9); do 
	Nrep=10${i}00
	n_sample=3${i}
	seed=${i}
	fst_test || exit 1
done
