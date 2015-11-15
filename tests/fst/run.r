rm(list = ls())
	source("tests/fst/hybridLambda_fst.r")
	prefix = "OUT"
	bjarki_fst = hybridLambdaFst( 100, 21 ) 
	fst = read.table("OUT_fst")$V1 
	fst = fst[!is.nan(fst)]
	hybridlambda_fst = c( mean(fst), var(fst))
	#bjarki_fst == c( mean(fst), var(fst))
	bjarki_fst
	hybridlambda_fst
	result_bool = (bjarki_fst - hybridlambda_fst) < 1e-6
	
	keyword = "Not "
	if ( result_bool[1] && result_bool[2] ) 
	keyword = ""
	
	print ( paste ( "Status ", keyword, "Ok", sep = ""))
	
	#(bjarki_fst - hybridlambda_fst)/ hybridlambda_fst < 1e-7
	#plot (ecdf(fst), add = TRUE, col = "blue") 
	
