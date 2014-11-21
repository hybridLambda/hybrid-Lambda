hybridLambdaFst <- function(N, n)  {
## N is number of files
## n is sample size per population
## assuming two populations same sample size
Fst <- c()
#prefix=scan("prefix",what="");
for( i in 1:N ){
  skra <- matrix( scan( paste( c( prefix, "seg_sites/site", i), collapse=""), what=c("",""), quiet=TRUE), 2*n, 2, byrow=T )
  ## see if file has mutations
  if( substring( skra[1,2], 1,1) != "A" ){
    ## at least one mutation in the tree
    ## computing within differences
    Hw <- Hb <- 0
    sites <- 1:nchar( skra[1,2] )
    ## first add to Hw for population A
    for( a1 in 1:(n-1)){
      for( a2 in (a1+1):n){
        seqA1 <- as.numeric( substring( skra[a1,2], sites, sites ) )
        seqA2 <- as.numeric( substring( skra[a2,2], sites, sites ) )
        seqB1 <- as.numeric( substring( skra[(a1+n),2], sites, sites ) )
        seqB2 <- as.numeric( substring( skra[(a2+n),2], sites, sites ) )
        Hw <- Hw  +   sum( abs( seqA1 - seqA2 ) )  +   sum( abs( seqB1 - seqB2 ) ) 
        }
#        print( paste("difference in pop A ",sum( abs( seqA1 - seqA2 ) ) ))
#        print( paste("difference in pop B ",sum( abs( seqB1 - seqB2 ) ) ))
    }
#    print (paste("Hw is ", Hw))
    ## now add to Hw for population B
    #for( a1 in (n+1):(2*n   -1)){
    #  for( a2 in (a1+1):(2*n)){
    #    seq1 <- as.numeric( substring( skra[a1,2], sites, sites ) )
    #    seq2 <- as.numeric( substring( skra[a2,2], sites, sites ) )
    #    Hw <- Hw  +   sum( abs( seq1 - seq2 ) ) }}
    #### done adding to Hw
    #### now adding to  Hb
    for( a1 in 1:n ){
      for( b1 in (n+1):(2*n) ){
        seq1 <- as.numeric( substring( skra[a1,2], sites, sites ) )
        seq2 <- as.numeric( substring( skra[b1,2], sites, sites ) )
        Hb <- Hb  +   sum( abs( seq1 - seq2 ) ) }}
    Fst <- c(Fst,  1 -  (Hw*n/(Hb*(n-1))) ) }
#    print (paste("Hw is ", Hb))

#print (1 -  (Hw*n/(Hb*(n-1))) )
}
##print( Fst )
plot (ecdf(Fst), col="red") 
return( c( mean( Fst[-1]), var( Fst[-1] ) ) )
}
