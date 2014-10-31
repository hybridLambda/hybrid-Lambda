/*
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010 -- 2014 Sha (Joe) Zhu
 * 
 * This file is part of hybrid-Lambda.
 * 
 * hybrid-Lambda is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "sim_gt.hpp"
#include "fst.hpp"
#include <cassert>
using namespace std;


//function(N, n)
  //{
    //// N is number of files
    //// n is sample size per population
    //// assuming two populations same sample size
    //Fst <- 0
    //for( i in 1:N ){
      //skra <- matrix( scan( paste( c("/home/knoppix/verk/Rvinna/cprograms/largesamplesize/jerome/dist/seg-sites/site", i), collapse=""), what=c("",""), quiet=TRUE), 2*n, 2, byrow=T )
      //// see if file has mutations
      //if( substring( skra[1,2], 1,1) != "A" ){
        //// at least one mutation in the tree
        //// computing within differences
        //Hw <- Hb <- 0
        //sites <- 1:nchar( skra[1,2] )
        //// first add to Hw for population A
        //for( a1 in 1:(n-1)){
          //for( a2 in (a1+1):n){
            //seqA1 <- as.numeric( substring( skra[a1,2], sites, sites ) )
            //seqA2 <- as.numeric( substring( skra[a2,2], sites, sites ) )
            //seqB1 <- as.numeric( substring( skra[(a1+n),2], sites, sites ) )
            //seqB2 <- as.numeric( substring( skra[(a2+n),2], sites, sites ) )
            //Hw <- Hw  +   sum( abs( seqA1 - seqA2 ) )  +   sum( abs( seqB1 - seqB2 ) ) }}
        //// now add to Hw for population B
        //#for( a1 in (n+1):(2*n   -1)){
        //#  for( a2 in (a1+1):(2*n)){
        //#    seq1 <- as.numeric( substring( skra[a1,2], sites, sites ) )
        //#    seq2 <- as.numeric( substring( skra[a2,2], sites, sites ) )
        //#    Hw <- Hw  +   sum( abs( seq1 - seq2 ) ) }}
        ////// done adding to Hw
        ////// now adding to  Hb
        //for( a1 in 1:n ){
          //for( b1 in (n+1):(2*n) ){
            //seq1 <- as.numeric( substring( skra[a1,2], sites, sites ) )
            //seq2 <- as.numeric( substring( skra[b1,2], sites, sites ) )
            //Hb <- Hb  +   sum( abs( seq1 - seq2 ) ) }}
        //Fst <- c(Fst,  1 -  (Hw*n/(Hb*(n-1))) ) }
    //}
    ////print( Fst )
    //return( c( mean( Fst[-1]), var( Fst[-1] ) ) )
  //}
