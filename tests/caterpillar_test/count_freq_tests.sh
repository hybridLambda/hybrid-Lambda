#!/bin/bash
# all the trees are caterpillar trees
numtree=100000

for (( numtax=3;numtax<=8;numtax++))
    do
    method=kingman
    hybrid-Lambda -spcu ${numtax}tax_caterpillar -num ${numtree} -o ${numtax}gt${method} -f 
    
    method=alpha
    alpha=1.4
    hybrid-Lambda -spcu ${numtax}tax_caterpillar -num ${numtree} -mm ${alpha} -o ${numtax}gt${method} -f 
    
    method=psi
    psi=0.1
    hybrid-Lambda -spcu ${numtax}tax_caterpillar -num ${numtree} -mm ${psi} -o ${numtax}gt${method} -f 
    done

