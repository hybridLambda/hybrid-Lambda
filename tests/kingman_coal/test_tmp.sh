#!/bin/bash

dir=test-demo
mkdir ${dir}
cd ${dir}
rm *pdf


rep=10000

## compare population sturture for a single population data
COMPAREFILE=compareDemo
rm ${COMPAREFILE}

#source chisq_r.src

source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	

#case 00
net=Pair
echo case00_${net} > current_case
rm ms* hybridLambda*
ms 2 ${rep} -T -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu "(A:0,B:0);" -S 1 1 -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num


rearrange_hybridLambdaout 2

foo


##case 0
#net=Onepop
#echo case0_${net} > current_case
#rm ms* hybridLambda*
#ms 5 ${rep} -T -t 10 > msout
#grep ";" msout > ms_gt
#cat msout | sample_stats > ms_stats
#hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
#hybrid-Lambda -spcu "(A:0,B:0);" -S 4 1 -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num


#rearrange_hybridLambdaout 5

#foo

#case 1
#net=5_tax_sp_nt1_para00_bl06
#echo case1_${net} > current_case
#rm ms* hybridLambda*
##ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -es 1.1 3 1.0 -ej 1.4 6 5 -es 1.4 3 1.0 -ej 1.7 7 5 -ej 1.7 3 1 -ej 2 5 1 -t 10 > msout
##grep ";" msout > ms_gt
#cat msout | sample_stats > ms_stats
#ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -ej 1.7 3 1 -ej 2 5 1 -t 10 > msout
#grep ";" msout > ms_gt
#cat msout | sample_stats > ms_stats
#hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
#hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

#rearrange_hybridLambdaout 5

#foo
