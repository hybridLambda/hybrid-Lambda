#!/bin/bash

dir=test-demo
mkdir ${dir}
cd ${dir}
rm *pdf


rep=100000

## compare population sturture for a single population data
COMPAREFILE=compareDemo
rm ${COMPAREFILE}

#source chisq_r.src

#source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	

#case 0
net=Onepop
echo case0_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu "(A:0,B:0);" -S 4 1 -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 1
net=5_tax_sp_nt1_para00_bl06
echo case1_${net} > current_case
rm ms* hybridLambda*
#ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -es 1.1 3 1.0 -ej 1.4 6 5 -es 1.4 3 1.0 -ej 1.7 7 5 -ej 1.7 3 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -ej 1.7 3 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 2
net=5_tax_sp_nt1_para02_bl06
echo case2_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -es 1.1 3 0.8 -ej 1.4 6 5 -es 1.4 3 0.8 -ej 1.7 7 5 -ej 1.7 3 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 3
net=5_tax_sp_nt1_para05_bl06
echo case3_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -es 1.1 3 0.5 -ej 1.4 6 5 -es 1.4 3 0.5 -ej 1.7 7 5 -ej 1.7 3 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 4
net=5_tax_sp_nt1_para08_bl06
echo case4_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -es 1.1 3 0.2 -ej 1.4 6 5 -es 1.4 3 0.2 -ej 1.7 7 5 -ej 1.7 3 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 5
net=5_tax_sp_nt1_para10_bl06
echo case5_${net} > current_case
rm ms* hybridLambda*
#ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -es 1.1 3 0.0 -ej 1.4 6 5 -es 1.4 3 0.0 -ej 1.7 7 5 -ej 1.7 3 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 0.8 2 3 -ej 1.4 3 5 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 6
net=5_tax_sp_nt1_para_bl100
echo case6_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 50.5 2 3 -es 100.5 3 0.5 -ej 150.5 6 5 -es 150.5 3 0.5 -ej 200.5 7 5 -ej 200.5 3 1 -ej 250.5 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 7
net=5_tax_sp_nt1_para_bl0
echo case7_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0 4 3 -ej 0 2 3 -es 0 3 0.5 -ej 0 6 5 -es 0 3 0.5 -ej 0 7 5 -ej 0 3 1 -ej 0 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

##ms 5 10000 -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej 1.05 2 3 -es 2.25 3 0.5 -ej 3.25 6 5 -es 3.3 3 0.5 -ej 4.8 7 5 -ej 4.95 3 1 -ej 5.45 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt

###########################################################################################################
#5taxnt2

#case 8
net=5_tax_sp_nt2_para00_bl06
echo case8_${net} > current_case
rm ms* hybridLambda*
#ms 5 10000 -T -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.0 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -ej 0.8 3 1 -ej 1.1 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 9
net=5_tax_sp_nt2_para02_bl06
echo case9_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.2 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case  10
net=5_tax_sp_nt2_para05_bl06
echo case10_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.5 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 11
net=5_tax_sp_nt2_para08_bl06
echo case11_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.8 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 12
net=5_tax_sp_nt2_para10_bl06
echo case12_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 1.0 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 13
net=5_tax_sp_nt2_para_bl0
echo case13_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0 4 3 -ej 0 2 3 -es 0 3 0.5 -ej 0 3 5 -ej 0 6 1 -ej 0 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 14
net=5_tax_sp_nt2_para_bl100
echo case14_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es .5 3 0.5 -ej 50.5 3 5 -ej 50.5 6 1 -ej 100.5 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo


#case 15
net=5_tax_sp_nt3_para00_bl06
echo case15_${net} > current_case
rm ms* hybridLambda*
#ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0 -ej 1.1 3 5 -es 1.1 6 0 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0 -ej 1.1 3 5 -es 1.1 6 0 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 16
net=5_tax_sp_nt3_para02_bl06
echo case16_${net} > current_case
rm ms* hybridLambda*

ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0.2 -ej 1.1 3 5 -es 1.1 6 0.2 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 17
net=5_tax_sp_nt3_para05_bl06
echo case17_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0.5 -ej 1.1 3 5 -es 1.1 6 0.5 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 18
net=5_tax_sp_nt3_para08_bl06
echo case18_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0.8 -ej 1.1 3 5 -es 1.1 6 0.8 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 19
net=5_tax_sp_nt3_para10_bl06
echo case19_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -ej 1.1 3 5 -ej 1.7 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 20
net=5_tax_sp_nt3_para_bl0
echo case20_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej 0 4 3 -ej 0 2 3 -es 0 3 0.5 -ej 0 3 5 -es 0 6 0.5 -ej 0 6 5 -ej 0 7 1 -ej 0 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 21
net=5_tax_sp_nt3_para_bl100
echo case21_${net} > current_case
rm ms* hybridLambda*
ms 5 ${rep} -T -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 50.5 3 0.5 -ej 100.5 3 5 -es 100.5 6 0.5 -ej 150.5 6 5 -ej 150.5 7 1 -ej 200.5 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo


###########################################################################################################
#6taxnt1
#ms 6 100000 -T -I 6 1 1 1 1 1 1 -ej 2.8 6 2 -ej 2.85 2 3 -ej 3.35 4 3 -es 4.35 3 0.5 -ej 5.405 3 5 -ej 5.355 7 1 -ej 5.455 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt

#case 22
net=6_tax_sp_nt1_para00_bl06
echo case22_${net} > current_case
rm ms* hybridLambda*
#ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.5 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -ej 1.7 3 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 23
net=6_tax_sp_nt1_para02_bl06
echo case23_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.2 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 24
net=6_tax_sp_nt1_para05_bl06
echo case24_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.5 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 25
net=6_tax_sp_nt1_para08_bl06
echo case25_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.8 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 26
net=6_tax_sp_nt1_para10_bl06
echo case26_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 1 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 27
net=6_tax_sp_nt1_para_bl0
echo case27_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0 6 2 -ej 0 2 3 -ej 0 4 3 -es 0 3 .5 -ej 0 3 5 -ej 0 7 1 -ej 0 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 28
net=6_tax_sp_nt1_para_bl100
echo case28_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 50.5 2 3 -ej 100.5 4 3 -es 150.5 3 .5 -ej 200.5 3 5 -ej 200.5 7 1 -ej 250.5 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo


###########################################################################################################
#6taxnt2
#ms 6 100000 -T -I 6 1 1 1 1 1 1 -ej 1.85 2 3 -ej 2.35 4 3 -es 2.35 3 0.5 -ej 1.15 6 5 -ej 2.855 5 3 -ej 3.605 7 1 -ej 3.705 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt

#case 29
net=6_tax_sp_nt2_para00_bl106
echo case29_${net} > current_case
rm ms* hybridLambda*
#ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.0 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -ej 0.8 6 5 -ej 1.1 3 1 -ej 1.4 5 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 30
net=6_tax_sp_nt2_para02_bl106
echo case30_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.2 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 31
net=6_tax_sp_nt2_para05_bl106
echo case31_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.5 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 32
net=6_tax_sp_nt2_para08_bl106
echo case32_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.8 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 33
net=6_tax_sp_nt2_para10_bl106
echo case33_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 1 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 34
net=6_tax_sp_nt2_para_bl0
echo case34_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0 2 3 -ej 0 4 3 -es 0 3 .5 -ej 0 6 5 -ej 0 5 3 -ej 0 7 1 -ej 0 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 35
net=6_tax_sp_nt2_para_bl100
echo case35_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 50.5 4 3 -es 50.5 3 .5 -ej 50.5 6 5 -ej 100.5 5 3 -ej 100.5 7 1 -ej 150.5 3 1 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

###########################################################################################################
#6taxnt3
#ms 6 100000 -T -I 6 1 1 1 1 1 1 -ej 0.85 2 3 -ej 0.85 6 3 -ej 0.85 4 3 -es 1.85 3 0.34 -ej 2.605 3 5 -ej 2.855 7 1 -ej 3.155 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt

#case 36
net=6_tax_sp_nt3_para00_bl06
echo case36_${net} > current_case
rm ms* hybridLambda*
#ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.0 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -ej 1.1 3 1 -ej 1.4 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 37
net=6_tax_sp_nt3_para02_bl06
echo case37_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.2 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 38
net=6_tax_sp_nt3_para05_bl06
echo case38_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.5 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 39
net=6_tax_sp_nt3_para08_bl06
echo case39_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.8 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 40
net=6_tax_sp_nt3_para10_bl06
echo case40_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 1 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 41
net=6_tax_sp_nt3_para_bl0
echo case41_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0 2 3 -ej 0 6 3 -ej 0 4 3 -es 0 3 .5 -ej 0 3 5 -ej 0 7 1 -ej 0 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 42
net=6_tax_sp_nt3_para_bl100
echo case42_${net} > current_case
rm ms* hybridLambda*
ms 6 ${rep} -T -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 50.5 3 .5 -ej 100.5 3 5 -ej 100.5 7 1 -ej 150.5 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

##########################################################################################################

#7taxnt1
#ms 7 100000 -T -I 7 1 1 1 1 1 1 1 -ej 0.05 6 7 -ej 0.25 7 2 -ej 0.8 2 3 -ej 2.45 4 3 -es 3.45 3 0.5 -ej 4.155 3 5 -ej 3.6 8 1 -ej 4.205 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt

#case 43
net=7_tax_sp_nt1_para00_bl06
echo case43_${net} > current_case
rm ms* hybridLambda*
#ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.0 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -ej 2 3 1 -ej 2.3 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 44
net=7_tax_sp_nt1_para02_bl06
echo case44_${net} > current_case
rm ms* hybridLambda*
ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.2 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 45
net=7_tax_sp_nt1_para05_bl06
echo case45_${net} > current_case
rm ms* hybridLambda*
ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.5 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 46
net=7_tax_sp_nt1_para08_bl06
echo case46_${net} > current_case
rm ms* hybridLambda*
ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.8 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 47
net=7_tax_sp_nt1_para10_bl06
echo case47_${net} > current_case
rm ms* hybridLambda*
ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 1.0 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 48
net=7_tax_sp_nt1_para_bl0
echo case48_${net} > current_case
rm ms* hybridLambda*
ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0 6 7 -ej 0 7 2 -ej 0 2 3 -ej 0 4 3 -es 0 3 0.8 -ej 0 3 5 -ej 0 8 1 -ej 0 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo

#case 49
net=7_tax_sp_nt1_para_bl100
echo case49_${net} > current_case
rm ms* hybridLambda*
ms 7 ${rep} -T -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 50.5 7 2 -ej 100.5 2 3 -ej 150.5 4 3 -es 200.5 3 0.5 -ej 250.5 3 5 -ej 250.5 8 1 -ej 300.5 1 5 | tail -n +4 | grep -v // | grep -v '^$' > ms_gt
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl

foo
