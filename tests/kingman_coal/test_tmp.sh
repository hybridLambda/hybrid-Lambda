#!/bin/bash

dir=test-demo
mkdir ${dir}
cd ${dir}
rm *pdf


rep=1000

## compare population sturture for a single population data
COMPAREFILE=compareDemo
rm ${COMPAREFILE}

#source chisq_r.src

source ../ks_r.src

source ../tmrca_r.src

source ../bl_r.src

source ../process_sample_stats.src	


#case 7
net=5_tax_sp_nt1_para_bl0
CURRENTCASE=case7_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0 4 3 -ej 0 2 3 -es 0 3 0.5 -ej 0 6 5 -es 0 3 0.5 -ej 0 7 5 -ej 0 3 1 -ej 0 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

##ms 5 10000 -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej 1.05 2 3 -es 2.25 3 0.5 -ej 3.25 6 5 -es 3.3 3 0.5 -ej 4.8 7 5 -ej 4.95 3 1 -ej 5.45 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats

###########################################################################################################
#5taxnt2

#case 8
net=5_tax_sp_nt2_para00_bl06
CURRENTCASE=case8_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
#ms 5 10000 -T -t 10 -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.0 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -ej 0.8 3 1 -ej 1.1 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 9
net=5_tax_sp_nt2_para02_bl06
CURRENTCASE=case9_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.2 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case  10
net=5_tax_sp_nt2_para05_bl06
CURRENTCASE=case10_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.5 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 11
net=5_tax_sp_nt2_para08_bl06
CURRENTCASE=case11_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 0.8 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 12
net=5_tax_sp_nt2_para10_bl06
CURRENTCASE=case12_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es 0.5 3 1.0 -ej 0.8 3 5 -ej 0.8 6 1 -ej 1.1 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 13
net=5_tax_sp_nt2_para_bl0
CURRENTCASE=case13_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0 4 3 -ej 0 2 3 -es 0 3 0.5 -ej 0 3 5 -ej 0 6 1 -ej 0 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 14
net=5_tax_sp_nt2_para_bl100
CURRENTCASE=case14_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0.5 4 3 -ej 0.5 2 3 -es .5 3 0.5 -ej 50.5 3 5 -ej 50.5 6 1 -ej 100.5 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo


#case 15
net=5_tax_sp_nt3_para00_bl06
CURRENTCASE=case15_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
#ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0 -ej 1.1 3 5 -es 1.1 6 0 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0 -ej 1.1 3 5 -es 1.1 6 0 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 16
net=5_tax_sp_nt3_para02_bl06
CURRENTCASE=case16_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*

ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0.2 -ej 1.1 3 5 -es 1.1 6 0.2 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 17
net=5_tax_sp_nt3_para05_bl06
CURRENTCASE=case17_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0.5 -ej 1.1 3 5 -es 1.1 6 0.5 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 18
net=5_tax_sp_nt3_para08_bl06
CURRENTCASE=case18_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 0.8 3 0.8 -ej 1.1 3 5 -es 1.1 6 0.8 -ej 1.4 6 5 -ej 1.4 7 1 -ej 1.7 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 19
net=5_tax_sp_nt3_para10_bl06
CURRENTCASE=case19_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -ej 1.1 3 5 -ej 1.7 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 20
net=5_tax_sp_nt3_para_bl0
CURRENTCASE=case20_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej 0 4 3 -ej 0 2 3 -es 0 3 0.5 -ej 0 3 5 -es 0 6 0.5 -ej 0 6 5 -ej 0 7 1 -ej 0 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo

#case 21
net=5_tax_sp_nt3_para_bl100
CURRENTCASE=case21_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 5 ${rep} -T -t 10 -I 5 1 1 1 1 1 -ej .5 4 3 -ej .5 2 3 -es 50.5 3 0.5 -ej 100.5 3 5 -es 100.5 6 0.5 -ej 150.5 6 5 -ej 150.5 7 1 -ej 200.5 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 5

foo


###########################################################################################################
#6taxnt1
#ms 6 100000 -T -t 10 -I 6 1 1 1 1 1 1 -ej 2.8 6 2 -ej 2.85 2 3 -ej 3.35 4 3 -es 4.35 3 0.5 -ej 5.405 3 5 -ej 5.355 7 1 -ej 5.455 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats

#case 22
net=6_tax_sp_nt1_para00_bl06
CURRENTCASE=case22_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
#ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.5 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -ej 1.7 3 1 -ej 2 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 23
net=6_tax_sp_nt1_para02_bl06
CURRENTCASE=case23_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.2 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 24
net=6_tax_sp_nt1_para05_bl06
CURRENTCASE=case24_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.5 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 25
net=6_tax_sp_nt1_para08_bl06
CURRENTCASE=case25_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 0.8 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 26
net=6_tax_sp_nt1_para10_bl06
CURRENTCASE=case26_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 0.8 2 3 -ej 1.1 4 3 -es 1.4 3 1 -ej 1.7 3 5 -ej 1.7 7 1 -ej 2 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 27
net=6_tax_sp_nt1_para_bl0
CURRENTCASE=case27_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0 6 2 -ej 0 2 3 -ej 0 4 3 -es 0 3 .5 -ej 0 3 5 -ej 0 7 1 -ej 0 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 28
net=6_tax_sp_nt1_para_bl100
CURRENTCASE=case28_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 6 2 -ej 50.5 2 3 -ej 100.5 4 3 -es 150.5 3 .5 -ej 200.5 3 5 -ej 200.5 7 1 -ej 250.5 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo


###########################################################################################################
#6taxnt2
#ms 6 100000 -T -t 10 -I 6 1 1 1 1 1 1 -ej 1.85 2 3 -ej 2.35 4 3 -es 2.35 3 0.5 -ej 1.15 6 5 -ej 2.855 5 3 -ej 3.605 7 1 -ej 3.705 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats

#case 29
net=6_tax_sp_nt2_para00_bl106
CURRENTCASE=case29_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
#ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.0 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -ej 0.8 6 5 -ej 1.1 3 1 -ej 1.4 5 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 30
net=6_tax_sp_nt2_para02_bl106
CURRENTCASE=case30_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.2 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 31
net=6_tax_sp_nt2_para05_bl106
CURRENTCASE=case31_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.5 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 32
net=6_tax_sp_nt2_para08_bl106
CURRENTCASE=case32_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 0.8 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 33
net=6_tax_sp_nt2_para10_bl106
CURRENTCASE=case33_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.8 4 3 -es 0.8 3 1 -ej 0.8 6 5 -ej 1.1 5 3 -ej 1.1 7 1 -ej 1.4 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 34
net=6_tax_sp_nt2_para_bl0
CURRENTCASE=case34_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0 2 3 -ej 0 4 3 -es 0 3 .5 -ej 0 6 5 -ej 0 5 3 -ej 0 7 1 -ej 0 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 35
net=6_tax_sp_nt2_para_bl100
CURRENTCASE=case35_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 50.5 4 3 -es 50.5 3 .5 -ej 50.5 6 5 -ej 100.5 5 3 -ej 100.5 7 1 -ej 150.5 3 1 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

###########################################################################################################
#6taxnt3
#ms 6 100000 -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.85 2 3 -ej 0.85 6 3 -ej 0.85 4 3 -es 1.85 3 0.34 -ej 2.605 3 5 -ej 2.855 7 1 -ej 3.155 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats

#case 36
net=6_tax_sp_nt3_para00_bl06
CURRENTCASE=case36_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
#ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.0 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -ej 1.1 3 1 -ej 1.4 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 37
net=6_tax_sp_nt3_para02_bl06
CURRENTCASE=case37_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.2 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 38
net=6_tax_sp_nt3_para05_bl06
CURRENTCASE=case38_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.5 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 39
net=6_tax_sp_nt3_para08_bl06
CURRENTCASE=case39_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 0.8 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 40
net=6_tax_sp_nt3_para10_bl06
CURRENTCASE=case40_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 0.8 3 1 -ej 1.1 3 5 -ej 1.1 7 1 -ej 1.4 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 41
net=6_tax_sp_nt3_para_bl0
CURRENTCASE=case41_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0 2 3 -ej 0 6 3 -ej 0 4 3 -es 0 3 .5 -ej 0 3 5 -ej 0 7 1 -ej 0 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

#case 42
net=6_tax_sp_nt3_para_bl100
CURRENTCASE=case42_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 6 ${rep} -T -t 10 -I 6 1 1 1 1 1 1 -ej 0.5 2 3 -ej 0.5 6 3 -ej 0.5 4 3 -es 50.5 3 .5 -ej 100.5 3 5 -ej 100.5 7 1 -ej 150.5 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 6

foo

##########################################################################################################

#7taxnt1
#ms 7 100000 -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.05 6 7 -ej 0.25 7 2 -ej 0.8 2 3 -ej 2.45 4 3 -es 3.45 3 0.5 -ej 4.155 3 5 -ej 3.6 8 1 -ej 4.205 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats

#case 43
net=7_tax_sp_nt1_para00_bl06
CURRENTCASE=case43_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
#ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.0 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -ej 2 3 1 -ej 2.3 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 7

foo

#case 44
net=7_tax_sp_nt1_para02_bl06
CURRENTCASE=case44_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.2 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 7

foo

#case 45
net=7_tax_sp_nt1_para05_bl06
CURRENTCASE=case45_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.5 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 7

foo

#case 46
net=7_tax_sp_nt1_para08_bl06
CURRENTCASE=case46_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 0.8 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 7

foo

#case 47
net=7_tax_sp_nt1_para10_bl06
CURRENTCASE=case47_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 0.8 7 2 -ej 1.1 2 3 -ej 1.4 4 3 -es 1.7 3 1.0 -ej 2 3 5 -ej 2 8 1 -ej 2.3 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 7

foo

#case 48
net=7_tax_sp_nt1_para_bl0
CURRENTCASE=case48_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0 6 7 -ej 0 7 2 -ej 0 2 3 -ej 0 4 3 -es 0 3 0.8 -ej 0 3 5 -ej 0 8 1 -ej 0 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 7

foo

#case 49
net=7_tax_sp_nt1_para_bl100
CURRENTCASE=case49_${net}
echo -e ${CURRENTCASE} > current_case
rm -rf ms* hybridLambda*
ms 7 ${rep} -T -t 10 -I 7 1 1 1 1 1 1 1 -ej 0.5 6 7 -ej 50.5 7 2 -ej 100.5 2 3 -ej 150.5 4 3 -es 200.5 3 0.5 -ej 250.5 3 5 -ej 250.5 8 1 -ej 300.5 1 5 -t 10 > msout
grep ";" msout > ms_gt
cat msout | sample_stats > ms_stats
hybrid-Lambda -gt ms_gt -o ms -tmrca -bl
hybrid-Lambda -spcu ../../../trees/${net} -num ${rep} -o hybridLambda -tmrca -bl -mu 0.0005 -seg -sim_Si_num

rearrange_hybridLambdaout 7

foo
