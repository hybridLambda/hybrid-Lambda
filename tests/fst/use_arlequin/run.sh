#!/bin/bash
arlequin=arlecore3513_32bit #location of the binary executable

rm pairwise_fst_summary
hybrid-Lambda -spcu "(A:100, B:100);" -seg -fst -S 30 30 -seed 101

rm -r current_arlequin.res
${arlequin} current_arlequin.arp fst_setting.ars > /dev/null
sed -e '1,/== Comparisons of pairs of population samples/d' current_arlequin.res/current_arlequin.htm | sed -e '/END OF RUN/,$d' > site1_arlequin.summary
head -18 site1_arlequin.summary | tail -1 >> pairwise_fst_summary

cat pairwise_fst_summary
cat OUT_fst

R CMD BATCH hybridLambda_fst.r
tail -10 hybridLambda_fst.r.Rout
