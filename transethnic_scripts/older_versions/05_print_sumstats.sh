for filename in theor_sumstats/*.out; do
sed '23q;d' $filename >> theor_sumstats/platelet_mcv_ldpred_sumstats.txt
done

for filename in theor_ld/*.out; do
sed '23q;d' $filename >> theor_ld/ld.txt
done


for filename in theor_sumstats_ct/*.out; do
grep -E "\[1\]" $filename >> theor_sumstats_ct/sumstats_ct.txt
done



#LDpred - WPC- Platelet - 21, 31
#LDpred - WPC - Platelet- 18