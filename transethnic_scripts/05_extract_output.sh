for filename in R-ss_ldpred_*.out; do
sed '23q;d' $filename >> data/theoretical/ldpred_sumstats.txt
done

for filename in R-ss_ss_ct_*.out; do
sed '23q;d' $filename >> data/theoretical/ct_sumstats.txt
done


for filename in R-ld_*.out; do
grep -hnr "fst LD|wpc LD" $filename >> data/theoretical/ld.txt
done



