cd /home/dyap/dyap_temp/ChIPseqAnalysis3/CX5461/EXPT-233

echo "Intersect all replicates..."

/home/dyap/bin/bedtools2/bin/multiIntersectBed  -header -i RAD51-IgG-10-7CX#1_c3.0_common.bed RAD51-IgG-10-7CX#2_c3.0_common.bed RAD51-IgG-10-7CX#3_c3.0_common.bed \
RAD51-IgG-nodrug_c3.0_common.bed > Expt.bed

echo "Intersect drug treated replicates..."

/home/dyap/bin/bedtools2/bin/multiIntersectBed  -header -i RAD51-IgG-10-7CX#1_c3.0_common.bed RAD51-IgG-10-7CX#2_c3.0_common.bed RAD51-IgG-10-7CX#3_c3.0_common.bed  \
 > RAD51-IgG-10-7CX_Reps.bed

echo "Intersect all replicate but only get the common regions..."

/home/dyap/bin/bedtools2/bin/multiIntersectBed  -header -i RAD51-IgG-10-7CX#1_c3.0_common.bed RAD51-IgG-10-7CX#2_c3.0_common.bed RAD51-IgG-10-7CX#3_c3.0_common.bed | \
awk -F"\t" '{ if ($4 == 3)  print $1,$2,$3 }' > RAD51-IgG-10-7CX_CommonRep.bed



