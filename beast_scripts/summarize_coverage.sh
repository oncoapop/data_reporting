echo "Please enter the full path of the directory were we want to analyze coverage:"
   read indir
echo "Please enter the full path of the directory were we want to store the output:"
   read outdir
cd $indir
awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5+$6 } END { for(i=1;i<=FNR;i++) print a[i] }' $(ls -1v variant_status/*.tsv) > $outdir"/coverage.tsv"
