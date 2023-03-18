#
# This script to run the full analysis if you specify an study and SNP:
#
# Please replace the two paths below with the locations
#  where you store the code and the data:
path_to_code=../Causal_genes_GWAS_loci_CAD
path_to_data=../data/GTEx_data

#  Choose study data file:
sample_study="Liv_3116"
chr=6
position=161108536
cp_suffix="_"$chr"_"$position

echo "# Begin Analysis:"
# Choose SNPs in the neighbourhood of the given SNP:
echo "# chromsome, position: " $chr", "$position
echo "# Choose SNPs in the neighbourhood:"
echo "#" python3 $path_to_code/Choose_SNPs.py  $path_to_data/"exp_"$sample_study".csv" $path_to_data/"out_"$sample_study".csv" $chr $position
python3 $path_to_code/Choose_SNPs.py  $path_to_data/"exp_"$sample_study".csv" $path_to_data/"out_"$sample_study".csv" $chr $position

# Data preparation:
echo "# Data preparation:"
echo "#" python3 $path_to_code/Data_preperation.py  $path_to_data/"exp_"$sample_study$cp_suffix".csv" $path_to_data/"out_"$sample_study$cp_suffix".csv"
python3 $path_to_code/Data_preperation.py  $path_to_data/"exp_"$sample_study$cp_suffix".csv" $path_to_data/"out_"$sample_study$cp_suffix".csv"

# Run MVMR_withoutLD:
echo "# Run MVMR_withoutLD:"
prep_suffix=$cp_suffix"_prepared"
echo "#" python3 $path_to_code/MVMR_withoutLD.py  $path_to_data/"exp_"$sample_study$prep_suffix".csv"
python3 $path_to_code/MVMR_withoutLD.py  $path_to_data/"exp_"$sample_study$prep_suffix".csv"

echo "# Done."
