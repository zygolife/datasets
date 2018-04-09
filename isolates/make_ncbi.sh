#PBS -l mem=256gb -j oe
module load perl
perl build_ncbi_zygo_taxa.pl
