#PBS -q js -j oe -l nodes=1:ppn=1,walltime=24:00:00 -N GAL_extract_pep
module load perl

N=$PBS_ARRAYID

#MAKE gbkfiles list with 
# make need additional curration
# ls genomes/GFF | grep -v v1 | grep -v v2 > genbank_genomes.dat
INFILE=genbank_genomes.dat
OUTDIR=genomes
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 echo "Need a PBS_ARRAYID or cmdline number"
 exit;
fi

if [ ! -f $INFILE ]; then
 ls $OUTDIR/GFF/*.gff3 | grep -v v1 | grep -v v2 > $INFILE
fi

line=`head -n $N $INFILE  | tail -n 1`
perl scripts/extract_pep_GAL_arrayjob.pl --gff $line --dna $OUTDIR/DNA
