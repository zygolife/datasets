#!env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my $debug = 0;
my $force = 0;
my $exe = '/opt/GAL/0.2.2/bin/gal_protein_sequence';
my $cdsexe = '/opt/GAL/0.2.2/bin/gal_CDS_sequence';
my $dir = "genomes";
my $gff_file;
my $dnadir = File::Spec->catfile($dir,"DNA");
my $pepdir = File::Spec->catfile($dir,"pep");
my $cdsdir = File::Spec->catfile($dir,"CDS");

GetOptions('force!'   => \$force,
	   'exe:s'    => \$exe,
	   'gff:s'    => \$gff_file,
	   'dna:s'    => \$dnadir,
	   'pep:s'    => \$pepdir,
	   'cds:s'    => \$cdsdir,
	   'v|debug!' => \$debug,
    );
if( ! $gff_file || ! $dnadir ||
    ! -f $gff_file || ! -d $dnadir ) {
    warn("must provide GFF file and DNA dir with --gff and --dna\n");
}
my (undef,undef,$file) = File::Spec->splitpath($gff_file);

if( $file =~ /(\S+)\.gff(3)?/ ) {
    my $base = $1;
    my ($genus,$species) = split(/\_/,$base);
    my $pref = substr($genus,0,1).substr($species,0,3);
    if ( ! -f "$pepdir/$base.aa.fasta" || $force ) {
	warn("attempting PEP for $base\n");
	`$exe $gff_file $dnadir/$base.fasta | perl -p -e 's/^>/>$pref|/' > $pepdir/$base.aa.fasta`;
    }
    if( ! -f "$cdsdir/$base.CDS.fasta" || $force ) {
	warn("attempting CDS for $base\n");
	`$cdsexe $gff_file $dnadir/$base.fasta | perl -p -e 's/^>/>$pref|/' > $cdsdir/$base.CDS.fasta`;
    }
} else {
    die "must provide GFF file with .gff[3] extension\n";
}
