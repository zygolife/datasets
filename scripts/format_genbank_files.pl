#!/usr/bin/perl 
use strict;
use warnings;
use Data::Dumper;
use File::Path;
use Getopt::Long;
use File::Spec;
use Text::CSV_XS qw(csv);
use Bio::SeqIO;
use IO::String;

my $force = 0;
my $debug = 0;

my $basedir = 'download';
my $targetdir = 'genomes';

GetOptions(
    'retmax:i'  => \$retmax,
    'f|force!'  => \$force, # force downloads even if file exists
    'v|d|debug|verbose!'   => \$debug,
    'f|force!'             => \$force,
    'b|basedir:s'          => \$basedir,
    't|targetdir:s'        => \$targetdir,);

my $ncbi_id_file = shift || 'genomes/zygo_genomes.csv';
if( ! -d $targetdir && ! -d $basedir ) {
    die "need a target and basedir\n";
}

my $csv = Text::CSV_XS->new ({ binary => 1, auto_diag => 1 });
open my $fh, "<:encoding(utf8)", $ncbi_id_file or die "$ncbi_id_file: $!";

my %orgs;
my %gbk_targets;
my $header = $csv->getline ($fh);
my %header;
my $x = 0;
for my $r ( @$header ) {
    $header{ $r } = $x++;
}

while (my $row = $csv->getline ($fh)) {
    next if( $row->[0] =~ /^\#/);
    $orgs{$row->[ $header{species}]} = { 'strain' => $row->[ $header{strain} ],
					 'family' => $row->[ $header{Family} ],
					 'source' => $row->[ $header{Source} ],
					 'accessions' => $row->[ $header{Accession}] };
    if( $row->[$header{Source} ] =~ /GB/ ) {
	# keep track of the species which we would like to get from GenBank
	$gbk_targets{$row->[ $header{species} ]} = 1;
    }
}


