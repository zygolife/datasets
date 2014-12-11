#!env perl
use strict;
use warnings;
use Text::CSV_XS qw(csv);
use XML::Simple;
use Data::Dumper;
use Getopt::Long;
use File::Spec;


my $url_base = 'http://genome.jgi.doe.gov';
my $outdir = 'download';

my $genome_file = 'genomes/zygo_genomes.csv';
# data fields
my @data_fields = qw(label filename md5 timestamp); # removed 'url' as it was redundant here
my $download_cmds = "jgi_download_curl.sh";
my $cookie_file = '.JGI_cookies';
my $debug = 0;
my $infile = 'fungi.xml';  # input file with XML
GetOptions(
    'v|debug|verbose!' => \$debug,
    'i|input:s'        => \$infile,
    'cookie:s'         => \$cookie_file,
    'g|genomes:s'      => \$genome_file,
    'o|dir|outdir:s'   => \$outdir,
    'download|script:s' => \$download_cmds,
    );

$infile = shift @ARGV unless defined $infile;

die "must provide an input file via -i or single argument" 
    if ! defined $infile;

#init the XML parser
my $xs = XML::Simple->new;

print join("\t", qw(Type Prefix LocalFile), @data_fields),"\n";
# Read whole file in memory as hash of hashes
my $csv = Text::CSV_XS->new ({ binary => 1, auto_diag => 1 });
open my $fh, "<:encoding(utf8)", $genome_file or die "$genome_file: $!";

my $parsed = $xs->XMLin($infile);
my $folder = $parsed->{folder}->{folder};

my %orgs;
my %jgi_targets;
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
					 'source' => $row->[ $header{Source} ] };
    if( $row->[$header{Source} ] eq 'JGI' ) {
	# keep track of the species which we would like to get from JGI 
	$jgi_targets{$row->[ $header{species} ]} = 1;
    }
}

open(my $curl_cmds => ">$download_cmds") || die $!;

my %sanity;
#warn(Dumper($fungi));
warn("types are: ", join(",", keys %$folder),"\n") if $debug;
my %data;
my $first = 1;
while( my ($type,$d) = each %$folder ) {
    if( $type eq 'Assembly' ) {
	my $asm = $d->{folder};
	while( my ($k,$n) = each %$asm ) {
	    if( $k eq 'Assembled scaffolds (unmasked)') {
		# sorting here on the prefix by extracting the prefix from 
		# $file->{url} which has '/' separating the path
		# so (split(/\//,$_->{url}))[1] gets the 2nd entry which 
		# is the prefix
		# this creates a new data structure
		# [ prefix, filehash (from the XML)]
		# now we can sort by the prefix -> hence the sort
		for my $file_obj ( 
		    sort { $a->[0] cmp $b->[0] }
		    map {[(split(/\//,$_->{url}))[1],$_]  } 
			       @{$n->{file}} ) {
		    my ($prefix,$file) = @$file_obj;
		    my $fullurl = sprintf("%s%s",$url_base,
					   $file->{url});
		    $file->{label} =~ s/(\s+var)(\s+)/$1.$2/;
		    my $filename = $file->{filename};
		    my @label_spl = split(/\s+/,$file->{label});
		    my $version = $label_spl[-1];
		    my $name;
		    if( $file->{label} =~ /\s+var(\.)?\s+/ ) {
			$name = join(" ", $label_spl[0],$label_spl[1],
				     $label_spl[2],$label_spl[3]);
		    } elsif ($file->{label} =~ /\s+f\.\s+sp\.\s+/ ) {
			$name = join(" ", $label_spl[0],$label_spl[1],
				     $label_spl[2],$label_spl[3],
				     $label_spl[4]);
		    } else {
			$name = join(" ", $label_spl[0],$label_spl[1]);
		    }		    
		    if( ! exists $orgs{$name} ) {
			#warn("cannot find '$name' in the query file\n");
			next;
		    }
		    next if $filename =~ /MitoScaffold|Mitochondria|mito\.scaffolds/;
		    if( $jgi_targets{$name} ) {
			my $family = $orgs{$name}->{family};
                        my $strain = $orgs{$name}->{strain};
			my $oname = $name . " ". $strain;
			$oname =~ s/\s+/_/g;
			$oname =~ s/\.//g;
#			$oname =~ s/\.$//g; # remove trailing periods for things that are like XX sp. 				
			my $oname_labeled;
			if( $version =~ s/(v\d+)(\.0)?/$1/ ) {
			    $oname_labeled = "$oname.$prefix.$version.assembly.fasta.gz";
			} else {
			    $oname_labeled = "$oname.$prefix.assembly.fasta.gz";
			}
			my $outfile = File::Spec->catfile($outdir,$family,$oname,$oname_labeled);

			print join("\t", "Genome",$prefix,
				   $outfile,
				   map { $file->{$_} || ''} @data_fields),"\n";
			push @{$sanity{$outfile}}, $filename;
			if( ! -f $outfile ) {
			    print $curl_cmds "curl $fullurl -b $cookie_file -o $outfile --create-dirs\n";
			} else {
	warn("skipping $outfile as it already exists\n");
			}
			$jgi_targets{$name} = 2;
		    }
		}
	    } 
	}
    } elsif( $type eq 'Annotation' ) {
	my $asm = $d->{folder};
	while( my ($k,$n) = each %$asm ) {
	    #warn("keys are $k\n");
	    if( $k eq 'Filtered Models ("best")' ) {
		for my $ftype ( qw(Proteins CDS Transcripts) ) {
		    my $f = $n->{'folder'}->{$ftype};
		    for my $file_obj ( 
			sort { $a->[0] cmp $b->[0] }
			map {[(split(/\//,$_->{url}))[1],$_]  } 
			@{$f->{file}} ) {
			my ($prefix,$file) = @$file_obj;			
			my $fullurl = sprintf("%s%s",$url_base,
					      $file->{url});
			$file->{label} =~ s/(\s+var)(\s+)/$1.$2/;
			my $filename = $file->{filename};

			next if( $filename =~ /\.tar.gz$/ || 
				 $filename =~ /ESTs|unsupported_short/ ||
				 $filename =~ /\.nt\.fasta/ );
# skip the compiled set of all gene modules (tar.gz) or the EST fastas
			
			my @label_spl = split(/\s+/,$file->{label});
			my $version = $label_spl[-1];
			my $name;
			if( $file->{label} =~ /\s+var(\.)?\s+/ ) {
			    $name = join(" ", $label_spl[0],$label_spl[1],
					 $label_spl[2],$label_spl[3]);
			} elsif ($file->{label} =~ /\s+f\.\s+sp\.\s+/ ) {
			    $name = join(" ", $label_spl[0],$label_spl[1],
					 $label_spl[2],$label_spl[3],
					 $label_spl[4]);
			} else {
			    $name = join(" ", $label_spl[0],$label_spl[1]);
			}
			warn("name is $name ftype is $ftype url is $fullurl\n") if $debug;
			if( ! exists $orgs{$name} ) {
			    #warn("cannot find '$name' in the query file\n");
			    next;
			}
			
			if( $jgi_targets{$name} ) {
			    my $family = $orgs{$name}->{family};
			    my $strain = $orgs{$name}->{strain};
			    my $oname = "$name $strain";
			    $oname =~ s/\s+/_/g;
			    $oname =~ s/\.$//; # remove trailing periods for things that are like XX sp. 
			    my $outfile = File::Spec->catdir($outdir,$family,$oname);

			    my $oname_labeled = "$oname.$prefix";
			    if( $version =~ s/(v\d+)(\.0)?/$1/ ) {
				$oname_labeled .= ".$version";
			    }
			    #warn("$filename for $oname_labeled\n");

			    if( $filename =~ /\.gff/ ) {
				$outfile = File::Spec->catfile($outfile,"$oname_labeled.gff3.gz");
				#warn("GFF outfile is $outfile\n");
			    } elsif( $filename =~ /\.aa\./ || 
				     $filename =~ /FilteredModels\d*\.aa|best_proteins|filtered_proteins|GeneCatalog\_?\d+\.proteins|_proteins/) 
			    {
				warn"$filename";
				$outfile = File::Spec->catfile($outfile,"$oname_labeled.aa.fasta.gz");
				if ( exists $sanity{$outfile} ) { 
				    warn("(AA) $filename second time around for $outfile (previous @{$sanity{$outfile}})\n");
				    next;
				}
			    } elsif( $filename =~ /[_\.](CDS|transcripts)/i ||
				     $filename =~ /best_transcripts|filtered_transcripts|_transcripts|FilteredModels\d*\.na|geneCatalog_CDS_\d+/){ 
				$outfile = File::Spec->catfile($outfile,"$oname_labeled.CDS.fasta.gz");
				if (exists $sanity{$outfile} ) {
				    warn("(CDS) $filename second time around for $outfile (previous @{$sanity{$outfile}}\n");
				    # deal with mis-placed or missing CDS sections that are listed as transcript sections
				    next;
				}
			    } else {
				warn("WARNING: unmatched filename is $filename\n");
				next;
			    }
			    
			    print join("\t",  $ftype,$prefix,$outfile,
				       map { $file->{$_} || ''} 
				       @data_fields),"\n";
			
			    push @{$sanity{$outfile}}, $filename;
			    if( ! -f $outfile ) {
				print $curl_cmds "curl $fullurl -b $cookie_file -o $outfile --create-dirs\n";
			    }
			    $jgi_targets{$name} = 2;
			} else {
			    # warn("skipping $name not an intending JGI to be the target source\n");
			}
		    }		    
		}
	    }
	}
    }
}

for my $outfile ( keys %sanity ) {
    if( scalar @{$sanity{$outfile}} > 1) {
	warn("$outfile seen more than once, will be overwritten\n");
    }
}
my @missing;
for my $name ( sort keys %jgi_targets ) {
    if( $jgi_targets{$name} == 1 ) {
	push @missing, $name;
    }
}
if( @missing ) {
    warn("Missing the following JGI targets: \n", join("\n", @missing),"\n");
}
