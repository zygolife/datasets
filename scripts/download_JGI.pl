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
#my $genome_file = 'lib/organisms.csv';
# data fields
my @data_fields = qw(label filename md5 timestamp); # removed 'url' as it was redundant here
my $download_cmds = "jgi_download_curl.sh";
my $cookie_file = '.JGI_cookies';
my $debug = 0;
my $infile;  # input file with XML
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

my %header;
my $row = $csv->getline ($fh);

my $i =0;
for my $c ( @$row ) {
    $header{$c} = $i++;
}

my $i =0;
for my $c ( @$row ) {
    $header{$c} = $i++;
}

while (my $row = $csv->getline ($fh)) {
    if ( ! defined $row->[$header{Source}] ) {
	warn("no source for: ",join(",",@$row),"\n");
    }
    my  $name = join(" ", $row->[$header{species}]);
    if( $row->[$header{Source}] eq 'JGI' ) {
	# keep track of the species which we would like to get from JGI
	$jgi_targets{$name} = 1;
    } else {
	next;
    }

    #$row->[$header{strain}]);
    if( exists $orgs{$name} ) {
	warn("writing over storage for $name \n");
    }
    $orgs{$name} = { 'strain' => $row->[$header{strain}],
		     'family' => $row->[$header{Family}],
		     'source' => $row->[$header{Source}],
		     'species' => $row->[$header{species}]};
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
		    my $filename = $file->{filename};
		    warn("prefix is $prefix and filename is $filename\n") if $debug;

		    my ($name,$version) = &standardize_name($file->{label});
	 	    
		    if( ! exists $orgs{$name} ) {
			warn("cannot find '$name' in the query file\n") if $debug;
			next;
		    }
		    warn("filename is $filename\n") if $debug; 
		    next if $filename =~ /MitoScaffold|Mitochondria|mito\.scaffolds/;
		    if( $jgi_targets{$name} ) {
			my $family = $orgs{$name}->{family};
			my $oname = $name;
			$oname =~ s/\s+/_/g;
			$oname =~ s/\.//g;
			if( my $str = $orgs{$name}->{strain} ) {
			    $str =~ s/\s+/_/g;
			    $oname .= "_" . $str;
			}
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
			$sanity{$outfile}->{ $filename }++;
			if( ! -f $outfile ) {
			    print $curl_cmds "curl $fullurl -b $cookie_file -o $outfile --create-dirs\n";
			}
			$jgi_targets{$name} = 2;
		    }
		}
	    } 
	}
    } elsif( $type eq 'Annotation' ) {
	my $asm = $d->{folder};
	while( my ($k,$n) = each %$asm ) {
	    warn("keys are $k\n");
	    if( $k eq 'Filtered Models ("best")' ) {
		for my $ftype ( qw(Proteins CDS Transcripts Genes) ) {
		    my $oftype = $ftype;
		    my $f = $n->{'folder'}->{$ftype};
		    for my $file_obj ( 
			sort { $a->[0] cmp $b->[0] }
			map {[(split(/\//,$_->{url}))[1],$_]  } 
			@{$f->{file}} ) {
			my ($prefix,$file) = @$file_obj;			
			my $filename = $file->{filename};
			
# skip the compiled set of all gene modules (tar.gz) or the EST fastas
			next if( $filename =~ /\.tar.gz$/ || 
				 $filename =~ /ESTs|unsupported_short/ ||
				 $filename =~ /\.nt\.fasta/ );
			
			my $fullurl = sprintf("%s%s",$url_base,
					      $file->{url});
			my ($name,$version) = &standardize_name($file->{label});
                        warn("Unfiltered filename is $filename\n");
			warn("name is $name ftype is $ftype url is $fullurl\n") if $debug;
			if( ! exists $orgs{$name} ) {
			    warn("cannot find '$name' in the query organisms file\n") if $debug;
			    next;
			}
			
			if( $jgi_targets{$name} ) {
			    my $family = $orgs{$name}->{family};
			    my $oname = $name;
			    $oname =~ s/\s+/_/g;
			    $oname =~ s/\.//g;
			    if( my $str = $orgs{$name}->{strain} ) {
				$str =~ s/\s+/_/g;
				$oname .= "_" . $str;
			    }
			    			    
			    my $outfile = File::Spec->catdir($outdir,$family,$oname);

			    my $oname_labeled = "$oname.$prefix";
			    if( $version =~ s/(v\d+)(\.0)?/$1/ ) {
				$oname_labeled .= ".$version";
			    }
			    warn("$filename for $oname_labeled\n");

			    if ( $filename =~ /\.tab(\.gz)?$/ ) { 
				warn("skipping $filename\n") if $debug;
				next;
			    } elsif( $filename =~ /\.gff/ ) {
				$outfile = File::Spec->catfile($outfile,"$oname_labeled.gff3.gz");
				warn("GFF outfile $outfile for $filename\n");
				if( $ftype ne 'Genes') {
				    warn("$filename misclassified file as $ftype not 'Genes'\n");
				    $oftype = 'Genes';
				}
			    } elsif( $filename =~ /\.aa\./ || 
				     $filename =~ /FilteredModels\d*\.(aa|proteins)|best_proteins|filtered_proteins|GeneCatalog\_?\d+\.proteins|_proteins/) 
			    {
				warn("AA outfile $outfile for $filename\n");
				$outfile = File::Spec->catfile($outfile,"$oname_labeled.aa.fasta.gz");
				$oftype = 'Proteins';
				if ( exists $sanity{$outfile} ) { 
				    warn("(AA) $filename second time around for $outfile (previous ",join(",",keys %{$sanity{$outfile}}),")\n");
				    next;
				}
			    } elsif( $filename =~ /[_\.](CDS|transcripts)/i ||
				     $filename =~ /best_transcripts|filtered_transcripts|_transcripts|FilteredModels\d*\.na|geneCatalog_CDS_\d+/){ 
				$outfile = File::Spec->catfile($outfile,"$oname_labeled.CDS.fasta.gz");
				$oftype = 'Transcripts';
				if (exists $sanity{$outfile} ) {
				    warn("(CDS) $filename second time around for $outfile (previous ",join(",",keys %{$sanity{$outfile}}),")\n");
				    # deal with mis-placed or missing CDS sections that are listed as transcript sections
				    next;
				}
			    } else {
				warn("WARNING: unmatched filename is $filename\n");
				next;
			    }
			    if( ! $sanity{$outfile}->{$filename}++ ) {
				print join("\t",  $oftype,$prefix,$outfile,
					   map { $file->{$_} || ''} 
					   @data_fields),"\n";
			    }
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
    if( scalar keys %{$sanity{$outfile}} > 1) {
	warn("$outfile seen more than once, will be overwritten ",join(" ",keys %{$sanity{$outfile}}),"\n");
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

# how to fix the names for each species with some inconsistencies
sub standardize_name {
    my $label = shift;
    $label =~ s/f\.sp\./f. sp./;
    $label =~ s/(\s+var)(\s+)/$1.$2/;
    $label =~ s/\s+$//;
    # special case
    $label =~ s/ \([^\)]+\)//;
    
    my @label_spl = split(/\s+/,$label);      
    my $version = $label_spl[-1];
    my $name;
    if( $label =~ /\s+var(\.)?\s+/ ) {
	$name = join(" ", $label_spl[0],$label_spl[1],
		     "var.",$label_spl[3]);
    } elsif ($label =~ /\s+f\.\s*sp\.\s+/ ) {
	$name = join(" ", $label_spl[0],$label_spl[1],"f. sp.",
		     $label_spl[4]);
	warn ">>>>$name\n";
	#$label_spl[4]); # we don't put strain in there
	#warn("DNA fsp -> $name\n") if $debug; 
    } else {
	$name = join(" ", $label_spl[0],$label_spl[1]);
    }
    ($name,$version);
}
