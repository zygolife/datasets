#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $debug =0;
my $debug_one = 0;
my $force = 0;

GetOptions('f|force!' => \$force,
	   'v|debug!' => \$debug,
	   'one!'     => \$debug_one,
    );
my $dir = shift || "download";
my $p_odir = shift || 'genomes/pep';
my $g_odir = shift || 'genomes/GFF';
my $d_odir = shift || 'genomes/DNA';
my $c_odir = shift || 'genomes/CDS';

for $dir ( $p_odir, $g_odir, $d_odir, $c_odir) {
 mkdir $dir;
}
opendir(DIR, $dir) || die $!;
open(my $lookup => ">lookup_prefix.dat") || die $!;
my %prefix;
for my $f ( readdir(DIR) ) {
    next if $f =~ /^\./;
    opendir(FAM,"$dir/$f") || die $!;
    for my $sp ( readdir(FAM) ) {	
	next if $sp =~ /^\./;
	next if ( ! -d "$dir/$f/$sp");
	warn("sp is $dir/$f/$sp\n") if $debug;
	opendir(SP,"$dir/$f/$sp") || die "$dir/$f/$sp $!";
	$sp =~ s/var\./var/;
	my ($stem1) = split(/\./,$sp);
#	my @n = split(/_/,$stem);
#	$stem =~ s/var\.?_\w+//;
	my ($g,$species,$rest) = split(/_/,$stem1,3);
	my $prefix = substr($g,0,1) . substr($species,0,3);
	if( $rest ) {
	    $rest =~ s/-/_/;
	    $rest =~ s/_#//;
	    my @n = split(/_/,$rest);
	    if( $n[0] eq 'var' ) {
		shift @n;
		shift @n;
	    }
	    $rest = join("_",@n);
	    $prefix .= "_$rest";
	}
	warn("g is $g species is $species -- $stem1 -- $prefix\n");
	$prefix{$stem1} = $prefix;
	for my $file ( readdir(SP) ) {
	    if( $file =~ /(\S+\.aa\.fasta).gz$/) {
		my $stem = $1;
		next if ! $force && -f "$p_odir/$stem";
		print("zcat $dir/$f/$sp/$file | ",
		      'perl -p -e \'s/^>(\w+)\|(\S+)\|(\d+)\|(\S+)/>$2|$2_$3 $4/\' > ',"$p_odir/$stem\n");
            } elsif ( $file =~ /(\S+\.CDS\.fasta).gz$/) {
		my $stem = $1;
                next if !$force && -f "$c_odir/$stem";
                print("zcat $dir/$f/$sp/$file | ",
                      'perl -p -e \'s/^>(\w+)\|(\S+)\|(\d+)\|(\S+)/>$2|$2_$3 $4/\' > ',"$c_odir/$stem\n");
	    } elsif ( $file =~ /(\S+\.assembly\.fasta).gz$/) {
		my $stem = $1;
		my (@all) = split(/\./,$stem);
		my $jgiprefix = $all[-4];
		$prefix{$stem1} = $jgiprefix;
		$stem =~ s/\.assembly//;
		next if !$force && -f "$d_odir/$stem";
		print("zcat $dir/$f/$sp/$file | perl -p -e 's/^>/>$jgiprefix|/' > $d_odir/$stem\n");
	    } elsif ( $file =~ /(\S+\.gff3).gz$/) {
		my $stem = $1;
		my (@all) = split(/\./,$stem);
		my $jgiprefix = $all[-3];
		$prefix{$stem1} = $jgiprefix;

		next if !$force && -f "$g_odir/$stem";
		
		open(my $in => "zcat $dir/$f/$sp/$file |") || die "cannot open $dir/$f/$sp/$file: $!";
		open(my $out => ">$g_odir/$stem")|| die "cannot open $g_odir/$stem: $!";
		while(<$in>) {
		    last if /##FASTA/;
		    if( ! /^\#/ ) {
			chomp;
			my @row = split(/\t/,$_);
			$row[0] = join("_",$prefix{$stem1}, $row[0]);
			my $last = pop @row;
			my (@order,%ninth);
			for my $ent ( split(/;/,$last) ) {
			    my ($id,$val) = split(/=/,$ent);
			    $ninth{$id} = $val;
			    push @order, $id;
			}
			if( exists $ninth{'Name'} ) {
			    my $val = $ninth{'Name'};
			    if( $val =~ /jgi\.p\|(\S+)\|(\d+)\|?$/ ) {
				$val = "$1|$1\_$2";
			    }
			    $ninth{'Name'} = $val;
			}
			push @row, join(";", map { sprintf("%s=%s", $_,
							   $ninth{$_}) }
					@order);
						       
			$_= join("\t",@row)."\n";
		    } 
		    print $out $_;
		}
		close($in);
		close($out);
	    }
	}
#	last if $debug_one;
    }
}

for my $stem ( sort keys %prefix ){
    print $lookup join("\t", $prefix{$stem}, $stem),"\n";
}
