#!/usr/bin/perl
use strict;
use warnings;
my $debug =0;
my $debug_one = 0;
my $force = 0;
my $dir = shift || "download";
my $p_odir = shift || 'genomes/pep';
my $g_odir = shift || 'genomes/GFF';
my $d_odir = shift || 'genomes/DNA';
my $c_odir = shift || 'genomes/CDS';

for $dir ( $p_odir, $g_odir, $d_odir, $c_odir) {
 mkdir $dir;
}
opendir(DIR, $dir) || die $!;
for my $f ( readdir(DIR) ) {
    next if $f =~ /^\./;
    opendir(FAM,"$dir/$f") || die $!;
    for my $sp ( readdir(FAM) ) {
	next if $sp =~ /^\./;
	next if ( ! -d "$dir/$f/$sp");
	warn("sp is $dir/$f/$sp\n") if $debug;
	opendir(SP,"$dir/$f/$sp") || die "$dir/$f/$sp $!";
	for my $file ( readdir(SP) ) {
	    if( $file =~ /(\S+\.aa\.fasta).gz$/) {
		my $stem = $1;
		next if -f "$p_odir/$stem";
		print("zcat $dir/$f/$sp/$file | ",
		      'perl -p -e \'s/^>(\w+)\|(\S+)\|(\d+)\|(\S+)/>$2|$2_$3 $4/\' > ',"$p_odir/$stem\n");
            } elsif ( $file =~ /(\S+\.CDS\.fasta).gz$/) {
		my $stem = $1;
                next if -f "$c_odir/$stem";
                print("zcat $dir/$f/$sp/$file | ",
                      'perl -p -e \'s/^>(\w+)\|(\S+)\|(\d+)\|(\S+)/>$2|$2_$3 $4/\' > ',"$c_odir/$stem\n");
	    } elsif ( $file =~ /(\S+\.assembly\.fasta).gz$/) {
		my $stem = $1;
		$stem =~ s/\.assembly//;
		next if -f "$d_odir/$stem";
		print("zcat $dir/$f/$sp/$file > $d_odir/$stem\n");
	    } elsif ( $file =~ /(\S+\.gff3).gz$/) {
		my $stem = $1;
		next if !$force && -f "$g_odir/$stem";
		open(my $in => "zcat $dir/$f/$sp/$file |") || die "cannot open $dir/$f/$sp/$file: $!";
		open(my $out => ">$g_odir/$stem")|| die "cannot open $g_odir/$stem: $!";
		while(<$in>) {
		    last if /##FASTA/;
		    if( ! /^\#/ ) {
			chomp;
			my @row = split(/\t/,$_);
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
