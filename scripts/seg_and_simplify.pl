#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my $ext = '.aa.fasta';
my $debug = 0;

GetOptions(
    'v|debug!'   => \$debug,
    'ext:s'      => \$ext,
    );
my $lookup = shift || 'lookup_prefix.dat';
my $dir = shift || ".";
my $odir = shift || ".";

open(L, $lookup) || die "cannot open $lookup: $!";
my %lk;
while(<L>) {
    my ($prefix,$stem) = split;
    warn("$stem => $prefix\n") if $debug;
    $lk{$stem} = $prefix;
}
opendir(DIR, $dir) || die "cannot open dir $dir: $!";
for my $file ( readdir(DIR) ) {
    next if $file =~ /^\./;
    next unless ( $file =~ /(\S+)\Q$ext\E/ );
    my $stem = $1;
    $stem =~ s/var\./var/;
    my ($base) = split(/\./,$stem);
    if( ! $lk{$base} ) {
	warn("cannot find $base in lookup file ($stem)\n");
	next;
    }
    print "pseg $dir/$file -q -z 1 > $odir/$lk{$base}.seg\n";
}
