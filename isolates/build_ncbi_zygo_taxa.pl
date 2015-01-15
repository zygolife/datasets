#!env perl
use strict;
use Bio::DB::Taxonomy;

mkdir 'idx' unless -d 'idx';
my $db = Bio::DB::Taxonomy->new(
    -source => 'flatfile',
    -directory => 'idx',
    -nodesfile => 'tmp/nodes.dmp',
    -namesfile => 'tmp/names.dmp'
);

my @zygoclades = qw(Entomophthoromycota Kickxellomycotina Zoopagomycotina
Mucoromycotina Mortierellomycotina Glomeromycota
Nephridiophagidae
Olpidiaceae
);
for my $toplin ( @zygoclades ) {
    my @taxonids = $db->get_taxonids($toplin);
    if( @taxonids > 1 ) {
	warn("multiple taxa have this name $toplin, ambiguous");
	next;
    } elsif ( ! @taxonids ) {
	warn("No Taxa with name $toplin\n");
	next;
    } else {
	my $node = $db->get_taxon(-name => $toplin);
	my @genus = &get_genus($db,$node);
	print "$toplin\n";
	for my $g ( @genus ) {
	    print join (",", @$g), "\n";
	}
    }
    last;
}
    
sub get_genus {
    my $db = shift;
    my $node = shift;
    return () if ! $node;
    if( $node->rank eq 'genus' ) {
	return [$node->ncbi_taxid, $node->name('scientific')];
    } 
    my @results;
    for my $n ( $db->each_Descendent($node) ) {
	push @results, &get_genus($db,$node);
    }
    @results;
}
