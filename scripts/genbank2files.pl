#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Bio::SeqFeature::Tools::Unflattener;

my $force = 0;
use constant MRNA => 'mRNA';
use constant GENE => 'gene';
use constant CDS  => 'CDS';
use constant EXON  => 'exon';
use constant MIN_LENGTH => 1_000;

my $write_fasta_separate = 1;
my $SRC = "genbank";

my $dir = "download";
my $odir = 'genomes';
my $debug = 0;
GetOptions(
    'v|debug|verbose!' => \$debug,
    'fasta!'  => \$write_fasta_separate,
    'force!'  => \$force,
    'd|dir:s' => \$dir,
    'o|out:s' => \$odir,
    's|src:s' => \$SRC,
    );

opendir(DIR, $dir) || die $!;
for my $d ( qw(GFF DNA) ) {
    mkdir("$odir/$d") if( ! -d "$odir/$d" );
}

my $ran_count = 0;
for my $family ( readdir(DIR) ) {    
    next if $family =~ /^\./;
    next unless -d "$dir/$family";
    opendir(SPECIES,"$dir/$family");    
    for my $species ( readdir(SPECIES) ) {
	next if $species =~ /^\./;
	my $gbkdir = "$dir/$family/$species/gbk";
	next unless -d $gbkdir;
	next if ( -f "$odir/GFF/$species.gff3" && ! $force);
	opendir(GBKFILES, $gbkdir);
	my $cdsnum = 1;
	my $first = 1;
	my (@ALL,$outfh,$outfa,@seqs);
	for my $file ( readdir(GBKFILES) ) {
	    next if $file =~ /^\./;
	    next unless $file =~ /(\S+)\.(gb[sk]|gbff)/;
	    $ran_count++;
	    if( $first ) {
		warn("species is $species\n");
		open($outfh, ">$odir/GFF/$species.gff3") || die $!;
		if( $write_fasta_separate ) {
		    $outfa = Bio::SeqIO->new(-format => 'fasta', 
					     -file   => ">$odir/DNA/$species.fasta");
		}
		print $outfh "##gff-version 3\n";
		print $outfh "#date ". localtime(time)."\n";
	    }
	    $first = 0;

	    my $fh;
	    if( $file =~ /\.gz$/) {
		open($fh => "zcat $gbkdir/$file |") || die $!;
	    } else {
		open($fh => "<$gbkdir/$file") || die $!;
	    }	    

	    #warn("$file\n");
	    my $seqio = Bio::SeqIO->new(-format => 'genbank',
					-fh     => $fh);

	    my %genes;
	    # generate an Unflattener object
	    my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;
	    $unflattener->error_threshold(1);
	    while( my $seq = $seqio->next_seq ) {
		warn("seq is ", $seq->display_id, " len=", $seq->length,"\n") if $debug;
		my @top_sfs = $seq->get_SeqFeatures;
		warn("n # features ", scalar @top_sfs, "\n") if $debug;
		
		if( $seq->length < MIN_LENGTH || ! defined $seq->seq ) {
		    # write each contig if no genes in it
#			grep { $_->primary_tag eq 'gene' } @top_sfs) {
		    next;
		}
		push @seqs, $seq;
		print $outfh  join("\t",
				   $seq->display_id,
				   'chromosome',
				   'scaffold',
				   1, $seq->length,
				   '.','.','.',
				   sprintf("ID=%s;Name=%s;Accession=%s.%d",
					   $seq->display_id, $seq->display_id, 
					   $seq->accession_number, 
					   $seq->seq_version)),"\n";

		# get top level unflattended SeqFeatureI objects
		eval {
		    $unflattener->unflatten_seq( -seq       => $seq,
						 -group_tag =>'locus_tag',
						 -use_magic => 1);
		}; 
		if( $@ ) {
			warn("error on seq $gbkdir/$file\n");
			die;
		}
		my $i = 0;
		foreach my $f (@top_sfs) {
		    next unless $f->primary_tag eq 'gene';
		    next if $f->has_tag('pseudo');
		    my $primarytag = $f->primary_tag;
		    # skip tRNAs?
		    next unless( $primarytag eq 'CDS' || 
				 $primarytag eq 'gene' || 
				 $primarytag eq 'mRNA');
		    my $genestr;
		    my ($min,$max,$strand) = ($f->start,$f->end, $f->strand);
		    my ($pname,%genexrefs, @genexrefs_a);

		    if( $f->has_tag('db_xref') ) { 
			for my $xref ( $f->get_tag_values('db_xref') ) {
			    my ($xref_src,$xref_id) = split(':',$xref);
			    $genexrefs{$xref_src} = $xref_id;
			    push @genexrefs_a, &escape($xref);
			} 
		    }
		    for my $ptag ( qw(locus_tag gene name) ) {
			if( $f->has_tag($ptag) ) {
			    ($pname) = $f->get_tag_values($ptag);
			    last;
			}
		    }
		    $pname = $genexrefs{GeneID} if( exists $genexrefs{'GeneID'} && 
						    ! defined $pname);
		    unless( defined $pname ) {
			warn("cannot find pname in ", $f->gff_string, "\n");
			last;
		    }

		    my %refs;
		    push @ALL, [$seq->display_id,
				$SRC,
				GENE, 
				$f->start,
				$f->end,
				'.',
				$f->strand < 0 ? '-' : '+',
				'.',
				{ 'ID' => "$pname.gene", 'Name'=>$pname} ];

		    if( @genexrefs_a ) {
			$ALL[-1]->[8]->{'Dbxref'} = join(",", @genexrefs_a);
		    }

		    my $mrnact = 1;
		    my @mrnas = $f->get_SeqFeatures;
		    for my $mRNA ( @mrnas ) {
			my $mrna_name = $pname;
			if( @mrnas > 1 ) {
			    $mrna_name = "$pname.$mrnact"; 
			}
			my $x = [$seq->display_id,
				 $SRC,
				 MRNA, 
				 $mRNA->start,
				 $mRNA->end,
				 '.',
				 $mRNA->strand < 0 ? '-' : '+',
				 '.',
				 { 'ID'     => $mrna_name,
				   'Name'   => "$pname.tr",
				   'Parent' => "$pname.gene"}];

			my %mRNAxref;
			my @m_xrefs;
			if( $mRNA->has_tag('db_xref') ) { 
			    for my $xref ( $mRNA->get_tag_values('db_xref') ) {
				my ($xref_src,$xref_id) = split(':',$xref);			    
				$mRNAxref{$xref_src} = $xref_id;
				push @m_xrefs, &escape($xref);
			    }
			}

			$x->[8]->{'Dbxref'} = join(",",@m_xrefs);
			my @note;
			if( $mRNA->has_tag('note') ) {
			    @note= $mRNA->get_tag_values('note');
			}
			if( $mRNA->has_tag('product') ) {
			    unshift @note, &escape( $mRNA->get_tag_values('product') );
			}		    

			$x->[8]->{'Note'} = join(",",&escape(@note));
			my @alias;
			for my $t ( qw(protein_id transcript_id synonym) ) {
			    if( $f->has_tag($t) ) {
				# warn("$t --> ", $f->get_tag_values($t), "\n");
				push @alias, &escape( $f->get_tag_values($t) );
			    }
			}
			if( @alias ) {
			    $x->[8]->{'Alias'} = join(",", @alias);
			}
			my (@exon_a,@cds_a);
			for my $CDS ( $mRNA->get_SeqFeatures ) {
			    if( $CDS->primary_tag eq 'CDS' ) {
				my (@cdsrefs_a,%cdsxrefs);
				if(! @cdsrefs_a && $CDS->has_tag('db_xref') ) {
				    for my $cxref ( $CDS->get_tag_values('db_xref') ) {
					my ($cxref_src,$cxref_id) = split(':',$cxref);
					$cdsxrefs{$cxref_src} = $cxref_id;
					push @cdsrefs_a, &escape($cxref);
				    } 
				}

				my $codon_start = 1;
				if( $CDS->has_tag('codon_start') ) {
				    ($codon_start) = $CDS->get_tag_values('codon_start');
				}
				my $cdsct = 1;
				for my $cdsloc ( sort { 
				    ($a->strand * $a->start) <=> 
					($b->strand * $b->start)} 
						 $CDS->location->each_Location ) {
				    my $id;
				    if( exists $cdsxrefs{'GI'} ) {
					$id = $cdsxrefs{'GI'};
				    }
				    if( ! defined $id ) {
					$id = sprintf("cds%05d",$cdsnum++);
				    } else {
					$id = sprintf("%s.cds%d",$id,$cdsct);
				    }

				    my ($cds_start,$cds_end) = ($cdsloc->start, 
								$cdsloc->end);
				    if( $codon_start > 1 && $cdsct == 1 ) {
					if( $cdsloc->strand < 0 ) {
					    $cds_end -= ($codon_start-1);
					} else {
					    $cds_start += ($codon_start-1);
					}
				    }
				    if( $cds_start == $cds_end ) {
					#    $cds_end+=1;
				    }
				    push @cds_a, [$seq->display_id,
						  $SRC,
						  CDS, 
						  $cds_start,
						  $cds_end,
						  '.',
						  $cdsloc->strand < 0 ? '-' : '+',
						  '.',
						  { 'ID' => $id,
						    'Parent' => $mrna_name}];		
				    if( @cdsrefs_a ) {
					$cds_a[-1]->[8]->{'Dbxref'} = join(",",@cdsrefs_a);
				    }

				    for my $copy ( qw(protein_id codon_start transl_table
						  inference translation) ) {
					if( ! exists $x->[8]->{$copy} && 
					    $CDS->has_tag($copy) ) {
					    ($x->[8]->{$copy}) = $CDS->get_tag_values($copy);
					}
				    }

				    $cdsct++;
				}
			    } else {
				push @exon_a, [$seq->display_id,
					       $SRC,
					       EXON, 
					       $CDS->start,
					       $CDS->end,
					       '.',
					       $CDS->strand < 0 ? '-' : '+',
					       '.',
					       { 'Parent' => $mrna_name}];
			    }
			}
			push @ALL, $x;
			push @ALL, @exon_a, @cds_a;
			$mrnact++;
		    }
		}		
	    }
	    if ( my @ps = $unflattener->get_problems ) {
		warn("Problems in file: $file \n");
		for my $p ( @ps ) {
		    warn("Severity: ",$p->[0], "\n",
			 $p->[1],"\n");
		}
	    }
	}
	warn("Total number of features are ", scalar @ALL,"\n");
	for my $feature ( @ALL ) {
	    my $lastcol = pop @$feature;
	    my @lastcol_n;
	    for my $field ( qw(ID Parent) ) {
		if( exists $lastcol->{$field} ) {
		    push @lastcol_n, "$field=".$lastcol->{$field};
		    delete $lastcol->{$field};
		}
	    }	
	    for my $k ( keys %{$lastcol}) {
		push @lastcol_n, sprintf("%s=%s",$k,$lastcol->{$k});
	    }
	    print $outfh join("\t", @$feature, join(";", @lastcol_n)),"\n";
	}
	if( $write_fasta_separate && $outfa ) {
	    $outfa->write_seq(@seqs);
	} elsif( $outfh && @ALL && @seqs) {
	    my $out = Bio::SeqIO->new(-format => 'fasta',
				      -fh     => $outfh);
	    $out->write_seq(@seqs);
	}
    }
    last if $debug && $ran_count;
}


sub uniqlst { 
    my %x;
    return grep { ! $x{$_}++}  @_;
}

sub escape {
    
    for my $value ( @_) {
	if(  defined $value && length($value) ) { 
	    if ($value =~ /[^a-zA-Z0-9\,\;\=\.:\%\^\*\$\@\!\+\_\?\-]/) {
		$value =~ s/\t/\\t/g;       # substitute tab and newline 
		# characters
		$value =~ s/\n/\\n/g;       # to their UNIX equivalents
		
# Unescaped quotes are not allowed in GFF3
#                   $value = '"' . $value . '"';
	    }
	    $value =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;
	} 
    }
    return @_;
}
