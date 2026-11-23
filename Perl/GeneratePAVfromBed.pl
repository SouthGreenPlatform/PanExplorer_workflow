#!/usr/bin/perl

use strict;

use Graph::Undirected;
use FindBin;
use lib $FindBin::Bin;
use GeneId qw(gff_is_eukaryotic gene_key_from_attributes);

my $strain_info = $ARGV[0];
my $bed_directory = $ARGV[1];
my $PAV = $ARGV[2];
my $threshold_coverage = $ARGV[3];

my %strains;
open(S,$strain_info);
while (my $line = <S>) {
	chomp $line;
	my ($id,$name) = split(/\t/,$line);
	$strains{$id} = $name;
}
close(S);

# BED of the genes projected on the graph for a genome (exact name first, see above)
sub bed_of {
	my ($genome) = @_;
	my $bed = "$bed_directory/$genome.gene_segments.bed";
	if (-e $bed){return $bed;}
	my @candidates = glob("$bed_directory/$genome.*.bed");
	return $candidates[0];
}

my %gene_lengths;
foreach my $id(keys(%strains)){

	# exact file name first: "ls $id*.gff" also matches other genomes whose name starts with $id
	my $gff = "$bed_directory/$id.gb.gff";    # prokaryotes (same file as in create_gene_paths)
	if (!-e $gff){$gff = "$bed_directory/$id.gff";}  # eukaryotes
	if (!-e $gff){
		my @candidates = glob("$bed_directory/$id*.gff");
		$gff = $candidates[0];
	}
	# same gene identifiers as in CreateGenePathsFromGFA.pl (see GeneId.pm)
	my $is_eukaryotic = gff_is_eukaryotic($gff);
	open(GFF,$gff);
	my $current_gene_length;
	while(my $line = <GFF>){
		chomp($line);
		my @infos = split(/\t/,$line);
		# gene lengths are computed from CDS only (summed over CDS lines sharing the same identifier);
		# "gene" features are not used anymore: with bp_genbank2gff3 their ID is the locus_tag,
		# which would add the gene length to the CDS length of the same key.
		if ($infos[2] eq "CDS"){
			my $gene = gene_key_from_attributes($infos[8], $is_eukaryotic);
			if (!defined $gene){next;}
			my $start = $infos[3];
			my $end = $infos[4];
			my $genelength = $end-$start;
			$gene_lengths{$id}{$gene} += $genelength;
		}
	}
	close(GFF);
}




my %genes;
my $graph = Graph::Undirected->new;
my $num = 0;
foreach my $id(keys(%strains)){

	my $file = bed_of($id);
	open(B,$file);
	while(my $line = <B>){
		chomp($line);
		my @infos = split(/\t/,$line);
		my $gene1 = $infos[3];
		$genes{"$id!$gene1"} = 1;
	}
	close(B);

	foreach my $id2(keys(%strains)){
		my $file2 = bed_of($id2);
		if (-e $file && -e $file2 && $file ne $file2){
			$num++;

   			system("awk {'if (\$3>10)print $_'} $file >$file.light");
      			system("awk {'if (\$3>10)print $_'} $file2 >$file2.light");
			#system("bedtools intersect -a $file -b $file2 -wo >$PAV.$num.intersect.out");
   			system("bedtools intersect -a $file.light -b $file2.light -wo >$PAV.$num.intersect.out");
			my %cumul_match;
			my %cumul_match2;
			open(INTER,"$PAV.$num.intersect.out");
			while(<INTER>){
				my $line = $_;
				$line =~s/\n//g;$line =~s/\r//g;
				my ($segment,$start_gene1,$end_gene1,$gene1,$segment,$start_gene2,$end_gene2,$gene2,$size_overlap) = split(/\t/,$line);
				$cumul_match{"$gene1;$gene2"}+=$size_overlap;
			}
			close(INTER);

			#unlink("$PAV.$num.intersect.out");

			foreach my $pair(keys(%cumul_match)){
				my $size1 = $cumul_match{$pair};
				my ($gene1,$gene2) = split(/;/,$pair);
				my $gene1_length = $gene_lengths{$id}{$gene1};
				my $gene2_length = $gene_lengths{$id2}{$gene2};
				# safety: unknown or null length (identifier not found in the GFF) -> no link
				if (!$gene1_length || !$gene2_length){next;}
				my $size_match_gene = $cumul_match{"$gene1;$gene2"};
				my $percentage_overlap_gene1 = ($size_match_gene/$gene1_length)*100;
				my $percentage_overlap_gene2 = ($size_match_gene/$gene2_length)*100;
				#print "$pair $size_match_gene $percentage_overlap_gene1 $percentage_overlap_gene2\n";
				if ($percentage_overlap_gene1 > $threshold_coverage && $percentage_overlap_gene2 > $threshold_coverage){
					$graph->add_edge("$id!$gene1","$id2!$gene2");
					delete($genes{"$id!$gene1"});
					delete($genes{"$id2!$gene2"});
				}
			}
		}
	}
}

open(OUT,">$PAV");
print OUT "ClutserID";
foreach my $id(sort keys(%strains)){
	my $name = $strains{$id};
	print OUT "\t".$name;
}
print OUT "\n";
my $clnum = 0;
my @cc = $graph->connected_components();
foreach my $component (@cc){
	$clnum++;
	print OUT $clnum;
        my @genes = @$component;
	my %h;
	foreach my $gene(@genes){
		my ($id,$genename) = split(/!/,$gene);
		$h{$id}.="$genename,";
	}	
	foreach my $id(sort keys(%strains)){
		my $ids = "-";
		if ($h{$id}){
			$ids = $h{$id};	
			chop($ids);
		}
		print OUT "\t".$ids;
	}
	print OUT "\n";
}
#######################
# add singletons
#######################
foreach my $gene(keys(%genes)){
	$clnum++;
	print OUT $clnum;
	my ($id,$genename) = split(/!/,$gene);
	my %h;
	$h{$id}.="$genename,";
	foreach my $id(sort keys(%strains)){
                my $ids = "-";
                if ($h{$id}){
                        $ids = $h{$id};
                        chop($ids);
                }
                print OUT "\t".$ids;
        }
        print OUT "\n";	
}
close(OUT);
