#!/usr/bin/perl

use strict;
use FindBin;
use lib $FindBin::Bin;
use GeneId qw(gff_is_eukaryotic gene_key_from_attributes);



my $gfa = $ARGV[0];
my $gff3 = $ARGV[1];
my $out = $ARGV[2];
# optional: name of the genome in the graph (PanSN prefix "sample#hap#chrom").
# When given, the genes of this GFF are only projected on the paths of this genome. This is required
# when several genomes share the same chromosome names (Chr1, Chr2... in eukaryotes): otherwise the
# coordinates of the genes would also be applied to the homonymous chromosomes of the other genomes.
my $sample_filter = $ARGV[3];

my %segment_is_covered_by_genes;
my %coordinates;
my %store_gene_paths;

my %allgenes;
my %coding_positions;
my %gene_length;
my @gene_order;

# Gene identifier: locus_tag for prokaryotic annotations (unique per gene, keeps multi-copy
# genes distinct), protein_id/Parent/ID for eukaryotic ones (see GeneId.pm).
my $is_eukaryotic = gff_is_eukaryotic($gff3);

open(F,$gff3);
while(my $line = <F>){
		chomp($line);
		$line =~s/\r//g;
		my @infos = split(/\t/,$line);
		#if ($infos[2] eq "gene" && /CITME_006g014440/){
		if ($infos[2] eq "CDS"){
			my $gene = gene_key_from_attributes($infos[8], $is_eukaryotic);
			if (!defined $gene){next;}
			my $chr = $infos[0];
			my $start = $infos[3];
			my $end = $infos[4];
			$coordinates{$chr}.= "$gene:$start-$end|";
			$allgenes{$gene}=1;
			for (my $i = $start; $i <= $end; $i++){
				# several genes can share a position (overlapping genes, isoforms): keep all of them
				my $current = $coding_positions{$chr}{$i};
				if (!defined $current){
					$coding_positions{$chr}{$i} = "$gene,";
				}
				elsif (index(",$current", ",$gene,") < 0){
					$coding_positions{$chr}{$i} .= "$gene,";
				}
			}
			# CDS split over several lines (join, frameshift) share the same identifier: sum their lengths
			push(@gene_order, $gene) if (!exists $gene_length{$gene});
			$gene_length{$gene} += $end-$start;
		}
}
close(F);

open(LENGTH,">$out.gene_length.txt");
foreach my $gene(@gene_order){
	print LENGTH "$gene	$gene_length{$gene}\n";
}
close(LENGTH);

# chromosome names of the paths may lack the version suffix (".1") of the GFF seqids, or the reverse
my %chr_alias;
foreach my $chr(keys(%coordinates)){
	my $short = $chr;
	$short =~s/\.\d+$//;
	$chr_alias{$short} = $chr if (!exists $chr_alias{$short});
}
my %chr_projected;

print "Start GFA parsing ($gfa)...\n";

my %segments_chr;
my %segments_of_gene;
my %segments;
#open(ZCAT,"zcat $gfa |");
open(GFA,$gfa);
while(my $line = <GFA>){
	chomp($line);
	if ($line =~/^S/){
		my ($type,$segment,$segment_sequence,$SN,$tag) = split(/\t/,$line);
		$segments{$segment} = $segment_sequence;
	}
	elsif ($line =~/^P/ or $line =~/^W/){

		my $path;
		my $chrom;
		my $path_sample;
		my @segments;
		if ($line =~/^W/){
	                my ($type,$sample,$strand,$chromW,$start,$length,$pathW) = split(/\t/,$line);
			$path = $pathW;
			$chrom = $chromW;
			$path_sample = $sample;
			@segments = split(/[\>\<]/,$path);
		}
		elsif ($line =~/^P/){
			my ($type,$strain_chrom,$pathP) = split(/\t/,$line);
			my ($stra,$haplo,$chromP) = split(/#/,$strain_chrom);
			if (!defined $chromP){$chromP = $strain_chrom;}  # path name without PanSN prefix
			$path = $pathP;
			$chrom = $chromP;
			$path_sample = $stra;
			@segments = split(/,/,$path);
		}
		if (defined $sample_filter && $sample_filter ne "" && $path_sample ne $sample_filter){next;}
		my $short_chrom = $chrom;
		$short_chrom =~s/\.\d+$//g;
		if (exists $coordinates{$chrom}){}
		elsif (exists $chr_alias{$short_chrom}){$chrom = $chr_alias{$short_chrom};}
		else {$chrom = $short_chrom;}
		if ($coordinates{$chrom}){
			$chr_projected{$chrom}++;
			my $end_segment = 0;
			my $start_segment = 0;
			my $nb_gene = 0;
			
			my %genes_located;
			my %first_segment_of_gene;
			print "$chrom: ". scalar @segments."\n";
			my $cumul = 0;
			LOOP_SEGMENT: foreach my $segment(@segments){
				if ($segment eq ""){next;}
				my $segment_including_strand = $segment;
				$segment =~s/\+//g;
				$segment =~s/\-//g;
				my $size = length($segments{$segment});
				$start_segment = $end_segment+1;
				$end_segment = $start_segment+$size-1;
				$cumul+=$size;
				my $current_gene;

				if ($segment_including_strand =~/\+/){
					my $pos_in_segment = 0;
					for (my $k = $start_segment; $k <= $end_segment; $k++){
						$pos_in_segment++;
						#print "$k\n";
					
						if ($coding_positions{$chrom}{$k}){
							my $info_gene = $coding_positions{$chrom}{$k};
							chop($info_gene);
							my @genes = split(/,/,$info_gene);
							foreach my $gene(@genes){
								$segments_of_gene{$gene} .= "s$segment:$pos_in_segment|";
							}
						}				
					}
				}
				elsif ($segment_including_strand =~/\-/){
					my $pos_in_segment = $size+1;
					for (my $k = $start_segment; $k <= $end_segment; $k++){
						$pos_in_segment--;
						if ($coding_positions{$chrom}{$k}){
							my $info_gene = $coding_positions{$chrom}{$k};
							chop($info_gene);
							my @genes = split(/,/,$info_gene);
							foreach my $gene(@genes){
								$segments_of_gene{$gene} .= "s$segment:$pos_in_segment|";
							}
						}
					}
				}
			
			}
			print "$chrom $cumul $end_segment\n";
		}
	}
}
close(GFA);

# summary: genes of the GFF that could not be projected on the graph
my $nb_genes = scalar keys(%allgenes);
my $nb_projected = scalar keys(%segments_of_gene);
print STDERR "[CreateGenePathsFromGFA] ".($sample_filter ? "$sample_filter: " : "")."$nb_projected / $nb_genes genes projected on the graph\n";
foreach my $chr(sort keys(%coordinates)){
	if (!$chr_projected{$chr}){
		print STDERR "[CreateGenePathsFromGFA] WARNING: no graph path found for GFF sequence '$chr'".($sample_filter ? " of $sample_filter" : "")."\n";
	}
	elsif ($chr_projected{$chr} > 1){
		print STDERR "[CreateGenePathsFromGFA] WARNING: GFF sequence '$chr' matches $chr_projected{$chr} graph paths: give the genome name as 4th argument\n";
	}
}

open(O,">$out");
open(BED,">$out.bed");
foreach my $gene(keys(%segments_of_gene)){
	my $segments = $segments_of_gene{$gene};
	my @positions = split(/\|/,$segments);
	my $previous_segment;
	
	#print O "$gene	$segments\n";
	print O "$gene\t";
	
	my $concat = "";
	foreach my $pos(@positions){
		
		my ($segment,$pos_in_segment) = split(/:/,$pos);
		if ($segment ne $previous_segment){
			$concat .= ">$segment:$pos_in_segment";
		}
		$concat .= "-$pos_in_segment";
		$previous_segment = $segment;
	}
	my @segments = split(/>/,$concat);
	foreach my $segment(@segments){
		if (!$segment){next;}
		my ($segmentname,$pos) = split(/:/,$segment);
		my @positions = split(/-/,$pos);
		my $start = $positions[0];
		my $end = $positions[$#positions];
		if ($start > $end){
			$start = $positions[$#positions];
			$end = $positions[0];
		}
		print O ">$segmentname:$start-$end";
		print BED "$segmentname	$start	$end	$gene\n";
	}
	print O "\n";
}
close(O);
close(BED);
