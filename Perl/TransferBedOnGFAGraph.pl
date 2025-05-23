#!/usr/bin/perl

use strict;



my $gfa = $ARGV[0];
my $bed = $ARGV[1];
my $out = $ARGV[2];

my %segment_is_covered_by_genes;
my %coordinates;
my %store_gene_paths;

my %allgenes;
my %coding_positions;
open(LENGTH,">$out.gene_length.txt");
open(F,$bed);
while(my $line = <F>){
		chomp($line);
		my @infos = split(/\t/,$line);
		my $chr = $infos[0];
		my $start = $infos[1];
		my $end = $infos[2];
		my $feature = $infos[3];
		$coordinates{$chr}.= "$feature:$start-$end|";
		$allgenes{$feature}=1;
		for (my $i = $start; $i <= $end; $i++){
			$coding_positions{$chr}{$i} = "$feature,";
		}
}
close(F);
close(LENGTH);

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
		my @segments;
		if ($line =~/^W/){
	                my ($type,$sample,$strand,$chromW,$start,$length,$pathW) = split(/\t/,$line);
			$path = $pathW;
			$chrom = $chromW;
			@segments = split(/[\>\<]/,$path);
		}
		elsif ($line =~/^P/){
			my ($type,$strain_chrom,$pathP) = split(/\t/,$line);
			my ($stra,$haplo,$chromP) = split(/#/,$strain_chrom);
			$path = $pathP;
			$chrom = $chromP;
			@segments = split(/,/,$path);
		}
		$chrom =~s/\.\d+$//g;
		my @genes = split(/\|/,$coordinates{$chrom});
		if ($coordinates{$chrom}){
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
